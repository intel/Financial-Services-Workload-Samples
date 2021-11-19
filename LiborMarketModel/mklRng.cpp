/*
Copyright (c) 2019-2020 
Intel Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#include "Timer.hpp"
#include "params.h"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <cassert>
#include <vector>
#include <numeric>
#include <omp.h>
#include <mkl_vsl.h>
#include <sstream>

inline void genGaussian(VSLStreamStatePtr stream, int num, float* z,
                        float mean, float stddev)
{
  vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, num, z, mean, stddev);
}

inline void genGaussian(VSLStreamStatePtr stream, int num, double* z,
                        double mean, double stddev)
{
  vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, num, z, mean, stddev);
}

// define this for timing the RNG on all threads individually
//#define RNG_TIME

// define this to use blocking in RNG
#define BLOCKING

const int BLOCKSZ = 1024;


// parallel random number generation using MKL and OpenMP
void generateRand(REAL* z, int num, REAL mean, REAL stddev, int seed)
{
#ifdef RNG_TIME
  int nT;
  Timer t2;
  t2.start();

  // just to fire up OMP and find number of threads used
#pragma omp parallel
  {
    int id = omp_get_thread_num();
    if (id == 0)
      nT = omp_get_num_threads();
  }
  t2.stop();
  std::cout << "thread start time: " << t2.getTime() << "\n";
  std::cout << "using " << nT << " threads\n";
  double *ctimes = new double[nT];
  double *dtimes = new double[nT];
  double *stimes = new double[nT];
  double *gentimes = new double[nT];
  t2.start();
#endif
  #pragma omp parallel
  {
    int id = omp_get_thread_num();
    int nThreads = omp_get_num_threads();
    int numThis = num / nThreads;
    int start = id * numThis;
    if (id == nThreads-1)
    {
      // adjust last thread to take all the rest
      numThis = num - start;
    }
#ifdef BLOCKING
    const int blockSz = BLOCKSZ;
#else
    int blockSz = numThis;
#endif
    int nblocks = numThis / blockSz;
    int ntail = numThis - nblocks * blockSz;

#ifdef RNG_TIME
    Timer t;
    t.start();
#endif
    VSLStreamStatePtr stream;
    vslNewStream(&stream, VSL_BRNG_MRG32K3A, seed);
#ifdef RNG_TIME
    t.stop();
    ctimes[id] = t.getTime();
    t.start();
#endif
    vslSkipAheadStream(stream, start);  // offset for current thread
#ifdef RNG_TIME
    t.stop();
    stimes[id] = t.getTime();
    t.start();
#endif
    int iblock = 0;
    for (; iblock < nblocks; ++iblock)
      genGaussian(stream, blockSz, &z[start + iblock*blockSz], mean, stddev);
    if (ntail)
    {
      genGaussian(stream, ntail, &z[start+iblock*blockSz], mean, stddev);
    }
#ifdef RNG_TIME
    t.stop();
    gentimes[id] = t.getTime();
    t.start();
#endif
    vslDeleteStream(&stream);
#ifdef RNG_TIME
    t.stop();
    dtimes[id] = t.getTime();
#endif
  }
#ifdef RNG_TIME
  t2.stop();
  std::cout << "core RNG time: " << t2.getTime() << "\n";
  std::cout << "id,create,skip,gen,destroy\n";
  for (int i = 0; i < nT; ++i)
  {
    std::cout << i << ","
        << ctimes[i] << ","
        << stimes[i] << ","
        << gentimes[i] << ","
        << dtimes[i] << "\n";
  }

  delete [] ctimes;
  delete [] dtimes;
  delete [] gentimes;
  delete [] stimes;
#endif
}

