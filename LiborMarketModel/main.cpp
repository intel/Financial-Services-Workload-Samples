
/*
Copyright (c) 2019-2020 
Intel Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

/** \file main.cpp
 * Prices a LIBOR swaption portfolio using a Monte-Carlo method
 *
 * Computes the portfolio value and Greek value using and adjoint approach
 *
 * The algorithm is taken from Mike Giles
 * (see http://people.maths.ox.ac.uk/gilesm/hpc).
 */

#include "params.h"
#include "liborSwaptionGreeks.hpp"
#include "Timer.hpp"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <ctime>
#include <mkl_vsl.h>
#include <omp.h>
#include <vector>
#include "mklRng.h"

// using 512K paths by default
#define NPATH 131072 //524288*8*2

int main(int argc, char** argv)
{ 
  /* Parse Command-line Parameters ************************************/
  // first one is number of paths
  int numPath = NPATH;

  if (argc > 1)
  {
    std::stringstream sstr;
    sstr << argv[1];
    sstr >> numPath;
  }

  std::cout << "Number of paths: " << numPath << std::endl;
  
  /* Algorithm ****************************************************/

  REAL result = 0;
  REAL Lb = 0;
  int nthreads = 1;
  nthreads = omp_get_max_threads();
  #pragma omp parallel
  {
    

  }
  VSLStreamStatePtr *streams;
  int streamCount = nthreads*numPath/SIMDVLEN;
   
  streams = (VSLStreamStatePtr *) malloc(sizeof(VSLStreamStatePtr)*streamCount);
  //@ int nb = NMAT*SIMDVLEN;
  const int num = NMAT*numPath;
  std::vector<REAL> z(num);
  REAL* ptrZ = &z[0];
  generateRand(ptrZ, num, 0.0, 1.0);

  /* Initialize RNG */
/*
  #pragma omp parallel for
  for( int j=0 ; j<streamCount ; j++)
  {
      vslNewStream( streams + j, VSL_BRNG_MRG32K3A , 12345 );
      vslSkipAheadStream(streams[j], nb*j);  // offset for current thread     
 
   }
*/

  std::cout << "Running... " << std::endl;

  double elapsed;
  {
    ScopedTimer  st(elapsed);

    liborSwaptionGreeks(&z[0],NMAT*numPath, 0, 1, 1234,
                        &result, &Lb, numPath);

  }

  std::cout << "Result: value=" << result
                << ",\tlambda=" << Lb << std::endl;
  std::cout << "Total time: " << elapsed << " ms" << std::endl;
  std::cout << "Paths/sec: " << numPath/elapsed << " KPaths/sec " << std::endl;
  free(streams);
  return EXIT_SUCCESS;
}

