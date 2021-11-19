/*

Copyright (c) 2019-2020 
Intel Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#pragma once

#include <mkl_vsl.h>

#define REAL_BYTES 8

#if REAL_BYTES == 8
    typedef double REAL;  // floating-point type used throughout the code. Code has only been tested with double
#else
    typedef float REAL;
#endif

#define NBUF 1024   // Number of paths per path generation thread work-unit - numPaths must be a multiple of NBUF
#define VECLEN 16   // Vector length for SIMD - used in Longstaff-Schwartz

const int numPaths = 131072*2*2;  // must be a multiple of NBUF
const  int numSteps = 256;

//func declarations
inline double det3 (double a, double b, double c,
                    double d, double e, double f,
                    double g, double h, double i);
inline REAL determinant(REAL x[][3]);
void invert3_arr(REAL a[][3], REAL x[][3], int numRows, int numCols);

void perform_sims(REAL **a, VSLStreamStatePtr *streams, int nthreads, REAL fwdFactor, REAL volSqrtdt);
