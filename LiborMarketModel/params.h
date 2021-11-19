/*
Copyright (c) 2019-2020 
Intel Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#ifndef LIBORSWAPTIONGREEKS_PARAMS_H_
#define LIBORSWAPTIONGREEK_PARAMS_H_
#ifdef USE_DOUBLE
typedef double REAL;
#define ONE 1.0
#define ZERO 0.0
#define HALF 0.5
#define HUNDRED 100.0
const int SIMDVLEN=16;
#else
// datatype to use
typedef float REAL;
#define ONE 1.0f
#define ZERO 0.0f
#define HALF 0.5f
#define HUNDRED 100.0f
const int SIMDVLEN=32;
#endif

#define NN 80          ///< number of forward LIBOR rates
#define NMAT 40        ///< random numbers needed to compute forward rates
#define L2_SIZE 3280   ///< NN*(NMAT+1) - stored forward/reverse rates interally
#define NOPT 15        ///< number of swaptions
#define SEED 1         ///< default seed for random number generator

#define DELTA ((REAL)0.25)
#define SQRT_DELTA ((REAL)0.5)
#define BLOCK 16

// maturities and swaprates for the swaptions in the portfolio
const int maturities[] = {4,4,4,8,8,8,20,20,20,28,28,28,40,40,40};
const REAL swaprates[] = {.045f,.05f,.055f,.045f,.05f,.055f,.045f,.05f,
                          .055f,.045f,.05f,.055f,.045f,.05f,.055f};


#endif
