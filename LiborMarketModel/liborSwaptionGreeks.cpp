/*
Copyright (c) 2019-2020 
Intel Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/


/** \file liborSwaptionGreeks.cpp
 * Implementation of Libor pricer
 */

#include <cstdlib> 
#include <cstdio> 
#include <cmath>
#include <string>
#include <cstring>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <omp.h>
#include "params.h"
#include "Timer.hpp"
#include "mklRng.h"

const int NLOOPS=10;//250;

/** Monte Carlo LIBOR path calculation.
 *
 * Forward path calculation storing data
 * for subsequent reverse path calculation
 */
static void path_calc_b1(REAL *z, REAL L[][SIMDVLEN], REAL L2[][SIMDVLEN], const REAL* lambda)
{
  int   i, n;
  REAL sqez[SIMDVLEN], lam[SIMDVLEN], con1[SIMDVLEN], v[SIMDVLEN], vrat[SIMDVLEN];
  REAL z1[SIMDVLEN];

  //memcpy(L2, L, NN*sizeof(REAL));
  #pragma omp simd
  for (int m =0; m<NN; m++)
  {
    for (int j=0; j<SIMDVLEN; j++) {
      L2[m][j] = L[m][j];
    }
    
  }
  
  for(n = 0; n < NMAT; n++)
  {

   # pragma omp simd  
     for (int j =0; j<SIMDVLEN; j++)
     {
       z1[j] = z[(j*NMAT)+n];
     
       sqez[j] = SQRT_DELTA * z1[j];
     
       v[j] =REAL(0); 
     }

    for (i=n+1; i<NN; i++)
    {
      REAL val = lambda[i-n-1];
  #pragma omp simd
      for (int j=0; j<SIMDVLEN; j++) {
	lam[j]  = val;
	con1[j] = DELTA * lam[j];
	v[j]   += con1[j] * L[i][j] / (ONE + DELTA * L[i][j]);
	vrat[j] = con1[j] * v[j] + lam[j] * (sqez[j] - (HALF * con1[j]));
	vrat[j] = std::exp(vrat[j]);
	L[i][j] = L[i][j] * vrat[j];

      // store these values for reverse path 
      L2[i+(n+1)*NN][j] = L[i][j];
      }


    }
  }
}


/* Reverse path calculation of DELTAs using stored data */
static void path_calc_b2(REAL L_b[][SIMDVLEN], const REAL L2[][SIMDVLEN], const REAL *lambda)
{
  int   i, n;
  REAL faci[SIMDVLEN] __attribute__((aligned(64))), v1[SIMDVLEN] __attribute__((aligned(64)));
  REAL lam[SIMDVLEN] __attribute__((aligned(64)));
  REAL L2_simd[SIMDVLEN] __attribute__((aligned(64)));
  REAL L2_simd1[SIMDVLEN] __attribute__((aligned(64)));

  for (n= NMAT-1; n >= 0; n--)
  {
    for (int j=0; j<SIMDVLEN; j++) {
      v1[j] = ZERO;
    }

    for (i = NN - 1; i > n; i--)
    {
      double val = lambda[i-n-1];
      #pragma omp simd
      for (int j=0; j<SIMDVLEN; j++) {
	lam[j] = val; 
	L2_simd[j] = L2[i+n*NN][j];
	L2_simd1[j] = L2[i+(n+1)*NN][j];
	v1[j]    += lam[j] * L2_simd1[j] * L_b[i][j];
	faci[j]   = DELTA / (ONE + DELTA * L2_simd[j]);
	L_b[i][j] = L_b[i][j] * L2_simd1[j] / L2_simd[j]
          + v1[j] * lam[j] * faci[j] * faci[j];
      }


    }
  }
}

/* calculate the portfolio value v */
static REAL* portfolio_b(REAL L[][SIMDVLEN], REAL v[SIMDVLEN]) 
{
  int   m, n;
  //static REAL v[SIMDVLEN]; 
  REAL b[SIMDVLEN], s[SIMDVLEN], swapval[SIMDVLEN];
  REAL B[NMAT][SIMDVLEN], S[NMAT][SIMDVLEN], B_b[NMAT][SIMDVLEN], S_b[NMAT][SIMDVLEN];

  for (int i=0; i<SIMDVLEN; i++) {
    b[i] = ONE;
    s[i] = ZERO;
    v[i] = ZERO; 
  }

  for (m = 0; m < NN - NMAT; m++)
  {
    n    = m + NMAT;
    for (int j=0; j<SIMDVLEN; j++) {
      b[j]    = b[j] / (ONE + DELTA * L[n][j]);
      s[j]    = s[j] + DELTA * b[j];
      B[m][j] = b[j];
      S[m][j] = s[j];
    }
    
  }

 

  //memset(B_b, 0, sizeof(B_b)*NMAT*SIMDVLEN);
  //memset(S_b, 0, sizeof(S_b)*NMAT*SIMDVLEN);
  for (int i =0; i < NMAT; i++)
  {
    for (int j=0; j<SIMDVLEN; j++) {
      B_b[i][j] = ZERO;
      S_b[i][j] = ZERO;
    }
   

  }

  for (n = 0; n < NOPT; n++)
  {
    m = maturities[n] - 1;
    for (int j=0; j<SIMDVLEN; j++) {
      swapval[j] = B[m][j] + swaprates[n]*S[m][j] - ONE;
      if (swapval[j] < ZERO) {
	  v[j]      += -HUNDRED * swapval[j];
	  S_b[m][j] += -HUNDRED * swaprates[n];
	  B_b[m][j] += -HUNDRED;
      }
    }
   
  }

  for (m = NN - NMAT - 1; m >= 0; m--)
  {
    n = m + NMAT;
    for ( int j=0; j<SIMDVLEN; j++) {
      B_b[m][j] += DELTA * S_b[m][j];
      L[n][j]  = -B_b[m][j] * B[m][j] * DELTA / (ONE + DELTA * L[n][j]);
      if (m > 0) {
	S_b[m-1][j] += S_b[m][j];
	B_b[m-1][j] += B_b[m][j] / (ONE + DELTA * L[n][j]);
      }
    }
   }

#pragma omp simd
    for (int j=0; j<SIMDVLEN; j++)
      b[j]=ONE;
    for (n = 0; n < NMAT; n++) {
      for (int j=0; j<SIMDVLEN; j++) {
	b[j] = b[j] / (ONE + DELTA * L[n][j]);
      }
    }

#pragma omp simd
    for ( int j=0; j<SIMDVLEN; j++)
      v[j]= b[j]*v[j];

    for (n = 0; n < NMAT; n++) {
      for (int j=0; j<SIMDVLEN; j++) {
	L[n][j] = -v[j] * DELTA / (ONE + DELTA * L[n][j]);
      }
    }

    for (n = NMAT; n < NN; n++) {
      for ( int j=0; j<SIMDVLEN; j++) {
	L[n][j] = b[j] * L[n][j];
      }
    }
  return v;
}

/*
/// Compute num normally distributed random numbers using the standard
/// rand() function and the Box-Mueller transform
/// Note: rand() is the worst RNG possible
static void generateRand(REAL* z, int num, REAL mean, REAL stddev)
{
  for (int n = 0; n < num/2; n++)
  {
    //two uniform random numbers 
    REAL x; 
    REAL y;
    REAL r, d;

    do {
      x = REAL(2) * ((REAL)rand() / (REAL)RAND_MAX) - 1;
      y = REAL(2) * ((REAL)rand() / (REAL)RAND_MAX) - 1;
      r = x*x + y*y;
    } while (r == REAL(0) || r >= REAL(1) );
    d = std::sqrt(-REAL(2) * std::log(r) / r);
 	      
    z[2*n] = x*d;
    z[2*n+1] = y*d;
  }
}
*/

/* -------------------------------------------------------- */

void liborSwaptionGreeks(REAL *z, int num, REAL mean, REAL std_dev,
                         int seed, REAL *d_v, REAL *d_Lb, int numPaths)
{
  REAL L2[L2_SIZE] __attribute__((aligned(64))); for (size_t i=0; i<L2_SIZE; i++) {L2[i]=0.0;} //@
  int path;
  Timer t,tt_v, tt_lb;
  double vTime=0, lbTime=0;
  //@ const int blockSize= NMAT*SIMDVLEN;
  
  t.start();
  //std::vector<REAL> z(num);
  //REAL *z;
  //z = (REAL *) _mm_malloc(sizeof(REAL)*num, 64);
  
  REAL* ptrZ = &z[0];
  t.stop();
  double vecSetupTime = t.getTime();
 
  t.start();
  //srand((unsigned)seed);
  generateRand(ptrZ, num, mean, std_dev);
  t.stop();
  double rngtime = t.getTime();
  //std::cout<< "RNG Done " << std::endl;

  
  t.start();
  // initialize constants
  REAL lambda[NN], L0[NN];
  std::fill(lambda, lambda+NN, REAL(0.2));
  std::fill(L0, L0+NN, REAL(0.051));
  t.stop();
  double initConstTime = t.getTime();

  t.start();
  double sumv = REAL(0);
  double sumlb = REAL(0);
for(int loops=0;loops<NLOOPS; loops++)
{
  //#pragma simd reduction(+:sumv) reduction(+:sumlb)
  #pragma omp parallel for reduction(+:sumv) reduction(+:sumlb)
  for (path=0; path<numPaths; path+=SIMDVLEN)
  {
   //#pragma simd reduction(+:sumv) reduction(+:sumlb)
   //for (int pathInn=0; pathInn < BLOCK; pathInn+=SIMDVLEN)
   {   
   //@ int tid = omp_get_thread_num();
    REAL L[NN][SIMDVLEN] __attribute__((aligned(64)));
    //REAL L_b[NN][SIMDVLEN];
    //REAL lambda1[NN][SIMDVLEN];
    REAL sumv_local[SIMDVLEN];
    //@ REAL sumlb_local[SIMDVLEN];
    REAL L2_local[L2_SIZE][SIMDVLEN] __attribute__((aligned(64)));
    //__attribute__(aligned(64)) REAL z[blockSize];

    //int streamAddr = (tid)*path/SIMDVLEN;
    //VSLStreamStatePtr stream = *(streams + streamAddr ); 
    //vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, blockSize, z, 0.0, 1.0);
    // intialize forward rates
    //memcpy(L, L0, NN*sizddeof(REAL));


    for(int i =0; i< NN; i++)
    {
      for (int j=0; j<SIMDVLEN; j++) {

	L[i][j]=L0[i];
      }

   
    }


    for (int i =0; i<L2_SIZE; i++)
    {
      for (int j=0; j<SIMDVLEN; j++) {

       L2_local[i][j] = L2[i];
      }

      /*       L2_local[i][:] = L2[i]; */
    }
    
    //REAL **L_b = &L[0];  //How do I use L to initialize L_b //This might be causing lambda error

    /* LIBOR path calculation and adjoint for Greeks */
    //tt_v.start();
    path_calc_b1(&ptrZ[path*NMAT], L, L2_local, lambda);
    portfolio_b(L,sumv_local);
    /*    sumv += __sec_reduce_add(sumv_local[:]);*/

    for (int j=0; j<SIMDVLEN; j++) {
      sumv += sumv_local[j];
    }
    //t_v.stop();
    //vTime += tt_v.getTime();
    
    //tt_lb.start();
    path_calc_b2(L, L2_local, lambda);
    //sumlb += L_b[NN-1];
    /*      sumlb += __sec_reduce_add(L[NN-1][:]); */
    for (int j=0; j<SIMDVLEN; j++) {
      sumlb += L[NN-1][j];
    }
    //tt_lb.stop();
    //lbTime += tt_lb.getTime();

    /* move pointer to start of next block */
    //ptrZ += (NMAT*SIMDVLEN);
   }
  }

}//Loops loop End

  *d_v = REAL(sumv / numPaths);
  *d_Lb = REAL(sumlb / numPaths);
  t.stop();
  double proctime = t.getTime();

  std::cout << "Vector setup time: " << vecSetupTime << std::endl;
  std::cout << "RNG time         : " << rngtime << std::endl;
  std::cout << "Init time        : " << initConstTime << std::endl;
  std::cout << "Proc time        : " << proctime << std::endl;
  std::cout << "vTime time       : " << vTime << std::endl;
  std::cout << "lbTime time      : " << lbTime << std::endl;
std::cout << "Number of Loops    : " << NLOOPS << std::endl;

}
