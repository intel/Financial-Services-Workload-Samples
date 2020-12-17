/*
Copyright (c) 2019-2020 
Intel Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#include <sys/time.h>
#include <time.h>

#include <stdlib.h>
#include <stdio.h>
#include <locale.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <tbb/scalable_allocator.h>

#ifndef __INTEL_COMPILER
#define __forceinline __attribute__((always_inline))
#endif

/******* VERSION *******/
#define MAJOR 1
#define MINOR 4

const int SIMDALIGN = 4096;

//Array size. Each kernel invocation calculates 2 options.
const int OPT_N = 65536;
const int  NUM_ITERATIONS = 32000;
#define CHUNKSIZE 256

const double      RISKFREE = 0.02;
const double    VOLATILITY = 0.30;
const double          HALF = 0.5;
const double RLOG2E =-RISKFREE*M_LOG2E; 
const double RVV =  (RISKFREE + 0.5 * VOLATILITY * VOLATILITY)/VOLATILITY;
const double RVLOG2E =  M_LN2 * (1.0/VOLATILITY);
const double VLOG2E =  M_LOG2E * VOLATILITY;
const double INV_LN2X2 = 1.0/(M_LN2 * 2.0);
const double VLOG2ERSQRT1_2 =  M_LOG2E * VOLATILITY / M_SQRT1_2;
const double RVVSQRT1_2  = RVV * M_SQRT1_2;
const double VOLATILITYSQRT1_2 = VOLATILITY * M_SQRT1_2;
const double SQRT1_2 = 0.70710678118654752440084436210485;

static double start_sec, end_sec; 
double second()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}

inline double RandDouble(double low, double high, unsigned int *seed){
    double t = (double)rand_r(seed) / (double)RAND_MAX;
    return (1.0 - t) * low + t * high;
}

#define inv_sqrt_2xPI 0.39894228040143270286
#define CNDF_C1 0.2316419
#define CNDF_C2 0.319381530
#define CNDF_C3 -0.356563782
#define CNDF_C4 1.781477937
#define CNDF_C5 -1.821255978
#define CNDF_C6 1.330274429


// Black-Scholes Reference
void BlackScholesRefImpl(
    double& callResult,
    double  Sf, //Stock price
    double  Xf, //Option strike
    double  Tf, //Option years
    double  Rf, //Riskless rate
    double  Vf  //Volatility rate
)
{
    // BSM Formula: https://www.nobelprize.org/prizes/economic-sciences/1997/press-release/
	double S = Sf, L = Xf, t = Tf, r = Rf, sigma = Vf;
	double N_d1 = 1. / 2. + 1. / 2. * erf(((log(S / L) + (r + 0.5 * sigma * sigma) * t) / (sigma * sqrt(t))) / sqrt(2.));
	double N_d2 = 1. / 2. + 1. / 2. * erf(((log(S / L) + (r - 0.5 * sigma * sigma) * t) / (sigma * sqrt(t))) / sqrt(2.));
	callResult = (S * N_d1 - L*exp(-r*t) * N_d2);
}


__forceinline void CNDF_C ( double * OutputX, double * InputX )
{
    double k2_2, k2_3, k2_4, k2_5;
    double tmp = 0.0;
    int flag=0;
    double    X = *InputX;
    
    if (X < 0.0) {
        X = -X;
        flag=1;
    }

    tmp = inv_sqrt_2xPI * exp2(-X*X*INV_LN2X2);

    double k2 = 1.0 / (1.0 + CNDF_C1 * X );
    k2_2 = k2*k2;
    k2_3 = k2_2*k2;
    k2_4 = k2_3*k2;
    k2_5 = k2_4*k2;
    
    *OutputX = 1.0 - ( tmp * ( (  CNDF_C2 * k2 ) +
                               ( ( CNDF_C3  * (k2_2) ) +
                                 ( CNDF_C4  * (k2_3) ) +
                                 ( CNDF_C5  * (k2_4) ) +
                                 ( CNDF_C6  * (k2_5) ) ) ) );
    if (flag)
        *OutputX = (1.0 - *OutputX);
        
    return;
}

int main(int argc, char* argv[])
{
	double sum_delta = 0.0, sum_ref = 0.0, L1norm = 0.0;
        unsigned int seed = 123;
	int verbose = 0; 
	if (argc > 2)
	{
		printf("usage: Black-Scholes <verbose> verbose = 1 for validtating result, the default is 0. \n");
        exit(1);
	}
	if (argc == 1)
	    verbose = 0;
	else if (argc == 2)
		verbose = atoi(argv[1]);
#ifdef _OPENMP
	int ThreadNum = omp_get_max_threads();
	omp_set_num_threads(ThreadNum); 
#else
	int ThreadNum = 1; 
#endif
	int OptPerThread = OPT_N / ThreadNum;
	int mem_size = sizeof(double) * OptPerThread; 
	setlocale(LC_ALL,"");
	printf("Double Precision Black-Scholes Formula, version %d.%d\n", MAJOR, MINOR);
	printf("Build Time        = %s %s\n", __DATE__, __TIME__);
	printf("Input Dataset     = %d\n", OPT_N);
	printf("Repetitions       = %d\n", NUM_ITERATIONS);
	printf("Chunk Size        = %d\n", CHUNKSIZE);
	printf("Worker Threads    = %d\n\n", ThreadNum);

	if (verbose)
		printf("Allocate and initialize memory on %d boundary,\n", SIMDALIGN);
#pragma omp parallel reduction(+ : sum_delta) reduction(+ : sum_ref)
{
#ifdef _OPENMP
	int threadID = omp_get_thread_num();
#else
	int threadID = 0; 
#endif
	double *CallResult = (double *)scalable_aligned_malloc(mem_size, SIMDALIGN);
	double *PutResult  = (double *)scalable_aligned_malloc(mem_size, SIMDALIGN);
	double *StockPrice    = (double *)scalable_aligned_malloc(mem_size, SIMDALIGN);
	double *OptionStrike  = (double *)scalable_aligned_malloc(mem_size, SIMDALIGN);
	double *OptionYears   = (double *)scalable_aligned_malloc(mem_size, SIMDALIGN);

	seed += threadID;
	for(int i = OptPerThread-1; i > -1 ; i--)
	{
		CallResult[i] = 0.0;
		PutResult[i]  = -1.0;
		StockPrice[i]    = RandDouble(5.0, 30.0, &seed);
		OptionStrike[i]  = RandDouble(1.0, 100.0, &seed);
		OptionYears[i]   = RandDouble(0.25, 10.0, &seed);
	}
#pragma omp barrier
	if (threadID == 0) {
		start_sec = second();
	}

	for(int i = 0; i < NUM_ITERATIONS; i++)
#pragma distribute_point
	    for (int chunkBase = 0; chunkBase < OptPerThread; chunkBase += CHUNKSIZE)
	    {
#pragma vector aligned
#pragma omp simd 
#pragma vector nontemporal (CallResult, PutResult)
		for(int opt = chunkBase; opt < (chunkBase+CHUNKSIZE); opt++)
		{
			double CNDD1, CNDD2; 
			double T = OptionYears[opt];
			double X = OptionStrike[opt];
			double S = StockPrice[opt];
			double sqrtT = sqrt(T);
			double d1 = (log(S / X) + (RISKFREE + 0.5 * VOLATILITY * VOLATILITY) * T) / (VOLATILITY * sqrtT);
			double d2 = d1 - VOLATILITY * sqrtT;
			CNDF_C (&CNDD1, &d1 );
			CNDF_C (&CNDD2, &d2 );
			double XexpRT = X*exp(-RISKFREE * T);
			double CallVal  = S * CNDD1 - XexpRT * CNDD2;
			double PutVal  = CallVal  +  XexpRT - S;
			CallResult[opt] = CallVal ;
			PutResult[opt] = PutVal ;
		}
	    }
#pragma omp barrier
	if (threadID == 0) {
        	end_sec = second();
	}

	if (verbose)
	{
       	double delta = 0.0, ref = 0.0, L1norm = 0.0;
		int max_index = 0;

		for(int i = 0; i < OptPerThread; i++)
		{
			double callReference = 0.0;
			BlackScholesRefImpl(callReference,
				                StockPrice[i],   //Stock price
				                OptionStrike[i], //Option strike
				                OptionYears[i],  //Option years
				                RISKFREE,        //Riskless rate
				                VOLATILITY       //Volatility rate
			);
			ref   = callReference;
			delta = fabs(callReference - CallResult[i]);
			sum_delta += delta;
			sum_ref   += fabs(ref);
		}
	}
	scalable_aligned_free(CallResult);
	scalable_aligned_free(PutResult);
	scalable_aligned_free(StockPrice);
	scalable_aligned_free(OptionStrike);
	scalable_aligned_free(OptionYears);
} //parallel section

	double sec = end_sec - start_sec;
	printf("=============================================\n");
	printf("Time Elapsed      = %5.4f sec\n", sec);
	printf("Throughput        = %5.4f GOptions/sec\n", (2.0f*NUM_ITERATIONS*OPT_N)/(1e9*sec));
	printf("=============================================\n");
	if (verbose)
	{
		L1norm = sum_delta / sum_ref;
		printf("L1 norm: %E\n", L1norm);
		printf((L1norm < 9e-5) ? "TEST PASSED\n" : "TEST FAILED\n");
	}
}
