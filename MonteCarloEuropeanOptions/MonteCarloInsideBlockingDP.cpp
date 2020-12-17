/*
Copyright (c) 2019-2020 
Intel Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/


#include <stdlib.h>
#include <string.h>
#include <cstdio>
#include <cstring>
#include <cstdarg>
#include <cerrno>
#include <ctime>
#include <tbb/scalable_allocator.h>
#include <omp.h>

#ifdef __INTEL_COMPILER
#include <mathimf.h>
#endif

#include <mkl_vsl.h>
#include <mkl_service.h>

/******* VERSION *******/
#define MAJOR 1
#define MINOR 2


#define MAX_THREADS 288
#define RANDSEED 777

#define NUM_ITERS 4

#include <stddef.h>
#include <sys/time.h>
double second()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}

static inline double RandFloat_T(double low, double high, unsigned int *seed){
    double t = (double)rand_r(seed) / (double)RAND_MAX;
    return (1.0 - t) * low + t * high;
}

static inline void die(const char *fmt, ...)
{
    va_list val;
    va_start(val, fmt);
    vfprintf(stderr, fmt, val);
    va_end(val);
    exit(EXIT_FAILURE);
}

static const int SIMDALIGN = 4096;  //alignment requirement SSE: 16, AVX: 32, LRBni: 64

static int         NTHREADS         = 0;
static int         OPT_PER_THREAD   = 0;
static int         OPT_N            = 0;
static int         RAND_N           = 0;
static int         RAND_BLOCK_LENGTH= 0;
static int         VERBOSE          = 0;
static double       F_RAND_N         = 0;
static double       INV_RAND_N       = 0;
static double       STDDEV_DENOM     = 0;
static double       CONFIDENCE_DENOM = 0;
static const double RISKFREE         = 0.06;
static const double VOLATILITY       = 0.10;

static const double  RLog2E	= -RISKFREE * M_LOG2E;
static const double MuLog2E	= M_LOG2E*(RISKFREE - 0.5 * VOLATILITY * VOLATILITY);
static const double VLog2E	= M_LOG2E * VOLATILITY;

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

static long long suffixed_atoll(const char *nptr)
{
    char   *mod = strdup(nptr);
    size_t  s   = strlen(mod);
    long long res;
    switch(mod[s-1])
    {
    case 'k':
        res = 1024;
        break;
    case 'K':
        res = 1000;
        break;
    case 'm':
        res = 1024*1024;
        break;
    case 'M':
        res = 1000000;
        break;
    case 'g':
        res = 1024*1024*1024;
        break;
    case 'G':
        res = 1000000000;
        break;
    default:
        res = 1;
    }
    mod[s]  = 0;
    res    *= atoll(mod);
    free(mod);
    return res;
}

int main(int argc, char* argv[])
{
    double
      mkl_sTime, mklTime= 0.0, sTime, eTime;

    double sum_delta  = 0.0;
    double sum_ref    = 0.0;
    double max_delta  = 0.0;
    double sumReserve = 0.0;

    const double SQRT_2 = sqrt(2.);

    if(argc < 5)
        die("Usage: %s <nthreads> <noptions> <path_length> <path_block_length> <verbose(default=0)> \n", argv[0]);
    NTHREADS = suffixed_atoll(argv[1]);
    if(NTHREADS < 1)
        die("NTHREADS < 1: (%d)\n", NTHREADS);
    else if(NTHREADS > 288)
        die("NTHREADS > MAX_THREADS: (%d)\n", NTHREADS);
    if(((NTHREADS % 4) != 0) && (NTHREADS != 1))
        die("NTHREADS really ought to be a multiple of 4! (%d)\n", NTHREADS);

    OPT_N = suffixed_atoll(argv[2]);
    OPT_PER_THREAD = OPT_N/NTHREADS;
    if(OPT_PER_THREAD < 1)
        die("OPTS_N should be at least 1! (%d)\n", OPT_N);


    RAND_N              = suffixed_atoll(argv[3]);
    if(RAND_N < 16)
        die("RAND_N should be at least 16 for SIMD's sake! (%d)\n", RAND_N);
    F_RAND_N            = (double)RAND_N;
    INV_RAND_N          = 1.0/(double)RAND_N;
    STDDEV_DENOM        = 1.0/((double)RAND_N * (double)(RAND_N - 1));
    CONFIDENCE_DENOM    = 1.96/sqrt((double)RAND_N);

    RAND_BLOCK_LENGTH  = suffixed_atoll(argv[4]);
    if(RAND_BLOCK_LENGTH < 16)
        die("RAND_BLOCK_LENGTH should be at least 16! (%d)\n", RAND_BLOCK_LENGTH);
    else if(RAND_BLOCK_LENGTH > RAND_N)
        die("RAND_BLOCK_LENGTH should be no more than RAND_N(=%d)! (%d)\n", RAND_N, RAND_BLOCK_LENGTH);

    if(RAND_BLOCK_LENGTH % 16 != 0)
        die("RAND_BLOCK_LENGTH should divisibly by 16! (/16 = %lf)\n", RAND_BLOCK_LENGTH/16.0);
    if(RAND_N % RAND_BLOCK_LENGTH != 0)
        die("RAND_BLOCK_LENGTH should evenly divide RAND_N! (RAND_N/BLOCK_LENGTH = %lf)\n", RAND_N/(double)RAND_BLOCK_LENGTH);

    if(argc == 6)
       VERBOSE = suffixed_atoll(argv[5]);


    printf("Double Precision Monte-Carlo European Option Pricing, version %d.%d\n\n", MAJOR, MINOR);
    printf("Build Time       = %s %s\n", __DATE__, __TIME__);
    printf("Path Length      = %d\n", RAND_N);
    printf("Number of Options= %d\n", OPT_N);
    printf("Block Size       = %d\n", RAND_BLOCK_LENGTH);
    printf("Worker Threads   = %d\n\n", NTHREADS);

    const int mem_size  = sizeof(double)*OPT_PER_THREAD;

#ifdef _OPENMP
    omp_set_num_threads(NTHREADS);
#else
    int NTHREADS = 1;
#endif

    // calculate the block number based on block size
    const int nblocks = RAND_N/RAND_BLOCK_LENGTH;

    for(int ITER=0; ITER<NUM_ITERS ;ITER++)
    {
    
#pragma omp parallel reduction(+ : sum_delta) reduction(+ : sum_ref) reduction(+ : sumReserve) reduction(max : max_delta)
        {
#ifdef _OPENMP
            int threadID = omp_get_thread_num();
#else
            int threadID = 0;
#endif
            unsigned int randseed = RANDSEED + threadID;
            srand(randseed);
            double *CallResultList     = (double *)scalable_aligned_malloc(mem_size, SIMDALIGN);
            double *CallConfidenceList = (double *)scalable_aligned_malloc(mem_size, SIMDALIGN);
            double *StockPriceList     = (double *)scalable_aligned_malloc(mem_size, SIMDALIGN);
            double *OptionStrikeList   = (double *)scalable_aligned_malloc(mem_size, SIMDALIGN);
            double *OptionYearsList    = (double *)scalable_aligned_malloc(mem_size, SIMDALIGN);
        
            float *samples;
            VSLStreamStatePtr rngStream;
    
            for(int i = 0; i < OPT_PER_THREAD; i++)
            {
                CallResultList[i]     = 0.0;
                CallConfidenceList[i] = 0.0;
                StockPriceList[i]     = RandFloat_T(5.0, 50.0, &randseed);
                OptionStrikeList[i]   = RandFloat_T(10.0, 25.0, &randseed);
                OptionYearsList[i]    = RandFloat_T(1.0, 5.0, &randseed);
            }
            samples = (float *)scalable_aligned_malloc(RAND_BLOCK_LENGTH * sizeof(float), SIMDALIGN);
        
            vslNewStream(&rngStream, VSL_BRNG_SFMT19937, RANDSEED + threadID); 
    
    #pragma omp barrier
            if (threadID == 0)
            {
                sTime = second();
            }
    
            for(int opt = 0; opt < OPT_PER_THREAD; opt++)
            {
                const double VBySqrtT = VLog2E * sqrt(OptionYearsList[opt]);
                const double sqrt_2_sigma = SQRT_2*VBySqrtT;
                const double MuByT    = MuLog2E * OptionYearsList[opt];
                const double Y        = StockPriceList[opt];
                const double Z        = OptionStrikeList[opt];
                    
                double v0 = 0.0;
                double v1 = 0.0;
                for(int block = 0; block < nblocks; ++block)
                {
                    vsRngGaussian (VSL_RNG_METHOD_GAUSSIAN_ICDF, rngStream, RAND_BLOCK_LENGTH, samples, MuByT, VBySqrtT);
    
    
    #pragma omp simd reduction(+:v0) reduction(+:v1)
    #pragma unroll(4)
                    for(int i=0; i < RAND_BLOCK_LENGTH; i++) 
                    {
                        double rngVal = samples[i];
                        double callValue  = Y * exp2(rngVal) - Z;
        
                        if (callValue>0.0)
                        {
                           v0 += callValue;
                           v1 += callValue * callValue;
                        }
                    }
                } //end of block
                const double  exprt    = exp2(RLog2E*OptionYearsList[opt]);
                CallResultList[opt]     = exprt * v0 * INV_RAND_N;
                const double  stdDev   = sqrt((F_RAND_N * v1 - v0 * v0) * STDDEV_DENOM);
                CallConfidenceList[opt] = (double)(exprt * stdDev * CONFIDENCE_DENOM);
            } //end of opt 
    
    #pragma omp barrier
            if (threadID == 0) {
                eTime = second();
            }
    
            if(VERBOSE){
                double delta = 0.0, ref = 0.0, L1norm = 0.0;
                int max_index = 0;
                double max_local  = 0.0;
                for(int i = 0; i < OPT_PER_THREAD; i++){
                    double callReference = 0.0;
                    BlackScholesRefImpl(callReference,
                                        StockPriceList[i],
                                        OptionStrikeList[i], 
                                        OptionYearsList[i],  
                                        RISKFREE, 
                                        VOLATILITY );
                    ref   = callReference;
                    delta = fabs(callReference - CallResultList[i]);
                    sum_delta += delta;
                    sum_ref   += fabs(ref);
                    if(delta > 1e-6)
                       sumReserve += CallConfidenceList[i] / delta;
                       max_local = delta>max_local? delta: max_local;
                }
                max_delta = max_local>max_delta? max_local: max_delta;
            }
            scalable_aligned_free(samples);
            scalable_aligned_free(CallResultList);
            scalable_aligned_free(CallConfidenceList);
            scalable_aligned_free(StockPriceList);
            scalable_aligned_free(OptionStrikeList);
            scalable_aligned_free(OptionYearsList);
        
            vslDeleteStream(&rngStream);
        }//end of parallel block
    }//end of iter loop
    
    if(VERBOSE){
        sumReserve          /= (double)OPT_N;
        const double L1norm  = sum_delta / sum_ref;
        printf("L1_Norm          = %E\n", L1norm);
        printf("Average RESERVE  = %f\n", sumReserve);
        printf("Max Error        = %E\n", max_delta);
        printf(sumReserve > 1.0f ? "Test passed\n" : "Test failed!\n");
    }
    
    printf("==========================================\n");
    printf("Time Elapsed = %lf\n", eTime-sTime);
    printf("Opt/sec      = %lf\n", OPT_N/(eTime-sTime));
    printf("==========================================\n");
    return 0;
}
