//=======================================================================================
//
// SAMPLE SOURCE CODE - SUBJECT TO THE TERMS OF SAMPLE CODE LICENSE AGREEMENT,
// http://software.intel.com/en-us/articles/intel-sample-source-code-license-agreement/
//
// Copyright 2020 Intel Corporation
//
// THIS FILE IS PROVIDED "AS IS" WITH NO WARRANTIES, EXPRESS OR IMPLIED, INCLUDING BUT
// NOT LIMITED TO ANY IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
// PURPOSE, NON-INFRINGEMENT OF INTELLECTUAL PROPERTY RIGHTS.
//
// ======================================================================================


#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include"mkl_rng.h"
#include"timer.h"

#define RANDSEED 42
const int NUM_LOOPS = 2;
typedef double tfloat;
typedef double tfloat_rng;

const int SIMSTEPS = 128;
const int VECLEN = 64;
const int SIMD_ALIGN = 128;
#if 1
const int OPT_N = 100;
const int RAND_N = 131072;
#else
const int OPT_N = 10;
const int RAND_N = 10;
#endif

static const tfloat  RISKFREE = 0.03;
static const tfloat  VOLATILITY = 0.20;
static const tfloat  RHO = -0.65;
static const tfloat  INV_RAND_N = 1.0 / RAND_N;
static const tfloat  F_RAND_N = static_cast<tfloat>(RAND_N);
static const tfloat STDDEV_DENOM = 1.0 / (F_RAND_N * (F_RAND_N - 1.0));
static const tfloat CONFIDENCE_DENOM = 1.0 / std::sqrt(F_RAND_N);

void* (*my_aligned_malloc)(size_t, int) = &MKL_malloc;
void (*my_aligned_free)(void*) = &MKL_free;

//Simulate Spot Path
inline void monteCarloByTimeStepKernel(tfloat* prevStepVal, tfloat dt, tfloat_rng* z, tfloat* volArray, tfloat* nextStepVal)
{
#pragma omp simd simdlen(VECLEN)
#pragma vector aligned
    for (int iv = 0; iv < VECLEN; iv++)
    {
        tfloat Vmax = std::max(volArray[iv], 0.0);
        tfloat MuByT = (RISKFREE - 0.5 * Vmax * Vmax) * dt;
        tfloat VBySqrtT = Vmax * std::sqrt(dt);
        nextStepVal[iv] = prevStepVal[iv] * std::exp(MuByT + VBySqrtT * z[iv]);
    }
}

//Simulate Vol Path
inline void volTimeStepKernel(tfloat* prevStepVol, tfloat dt, tfloat_rng* z, tfloat kappa, tfloat theta, tfloat xi, tfloat* nextStepVol)
{
#pragma omp simd simdlen(VECLEN)
#pragma vector aligned
    for (int iv = 0; iv < VECLEN; iv++)
    {
        tfloat Vmax = std::max(prevStepVol[iv], 0.0);
        nextStepVol[iv] = prevStepVol[iv] + kappa * dt * (theta - Vmax) + xi * std::sqrt(Vmax * dt) * z[iv];
    }
}

void MonteCarloTimeStep(
    tfloat* h_CallResult,
    tfloat* h_CallConfidence,
    tfloat* S,
    tfloat* X,
    tfloat* T,
    tfloat* initialVol,
    tfloat* kappaVals,
    tfloat* thetaVals,
    tfloat* xiVals
)
{

    //Size of Random Numbers Needed
    int rngNums = VECLEN * SIMSTEPS;
    const int mem_size = sizeof(tfloat_rng) * rngNums;
    mkl_rng<tfloat_rng> stream1(VSL_BRNG_MT19937, RANDSEED, VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, 0.0, 1.0, mem_size, SIMD_ALIGN);
    mkl_rng<tfloat_rng> stream2(VSL_BRNG_MT19937, RANDSEED+1, VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, 0.0, 1.0, mem_size, SIMD_ALIGN);

    int location;

    double time = 0.0;
    double time_copy = 0.0;
    double time_steps = 0.0;
    double time_init = 0.0;
    double time_RNG = 0.0;

    __declspec(align(64)) tfloat simStepResult[SIMSTEPS + 1][VECLEN];
    __declspec(align(64)) tfloat volArray[SIMSTEPS + 1][VECLEN];     //Simulated Volatilities
    __declspec(align(64)) tfloat avgMean[VECLEN];
    __declspec(align(64)) tfloat callValue[VECLEN];

    timer<timer_enabled> tt;
    timer<timer_enabled> ttRng;

    tt.start();
    for (int iter = 0; iter < NUM_LOOPS; iter++)
    {
        for (int opt = 0; opt < OPT_N; opt++)
        {
            tfloat Sval = S[opt];
            tfloat Xval = X[opt];
            tfloat val = 0.0, val2 = 0.0;
            tfloat dt = T[opt] / SIMSTEPS;

            tfloat SqrtT = std::sqrt(dt);
            tfloat sqrtRhoSq = std::sqrt(1 - RHO * RHO);

            //Volatility simulation related
            tfloat v0 = initialVol[opt];
            tfloat kappa = kappaVals[opt];
            tfloat theta = thetaVals[opt];
            tfloat xi = xiVals[opt];

            for (int pos = 0; pos < RAND_N; pos += VECLEN)
            {

                ttRng.start();
                tfloat_rng* l_Random1 = stream1.get_gaussian(rngNums);
                tfloat_rng* l_Random2 = stream2.get_gaussian(rngNums);
                ttRng.stop();
                time_RNG += ttRng.duration();

                //#pragma omp simd
                for (int iv = 0; iv < VECLEN; iv++)
                {
                    simStepResult[0][iv] = Sval;
                    avgMean[iv] = 0.0;
                    callValue[iv] = 0.0;
                    volArray[0][iv] = v0;
                }

                //Simulate volatilities
                for (int simStep = 0; simStep < SIMSTEPS; simStep++)
                {
                    location = simStep * VECLEN;
                    //Correlate two random numbers
#pragma omp simd simdlen(VECLEN)
#pragma vector aligned
                    for (int i = location; i < location + VECLEN; i++)
                    {
                        l_Random2[i] = RHO * l_Random1[i] + sqrtRhoSq * l_Random2[i];
                    }
                    volTimeStepKernel(volArray[simStep], dt, &l_Random2[location], kappa, theta, xi, volArray[simStep + 1]);
                }

                //Simulate spot
                for (int simStep = 0; simStep < SIMSTEPS; simStep++)
                {
                    location = simStep * VECLEN;
                    monteCarloByTimeStepKernel(simStepResult[simStep], dt, &l_Random1[location], volArray[simStep], simStepResult[simStep + 1]);

                    //#pragma omp simd reduction(+:avgMean)
                    for (int iv = 0; iv < VECLEN; iv++)
                    {
                        avgMean[iv] += simStepResult[simStep + 1][iv];
                    }

                }

                //Use Arithmetic Mean
                //#pragma omp simd
                for (int iv = 0; iv < VECLEN; iv++)
                {
                    //avgMean[iv] = avgMean[iv]/SIMSTEPS;
                    callValue[iv] = std::max(((avgMean[iv] / SIMSTEPS) - Xval), 0.0);
                }

                tfloat tmpVal = 0.0;

                //#pragma omp simd reduction(+:tmpVal)
                for (int iv = 0; iv < VECLEN; iv++)
                {
                    tmpVal += callValue[iv];
                }
                val += tmpVal / VECLEN;
                val2 += val * val;

            }


            tfloat exprt = std::exp(-RISKFREE * T[opt]);
            h_CallResult[opt] = exprt * val * INV_RAND_N * VECLEN;
            tfloat stdDev = std::sqrt((F_RAND_N * val2 - val * val) * STDDEV_DENOM);

            h_CallConfidence[opt] = (exprt * 1.96 * stdDev * CONFIDENCE_DENOM);

        } //end of for
    }

    tt.stop();

    time = tt.duration() / (NUM_LOOPS);

    printf("Completed RNG %8.4f seconds \n", time_RNG / NUM_LOOPS);
    printf("Completed MC by timestep in %8.4f seconds \n", time);
    printf("Computation rate MC by time Step - %8.4f.\n", OPT_N / time);
    fprintf(stderr, "%8.4f\n", OPT_N / time);

}

int main(int argc, char* argv[])
{
    tfloat
        * CallResultParallel,
        * CallConfidence,
        * StockPrice,
        * OptionStrike,
        * OptionYears,
        //For heston Volatility
        * initialVol,
        * kappa,
        * theta,
        * xi;

    int
        i, mem_size, rand_size;

    //const int RAND_N = 1 << 18;

    mem_size = sizeof(tfloat) * OPT_N;
    rand_size = sizeof(tfloat) * RAND_N;

    CallResultParallel = (tfloat*)my_aligned_malloc(mem_size, SIMD_ALIGN);
    CallConfidence = (tfloat*)my_aligned_malloc(mem_size, SIMD_ALIGN);
    StockPrice = (tfloat*)my_aligned_malloc(mem_size, SIMD_ALIGN);
    OptionStrike = (tfloat*)my_aligned_malloc(mem_size, SIMD_ALIGN);
    OptionYears = (tfloat*)my_aligned_malloc(mem_size, SIMD_ALIGN);

    kappa = (tfloat*)my_aligned_malloc(mem_size, SIMD_ALIGN);
    theta = (tfloat*)my_aligned_malloc(mem_size, SIMD_ALIGN);
    xi = (tfloat*)my_aligned_malloc(mem_size, SIMD_ALIGN);
    initialVol = (tfloat*)my_aligned_malloc(mem_size, SIMD_ALIGN);

    auto get_rand = [](tfloat low, tfloat  high) {
        tfloat t = static_cast<tfloat>(std::rand()) / static_cast<tfloat>(RAND_MAX);
        return (tfloat(1.0) - t) * low + t * high;
    };

    for (i = 0; i < OPT_N; i++)
    {
        CallResultParallel[i] = 0.0;
        CallConfidence[i] = -1.0;
        StockPrice[i] = get_rand(5.0, 50.0);
        OptionStrike[i] = get_rand(10., 25.0);
        OptionYears[i] = get_rand(1.0, 5.0);
        initialVol[i] = get_rand(0.01, 0.015);
        kappa[i] = get_rand(6.0, 7.0);
        theta[i] = get_rand(0.02, 0.03);
        xi[i] = get_rand(0.45, 0.6);
    }

    printf("Pricing Asian Options under Heston Model ... \n");

    MonteCarloTimeStep(
        CallResultParallel,
        CallConfidence,
        StockPrice,
        OptionStrike,
        OptionYears,
        initialVol,
        kappa,
        theta,
        xi);


    my_aligned_free(CallResultParallel);
    my_aligned_free(CallConfidence);
    my_aligned_free(StockPrice);
    my_aligned_free(OptionStrike);
    my_aligned_free(OptionYears);
    my_aligned_free(kappa);
    my_aligned_free(theta);
    my_aligned_free(xi);
    my_aligned_free(initialVol);
    return 0;
}
