/*
Copyright (c) 2019-2020 
Intel Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/


#include <math.h>
#include <stdio.h>
#include <typeinfo>
#if _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <chrono>

#if __INTEL_COMPILER || _MSC_VER
#define FORCEINLINE __forceinline
#define ALIGN(x) __declspec(align(x))
#else
#define FORCEINLINE 
#define ALIGN(x)
#endif

#include "binomial.h"

template<> void binomial<TPRECISION>::initialize() {}

template<typename DataType>
void binomial<DataType>::run()
{
#if _OPENMP
	kmp_set_defaults("KMP_AFFINITY=scatter,granularity=thread");
#endif

#pragma omp parallel
	{
	}

	printf("%s Precision Binomial Option Pricing, version %d.%d\n", isDP() ? "Double" : "Single", MAJOR, MINOR);
	printf("Compiler Version  = %d\n", __INTEL_COMPILER / 100);
	printf("Release Update    = %d\n", __INTEL_COMPILER_UPDATE);
	printf("Build Time        = %s %s\n", __DATE__, __TIME__);
	printf("Input Dataset     = %d\n", optN);

	auto body = [this](
		binomial_vector& h_CallResult,
		const binomial_vector&  S,
		const binomial_vector&  X,
		const binomial_vector&  T,
		const DataType Riskfree,
		const DataType Volatility,
		const int optN,
		const int timestepN
		)
	{
		printf("\nUsing %d OpenMP thread(s)%s.\n", omp_get_max_threads(), BINOMIAL_USE_OMP_SIMD ? " and SIMD Vectorization" : "");
		printf("Pricing %d Option with time step of %d.\n", optN, timestepN);
		auto sTime = std::chrono::system_clock::now();
		auto NUM_STEPS = timestepN;
#pragma omp parallel for
		for (int opt = 0; opt < optN; opt++)
		{
			ALIGN(1024) DataType Call[NUM_STEPS + 1];
			const DataType       Sx = S[opt];
			const DataType       Xx = X[opt];
			const DataType       Tx = T[opt];
			const DataType      dt = Tx / static_cast<DataType>(NUM_STEPS);
			const DataType     vDt = Volatility * std::sqrt(dt);
			const DataType     rDt = Riskfree * dt;
			const DataType      If = std::exp(rDt);
			const DataType      Df = std::exp(-rDt);
			const DataType       u = std::exp(vDt);
			const DataType       d = std::exp(-vDt);
			const DataType      pu = (If - d) / (u - d);
			const DataType      pd = 1.0f - pu;
			const DataType  puByDf = pu * Df;
			const DataType  pdByDf = pd * Df;
			const DataType mul_c = vDt* static_cast<DataType>(2.0);
			DataType id = vDt*static_cast<DataType>(-NUM_STEPS);
#if BINOMIAL_USE_OMP_SIMD
#pragma omp simd
#endif
#pragma unroll(8)
			for (int i = 0; i <= NUM_STEPS; i++)
			{
				DataType d = Sx * std::exp(id) - Xx;
				Call[i] = (d > 0) ? d : 0;
				id += mul_c;
			}
			// Start at the final tree time step nodes(leaves) and walk backwards
			// to calculate the call option price. 
			for (int i = NUM_STEPS; i > 0; i--)
#if BINOMIAL_USE_OMP_SIMD
#pragma omp simd
#endif
#pragma unroll(8)
				for (int j = 0; j <= i - 1; j++)
					Call[j] = puByDf * Call[j + 1] + pdByDf * Call[j];
				h_CallResult[opt] = Call[0];
		}
		auto eTime = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = eTime - sTime;
		printf("Completed in %10.5f seconds. Options per second: %10.5f\n", elapsed_seconds.count(), static_cast<double>(optN) / (elapsed_seconds.count()));
		printf("Time Elapsed =  %10.5f seconds\n", elapsed_seconds.count());
	};

	body(h_CallResult, h_StockPrice, h_OptionStrike, h_OptionYears, Riskfree, Volatility, optN, stepsN);
}

template void binomial<TPRECISION>::initialize();
template void binomial<TPRECISION>::run();

