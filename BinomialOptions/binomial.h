/*
Copyright (c) 2019-2020 
Intel Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#ifndef __BINOMIAL_H__ 
#define __BINOMIAL_H__ 

#ifndef TPRECISION 
#define TPRECISION double
#endif

#ifndef VERBOSE 
#define VERBOSE    false
#endif

/******* VERSION *******/
#define MAJOR 1
#define MINOR 1

#ifndef ENABLE_CPU 
#define ENABLE_CPU 1
#endif

#ifndef ENABLE_OPENCL 
#define ENABLE_OPENCL    0
#endif

#ifndef ENABLE_SYCL 
#define ENABLE_SYCL 0
#endif

const float VOLATILITY = 0.10f;
const float RISKFREE = 0.06f;
const int SIMDALIGN = 32;

const int   NUM_STEPS = 2048;
const int   OPT_N = 48000;

#ifndef __INTEL_COMPILER 
#define __INTEL_COMPILER    0
#endif

#ifndef __INTEL_COMPILER_UPDATE
#define __INTEL_COMPILER_UPDATE    0
#endif

#ifndef BINOMIAL_USE_TBB_ALLOCATOR
#define BINOMIAL_USE_TBB_ALLOCATOR    1
#endif

#ifndef BINOMIAL_USE_OMP_SIMD
#define BINOMIAL_USE_OMP_SIMD    1
#endif

#ifndef BINOMIAL_USE_ACC_THREADS
#define BINOMIAL_USE_ACC_THREADS    0
#endif

#include <vector>
#include <cstdlib>
#if BINOMIAL_USE_TBB_ALLOCATOR
#include <tbb/cache_aligned_allocator.h>
#endif

template<typename DataType = double>
class binomial
{
public:
	binomial(bool, int, int, DataType, DataType);
	~binomial();
	void run();
	void check(const bool);
	inline DataType rand(DataType low, DataType high) {
		DataType t = static_cast<DataType>(std::rand()) / static_cast<DataType>(RAND_MAX);
		return (DataType(1.0) - t) * low + t * high;
	}

	inline bool isDP()
	{
		return (int(sizeof(DataType)) > 4);
	}

	typedef
#if BINOMIAL_USE_TBB_ALLOCATOR
		std::vector<DataType, tbb::cache_aligned_allocator<DataType> >
#else
		std::vector<DataType>
#endif
		binomial_vector;

	binomial_vector h_CallResult;
	binomial_vector h_StockPrice;
	binomial_vector h_OptionStrike;
	binomial_vector h_OptionYears;

	bool verbose;
	int optN;
	int stepsN;
	DataType Riskfree;
	DataType Volatility;

private:
	binomial();
	void initialize();
};

#endif
