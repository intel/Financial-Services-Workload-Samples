/*
Copyright (c) 2019-2020 
Intel Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#include "binomial.h"

#include<cmath>
#include<limits>
#include<cstdio>
#include<algorithm>


// Black-Scholes Reference
void BlackScholesRefImpl(
    double& callResult,
    double Sf, //Stock price
    double Xf, //Option strike
    double Tf, //Option years
    double Rf, //Riskless rate
    double Vf  //Volatility rate
)
{
    // BSM Formula: https://www.nobelprize.org/prizes/economic-sciences/1997/press-release/
	double S = Sf, L = Xf, t = Tf, r = Rf, sigma = Vf;
	double N_d1 = 1. / 2. + 1. / 2. * erf(((log(S / L) + (r + 0.5 * sigma * sigma) * t) / (sigma * sqrt(t))) / sqrt(2.));
	double N_d2 = 1. / 2. + 1. / 2. * erf(((log(S / L) + (r - 0.5 * sigma * sigma) * t) / (sigma * sqrt(t))) / sqrt(2.));
	callResult = (S * N_d1 - L*exp(-r*t) * N_d2);
}

template <class T>
void BlackScholesCPUTemplate(double *h_CallResult, T *h_StockPrice, T *h_OptionStrike, T *h_OptionYears, T Riskfree, T Volatility, int optN)
{
	for (int opt = 0; opt < optN; opt++)
		BlackScholesBodyCPU(h_CallResult[opt], h_StockPrice[opt], h_OptionStrike[opt], h_OptionYears[opt], Riskfree, Volatility);
}

template <class Basetype>
void BinomialPragmaSIMD(
	Basetype *h_CallResult,
	Basetype *S,
	Basetype *X,
	Basetype *T,
	Basetype Riskfree,
	Basetype Volatility,
	int optN,
	int timestepN
);

template<> binomial<>::~binomial() {}

int main(int argc, char** argv)
{
	binomial<TPRECISION> test(VERBOSE, OPT_N, NUM_STEPS, RISKFREE, VOLATILITY);
	test.run();
	test.check(test.verbose);
	return 0;
}

template<typename DataType>
void binomial<DataType>::check(const bool verbose)
{
	if (verbose)
	{
		std::printf("Creating the reference result...\n");
		//Calculate options values on CPU
		std::vector<double> h_CallResultCPU(optN);

		for (int opt = 0; opt < optN; opt++)
			BlackScholesRefImpl(h_CallResultCPU[opt], h_StockPrice[opt], h_OptionStrike[opt], h_OptionYears[opt], Riskfree, Volatility);

		std::printf("Benchmarking the reference result...\n");

		double sum_delta = 0,
			sum_ref = 0,
			max_delta = 0.0,
			sumReserve = 0,
			errorVal;

		for (auto i = 0; i < optN; i++) {
			auto ref = h_CallResultCPU[i];
			auto delta = fabs(h_CallResultCPU[i] - h_CallResult[i]);
			if (delta > max_delta)
			{
				max_delta = delta;
			}
			sum_delta += delta;
			sum_ref += fabs(ref);
		}
		sumReserve /= (double)optN;
		if (sum_ref > 1E-5)
			printf("L1 norm: %E\n", errorVal = sum_delta / sum_ref);
		else
			printf("Avg. diff: %E\n", errorVal = sum_delta / (double)optN);
		printf((errorVal < 5e-4) ? "TEST PASSED\n" : "TEST FAILED\n");

	}
}

template<typename DataType>
binomial<DataType>::binomial(bool _verbose, int _optN, int _stepsN, DataType _Riskfree, DataType _Volatility)
	:verbose(_verbose), optN(_optN), stepsN(_stepsN), Riskfree(_Riskfree), Volatility(_Volatility)
{
	h_CallResult.resize(optN);
	h_StockPrice.resize(optN);
	h_OptionStrike.resize(optN);
	h_OptionYears.resize(optN);

	std::fill(h_CallResult.begin(), h_CallResult.end(), 0.0);
	std::generate(h_StockPrice.begin(), h_StockPrice.end(), [this]()->DataType { return this->rand(static_cast<DataType>(5.0), static_cast<DataType>(50.0)); });
	std::generate(h_OptionStrike.begin(), h_OptionStrike.end(), [this]()->DataType { return this->rand(static_cast<DataType>(10.0), static_cast<DataType>(25.0)); });
	std::generate(h_OptionYears.begin(), h_OptionYears.end(), [this]()->DataType { return this->rand(static_cast<DataType>(1.0), static_cast<DataType>(5.0)); });
	initialize();
}

template binomial<TPRECISION>::binomial(bool, int, int, TPRECISION, TPRECISION);
template void binomial<TPRECISION>::check(const bool);

