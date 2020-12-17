#Copyright (c) Intel Corporation 2019-2020

#Intel Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#!/bin/bash

CXX=icpc

CXXFLAGS=-O3 -g -restrict -fimf-domain-exclusion=31 -fimf-accuracy-bits=11 -no-prec-div -no-prec-sqrt -fno-alias -qopt-report=5 -qopt-report-phase=vec


LDFLAGS=-lpthread -ltbbmalloc

all: binomial.avx512

binomial.avx512: binomial_cpu.cpp binomial_main.cpp binomial.h
	$(CXX) -std=c++14 -tbb -xCORE-AVX512 -qopenmp -qopt-report-file=binomial.avx512.optrpt $(CXXFLAGS) -o $@ $(LDFLAGS) binomial_main.cpp binomial_cpu.cpp -qopt-zmm-usage=high


clean	:
	rm -f binomial.avx512 *.optrpt
