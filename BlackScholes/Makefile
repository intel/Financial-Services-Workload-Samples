#Copyright (c) Intel Corporation 2019-2020

#Intel Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#!/bin/bash

CXX=icpc
COMPILER_VERSION:="$(CXX)-$(shell $(CXX) --version | head -n1 | cut -d' ' -f4)"

CXXFLAGS=-g -O3 -qopenmp  -fimf-precision=low -restrict -qopt-report=5  -fimf-domain-exclusion=31 -no-prec-div  -no-prec-sqrt -fno-alias \
         -DCOMPILER_VERSION=\"$(COMPILER_VERSION)\" -qopt-report-phase=loop -qopt-report-phase=vec -vecabi=cmdtarget 

LDFLAGS= -ltbbmalloc
all: BlackScholesDP.avx512 

BlackScholesDP.avx512: BlackScholesDP.cpp
	$(CXX) -xCORE-AVX512 -qopt-zmm-usage=high $(CXXFLAGS) $(LDFLAGS) -qopt-report-file=DP_AVX512.optrpt  -o $@ $^

clean   :
	rm -f *optrpt BlackScholesDP.avx512   
