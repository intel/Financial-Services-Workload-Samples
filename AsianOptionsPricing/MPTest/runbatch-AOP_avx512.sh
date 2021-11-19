#Copyright (c) 2019-2020 
#Intel Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#!/bin/bash


rm -rf inst*

#export PR=`lscpu | grep '^CPU(s)' | cut -d':' -f 2 | sed -e 's/^[ \t]*//'`

for i in `seq 0 $(($PR - 1))`
do

   mkdir inst$i
done

for i in `seq 0 $(($PR - 1))`
do
   cp asian-opt-avx512.bm_rng.double.new inst$i/. 
done

export OMP_NUM_THREADS=1
export MKL_ENABLE_INSTRUCTIONS=AVX512
#export MKL_ENABLE_INSTRUCTIONS=SSE4_2



#EXEFILE=$1

for i in `seq 0 $(($PR - 1))`
do
   taskset -c $i inst$i/asian-opt-avx512.bm_rng.double.new > inst$i/results$i.txt 2>&1 &
done
