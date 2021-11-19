COMMON_FLAGS_AVX512=" -std=c++14 -fimf-precision=low -no-prec-div -no-prec-sqrt -vecabi=cmdtarget -qopt-report=5 -qopt-report-phase=vec  -Wall -mkl -xCORE-AVX512 -qopt-zmm-usage=high -fimf-domain-exclusion=31 -fno-alias -qopt-assume-safe-padding -restrict -Bstatic -Wl,--start-group -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -Wl,--end-group -Bdynamic -lpthread "

icc -O3 -g -o asian-opt-avx512.bm_rng.double.new  $COMMON_FLAGS_AVX512 MCAsianOptions.cpp &


COMMON_FLAGS_AVX2=" -std=c++14  -fimf-precision=low -no-prec-div -no-prec-sqrt -vecabi=cmdtarget -qopt-report=5 -qopt-report-phase=vec  -Wall -mkl -mavx2 -fimf-domain-exclusion=31 -fno-alias -qopt-assume-safe-padding -restrict -Bstatic -Wl,--start-group -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -Wl,--end-group -Bdynamic -lpthread  "


icc -O3 -g -o asian-opt-avx2.bm_rng.double.new -DRNGBM=1 $COMMON_FLAGS_AVX2 MCAsianOptions.cpp &


COMMON_FLAGS_NoVEC=" -std=c++14 -no-vec -no-simd -qno-openmp-simd -fimf-precision=low -no-prec-div -no-prec-sqrt -vecabi=cmdtarget -qopt-report=5 -qopt-report-phase=vec  -Wall -mkl -fimf-domain-exclusion=31 -fno-alias -qopt-assume-safe-padding -restrict -Bstatic -Wl,--start-group -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -Wl,--end-group -Bdynamic -lpthread  "


icc -O3 -g -o asian-opt-NoVEC.bm_rng.double.new -DRNGBM=1 $COMMON_FLAGS_NoVEC MCAsianOptions.cpp &



wait

