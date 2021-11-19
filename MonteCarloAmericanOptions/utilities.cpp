/*

Copyright (c) 2019-2020 
Intel Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#include "utilities.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <mkl_vsl.h>
#include <omp.h>

#define FUZZ 0.00000001

/*static void generateRand_vanilla(REAL* z, int num, REAL mean, REAL stddev)
{
    for (int n = 0; n < num/2; n++)
    {
        // two uniform random numbers 
        REAL x; 
        REAL y;
        REAL r, d;

        do {
            x = REAL(2) * ((REAL)rand() / (REAL)RAND_MAX) - 1;
            y = REAL(2) * ((REAL)rand() / (REAL)RAND_MAX) - 1;
            r = x*x + y*y;

        } while (r == REAL(0) || r >= REAL(1) );

        d = std::sqrt(-REAL(2) * std::log(r) / r);
        z[2*n] = x*d;
        z[2*n+1] = y*d;
    }
}
*/
void perform_sims(REAL **a, VSLStreamStatePtr *streams, int nthreads, REAL fwdFactor, REAL volSqrtdt)
{
    int nb = numPaths/NBUF;
    int ns_blocks =1;

    int outerSteps = numSteps/ns_blocks;

    #pragma omp parallel for
    for (int ip=0; ip <nb; ip++)
    {
        int tid = omp_get_thread_num();
        __attribute__((aligned(64))) REAL rn[NBUF];
        int block_size = (ip != (nb-1)) ? NBUF:(numPaths-(nb-1)*NBUF);
 
        for (int is =1; is <=outerSteps; is++)
        {
            for (int iss =0; iss<ns_blocks; iss++)
            {
                int s = (is-1)*ns_blocks + iss+1;
                REAL *ap = &a[s][0];
                REAL *ap1 = &a[s-1][0];

                int streamAddr = (tid)*numSteps + (s-1);
                //std::cout << "streamAddr " << streamAddr << " s = " << s <<std::endl;

                VSLStreamStatePtr stream = *(streams + streamAddr ); 

                #if REAL_BYTES == 8
                    vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, block_size, rn, 0.0, 1.0);
                #else
                    vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, block_size, rn, 0.0, 1.0);
                #endif

                // Alternative RNG - if using this, remove all of the VSL calls. 
                //generateRand_vanilla(rn,block_size,0.0,1.0);

                #pragma omp simd
                for (int ipp = 0; ipp<block_size; ipp++)
                {
                    int p = NBUF*ip+ ipp;
                    REAL z = rn[ipp];
                    ap[p] = ap1[p] * exp2(fwdFactor+ volSqrtdt*z);
                }
            }
        }
    }
}

inline double det3 (double a, double b, double c,
                    double d, double e, double f,
                    double g, double h, double i)
{ return a*e*i + d*h*c + g*b*f - g*e*c - d*b*i - a*h*f; }

// determinant of a 3x3 matrix (for inversion)
//inline REAL determinant(REAL **x)

inline REAL determinant(REAL x[][3])
{
    REAL det;
    det = det3(x[0][0], x[0][1], x[0][2],
               x[1][0], x[1][1], x[1][2],
               x[2][0], x[2][1], x[2][2]);
    return det;
}

void invert3_arr(REAL a[][3], REAL x[][3], int numRows, int numCols)
{
    REAL det  = determinant(x);

    if (std::abs(det)<FUZZ)
        std::cout << "Trying to invert non invertible matrix" << std::endl;

    a[0][0] = (x[1][1]*x[2][2] - x[2][1]*x[1][2])/det;
    a[0][1] = (x[0][2]*x[2][1] - x[0][1]*x[2][2])/det;
    a[0][2] = (x[0][1]*x[1][2] - x[0][2]*x[1][1])/det;
    a[1][0] = (x[1][2]*x[2][0] - x[1][0]*x[2][2])/det;
    a[1][1] = (x[0][0]*x[2][2] - x[0][2]*x[2][0])/det;
    a[1][2] = (x[0][2]*x[1][0] - x[0][0]*x[1][2])/det;
    a[2][0] = (x[1][0]*x[2][1] - x[1][1]*x[2][0])/det;
    a[2][1] = (x[0][1]*x[2][0] - x[0][0]*x[2][1])/det;
    a[2][2] = (x[0][0]*x[1][1] - x[0][1]*x[1][0])/det;
}
