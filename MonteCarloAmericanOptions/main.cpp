/*
Copyright (c) 2019-2020 
Intel Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

This application prices American Options using Longstaff Schwartz Method
Ref: https://people.math.ethz.ch/~hjfurrer/teaching/LongstaffSchwartzAmericanOptionsLeastSquareMonteCarlo.pdf
This file also contains problem-size definitions

*/

#include "utilities.h"

#include "chrono"
#include <ctime>

#include "Timer.hpp"
#include <cmath>
#include <cstdlib>
#include <malloc.h>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <mkl_vsl.h>
#include <omp.h>

#define progLoops 50//*120-power

inline int GetStripInterval( int& pfirst, int& plast, int numPathStrips ) 
{
   // Figure out which strips current thread computes.
    int nthread = omp_get_num_threads();
    int chunksize = (numPathStrips+nthread-1)/nthread;
    int id = omp_get_thread_num();
    pfirst = id*chunksize;
    plast = std::min(pfirst+chunksize,numPathStrips);
    return id;
}


int main(int argc, char** argv)
{    
    std::cout << "numPaths = " << numPaths << ", numSteps = " << numSteps << std::endl;
        
    /* ------------- option parameters      */
    REAL spot = 100.0;
    REAL strike = 110.0;
    REAL vol = 0.25;
    REAL maturity = 2.0;
    
    Timer tt;
    tt.start();
    
    // pre-compute some often-used values
    REAL dt = maturity / REAL(numSteps);
    REAL volSqrtdt = vol*std::sqrt(dt);
    REAL fwdFactor = vol*vol*dt*-0.5;
    double vsqrtdt_log2e = volSqrtdt * M_LOG2E;
    double fwdFactor_log2e = fwdFactor * M_LOG2E;
    REAL avgPathFactor = REAL(1.0) / REAL(numPaths);
        
    REAL *cashFlowPut1;
    cashFlowPut1 = (REAL*) memalign(64,sizeof(REAL)*(numPaths));
    
    #pragma omp parallel for 
    for(int i=0; i<numPaths; i++)
        cashFlowPut1[i] = 0.0;
  
    int nStepsPlusOne = numSteps+1;
    // now go back in time using regression
  
    REAL **asset;
    asset = (REAL**) memalign(64,sizeof(REAL*)*(nStepsPlusOne));
    for(int i =0; i<nStepsPlusOne; i++)
        asset[i] = (REAL *) memalign(64,sizeof(REAL)*(numPaths));
    
    int nthreads = 1;
    nthreads = omp_get_max_threads();
  
    // RK
  
    for (int i =0; i<nStepsPlusOne; i++)
    {
        #pragma omp parallel for //simd 
        for(int j =0; j<numPaths; j++)
            asset[i][j] = spot;
    }
    
    VSLStreamStatePtr *streams;
    int nb = numPaths/NBUF;
    
    int streamCount = nthreads*numSteps;

    streams = (VSLStreamStatePtr *) malloc(sizeof(VSLStreamStatePtr)*streamCount);

    /* Initialize RNG */
    #pragma omp parallel for
    for (int i =0; i < nthreads; i++)
    {
        for( int ij=0 ; ij<numSteps ; ij++)
        {
            int j = numSteps*i + ij;
            vslNewStream( streams + j, VSL_BRNG_MRG32K3A , 12345 );
            vslSkipAheadStream(streams[j], nb*j);  // offset for current thread
        }
    }

double initAssetTime=0, computeTime=0, elapsed=0.0, reginitTime = 0.0, applyTime = 0.0; 
	REAL amerpayoffput = 0.0, finalputpayoff=0.0;
	tt.stop();
	std::cout << "Number of Loops = " << progLoops << std::endl;
    std::cout << "nThreads  " << nthreads << std::endl;  
    std::cout << "Init time: " << tt.getTime() << " ms" << std::endl;  
    
    for (int pLoop=0; pLoop<progLoops+1; pLoop++)
    {
        tt.start();
        Timer ttt ;
    
        ttt.start();
        perform_sims(asset,streams,nthreads,fwdFactor_log2e,vsqrtdt_log2e);
        ttt.stop();

        if (pLoop >=1)
            initAssetTime += ttt.getTime();
        //std::cout << "sims done" << std::endl;
        
        REAL betaPut0=0.0, betaPut1=0.0, betaPut2=0.0;
        int numPathStrips = numPaths/VECLEN;
  
        for (int s = numSteps-1; s > 0; s--)
        {
            // these sums compute the 3x3 regression matrix elements on the fly,
            // i.e., sum(1), sum(assetval), sum(assetval**2), sum(assetval**3)
            // for all paths where assetval <= strike
            REAL sum1=0.0, sum2=0.0, sum3=0.0, sum4=0.0;
            REAL sum0i = 0;  // use integer for first sum
            REAL x0=0, x1=0, x2=0;
            Timer t1, t2;
            REAL *ap = &asset[s][0];
            
            t2.start();
            
            #pragma omp parallel reduction(+:sum0i,sum1,sum2,sum3,sum4,x0,x1,x2)
            // this block of code is a thread body
            {
                int pfirst, plast;
                GetStripInterval( pfirst, plast, numPathStrips );
                __attribute__((aligned(64))) REAL sum0i_local[VECLEN];
                __attribute__((aligned(64))) REAL sum1_local[VECLEN];
                __attribute__((aligned(64))) REAL sum2_local[VECLEN];
                __attribute__((aligned(64))) REAL sum3_local[VECLEN];
                __attribute__((aligned(64))) REAL sum4_local[VECLEN];
                __attribute__((aligned(64))) REAL x0_local[VECLEN];
                __attribute__((aligned(64))) REAL x1_local[VECLEN];
                __attribute__((aligned(64))) REAL x2_local[VECLEN];
                __attribute__((aligned(64))) REAL assetval[VECLEN];
        
                for (int in = 0; in < VECLEN; in++) 
                {
                    sum0i_local[in] = 0;
                    sum1_local[in] = 0.0;
                    sum2_local[in] = 0.0;
                    sum3_local[in] = 0.0;
                    sum4_local[in] = 0.0;
                    x0_local[in] = x1_local[in] = x2_local[in] = 0.0;
                }       
            
                for( int ip=pfirst; ip<plast; ++ip ) 
                {       
                    int p = VECLEN*ip;
            
                    for (int i=0; i< VECLEN; i++)
                        assetval[i] = ap[p+i];
        
                    if (s < numSteps-1)
                    {
                        for(int in=0;in<VECLEN;in++)
                        {
                            if (assetval[in] <= strike)
                            {
                                REAL regImpliedCashFlow;
                                REAL payoff;

                                regImpliedCashFlow = betaPut0 + betaPut1*assetval[in] + betaPut2*(assetval[in]*assetval[in]);
                                payoff = strike-assetval[in];
                                if( payoff>regImpliedCashFlow )
                                {
                                    for (int i=0; i<VECLEN; i++)
                                        cashFlowPut1[p+in] = payoff;
                                }
                            }
                        }
                    }
                    else
                    {
                        for (int i=0; i<VECLEN; i++)
                            cashFlowPut1[p+i] = std::max(strike - assetval[i], REAL(0)) * avgPathFactor;
                    }
        
                    if( s > 0 ) 
                    {
                        for(int in=0; in<VECLEN; in++)
                        {
                            if (assetval[in] <= strike) 
                            {
                                REAL assetval2;
                                assetval2 = assetval[in]*assetval[in];
                                sum0i_local[in] += 1;
                                sum1_local[in] += assetval[in];
                                sum2_local[in] += assetval2;
                                sum3_local[in] += assetval2 * assetval[in];
                                sum4_local[in] += assetval2 * assetval2;
        
                                x0_local[in] += cashFlowPut1[p+in];
                                x1_local[in] += cashFlowPut1[p+in]*assetval[in];
                                x2_local[in] += cashFlowPut1[p+in]*assetval2;
                            }
                        }
                    }
                } //End of ip loop
            
                for(int in=0; in<VECLEN; in++)
                {    
                    sum0i +=sum0i_local[in];
                    sum1 += sum1_local[in];
                    sum2 += sum2_local[in];
                    sum3 += sum3_local[in];
                    sum4 += sum4_local[in];
                    x0 += x0_local[in];
                    x1 += x1_local[in];
                    x2 += x2_local[in];
                }
            } //End of parallel region
            
            t2.stop();
            if (pLoop >=1)
                reginitTime += t2.getTime();
            
            if( s==0 ) 
                break;
    
            //compute the regression
            REAL ridgecorrect = 0.01;
            REAL sum0 = REAL(sum0i);
            sum0+=ridgecorrect; sum1+=ridgecorrect; sum2+=ridgecorrect;sum3+=ridgecorrect;sum4+=ridgecorrect;
            REAL sqrtest1[3][3];
            sqrtest1[0][0] = sum0; sqrtest1[0][1] = sum1; sqrtest1[0][2] = sum2;
            sqrtest1[1][0] = sum1; sqrtest1[1][1] = sum2; sqrtest1[1][2] = sum3;
            sqrtest1[2][0] = sum2; sqrtest1[2][1] = sum3; sqrtest1[2][2] = sum4;
            REAL invsqr1[3][3];
            invert3_arr(invsqr1, sqrtest1, 3,3);
            //std::cout << "invert done " << std::endl;
        
            //invsqr.transpose(); -> it's symmetric
            // the inverse of a symmetric matrix is a symmetric matrix, so we can use
            // symmetric matrix multiplication here
            t2.start();
            betaPut0 = x0*invsqr1[0][0] + x1*invsqr1[1][0] + x2*invsqr1[2][0];
            betaPut1 = x0*invsqr1[0][1] + x1*invsqr1[1][1] + x2*invsqr1[2][1];
            betaPut2 = x0*invsqr1[0][2] + x1*invsqr1[1][2] + x2*invsqr1[2][2];
            
            t2.stop();
            if (pLoop >=1) 
            applyTime += t2.getTime();
        }  // time steps loop

        // final results
        amerpayoffput = 0;
        finalputpayoff = 0.0;
        #pragma omp parallel for reduction(+:amerpayoffput,finalputpayoff)
        for(int i=0; i < numPaths; i++)
        {
            amerpayoffput += cashFlowPut1[i];
            finalputpayoff += std::max(strike - asset[numSteps][i], 0.0);
        }
        
        amerpayoffput /= numPaths;
        finalputpayoff /= numPaths;
        
        tt.stop();
        if (pLoop >=1) 
            computeTime += tt.getTime();

    } //End of progLoop
    
    elapsed = computeTime/progLoops;
    std::cout << "init asset time = " << initAssetTime/progLoops << std::endl;
    std::cout << "reginit Time = " << reginitTime/progLoops << std::endl;
    std::cout << "apply Time   = " << applyTime/progLoops << std::endl;
std::cout << "Computation time per loop: " << elapsed << " ms" << std::endl;
    std::cout << "Paths/sec: " << numPaths/elapsed << " Kpaths/sec" << std::endl;
        
    std::cout <<std::setprecision(9) << "European Put Option Price: "<< finalputpayoff << std::endl;
    std::cout <<std::setprecision(9) << "American Put Option Price: "<< amerpayoffput << std::endl;

}



