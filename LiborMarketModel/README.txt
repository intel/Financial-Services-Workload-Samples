FINANCIAL SAMPLE WORKLODAS:
===========================

The FSI customers deploy their applications in a multi-process (one process per core) environment. So, all the workload in this sample will execute in Multi-processes/Instances mode. Each process will execute single thread and inter process communication and resources sharing between the threads wouldn't exist while working on separate instruments (i.e. Options). In the package the compilation will build AVX512 binaries from source code. And Intel C++ Compiler is needed to compile the sources.

BLACK SCHOLES OPTION PRICING:
==============================

Black Scholes is a popular mathematical model used in finance for European option valuation. This is a double precision version. This benchmark prices call and put options
using Black Scholes formula. Large set of options are fed via an input array and results are written into call and put arrays. Exp function is used for simulating Cumulative
e Normal Distribution (CND). This is computing bound, double precision workload and benefitted from Turbo. Runs in HT OFF mode.

BINOMIAL OPTIONS PRICICNG:
==========================

Binomial Options Pricing is a lattice-based approach (Cox, Ross and Rubenstein method) that uses a discrete-time model of the varying price over time of the underlying fina
ncial instrument (a European call option). At every time step, the value of stock “S” can either go up by u*S or go down by “v*S” (0<v<1<u). The values of “u” and “v” are c
onstant for every time step. This is compute bound, double precision workload and benefitted from Turbo and SMT. Runs in HT ON mode.

MONTECARLO EUROPEAN OPTIOSN PRICING:
===================================

Monte Carlo European options is a numerical method that uses statistical sampling techniques to approximate solutions to quantitative problems. In Computational Finance, Monte Carlo algorithms are used to calculate the value of an option with multiple sources of uncertainties and random features, such as changing interest rates, stock prices
or exchange rates, etc. to evaluate complex instruments, portfolios, and investments. This is compute bound, double precision workload and benefitted from Turbo and SMT. Runs in HT ON mode.

MONTECARLO AMERICAN OPTIOSN PRICING:
====================================
This benchmark prices an American put option by Monte Carlo Simulation using Long Staff Schwartz
Least Square Approach (LS). Volatility and Interest rates are constants. Time steps and paths loops were interchanged to achieve maximum performance. This is compute bound,
 double precision workload and benefitted from Turbo and SMT. Runs in HT ON mode.


LIBOR MARKET MODEL:
===================
This benchmark prices a portfolio of LIBOR swaptions on a LIBOR Market Model using a Monte Carlo simulation. It also simultaneously calculates delta Greeks using path-wise
Adjoint Differentiation (AD). In each path, the LIBOR forward rates are generated randomly at all required maturities following the LIBOR Market Model. This is compute bound, double precision workload and benefitted from Turbo and SMT. Runs in HT ON mode.

BASIC INSTRUCTIONS HOW TO USE THE SAMEPLE
==========================================
1. untar the FSI_Sample_Workloads.tar file. There will be 5 workload folders.
2. Each Workload Folder has sources and MPTest directory. The binaries need to be executed in MPTEST Directory.
3. To compile and build Intel C++ Compiler environment needed to be set. After compilation of the sources, the built binary need to be copy in respective MPTest folder. MPTest folder has 3 shell script. clean.sh will clean the folders and delete any results from previous runs. to run the binary, the binary which has been copied in the MPTest folder need to be invoked with runbatch.sh script. with top command completions for the test run needed to be confirmed. once execution of the binary completed getresult.sh need to be invoked with a file name which will collect the run result.

EXAMPLE RUN:
=============

[user@SUT FSI_SAMPEL_WORKLOADS] # ls
BinomialOptions  BlackScholes  LiborMarketModel  MonteCarloAmericanOptions  MonteCarloEuropeanOptions
[user@SUT FSI_SAMPEL_WORKLOADS]# cd MonteCarloEuropeanOptions
[user@SUT FSI_SAMPEL_WORKLOADS MonteCarloEuropeanOptions]# ls
Makefile  MonteCarloInsideBlockingDP.cpp  MPTest  REDME.txt
[user@SUT FSI_SAMPEL_WORKLOADS MonteCarloEuropeanOptions]# ls MPTest
clean.sh  getresults.sh  runbatch.sh
[user@SUT FSI_SAMPEL_WORKLOADS MonteCarloEuropeanOptions]#source /opt/intel/compilers_and_libraries_........../linux/bin/compilervars.sh intel64
[user@SUT FSI_SAMPEL_WORKLOADS MonteCarloEuropeanOptions]#export PR="Number OF Cores in the system reported by lscpu"
//put numerical value of the cores. For BlackScholes HT need to OFF state and and for other WLs require HT ON
[user@SUT FSI_SAMPEL_WORKLOADS MonteCarloEuropeanOptions]#make
[user@SUT FSI_SAMPEL_WORKLOADS MonteCarloEuropeanOptions]#ls
Makefile  MonteCarloInsideBlockingDP.avx512  MonteCarloInsideBlockingDP.cpp  MonteCarloInsideBlockingDP.optrpt  MPTest  REDME.txt
[user@SUT FSI_SAMPEL_WORKLOADS MonteCarloEuropeanOptions]#cp MonteCarloInsideBlockingDP.avx512 MPTest
[user@SUT FSI_SAMPEL_WORKLOADS MonteCarloEuropeanOptions]#cd MPTest
[user@SUT FSI_SAMPEL_WORKLOADS MonteCarloEuropeanOptions]#ls MPTEST
clean.sh  getresults.sh  MonteCarloInsideBlockingDP.avx512  runbatch.sh

[user@SUT FSI_SAMPEL_WORKLOADS MonteCarloEuropeanOptions]#./runbatch.sh MonteCarloInsideBlockingDP.avx512

// Need to make sure with top command that all the precesses started and ended

[user@SUT FSI_SAMPEL_WORKLOADS MonteCarloEuropeanOptions]# top
top - 14:04:25 up 3 days,  1:26,  3 users,  load average: 4.57, 1.01, 0.47
Tasks: 996 total,  97 running, 899 sleeping,   0 stopped,   0 zombie
%Cpu(s): 99.8 us,  0.2 sy,  0.0 ni,  0.0 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
KiB Mem : 19651051+total, 18213664+free,  2937204 used, 11436676 buff/cache
KiB Swap: 65535996 total, 65535996 free,        0 used. 19252812+avail Mem

   PID USER      PR  NI    VIRT    RES    SHR S  %CPU %MEM     TIME+ COMMAND
303188 root      20   0   34948   3048   1892 R  87.5  0.0   0:04.88 MonteCarloInsid
303189 root      20   0   34948   3052   1892 R  87.5  0.0   0:04.88 MonteCarloInsid
303190 root      20   0   34948   3052   1892 R  87.5  0.0   0:04.88 MonteCarloInsid
303191 root      20   0   34948   3052   1892 R  87.5  0.0   0:04.87 MonteCarloInsid
303192 root      20   0   34948   3052   1892 R  87.5  0.0   0:04.88 MonteCarloInsid
303193 root      20   0   34948   3048   1892 R  87.5  0.0   0:04.88 MonteCarloInsid
303194 root      20   0   34948   3048   1892 R  87.5  0.0   0:04.88 MonteCarloInsid
303195 root      20   0   34948   3052   1892 R  87.5  0.0   0:04.88 MonteCarloInsid
303196 root      20   0   34948   3048   1892 R  87.5  0.0   0:04.88 MonteCarloInsid
303197 root      20   0   34948   3048   1892 R  87.5  0.0   0:04.88 MonteCarloInsid
303198 root      20   0   34948   3052   1892 R  87.5  0.0   0:04.88 MonteCarloInsid

//getting results: invoke getresult.sh script passing with result file name. All data will be dumped in the file and main result will be displayed (execution time in sec/ms
)

[user@SUT FSI_SAMPEL_WORKLOADS MonteCarloEuropeanOptions]#./getresult result
[root@wsl-skl2s-02 MPTest]# ./getresults.sh result
7.00145
[user@SUT FSI_SAMPEL_WORKLOADS MonteCarloEuropeanOptions]#ls
clean.sh       inst12  inst18  inst23  inst29  inst34  inst4   inst45  inst50  inst56  inst61  inst67  inst72  inst78  inst83  inst89  inst94
getresults.sh  inst13  inst19  inst24  inst3   inst35  inst40  inst46  inst51  inst57  inst62  inst68  inst73  inst79  inst84  inst9   inst95
inst0          inst14  inst2   inst25  inst30  inst36  inst41  inst47  inst52  inst58  inst63  inst69  inst74  inst8   inst85  inst90  MonteCarloInsideBlockingDP.avx512
inst1          inst15  inst20  inst26  inst31  inst37  inst42  inst48  inst53  inst59  inst64  inst7   inst75  inst80  inst86  inst91  result.txt
inst10         inst16  inst21  inst27  inst32  inst38  inst43  inst49  inst54  inst6   inst65  inst70  inst76  inst81  inst87  inst92  runbatch.sh
inst11         inst17  inst22  inst28  inst33  inst39  inst44  inst5   inst55  inst60  inst66  inst71  inst77  inst82  inst88  inst93
[user@SUT FSI_SAMPEL_WORKLOADS MonteCarloEuropeanOptions]#

Claculations of the performances: number of options*Number of processes/ Average Time( reported on the screen via running getresults.sh. Unit: Options/sec

MonteCarloEuropeanOptions Throughput =Total Number of Options/Average Time needs to complete all processes 
= (Number of Options * NLoops(Iteration)*Nprocesses) / Time"	getresult calculates average time to complete all processes to complete

MonteCarloAmericanOptions Throughput =Total Number of Paths/Total Time 
= (Number of Path * NLoops*Nprocesses) /(Elapesd Time *NLoops) 
= Number of Path * Nprocesses/Elapsed Time Kpaths/Sec"	getresult.sh calculates average elapsed time which is total time / NLoops

LiborMarketModel Throughput =Total Number of Paths/Total Time 
= (Number of Path * NLoops*NProcesses) /Total Time"	getresult.sh calculates average Total time which is total time / NLoops
	
		
BinomialOptions Throughput =Total Number of Options/Average Time needs to complete all processes 
= (Number of Options * NLoops(Iteration)*NProcesses) / Time"	getresults.sh calculates avrg time tto complete all processes to complete
		
BlackScholes Throughput =Total Number of Options/Average Time needs to complete all processes 
= (Number of Options * Nloops(Iteration)*Nprocesses) / Time"	getresult.sh calculates avrg time tto complete all processes to complete

