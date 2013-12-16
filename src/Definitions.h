#ifndef defs_H_INCLUDED
#define defs_H_INCLUDED

#ifdef __linux__
  #define LINUX
#else
  #define WINDOWS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>       /* time */
#include <cstring>
#include <math.h>
#include <omp.h>



#ifdef WINDOWS
    #include "lapacke.h"
    #include <windows.h>
    #include <cblas.h>
    #define cpu_freq 4.6
#else
    #include "mkl.h"
    //#include "mpi.h"
    //#define cpu_freq 3.0
    #define cpu_freq 3.0
    //#include "include/lapacke.h"

#endif

#ifdef __INTEL_MKL__
    #define blas_set_num_threads(n) mkl_set_num_threads(n)
#else
    #define blas_set_num_threads(n) goto_set_num_threads(n)
#endif

#include <unistd.h>
#include <pthread.h>
#include <limits.h>
#include <queue>
#include <iostream>


//!SETTINGS
#define EXTENDEDTEST 0
#define PRINT 0 //&& 1EXTENDEDTEST
#define ForceCheck 0 || EXTENDEDTEST
#define OUTPUT 0

#define STORAGE_TYPE LAPACK_COL_MAJOR
#define type_precision float

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define _1MB 1024*1024
#define _10MB 10*_1MB
#define _1GB 1024*1024*1024

//!for CPU speed!

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
typedef unsigned __int64 usCount;
static usCount GetUsCount()
{
    static LARGE_INTEGER ticksPerSec;
    static double scalefactor;
    LARGE_INTEGER val;
    if(!scalefactor)
    {
        if(QueryPerformanceFrequency(&ticksPerSec))
            scalefactor=ticksPerSec.QuadPart/1000000000000.0;
        else
            scalefactor=1;
    }
    if(!QueryPerformanceCounter(&val))
        return (usCount) GetTickCount() * 1000000000;
    return (usCount) (val.QuadPart/scalefactor);
}
#else
#include <sys/time.h>
#include <time.h>
#include <sched.h>
typedef unsigned long long usCount;
static usCount GetUsCount()
{
#ifdef CLOCK_MONOTONIC
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ((usCount) ts.tv_sec*1000000000000LL)+ts.tv_nsec*1000LL;
#else
    struct timeval tv;
    gettimeofday(&tv, 0);
    return ((usCount) tv.tv_sec*1000000000000LL)+tv.tv_usec*1000000LL;
#endif
}
#endif
static usCount usCountOverhead, CPUClockSpeed;
#ifdef __GNUC__
#include "x86intrin.h"
#define __rdtsc() __builtin_ia32_rdtsc()
#endif
static usCount GetClockSpeed()
{
  int n;
  usCount start, end, start_tsc, end_tsc;
  if(!usCountOverhead)
  {
    usCount foo=0;
    start=GetUsCount();
    for(n=0; n<1000000; n++)
    {
      foo+=GetUsCount();
    }
    end=GetUsCount();
    usCountOverhead=(end-start)/n;
  }
  start=GetUsCount();
  start_tsc=__rdtsc();
  for(n=0; n<1000; n++)
#ifdef WIN32
    Sleep(0);
#else
    sched_yield();
#endif
  end_tsc=__rdtsc();
  end=GetUsCount();
  return (usCount)((1000000000000.0*(end_tsc-start_tsc))/(end-start-usCountOverhead));
}


using namespace std;

struct Settings
{

    int n;
    int t;
    int m;
    int l;
    int r;
    int p;
    int tb;
    int mb;
    int id;

    int threads;

    bool use_fake_files;

    string fnameAL;
    string fnameAR;
    string fnameY;
    string fnameOutB;

    bool doublefileType;

};

struct Outputs
{

    double duration;
    double gflops;

    double acc_RTL_QLY;
    double acc_loady;
    double acc_loadxr;
    double acc_pre;
    double acc_gemm;
    double acc_b;
    double firstloop;


};


#endif // UTILITY_H_INCLUDED


