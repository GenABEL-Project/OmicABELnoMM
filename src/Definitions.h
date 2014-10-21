#ifndef defs_H_INCLUDED
#define defs_H_INCLUDED



#include <limits>
#include <list>


#ifdef __linux__
    #define LINUX
#else
    #define WINDOWS
#endif

#include <stdint.h>
#include <unistd.h>
#include <limits.h>
#include <queue>
#include <map>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>       /* time */
#include <cstring>
#include <math.h>
#include <omp.h>
#include <pthread.h>
#include <utility>
#include <iomanip>
#include <mpi.h>

#ifdef WINDOWS
    #include <windows.h>
#else

#endif

//!For intel use propetary MKL, it will be preferred over others
#ifdef __INTEL_MKL__
    #include <boost/math/distributions/normal.hpp>
    #pragma message("MKL will Probably NOT compile")
    #include "mkl.h"
    #include "cblas.h"
    #include <lapacke.h>
    #define blas_set_num_threads(n) mkl_set_num_threads(n)
    #define STORAGE_TYPE LAPACK_COL_MAJOR
#else

    //!For AMD systems use the proper ACML library, preferred over openblas ON AMD
    #ifdef _acml_
        #pragma message("Compiled with AMD ACML")
        #define blas_set_num_threads(n) omp_set_num_threads(n)

        #include <acml.h>

        #define lapack_int int

        #define CblasTrans 'T'
        #define CblasNoTrans 'N'
        #define CblasUpper 'U'
        #define CblasColMajor 1


        #define STORAGE_TYPE CblasColMajor

         #define cblas_sdot sdot
        #define cblas_snrm2 snrm2
        #define cblas_saxpy saxpy

        #ifndef BLASdefs_H_INCLUDED
        #define BLASdefs_H_INCLUDED

        inline  void cblas_sgemm(int storage, char transa, char transb, int m, int n, int k, float alpha, float *a, int lda, float *b, int ldb, float beta, float *c, int ldc)
        {
           sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        }

        inline  lapack_int LAPACKE_sposv( int matrix_order, char uplo, lapack_int n,  lapack_int nrhs, float* a, lapack_int lda, float* b, lapack_int ldb )
        {
            int info;
            sposv( uplo, n,  nrhs, a, lda, b,  ldb, &info);
            return info;
        }


        inline  void cblas_ssyrk(int Order, char uplo, char Trans,
		 int N, int K, float alpha, float *A, int lda,  float beta, float *C, int ldc)
        {
            ssyrk(uplo, Trans, N, K, alpha, A, lda, beta, C, ldc);
        }

        inline  lapack_int LAPACKE_sgels( int matrix_order, char trans, lapack_int m, lapack_int n, lapack_int nrhs, float* a,  lapack_int lda, float* b, lapack_int ldb )
        {
            int info;
            sgels(trans, m, n, nrhs, a, lda, b, ldb,&info);
            return info;
        }

        inline  lapack_int LAPACKE_spotrf( int matrix_order, char trans, lapack_int n, float* a,  lapack_int lda )
        {
            int info;
            sgels(trans, n, nrhs, a, lda,&info);
            return info;
        }

        inline  lapack_int LAPACKE_spotri( int matrix_order, char trans,lapack_int n, float* a,  lapack_int lda)
        {
            int info;
            sgels(trans, n, a, lda,&info);
            return info;
        }

        #endif



        #include <boost/math/distributions/normal.hpp>

    #else

        //!IF MKL is not present on INTEL, use openblas
        #ifdef _openblas_
            #include <boost/math/distributions/normal.hpp>
            #pragma message("Compiled with OPENBLAS")
            #define STORAGE_TYPE LAPACK_COL_MAJOR
            #include "cblas.h"
            #include <lapacke.h>
            extern "C" void openblas_set_num_threads(int num_threads);
            #define blas_set_num_threads(n) openblas_set_num_threads(n)
        #endif

    #endif

#endif



//!SETTINGS

#define EXTENDEDTEST 0
#define PRINT 1 //&& 1EXTENDEDTEST

#define OUTPUT 0


#define type_precision float

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define _1MB 1024*1024
#define _10MB 10*_1MB
#define _1GB 1024*1024*1024

#define one_oversqrt2 0.7071067811865475244008443621048490392848359376884740



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
    int mpi_id;
    int mpi_num_threads;

    bool use_fake_files;



    bool storeBin;

    float minPdisp;
    float minPstore;
    float minR2disp;
    float minR2store;
    bool storePInd;
    bool disp_cov;
    bool dosages;
    int model;//recessive additive dominant etc
    int modelR;

    string fnameAL;
    string fnameAR;
    string fnameY;
    string fnameOutFiles;
    string fname_excludelist;
    string fname_dosages;

    bool use_interactions;
    bool keep_depVar;
    bool use_multiple_interaction_sets;
    string fname_interactions;
    int ini_IDX_interactions;
    int end_IDX_interactions;

    bool doublefileType;

    bool ForceCheck;

    int limit_m;
    int limit_t;
    int limit_n;
};


struct Outputs
{

    double duration;
    double gflops;

    double acc_real_loady;
    double acc_real_loadxr;
    double acc_real_storeb;

    double acc_real_innerloops;

    double acc_loady;
    double acc_loadxr;
    double acc_storeb;

    double acc_other;

    double acc_gemm;

    double acc_firstAR;
    double acc_firstY;


    double acc_stl;

    double acc_sbr;
    double acc_str;
    double acc_solve;

    double acc_scorrect;
    double acc_inner_scorrect;

    double acc_stats;

    unsigned int total_sig_results;


};



struct result
{
    float T;
    long double P;
    float SE;
    float B;
    uint8_t AL_name_idx;
    int AR_name_idx;
};

class resultH
{
    public:
    int Y_name_idx;
    int ARoffset;
    uint16_t nUsed;
    float nUsedPct;
    float R2;
    list< result > res_p;
};




#endif // UTILITY_H_INCLUDED
