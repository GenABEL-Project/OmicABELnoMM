#ifndef UTILITY_H_INCLUDED
#define UTILITY_H_INCLUDED

#include "Definitions.h"
#include <string.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <sys/time.h>
#include <vector>
#include <list>
#include <bitset>
//


double gemm_flops(double m, double n, double k,int sum);
void myassert(int cond,string msg);
void myassert(int cond,string msg,int var);
type_precision* random_vec(int size);
void re_random_vec(type_precision* vec, int size);
inline void copy_vec(type_precision*old, type_precision* new_vec, int size)
{
    memcpy( (type_precision*)new_vec, (type_precision*)old, size * sizeof(type_precision) );
}
type_precision* replicate_vec(type_precision*old, int size);
void matlab_print_matrix(string name,int m,int n,type_precision* A);
void cpu_benchmark(int n,int samples, double &duration, double &gflops);

void re_random_vec_nan(type_precision* vec, int size);

void replace_with_zeros(list<long int>* indexs, type_precision* vec, int n, int r,int block_count);
void replace_nans(list<long int>* indexs, int blocksize, type_precision* vec, int rows , int cols);
void replace_nans_avgs(int vec_blocksize, type_precision* vec, int rows , int cols,list<long int>* ar_nan_idxs);


void float32(float* out, const uint16_t in);
void float16(uint16_t &out, const float in);

float getTvalue(float pval);

#define cputime_type struct timeval

    #define get_ticks(t)\
    {\
          gettimeofday(&t, NULL);\
    }

    #define ticks2sec(end_time, start_time)\
          ((int)(end_time.tv_sec  - start_time.tv_sec) * 1e6 +\
                       (int)(end_time.tv_usec - start_time.tv_usec))/1000000.0\


#endif // UTILITY_H_INCLUDED
