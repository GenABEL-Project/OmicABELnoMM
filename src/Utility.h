#ifndef UTILITY_H_INCLUDED
#define UTILITY_H_INCLUDED

#include "Definitions.h"
#include <string.h>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <vector>
#include <list>

double gemm_flops(double m, double n, double k,int sum);
void assert(int cond,char msg[]);
type_precision* random_vec(int size);
void re_random_vec(type_precision* vec, int size);
void copy_vec(type_precision*old, type_precision* new_vec, int size);
type_precision* replicate_vec(type_precision*old, int size);
void matlab_print_matrix(char name[],int m,int n,type_precision* A);
void cpu_benchmark(int n,int samples, int cpu_frequency, double &duration, double &gflops);
int disk_benchmark(int n, int min_size, int max_size, int setps);

void re_random_vec_nan(type_precision* vec, int size);

void replace_with_zeros(list<long int>* indexs, type_precision* vec, int n, int r,int block_count);
int replace_nans(list<long int>* indexs, type_precision* vec, int rows , int cols);
//int replace_nans(list<long int> &indexs, type_precision* old, type_precision* new_vec, int size);
#ifdef WINDOWS

    #define cputime_type double

    #define get_ticks(var) var=clock()

    #define ticks2sec(ticks,ghz) (double)(ticks)/(CLOCKS_PER_SEC)


#else

    #define cputime_type unsigned long

    #define get_ticks(var) { \
           unsigned int __a, __d; \
           asm volatile("rdtsc" : "=a" (__a), "=d" (__d)); \
           var = ((unsigned long) __a) | (((unsigned long) __d) << 32); \
       } while(0)

    #define ticks2sec(ticks,ghz) ((double)(ticks)/((double)(ghz)*1000.0*1000.0*1000.0))


#endif




#endif // UTILITY_H_INCLUDED
