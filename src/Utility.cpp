#include "Utility.h"

double gemm_flops(double m, double n, double k,int sum)
{
    return  (double)((2.0*m/1000.0*n/1000.0*k/1000.0+sum*m/1000.0*k/1000.0/1000.0));
}



void assert(int cond,string msg)
{
    if(cond < 0)
    {
        cout <<"\nCondition violated on param "<< cond << ": ";
        cout << msg << endl;

        exit(cond);
    }
    else
    {
        if(cond > 0)
        {

//            printf("\n%%WARINING %d in ",cond);
//            printf(msg);
//            printf("\n");
        }
    }

}


type_precision* random_vec(int size)
{

    type_precision* vec = (type_precision*)malloc(size*sizeof(type_precision));
    if(vec==0)
    {
        printf("\nNot enough RAM! %dMB\n",size*sizeof(type_precision)/1024/1024 );
        system("pause");
        exit(1);
    }
    for( int i = 0; i < size; i++)
    {
        vec[i] = (type_precision)rand() / (type_precision)RAND_MAX;
    }
    return vec;
}

void re_random_vec(type_precision* vec, int size)
{


    for( int i = 0; i < size; i++)
    {
        vec[i] = (type_precision)rand() / (type_precision)RAND_MAX;
    }
}

void re_random_vec_nan(type_precision* vec, int size)
{
    //int i;

//    for( i = 0; i < size; i++)
//    {
//        if((type_precision)rand() / (type_precision)RAND_MAX > 0.5)
//            vec[i] = nanf("");
//    }
}

//no allocation!
void copy_vec(type_precision*old, type_precision* new_vec, int size)
{
    memcpy( (type_precision*)new_vec, (type_precision*)old, size * sizeof(type_precision) );
}

type_precision* replicate_vec(type_precision*old, int size)
{

    type_precision* vec = (type_precision*)malloc(size*sizeof(type_precision));
    if(vec==0)
    {
        printf("\nNot enough RAM! %dMB\n",size*sizeof(type_precision)/1024/1024 );
        system("pause");
        exit(1);
    }

    memcpy( (type_precision*)vec, (type_precision*)old, size * sizeof(type_precision) );

    return vec;
}

int replace_nans(list<long int>* indexs, type_precision* vec, int rows , int cols)
{
    int count_nans = 0;
    //#pragma omp parallel default(shared) reduction(+:count_nans)
    {

        //#pragma omp for nowait  schedule(static)
        for( int i = 0; i < cols; i++)//go thru all columns
        {

             for( int j = 0; j < rows; j++)//move over the rows of this column
             {
                int idx = i*rows+j;
                if(isnan( vec[idx] ))
                {
                    vec[idx] = 0;
                    if(indexs)
                    {
                        //this col had this rows with nans
                        indexs[i].push_back(j);
                    }

                    count_nans++;
                }

             }

        }
    }
    return count_nans;
}

void replace_with_zeros(list<long int>* indexs, type_precision* vec, int n, int r,int block_count)
{
    if(indexs)
    {


        //#pragma omp parallel default(shared)
        {
            //#pragma omp for nowait  schedule(dynamic)
            for( int i = 0; i < block_count; i++)
            {
                int idx;
               for( int j = 0; j < r; j++)
                {
                    for (list<long int>::iterator it = indexs->begin(); it != indexs->end(); it++)
                    {
                        idx = i*r*n+j*n+(*it);
                        vec[idx] = 0;
                    }
                }
            }
        }
    }
}

void matlab_print_matrix(string name,int m,int n,type_precision* A)
{//fix for row major
    if(PRINT)
    {
        cout << endl << name << " = [\n";
        int i;
        int j;
        int index=0;
        for(j=0;j<m;j++)
        {
          for(i=0;i<n;i++)
          {
              if(i != 0)
                printf(",\t");
                index = j+i*m;

             printf("%.6g",A[index]);
          }
          cout << " ; \n";
        }
        cout << "];\n";

    }
}


void cpu_benchmark(int n,int samples, int cpu_frequency, double &duration, double &GFLOPS)
{
    type_precision* A = new type_precision[n*n];
    type_precision* B = new type_precision[n*n];
    type_precision* C = new type_precision[n*n];


    cputime_type start_tick, end_tick;
    duration = 9999999999.0;
    int b = 0;

    cout << "\n%%Performing CPU GEMM Benchmark" << endl;

    for(int i = 0; i < samples; i++)
    {
        if(samples < 10 || i%(samples/10)==0)
            cout << "%" << flush;

        re_random_vec(A,n*n);
        re_random_vec(B,n*n);
        re_random_vec(C,n*n);

        get_ticks(start_tick);
        cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,n,n,n,1.0,A,n,B,n,1.0,C,n);
        get_ticks(end_tick);
        duration = min(duration,(double)(ticks2sec(end_tick-start_tick,cpu_freq)));
        int a = 0;
        for(int j = 0; j < n*n ; j++)
        {
            a += A[j]+B[j]+C[j];
        }
        b += a;
    }
    //!2nnn - nn + 2nn (from+c)
    GFLOPS = gemm_flops(n,n,n,0);


    cout << b;

    delete []A;
    delete []B;
    delete []C;

}

