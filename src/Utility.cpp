#include "Utility.h"


double gemm_flops(double m, double n, double k, int sum)
{
    return  (double)((2.0*m / 1000.0*n / 1000.0*k / 1000.0 +
                      sum*m / 1000.0*k / 1000.0 / 1000.0));
}


void myassert(int cond, string msg, int var)
{
    if (!cond)
    {
        cout << msg << " "<< var << endl;
        exit(1);
    }
}


void myassert(int cond, string msg)
{
    if (!cond)
    {
        cout << msg << endl;
        exit(1);
    }
}





void re_random_vec(type_precision* vec, int size)
{
    for (int i = 0; i < size; i++)
    {
        vec[i] = (type_precision)rand() / (type_precision)RAND_MAX;
    }
}





type_precision* replicate_vec(type_precision* old, int size)
{
    type_precision* vec =
        (type_precision*)malloc(size * sizeof(type_precision));

    if (vec==0)
    {
        cout << "\nNot enough RAM! "
             << (int)(size*sizeof(type_precision)/1024/1024) << "MB\n";
        //system("pause");
        exit(1);
    }

    memcpy( (type_precision*)vec,
            (type_precision*)old,
            size * sizeof(type_precision) );

    return vec;
}


void replace_nans_avgs(int vec_blocksize, type_precision* vec,
                       int rows , int cols, list<long int>* ar_nan_idxs)
{
    type_precision* ones = new type_precision[rows];

    for (int i = 0; i < rows; i++)
    {
        ones[i] = 1.0;
    }

    for (int i = 0; i < vec_blocksize * cols; i++)
    {
        type_precision sum = cblas_sdot(rows, ones, 1, &vec[i * rows], 1);
        type_precision avg = sum / (rows - ar_nan_idxs[i].size());
        for (list<long int>::iterator it = ar_nan_idxs[i].begin();
             it != ar_nan_idxs[i].end();
             it++)
        {
            int idx = i * rows + (*it);
            vec[idx] = avg;
        }
    }

    delete []ones;
}


void replace_nans(list<long int>* indexs_vec, int vec_blocksize,
                  type_precision* vec, int rows , int cols)
{
    //#pragma omp parallel default(shared)
    if (indexs_vec)
    {
        for (int k = 0; k < vec_blocksize; k++)
        {
            indexs_vec[k].clear();
        }
    }

    for (int k = 0; k < vec_blocksize; k++)
    {
        //#pragma omp for nowait  schedule(static)
        for (int i = 0; i < cols; i++)      // go through all columns
        {
             for (int j = 0; j < rows; j++) // move over the rows of this column
             {
                int idx = k * cols * rows + i * rows + j;

                //cout << idx;
                if (idx >= rows * cols * vec_blocksize)
                {
                    exit(1);
                }

                if (isnan(vec[idx] ))
                {
                    vec[idx] = 0;
                    if (indexs_vec)
                    {
                        //this col had this rows with nans
                        indexs_vec[k].push_back(j);
                    }
                }
             }
        }
        if (indexs_vec)
        {
            indexs_vec[k].sort();
            indexs_vec[k].unique();
        }
    }

    if (cols > 1)
    {
        replace_with_zeros(indexs_vec, vec, rows, cols, vec_blocksize);
    }
}


void replace_with_zeros(list<long int>* indexs, type_precision* vec,
                        int n, int cols, int vec_block_count)
{
    if (indexs)
    {
        for (int i = 0; i < vec_block_count; i++)
        {
            int idx;
            for (int j = 0; j < cols; j++)
            {
                for (list<long int>::iterator it = indexs->begin();
                     it != indexs->end();
                     it++)
                {
                    idx = i * cols * n + j * n + (*it);
                    vec[idx] = 0;
                }
            }
        }
    }
}


void matlab_print_matrix(string name, int m, int n, type_precision* A)
{//fix for row major
    if (PRINT)
    {
        cout << endl << name << " = [\n";
        int i;
        int j;
        int index=0;
        for (j=0;j<m;j++)
        {
          for (i=0;i<n;i++)
          {
              if (i != 0)
              {
                  printf(", \t");
              }
              index = j + i * m;
              if (A[index] == 0.0)
              {
                  cout << " 0.000";
              }
              else
              {
                  if (A[index] < 0.0)
                  {
                      cout << std::setprecision(2) << A[index];
                  }
                  else
                  {
                      cout << " "<<std::setprecision(2) << A[index];
                  }
              }
          }
          cout << " ; \n";
        }
        cout << "];\n";
    }
}

float getTvalue(float pval)
{
//    cout << pval << endl;
    float T = boost::math::erf_inv((long double)(1.0 - pval)) * sqrt(2.0);
//    cout << T << endl;
    return T;
}


void float32(float* out, const uint16_t in)
{
    uint32_t t1 = in;
    t1 <<= 16;        // convert to float again
    t1 &= 0xffff0000; // keep sign and exponent only
    *((uint32_t*)out) =  t1;
}


void float16(uint16_t &out, const float in)
{
    uint32_t* intp = (uint32_t*)&(in);
    uint32_t inu = *(intp);
    inu >>= 16; // left most 16 are kept after next truncation
    out = inu;  // truncates zeroes away
}
