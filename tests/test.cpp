//export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_NESTED=true
#include <unistd.h>
#include <getopt.h>



#include "../src/Definitions.h"
#include "../src/Algorithm.h"



int main(int argc, char *argv[] )
{
    struct Settings params;

    params.ForceCheck = true;


    //!default params
    params.r = 1;
    params.threads = 1;

    omp_set_num_threads(params.threads);
    blas_set_num_threads(params.threads);

    Algorithm alg;

    params.use_fake_files = true;
    int iters = 10;
    int max_threads = 2;


    for (int th = 0; th < max_threads; th++)
    {
        params.threads = th;

        params.n=10; params.l=4;  params.r=1;
        params.t=16; params.tb=1; params.m=16; params.mb=1;

        struct Outputs out = {0};
        for (int i = 0; i < iters; i++)
        {
            alg.solve(params, out, P_NEQ_B_OPT_MD);
        }
        cout << endl;
        /******************************/
        params.n=10; params.l=4;  params.r=2;
        params.t=16; params.tb=4; params.m=16; params.mb=4;


        for (int i = 0; i < iters; i++)
        {
            alg.solve(params, out, P_NEQ_B_OPT_MD);
        }
        cout << endl;
        /******************************/
        params.n=10; params.l=4;  params.r=2;
        params.t=16; params.tb=5; params.m=16; params.mb=3;


        for (int i = 0; i < iters; i++)
        {
            alg.solve(params, out, P_NEQ_B_OPT_MD);
        }
        cout << endl;
        /******************************/
        params.n=10; params.l=4;  params.r=2;
        params.t=4; params.tb=4; params.m=4; params.mb=4;


        for (int i = 0; i < iters; i++)
        {
            alg.solve(params, out, P_NEQ_B_OPT_MD);
        }
        cout << endl;

    }

    cout << "\nTest finished succesfully\n";




    return 0;
}
