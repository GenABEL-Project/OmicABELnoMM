//export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_NESTED=true
#include <unistd.h>
#include <getopt.h>



#include "../src/Definitions.h"
#include "../src/Algorithm.h"

void print_output(struct Outputs out)
{
    std::cout << std::fixed;


    cout << "\nDuration:\t\t" << out.duration <<" \t"<<  out.duration/out.duration*100 << "\n";

    double missing = (out.duration - (out.acc_loady+out.acc_loadxr
                                                            +out.acc_real_innerloops+out.acc_gemm));

    cout << "Missing Sec:\t\t" << missing <<" \t"<<  missing/out.duration*100 << "\n";

    missing = (out.acc_real_innerloops - (out.acc_solve+out.acc_sbr +out.acc_stl+out.acc_str+out.acc_other));

    cout << endl;

    cout << "GEMM Sec:\t\t" << out.acc_gemm <<" \t"<<  out.acc_gemm/out.duration*100 << "\n";

    cout << "InnerLoops  Sec:\t" << out.acc_real_innerloops <<" \t"<<  out.acc_real_innerloops/out.duration*100 << "\n";
    cout << "\tStl Sec:\t" << out.acc_stl <<" \t"<<  out.acc_stl/out.duration*100 << "\n";
    cout << "\tStr Sec:\t" << out.acc_str <<"\t"<<  out.acc_str/out.duration*100 << "\n";
    cout << "\tSbr Sec:\t" << out.acc_sbr <<"\t"<<  out.acc_sbr/out.duration*100 << "\n";
    cout << "\tSolve Sec:\t" << out.acc_solve <<" \t"<<  out.acc_solve/out.duration*100 << "\n";
    cout << "\tOther Sec:\t" << out.acc_other <<" \t"<<  out.acc_other/out.duration*100 << "\n";
    cout << "\tMissing Sec:\t" << missing <<"\t"<<  missing/out.duration*100 << "\n";

    //cout << "RQy Sec:\t" << out.acc_real_innerloops-out.acc_gemm <<" \t"<<  (out.acc_real_innerloops-out.acc_gemm)/out.duration*100 << "\n";
    cout << endl;
    cout << "Betas Sec:\t\t" << out.acc_storeb <<" \t"<<  out.acc_storeb/out.duration*100 << "\n";
    cout << "LoadXR Sec:\t\t" << out.acc_loadxr <<" \t"<<  out.acc_loadxr/out.duration*100 << "\n";
    cout << "LoadY Sec:\t\t" << out.acc_loady <<" \t"<<  out.acc_loady/out.duration*100 << "\n";
    //cout << endl;
//    cout << "FirstXR Sec:\t" << out.acc_firstAR <<" \t"<<  out.acc_firstAR/out.duration*100 << "\n";
//    cout << "FirstY Sec:\t" << out.acc_firstY <<" \t"<<  out.acc_firstY/out.duration*100 << "\n";
    //cout << endl;


//    cout << "Real StoreB Sec:\t" << out.acc_real_storeb <<" \t"<<  out.acc_real_storeb/out.duration*100 << "\n";
//    cout << "Real LoadXR Sec:\t" << out.acc_real_loadxr <<" \t"<<  out.acc_real_loadxr/out.duration*100 << "\n";
//    cout << "Real LoadY  Sec:\t" << out.acc_real_loady <<" \t"<<  out.acc_real_loady/out.duration*100 << "\n";




}

int main(int argc, char *argv[] )
{
    struct Settings params;


    params.ForceCheck = false;

    //!default params
    params.r = 1;
    params.threads = 1;



    omp_set_nested(false);
    omp_set_dynamic(false);

    Algorithm alg;

    params.use_fake_files = true;
    int iters = 10;
    int max_threads = 2;

    params.threads = max_threads;

    omp_set_num_threads(params.threads);
    blas_set_num_threads(params.threads);

    cout << "Perfromance Test\n" << flush;

    double duration, gemmgflops, gemm_gflopsPsec;
    cpu_benchmark((2*1024.0),2,duration,gemmgflops);
    gemm_gflopsPsec = gemmgflops/duration;
    cout << "\nGEMM GFLOPS/s " << gemm_gflopsPsec << endl;

    struct Outputs out2 = {0};
    params.n=2000; params.l=5;  params.r=1;
    params.t=2000; params.tb=1000; params.m=2000; params.mb=100;
    alg.solve(params, out2, P_NEQ_B_OPT_MD);

    cout <<endl<< "Duration:"<< out2.duration << endl;
    cout <<endl<< "GFLOPS:"<< out2.gflops << endl;
    cout << "GFLOPS/s: " << (out2.gflops/out2.duration) << endl;
    cout  <<"Perf: " << (out2.gflops/out2.duration)/gemm_gflopsPsec << endl ;

    print_output(out2);


    cout << "\nDone\n";
    cout << "\nMisc Tests\n" << flush;
    params.ForceCheck = true;

    for(int th = 0; th < max_threads; th++)
    {
        cout << "****" << flush;
    }

    cout << "*\n*" << flush;

    max_threads = 2;

    for (int th = 1; th < max_threads+1; th++)
    {

        omp_set_num_threads(th);
        blas_set_num_threads(th);
        params.threads = th;

        params.n=10; params.l=4;  params.r=1;
        params.t=16; params.tb=1; params.m=16; params.mb=1;



        struct Outputs out = {0};
        for (int i = 0; i < iters; i++)
        {
            alg.solve(params, out, P_NEQ_B_OPT_MD);
        }
        cout << "*" << flush;
        /******************************/
        params.n=10; params.l=4;  params.r=2;
        params.t=16; params.tb=4; params.m=16; params.mb=4;


        for (int i = 0; i < iters; i++)
        {
            alg.solve(params, out, P_NEQ_B_OPT_MD);
        }
        cout << "*" << flush;
        /******************************/
        params.n=10; params.l=4;  params.r=2;
        params.t=16; params.tb=5; params.m=16; params.mb=3;


        for (int i = 0; i < iters; i++)
        {
            alg.solve(params, out, P_NEQ_B_OPT_MD);
        }
        cout << "*" << flush;
        /******************************/
        params.n=10; params.l=4;  params.r=2;
        params.t=4; params.tb=4; params.m=4; params.mb=4;


        for (int i = 0; i < iters; i++)
        {
            alg.solve(params, out, P_NEQ_B_OPT_MD);
        }
        cout << "*" << flush;

    }

    cout << "\nTest finished succesfully\n";




    return 0;
}
