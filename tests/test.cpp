#include <unistd.h>
#include <getopt.h>



#include "../src/Definitions.h"
#include "../src/Algorithm.h"

void print_output(struct Outputs out, float gemm_gflopsPsec)
{
    std::cout << std::fixed;


    cout << "\nDuration:\t\t" << out.duration <<" \t"<<  out.duration/out.duration*100 << "\n";

    double missing = (out.duration - (out.acc_loady+out.acc_loadxr
                                                            +out.acc_real_innerloops+out.acc_gemm));

    cout << "Missing Sec:\t\t" << missing <<" \t"<<  missing/out.duration*100 << "\n";

    missing = (out.acc_real_innerloops - (out.acc_solve+out.acc_sbr +out.acc_stl+out.acc_str+out.acc_other+out.acc_stats+out.acc_storeb));

    cout << endl;

    cout << "GEMM Sec:\t\t" << out.acc_gemm <<" \t"<<  out.acc_gemm/out.duration*100 << "\n";

    cout << "InnerLoops  Sec:\t" << out.acc_real_innerloops <<" \t"<<  out.acc_real_innerloops/out.duration*100 << "\n";
    cout << "Stl Sec:\t\t" << out.acc_stl <<" \t"<<  out.acc_stl/out.duration*100 << "\n";
    cout << "Str Sec:\t\t" << out.acc_str <<"\t"<<  out.acc_str/out.duration*100 << "\n";
    cout << "Sbr Sec:\t\t" << out.acc_sbr <<"\t"<<  out.acc_sbr/out.duration*100 << "\n";
    cout << "Solve Sec:\t\t" << out.acc_solve <<" \t"<<  out.acc_solve/out.duration*100 << "\n";
    cout << "Other Sec:\t\t" << out.acc_other <<" \t"<<  out.acc_other/out.duration*100 << "\n";
    cout << "Stats Sec:\t\t" <<  out.acc_stats <<"\t"<<  out.acc_stats/out.duration*100 << "\n";
    cout << "Missing Sec:\t\t" << missing <<"\t"<<  missing/out.duration*100 << "\n";

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

    cout << endl<< "Duration:"<< out.duration << endl;
    cout << endl<< "GFLOPS:"<< out.gflops << endl;
    cout << "GFLOPS/s: " << (out.gflops/out.duration) << endl;
    cout  <<"Perf: " << (out.gflops/out.duration)/gemm_gflopsPsec << endl ;




}

int main(int argc, char *argv[] )
{
    struct Settings params;



    omp_set_nested(false);
    omp_set_dynamic(false);

    Algorithm alg;


    //!default params
    alg.applyDefaultParams(params);

    params.use_fake_files = true;



    int max_threads = 2;

    params.threads = max_threads;

    omp_set_num_threads(params.threads);
    blas_set_num_threads(params.threads);

    cout << "Perfromance Test\n" << flush;

    double duration, gemmgflops, gemm_gflopsPsec;
    cpu_benchmark((1*1024.0),3,duration,gemmgflops);
    gemm_gflopsPsec = gemmgflops/duration;
    cout << "\nGEMM GFLOPS/s " << gemm_gflopsPsec << endl;

    struct Outputs out2 = {0};
    int factor = 10;
    params.n=4000; params.l=5;  params.r=1;
    params.t=100*factor; params.tb=min(1000,100*factor); params.m=100*factor; params.mb=min(200,100*factor);
    alg.solve(params, out2, P_NEQ_B_OPT_MD);

    print_output(out2, gemm_gflopsPsec);


    cout << "\nDone\n";
    cout << "\nMisc Tests\n" << flush;

    //!default params
    alg.applyDefaultParams(params);

    params.minPstore = 0.1;
    params.minPdisp = 0.05;
    params.minR2store = 0.001;
    params.minR2disp = 0.001;


    params.r = 2;
    params.fnameAL="examples/XL";
    params.fnameAR="examples/XR";
    params.fnameY="examples/Y";
    params.fnameOutFiles="resultsSig";


    for(int th = 0; th < max_threads; th++)
    {
        cout << "****" << flush;
    }
    struct Outputs out23 = {0};
    cout << "**\n*" << flush;
    alg.solve(params, out23, P_NEQ_B_OPT_MD);

    print_output(out23, gemm_gflopsPsec);


    cout << "*" << flush;

    alg.applyDefaultParams(params);
    params.use_fake_files = true;
    params.ForceCheck = true;

    max_threads = 2;
    int iters = 10;

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
        params.t=16; params.tb=3; params.m=16; params.mb=5;


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
