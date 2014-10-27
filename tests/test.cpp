#include <unistd.h>
#include <getopt.h>
#include <string>
#include <vector>

using namespace std;


#include "../src/Definitions.h"
#include "../src/Algorithm.h"

void print_params(struct Settings params)
{
    cout << "Threads:" << params.threads << " ";
    cout << "n:" << params.n << " ";
    cout << "l:" << params.l << " ";
    cout << "r:" << params.r << " ";
    cout << "m:" << params.m << " ";
    cout << "mb:" << params.mb << " ";
    cout << "t:" << params.t << " ";
    cout << "tb:" << params.tb << endl;
}

void print_output(struct Outputs out, float gemm_gflopsPsec)
{
    std::cout << std::fixed;

    cout << "\nTotal Results: \t" << out.total_sig_results << endl;

    cout << "\nDuration:\t\t" << out.duration << " \t"
         <<  out.duration/out.duration * 100 << "\n";

    double missing = (out.duration - (out.acc_loady+out.acc_loadxr
                                      + out.acc_real_innerloops
                                      + out.acc_gemm + out.acc_scorrect));

    cout << "Missing Sec:\t\t" << missing << " \t"
         <<  missing/out.duration * 100 << "\n";

    missing = (out.acc_real_innerloops -
               (out.acc_solve + out.acc_sbr + out.acc_stl + out.acc_str
                + out.acc_other + out.acc_stats + out.acc_storeb
                + out.acc_inner_scorrect));


    cout << endl;

    cout << "GEMM Sec:\t\t" << out.acc_gemm << " \t" <<
        out.acc_gemm/out.duration * 100 << "\n";
    cout << "SGenCorr Sec:\t\t" << out.acc_scorrect << "\t"
         <<  out.acc_scorrect/out.duration * 100 << "\n";

    cout << "InnerLoops  Sec:\t" << out.acc_real_innerloops << " \t"
         <<  out.acc_real_innerloops/out.duration * 100 << "\n";
    cout << "Stl Sec:\t\t" << out.acc_stl << " \t"
         << out.acc_stl/out.duration * 100 << "\n";
    cout << "Str Sec:\t\t" << out.acc_str << "\t"
         <<  out.acc_str/out.duration * 100 << "\n";
    cout << "Sbr Sec:\t\t" << out.acc_sbr << "\t"
         <<  out.acc_sbr/out.duration * 100 << "\n";
    cout << "Solve Sec:\t\t" << out.acc_solve << " \t"
         << out.acc_solve/out.duration * 100 << "\n";
    cout << "Other Sec:\t\t" << out.acc_other << " \t"
         << out.acc_other/out.duration * 100 << "\n";
    cout << "Stats Sec:\t\t" <<  out.acc_stats << "\t"
         << out.acc_stats/out.duration * 100 << "\n";
    cout << "SCorr Sec:\t\t" <<  out.acc_inner_scorrect << "\t"
         <<  out.acc_inner_scorrect/out.duration * 100 << "\n";

    cout << "Missing Sec:\t\t" << missing << "\t"
         <<  missing/out.duration * 100 << "\n";

    // cout << "RQy Sec:\t" << out.acc_real_innerloops-out.acc_gemm
    // << " \t" <<  (out.acc_real_innerloops-out.acc_gemm)
    // / out.duration * 100 << "\n";
    cout << endl;
    cout << "Betas Sec:\t\t" << out.acc_storeb << " \t"
         <<  out.acc_storeb/out.duration * 100 << "\n";
    cout << "LoadXR Sec:\t\t" << out.acc_loadxr << " \t"
         << out.acc_loadxr/out.duration * 100 << "\n";
    cout << "LoadY Sec:\t\t" << out.acc_loady << " \t"
         << out.acc_loady/out.duration * 100 << "\n";
    // cout << endl;
    // cout << "FirstXR Sec:\t" << out.acc_firstAR << " \t"
    // <<  out.acc_firstAR/out.duration * 100 << "\n";
    // cout << "FirstY Sec:\t" << out.acc_firstY << " \t"
    // <<  out.acc_firstY/out.duration * 100 << "\n";
    // cout << endl;


    // cout << "Real StoreB Sec:\t" << out.acc_real_storeb << " \t"
    // <<  out.acc_real_storeb/out.duration * 100 << "\n";
    // cout << "Real LoadXR Sec:\t" << out.acc_real_loadxr << " \t"
    // <<  out.acc_real_loadxr/out.duration * 100 << "\n";
    // cout << "Real LoadY  Sec:\t" << out.acc_real_loady << " \t"
    // <<  out.acc_real_loady/out.duration * 100 << "\n";

    cout << endl << "Duration:" << out.duration << endl;
    cout << endl << "GFLOPS:" << out.gflops << endl;
    cout << "GFLOPS/s: " << (out.gflops/out.duration) << endl;
    cout << "Perf: " << (out.gflops/out.duration)/gemm_gflopsPsec << endl;
}


int main(int argc, char *argv[])
{
    struct Settings params;

    //!MPI
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //printf("SIZE = %d RANK = %d\n",size,rank);


    omp_set_nested(false);
    omp_set_dynamic(false);

    Algorithm alg;

    int max_threads = 1;

    params.threads = max_threads;

    omp_set_num_threads(params.threads);
    blas_set_num_threads(params.threads);

    cout << "Perfromance Test\n" << flush;

    double duration, gemmgflops, gemm_gflopsPsec;
    cpu_benchmark((1024.0), 5, duration, gemmgflops);
    gemm_gflopsPsec = gemmgflops/duration;
    cout << "\nGEMM GFLOPS/s " << gemm_gflopsPsec << endl;

    struct Outputs out = {0};
    //    params.use_fake_files = false;
    //    !default params
    //    alg.applyDefaultParams(params);
    //    int factor = 0;//need to repair fakefiles?
    //    params.n=2000; params.l=3;  params.r=1;
    //    params.t=1000; params.tb=min(1000, params.t); params.m=1000;
    //    params.mb=min(200, params.m);
    //    print_params(params);
    //    alg.solve(params, out, P_NEQ_B_OPT_MD);
    //    print_output(out, gemm_gflopsPsec);
    //
    //
    //    cout << "\nDone\n";


    //!-------------------------------------
    //!default params
    alg.applyDefaultParams(params);
    params.threads = max_threads;

    params.minPstore = 0.1;
    params.minPdisp = 0.00005;
    params.minR2store = 0.00001;
    params.minR2disp = 0.000001;
    params.disp_cov = true;

    params.mpi_id = rank;
    params.mpi_num_threads = size;


    params.r = 2;
    params.fnameOutFiles = "examples/results/normal";
    params.fnameAL = "examples/XL";
    params.fnameAR = "examples/XR";
    params.fnameY = "examples/Y";

    params.limit_t = 1000;
    params.limit_m = params.limit_t * params.r;

    alg.solve(params, out, P_NEQ_B_OPT_MD);
    print_output(out, gemm_gflopsPsec);

    //!MPI
    MPI_Finalize();

    // exit(0);

    cout << "\nInteraction Tests\n" << flush;
    //!-------------------------------------
    //!-------------------------------------
    //!interactions
    params.fnameOutFiles = "examples/results/normal";
    params.fnameAL = "examples/interactions/XL";
    params.fnameAR = "examples/interactions/XR";
    params.fnameY = "examples/interactions/Y";

    params.use_interactions = true;
    params.r = 2;
    params.limit_t = 50;
    params.limit_m = params.limit_t * params.r;
    // params.limit_n = 1000;
    params.fname_excludelist = "examples/exclude_individuals.txt";
    string source_path = "examples/interactions/INT";
    string out_path[] = {"examples/results/single_inter_",
                         "examples/results/multi_inter_"};
    string out_keep[] = {"", "keep_"};
    string dosagename[] = {"", "add_", "dom_", "res_", "mylin_", "myadd_"};
    int dosage_r[] = {1, 2, 2, 2, 1, 2};
    //!--------------No Model-----------------------
    string test_names[] = {"0", "1", "R", "R1", "R0", "10", "11", "011", "1RR"};
    string keep_test_names[] = {"R", "R0", "RR"};
    uint32_t exected_res_single[] = {0, (uint32_t)params.limit_t, 0,
                                     (uint32_t)params.limit_t, 0, 0, 0,
                                     0, (uint32_t)params.limit_t};
    uint32_t exected_res_mult[] = {0, (uint32_t)params.limit_t, 0,
                                   (uint32_t)params.limit_t, 0,
                                   (uint32_t)params.limit_t,
                                   2 *(uint32_t)params.limit_t,
                                   (uint32_t)2 *params.limit_t,
                                   (uint32_t)params.limit_t};
    uint32_t exected_res_single_keep[] = {(uint32_t)params.limit_t, 0,
                                          (uint32_t)params.limit_t};
    uint32_t exected_res_mult_keep[] = {(uint32_t)params.limit_t,
                                        (uint32_t)params.limit_t,
                                        (uint32_t)2 *params.limit_t};


    for (int j = 0; j < 2; j++ )
    {
        params.use_multiple_interaction_sets = j;
        for (int k = 0; k < 2; k++)
        {
            params.keep_depVar = k;
            for (int d = 0; d < 6; d++)
            {
                params.dosages = d;
                params.model = d - 1;
                params.r=dosage_r[d];
                params.fname_dosages = "examples/dosages_";
                params.fname_dosages += std::to_string(params.r);
                params.fname_dosages += ".txt";

                if (params.keep_depVar)
                {
                    for (int i = 0; i < 3; i++ )
                    {
                        params.fnameOutFiles = out_path[j] +
                            out_keep[k] + dosagename[d] + keep_test_names[i];
                        cout << params.fnameOutFiles << endl;
                        params.fname_interactions = source_path +
                            keep_test_names[i];
                        memset(&out, 0, sizeof(out));
                        alg.solve(params, out, P_NEQ_B_OPT_MD);
                        cout << "\nTotal Results: \t"
                             << out.total_sig_results << endl;
                        cout << "=========================================\n";
                        if (d != 3)
                        {
                            if (j == 0)
                            {
                                myassert(exected_res_single_keep[i] <=
                                         out.total_sig_results,
                                         "Unexpected number of Results, check manually!");
                            }
                            else
                            {
                                myassert(exected_res_mult_keep[i] <=
                                         out.total_sig_results,
                                         "Unexpected number of Results, check manually!");
                            }
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < 9; i++ )
                    {
                        //cout << "=========================================\n";
                        params.fnameOutFiles=out_path[j] + out_keep[k] +
                            dosagename[d] + test_names[i];
                        cout << params.fnameOutFiles << endl;
                        params.fname_interactions=source_path + test_names[i];

                        memset(&out, 0, sizeof(out));
                        alg.solve(params, out, P_NEQ_B_OPT_MD);
                        cout << "\nTotal Results: \t"
                             << out.total_sig_results << endl;
                        cout << "=========================================\n";

                        if (d != 3)
                        {
                            if (j == 0)
                                myassert(exected_res_single[i] <=
                                         out.total_sig_results,
                                         "Unexpected number of Results, check manually!");
                            else
                                myassert(exected_res_mult[i] <=
                                         out.total_sig_results,
                                         "Unexpected number of Results, check manually!");
                        }
                    }
                }
                cout << "========================================\n";
            }
            cout << "===============KEEP=============\n";
        }
        cout << "==================MULTI=========================\n";
    }
    //!-------------------------------------
    //!interactions END
    //!-------------------------------------

    cout << "\nDosage Tests\n" << flush;
    //!-------------------------------------
    //!dosages alone
    //!-------------------------------------
    string out_dos_path[] = {"examples/results/dosages_",
                             "examples/results/dosages_excl_"};

    for (int j = 0; j < 2; j++)
    {
        params.fname_excludelist = "";
        if (j)
            params.fname_excludelist = "examples/exclude_individuals.txt";

        for (int d = 0; d < 6; d++)
        {
            params.dosages = d;
            params.model = d-1;
            params.r = dosage_r[d];
            params.fname_dosages = "examples/dosages_";
            params.fname_dosages += std::to_string(params.r);
            params.fname_dosages += ".txt";

            params.fnameOutFiles = out_dos_path[j] + dosagename[d];
            cout << params.fnameOutFiles << endl;

            memset(&out, 0, sizeof(out));
            alg.solve(params, out, P_NEQ_B_OPT_MD);
            cout << "\nTotal Results: \t" << out.total_sig_results << endl;
            cout << "=========================================\n";

            if (d != 3)
            {
                myassert((uint32_t)params.limit_t <= out.total_sig_results,
                         "Unexpected number of Results, check manually!");
            }
        }
    }


    exit(0);
    //!-------------------------------------
    alg.applyDefaultParams(params);
    memset(&out, 0, sizeof(out));
    params.use_fake_files = true;
    params.ForceCheck = true;

    max_threads = 2;
    int iters = 10;

    cout << "Misc tests" << endl;

    for (int th = 1; th < max_threads + 1; th++)
    {
        omp_set_num_threads(th);
        blas_set_num_threads(th);
        params.threads = th;

        params.n=10; params.l=4;  params.r=1;
        params.t=16; params.tb=1; params.m=16; params.mb=1;

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
