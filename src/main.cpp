//export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_NESTED=true

#include <boost/math/special_functions/erf.hpp>

#include <unistd.h>
#include <getopt.h>

#include "Definitions.h"
#include "Algorithm.h"





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


void parse_params(int argc, char *argv[], struct Settings &params )
{
    int c;

    unsigned pos;

    bool phe = false;
    bool snp = false;
    bool cov = false;
    bool bout = false;


    while (1)
    {
        static struct option long_options[] =
        {   //-c XL --snp XR -p Y -b B -r 2 -t 2
            // These options set a flag.
            //{"", no_argument,       &verbose_flag, 1},
            //{"",   no_argument,       &verbose_flag, 0},
            // These options don't set a flag.
            //We distinguish them by their indices.
            {"phe",    required_argument, 0, 'p'},
            {"geno",   required_argument, 0, 'g'},
            {"cov",    required_argument, 0, 'c'},
            {"out",    required_argument, 0, 'o'},
            {"ngpred", required_argument, 0, 'n'},//r
            {"thr",    required_argument, 0, 't'},
            {"mem",    required_argument, 0, 'm'},
            {"excl",    required_argument, 0, 'x'},
            {"sigth",    required_argument, 0, 's'},//pvalthreashold
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here.
        int option_index = 0;

        c = getopt_long(argc, argv, "c:p:g:o:n:t:m:x:", long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        if (!optarg)
        {
            cout << "\nerror with argument parameter " << c << endl;
            exit(1);
        }

        switch (c)
        {
        case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
                break;
            printf("Long option %s", long_options[option_index].name);
            if (optarg)
                printf (" with arg %s", optarg);
            printf("\n");
            break;

        case 'p':
            phe = !phe;
            params.fnameY = string(optarg);

            pos = string(optarg).find(".");
            if (pos != string::npos)
                params.fnameY = string(optarg).substr(0, pos);

            cout << "using -p with phenotypes from file " << optarg << endl;
            break;

        case 'g':
            snp = !snp;
            params.fnameAR = string(optarg);

           pos = string(optarg).find(".");
            if (pos != string::npos)
                params.fnameAR = string(optarg).substr(0, pos);

            cout << "using -g with genotype data from file " << optarg << endl;
            break;

         case 'n':
            params.r = atoi(optarg);
            params.r = max(params.r, 1);
            cout << "using -n with value " << params.r << endl;
            break;

        case 'c':
            cov = !cov;
            params.fnameAL = string(optarg);

            pos = string(optarg).find(".");
            if (pos != string::npos)
                params.fnameAL = string(optarg).substr(0, pos);

            cout << "using -c with covariates from file " << optarg << endl;
            break;

        case 'o':
            bout = !bout;
            params.fnameOutFiles = string(optarg);

            pos = string(optarg).find(".");
            if (pos != string::npos)
                params.fnameOutFiles = string(optarg).substr(0, pos);

            cout << "using -o for output files " << optarg << endl;
            break;

        case 't':
            params.threads = atoi(optarg);
            params.threads = min(max(params.threads, 1), 32);
            cout << "using -t with value " << params.threads << endl;
            break;

         case 'x':
            params.fname_excludelist = string(optarg);

            cout << "Excluding ids on " << params.fname_excludelist << endl;
            break;


         case 's':
            params.sig_threshold = atof(optarg);

            cout << "Significance will be set at P-val's below " << params.sig_threshold << endl;
            break;

        case '?':
            /* getopt_long already printed an error message. */
            break;

        default:
            abort();
        }
    }

    /* Print any remaining command line arguments (not options). */
    if (optind < argc)
    {
        printf("non-option ARGV-elements: ");
        while (optind < argc)
            printf("%s ", argv[optind++]);

        putchar('\n');
    }

    if (!bout || !phe || !cov || !snp)
    {
        cout << "Missing arguments! p:" << phe
             << " g:" << snp
             << " c:" << cov
             << " o:"<< bout
             << endl;
        exit(1);
    }
}


int main(int argc, char *argv[] )
{
    struct Settings params;
    params.ForceCheck = false;

//    cout << "Using Arguments: ";
//    for (int i = 0; i < argc; i++)
//    {
//        cout  << argv[i] << " ";
//    }
//    cout << endl;

    //!default params
    params.r = 1;
    params.threads = 1;

    params.fname_excludelist = "";

    parse_params(argc, argv, params);


    //omp_set_nested(false);
    //mkl_set_dynamic(false);
    //omp_set_dynamic(false);

    omp_set_num_threads(params.threads);
    blas_set_num_threads(params.threads);

    params.use_fake_files = false;

    //params.use_fake_files = true;
    if (params.use_fake_files)
    {
        int multiplier = 1000;

        params.n = 4 * multiplier;
        params.l = 3;
        params.r = 1;

        params.t  = 1 * multiplier;
        params.tb = 1 * multiplier;

        params.m  = 1 * multiplier;
        params.mb = 1 * multiplier;
    }


    Algorithm alg;
    struct Outputs out = {0};
    if (params.use_fake_files)
    {
            alg.solve(params, out, P_NEQ_B_OPT_MD);
    }
    else
    {
        alg.solve(params, out, P_NEQ_B_OPT_MD);
    }

    cout << endl;

    return 0;
}
