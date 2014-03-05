//export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_NESTED=true
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
            {"ngpred", required_argument, 0, 'n'},
            {"thr",    required_argument, 0, 't'},
            {"mem",    required_argument, 0, 'm'},
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here.
        int option_index = 0;

        c = getopt_long(argc, argv, "c:p:g:o:n:t:m:", long_options, &option_index);


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
            params.fnameOutB = string(optarg);

            pos = string(optarg).find(".");
            if (pos != string::npos)
                params.fnameOutB = string(optarg).substr(0, pos);

            cout << "using -o with output file " << optarg << endl;
            break;

        case 't':
            params.threads = atoi(optarg);
            params.threads = min(max(params.threads, 1), 32);
            cout << "using -t with value " << params.threads << endl;
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

    parse_params(argc, argv, params);


    //omp_set_nested(false);
    //mkl_set_dynamic(false);
    //omp_set_dynamic(false);

    omp_set_num_threads(params.threads);
    blas_set_num_threads(params.threads);

    params.use_fake_files = false;

//    if (params.use_fake_files)
//    {
//        int multiplier = 1024;
//
//        params.n = 20 * multiplier;
//        params.l = 8;
//        params.r = 2;
//
//        params.t  = (16 * multiplier, 16 * multiplier);
//        params.tb = (4 * multiplier,4 * multiplier);
//
//        params.m  = (4 * multiplier,4 * multiplier);
//        params.mb = (4 * multiplier,4 * multiplier);
//    }


    Algorithm alg;
    struct Outputs out = {0};
    if (params.use_fake_files)
    {
        for (int i = 0; i < 100; i++)
            alg.solve(params, out, P_NEQ_B_OPT_MD);
    }
    else
    {
        alg.solve(params, out, P_NEQ_B_OPT_MD);
    }

    return 0;
}
