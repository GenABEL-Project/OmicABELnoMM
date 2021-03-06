//export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_NESTED=true

#include <boost/math/special_functions/erf.hpp>
#include <boost/math/distributions/normal.hpp>

#include <unistd.h>
#include <getopt.h>

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include "Definitions.h"
#include "Algorithm.h"

bool help_request = false;

string helpcmd = "usage: omicabelnomm -c <path/fname> --geno <path/fname> "
    "-p <path/fname> -o <path/fname> -x <path/fname> -n <#SNPcols> "
    "-t <#CPUs> -d <0.0~1.0> -r <-10.0~1.0> -b -s <0.0~1.0> "
    "-e <-10.0~1.0> -i -f";

string helpcmd_expl =
    PACKAGE_NAME " v" PACKAGE_VERSION " \n\t"
"Required: \n\t"
"-p --phe    \t <path/filename> to the inputs containing phenotypes. \n\t"
"-g --geno   \t <path/filename> to the inputs containing genotypes. \n\t"
"-c --cov    \t <path/filename> to the inputs containing covariates. \n\t"
"-o --out    \t <path/filename> to store the output to (used for all .txt"
    " and .ibin & .dbin). \n\n"
"Optional: \n\t"
"-n --ngpred \t <#SNPcols> Number of columns in the geno file that represent"
    " a single SNP. \n\t"
"-t --thr    \t <#CPUs> Number of computing threads to use to speed up"
    " computations. Recommended is 4-8 per node (see MPI). \n\t"
"-x --excl   \t <path/filename> file containing list of individuals to"
    " exclude from input files, (see example file). \n\t"
"-d --pdisp  \t <0.0~1.0> Value to use as maximum threshold for"
    " significance.\n\t"
"\t\t Results with P-values BELOW this threshold will be displayed in"
    " the output .txt file. \n\t"
"-r --rdisp  \t <-10.0~1.0> Value to use as minimum threshold for R^2. \n\t"
"\t\t Results with R^2-values ABOVE this threshold will be displayed in"
    " the output .txt file. \n\t"
"-b --stobin \t Flag that forces to ALSO store results in a smaller binary"
    " format (*.ibin & *.dbin). \n\t"
"-s --psto   \t <0.0~1.0> Results with P-values BELOW this threshold will"
    " be displayed in the output binary files. \n\t"
"-e --rsto   \t <-10.0~1.0> Results with R^2-values ABOVE this threshold"
    " will be stored in the output binary files. \n\t"
"-i --fdcov  \t Flag that forces to include covariates (when its genotype"
    " is significant) as part of the results stored. \n\t"
"-f --fdgen  \t Flag that forces to consider all included results (causes"
    " the analysis to ignore ALL threshold values). \n\t"
"-j --additive  \t Flag that runs the analysis with an Additive Model"
    " with (2*AA, 1*AB, 0*BB) effects. \n\t"
"-k --dominant  \t Flag that runs the analysis with a Dominant Model"
    " with (1*AA, 1*AB, 0*BB) effects. \n\t"
"-l --recessive \t Flag that runs the analysis with a Recessive Model"
    " with (1*AA, 0*AB, 0*BB) effects. \n\t"
"-z --mylinear \t <path/filename> to read Factors 'f_i' for a Custom"
    " Linear Model with f1*X1, f2*X2, f3*X3...fn*X_ngpred as effects,\n\t"
"              \t each column of each independent variable will be"
    " multiplied with the specified factors. \n\t"
"              \t Formula: y~alpha*cov + beta_1*f1*X1 + beta_2*f2*X2"
    " +...+ beta_n*fn*Xn, (see example files!). \n\t"
"-y --myaddit  \t <path/filename> to read Factors 'f_i' for a Custom"
    " Additive Model with (f1*X1, f2*X2, f3*X3...fn*X_ngpred) as effects,\n\t"
"              \t each column of each independent variable will be"
    " multiplied with the specified factors and then added together. \n\t"
"              \t Formula: y~alpha*cov + beta*(f1*X1 + f2*X2 +...+"
    " fn*Xn), (see example files!). \n\t"
"-v --simpleinter <path/filename> to read the interactions from; for"
    " single analysis using multiple interactions. \n\t"
"-w --multinter \t <path/filename> to read the interactions from; for"
    " multiple analysis using single interaction per analysis. \n\t"
"-u --keepinter \t Flag that sets if the interaction analysis chosen is"
    " to too keep the dependent variable X.\n\t"
"            \t If set, Formula: y~alpha*cov + beta_1*INT*X + beta_2*X,"
    " (see example files!). \n\t"
"            \t Default not set, Formula: y~alpha*cov + beta_1*INT*X,"
    " (see example files!). \n\t\n\t"
"            \t Support for MPI is available. Simply use mpirun -np"
    " <#nodes> omicabelnomm <params> on an Open-MPI enabled computer/cluster.\n\t"
"            \t Recommended is to use MPI when dealing with problems"
    " with over 2000 genotypes, at a rate of 1 node per 2000 genotypes.\n";


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


void parse_params(int argc, char *argv[], struct Settings &params)
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
            {"ngpred", required_argument, 0, 'n'}, //r
            {"thr",    required_argument, 0, 't'},
            {"mem",    required_argument, 0, 'm'},
            {"excl",    required_argument, 0, 'x'},
            {"stobin",    no_argument, 0, 'b'},
            {"pdisp",    required_argument, 0, 'd'}, //
            {"psto",    required_argument, 0, 's'}, //
            {"rdisp",    required_argument, 0, 'r'}, //
            {"rsto",    required_argument, 0, 'e'}, //
            {"fdcov",    no_argument, 0, 'i'}, //
            {"fdgen",    no_argument, 0, 'f'}, //
            {"additive",    no_argument, 0, 'j'}, //
            {"dominant",   no_argument, 0, 'k'}, //
            {"recessive",    no_argument, 0, 'l'}, //
            {"mylinear",    required_argument, 0, 'z'}, //
            {"myaddit",    required_argument, 0, 'y'}, //
            {"simpleinter",    required_argument, 0, 'v'}, //
            {"multinter",    required_argument, 0, 'w'}, //
            {"keepinter",    no_argument, 0, 'u'}, //
            {"help",    no_argument, 0, 'h'}, //
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here.
        int option_index = 0;

        c = getopt_long(argc, argv,
                        "c:p:g:o:n:t:m:x:d:s:r:e:z:y:v:w:ufibhjkl",
                        long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        if (!optarg && (c != 'f' && c != 'b' && c != 'i' && c != 'h'))
        {
            if (params.mpi_id == 0) {
                cerr << "\nerror with argument parameter " <<
                    static_cast<char>(c) << endl;
            }
            exit(1);
        }

        switch (c)
        {
        case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
            {
                break;
            }
            printf("Long option %s", long_options[option_index].name);
            if (optarg)
            {
                printf(" with arg %s", optarg);
            }
            printf("\n");
            break;

        case 'p':
            phe = !phe;
            params.fnameY = string(optarg);

            pos = string(optarg).find_last_of(".");
            if (pos != string::npos)
            {
                params.fnameY = string(optarg).substr(0, pos);
            }
            if (params.mpi_id == 0)
                cout << "-p Reading phenotypes from file " << optarg << endl;
            break;

        case 'g':
            snp = !snp;
            params.fnameAR = string(optarg);

            pos = string(optarg).find_last_of(".");
            if (pos != string::npos)
                params.fnameAR = string(optarg).substr(0, pos);
            if (params.mpi_id == 0)
            {
                cout << "-g Reading with genotype data from file "
                     << optarg << endl;
            }
            break;

         case 'n':
            params.r = atoi(optarg);
            params.r = max(params.r, 1);
            if (params.mpi_id == 0)
                cout << "-n Using columns per snp as " << params.r << endl;
            break;

        case 'c':
            cov = !cov;
            params.fnameAL = string(optarg);

            pos = string(optarg).find_last_of(".");
            if (pos != string::npos)
                params.fnameAL = string(optarg).substr(0, pos);
            if (params.mpi_id == 0)
                cout << "-c Reading covariates from file " << optarg << endl;
            break;

        case 'o':
            bout = !bout;
            params.fnameOutFiles = string(optarg);

            pos = string(optarg).find_last_of(".");
            if (pos != string::npos)
                params.fnameOutFiles = string(optarg).substr(0, pos);
            if (params.mpi_id == 0)
                cout << "-o Writing output files to " << optarg << endl;
            break;

        case 't':
            params.threads = atoi(optarg);
            params.threads = min(max(params.threads, 1), 32);
            if (params.mpi_id == 0) {
                cout << "-t Using " << params.threads
                     << " CPU threads for parallelism" << endl;
            }
            break;

         case 'x':
            params.fname_excludelist = string(optarg);
            if (params.mpi_id == 0) {
                cout << "-x Excluding ids on "
                     << params.fname_excludelist
                     << endl;
            }
            break;


         case 'd':
            params.minPdisp = fabs(atof(optarg));
            if (params.mpi_id == 0)
            {
                cout << "-d Significance to display in .txt "
                     << "will be P-val's below " << params.minPdisp << endl;
            }
            break;

        case 's':
            params.minPstore = fabs(atof(optarg));
            if (params.mpi_id == 0)
            {
                cout << "-s Significance to store in .bin will be "
                     << "P-val's below " << params.minPstore << endl;
            }
            break;

        case 'r':
            params.minR2disp = (atof(optarg));
            if (params.mpi_id == 0)
                {
                    cout << "-r Minimum R2 to display in .txt will be above "
                         << params.minR2disp << endl;
                }
            break;

        case 'e':
            params.minR2store =  (atof(optarg));
            if (params.mpi_id == 0) {
                cout << "-e Minimum R2 to store in .bin will be above "
                     << params.minR2store << endl;
            }
            break;

        case 'i':
            params.disp_cov = true;
            if (params.mpi_id == 0) {
                cout << "-i Covariate results (significant) will be "
                     << "included in results whenever their respective "
                     << "snps are also significant" << endl;
            }
            break;

        case 'f':
            params.storePInd = true;
            if (params.mpi_id == 0) {
                cout << "-f Forcing all included results to be considered"
                     << " independently of max P-val or min R2. (SLOW!)"
                     << endl;
            }
            break;

        case 'j':
            params.model = 0;
            params.dosages = true;
            if (params.mpi_id == 0) {
                cout << "-j Using Additive Model with (2*AA, 1*AB, 0*BB)"
                     << " effects" << endl;
            }
            break;

        case 'k':
            params.model = 1;
            params.dosages = true;
            if (params.mpi_id == 0) {
                cout << "-k Using Dominant Model with (1*AA, 1*AB, 0*BB)"
                     << " effects" << endl;
            }
            break;

        case 'l':
            params.model = 2;
            params.dosages = true;
            if (params.mpi_id == 0) {
                cout << "-j Using Recessive Model with (1*AA, 0*AB, 0*BB)"
                     << " effects" << endl;
            }
            break;

        case 'z':
            params.model = 3;
            params.dosages = true;
            if (params.mpi_id == 0) {
                cout << "-z Using Custom Linear Model with parameters "
                     << "read from the file " << params.fname_dosages << endl;
            }
            break;

        case 'y':
            params.model = 4;
            params.dosages = true;
            if (params.mpi_id == 0) {
                cout << "-z Using Custom Additive Model with parameters"
                     << " read from the file " << params.fname_dosages
                     << endl;
            }
            break;

        case 'b':
            params.storeBin = true;

            if (params.mpi_id == 0) {
                cout << "-b Results will be stored in binary format too"
                     << endl;
            }
            break;

        case 'v':
            params.fname_interactions =  (atof(optarg));
            params.use_interactions = true;
            params.keep_depVar = true;
            params.use_multiple_interaction_sets = false;
            if (params.mpi_id == 0) {
                cout << "-v File containing single interactions "
                     << params.fname_interactions << endl;
            }
            break;

        case 'w':
            params.fname_interactions =  (atof(optarg));
            params.use_interactions = true;
            params.keep_depVar = true;
            params.use_multiple_interaction_sets = true;
            if (params.mpi_id == 0) {
                cout << "-w File containing multiple interactions "
                     << params.fname_interactions << endl;
            }
            break;

        case 'u':
            params.keep_depVar = true;
            if (params.mpi_id == 0) {
                cout << "-u Keeping independent variable for interaction"
                     << " analysis " << endl;
            }
            break;

        case 'h':
            cout << helpcmd << endl;
            cout << helpcmd_expl << endl;
            help_request = true;
            exit(1);
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
        {
            printf("%s ", argv[optind++]);
        }

        putchar('\n');
    }

    //required arguments
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
    Algorithm alg;

    //!default params
    alg.applyDefaultParams(params);

    //!MPI
    #ifdef USE_MPI
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //printf("SIZE = %d RANK = %d\n", size,rank);
    params.mpi_id = rank;
    params.mpi_num_threads = size;
    #endif



    parse_params(argc, argv, params);


    //omp_set_nested(false);
    //mkl_set_dynamic(false);
    //omp_set_dynamic(false);

    omp_set_num_threads(params.threads);
    blas_set_num_threads(params.threads);
    if (params.mpi_id == 0)
    {
        cout << PACKAGE_NAME " v" PACKAGE_VERSION "\n";
    }



    if (!help_request)
    {
        struct Outputs out = {0};
        if (params.use_fake_files)
        {
            alg.solve(params, out, P_NEQ_B_OPT_MD);
        }
        else
        {
            alg.solve(params, out, P_NEQ_B_OPT_MD);
        }
        cout << endl << "Number of Significant Results: "
             << out.total_sig_results <<  endl;
    }

    cout << endl;

    //!MPI
    #ifdef USE_MPI
    MPI_Finalize();
    #endif

    return 0;
}
