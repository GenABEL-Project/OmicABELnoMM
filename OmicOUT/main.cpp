#include <boost/math/special_functions/erf.hpp>

#include <fstream>
#include <stdint.h>
#include <unistd.h>
#include <limits.h>
#include <queue>
#include <list>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cstring>
#include <string>
#include <math.h>
#include <getopt.h>


using namespace std;

struct result
{
    float T;
    long double P;
    float SE;
    float B;
    int AL_name_idx;
    int AR_name_idx;
};

class resultH
{
    public:
    int Y_name_idx;
    int nUsed;
    float nUsedPct;
    float R2;
    list< result > res_p;
};

void float32(float* out, const uint16_t in)
{
    uint32_t t1 = in;
    t1 <<= 16;//convert to float again
    t1 &= 0xffff0000;//keep sign and exponent only
    *((uint32_t*)out) =  t1;
}

void statistics(float b,  float SE, int AL_name_idx, int AR_name_idx, resultH &res, double min_p_disp)
{
        long double ptemp;

        long double t;
        long double one_oversqrt2 = 0.7071067811865475244008443621048490392848359376884740;
        float t1;

        t=fabs(b/SE);

        ptemp=erfc(t*one_oversqrt2);

        }

        result rdata;

        rdata.B = b;
        rdata.T = t;
        rdata.P = ptemp;
        rdata.SE = SE;

        rdata.AR_name_idx = AR_name_idx;
        rdata.AL_name_idx = AL_name_idx;

        if((double)ptemp < min_p_disp)
        {
            res.res_p.push_back(rdata);

        }

}


void parse_params(int argc, char *argv[], string &source, string &dest, bool &forceAllRes, bool &includeCov, float &minP, float &minR2 )
{
    int c;

    unsigned pos;

    bool out = false;
    bool bin = false;
    bool pval = false;



    while (1)
    {
        static struct option long_options[] =
        {   //-c XL --snp XR -p Y -b B -r 2 -t 2
            // These options set a flag.
            //{"", no_argument,       &verbose_flag, 1},
            //{"",   no_argument,       &verbose_flag, 0},
            // These options don't set a flag.
            //We distinguish them by their indices.
            {"binsource",    required_argument, 0, 'b'},
            {"out",    required_argument, 0, 'o'},

            {"pdisp",    required_argument, 0, 'd'},//
            {"rdisp",    required_argument, 0, 'r'},//
            {"fdcov",    no_argument, 0, 'i'},//
            {"fdgen",    no_argument, 0, 'f'},//
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here.
        int option_index = 0;

        c = getopt_long(argc, argv, "b:o:d:r:fi", long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        if (!optarg && (c != 'f' && c != 'i'))
        {
            cout << "\nerror with argument parameter " << (char)c << endl;
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


        case 'o':
            out = !out;
            dest = string(optarg);

//            dest = string(optarg).find(".");
//            if (pos != string::npos)
//                dest = string(optarg).substr(0, pos);

            cout << "-o Writing output text file to " << dest << endl;
            break;

        case 'b':
            bin = !bin;
            source = string(optarg);

            pos = string(optarg).find(".");
            if (pos != string::npos)
                source = string(optarg).substr(0, pos);

            cout << "-o Reading input from file " << source << endl;
            break;


         case 'd':
            minP = fabs(atof(optarg));
            pval = true;

            cout << "-d Significance to display in .txt will be P-val's below " << minP << endl;
            break;

        case 'r':
            minR2 = (atof(optarg));

            cout << "-r Minimum R2 to display in .txt will be above " << minR2  << endl;
            break;

        case 'i':
            includeCov = true;

            cout << "-i Covariate results will be included in results" << endl;
            break;

        case 'f':
            forceAllRes = true;

            cout << "-f Forcing all includec results to be considered independently of max P-val or min R2. (Overrides -t and -d, Big file!)"<< endl;
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

    //required arguments
    if (!bin || !out || !pval)
    {
        cout << "Missing arguments!:"
             << " o:"<< out << " b:" << bin << " d:" << pval
             << endl;
        exit(1);
    }
}

int main(int argc, char *argv[] )
{
        string source;
        string dest;
        bool forceAllRes = false;
        bool includeCov = false;
        float minP;
        float minR2 = -100.0;

        parse_params(argc,argv,source,dest,forceAllRes,includeCov,minP,minR2);

        ifstream fp_allResults;
        fp_allResults.open((source + ".dbin").c_str(),ios::in | ios::binary );
        if(fp_allResults == 0)
        {
            cout << "Error reading File "<< (source + ".dbin") << endl;
            exit(1);
        }
        ifstream fp_InfoResults;
        fp_InfoResults.open((source + ".ibin").c_str(),ios::in | ios::binary );
        if(fp_InfoResults == 0)
        {
            cout << "Error Creating File " << (source + ".ibin") << endl;
            exit(1);
        }
        ofstream fp_sigResults;
        fp_sigResults.open((dest + ".txt").c_str(),ios::out | ios::trunc);
        if(fp_sigResults == 0)
        {
            cout << "Error Creating File "<< (dest + ".txt") << endl;
            exit(1);
        }


        int n,l,r,p,t,m,namelength;
        bool disp_cov, storePInd;

        bool new_disp_cov = includeCov;

        double min_p_disp = minP;

        if(forceAllRes)
        {
            min_p_disp = 2.0;
            minR2 = -100.0;
        }

        bool only_y_selected = false;
        list < string > y_selected;//TODO a list of selceted Y to display, onyl these values are meant to be displayed

        char* ALnames;
        char* ARnames;
        char* Ynames;

        //info
        {

            fp_InfoResults.read( (char*)&n,sizeof(int));
            fp_InfoResults.read( (char*)&l,sizeof(int));
            fp_InfoResults.read( (char*)&r,sizeof(int));
            fp_InfoResults.read( (char*)&t,sizeof(int));
            fp_InfoResults.read( (char*)&m,sizeof(int));
            fp_InfoResults.read( (char*)&namelength,sizeof(int));
            fp_InfoResults.read( (char*)&disp_cov,sizeof(bool));
            fp_InfoResults.read( (char*)&storePInd,sizeof(bool));

            new_disp_cov = new_disp_cov && disp_cov;

            if(!(n >0 && l > 1 && r > 0 && m > 0 && t > 0 && namelength > 0))
            {
                cout << "Input data had invalid info data!" << endl;
                exit(1);
            }

            if(l>1)
            {
                ALnames = new char[namelength*(l-1)];
                fp_InfoResults.read( ALnames,namelength*(l-1));
            }

            ARnames = new char[namelength*(r*m)];
            Ynames = new char[namelength*(t)];



            fp_InfoResults.read( ARnames,namelength*r*m);
            fp_InfoResults.read( Ynames,namelength*t);

        }

        p=l+r;



        //if(storePInd)
        {
            int h_max = r;
            if(disp_cov)
                h_max+=l-1;

            uint16_t n_Used;
            float R2;
            uint16_t R2loaded;

            float B;
            float SE;


            uint8_t AL_name_idx;
            uint16_t AR_name_idx;
            int Y_name_idx;

            int size_block = m;

            fp_allResults.seekg (0, fp_allResults.end);
            unsigned long int max_to_read = fp_allResults.tellg();
            fp_allResults.seekg (0, fp_allResults.beg);



            unsigned long int read=0;

            for(int j=0; j < t && read < max_to_read; j++)
            {
                Y_name_idx = j;

                int ARoffset=0;

                if(!storePInd)
                {
                     fp_allResults.read((char*)&Y_name_idx,sizeof(int));
                     fp_allResults.read((char*)&size_block,sizeof(int));
                     fp_allResults.read((char*)&ARoffset,sizeof(int));
                     read+=sizeof(int)*3;
                }

                if(only_y_selected)
                {
                    list < string >::iterator iter = find( y_selected.begin(),y_selected.end() , string(&Ynames[Y_name_idx]));
                    if( iter == y_selected.end() )
                    {
                        size_block = 0;//this means that the entire list of X for this y wont be even considered int he for loop below
                    }
                }

                for(int i=0; i < size_block; i++)
                {
                    fp_allResults.read((char*)&n_Used,sizeof(uint16_t));
                    fp_allResults.read((char*)&R2loaded,sizeof(uint16_t));
                    read+=sizeof(uint16_t)*2;

                    float32(&R2, R2loaded);

                    uint8_t size_p = (uint8_t)h_max;

                    if(!storePInd)
                    {
                        fp_allResults.read((char*)&size_p,sizeof(uint8_t));
                        read+=sizeof(uint8_t);
                    }

                    resultH res;
                    res.R2 = R2;
                    res.nUsed = n_Used;
                    res.Y_name_idx = Y_name_idx;
                    res.nUsedPct = (float)n_Used/(float)n;
                    int p_max = (int) size_p;


                    for(int h=0; h < p_max; h++)//intercept was not stored so size_p will be at most p-1
                    {

                        if(disp_cov && !storePInd)
                        {
                            fp_allResults.read((char*) &AL_name_idx,sizeof(uint8_t));
                            read+=sizeof(uint8_t);
                        }
                        if(!storePInd)
                        {
                            fp_allResults.read((char*) &AR_name_idx,sizeof(uint16_t));
                            read+=sizeof(uint16_t);
                            AR_name_idx += ARoffset;
                        }
                        else
                        {
                            AR_name_idx = i*r+h;
                            if(disp_cov)
                            {
                                AR_name_idx-=l-1;
                                if(h < (l-1))
                                {
                                    AR_name_idx = i*r;
                                    AL_name_idx = h;
                                }
                                else
                                {
                                    AL_name_idx = -1;//255 in ubyte form
                                }

                            }


                        }

                        fp_allResults.read((char*)&B,sizeof(float));
                        fp_allResults.read((char*)&SE,sizeof(float));

                        read+=sizeof(float)*2;

                        bool skip_cov = (!new_disp_cov && AL_name_idx != 255);

                        if( !disp_cov || !skip_cov )//-1;//255 in ubyte form
                            statistics(B,SE,AL_name_idx,AR_name_idx,res, min_p_disp);
                    }


                    for (list<  result  >::iterator it2=res.res_p.begin(); it2 != res.res_p.end(); ++it2)
                    {
                        if(res.R2 > minR2)
                        {
                            result res_p = *it2;

                            string Aname = " ";

                            if(new_disp_cov && disp_cov && res_p.AL_name_idx != 255 )
                            {
                                Aname = string(&ALnames[namelength*(res_p.AL_name_idx)]) + ":";
                            }

                            Aname += string(&ARnames[namelength*res_p.AR_name_idx]);


                            fp_sigResults << string(&Ynames[namelength*res.Y_name_idx]) << "\t" << Aname << "\t"
                                    << std::fixed << std::setprecision(2) << res.nUsed << "\t" << res.nUsedPct << "\t"
                                    << std::setprecision(std :: numeric_limits<float> :: digits10 + 1) << res_p.B  << "\t" << std::fixed << std::setprecision(4) << res.R2
                                          << "\t" << std::setprecision(std :: numeric_limits<float> :: digits10 + 1) << res_p.SE   << "\t" << res_p.T
                                                << "\t" << std::scientific  << res_p.P << std::fixed <<  std::setprecision(std :: numeric_limits<float> :: digits10 + 1) << endl;
                        }
                    }




                }
            }
        }












}
