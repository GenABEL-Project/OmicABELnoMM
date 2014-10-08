#include "AIOwrapper.h"



AIOwrapper::AIOwrapper()
{
    Fhandler = &FHandler;
    io_overhead = "*";
}

AIOwrapper::~AIOwrapper()
{

}


void AIOwrapper::initialize(struct Settings &params)
{


    pthread_mutex_init( &(FHandler.m_more), NULL);
    pthread_mutex_init( &(FHandler.m_read), NULL);
    pthread_mutex_init( &(FHandler.m_buff_upd), NULL);
    pthread_mutex_init( &(FHandler.out_buff_upd), NULL);

    pthread_cond_init(&(FHandler.condition_more), NULL);
    pthread_cond_init(&(FHandler.condition_read), NULL);

    pthread_barrier_init(&(FHandler.finalize_barrier),NULL,2);


    Fhandler->fakefiles = params.use_fake_files;


    Fhandler->use_dosages = params.dosages;
    if(params.dosages && params.model ==-1)
    {
        cout << "Requested dosages model wihtout a valid model!" << endl;
        exit(1);
    }

    Fhandler->not_done = true;
    Fhandler->model = params.model;
    Fhandler->fname_dosages = params.fname_dosages;
    Fhandler->dosage_skip = 0;
    Fhandler->fileR = 0;
    Fhandler->add_dosages = false;



    if(!Fhandler->fakefiles)
    {
        Fhandler->fnameAL = params.fnameAL;
        Fhandler->fnameAR = params.fnameAR;
        Fhandler->fnameY = params.fnameY;
        Fhandler->fnameOutFiles = params.fnameOutFiles;

        Fhandler->storeBin = params.storeBin;

        Fhandler->disp_cov = params.disp_cov;
        Fhandler->storePInd = params.storePInd;

        Fhandler->min_p_disp = params.minPdisp;
        Fhandler->min_R2_disp = params.minR2disp;


        Yfvi  = load_databel_fvi( (Fhandler->fnameY+".fvi").c_str() );
        ALfvi = load_databel_fvi( (Fhandler->fnameAL+".fvi").c_str() );
        ARfvi = load_databel_fvi( (Fhandler->fnameAR+".fvi").c_str() );

        if(ALfvi->fvi_header.numObservations !=  ARfvi->fvi_header.numObservations &&  ARfvi->fvi_header.numObservations != Yfvi->fvi_header.numObservations
              &&  ARfvi->fvi_header.numObservations != Yfvi->fvi_header.numObservations   )
        {
            cout << "The number of elements/individuals from the input files do not coincide with each other! See:" << endl;
            cout << Fhandler->fnameY << ":" <<  Yfvi->fvi_header.numObservations << endl;
            cout << Fhandler->fnameAL << ":" << ALfvi->fvi_header.numObservations  << endl;
            cout << Fhandler->fnameAR << ":" << ARfvi->fvi_header.numObservations << endl;
            exit(1);
        }

        params.n = min((int)(ALfvi->fvi_header.numObservations),params.limit_n);
        Fhandler->fileN = params.n;
        realN = params.n;


        params.m = min((int)(ARfvi->fvi_header.numVariables/params.r),params.limit_m/params.r);

        #ifdef DEBUG
        cout << "calcM:" << params.m << endl;
        #endif

        if((int)(ARfvi->fvi_header.numVariables) % params.r != 0)
        {
            cout << "Number of elements in "<< (Fhandler->fnameAR+".fvi") <<":"<<  (int)(ARfvi->fvi_header.numVariables)
                        << " is not consistent with the value of --ngpred (-n)" << params.r << " provided" << endl;
            exit(1);
        }

        Fhandler->fileM = params.m;


        Fhandler->fileR = params.r;
        Fhandler->modelR = Fhandler->fileR;

        params.t = min((int)(Yfvi->fvi_header.numVariables),params.limit_t);
        params.l = ALfvi->fvi_header.numVariables;


        //!check id names order
        string  YidNames;
        string  ALidName;
        string  ARidName;
        string  INTidName;

        int yname_idx=0;//starting idx for names on ALfvi->data
        for(int i = 0; i < params.n; i++)
        {
            YidNames = string(&(Yfvi->fvi_data[yname_idx]));
            ALidName = string(&(ALfvi->fvi_data[yname_idx]));
            ARidName = string(&(ARfvi->fvi_data[yname_idx]));
            if(YidNames.compare(ALidName))
            {
                cout << "Warning, ID names between -p file and -c file do not coincide! The files must have the same meaningful order" << endl;
                i = params.n;//exit
            }
            if(YidNames.compare(ARidName))
            {
                cout << "Warning, ID names between -p file and -g file do not coincide! The files must have the same meaningful order" << endl;
                i = params.n;//exit
            }
            if(ARidName.compare(ALidName))
            {
                cout << "Warning, ID names between -g file and -c file do not coincide! The files must have the same meaningful order" << endl;
                i = params.n;//exit
            }
            yname_idx += Yfvi->fvi_header.namelength;


        }


        for(int i = 0; i < params.t; i++)
        {
            Fhandler->Ynames.push_back(string(&(Yfvi->fvi_data[yname_idx])));
            yname_idx += Yfvi->fvi_header.namelength;
        }


        int Aname_idx=params.n*ALfvi->fvi_header.namelength;
        for(int i = 0; i < params.l; i++)
        {
            Fhandler->ALnames.push_back(string(&(ALfvi->fvi_data[Aname_idx])));
            Aname_idx += ALfvi->fvi_header.namelength;

//            cout << i << ": " << string(&(ARfvi->fvi_data[Aname_idx])) << "\t";
//            cout << i << ": " << ALnames[i] << "\t";
        }



        int opt_tb = 1000;
        int opt_mb = 500;

        params.mb = min(params.m, opt_mb);
        params.tb = min(params.t, opt_tb);





    }
    else
    {
        //other params come from outside
    }

    params.mb = min(params.m,params.mb);
    params.tb = min(params.t,params.tb);

    Fhandler->Ar_file_blocksize = params.mb;


    //!excludeLIST
    int excl_count = 0;
    int Almissings = 0;
    Fhandler->excl_List = new list< pair<int,int> >();






    (Fhandler->excl_List)->push_back( make_pair(0,params.n) );

    if(params.fname_excludelist.size()!=0)
    {
        read_excludeList( Fhandler->excl_List,excl_count,params.n,params.fname_excludelist);
//         for (list<  pair<int,int>  >::iterator it= (Fhandler->excl_List)->begin(); it !=  (Fhandler->excl_List)->end(); ++it)
//         {
//            cout << it->first << "-" << it->second << " | ";
//         }
    }

    if(!Fhandler->fakefiles)
    {

        removeALmissings(Fhandler->excl_List,params,Almissings);

    }


    //!DOSAGES
    if(params.dosages)
    {

        Fhandler->ArDosage = new float[Fhandler->fileR*params.n];
        Fhandler->dosages = new float[Fhandler->fileR];
        Fhandler->modelR = 1;


        switch (Fhandler->model)
        {
        case -1://nomodel

        break;
        case 0://add
            if(Fhandler->fileR != 3 && Fhandler->fileR != 2)
            {
                cout << "The amount of columns per Independent Variable (--ngpred) is not 2 or 3 for a valid Additive Model!" << endl;
                exit(1);
            }
            Fhandler->dosages[0] = 2;Fhandler->dosages[1] = 1;//Fhandler->dosages[2] = 0;
            if(Fhandler->fileR == 3)
                Fhandler->dosage_skip = 1;//last columns is irrelevant
            params.r = 1;
            Fhandler->add_dosages = true;
        break;
        case 1://dom
            if(Fhandler->fileR != 3 && Fhandler->fileR != 2)
            {
                cout << "The amount of columns per Independent Variable (--ngpred) is not 2 or 3 for a valid Dominant Model!" << endl;
                exit(1);
            }
            Fhandler->dosages[0] = 1;Fhandler->dosages[1] = 1;//Fhandler->dosages[2] = 0;
            if(Fhandler->fileR == 3)
                Fhandler->dosage_skip = 1;//last columns is irrelevant
            params.r = 1;
            Fhandler->add_dosages = true;
        break;
        case 2://rec
            if(Fhandler->fileR != 3 && Fhandler->fileR != 2)
            {
                cout << "The amount of columns per Independent Variable (--ngpred) is not 2 or 3 for a valid Recessive Model!" << endl;
                exit(1);
            }
            Fhandler->dosages[0] = 1;//Fhandler->dosages[1] = 0;Fhandler->dosages[2] = 0;
            Fhandler->dosage_skip = 2;
            if(Fhandler->fileR == 3)
                (Fhandler->dosage_skip)--;//last 2 columns is irrelevant
            params.r = 1;
            Fhandler->add_dosages = true;
        break;
        case 3://linear

            read_dosages(params.fname_dosages,Fhandler->fileR,Fhandler->dosages);
            Fhandler->add_dosages = false;
            Fhandler->dosage_skip = 0 ;
            Fhandler->modelR = params.r;
        break;
        case 4://additive
            read_dosages(params.fname_dosages,Fhandler->fileR,Fhandler->dosages);
            params.r = 1;
            Fhandler->add_dosages = true;
        break;
        }
    }

    //!INTERACTIONS

    Fhandler->use_interactions = params.use_interactions;
    Fhandler->keep_depVar = params.keep_depVar;
    Fhandler->use_multiple_interaction_sets = params.use_multiple_interaction_sets;
    Fhandler->fname_interactions = params.fname_interactions;

    if(Fhandler->use_interactions && Fhandler->modelR > 1)
    {
        if(Fhandler->use_dosages)
        {
            cout << "Cannot use Interactions with the current dosage model! #Of factors can be at most 1!" << endl;
        }

        {
            cout << "Cannot use Interactions with --ngpred set to: " << params.r << ". Required is 1!" <<  endl;
            cout << "Interactions can be only of the form y~cov+V1*X+...+Vi*X, (only with --ngpred=1 or models excluding(--mylinear)) !" << endl;
        }
        exit(1);
    }


    int interaction_missings = 0;
    if( Fhandler->use_interactions)
    {

        //cout << Fhandler->fname_interactions;
        Ifvi  = load_databel_fvi( (Fhandler->fname_interactions+".fvi").c_str() );

        if(Ifvi->fvi_header.numObservations !=  ARfvi->fvi_header.numObservations  )
        {
            cout << "The number of elements/individuals from the input files do not coincide with each other! See:" << endl;
            cout << Fhandler->fnameY << ":" <<  Ifvi->fvi_header.numObservations << endl;
            cout << Fhandler->fnameAR << ":" << ARfvi->fvi_header.numObservations << endl;
            exit(1);
        }

        int name_idx= 0;
        for(int i = 0; i < params.n; i++)
        {
            string ARidName = string(&(ARfvi->fvi_data[name_idx]));
            string INTidName = string(&(Ifvi->fvi_data[name_idx]));
            if(ARidName.compare(INTidName))
            {
                cout << "Warning, ID names between -g file and interactions file do not coincide! The files must have the same meaningful order!" << endl;
                i = params.n;//exit
            }

            name_idx += Yfvi->fvi_header.namelength;
        }


        if(params.end_IDX_interactions != -1 || params.ini_IDX_interactions != -1)
        {
            Fhandler->numInter = params.end_IDX_interactions - params.ini_IDX_interactions+1;
            Fhandler->ini_IDX_interactions =  params.ini_IDX_interactions;
            Fhandler->end_IDX_interactions = params.end_IDX_interactions;
        }
        else
        {
            Fhandler->numInter = Ifvi->fvi_header.numVariables;
            Fhandler->ini_IDX_interactions = 0;
            Fhandler->end_IDX_interactions = Fhandler->numInter;
        }

        if((unsigned int)(Fhandler->numInter) > Ifvi->fvi_header.numVariables)
        {
            cout << "Amount of data requested to be read from the Interactions file exceeds its contents!" << endl;
            exit(1);
        }

        if((unsigned int)Fhandler->ini_IDX_interactions > Ifvi->fvi_header.numVariables || (unsigned int)Fhandler->end_IDX_interactions > Ifvi->fvi_header.numVariables)
        {
            cout << "The specified range for data to be read from the Interactions file exceeds its content dimensions!" << endl;
            exit(1);
        }

        if(Fhandler->use_multiple_interaction_sets)
        {
            cout << "rm:"<< params.m<<  " rmb:" << params.mb << endl;

            params.m *= Fhandler->numInter;
            Fhandler->Ar_file_blocksize = params.mb/Fhandler->numInter;

            if(Fhandler->Ar_file_blocksize < 1)
                Fhandler->Ar_file_blocksize = 1;

            params.mb = (Fhandler->Ar_file_blocksize*Fhandler->numInter);
            #ifdef DEBUG
            cout << " lmb:"<< params.mb << " fmb:" << Fhandler->Ar_file_blocksize << " lr:" << params.r << " fr:" << Fhandler->fileR << endl << flush;
            #endif
        }
        else
        {
            params.r = Fhandler->numInter;
            Fhandler->Ar_file_blocksize = params.mb;
            #ifdef DEBUG
            cout << " lmb:"<< params.mb << " fmb:" << Fhandler->Ar_file_blocksize << " lr:" << params.r << " fr:" << Fhandler->fileR << endl << flush;
            #endif
        }

        if(Fhandler->keep_depVar)
        {
            params.r += 1;
        }
        #ifdef DEBUG
        cout << Fhandler->Ar_file_blocksize << "!" << endl << flush;
        #endif




        ifstream fp_I;
        fp_I.open ((Fhandler->fname_interactions+".fvd").c_str(), ios::in | ios::binary);
        if(!fp_I.is_open())
        {
            cout << "Error opening Interactions File " << Fhandler->fname_interactions << endl;
            exit(1);
        }


        Fhandler->interaction_data = new float[Fhandler->numInter*Fhandler->fileN];

        int ini_filepos = Fhandler->ini_IDX_interactions*Fhandler->fileN;
        fp_I.seekg(  ini_filepos*sizeof(type_precision) , ios::beg );
        if(fp_I.fail())
        {
            cout << "Error reading Interactions File while positioning! " <<  ini_filepos << endl;
            exit(1);
        }

        //!read interactions
        list< pair<int,int> >* excl_List = Fhandler->excl_List;



        int chunk_size_buff = Fhandler->numInter*Fhandler->fileN;

        fp_I.read ((char*)(Fhandler->interaction_data),sizeof(type_precision)*chunk_size_buff);



        for (int i=0; i < Fhandler->numInter; i++)
        {
            for (int j=0; j < Fhandler->fileN; j++)
            {
                if(isnan(Fhandler->interaction_data[i*Fhandler->fileN+j]))
                {
                    //cout << "fmi"<<" ";
                    bool removed  = splitpair(j,excl_List,Fhandler->fileN );
                    if(removed)
                    {
                        interaction_missings++;
                    }

                }
            }
        }
        if(interaction_missings)
            cout << "Excluding " << interaction_missings << " from Interaction missings" << endl;
        #ifdef DEBUG
        cout << params.mb << " " << params.r << " " << params.m << endl;
        #endif

    }


//    for (list<  pair<int,int>  >::iterator it=Fhandler->excl_List->begin(); it != Fhandler->excl_List->end(); ++it)
//    {
//        cout << it->first <<":" << it->second <<" ";
//    }
    #ifdef DEBUG
    cout << Fhandler->excl_List->size() << endl;
    #endif



    params.n -= (excl_count + Almissings + interaction_missings);

    params.p = params.l + params.r;

    if(!Fhandler->fakefiles)
    {

        //!******nameGATHER****************

        int Aname_idx=Fhandler->fileN*ARfvi->fvi_header.namelength;//skip the names of the rows
        if(Fhandler->use_dosages && Fhandler->add_dosages)
        {
            for(int i = 0; i < Fhandler->fileM; i++)
            {
                if(Fhandler->use_interactions && !Fhandler->use_multiple_interaction_sets)
                {
                    for(int ii = 0; ii < Fhandler->numInter; ii++)
                    {

                        string interaction_name = "I";
                        interaction_name += std::to_string(ii+1);
                        interaction_name += "x";
                        interaction_name += string(&(ARfvi->fvi_data[Aname_idx]));
                        Fhandler->ARnames.push_back(interaction_name);
                    }
                }
                if( !Fhandler->use_interactions || Fhandler->keep_depVar ||  Fhandler->use_multiple_interaction_sets)
                {
                    Fhandler->ARnames.push_back(string(&(ARfvi->fvi_data[Aname_idx])));
                }

                Aname_idx += ARfvi->fvi_header.namelength*Fhandler->fileR;
            }
        }
        else
        {
            for(int i = 0; i < Fhandler->fileM*Fhandler->fileR; i++)
            {
                if(Fhandler->use_interactions && !Fhandler->use_multiple_interaction_sets)
                {
                    for(int ii = 0; ii < Fhandler->numInter; ii++)
                    {

                        string interaction_name = "I";
                        interaction_name += std::to_string(ii+1);
                        interaction_name += "x";
                        interaction_name += string(&(ARfvi->fvi_data[Aname_idx]));
                        Fhandler->ARnames.push_back(interaction_name);
                       // cout << interaction_name<< endl;
                    }
                }
                if( Fhandler->use_multiple_interaction_sets || !Fhandler->use_interactions || Fhandler->keep_depVar)
                {
                    Fhandler->ARnames.push_back(string(&(ARfvi->fvi_data[Aname_idx])));
                   // cout << string(&(ARfvi->fvi_data[Aname_idx])) << endl;
                }
                Aname_idx += ARfvi->fvi_header.namelength;
            }

        }

        if(Fhandler->use_interactions && Fhandler->use_multiple_interaction_sets)
        {
            int ar_name_restamount = Fhandler->fileM;
            vector<string> original_names = std::move(Fhandler->ARnames);

            //cout <<"ar_name_restamount:" <<ar_name_restamount <<" " << Fhandler->Ar_file_blocksize << flush;

            int ini_idx=0;
            while(ar_name_restamount > 0 )
            {

                int ar_name_blocksize = Fhandler->Ar_file_blocksize;
                if(ar_name_restamount < ar_name_blocksize)
                {
                    ar_name_blocksize = ar_name_restamount;
                }

                for(int ii = 0; ii < Fhandler->numInter; ii++)
                {
                    string interaction_name = "I";
                    interaction_name += std::to_string((ii+1));
                    interaction_name += "_";
                    int name_idx = ini_idx;
                    for(int i = 0; i < ar_name_blocksize; i++)
                    {
                        Fhandler->ARnames.push_back(interaction_name+original_names[name_idx]);
                       // cout <<interaction_name+original_names[name_idx] << endl;
                        if( Fhandler->keep_depVar)
                        {
                            string interaction_name_orig= interaction_name+"o_";
                            Fhandler->ARnames.push_back(interaction_name_orig+original_names[name_idx]);
                           // cout << interaction_name_orig+original_names[name_idx] << endl;
                        }
                        name_idx++;
                    }
                }

                ini_idx += ar_name_blocksize;

                ar_name_restamount -= ar_name_blocksize;
            }
        }

        //cout << "ARnames:" << Fhandler->ARnames.size() << endl;


        //!******endofnameGATHER****************

        if(params.l > 255 || params.p > 255 || params.n > 65535)//can remove if fixed for output files
        {
            cout << "Warning, output binary format does not yet support current problem sizes for the provided p, l, r and n." << endl;
            cout << "Omitting outputfile." << endl;
        }

        if(params.storeBin)
        {
            //!write info file for results
            ofstream fp_InfoResults;

            fp_InfoResults.open((Fhandler->fnameOutFiles + "_sto.ibin").c_str(),ios::out | ios::binary | ios::trunc);
            if(fp_InfoResults == 0)
            {
                cout << "Error Creating File: "  << (Fhandler->fnameOutFiles + "_sto.ibin") << endl;
                exit(1);
            }
            fp_InfoResults.write( (char*)&params.n,sizeof(int));
            fp_InfoResults.write( (char*)&params.l,sizeof(int));
            fp_InfoResults.write( (char*)&params.r,sizeof(int));
            fp_InfoResults.write( (char*)&params.t,sizeof(int));
            fp_InfoResults.write( (char*)&params.m,sizeof(int));
            fp_InfoResults.write( (char*)&(ARfvi->fvi_header.namelength),sizeof(int));
            fp_InfoResults.write( (char*)&params.disp_cov,sizeof(bool));
            fp_InfoResults.write( (char*)&params.storePInd,sizeof(bool));


            Aname_idx=(Fhandler->fileN+1)*ALfvi->fvi_header.namelength;//skip the names of the rows + intercept(+1)
            fp_InfoResults.write( (char*)&ALfvi->fvi_data[Aname_idx],ALfvi->fvi_header.namelength*(params.l-1)*sizeof(char));

            Aname_idx=Fhandler->fileN*ARfvi->fvi_header.namelength;//skip the names of the rows
            if(Fhandler->use_dosages && Fhandler->add_dosages)
            {
                for(int i = 0; i < params.m; i++)
                {
                    fp_InfoResults.write( (char*)&ARfvi->fvi_data[Aname_idx],ARfvi->fvi_header.namelength*sizeof(char));
                    Aname_idx += Fhandler->fileR*ARfvi->fvi_header.namelength*sizeof(char);
                }
            }
            else
            {
                fp_InfoResults.write( (char*)&ARfvi->fvi_data[Aname_idx],Fhandler->fileR*params.m*ARfvi->fvi_header.namelength*sizeof(char));
            }

            int Yname_idx=Fhandler->fileN*Yfvi->fvi_header.namelength;//skip the names of the rows
            fp_InfoResults.write( (char*)&Yfvi->fvi_data[Yname_idx],Yfvi->fvi_header.namelength*params.t*sizeof(char));

        }


    }


    //check p < n!!!!
    if(params.p >= params.n)
     {
        cout << "Error! The amount of unknows (individuals|elements) must be greater than total of variables (#covariates+#ind.variables)!" << endl;
        exit(1);
    }

    //block size to keep mem under 1 gigabyte
//    int opt_block = params.n/(4*1000^3)*(1/(2*params.r));
//    int opt_tb = max(4*2000,opt_block);
//    int opt_mb = max(2000,opt_block);
//
    #ifdef DEBUG
    cout << "m:" << params.m << " lmb:" << params.mb << " ln:"<<  params.n <<  " lr:" << params.r << endl ;
    #endif

    prepare_AL(params.l,params.n);
    prepare_AR(  params.mb,  params.n,  Fhandler->fileM,  params.r);
    prepare_OutFiles();
    prepare_Y(params.tb, params.n, params.t);



    //

}

void AIOwrapper::generate_singleinteraction_multset(float* interaction_data, int cols_indp_data, float* ARorig, int ar_orig_block_size, int n, float* dest, bool keep_original )
{
    int new_r = 1;
    if(keep_original)
    {
        new_r += 1;
    }

    //cout << ar_orig_block_size << ":" << new_r << ":" << cols_indp_data << ":" << n << endl;
    for(int h = 0; h < cols_indp_data; h++)
    {
        for(int i = 0; i < ar_orig_block_size; i++)
        {
            for(int k = 0; k < n; k++)
            {
                float ak = ARorig[i*n+k];
                float idta = interaction_data[h*n+k];

                dest[h*ar_orig_block_size*n*new_r+i*n*new_r+k] = idta*ak;

            }
            if(keep_original)
            {
                for(int k = 0; k < n; k++)
                {
                    dest[h*n*ar_orig_block_size*new_r+i*new_r*n+1*n+k] = ARorig[i*n+k];
                }
            }
        }
    }
}

void AIOwrapper::generate_multinteraction_singleset(float* interaction_data, int cols_data, float* AR, int ar_block_size, int n, float* dest, bool keep_original  )
{//colsdata =1  and dest


    int new_r = cols_data;
    if(keep_original)
    {
        new_r = cols_data+1;
    }

    //cout << ar_block_size << ":" << new_r << ":" << cols_data << ":" << n << endl;
    for(int i = 0; i < ar_block_size; i++)
    {
        for(int h = 0; h < cols_data; h++)
        {
            for(int k = 0; k < n; k++)
            {
                dest[i*new_r*n+h*n+k] = interaction_data[h*n+k]*AR[i*n+k];
//                if(h==2)
//                    dest[i*new_r*n+h*n+k] = interaction_data[h*n+k];
            }
            //cout << dest[i*new_r*n+h*n] << ":" << AR[i*n] << " ";
        }
        if(keep_original)
        {
            for(int k = 0; k < n; k++)
            {
                dest[i*(new_r)*n+cols_data*n+k] = AR[i*n+k];
            }
            //cout << dest[i*(new_r)*n+cols_data*n] << ":0";
        }
    }
}

void AIOwrapper::finalize()
{

    Fhandler->not_done = false;
    pthread_mutex_lock(&(Fhandler->m_more));
    pthread_cond_signal( &(Fhandler->condition_more ));
    pthread_mutex_unlock(&(Fhandler->m_more));

    pthread_barrier_wait(&(Fhandler->finalize_barrier));

    pthread_join( (Fhandler->iothread), NULL);

    finalize_Y();
    finalize_AR();
    finalize_AL();
    finalize_OutFiles();

    pthread_attr_destroy(&(Fhandler->attr));

    pthread_mutex_destroy(&(Fhandler->m_more));
    pthread_cond_destroy(&(Fhandler->condition_more));

    pthread_mutex_destroy(&(Fhandler->m_read));
    pthread_cond_destroy(&(Fhandler->condition_read));

    delete Fhandler->excl_List;
    if(Fhandler->use_dosages)
    {
            delete [](Fhandler->ArDosage);
            delete [](Fhandler->dosages);
    }






}




void AIOwrapper::finalize_OutFiles()
{

}


void* AIOwrapper::async_io( void *ptr )
{
    //cout << "async_io\n" << flush;
    type_fileh* Fhandler = (type_fileh *)ptr;
    int size_buff,tmp_y_blockSize,tmp_ar_blockSize;

    struct timespec timeToWait;
    FILE*  fp_Y;
    ifstream fp_Ar;
    ifstream fp_I;

    ofstream fp_sigResults;
    ofstream fp_allResults;

    int y_file_pos = 0;
    int ar_file_pos = 0;


    int seq_count;
    int max_secuential_write_count= 10;


    if(!Fhandler->fakefiles)
    {
        fp_Y = fopen((Fhandler->fnameY+".fvd").c_str(), "rb");
        if(fp_Y == 0)
        {
            cout << "Error opening File Y " << Fhandler->fnameY << endl;
            exit(1);
        }
        fp_Ar.open ((Fhandler->fnameAR+".fvd").c_str(), ios::in | ios::binary);

        if(!fp_Ar.is_open())
        {
            cout << "Error opening File Xr " << Fhandler->fnameAR << endl;
            exit(1);
        }


        fp_sigResults.open((Fhandler->fnameOutFiles + "_dis.txt").c_str(),ios::out | ios::trunc);
        if(fp_sigResults == 0)
        {
            cout << "Error Creating File " << (Fhandler->fnameOutFiles + "_dis.txt") << endl;
            exit(1);
        }

        fp_sigResults << "Phe\tSNP\tn\tnPct\tB\tR2\tSE\tT\tP" << endl;

        //float *a= new float [Fhandler->fileM*Fhandler->fileM*Fhandler->fileM];

        if(Fhandler->use_interactions)
        {
            //!************************************************
            int ini_filepos = Fhandler->ini_IDX_interactions*Fhandler->fileN;
            int chunk_size_buff = 0;


            fp_I.open ((Fhandler->fname_interactions+".fvd").c_str(), ios::in | ios::binary);
            if(!fp_I.is_open())
            {
                cout << "Error opening Interactions File in AIO" << Fhandler->fname_interactions << endl;
                exit(1);
            }

            fp_I.seekg(  ini_filepos*sizeof(type_precision) , ios::beg );
            if(fp_I.fail())
            {
                cout << "Error reading Interactions File while positioning in AIO! " <<  ini_filepos << endl;
                exit(1);
            }
            list< pair<int,int> >* excl_List = Fhandler->excl_List;

            int buff_pos=0;
            int file_pos;



            for (int i=Fhandler->ini_IDX_interactions; i < Fhandler->end_IDX_interactions; i++)
            {
                for (list<  pair<int,int>  >::iterator it=excl_List->begin(); it != excl_List->end(); ++it)
                {
                    file_pos = i*Fhandler->fileN+ it->first;
                    fp_I.seekg(  file_pos*sizeof(type_precision) , ios::beg );
                    chunk_size_buff = it->second;

                    fp_I.read ((char*)&(Fhandler->interaction_data[buff_pos]),sizeof(type_precision)*chunk_size_buff);
                    if(fp_I.fail())
                    {
                        cout << "Error reading Interaction File! in AIO, size:"<<  chunk_size_buff<< " pos:" << file_pos <<endl;
                        exit(1);
                    }
                    buff_pos += chunk_size_buff;
                }
            }

//            replace_nans(ar_nan_idxs, Fhandler->numInter, backupAR, Fhandler->n , 1);
//            replace_nans_avgs(Fhandler->numInter, backupAR, Fhandler->n, r, ar_nan_idxs);

          //  matlab_print_matrix("int",Fhandler->n,Fhandler->numInter,Fhandler->interaction_data);



        }

        //!*************************************************

        if(Fhandler->storeBin)
        {

            fp_allResults.open((Fhandler->fnameOutFiles + "_sto.dbin").c_str(),ios::out | ios::binary | ios::trunc);
            if(fp_allResults == 0)
            {
                cout << "Error Creating File "<< (Fhandler->fnameOutFiles + "_sto.bin") << endl;
                exit(1);
            }

        }



    }
    else
    {
        //cout << "\nPreping files\n" << flush;
        fp_Y = fopen("tempY.bin", "w+b");
        if(fp_Y == 0)
        {
            cout << "Error creating temp File Y " << endl;
            exit(1);
        }
        type_precision* tempbuff1 = new type_precision[Fhandler->n*Fhandler->y_blockSize];
        fwrite(tempbuff1, sizeof(type_precision), Fhandler->n*Fhandler->y_blockSize, fp_Y);

        delete []tempbuff1;


        ofstream fp_Art;
        fp_Art.open ("tempAR.bin",ios::out | ios::binary );
        if(!fp_Art.is_open())
        {
            cout << "Error creating temp File AR "  << endl;
            exit(1);
        }
        type_precision* tempbuff2 = new type_precision[Fhandler->n*Fhandler->Ar_blockSize*Fhandler->r];
        fp_Art.write((char*)tempbuff2, sizeof(type_precision)*Fhandler->n*Fhandler->Ar_blockSize*Fhandler->r);
        fp_Art.close();
        fp_Ar.open ("tempAR.bin",ios::in | ios::binary );

        delete []tempbuff2;


        //cout << "\nEnd preping files\n" << flush;

    }

    //pthread_barrier_wait(&(Fhandler->finalize_barrier));//for testing only


    bool Local_not_done = true;
    Fhandler->reset_wait = false;
    #ifdef DEBUG
    cout << " m:" << Fhandler->Ar_to_readSize << " mb:" <<   Fhandler->Ar_file_blocksize << endl<<flush ;
    #endif
    while(Local_not_done)
    {

        while(!Fhandler->empty_buffers.empty() && Fhandler->y_to_readSize && Fhandler->not_done)
        {

            tmp_y_blockSize = Fhandler->y_blockSize;
            if(Fhandler->y_to_readSize < Fhandler->y_blockSize)
                tmp_y_blockSize = Fhandler->y_to_readSize;

            Fhandler->y_to_readSize -= tmp_y_blockSize;
            size_buff = Fhandler->n * tmp_y_blockSize;



            pthread_mutex_lock(&(Fhandler->m_buff_upd));


            type_buffElement* tobeFilled = Fhandler->empty_buffers.front();
            Fhandler->empty_buffers.pop();


            tobeFilled->size = tmp_y_blockSize;
            //cout << "tbz:" << tmp_y_blockSize << " " << flush;

            if(Fhandler->fakefiles)
            {
                fseek ( fp_Y , 0 , SEEK_SET );
                size_t result = fread (tobeFilled->buff,sizeof(type_precision),size_buff,fp_Y);
                result++;
                int old_seed = Fhandler->seed;
                srand (old_seed);
                re_random_vec(tobeFilled->buff, size_buff );
                re_random_vec_nan(tobeFilled->buff, size_buff );
                Fhandler->seed += 75;
            }
            else
            {
                list< pair<int,int> >* excl_List = Fhandler->excl_List;


                int chunk_size_buff;
                int buff_pos=0;
                int file_block_offset;


                for(int i = 0; i < tmp_y_blockSize; i++)
                {
                    for (list<  pair<int,int>  >::iterator it=excl_List->begin(); it != excl_List->end(); ++it)
                    {
                        file_block_offset = Fhandler->fileN*i + it->first;
                        fseek ( fp_Y , (y_file_pos+file_block_offset)*sizeof(type_precision) , SEEK_SET );
                        chunk_size_buff = it->second;

                        size_t result = fread (&tobeFilled->buff[buff_pos],sizeof(type_precision),chunk_size_buff,fp_Y); result++;
                        buff_pos += chunk_size_buff;


                    }
                }

                y_file_pos += tmp_y_blockSize*Fhandler->fileN;


                if(Fhandler->y_to_readSize <= 0)
                {
                    fseek ( fp_Y , 0 , SEEK_SET );
                    //cout << "resetY" << endl;
                }
            }



            Fhandler->full_buffers.push(tobeFilled);

            pthread_mutex_unlock(&(Fhandler->m_buff_upd));

            pthread_mutex_lock(&(Fhandler->m_read));
            pthread_cond_signal( &(Fhandler->condition_read ));
            pthread_mutex_unlock(&(Fhandler->m_read));

        }

        while(!Fhandler->ar_empty_buffers.empty() && Fhandler->Ar_to_readSize && Fhandler->not_done)
        {

            tmp_ar_blockSize = Fhandler->Ar_file_blocksize;
            if(Fhandler->Ar_to_readSize < Fhandler->Ar_file_blocksize)
                tmp_ar_blockSize = Fhandler->Ar_to_readSize;


            #ifdef DEBUG
            cout << "tmp_ar_blockSize:" << tmp_ar_blockSize << " " ;
            #endif

            Fhandler->Ar_to_readSize -= tmp_ar_blockSize;
            size_buff = Fhandler->n * tmp_ar_blockSize*Fhandler->r;

            pthread_mutex_lock(&(Fhandler->m_buff_upd));
            type_buffElement* tobeFilled = Fhandler->ar_empty_buffers.front();
            Fhandler->ar_empty_buffers.pop();


            tobeFilled->size = tmp_ar_blockSize;
            #ifdef DEBUG
            cout << tmp_ar_blockSize << " " << Fhandler->Ar_file_blocksize << " " ;
            #endif

            if(Fhandler->fakefiles)
            {
                fp_Ar.seekg ( 0 ,  ios::beg  );
                fp_Ar.read ((char*)tobeFilled->buff,sizeof(type_precision)*size_buff);

                re_random_vec(tobeFilled->buff , Fhandler->n * tmp_ar_blockSize*Fhandler->r );
                re_random_vec_nan(tobeFilled->buff , Fhandler->n * tmp_ar_blockSize*Fhandler->r );

            }
            else
            {
                #ifdef DEBUG
                cout << " "<< Fhandler->use_dosages << " " <<  Fhandler->add_dosages <<" " <<  Fhandler->model << " " << Fhandler->fileR << " " << Fhandler->dosage_skip << endl;
                #endif

                list< pair<int,int> >* excl_List = Fhandler->excl_List;

                int chunk_size_buff;
                int buff_pos = 0;

                int file_block_offset = 0;



                float * file_data = tobeFilled->buff;
                if(Fhandler->use_interactions)
                {
                    file_data = new float[tmp_ar_blockSize*Fhandler->fileR*Fhandler->n];
                }

                if(Fhandler->use_dosages)
                {
                    //cout << "loading " << endl << flush;

                    for(int i = 0; i < tmp_ar_blockSize; i++)
                    {
                        buff_pos = 0;
                        for(int ii = 0; ii < (Fhandler->fileR-Fhandler->dosage_skip); ii++)
                        {

                            for (list<  pair<int,int>  >::iterator it=excl_List->begin(); it != excl_List->end(); ++it)
                            {
                                file_block_offset = i*(Fhandler->fileN*Fhandler->fileR)+ ii*Fhandler->fileN + it->first;
                                fp_Ar.seekg(  (ar_file_pos+file_block_offset)*sizeof(type_precision) , ios::beg );
                                if(fp_Ar.fail())
                                {
                                    cout << "Error reading AR File while positioning when using dosages! " <<  (ar_file_pos+file_block_offset) << endl;
                                    exit(1);
                                }

                                chunk_size_buff = it->second;

                                fp_Ar.read ((char*)&(Fhandler->ArDosage[buff_pos]),sizeof(type_precision)*chunk_size_buff);
                                if(fp_Ar.fail())
                                {
                                    cout << "Error reading AR File when using dosages" <<  (ar_file_pos+file_block_offset) << endl;
                                    exit(1);
                                }

                                buff_pos += chunk_size_buff;
                            }
                        }

                        if(Fhandler->add_dosages)
                        {
                            //cout << "adding ";
                            cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                Fhandler->n, 1, Fhandler->fileR-Fhandler->dosage_skip, 1.0, Fhandler->ArDosage, Fhandler->n, Fhandler->dosages,Fhandler->fileR-Fhandler->dosage_skip ,
                                    0.0, &(file_data[i*Fhandler->n]), Fhandler->n);
                        }
                        else
                        {
                            //cout << "mult ";
                            for(int ii = 0; ii < Fhandler->fileR; ii++)
                            {
                                for(int k=0; k < Fhandler->n; k++)
                                {
                                   file_data[i*Fhandler->n*Fhandler->fileR + ii*Fhandler->n+k] = Fhandler->ArDosage[Fhandler->n*ii+k] * Fhandler->dosages[ii];
                                }
                            }
                        }


                    }
                    ar_file_pos += tmp_ar_blockSize*Fhandler->fileN*Fhandler->fileR;


                    //cout << "endloading " << endl << flush;

                }
                else
                {
                    //cout << tmp_ar_blockSize*Fhandler->fileR << " " << Fhandler->n << flush;
                    for(int i = 0; i < tmp_ar_blockSize*Fhandler->fileR; i++)
                    {
                       // cout << "    b:" << i  << ":" << tmp_ar_blockSize*Fhandler->fileR << ":" << excl_List->size() << "       " << flush;
                        //int idx = 0;
                        for (list<  pair<int,int>  >::iterator it=excl_List->begin(); it != excl_List->end(); ++it)
                        {
                            //idx++;
                            //cout << " "<< idx  << " "<< flush;

                            file_block_offset = Fhandler->fileN*i + it->first;


                            fp_Ar.seekg((ar_file_pos+file_block_offset)*sizeof(type_precision), ios::beg);
                            if(fp_Ar.fail())
                            {
                                cout << "Error reading AR File while positioning! " <<  (ar_file_pos+file_block_offset) << endl;

                                exit(1);
                            }

                            chunk_size_buff = it->second;


//                            if(buff_pos >= Fhandler->Ar_file_blocksize*Fhandler->fileR*Fhandler->n)
//                            {
                            //cout << " "<< buff_pos<< ":" << Fhandler->Ar_file_blocksize*Fhandler->fileR*Fhandler->n << ":"<< chunk_size_buff << " "<< flush;
                                //cout << buff_pos<< ":" << chunk_size_buff << " "<< flush;
                               // exit(1);
                            //}

                            fp_Ar.read ((char*)&(file_data[buff_pos]),sizeof(type_precision)*chunk_size_buff);// result++;
                            if(fp_Ar.fail())
                            {

                                cout << "Error reading AR File! " <<  (ar_file_pos+file_block_offset) << ":" << buff_pos << ":" << chunk_size_buff<<  endl;
                                if((ar_file_pos+file_block_offset) > Fhandler->fileM*Fhandler->fileR*Fhandler->fileN)
                                {
                                    cout << "File pos exceeds file size!" << endl << flush;

                                }

                                if(ar_file_pos > Fhandler->fileM*Fhandler->fileR*Fhandler->fileN)
                                {
                                    cout << "File pos counter exceeds file size!" << endl << flush;

                                }
                                exit(1);
                            }
                            buff_pos += chunk_size_buff;
                            //cout << "|" << buff_pos<< flush;


                        }
                    }
                    ar_file_pos += tmp_ar_blockSize*Fhandler->fileN*Fhandler->fileR;


                }

                if(Fhandler->use_interactions)
                {
                    //cout << "int" << endl << flush;

                    //matlab_print_matrix("orig",Fhandler->n,tmp_ar_blockSize,file_data);

                    if(Fhandler->modelR != 1)
                    {
                        cout << "unexpected error where Fhandler->modelR !=1 with value of " <<Fhandler->modelR << endl;
                        exit(1);
                    }
                    if(Fhandler->use_multiple_interaction_sets)
                    {
                        generate_singleinteraction_multset(Fhandler->interaction_data,Fhandler->numInter,file_data,tmp_ar_blockSize,Fhandler->n,tobeFilled->buff,Fhandler->keep_depVar);
                        tobeFilled->size = tmp_ar_blockSize*Fhandler->numInter;
                    }
                    else
                    {
                        generate_multinteraction_singleset(Fhandler->interaction_data,Fhandler->numInter,file_data,tmp_ar_blockSize,Fhandler->n,tobeFilled->buff,Fhandler->keep_depVar);

                    }




                    //matlab_print_matrix("exp",Fhandler->n,tmp_ar_blockSize*Fhandler->numInter,tobeFilled->buff);

                    delete []file_data;
                }








            }

            if(Fhandler->Ar_to_readSize <= 0)
            {
                Fhandler->Ar_to_readSize = Fhandler->Ar_Amount;
                ar_file_pos = 0;
                fp_Ar.seekg(0, ios::beg);

                if(fp_Ar.fail())
                {
                    cout << "Error reading AR file again!" << flush <<  endl;
                    exit(1);
                }
                #ifdef DEBUG
                cout << "reseting ar" << flush;
                #endif
            }

            Fhandler->ar_full_buffers.push(tobeFilled);
            //  cout << "\nStoring " << tobeFilled << endl;
            pthread_mutex_unlock(&(Fhandler->m_buff_upd));

            pthread_mutex_lock(&(Fhandler->m_read));
            pthread_cond_signal( &(Fhandler->condition_read ));
            pthread_mutex_unlock(&(Fhandler->m_read));

        }
        //B write

        while(!Fhandler->write_full_buffers.empty() && Local_not_done && seq_count < max_secuential_write_count)
        {
            seq_count++;
            //cout << "S" << " ";

            pthread_mutex_lock(&(Fhandler->out_buff_upd));
            list < resultH >* tobeWritten = Fhandler->write_full_buffers.front();
            Fhandler->write_full_buffers.pop();


            if(Fhandler->fakefiles)
            {

            }

            if(!Fhandler->fakefiles && !tobeWritten->empty())
            {

                if(!Fhandler->storePInd && Fhandler->storeBin)
                {
                    resultH first = tobeWritten->front();
                    fp_allResults.write( (char*)&first.Y_name_idx,sizeof(int));

                    int size_block = tobeWritten->size();
                    fp_allResults.write( (char*)&size_block,sizeof(int));
                    fp_allResults.write( (char*)&first.ARoffset,sizeof(int));
                }


                for (list<  resultH  >::iterator it=tobeWritten->begin(); it != tobeWritten->end(); ++it)
                {
                    resultH current = *it;
                    uint8_t size_p = current.res_p.size();

                    uint16_t R2=-166;
                    float16(R2, current.R2);

                    //cout << current.R2 << ":";
//                    current.R2=-1.0;
//
//                    float32(&(current.R2),R2);
                    //cout << current.R2 << "\t";


                    if(Fhandler->storeBin)
                    {
                        fp_allResults.write( (char*)&current.nUsed,sizeof(uint16_t));
                        fp_allResults.write( (char*)&R2,sizeof(uint16_t));
                    }

                    if(!Fhandler->storePInd && Fhandler->storeBin)
                    {
                        fp_allResults.write( (char*)&size_p,sizeof(uint8_t));
                    }

                    for (list<  result  >::iterator it2=current.res_p.begin(); it2 != current.res_p.end(); ++it2)
                    {
                        result res_p = *it2;
                        if((res_p.P <= Fhandler->min_p_disp && current.R2 >= Fhandler->min_R2_disp )|| Fhandler->storePInd)
                        {
                            Fhandler->io_overhead = "!";//let the user know a result was found
                            string Aname = " ";

                            if(res_p.AL_name_idx < Fhandler->l)
                            {
                                //cout << (int)res_p.AL_name_idx << " " << Fhandler->l << flush;
                                Aname = Fhandler->ALnames[res_p.AL_name_idx] + ":";
                            }

                            Aname += Fhandler->ARnames[res_p.AR_name_idx+current.ARoffset];


                            //cout << res_p.AR_name_idx+ current.ARoffset<< " " <<flush;
                            fp_sigResults << Fhandler->Ynames[current.Y_name_idx] << "\t" << Aname << "\t" << current.nUsed << "\t" <<   std::fixed << std::setprecision(2) << current.nUsedPct << std::setprecision(std :: numeric_limits<float> :: digits10 + 1)
                                    << "\t" << res_p.B  << "\t" << std::setprecision(4) << current.R2 << std::setprecision(std :: numeric_limits<float> :: digits10 + 1)
                                          << "\t" << res_p.SE   << "\t" << res_p.T   << "\t" << std::scientific  << res_p.P <<  std::fixed << std::setprecision(std :: numeric_limits<float> :: digits10 + 1)  << endl;
//                                     cout <<  std::setprecision(std :: numeric_limits<float> :: digits10 + 1) << Fhandler->Ynames[current.Y_name_idx] << "\t" << Aname << "\t" << current.nUsed << "\t" << current.nUsedPct
//                                    << "\t" << res_p.B  << "\t" << current.R2
//                                          << "\t" << res_p.SE   << "\t" << res_p.T   << "\t" << std::scientific  << res_p.P << std::fixed  << endl;


                        }

                        if(Fhandler->storeBin)
                        {
                            if(Fhandler->disp_cov && !Fhandler->storePInd)
                            {
                                fp_allResults.write( (char*)&res_p.AL_name_idx,sizeof(uint8_t));
                            }
                            if(!Fhandler->storePInd)
                            {
                                fp_allResults.write( (char*)&res_p.AR_name_idx,sizeof(uint16_t));
                            }
                            fp_allResults.write( (char*)&res_p.B,sizeof(float));
                            fp_allResults.write( (char*)&res_p.SE,sizeof(float));
                        }

                    }


                }
            }



            Fhandler->write_empty_buffers.push(tobeWritten);
            //  cout << "\nStoring " << tobeWritten << endl;

            //cout << "Done Storing" << endl;





            pthread_mutex_unlock(&(Fhandler->out_buff_upd));



            pthread_mutex_lock(&(Fhandler->m_read));
            pthread_cond_signal( &(Fhandler->condition_read ));
            pthread_mutex_unlock(&(Fhandler->m_read));

        }

        seq_count = 0;


        if(!Fhandler->not_done && Fhandler->write_full_buffers.empty())
                {Local_not_done = false;}


#ifdef WINDOWS
        SYSTEMTIME time;
        GetSystemTime(&time);

        timeToWait.tv_sec = time.wSecond + 500/1000;
        long int morenanos = (500%1000)*1000000;
        timeToWait.tv_nsec = time.wMilliseconds*1000 + morenanos ;
#else
        clock_gettime(CLOCK_REALTIME, &timeToWait);
        //timeToWait.tv_nsec += 10000000;
        timeToWait.tv_nsec += 10000;
#endif

        pthread_mutex_lock(&(Fhandler->m_more));
        pthread_cond_timedwait( &(Fhandler->condition_more), &(Fhandler->m_more), &timeToWait );
        pthread_mutex_unlock( &(Fhandler->m_more ));

        pthread_mutex_lock(&(Fhandler->m_read));
        pthread_cond_signal( &(Fhandler->condition_read ));
        pthread_mutex_unlock(&(Fhandler->m_read));


    }



    //barrier
    pthread_barrier_wait(&(Fhandler->finalize_barrier));

    {
    type_buffElement* tmp;

    if(Fhandler->currentReadBuff)
    {
        Fhandler->full_buffers.push(Fhandler->currentReadBuff);
        Fhandler->currentReadBuff=0;
    }

    while(!Fhandler->full_buffers.empty())
    {
       tmp= Fhandler->full_buffers.front();
       Fhandler->full_buffers.pop();
       delete []tmp->buff;
       delete tmp;
    }

    while(!Fhandler->empty_buffers.empty())
    {
       tmp= Fhandler->empty_buffers.front();
       Fhandler->empty_buffers.pop();
        delete []tmp->buff;
        delete tmp;

    }

    if(Fhandler->Ar_currentReadBuff)
    {
        Fhandler->ar_full_buffers.push(Fhandler->Ar_currentReadBuff);
        Fhandler->Ar_currentReadBuff=0;
    }

    while(!Fhandler->ar_full_buffers.empty())
    {
       tmp= Fhandler->ar_full_buffers.front();
       Fhandler->ar_full_buffers.pop();
       delete []tmp->buff;
       delete tmp;
    }

    while(!Fhandler->ar_empty_buffers.empty())
    {
       tmp= Fhandler->ar_empty_buffers.front();
       Fhandler->ar_empty_buffers.pop();
       delete []tmp->buff;
       delete tmp;
    }

    while(!Fhandler->write_full_buffers.empty())
    {
       list < resultH >* tmp2= Fhandler->write_full_buffers.front();
       Fhandler->write_full_buffers.pop();
       delete tmp2;
    }

    while(!Fhandler->write_empty_buffers.empty())
    {
       list < resultH >* tmp2= Fhandler->write_empty_buffers.front();
       Fhandler->write_empty_buffers.pop();
       delete tmp2;
    }


    if(Fhandler->use_interactions)
    {
        delete [](Fhandler->interaction_data);
        fp_I.close();
    }

    }



        fclose(fp_Y);
        fp_Ar.close();

        fp_sigResults.close();
        fp_allResults.close();



        return 0;



}

void AIOwrapper::load_ARblock(type_precision** Ar, int &Ar_blockSize)
{

    //int status;
   // int createstatus = 0;
    //cout<<"^";

    while(Fhandler->ar_full_buffers.empty())
    {
        pthread_mutex_lock(&(Fhandler->m_more));
        pthread_cond_signal( &(Fhandler->condition_more ));
        pthread_mutex_unlock(&(Fhandler->m_more));

        io_overhead = "#";

        pthread_mutex_lock(&(Fhandler->m_read));
        pthread_cond_wait( &(Fhandler->condition_read), &(Fhandler->m_read ));
        pthread_mutex_unlock(&(Fhandler->m_read));
    }


    //!read new rdy buffer
    pthread_mutex_lock(&(Fhandler->m_buff_upd));

    if(Fhandler->Ar_currentReadBuff)
    {
        Fhandler->ar_empty_buffers.push(Fhandler->Ar_currentReadBuff);
    }

    Fhandler->Ar_currentReadBuff = Fhandler->ar_full_buffers.front();
    Fhandler->ar_full_buffers.pop();

    //cout << "\nReading " << Fhandler->Ar_currentReadBuff << endl;
    Fhandler->Ar = Fhandler->Ar_currentReadBuff->buff;
    Ar_blockSize = Fhandler->Ar_currentReadBuff->size;


    pthread_mutex_unlock(&(Fhandler->m_buff_upd));



    pthread_mutex_lock(&(Fhandler->m_more));
    pthread_cond_signal( &(Fhandler->condition_more ));
    pthread_mutex_unlock(&(Fhandler->m_more));


    (*Ar) = Fhandler->Ar;






}

void AIOwrapper::load_Yblock(type_precision** Y, int &y_blockSize)
{

    //int status;
    //int createstatus = 0;

    while(Fhandler->full_buffers.empty())
    {

        pthread_mutex_lock(&(Fhandler->m_more));
        pthread_cond_signal( &(Fhandler->condition_more ));
        pthread_mutex_unlock(&(Fhandler->m_more));

        io_overhead = "$";

        pthread_mutex_lock(&(Fhandler->m_read));
        pthread_cond_wait( &(Fhandler->condition_read), &(Fhandler->m_read ));
        pthread_mutex_unlock(&(Fhandler->m_read));

    }

    //!read new rdy buffer
    pthread_mutex_lock(&(Fhandler->m_buff_upd));


        if(Fhandler->currentReadBuff)
        {
            Fhandler->empty_buffers.push(Fhandler->currentReadBuff);
        }
        Fhandler->currentReadBuff = Fhandler->full_buffers.front();
        Fhandler->full_buffers.pop();


    Fhandler->Yb = Fhandler->currentReadBuff->buff;
    y_blockSize = Fhandler->currentReadBuff->size;

    (*Y) = Fhandler->Yb;



    pthread_mutex_unlock(&(Fhandler->m_buff_upd));



    pthread_mutex_lock(&(Fhandler->m_more));
    pthread_cond_signal( &(Fhandler->condition_more ));
    pthread_mutex_unlock(&(Fhandler->m_more));

    //matlab_print_matrix("Y",Fhandler->n,y_blockSize,*Y);


}


void AIOwrapper::prepare_Y(int y_blockSize, int n, int totalY)
{
    //for fake files

    Fhandler->seed = 1337;
    srand (Fhandler->seed);

    Fhandler->y_blockSize = y_blockSize;

    Fhandler->n= n;
    Fhandler->Y_Amount=totalY;
    Fhandler->y_to_readSize = Fhandler->Y_Amount;
    Fhandler->buff_count = max(3,1+(totalY+ y_blockSize - 1)/y_blockSize) ;
    //cout << "buffcount " << Fhandler->buff_count;


    Fhandler->currentReadBuff = 0;
    type_buffElement* tmp;

    for(int i = 0; i< Fhandler->buff_count  ; i++)
    {
        tmp = new type_buffElement();
        tmp->buff = new type_precision[Fhandler->n*Fhandler->y_blockSize];
        tmp->size = y_blockSize;
        Fhandler->empty_buffers.push(tmp);
        Fhandler->Yb = tmp->buff;
    }




    pthread_mutex_init(&(Fhandler->m_buff_upd), NULL);
    pthread_mutex_init(&(Fhandler->out_buff_upd), NULL);
    pthread_mutex_init(&(Fhandler->m_more), NULL);
    pthread_mutex_init(&(Fhandler->m_read), NULL);

    pthread_attr_init(&(Fhandler->attr));
    pthread_attr_setdetachstate(&(Fhandler->attr), PTHREAD_CREATE_JOINABLE);

    pthread_create( &(Fhandler->iothread),&(Fhandler->attr), AIOwrapper::async_io, (void*)Fhandler);



}

void AIOwrapper::getCurrentWriteBuffers(list < resultH >* &sigResults)
{
    while(Fhandler->write_empty_buffers.empty())
    {
        pthread_mutex_lock(&(Fhandler->m_more));
        pthread_cond_signal( &(Fhandler->condition_more ));
        pthread_mutex_unlock(&(Fhandler->m_more));

        io_overhead = "W";

        pthread_mutex_lock(&(Fhandler->m_read));
        pthread_cond_wait( &(Fhandler->condition_read), &(Fhandler->m_read ));
        pthread_mutex_unlock(&(Fhandler->m_read));
    }
    pthread_mutex_lock(&(Fhandler->out_buff_upd));

    sigResults = Fhandler->write_empty_buffers.front();
    Fhandler->write_empty_buffers.pop();

    pthread_mutex_unlock(&(Fhandler->out_buff_upd));


    sigResults->clear();


    pthread_mutex_lock(&(Fhandler->m_more));
    pthread_cond_signal( &(Fhandler->condition_more ));
    pthread_mutex_unlock(&(Fhandler->m_more));
}

void AIOwrapper::write_OutFiles(list < resultH >* &sigResults)
{

    pthread_mutex_lock(&(Fhandler->out_buff_upd));


    Fhandler->write_full_buffers.push(sigResults);
    sigResults = 0;
    //cout << Fhandler->write_full_buffers.size() << endl;

    if(Fhandler->io_overhead.compare(""))
    {
        io_overhead = Fhandler->io_overhead;
    }
    else
    {
        Fhandler->io_overhead = "";
    }

    pthread_mutex_unlock(&(Fhandler->out_buff_upd));


    pthread_mutex_lock(&(Fhandler->m_more));
    pthread_cond_signal( &(Fhandler->condition_more ));
    pthread_mutex_unlock(&(Fhandler->m_more));
}





void AIOwrapper::prepare_OutFiles()
{


    int buff_count = 500000;

    list < resultH >* tmp;


    for(int i = 0; i< buff_count  ; i++)
    {
        tmp = new list < resultH >();
        Fhandler->write_empty_buffers.push(tmp);
    }
    Fhandler->currentWriteBuff = Fhandler->write_empty_buffers.front();
    Fhandler->write_empty_buffers.pop();


}





void AIOwrapper::reset_Y()
{
    //void *status;

    Fhandler->seed = 1337;

    cout << "ry" << flush;

    Fhandler->reset_wait = true;
    pthread_barrier_wait(&(Fhandler->finalize_barrier));

    pthread_mutex_lock(&(Fhandler->m_buff_upd));
    Fhandler->y_to_readSize = Fhandler->Y_Amount;

    if(Fhandler->currentReadBuff)
    {
        Fhandler->full_buffers.push(Fhandler->currentReadBuff);
        Fhandler->currentReadBuff=0;
    }

    while(!Fhandler->full_buffers.empty())
    {
        Fhandler->empty_buffers.push(Fhandler->full_buffers.front());
        for( int i = 0; i < Fhandler->n*Fhandler->y_blockSize; i++)
        {
            ((Fhandler->full_buffers.front())->buff)[i] = 0;
        }
        Fhandler->full_buffers.pop();
    }
    pthread_mutex_unlock(&(Fhandler->m_buff_upd));

    Fhandler->reset_wait = false;

    pthread_barrier_wait(&(Fhandler->finalize_barrier));

    pthread_mutex_lock(&(Fhandler->m_more));
    pthread_cond_signal( &(Fhandler->condition_more ));
    pthread_mutex_unlock(&(Fhandler->m_more));


}

void AIOwrapper::reset_AR()
{

}

void AIOwrapper::finalize_Y()
{


}

void AIOwrapper::prepare_AR( int desired_blockSize, int n, int totalRfile, int columnsAR)
{

    Fhandler->Ar = new type_precision[desired_blockSize*columnsAR*n];

    Fhandler->Ar_blockSize = desired_blockSize;


    Fhandler->r = columnsAR;
    Fhandler->Ar_Amount = totalRfile;
    Fhandler->Ar_to_readSize = Fhandler->Ar_Amount;
    #ifdef DEBUG
    cout << "pr:" << Fhandler->r << " pm:" << Fhandler->Ar_to_readSize << endl << flush;
    #endif



    int buff_count = max(10,4*(Fhandler->Ar_Amount+desired_blockSize-1)/desired_blockSize+1);

    #ifdef DEBUG
    cout << "buff_count" << buff_count << endl  << flush;
    #endif

    Fhandler->Ar_currentReadBuff = 0;
    type_buffElement* tmp;

    for(int i = 0; i< buff_count  ; i++)
    {
        tmp = new type_buffElement();
        tmp->buff = new type_precision[n*desired_blockSize*columnsAR];
        tmp->size = desired_blockSize;
//        memset(tmp->buff,-9,n*desired_blockSize*columnsAR);

        Fhandler->ar_empty_buffers.push(tmp);
        Fhandler->Ar = tmp->buff;
    }


}


void AIOwrapper::finalize_AR()
{

}

void AIOwrapper::removeALmissings(list< pair<int,int> >* excl_List, struct Settings params, int &Almissings)
{
    float* tempAL = new float[params.n*params.l];

    FILE *fp;
    fp = fopen((Fhandler->fnameAL+".fvd").c_str(), "rb");
    if(fp == 0)
    {
        cout << "Error Reading File " << Fhandler->fnameAL << endl;
        exit(1);
    }
    size_t result = fread (tempAL,sizeof(type_precision),params.n*params.l,fp); result++;
    fclose(fp);



    Almissings=0;

    for (int h=0; h < params.l; h++)
    {
        for (int i=0; i < params.n; i++)
        {
            if(isnan(tempAL[h*params.n+i]))
            {
                bool removed  = splitpair(i,excl_List,params.n );
                if(removed)
                {
                    Almissings++;
                }

            }
        }
    }


    if(Almissings > 0)
        cout << "Excluding " << Almissings << " from covariate missings" << endl;

    delete tempAL;


}





bool AIOwrapper::splitpair(int value, list< pair<int,int> >* excl_List, int n)
{


    int buffsize = 0;
    bool removed = false;

    //cout << excl_List->size();
    for (list<  pair<int,int>  >::iterator it=excl_List->begin(); it != excl_List->end(); ++it)
    {
        int inipos = it->first;
        int endpos = it->first+it->second-1;
        int origbuffsize = it->second;
        //cout << inipos<< "-" << endpos << ":" << origbuffsize << " ";


        if(value <=  endpos)
        {
           // cout << "(" << inipos<< ":" << endpos << ") ";

            if(value > inipos && value < endpos)
            {
                buffsize = value-inipos;
                excl_List->insert (it,make_pair(inipos,buffsize));
                buffsize = origbuffsize-buffsize-1;
                excl_List->insert (it,make_pair(value+1,buffsize));
//                cout << inipos<< ":" << endpos << " ";
//                it--;it--;
//                cout << it->first<< ":" << it->second << " ";
//                 it++;
//                cout << it->first<< ":" << it->second << " | ";
//                it++;

                excl_List->remove((*it));

                it = excl_List->end();
                removed = true;

            }
            else
            {
                if(value == inipos && inipos != endpos)
                {
                    excl_List->insert (it,make_pair(min(inipos+1,n),origbuffsize-1));
//                    cout << inipos<< ":" << endpos << " ";
//                    it--;
//                    cout << it->first<< ":" << it->second << " | ";
//                    it++;

                    excl_List->remove((*it));
                    it = excl_List->end();
                    removed = true;
                }
                else
                {
                     if(value == endpos && endpos != inipos)
                     {
                        excl_List->insert (it,make_pair(inipos,origbuffsize-1));
//                        cout << inipos<< ":" << endpos << " ";
//                        it--;
//                        cout << it->first<< ":" << it->second << " | ";
//                        it++;

                        excl_List->remove((*it));
                        it = excl_List->end();
                        removed = true;
                     }
                     else
                     {
                        if(value == inipos && value == endpos && endpos == inipos)
                        {
                            excl_List->remove((*it));
                            it = excl_List->end();
                            removed = true;
                        }
                     }
                }
            }
        }
    }

    return removed;
}

void AIOwrapper::load_AL(type_precision** AL)
{

    if(Fhandler->fakefiles)
    {
        FILE *fp;
        fp = fopen("tempAL.bin", "rb");
        if(fp == 0)
        {
            cout << "Error Reading File tempAL.bin" << endl;
            exit(1);
        }

        size_t result = fread (Fhandler->AL,sizeof(type_precision),Fhandler->l*Fhandler->n,fp);
        result++;
        fclose(fp);
        srand(22);
        re_random_vec(Fhandler->AL,Fhandler->n*Fhandler->l);
        re_random_vec_nan(Fhandler->AL,Fhandler->n*Fhandler->l);
        (*AL) = Fhandler->AL;
    }
    else
    {
        FILE *fp;
        fp = fopen((Fhandler->fnameAL+".fvd").c_str(), "rb");
        if(fp == 0)
        {
            cout << "Error Reading File " << Fhandler->fnameAL << endl;
            exit(1);
        }

        list< pair<int,int> >* excl_List = Fhandler->excl_List;

        int chunk_size_buff;
        int buff_pos=0;
        int file_pos;

        for (int i=0; i < Fhandler->l; i++)
        {
            for (list<  pair<int,int>  >::iterator it=excl_List->begin(); it != excl_List->end(); ++it)
            {

                file_pos = i*Fhandler->fileN+ it->first;
                fseek ( fp , file_pos*sizeof(type_precision) , SEEK_SET );
                chunk_size_buff = it->second;

                size_t result = fread (&(Fhandler->AL[buff_pos]),sizeof(type_precision),chunk_size_buff,fp); result++;
                buff_pos += chunk_size_buff;
            }
        }

        //cout << Fhandler->n;


//        size_t result = fread (Fhandler->AL,sizeof(type_precision),Fhandler->l*Fhandler->n,fp);
//
//        result++;
        fclose(fp);
    }



    (*AL) = Fhandler->AL;

    //matlab_print_matrix("AL",Fhandler->n,Fhandler->l,*AL);
}

void AIOwrapper::prepare_AL( int columnsAL, int n)
{

    Fhandler->AL = new type_precision[columnsAL*n];
    Fhandler->l=columnsAL;
    if(Fhandler->fakefiles)
    {
        FILE* fp_AL = fopen("tempAL.bin", "w+b");
        if(fp_AL == 0)
        {
            cout << "Error creating temp File AL "<< endl;
            exit(1);
        }
        fwrite(Fhandler->AL, sizeof(type_precision), n*columnsAL, fp_AL);
        fclose(fp_AL);
    }
}

void AIOwrapper::finalize_AL()
{
    delete []Fhandler->AL;
}


void AIOwrapper::read_excludeList(list< pair<int,int> >* excl, int &excl_count, int n, string fname_excludeList)
{

    ifstream fp_exL(fname_excludeList.c_str());
    if(fp_exL == 0)
    {
        cout << "Error reading exclude list file."<< endl;
        exit(1);
    }


    string line;
    excl_count = 0;

    int first,second;

     cout << "Excluding Ids: \n";

      while (std::getline(fp_exL, line) && second < n )
      {
            std::istringstream iss(line);

            iss >> first;

            iss >> second;
            cout << first << "-" << second << ", ";
            if(first > second)
            {
                cout << "\nPlease provide ordered pairs!\n";
                exit( 1 );
            }
            if(first < 1 || second < 1)
            {
                cout << "\nPlease provide valid ositive indixes!\n";
                exit( 1 );
            }

            for(int i = first-1; i <= second-1;i++ )
            {

                bool removed  = splitpair(i,excl, n );
                if(removed)
                {
                    excl_count++;
                }
            }


      }


    if(excl_count >= n)
    {
		cout << "\nExclusion List excluded all data!\n";
		cout << "Total Ids: " << n << "\nExcluded Ids: " << excl_count << endl;
		exit( 1 );
	}
	 cout << "Excluded: "  << excl_count << " Using: " << n-excl_count<< endl;



}


void AIOwrapper::read_dosages(string fname_dosages, int expected_count, float* vec)
{
    ifstream fp_dos(fname_dosages.c_str());
    if(fp_dos == 0)
    {
        cout << "Error reading dosages file."<< endl;
        exit(1);
    }
    int i;
    for (i=0; i < expected_count && !fp_dos.eof(); i++)
    {
       fp_dos >> vec[i];
       //cout << vec[i];
    }
    if(i!=expected_count)
    {
        cout << "The number of factors for the dosage model! does not coincide with the source data #facctors:" << i << ", expected:" <<expected_count << endl;
        exit(1);
    }

}


void AIOwrapper::free_databel_fvi( struct databel_fvi **fvi )
{
	free ((*fvi)->fvi_data);
	free (*fvi);
	*fvi = NULL;
}

FILE * AIOwrapper::fgls_fopen( const char *path, const char *mode )
{
	FILE * f;
	//char err[100];

	f = fopen( path, mode );
	if ( f == NULL )
	{
		cout << "\nError reading file: " << path << endl;
		exit( 1 );
	}
	return f;
}

void * AIOwrapper::fgls_malloc_impl( const char* file, long line, size_t size )
{
    void *buf;

    if ( (buf = malloc( size )) == NULL ) {
        cout<< "\nCouldn't allocate %ld bytes of memory in %s:%ld\n";
        exit(1);
    }

    return buf;
}

struct databel_fvi * AIOwrapper::load_databel_fvi( const char *path )
{
	FILE *f;
	databel_fvi *fvi;
	size_t data_size;

	f = fgls_fopen( path, "r" );

	fvi = (databel_fvi*) fgls_malloc( sizeof(databel_fvi) );
	// Header
	size_t result = fread( &fvi->fvi_header, sizeof(databel_fvi_header), 1, f );
	result++;
	// Labels
	data_size = (fvi->fvi_header.numVariables +fvi->fvi_header.numObservations ) *
		         fvi->fvi_header.namelength * sizeof(char);
	fvi->fvi_data = (char *) fgls_malloc ( data_size );
	// Load labels
	result = fread( fvi->fvi_data, 1, data_size, f );
	result++;

	fclose( f );

	return fvi;
}
