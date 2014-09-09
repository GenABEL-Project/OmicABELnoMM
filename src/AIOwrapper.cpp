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



        params.n = ALfvi->fvi_header.numObservations;
        Fhandler->fileN = params.n;
        Fhandler->fileR = params.r;
        params.m = ARfvi->fvi_header.numVariables/params.r;
        params.t = Yfvi->fvi_header.numVariables;
        params.l = ALfvi->fvi_header.numVariables;

        int yname_idx=0;//starting idx for names on ALfvi->data
        for(int i = 0; i < params.n; i++)
        {
            //Nnames.push_back(string(&(Yfvi->fvi_data[yname_idx])));
            yname_idx += Yfvi->fvi_header.namelength;

            //cout << i << ": " << string(&(Yfvi->fvi_data[yname_idx])) << "\t";
            //cout << i << ": " << Ynames[i] << "\t";
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
        int opt_mb = 100;

        params.mb = min(params.m, opt_mb);
        params.tb = min(params.t, opt_tb);



    }
    else
    {
        //other params come from outside
    }

    //params.fname_excludelist = "exclfile.txt";
    int excl_count = 0;
    int Almissings = 0;
    Fhandler->excl_List = new list< pair<int,int> >();



    if(params.fname_excludelist.size()==0)
    {
        (Fhandler->excl_List)->push_back( make_pair(0,params.n) );
    }
    else
    {
        read_excludeList( Fhandler->excl_List,excl_count,params.n,params.fname_excludelist);
    }

    if(!Fhandler->fakefiles)
    {

        removeALmissings(Fhandler->excl_List,params,Almissings);

    }

    params.n -= (excl_count + Almissings);

    if(params.dosages)
    {

        Fhandler->ArDosage = new float[Fhandler->fileR*params.n];
        Fhandler->dosages = new float[Fhandler->fileR];


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
        break;
        case 4://additive
            read_dosages(params.fname_dosages,Fhandler->fileR,Fhandler->dosages);
            params.r = 1;
            Fhandler->add_dosages = true;
        break;
        }
    }





    params.p = params.l + params.r;

    if(!Fhandler->fakefiles)
    {

        int Aname_idx=params.n*ARfvi->fvi_header.namelength;//skip the names of the rows
        if(Fhandler->use_dosages && Fhandler->add_dosages)
        {
            for(int i = 0; i < params.m; i++)
            {
                Fhandler->ARnames.push_back(string(&(ARfvi->fvi_data[Aname_idx])));
                Aname_idx += ARfvi->fvi_header.namelength*Fhandler->fileR;
            }
        }
        else
        {
            for(int i = 0; i < params.m*params.r; i++)
            {
                Fhandler->ARnames.push_back(string(&(ARfvi->fvi_data[Aname_idx])));
                Aname_idx += ARfvi->fvi_header.namelength;
            }
        }

        if(params.l > 255 || params.p > 255 || params.n > 65535)//can remove if fixed for output files
        {
            cout << "Warning, output binary format does not yet support current problem sizes for the provided p, l, r and n." << endl;
            cout << "Omitting outputfile." << endl;
        }


        //!write info file for results
        ofstream fp_InfoResults;

        fp_InfoResults.open((Fhandler->fnameOutFiles + "_sto.ibin").c_str(),ios::out | ios::binary | ios::trunc);
        if(fp_InfoResults == 0)
        {
            cout << "Error Creating File InfoResults.bin" << endl;
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


        Aname_idx=(params.n+1)*ALfvi->fvi_header.namelength;//skip the names of the rows + intercept
        fp_InfoResults.write( (char*)&ALfvi->fvi_data[Aname_idx],ALfvi->fvi_header.namelength*(params.l-1)*sizeof(char));

        Aname_idx=params.n*ARfvi->fvi_header.namelength;//skip the names of the rows
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

        int Yname_idx=params.n*Yfvi->fvi_header.namelength;//skip the names of the rows
        fp_InfoResults.write( (char*)&Yfvi->fvi_data[Yname_idx],Yfvi->fvi_header.namelength*params.t*sizeof(char));


    }




    //block size to keep mem under 1 gigabyte
//    int opt_block = params.n/(4*1000^3)*(1/(2*params.r));
//    int opt_tb = max(4*2000,opt_block);
//    int opt_mb = max(2000,opt_block);
//
    params.mb = min(params.m,params.mb);
    params.tb = min(params.t,params.tb);

    prepare_AL(params.l,params.n);
    prepare_AR(  params.mb,  params.n,  params.m,  params.r);
    prepare_OutFiles(params.mb, params.l+params.r);
    prepare_Y(params.tb, params.n, params.t);



    //

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
    FILE*  fp_Ar;
    ofstream fp_sigResults;
    ofstream fp_allResults;

    if(!Fhandler->fakefiles)
    {
        fp_Y = fopen((Fhandler->fnameY+".fvd").c_str(), "rb");
        if(fp_Y == 0)
        {
            cout << "Error Reading File Y " << Fhandler->fnameY << endl;
            exit(1);
        }

        fp_Ar = fopen((Fhandler->fnameAR+".fvd").c_str(), "rb");
        if(fp_Ar == 0)
        {
            cout << "Error Reading File Xr " << Fhandler->fnameAR << endl;
            exit(1);
        }


        fp_sigResults.open((Fhandler->fnameOutFiles + "_dis.txt").c_str(),ios::out | ios::trunc);
        if(fp_sigResults == 0)
        {
            cout << "Error Creating File " << (Fhandler->fnameOutFiles + "_dis.txt") << endl;
            exit(1);
        }

        fp_sigResults << "Phe\tSNP\tn\tnPct\tB\tR2\tSE\tT\tP" << endl;

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
        //fclose(fp_Y);
        delete []tempbuff1;


        fp_Ar = fopen("tempAR.bin", "w+b");
        if(fp_Ar == 0)
        {
            cout << "Error creating temp File AR "  << endl;
            exit(1);
        }
        type_precision* tempbuff2 = new type_precision[Fhandler->n*Fhandler->Ar_blockSize*Fhandler->r];
        fwrite(tempbuff2, sizeof(type_precision), Fhandler->n*Fhandler->Ar_blockSize*Fhandler->r, fp_Ar);
        //fclose(fp_Ar);
        delete []tempbuff2;


        //cout << "\nEnd preping files\n" << flush;

    }

    //pthread_barrier_wait(&(Fhandler->finalize_barrier));//for testing only


    bool Local_not_done = true;
    Fhandler->reset_wait = false;


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
                int file_pos;

                for(int i = 0; i < tmp_y_blockSize; i++)
                {
                    for (list<  pair<int,int>  >::iterator it=excl_List->begin(); it != excl_List->end(); ++it)
                    {
                        file_pos = Fhandler->fileN*i + it->first;
                        fseek ( fp_Y , file_pos*sizeof(type_precision) , SEEK_SET );
                        chunk_size_buff = it->second;

                        size_t result = fread (&tobeFilled->buff[buff_pos],sizeof(type_precision),chunk_size_buff,fp_Y); result++;
                        buff_pos += chunk_size_buff;


                    }
                }


                if(Fhandler->y_to_readSize <= 0)
                {
                    fseek ( fp_Y , 0 , SEEK_SET );
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




            tmp_ar_blockSize = Fhandler->Ar_blockSize;
            if(Fhandler->Ar_to_readSize < Fhandler->Ar_blockSize)
                tmp_ar_blockSize = Fhandler->Ar_to_readSize;

            Fhandler->Ar_to_readSize -= tmp_ar_blockSize;
            size_buff = Fhandler->n * tmp_ar_blockSize*Fhandler->r;

            pthread_mutex_lock(&(Fhandler->m_buff_upd));
            type_buffElement* tobeFilled = Fhandler->ar_empty_buffers.front();
            Fhandler->ar_empty_buffers.pop();


            tobeFilled->size = tmp_ar_blockSize;

            if(Fhandler->fakefiles)
            {
                fseek ( fp_Ar , 0 , SEEK_SET );
                size_t result = fread (tobeFilled->buff,sizeof(type_precision),size_buff,fp_Ar); result++;

                re_random_vec(tobeFilled->buff , Fhandler->n * tmp_ar_blockSize*Fhandler->r );
                re_random_vec_nan(tobeFilled->buff , Fhandler->n * tmp_ar_blockSize*Fhandler->r );

            }
            else
            {
                //cout << " "<< Fhandler->use_dosages << " " <<  Fhandler->add_dosages <<" " <<  Fhandler->model << " " << Fhandler->fileR << " " << Fhandler->dosage_skip << endl;


                list< pair<int,int> >* excl_List = Fhandler->excl_List;

                int chunk_size_buff;

                int file_pos;
                int buff_pos = 0;

                if(Fhandler->use_dosages)
                {


                    for(int i = 0; i < tmp_ar_blockSize; i++)
                    {
                        buff_pos=0;
                        for(int ii = 0; ii < (Fhandler->fileR-Fhandler->dosage_skip); ii++)
                        {

                            for (list<  pair<int,int>  >::iterator it=excl_List->begin(); it != excl_List->end(); ++it)
                            {
                                file_pos = i*(Fhandler->fileN*Fhandler->fileR)+ ii*Fhandler->fileN + it->first;
                                fseek ( fp_Ar , file_pos*sizeof(type_precision) , SEEK_SET );

                                chunk_size_buff = it->second;

                                size_t result = fread (&(Fhandler->ArDosage[buff_pos]),sizeof(type_precision),chunk_size_buff,fp_Ar); result++;
                                buff_pos += chunk_size_buff;
                            }
                        }

                        if(Fhandler->add_dosages)
                        {
                            //cout << "adding ";
                            cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                Fhandler->n, 1, Fhandler->fileR-Fhandler->dosage_skip, 1.0, Fhandler->ArDosage, Fhandler->n, Fhandler->dosages,Fhandler->fileR-Fhandler->dosage_skip ,
                                    0.0, &(tobeFilled->buff[i*Fhandler->n]), Fhandler->n);
                        }
                        else
                        {
                            //cout << "mult ";
                            for(int ii = 0; ii < Fhandler->fileR; ii++)
                            {
                                for(int k=0; k < Fhandler->n; k++)
                                {
                                   tobeFilled->buff[i*Fhandler->n*Fhandler->r + ii*Fhandler->n+k] = Fhandler->ArDosage[Fhandler->n*ii+k] * Fhandler->dosages[ii];
                                }
                            }
                        }

                    }



                }
                else
                {
                    for(int i = 0; i < tmp_ar_blockSize*Fhandler->r; i++)
                    {
                        for (list<  pair<int,int>  >::iterator it=excl_List->begin(); it != excl_List->end(); ++it)
                        {
                            file_pos = Fhandler->fileN*i + it->first;
                            fseek ( fp_Ar , file_pos*sizeof(type_precision) , SEEK_SET );

                            chunk_size_buff = it->second;
                            size_t result = fread (&tobeFilled->buff[buff_pos],sizeof(type_precision),chunk_size_buff,fp_Ar); result++;
                            buff_pos += chunk_size_buff;

                        }
                    }
                }




            }

            if(Fhandler->Ar_to_readSize <= 0)
            {
                Fhandler->Ar_to_readSize = Fhandler->Ar_Amount;
                fseek ( fp_Ar , 0 , SEEK_SET );
            }

            Fhandler->ar_full_buffers.push(tobeFilled);
            //  cout << "\nStoring " << tobeFilled << endl;
            pthread_mutex_unlock(&(Fhandler->m_buff_upd));

            pthread_mutex_lock(&(Fhandler->m_read));
            pthread_cond_signal( &(Fhandler->condition_read ));
            pthread_mutex_unlock(&(Fhandler->m_read));

        }
        //B write

        while(!Fhandler->write_full_buffers.empty() && Local_not_done)
        {

            pthread_mutex_lock(&(Fhandler->m_buff_upd));
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
                        if(res_p.P <= Fhandler->min_p_disp && current.R2 >= Fhandler->min_R2_disp)
                        {
                            string Aname = " ";

                            if(res_p.AL_name_idx < Fhandler->l)
                            {
                                //cout << (int)res_p.AL_name_idx << " " << Fhandler->l << flush;
                                Aname = Fhandler->ALnames[res_p.AL_name_idx] + ":";
                            }

                            Aname += Fhandler->ARnames[res_p.AR_name_idx+current.ARoffset];

                            //cout << Aname << " " <<flush;
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





            pthread_mutex_unlock(&(Fhandler->m_buff_upd));



            pthread_mutex_lock(&(Fhandler->m_read));
            pthread_cond_signal( &(Fhandler->condition_read ));
            pthread_mutex_unlock(&(Fhandler->m_read));

        }


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
        timeToWait.tv_nsec += 150;
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

    }


    //

    //pthread_exit(NULL);
        fclose(fp_Y);
        fclose(fp_Ar);
        fp_sigResults.close();
        fp_allResults.close();

        //cout << "\nexited io\n";

        return 0;

//
//            //!induce realistic fileread delay

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



//secueantial
//    int tmp_ar_blockSize = Fhandler->Ar_blockSize;
//    if(Fhandler->Ar_to_readSize < Fhandler->Ar_blockSize)
//        tmp_ar_blockSize = Fhandler->Ar_to_readSize;
//
//    Fhandler->Ar_to_readSize -= tmp_ar_blockSize;
//    int size_buff = Fhandler->n * tmp_ar_blockSize * Fhandler->r;
//    re_random_vec(Fhandler->Ar,size_buff);
//Ar_blockSize=tmp_ar_blockSize;


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

        io_overhead = "!";

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
    Fhandler->buff_count = min(3,(totalY+ y_blockSize - 1)/y_blockSize) ;
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
    pthread_mutex_lock(&(Fhandler->m_buff_upd));

    sigResults = Fhandler->write_empty_buffers.front();
    Fhandler->write_empty_buffers.pop();

    pthread_mutex_unlock(&(Fhandler->m_buff_upd));


    pthread_mutex_lock(&(Fhandler->m_more));
    pthread_cond_signal( &(Fhandler->condition_more ));
    pthread_mutex_unlock(&(Fhandler->m_more));
}

void AIOwrapper::write_OutFiles(list < resultH >* &sigResults)
{

    pthread_mutex_lock(&(Fhandler->m_buff_upd));


    Fhandler->write_full_buffers.push(sigResults);
    sigResults = 0;
    //cout << Fhandler->write_full_buffers.size() << endl;

    pthread_mutex_unlock(&(Fhandler->m_buff_upd));


    pthread_mutex_lock(&(Fhandler->m_more));
    pthread_cond_signal( &(Fhandler->condition_more ));
    pthread_mutex_unlock(&(Fhandler->m_more));
}





void AIOwrapper::prepare_OutFiles(int max_b_blockSize, int p)
{

    Fhandler->max_b_blockSize = max_b_blockSize;//useless?
    Fhandler->p=p;//useless?
    int buff_count = 4;

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

void AIOwrapper::prepare_AR( int desired_blockSize, int n, int totalR, int columnsAR)
{

    Fhandler->Ar = new type_precision[desired_blockSize*columnsAR*n];

    Fhandler->Ar_blockSize = desired_blockSize;
    Fhandler->r = columnsAR;
    Fhandler->Ar_Amount = totalR;
    Fhandler->Ar_to_readSize = Fhandler->Ar_Amount;

    int buff_count = 4;

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

void AIOwrapper::removeALmissings(list< pair<int,int> >* excl_List,struct Settings params, int &Almissings)
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

    int prev=0;
    list <int> missings;

    for (int h=0; h < params.l; h++)
    {
        for (int i=0; i < params.n; i++)
        {
            if(isnan(tempAL[h*params.n+i]))
            {
                splitpair(i,excl_List,params );
                missings.push_back(i);
            }
        }
    }

    missings.sort();
    missings.unique();
    Almissings = missings.size();

    delete tempAL;


}

void AIOwrapper::splitpair(int value, list< pair<int,int> >* excl_List,struct Settings params)
{


    int buffsize = 0;

    for (list<  pair<int,int>  >::iterator it=excl_List->begin(); it != excl_List->end(); ++it)
    {
        int inipos = it->first;
        int endpos = it->first+it->second-1;
        int origbuffsize = it->second;

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
            }
            else
            {
                if(value == inipos && inipos != endpos)
                {
                    excl_List->insert (it,make_pair(min(inipos+1,params.n),origbuffsize-1));
//                    cout << inipos<< ":" << endpos << " ";
//                    it--;
//                    cout << it->first<< ":" << it->second << " | ";
//                    it++;

                    excl_List->remove((*it));
                    it = excl_List->end();
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
                     }
                     else
                     {
                        if(value == inipos && value == endpos && endpos == inipos)
                        {
                            excl_List->remove((*it));
                            it = excl_List->end();
                        }
                     }
                }
            }
        }
    }
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


void AIOwrapper::read_excludeList(list< pair<int,int> >* excl, int &excl_count, int max_excl, string fname_excludeList)
{

    ifstream fp_exL(fname_excludeList.c_str());
    if(fp_exL == 0)
    {
        cout << "Error reading exclude list file."<< endl;
        exit(1);
    }


    string line;
    excl_count = 0;
    bool early_EOF;
    int first,second, prev_second;

     cout << "Excluding Ids: \n";

    std::getline(fp_exL, line);
    std::istringstream iss(line);
    iss >> first;
    early_EOF = iss.eof();
    iss >> prev_second;

    if(prev_second < first || early_EOF)
            prev_second = first;

    second = prev_second;


    if(first > max_excl)
    {
        excl->push_back( make_pair(0,max_excl) );
        cout << "\nNothing to Exclude!\n";
    }
    else
    {
        if(second > max_excl)
        {
                excl_count += max_excl-first+1;
                excl->push_back( make_pair(0,first-1) );
                cout << first << "-" << second << ", ";
        }
        else
        {

            cout << first << "-" << second << ", ";
            excl->push_back( make_pair(0,first - 1) );
            excl_count += second-(first-1);

            while (std::getline(fp_exL, line) && second < max_excl )
            {
                std::istringstream iss(line);

                iss >> first;
                early_EOF = iss.eof();
                iss >> second;
                if(prev_second >= first)
                {
                    cout << "\nPlease give an ordered Exlusion List!\n";
                    cout << "?? " << prev_second << "\n" << first << " ??"<< endl;
                    exit( 1 );
                }


                if(second < first || early_EOF )
                    second = first;


                if(second > max_excl)
                    excl_count += max_excl-first+1;
                else if(first < max_excl)
                    excl_count += second-(first-1);



                cout << first << "-" << second << ", ";

                excl->push_back( make_pair(prev_second,first - prev_second - 1) );


                prev_second = second;


            }
        }
    }




    if(excl_count >= max_excl)
    {
		cout << "\nExclusion List excluded all data!\n";
		cout << "Total Ids: " << max_excl << "\nExcluded Ids: " << excl_count << endl;
		exit( 1 );
	}
	 cout << "Excluded: "  << excl_count << " Using: " << max_excl-excl_count<< endl;



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
        cout << "not enough factor for the dosage model! required " << expected_count << endl;
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
		cout << "\nerror on fgls_fopen\n";
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
