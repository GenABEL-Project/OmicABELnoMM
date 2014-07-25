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






    if(!Fhandler->fakefiles)
    {
        Fhandler->fnameAL = params.fnameAL;
        Fhandler->fnameAR = params.fnameAR;
        Fhandler->fnameY = params.fnameY;
        Fhandler->fnameOutFiles = params.fnameOutFiles;


        Yfvi  = load_databel_fvi( (Fhandler->fnameY+".fvi").c_str() );
        ALfvi = load_databel_fvi( (Fhandler->fnameAL+".fvi").c_str() );
        ARfvi = load_databel_fvi( (Fhandler->fnameAR+".fvi").c_str() );
        params.n = ALfvi->fvi_header.numObservations;
        params.m = ARfvi->fvi_header.numVariables/params.r;
        params.t = Yfvi->fvi_header.numVariables;
        params.l = ALfvi->fvi_header.numVariables;

        int opt_tb = 1000;
        int opt_mb = 1000;

        params.mb = min(params.m, opt_tb);
        params.tb = min(params.t, opt_mb);

    }
    else
    {

    }

    //params.fname_excludelist = "exclfile.txt";
    int excl_count = 0;
    Fhandler->excl_List = new list< pair<int,int> >();

    if(params.fname_excludelist.size()==0)
    {
        (Fhandler->excl_List)->push_back( make_pair(0,params.n) );
    }
    else
    {
        read_excludeList( Fhandler->excl_List,excl_count,params.n,params.fname_excludelist);
    }

    params.n -= excl_count;

    params.p = params.l + params.r;


    //block size to keep mem under 1 gigabyte
//    int opt_block = params.n/(4*1000^3)*(1/(2*params.r));
//    int opt_tb = max(4*2000,opt_block);
//    int opt_mb = max(2000,opt_block);
//
//    params.mb = min(params.m,opt_tb);
//    params.tb = min(params.t,opt_mb);

    prepare_AL(params.l,params.n);
    prepare_AR(  params.mb,  params.n,  params.m,  params.r);
    prepare_OutFiles(params.mb, params.l+params.r);
    prepare_Y(params.tb, params.n, params.t);






}

void AIOwrapper::finalize()
{
    //cout << "f";
    //void *status;

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
    FILE*  fp_B;
    FILE*  fp_R;
    FILE*  fp_SD2;
    FILE*  fp_P;
    FILE*  fp_Ar;
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

        fp_B = fopen((Fhandler->fnameOutFiles+"_B.fvd").c_str(), "w+b");
        if(fp_B == 0)
        {
            cout << "Error Opening File B " << Fhandler->fnameOutFiles << "_B" << endl;
            exit(1);
        }
        fp_R = fopen((Fhandler->fnameOutFiles+"_R.fvd").c_str(), "w+b");
        if(fp_R == 0)
        {
            cout << "Error Opening File R " << Fhandler->fnameOutFiles << "_R" << endl;
            exit(1);
        }
        fp_SD2 = fopen((Fhandler->fnameOutFiles+"_SD2.fvd").c_str(), "w+b");
        if(fp_SD2 == 0)
        {
            cout << "Error Opening File SD2 " << Fhandler->fnameOutFiles << "_SD2" << endl;
            exit(1);
        }
        fp_P = fopen((Fhandler->fnameOutFiles+"_P.fvd").c_str(), "w+b");
        if(fp_P == 0)
        {
            cout << "Error Opening File P " << Fhandler->fnameOutFiles << "_P" << endl;
            exit(1);
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

        fp_B = fopen("tempB.bin", "w+b");
        if(fp_B == 0)
        {
            cout << "Error setting up temp File B " << endl;
            exit(1);
        }
        fp_R = fopen("tempR.bin", "w+b");
        if(fp_R == 0)
        {
            cout << "Error setting up temp File R " << endl;
            exit(1);
        }
        fp_SD2 = fopen("tempSD2.bin", "w+b");
        if(fp_SD2 == 0)
        {
            cout << "Error setting up temp File SD2 " << endl;
            exit(1);
        }
        fp_P = fopen("tempP.bin", "w+b");
        if(fp_P == 0)
        {
            cout << "Error setting up temp File P " << endl;
            exit(1);
        }
        //cout << "\nEnd preping files\n" << flush;

    }

    //pthread_barrier_wait(&(Fhandler->finalize_barrier));//for testing only


    Fhandler->not_done = true;
    Fhandler->reset_wait = false;


    while(Fhandler->not_done)
    {

        while(!Fhandler->empty_buffers.empty() && Fhandler->y_to_readSize)
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
                        file_pos = Fhandler->n*i + it->first;
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

        while(!Fhandler->ar_empty_buffers.empty() && Fhandler->Ar_to_readSize )
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

                list< pair<int,int> >* excl_List = Fhandler->excl_List;

                int chunk_size_buff;
                int buff_pos=0;
                int file_pos;

                for(int i = 0; i < tmp_ar_blockSize*Fhandler->r; i++)
                {
                    for (list<  pair<int,int>  >::iterator it=excl_List->begin(); it != excl_List->end(); ++it)
                    {
                        file_pos = Fhandler->n*i + it->first;
                        fseek ( fp_Ar , file_pos*sizeof(type_precision) , SEEK_SET );

                        chunk_size_buff = it->second;
                        size_t result = fread (&tobeFilled->buff[buff_pos],sizeof(type_precision),chunk_size_buff,fp_Ar); result++;
                        buff_pos += chunk_size_buff;


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

        while(!Fhandler->write_full_buffers.empty())
        {


            pthread_mutex_lock(&(Fhandler->m_buff_upd));
            type_buffElement* tobeWritten = Fhandler->write_full_buffers.front();
            Fhandler->write_full_buffers.pop();
            int size = Fhandler->p*Fhandler->b_blockSize;

            if(Fhandler->fakefiles)
            {
                fseek ( fp_B , 0 , SEEK_SET );
                fseek ( fp_R , 0 , SEEK_SET );
                fseek ( fp_SD2 , 0 , SEEK_SET );
                fseek ( fp_P , 0 , SEEK_SET );
            }
            fwrite (&(tobeWritten->buff[0]),sizeof(type_precision),size,fp_B);
            fwrite (&(tobeWritten->buff[Fhandler->max_b_blockSize*Fhandler->p]),sizeof(type_precision),Fhandler->b_blockSize,fp_R);
            fwrite (&(tobeWritten->buff[Fhandler->max_b_blockSize*(Fhandler->p+1)]),sizeof(type_precision),Fhandler->b_blockSize,fp_SD2);
            fwrite (&(tobeWritten->buff[Fhandler->max_b_blockSize*(Fhandler->p+2)]),sizeof(type_precision),size,fp_P);


            Fhandler->write_empty_buffers.push(tobeWritten);
            //  cout << "\nStoring " << tobeWritten << endl;
            pthread_mutex_unlock(&(Fhandler->m_buff_upd));

            pthread_mutex_lock(&(Fhandler->m_read));
            pthread_cond_signal( &(Fhandler->condition_read ));
            pthread_mutex_unlock(&(Fhandler->m_read));

        }




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

//        if(Fhandler->reset_wait)
//        {
//            pthread_barrier_wait(&(Fhandler->finalize_barrier));
//            //wait for main thread to reset everything
//
//            pthread_mutex_lock(&(Fhandler->m_buff_upd));
//            Fhandler->Ar_to_readSize = Fhandler->Ar_Amount;
//
//            if(Fhandler->Ar_currentReadBuff)
//            {
//                Fhandler->ar_full_buffers.push(Fhandler->Ar_currentReadBuff);
//                Fhandler->Ar_currentReadBuff=0;
//            }
//            while(!Fhandler->ar_full_buffers.empty())
//            {
//                Fhandler->ar_empty_buffers.push(Fhandler->ar_full_buffers.front());
//                Fhandler->ar_full_buffers.pop();
//            }
//            pthread_mutex_unlock(&(Fhandler->m_buff_upd));
//
//            Fhandler->reset_wait = false;
//
//
//            pthread_barrier_wait(&(Fhandler->finalize_barrier));
//        }


    }
    //cout << "k" << flush;
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
       tmp= Fhandler->write_full_buffers.front();
       Fhandler->write_full_buffers.pop();
       delete []tmp->buff;
       delete tmp;
    }

    while(!Fhandler->write_empty_buffers.empty())
    {
       tmp= Fhandler->write_empty_buffers.front();
       Fhandler->write_empty_buffers.pop();
       delete []tmp->buff;
       delete tmp;
    }
    }


    //

    //pthread_exit(NULL);
        fclose(fp_Y);
        fclose(fp_Ar);
        fclose(fp_B);
        fclose(fp_R);
        fclose(fp_SD2);
        fclose(fp_P);

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

void AIOwrapper::getCurrentWriteBuffers(type_precision* &B,type_precision* &R,type_precision* &SD2,type_precision* &P)
{
    B = &(Fhandler->currentWriteBuff->buff[0]);
    R = &(Fhandler->currentWriteBuff->buff[Fhandler->max_b_blockSize*Fhandler->p]);
    SD2 = &(Fhandler->currentWriteBuff->buff[Fhandler->max_b_blockSize*(Fhandler->p+1)]);
    P = &(Fhandler->currentWriteBuff->buff[Fhandler->max_b_blockSize*(Fhandler->p+2)]);
}

void AIOwrapper::write_OutFiles(type_precision* &B,type_precision* &R,type_precision* &SD2,type_precision* &P,  int blockSize)
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


    Fhandler->write_full_buffers.push(Fhandler->currentWriteBuff);
    Fhandler->b_blockSize = blockSize;


    Fhandler->currentWriteBuff = Fhandler->write_empty_buffers.front();
    Fhandler->write_empty_buffers.pop();

    B = &(Fhandler->currentWriteBuff->buff[0]);
    R = &(Fhandler->currentWriteBuff->buff[Fhandler->b_blockSize*Fhandler->p]);
    SD2 = &(Fhandler->currentWriteBuff->buff[Fhandler->b_blockSize*(Fhandler->p+1)]);
    P = &(Fhandler->currentWriteBuff->buff[Fhandler->b_blockSize*(Fhandler->p+2)]);


    pthread_mutex_unlock(&(Fhandler->m_buff_upd));


    pthread_mutex_lock(&(Fhandler->m_more));
    pthread_cond_signal( &(Fhandler->condition_more ));
    pthread_mutex_unlock(&(Fhandler->m_more));
}





void AIOwrapper::prepare_OutFiles(int max_b_blockSize, int p)
{

    Fhandler->max_b_blockSize = max_b_blockSize;
    Fhandler->p=p;
    int buff_count = 4;

    type_buffElement* tmp;


    for(int i = 0; i< buff_count  ; i++)
    {
        tmp = new type_buffElement();
        tmp->buff = new type_precision[Fhandler->max_b_blockSize*(2*Fhandler->p+2)];
        tmp->size = max_b_blockSize;
        Fhandler->write_empty_buffers.push(tmp);
    }
    Fhandler->currentWriteBuff = Fhandler->write_empty_buffers.front();
    Fhandler->write_empty_buffers.pop();


}


 void AIOwrapper::write_significantValues(int Y, int X_R, float R, float SD2, float P)
 {

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
    //void *status;


    //cout << "ra" << flush;

//    Fhandler->reset_wait = true;
//    pthread_barrier_wait(&(Fhandler->finalize_barrier));
//
////    pthread_mutex_lock(&(Fhandler->m_buff_upd));
////    Fhandler->Ar_to_readSize = Fhandler->Ar_Amount;
////
////    if(Fhandler->Ar_currentReadBuff)
////    {
////        Fhandler->ar_full_buffers.push(Fhandler->Ar_currentReadBuff);
////        Fhandler->Ar_currentReadBuff=0;
////    }
////
////    while(!Fhandler->ar_full_buffers.empty())
////    {
////        Fhandler->ar_empty_buffers.push(Fhandler->ar_full_buffers.front());
//////        for( int i = 0; i < Fhandler->n*Fhandler->r*Fhandler->Ar_blockSize; i++)
//////        {
//////            ((Fhandler->ar_full_buffers.front())->buff)[i] = 0;
//////        }
////        Fhandler->ar_full_buffers.pop();
////    }
////    pthread_mutex_unlock(&(Fhandler->m_buff_upd));
////
////    Fhandler->reset_wait = false;
//
//    pthread_barrier_wait(&(Fhandler->finalize_barrier));
//
//    pthread_mutex_lock(&(Fhandler->m_more));
//    pthread_cond_signal( &(Fhandler->condition_more ));
//    pthread_mutex_unlock(&(Fhandler->m_more));


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
                file_pos = i*Fhandler->n+ it->first;
                fseek ( fp , file_pos*sizeof(type_precision) , SEEK_SET );
                chunk_size_buff = it->second;

                size_t result = fread (&(Fhandler->AL[buff_pos]),sizeof(type_precision),chunk_size_buff,fp); result++;
                buff_pos += chunk_size_buff;
            }
        }


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
