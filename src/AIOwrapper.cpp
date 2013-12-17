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

    pthread_barrier_init(&(FHandler.finalize_barrier), NULL, 2);


    Fhandler->fakefiles = params.use_fake_files;

    databel_fvi* Yfvi;
    databel_fvi* ALfvi;
    databel_fvi* ARfvi;

    if (!Fhandler->fakefiles)
    {
        Fhandler->fnameAL   = params.fnameAL;
        Fhandler->fnameAR   = params.fnameAR;
        Fhandler->fnameY    = params.fnameY;
        Fhandler->fnameOutB = params.fnameOutB;

        databel_fvi* Yfvi  = load_databel_fvi( (Fhandler->fnameY+".fvi").c_str() );
        databel_fvi* ALfvi = load_databel_fvi( (Fhandler->fnameAL+".fvi").c_str() );
        databel_fvi* ARfvi = load_databel_fvi( (Fhandler->fnameAR+".fvi").c_str() );
        params.n = ALfvi->fvi_header.numObservations;
        params.m = ARfvi->fvi_header.numVariables / params.r;
        params.t = Yfvi->fvi_header.numVariables;
        params.l = ALfvi->fvi_header.numVariables;

        //block size to keep mem under 1 gigabyte
        int opt_block = params.n / (4*1000^3) * (1/(2*params.r));
        int opt_tb = 1000;
        int opt_mb = 1000;

        params.mb = min(params.m, opt_tb);
        params.tb = min(params.t, opt_mb);
    }
    else
    {
    }

    params.p = params.l + params.r;


    //block size to keep mem under 1 gigabyte
//    int opt_block = params.n/(4*1000^3)*(1/(2*params.r));
//    int opt_tb = max(4*2000, opt_block);
//    int opt_mb = max(2000, opt_block);
//
//    params.mb = min(params.m, opt_tb);
//    params.tb = min(params.t, opt_mb);

    prepare_B();
    prepare_AL(params.l, params.n);
    prepare_AR(params.mb, params.n, params.m, params.r);
    prepare_Y(params.tb, params.n, params.t);
}


void AIOwrapper::finalize()
{
    //cout << "f";
    void *status;

    Fhandler->not_done = false;
    pthread_mutex_lock(&(Fhandler->m_more));
    pthread_cond_signal( &(Fhandler->condition_more ));
    pthread_mutex_unlock(&(Fhandler->m_more));

    pthread_barrier_wait(&(Fhandler->finalize_barrier));
    finalize_Y();
    finalize_AR();
    finalize_AL();
    finalize_B();

    pthread_attr_destroy(&(Fhandler->attr));

    pthread_mutex_destroy(&(Fhandler->m_more));
    pthread_cond_destroy(&(Fhandler->condition_more));

    pthread_mutex_destroy(&(Fhandler->m_read));
    pthread_cond_destroy(&(Fhandler->condition_read));
}


void AIOwrapper::prepare_B()
{
    Fhandler->fp_B = fopen((Fhandler->fnameOutB + ".fvd").c_str(), "w+b");
    if (Fhandler->fp_B == 0)
    {
        cout << "Error Opening File B " << Fhandler->fnameOutB << endl;
        exit(1);
    }
}


void AIOwrapper::finalize_B()
{
    fclose(Fhandler->fp_B);
}


void* AIOwrapper::async_io( void *ptr )
{
    //cout << "async_io\n" << flush;
    type_fileh* Fhandler = (type_fileh *)ptr;
    int size_buff, tmp_y_blockSize, tmp_ar_blockSize;

    struct timespec timeToWait;
    FILE*  fp_Y;

    FILE*  fp_Ar;
    if (!Fhandler->fakefiles)
    {
        fp_Y = fopen((Fhandler->fnameY+".fvd").c_str(), "rb");
        if (fp_Y == 0)
        {
            cout << "Error Reading File Y " << Fhandler->fnameY << endl;
            exit(1);
        }

        fp_Ar = fopen((Fhandler->fnameAR+".fvd").c_str(), "rb");
        if (fp_Ar == 0)
        {
            cout << "Error Reading File Xr " << Fhandler->fnameAR << endl;
            exit(1);
        }
    }
    else
    {
//        fp_Y = fopen("tempY.bin", "w+b");
//        if (fp_Y == 0)
//        {
//            cout << "Error creating temp File Y " << Fhandler->fnameY << endl;
//            exit(1);
//        }
//        fwrite(Fhandler->Yb, sizeof(type_precision), Fhandler->n*Fhandler->y_blockSize, fp_Y);
//        fclose(fp_Y);
//
//
//        fp_Ar = fopen("tempAR.bin", "w+b");
//        if (fp_Ar == 0)
//        {
//            cout << "Error creating temp File Y " << Fhandler->fnameY << endl;
//            exit(1);
//        }
//        fwrite(Fhandler->Yb, sizeof(type_precision), Fhandler->n*Fhandler->y_blockSize, fp_Ar);
//        fclose(fp_Ar);
//
//        fp_Y = fopen("tempY.bin", "rb");
//        if (fp_Y == 0)
//        {
//            cout << "Error Reading File Y " << Fhandler->fnameY << endl;
//            exit(1);
//        }
//
//        fp_Ar = fopen("tempAR.bin", "rb");
//        if (fp_Ar == 0)
//        {
//            cout << "Error Reading File Xr " << Fhandler->fnameAR << endl;
//            exit(1);
//        }

    }


    Fhandler->not_done = true;
    Fhandler->reset_wait = false;

    //cout << "c" << flush;

    while (Fhandler->not_done)
    {
        while (!Fhandler->empty_buffers.empty() && Fhandler->y_to_readSize)
        {
            tmp_y_blockSize = Fhandler->y_blockSize;
            if (Fhandler->y_to_readSize < Fhandler->y_blockSize)
                tmp_y_blockSize = Fhandler->y_to_readSize;

            Fhandler->y_to_readSize -= tmp_y_blockSize;
            size_buff = Fhandler->n * tmp_y_blockSize;
            //cout << Fhandler->y_to_readSize << endl;

            pthread_mutex_lock(&(Fhandler->m_buff_upd));
            //cout << " pre;" << Fhandler->full_buffers.size() << ";" << Fhandler->empty_buffers.size() << endl;
            type_buffElement* tobeFilled = Fhandler->empty_buffers.front();
            Fhandler->empty_buffers.pop();
            //pthread_mutex_unlock(&(Fhandler->m_buff_upd));

            tobeFilled->size = tmp_y_blockSize;

            if (Fhandler->fakefiles)
            {
                //fread (tobeFilled->buff, sizeof(type_precision), size, fp_Y);
                //fseek ( fp_Y , 0 , SEEK_SET );
                int old_seed = Fhandler->seed;
                srand (old_seed);
                re_random_vec(tobeFilled->buff, size_buff );
                re_random_vec_nan(tobeFilled->buff, size_buff );
                Fhandler->seed += 75;
            }
            else
            {
                fread (tobeFilled->buff, sizeof(type_precision), size_buff, fp_Y);
                if (Fhandler->y_to_readSize <= 0)
                {
                    fseek ( fp_Y , 0 , SEEK_SET );
                }
            }


            //pthread_mutex_lock(&(Fhandler->m_buff_upd));
            Fhandler->full_buffers.push(tobeFilled);
            //  cout << "\nStoring " << tobeFilled << endl;
            //cout << " post;" << Fhandler->full_buffers.size() << ";" << Fhandler->empty_buffers.size() << endl;
            pthread_mutex_unlock(&(Fhandler->m_buff_upd));

            pthread_mutex_lock(&(Fhandler->m_read));
            pthread_cond_signal( &(Fhandler->condition_read ));
            pthread_mutex_unlock(&(Fhandler->m_read));
        }

        while (!Fhandler->ar_empty_buffers.empty() && Fhandler->Ar_to_readSize)
        {
            tmp_ar_blockSize = Fhandler->Ar_blockSize;
            if (Fhandler->Ar_to_readSize < Fhandler->Ar_blockSize)
                tmp_ar_blockSize = Fhandler->Ar_to_readSize;

            Fhandler->Ar_to_readSize -= tmp_ar_blockSize;
            size_buff = Fhandler->n * tmp_ar_blockSize*Fhandler->r;

            pthread_mutex_lock(&(Fhandler->m_buff_upd));
            type_buffElement* tobeFilled = Fhandler->ar_empty_buffers.front();
            Fhandler->ar_empty_buffers.pop();


            tobeFilled->size = tmp_ar_blockSize;

            if (Fhandler->fakefiles)
            {
                //fread (tobeFilled->buff, sizeof(type_precision), size, fp_Ar);
                //fseek ( fp_Ar , 0 , SEEK_SET );
                re_random_vec(tobeFilled->buff,
                              Fhandler->n * tmp_ar_blockSize*Fhandler->r );
                re_random_vec_nan(tobeFilled->buff,
                                  Fhandler->n * tmp_ar_blockSize*Fhandler->r );
            }
            else
            {
                fread (tobeFilled->buff, sizeof(type_precision), size_buff, fp_Ar);
                if (Fhandler->Ar_to_readSize <= 0)
                {
                    fseek ( fp_Ar , 0 , SEEK_SET );
                }
            }

            Fhandler->ar_full_buffers.push(tobeFilled);
            //  cout << "\nStoring " << tobeFilled << endl;
            pthread_mutex_unlock(&(Fhandler->m_buff_upd));

            pthread_mutex_lock(&(Fhandler->m_read));
            pthread_cond_signal( &(Fhandler->condition_read ));
            pthread_mutex_unlock(&(Fhandler->m_read));
        }




#ifdef WINDOWS
        SYSTEMTIME time;
        GetSystemTime(&time);

        timeToWait.tv_sec = time.wSecond + 500/1000;
        long int morenanos = (500%1000) * 1000000;
        timeToWait.tv_nsec = time.wMilliseconds * 1000 + morenanos ;
#else
        clock_gettime(CLOCK_REALTIME, &timeToWait);
        timeToWait.tv_nsec += 150;
#endif

        pthread_mutex_lock(&(Fhandler->m_more));
        pthread_cond_timedwait( &(Fhandler->condition_more),
                                &(Fhandler->m_more), &timeToWait );
        pthread_mutex_unlock( &(Fhandler->m_more ));

        pthread_mutex_lock( &(Fhandler->m_read) );
        pthread_cond_signal( &(Fhandler->condition_read ));
        pthread_mutex_unlock( &(Fhandler->m_read) );

        if (Fhandler->reset_wait)
        {
            pthread_barrier_wait(&(Fhandler->finalize_barrier));
            //wait for main thread to reset everything
            pthread_barrier_wait(&(Fhandler->finalize_barrier));
        }
    }
    //cout << "k" << flush;
    //barrier
    pthread_barrier_wait(&(Fhandler->finalize_barrier));

    type_buffElement* tmp;
    while (!Fhandler->full_buffers.empty())
    {
       tmp = Fhandler->full_buffers.front();
       Fhandler->full_buffers.pop();
       delete[] tmp->buff;
       delete tmp;
    }

    while (!Fhandler->empty_buffers.empty())
    {
       tmp = Fhandler->empty_buffers.front();
       Fhandler->empty_buffers.pop();
       delete[] tmp->buff;
       delete tmp;
    }

    while (!Fhandler->ar_full_buffers.empty())
    {
       tmp = Fhandler->ar_full_buffers.front();
       Fhandler->ar_full_buffers.pop();
       delete tmp->buff;
       delete tmp;
    }

    while (!Fhandler->ar_empty_buffers.empty())
    {
       tmp = Fhandler->ar_empty_buffers.front();
       Fhandler->ar_empty_buffers.pop();
       delete[] tmp->buff;
       delete tmp;
    }

    //

    //pthread_exit(NULL);

    if (!Fhandler->fakefiles)
    {
        fclose(fp_Y);
        fclose(fp_Ar);
    }
    else
    {
        fclose(fp_Y);
        fclose(fp_Ar);
    }
//
//            //!induce realistic fileread delay
}


void AIOwrapper::load_ARblock(type_precision** Ar, int &Ar_blockSize)
{
    int status;
    int createstatus = 0;
    //cout<<"^";

    while (Fhandler->ar_full_buffers.empty())
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
    if (Fhandler->Ar_currentReadBuff)
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
//    if (Fhandler->Ar_to_readSize < Fhandler->Ar_blockSize)
//        tmp_ar_blockSize = Fhandler->Ar_to_readSize;
//
//    Fhandler->Ar_to_readSize -= tmp_ar_blockSize;
//    int size_buff = Fhandler->n * tmp_ar_blockSize * Fhandler->r;
//    re_random_vec(Fhandler->Ar, size_buff);
//Ar_blockSize=tmp_ar_blockSize;

    (*Ar) = Fhandler->Ar;
}


void AIOwrapper::load_Yblock(type_precision** Y, int &y_blockSize)
{
    int status;
    int createstatus = 0;

    while (Fhandler->full_buffers.empty())
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
    //cout << " pre," << Fhandler->full_buffers.size() << ";" << Fhandler->empty_buffers.size() << endl;

    if (Fhandler->currentReadBuff)
    {
        //memset(Fhandler->currentReadBuff->buff, 0, y_blockSize);
        Fhandler->empty_buffers.push(Fhandler->currentReadBuff);
    }
    Fhandler->currentReadBuff = Fhandler->full_buffers.front();
    Fhandler->full_buffers.pop();

    //cout << "\nReading " << Fhandler->currentReadBuff << endl;
    Fhandler->Yb = Fhandler->currentReadBuff->buff;
    y_blockSize = Fhandler->currentReadBuff->size;

    (*Y) = Fhandler->Yb;

     //cout << " post," << Fhandler->full_buffers.size() << ";" << Fhandler->empty_buffers.size() << endl;

    pthread_mutex_unlock(&(Fhandler->m_buff_upd));

    pthread_mutex_lock(&(Fhandler->m_more));
    pthread_cond_signal( &(Fhandler->condition_more ));
    pthread_mutex_unlock(&(Fhandler->m_more));

    //matlab_print_matrix("Y", Fhandler->n, y_blockSize,*Y);
}


void AIOwrapper::write_B(type_precision* B, int p, int blockSize)
{
    fwrite(B , sizeof(type_precision), p*blockSize, Fhandler->fp_B);
}


void AIOwrapper::prepare_Y(int y_blockSize, int n, int totalY)
{
    //for fake files

    Fhandler->seed = 1337;
    srand(Fhandler->seed);

    Fhandler->y_blockSize = y_blockSize;

    Fhandler->n = n;
    Fhandler->Y_Amount = totalY;
    Fhandler->y_to_readSize = Fhandler->Y_Amount;
    Fhandler->buff_count = min(2,(totalY + y_blockSize - 1)/y_blockSize) ;
    //cout << "buffcount " << Fhandler->buff_count;


    Fhandler->currentReadBuff = 0;
    type_buffElement* tmp;

    for (int i = 0; i < Fhandler->buff_count; i++)
    {
        tmp = new type_buffElement();
        tmp->buff = new type_precision[Fhandler->n * Fhandler->y_blockSize];
        tmp->size = y_blockSize;
//        for ( int i = 0; i < Fhandler->n*Fhandler->y_blockSize; i++)
//        {
//            (tmp->buff)[i] = 0;
//        }
        Fhandler->empty_buffers.push(tmp);
    }
    Fhandler->Yb = tmp->buff;



    pthread_mutex_init(&(Fhandler->m_buff_upd), NULL);
    pthread_mutex_init(&(Fhandler->m_more), NULL);
    pthread_mutex_init(&(Fhandler->m_read), NULL);

    pthread_attr_init(&(Fhandler->attr));
    pthread_attr_setdetachstate(&(Fhandler->attr), PTHREAD_CREATE_JOINABLE);

    pthread_create( &(Fhandler->iothread), &(Fhandler->attr),
                    AIOwrapper::async_io, (void*)Fhandler);
}


void AIOwrapper::reset_Y()
{
    void *status;

    Fhandler->seed = 1337;

    //cout << "ry" << flush;

    Fhandler->reset_wait = true;
    pthread_barrier_wait(&(Fhandler->finalize_barrier));

    pthread_mutex_lock(&(Fhandler->m_buff_upd));
    Fhandler->y_to_readSize = Fhandler->Y_Amount;

    if (Fhandler->currentReadBuff)
    {
        Fhandler->full_buffers.push(Fhandler->currentReadBuff);
        Fhandler->currentReadBuff=0;
    }

    while (!Fhandler->full_buffers.empty())
    {
        Fhandler->empty_buffers.push(Fhandler->full_buffers.front());
        for ( int i = 0; i < Fhandler->n*Fhandler->y_blockSize; i++)
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
    void *status;


    //cout << "ra" << flush;

    Fhandler->reset_wait = true;
    pthread_barrier_wait(&(Fhandler->finalize_barrier));

    pthread_mutex_lock(&(Fhandler->m_buff_upd));
    Fhandler->Ar_to_readSize = Fhandler->Ar_Amount;

    if (Fhandler->Ar_currentReadBuff)
    {
        Fhandler->ar_full_buffers.push(Fhandler->Ar_currentReadBuff);
        Fhandler->Ar_currentReadBuff=0;
    }

    while (!Fhandler->ar_full_buffers.empty())
    {
        Fhandler->ar_empty_buffers.push(Fhandler->ar_full_buffers.front());
//        for ( int i = 0; i < Fhandler->n*Fhandler->r*Fhandler->Ar_blockSize; i++)
//        {
//            ((Fhandler->ar_full_buffers.front())->buff)[i] = 0;
//        }
        Fhandler->ar_full_buffers.pop();
    }
    pthread_mutex_unlock(&(Fhandler->m_buff_upd));

    Fhandler->reset_wait = false;

    pthread_barrier_wait(&(Fhandler->finalize_barrier));

    pthread_mutex_lock(&(Fhandler->m_more));
    pthread_cond_signal( &(Fhandler->condition_more ));
    pthread_mutex_unlock(&(Fhandler->m_more));
}


void AIOwrapper::finalize_Y()
{
}


void AIOwrapper::prepare_AR( int desired_blockSize, int n, int totalR, int columnsAR)
{
    Fhandler->Ar = new type_precision[desired_blockSize * columnsAR * n];
    Fhandler->Ar_blockSize = desired_blockSize;
    Fhandler->r = columnsAR;
    Fhandler->Ar_Amount = totalR;
    Fhandler->Ar_to_readSize = Fhandler->Ar_Amount;

    int buff_count = min(2, (totalR+ desired_blockSize - 1)/desired_blockSize);

    Fhandler->Ar_currentReadBuff = 0;
    type_buffElement* tmp;

    for (int i = 0; i < buff_count  ; i++)
    {
        tmp = new type_buffElement();
        tmp->buff = new type_precision[n * desired_blockSize * columnsAR];
        tmp->size = desired_blockSize;
//        memset(tmp->buff,-9, n*desired_blockSize*columnsAR);

        Fhandler->ar_empty_buffers.push(tmp);
    }
    Fhandler->Ar = tmp->buff;
}


void AIOwrapper::finalize_AR()
{
}


void AIOwrapper::load_AL(type_precision** AL)
{
    if (Fhandler->fakefiles)
    {
        re_random_vec(Fhandler->AL, Fhandler->n*Fhandler->l);
        re_random_vec_nan(Fhandler->AL, Fhandler->n*Fhandler->l);
        (*AL) = Fhandler->AL;
    }
    else
    {
        FILE *fp;
        fp = fopen((Fhandler->fnameAL+".fvd").c_str(), "rb");
        if (fp == 0)
        {
            cout << "Error Reading File " << Fhandler->fnameAL << endl;
            exit(1);
        }

        fread (Fhandler->AL, sizeof(type_precision), Fhandler->l*Fhandler->n, fp);
        fclose(fp);
    }



    (*AL) = Fhandler->AL;

    //matlab_print_matrix("AL", Fhandler->n, Fhandler->l,*AL);
}


void AIOwrapper::prepare_AL( int columnsAL, int n)
{
    Fhandler->AL = new type_precision[columnsAL * n];
    Fhandler->l = columnsAL;
}


void AIOwrapper::finalize_AL()
{
    delete []Fhandler->AL;
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
        char err[100];

        f = fopen(path, mode );
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

    if ( (buf = malloc(size)) == NULL ) {
        cout << "\nCouldn't allocate %ld bytes of memory in %s:%ld\n";
        exit(1);
    }
    return buf;
}


struct databel_fvi * AIOwrapper::load_databel_fvi( const char *path )
{
        FILE *f;
        databel_fvi *fvi;
        size_t data_size;

        f = fgls_fopen(path, "r");

        fvi = (databel_fvi*) fgls_malloc( sizeof(databel_fvi) );
        // Header
        fread( &fvi->fvi_header, sizeof(databel_fvi_header), 1, f );
        // Labels
        data_size = (fvi->fvi_header.numVariables + fvi->fvi_header.numObservations) *
            fvi->fvi_header.namelength * sizeof(char);
        fvi->fvi_data = (char *) fgls_malloc(data_size);
        // Load labels
        fread(fvi->fvi_data, 1, data_size, f);

        fclose(f);

        return fvi;
}
