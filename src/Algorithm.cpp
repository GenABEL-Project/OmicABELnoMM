#include "Algorithm.h"

Algorithm::Algorithm()
{
    //ctor
}

Algorithm::~Algorithm()
{
    //dtor
}


void Algorithm::solve(struct Settings params, struct Outputs &out, int type)
{
    switch (type)
    {
        case FULL_NEQ:
                fullNEQ(params,out);
            break;
        case P_NEQ:
                partialNEQ(params,out);
            break;
        case P_NEQ_B_OPT:
                partialNEQ_Blocked_STL(params,out);
            break;
        case FULL_QR:
                fullQR(params,out);
            break;
        case P_QR:
                partialQR(params,out);
            break;
        case P_QR_B_OPT:
                partialQR_Blocked_Rtl(params,out);
            break;
        case P_NEQ_B_OPT_MD:
                partialNEQ_Blocked_STL_MD(params,out);
            break;

        default:
            break;
    }

}

void Algorithm::extract_subMatrix(type_precision* source, type_precision* dest,int dim1_source, int dim2_source,
                                                            int dim1_ini,int dim1_end,int dim2_ini,int dim2_end)
{

    int i,j,idx=0;
    int size, source_ini;
    for(i = dim2_ini; i<dim2_end; i++)
    {
        j = dim1_ini;
        source_ini = i*dim1_source+j;
        size = dim1_end-dim1_ini;
        memcpy( (type_precision*)&dest[idx], (type_precision*)&source[source_ini], size * sizeof(type_precision) );
//        for(j = dim1_ini; j<dim1_end; j++)
//        {
//            dest[idx] = source[i*dim1_source+j];
//            idx++;
//        }
        idx += size;
    }

}

void Algorithm::prepare_Bfinal(type_precision* bfinal, type_precision* top,type_precision* bot,int dim1_b, int dim2_b,int dim1_b_bot)
{
    //memcpy are faster version of the fors
    int i,k,w,top_idx,bot_idx,max = dim1_b*dim2_b;
    int size;
    top_idx = 0;
    bot_idx = 0;
    for(k = 0; k < dim2_b; k++)
    {
        size = k*dim1_b+(dim1_b-dim1_b_bot)-(k*dim1_b);
        memcpy( (type_precision*)&bfinal[k*dim1_b], (type_precision*)&top[top_idx], size * sizeof(type_precision) );
//        for(i = k*dim1_b; i < k*dim1_b+(dim1_b-dim1_b_bot); i++)
//        {
//            bfinal[i] = top[top_idx];
//            top_idx++;
//        }
        top_idx += size;
        i = k*dim1_b + size;
        w=i;

        size = w+dim1_b_bot - w;
        memcpy( (type_precision*)&bfinal[w], (type_precision*)&bot[bot_idx], size * sizeof(type_precision) );
//        for(i = w; i < w+dim1_b_bot; i++)
//        {
//            bfinal[i] = bot[bot_idx];
//            bot_idx++;
//        }
        bot_idx += size;
    }

}

void Algorithm::prepare_QY(type_precision* qy, type_precision* top,type_precision* bot,int dim1_QY, int dim2_QY,int dim1_qy_bot,int bot_blocks )
{

    int i,k,w,top_idx,bot_idx,max = dim1_QY*dim2_QY;
    top_idx = 0;
    bot_idx = 0;
    for(k = 0; k < dim2_QY; k++)
    {
        for(i = k*dim1_QY; i < (k+1)*dim1_QY-dim1_qy_bot; i++)
        {
            qy[i] = top[top_idx];
            top_idx++;
        }
        w=i;

        for(i = w; i < w+dim1_qy_bot; i++)
        {
            qy[i] = bot[bot_idx];
            bot_idx++;
        }
        bot_idx+=(bot_blocks-1)*dim1_qy_bot;
    }

}

type_precision* Algorithm::extract_R(type_precision* A,int dim1_A, int dim2_A)
{
    type_precision* R = (type_precision*)calloc(dim2_A*dim2_A,sizeof(type_precision));
    int i,j;

    int R_idx=0;

    for(i = 0; i < dim2_A; i++)
    {
        for(j = 0; j <= i; j++)
        {
            R[R_idx] = A[j+i*dim1_A];
            R_idx++;
        }
        R_idx = dim2_A*(i+1);
    }
    return R;
}

type_precision* Algorithm::prepare_R(type_precision* RL,int dim1_A, int dim2_AL,int dim2_AR)
{
    int R_dim = (dim2_AR+dim2_AL);
    type_precision* R = new type_precision[R_dim*R_dim];

    int i,j;

    int RL_idx=0;
    int R_idx=0;

    for(i = 0; i < dim2_AL; i++)
    {
        for(j = 0; j <= i; j++)
        {
            RL_idx = i*dim1_A + j;
            R_idx = i*R_dim + j;
            R[R_idx] = RL[RL_idx];

        }

    }
    return R;
}

void Algorithm::update_R(type_precision* R, type_precision* topRr, type_precision* botRr,int dim1, int dim2, int r)
{
    int i,j,w;
    int max = dim1*dim2;
    int rtr_idx=0;
    int rbr_idx=0;
    for( j = r; j > 0 ; j--)
    {
        for(i = max-dim1*j; i < max-dim1*j+dim2-r; i++)
        {
            R[i] = topRr[rtr_idx];
            rtr_idx++;
        }
        w = i;
        for(i = w; i < w+r; i++)
        {
            R[i] = botRr[rbr_idx];
            rbr_idx++;
        }

    }

}


void Algorithm::build_S(type_precision* S,type_precision* Stl,type_precision* Str,type_precision* Sbr,int l,int r)
{
    int Sidx;
    int p = l+r;
    for(int i = 0; i < l; i++)
    {
        Sidx = i*p;
        for(int j= 0; j <= i; j++)
        {
            S[Sidx] = Stl[j+i*l];
            Sidx++;
        }
    }
    for(int i = 0; i < r; i++)
    {
        Sidx = l*p+p*i;
        for(int j= 0; j < l; j++)
        {
            S[Sidx] = Str[j+i*l];
            Sidx++;
        }
    }

    for(int i = 0; i < r; i++)
    {
        Sidx = l*p+l+p*i;
        for(int j= 0; j <= i; j++)
        {
            S[Sidx] = Sbr[j+i*r];
            Sidx++;
        }
    }



}

void Algorithm::check_result(type_precision* AL, type_precision* AR,int rowsA,int colsA, int rhs,int colsAR,
                                                    type_precision* y, type_precision* res)
{
    type_precision* A = (type_precision*)malloc(rowsA*colsA*sizeof(type_precision));

    int i,ar_idx=0;
    for(i = 0; i <rowsA*(colsA-colsAR) ; i++)
    {
        A[i] =  AL[i];
    }

    for(i = rowsA*(colsA-colsAR); i <rowsA*colsA ; i++)
    {
        A[i] =  AR[ar_idx];
        ar_idx++;
    }



    type_precision* ynew = replicate_vec(y,rowsA*rhs);
    type_precision* new_sol = (type_precision*)malloc(colsA*rhs*sizeof(type_precision));

    lapack_int info = LAPACKE_sgels(STORAGE_TYPE,'N',rowsA,colsA,rhs,A,rowsA,ynew,rowsA);
    assert(info == 0,"Error Check");





    int index=0;
    int index_new = 0;
    for(i = 0; i < rhs; i++)
    {
         copy_vec( &ynew[index], &new_sol[index_new], colsA);
         index += rowsA;
         index_new += colsA;
    }


//    if(PRINT)
//        printf("\n Btop=(Rtl\\(Ql'*Y))-(Rtl\\Rtr)*(Rbr\\(Qr'*Y)); \n [Btop ; Rbr\\Qr'*Y] - bcomputed \n");

    cblas_saxpy(rhs*colsA, -1.0, res, 1,new_sol,1);
    type_precision u_norm=cblas_snrm2(rhs*colsA,new_sol,1);
    //
    if(abs(u_norm) >= 0.0001 || isnan(u_norm))
    {


        fflush(stdout);
        matlab_print_matrix("AL",rowsA,colsA-colsAR,AL);
        matlab_print_matrix("AR",rowsA,colsAR,AR);
        matlab_print_matrix("Y",rowsA,rhs,y);
        printf("\nA = [AL AR]; [Q,R] = qr(A,0); rr=R\\(Q'*Y)\n");
        matlab_print_matrix("bcomputed",colsA,rhs,res);
        matlab_print_matrix("newsol",colsA,rhs,ynew);
        printf("\n%%\tnrom: %0.2g", u_norm);
         exit(1);
    }
    else
    {
        matlab_print_matrix("bcomputed",colsA,rhs,res);
        matlab_print_matrix("newsol",colsA,rhs,ynew);
        printf("\n%%\tnrom: %0.2g", u_norm);
    }



    //cout << "\t**************";


    free(ynew);
    free(new_sol);
    free(A);






}

///////////////////////////////



void Algorithm::partialNEQ_Blocked_STL_MD(struct Settings params, struct Outputs &out)
{
    int max_threads = params.threads;


    srand (time(NULL));

    blas_set_num_threads(max_threads);


    type_precision *Ytemp;
    lapack_int info,n,lda,ldy,l,r,k,p;

    int i,j,w;
    int m;

    cputime_type start_tick,start_tick2, end_tick;

    AIOfile.initialize(params);

    n = params.n; l=params.l; r=params.r; p = l+r;

    int y_amount = params.t;
    int y_block_size = params.tb;//kk

    int a_amount = params.m;
    int a_block_size = params.mb;

    int a_iters = (a_amount + a_block_size - 1) / a_block_size;

    int y_iters = (y_amount + y_block_size - 1) / y_block_size;


    lda = n; ldy = n;
    k = p;

    for(int j = 0; j < y_iters; j++)
    {

        if(y_iters >= 40 && (i%(y_iters/40))==0)
        {
            cout << "*" << flush;

        }
        for(int i = 0; i < a_iters; i++)
        {
            for(int jj = 0; jj < y_block_size; jj++)
            {
                if(y_iters < 40 && (y_block_size < 3 || (jj%(y_block_size/3))==0))
                {
                    cout << "*" << flush;

                }
            }
        }


    }
    cout << endl;


    type_precision Stl[l*l];
    type_precision Str[l*r*a_block_size];


    type_precision* Ay_top = new type_precision[l*y_amount];
    type_precision* Ay_bot = new type_precision[y_block_size*a_block_size*r];

    type_precision* AR = new type_precision[n*r*a_block_size*1];
    type_precision* AL = new type_precision[n*l*1];
    type_precision* backupAR;// = new type_precision[n*r*a_block_size];
    type_precision* backupAL;// = new type_precision[n*l];


    AIOfile.load_AL(&backupAL);
    int total_al_nans = replace_nans(0,backupAL, n,l);
    copy_vec(backupAL,AL,n*l);


    type_precision* Y;

    //printf("\n\n%%Computations\n%%");


    get_ticks(start_tick);

    for(int j = 0; j < y_iters; j++)
    {

        if(y_iters >= 40 && (i%(y_iters/40))==0)
        {
            cout << AIOfile.io_overhead << flush;
            AIOfile.io_overhead = "*";
        }

        AIOfile.load_Yblock(&Y,y_block_size);

        list<long int> y_nan_idxs[y_block_size];

        int total_y_nans = replace_nans(&y_nan_idxs[0],Y,n,y_block_size);


        //!Ay_top=AL'*Y
        cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,
                                    l,y_block_size,n,1.0,AL,n,Y,n,0.0,&Ay_top[j*l*y_block_size],l);



        for(int i = 0; i < a_iters; i++)
        {

            AIOfile.load_ARblock(&backupAR, a_block_size);
            int total_ar_nans = replace_nans(0,backupAR, n,a_block_size*r);
            copy_vec(backupAR,AR,n*r*a_block_size);

            get_ticks(start_tick2);

            get_ticks(end_tick);
            out.acc_pre += ticks2sec(end_tick-start_tick2,cpu_freq);

            get_ticks(start_tick2);
            //!Ay_bot=AR'*Y
            cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,
                                r*a_block_size,y_block_size,n,1.0,AR,n,Y,n,0.0,Ay_bot,r*a_block_size);

            get_ticks(end_tick);
            out.firstloop += ticks2sec(end_tick-start_tick2,cpu_freq);

            #pragma omp parallel default(shared)
            {

                for(int jj = 0; jj < y_block_size; jj++)
                {
                    int thread_id = 0*omp_get_thread_num();
                    int aL_idx = thread_id*l*n;
                    int aR_idx = thread_id*r*n*a_block_size;

                    if(y_iters < 40 && (y_block_size < 3 || (jj%(y_block_size/3))==0))
                    {
                        cout << AIOfile.io_overhead << flush;
                        AIOfile.io_overhead = "*";
                    }

                    copy_vec(backupAL,&AL[aL_idx],n*l);

                    replace_with_zeros(&y_nan_idxs[jj],&AL[aL_idx], n,l,1);
                    //!Generate Stl
                    cblas_ssyrk(CblasColMajor,CblasUpper,CblasTrans,l,n,1.0,&AL[aL_idx],lda,0.0,Stl,l);

                    copy_vec(backupAR,&AR[aR_idx],n*r*a_block_size);
                    replace_with_zeros(&y_nan_idxs[jj],&AR[aR_idx], n,r,a_block_size);

                    get_ticks(start_tick2);
                    //!Generate Str
                    cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,
                                            l,r*a_block_size,n,1.0,&AL[aL_idx],n,&AR[aR_idx],n,0.0,Str,l);
                    get_ticks(end_tick);
                    out.acc_RTL_QLY += ticks2sec(end_tick-start_tick2,cpu_freq);


                    type_precision Sbr[r*r*a_block_size];
                    type_precision Ay[p*a_block_size];

                    #pragma omp for nowait  schedule(dynamic)
                    for(int ii= 0; ii < a_block_size ; ii++)
                    {
                        get_ticks(start_tick2);
                        //!Generate Sbr
                        cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,
                                                r,r,n,1.0,&AR[aR_idx+ii*r*n],n,&AR[aR_idx+ii*r*n],n,0.0,&Sbr[ii*r*r],r);
                        get_ticks(end_tick);
                        out.acc_gemm += ticks2sec(end_tick-start_tick2,cpu_freq);

                        get_ticks(start_tick2);

                        copy_vec(&Ay_top[j*l*y_block_size+jj*l],&Ay[ii*p],l);
                        copy_vec(&Ay_bot[ii*r + jj*r*a_block_size],&Ay[l+ii*p],r);


                        type_precision* B = Ay;
                        type_precision S[p*p];

                        //!Rebuild S
                        build_S(S,Stl,&Str[ii*r*l],&Sbr[ii*r*r],l,r);
                        //matlab_print_matrix("S",p,p,S);

                        get_ticks(end_tick);
                        out.acc_b += ticks2sec(end_tick-start_tick2,cpu_freq);

                        get_ticks(start_tick2);

                        //!b=S\Ay
                        info = LAPACKE_sposv(STORAGE_TYPE,'U',p,1,S,p,&Ay[ii*p],p);

                        get_ticks(end_tick);
                        out.acc_loadxr += ticks2sec(end_tick-start_tick2,cpu_freq);
                        assert(info == 0,"POSV");

                        if( ForceCheck)
                        {
                            #pragma omp critical
                            {
                                check_result(AL,&AR[aR_idx+ii*r*n],n,p,1,r,&Y[jj*n],&Ay[ii*p]);
                            }
                        }
                    }

                    AIOfile.write_B(Ay,p,a_block_size);
                }

            }


        }

         AIOfile.reset_AR();

    }

    get_ticks(end_tick);
    out.duration = ticks2sec(end_tick-start_tick,cpu_freq);
    //out.gflops = a_amount/1000.0*(n*l*r+n*r*r+y_amount*(n*r+p*p*p)))/1000.0/1000.0;
    out.gflops = y_amount*(gemm_flops(l,n,1,0) + a_amount*(gemm_flops(r,n,1,0) +
                                gemm_flops(l,n,l,0)+ gemm_flops(l,n,r,0)+ gemm_flops(r,n,r,0) + (p*p*p/3.0)/1000.0/1000.0/1000.0));


    AIOfile.finalize();

    delete []Ay_top;
    delete []Ay_bot;
    delete []AR;
    delete []AL;
    //delete []backupAL;
    //delete []backupAR;

}


//!*************************************************!//

void Algorithm::partialNEQ_Blocked_STL(struct Settings params, struct Outputs &out)
{
    int max_threads = params.threads;


    srand (time(NULL));

    blas_set_num_threads(max_threads);


    type_precision *Ytemp;
    lapack_int info,n,lda,ldy,l,r,k,p;

    int i,j,w;
    int m;

    cputime_type start_tick,start_tick2, end_tick;

    AIOfile.initialize(params);

    n = params.n; l=params.l; r=params.r; p = l+r;

    int y_amount = params.m;
    int y_block_size = params.mb;//kk

    int a_amount = params.t;
    int a_block_size = params.tb;

    int a_iters = (a_amount + a_block_size - 1) / a_block_size;

    lda = n; ldy = n;
    k = p;



    type_precision* backupAL;

    AIOfile.load_AL(&backupAL);


    type_precision Stl[l*l];
    type_precision Str[l*r*a_block_size];


    type_precision* Ay_top = new type_precision[l*y_amount];
    type_precision* Ay_bot = new type_precision[y_block_size*a_block_size*r];
    type_precision* A = new type_precision[n*p];
    type_precision* AR = new type_precision[n*r*a_block_size];
    type_precision* AL = new type_precision[n*l];

    copy_vec(backupAL,AL,n*l);
    copy_vec(backupAL,A,n*l);


    int y_iters = y_amount/y_block_size;

    type_precision* Y;
    type_precision* backupAR;

    //printf("\n\n%%Computations\n%%");


    get_ticks(start_tick);


    //!Generate S
    cblas_ssyrk(CblasColMajor,CblasUpper,CblasTrans,l,n,1.0,backupAL,lda,0.0,Stl,l);

    for(j = 0; j < y_iters; j++)
    {
        AIOfile.load_Yblock(&Y,y_block_size);
        //!Ay_top=AL'*Y
        cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,
                                    l,y_block_size,n,1.0,AL,n,Y,n,0.0,&Ay_top[j*l*y_block_size],l);
    }


    for(i = 0; i < a_iters; i++)
    {

        if(a_iters > 10 && (i%(a_iters/10))==0)
        {
            cout << "%" << flush;
        }


        AIOfile.load_ARblock(&backupAR, a_block_size);
        AIOfile.reset_Y();
        copy_vec(backupAR,AR,n*r*a_block_size);

        //matlab_print_matrix("A",n,p,A);
        //!Generate Str
        cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,
                                l,r*a_block_size,n,1.0,AL,n,AR,n,0.0,Str,l);

        type_precision Sbr[r*r*a_block_size];

        #pragma omp parallel for default(shared) schedule(static)
        for(int ii= 0; ii < a_block_size ; ii++)
        {
            //!Generate Sbr
            cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,
                                    r,r,n,1.0,&AR[ii*r*n],n,&AR[ii*r*n],n,0.0,&Sbr[ii*r*r],r);
        }

        type_precision Ay[y_block_size*p];
        type_precision S[p*p];


        for(j = 0; j < y_iters; j++)
        {
            if(a_iters < 10 && (y_iters < 10 || (j%(y_iters/10))==0))
            {
                cout << "%" << flush;
            }


            AIOfile.load_Yblock(&Y,y_block_size);

            //!Ay_bot=AR'*Y
            cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,
                                    r*a_block_size,y_block_size,n,1.0,AR,n,Y,n,0.0,Ay_bot,r*a_block_size);

            #pragma omp parallel for private(S,Ay) default(shared) schedule(static)
            for(int ii= 0; ii < a_block_size ; ii++)
            {


                //matlab_print_matrix("Y",n,y_block_size,Y);

                //!Rebuild AY

                for(int jj = 0; jj < y_block_size; jj++)
                {
                    copy_vec(&Ay_top[j*l*y_block_size+jj*l],&Ay[jj*p],l);
                    copy_vec(&Ay_bot[ii*r + jj*r*a_block_size],&Ay[l+jj*p],r);
                }

                type_precision* B = Ay;
                //matlab_print_matrix("Ay_top",l,y_block_size,&Ay_top[j*l*y_block_size]);
                //matlab_print_matrix("Ay_bot",r,y_block_size,Ay_bot);
                //matlab_print_matrix("Ay",p,y_block_size,Ay);

                //!Rebuild S
                build_S(S,Stl,&Str[ii*r*l],&Sbr[ii*r*r],l,r);
                //matlab_print_matrix("S",p,p,S);

                //!b=S\Ay
                info = LAPACKE_sposv(STORAGE_TYPE,'U',p,y_block_size,S,p,Ay,p);
                //assert(info == 0,"POSV");


                if( ForceCheck)
                {
                    check_result(backupAL,&AR[ii*r*n],n,p,y_block_size,r,Y,B);
                }
            }

        }



    }

    get_ticks(end_tick);
    out.duration = ticks2sec(end_tick-start_tick,cpu_freq);
    //out.gflops = a_amount/1000.0*(n*l*r+n*r*r+y_amount*(n*r+p*p*p)))/1000.0/1000.0;
    out.gflops = gemm_flops(l,n,l,0) + gemm_flops(l,n,y_amount,0) +
                        a_amount*(gemm_flops(l,n,r,0) + gemm_flops(r,n,r,0) +
                                    gemm_flops(r,n,y_amount,0)+(y_amount/1000.0)*(p*p*p/3.0)/1000.0/1000.0);


    AIOfile.finalize();

    delete []Ay_top;
    delete []Ay_bot;
    delete []A;
    delete []AL;
    delete []AR;

}





void Algorithm::partialNEQ(struct Settings params, struct Outputs &out)
{
    int max_threads = params.threads;


    srand (time(NULL));


    blas_set_num_threads(max_threads);


    type_precision *Ytemp;
    lapack_int info,n,lda,ldy,l,r,k,p;

    int i,j,w;
    int m;

    cputime_type start_tick,start_tick2, end_tick;
    AIOfile.initialize(params);

    n = params.n; l=params.l; r=params.r; p = l+r;

    int y_amount = params.m;
    int y_block_size = params.mb;//kk

    int a_amount = params.t;
    int a_block_size = 1;

    int a_iters = a_amount/a_block_size;

    lda = n; ldy = n;
    k = p;



    type_precision* backupAL;

    AIOfile.load_AL(&backupAL);

    type_precision* S = new type_precision[p*p];
    type_precision Stl[l*l];
    type_precision Str[l*r];
    type_precision Sbr[r*r];
    type_precision* A = new type_precision[n*p];
    type_precision* AL = A;
    type_precision* AR = &A[l*n];

    copy_vec(backupAL,AL,n*l);


    int y_iters = y_amount/y_block_size;

    type_precision* Y;
    type_precision* backupAR;

    //printf("\n\n%%Computations\n%%");

    get_ticks(start_tick);

    //!Generate S
    cblas_ssyrk(CblasColMajor,CblasUpper,CblasTrans,l,n,1.0,backupAL,lda,0.0,Stl,l);

    for(i = 0; i < a_iters; i++)
    {

        if(a_iters < 10 || (i%(a_iters/10))==0)
        {
            cout << "%" << flush;
        }

        AIOfile.load_ARblock(&backupAR, a_block_size);
        AIOfile.reset_Y();
        copy_vec(backupAR,AR,n*r);

        //matlab_print_matrix("A",n,p,A);
        //!Generate Str
        cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,
                                l,r,n,1.0,AL,n,AR,n,0.0,Str,l);

        //!Generate Sbr
        cblas_ssyrk(CblasColMajor,CblasUpper,CblasTrans,r,n,1.0,backupAR,lda,0.0,Sbr,r);



        //matlab_print_matrix("S",p,p,S);

        for(j = 0; j < y_iters; j++)
        {

            build_S(S,Stl,Str,Sbr,l,r);

            AIOfile.load_Yblock(&Y,y_block_size);
            type_precision Ay[y_block_size*p];

            //matlab_print_matrix("Y",n,y_block_size,Y);

            //!Ay=A'*Y
            cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,
                                p,y_block_size,n,1.0,A,n,Y,n,0.0,Ay,p);

            type_precision* B = Ay;
            //matlab_print_matrix("Ay",p,y_block_size,Ay);

            //!b=S\qy
            info = LAPACKE_sposv(STORAGE_TYPE,'U',p,y_block_size,S,p,Ay,p);
            assert(info == 0,"POSV");

            //matlab_print_matrix("B",p,y_block_size,B);
            //exit(0);
            if( ForceCheck)
            {
                check_result(backupAL,backupAR,n,p,y_block_size,r,Y,B);
            }
        }


    }

        get_ticks(end_tick);
    out.duration = ticks2sec(end_tick-start_tick,cpu_freq);
    out.gflops = 2.0*(n*l*l/1000.0+a_amount/1000.0*(n*l*r+n*r*r+y_amount*(p*n+p*p*p)))/1000.0/1000.0;



    AIOfile.finalize();

    delete []A;
    delete []S;

}


void Algorithm::fullNEQ(struct Settings params, struct Outputs &out)
{
    int max_threads = params.threads;


    srand (time(NULL));


    blas_set_num_threads(max_threads);


    type_precision *Ytemp;
    lapack_int info,n,lda,ldy,l,r,k,p;

    int i,j,w;
    int m;

    cputime_type start_tick,start_tick2, end_tick;

    AIOfile.initialize(params);

    n = params.n; l=params.l; r=params.r; p = l+r;

    int y_amount = params.m;
    int y_block_size = params.mb;//kk

    int a_amount = params.t;
    int a_block_size = 1;

    int a_iters = a_amount/a_block_size;

    lda = n; ldy = n;
    k = p;



    type_precision* backupAL;

    AIOfile.load_AL(&backupAL);

    type_precision* S = new type_precision[p*p];
    type_precision* Stemp = new type_precision[p*p];
    type_precision* A = new type_precision[n*p];
    type_precision* AL = A;
    type_precision* AR = &A[l*n];

    int y_iters = y_amount/y_block_size;

    type_precision* Y;
    type_precision* backupAR;

    //printf("\n\n%%Computations\n%%");

    get_ticks(start_tick);

    for(i = 0; i < a_iters; i++)
    {

        if(a_iters < 10 || (i%(a_iters/10))==0)
        {
            cout << "%" << flush;
        }

        copy_vec(backupAL,AL,n*l);

        AIOfile.load_ARblock(&backupAR, a_block_size);
        AIOfile.reset_Y();
        copy_vec(backupAR,AR,n*r);

        //matlab_print_matrix("A",n,p,A);
        //!Generate S
        cblas_ssyrk(CblasColMajor,CblasUpper,CblasTrans,p,n,1.0,A,lda,0.0,S,p);
        //matlab_print_matrix("S",p,p,S);

        for(j = 0; j < y_iters; j++)
        {

            AIOfile.load_Yblock(&Y,y_block_size);
            type_precision Ay[y_block_size*p];
            copy_vec(S,Stemp,p*p);

            //matlab_print_matrix("Y",n,y_block_size,Y);

            //!Ay=A'*Y
            cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,
                                p,y_block_size,n,1.0,A,n,Y,n,0.0,Ay,p);

            type_precision* B = Ay;
            //matlab_print_matrix("Ay",p,y_block_size,Ay);

            //!b=S\qy
            info = LAPACKE_sposv(STORAGE_TYPE,'U',p,y_block_size,Stemp,p,Ay,p);
            assert(info == 0,"POSV");

            //matlab_print_matrix("B",p,y_block_size,B);
            //exit(0);
            if( ForceCheck)
            {
                check_result(backupAL,backupAR,n,p,y_block_size,r,Y,B);
            }
        }


    }

    get_ticks(end_tick);
    out.duration = ticks2sec(end_tick-start_tick,cpu_freq);
    out.gflops = 2.0*a_amount*(n*p*p/1000.0+y_amount/1000.0*(p*n+p*p*p))/1000.0/1000.0;


    AIOfile.finalize();


    free(AL);
    delete []A;
    delete []S;
    delete []Stemp;
}



//!*************************************************!//

void Algorithm::fullQR(struct Settings params, struct Outputs &out)
{
    int max_threads = params.threads;


    srand (time(NULL));

    blas_set_num_threads(max_threads);



    cputime_type start_tick,start_tick2, end_tick;
    double acc_gemm=0;
    double acc_trsm=0;
    double acc_other = 0;
    double acc_atemp=0;
    double acc_rtr=0;
    double acc_rqr=0;


    type_precision *Ytemp;
    lapack_int info,n,lda,ldy,l,r,k,p;

    int i,j,w;
    int m;

    AIOfile.initialize(params);

    n = params.n; l=params.l; r=params.r; p = l+r;

    int y_amount = params.m;
    int y_block_size = 1;//kk

    int a_amount = params.t;
    int a_block_size = 1;

    int a_iters = a_amount/a_block_size;

    lda = n; ldy = n;
    k = p;



    type_precision* backupAL;

    AIOfile.load_AL(&backupAL);

    type_precision* Q = new type_precision[n*p];
    type_precision* A = Q;
    type_precision* AL = A;
    type_precision* AR = &A[l*n];

    type_precision tau[k];


    int y_iters = y_amount/y_block_size;


    type_precision* Y;
    type_precision* backupAR;



    //printf("\n\n%%Computations\n%%");

    get_ticks(start_tick);
    for(i = 0; i < a_iters; i++)
    {


        if(a_iters < 10 || (i%(a_iters/10))==0)
        {
            cout << "%" << flush;
        }

        copy_vec(backupAL,AL,n*l);

        AIOfile.load_ARblock(&backupAR, a_block_size);
        AIOfile.reset_Y();
        copy_vec(backupAR,AR,n*r);

        //!Generate R
        info = LAPACKE_sgeqrf(STORAGE_TYPE,n,p,A,lda,tau);
        assert(info == 0,"QR decomp");
        type_precision* R = extract_R(A,n,p);

        //!generate Q
        info = LAPACKE_sorgqr(STORAGE_TYPE,n,p,k,Q,lda,tau);
        assert(info == 0,"Q form");



//        matlab_print_matrix("ALL",n,l,backupAL);
//        matlab_print_matrix("ARR",n,r,backupAR);

        for(j = 0; j < y_iters; j++)
        {

            //matlab_print_matrix("RR",p,p,R);
            //matlab_print_matrix("QQ",n,p,Q);

            AIOfile.load_Yblock(&Y,y_block_size);
            type_precision Qy[y_block_size*p];

            //!qy=Q'*Y
            cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,
                                p,y_block_size,n,1.0,Q,n,Y,n,0.0,Qy,p);

            type_precision* B = Qy;

            //!b=R\qy
            cblas_strsm(CblasColMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,
                                                p,y_block_size,1.0,R,p,Qy,p);

            if( ForceCheck)
            {
                check_result(backupAL,backupAR,n,p,y_block_size,r,Y,B);
            }
        }

        free(R);


    }

    get_ticks(end_tick);
    out.duration = ticks2sec(end_tick-start_tick,cpu_freq);
    out.gflops = 2.0*a_amount*(2.0*n*p*p/1000.0+y_amount/1000.0*(p*n+p*p*p))/1000.0/1000.0;


    AIOfile.finalize();


    free(AL);
    delete []Q;
}

void Algorithm::partialQR(struct Settings params, struct Outputs &out)
{
    int max_threads = params.threads;


    srand (time(NULL));

    blas_set_num_threads(max_threads);



    cputime_type start_tick,start_tick2, end_tick;
    double acc_gemm=0;
    double acc_trsm=0;
    double acc_other = 0;
    double acc_atemp=0;
    double acc_rtr=0;
    double acc_rqr=0;


    type_precision *RL, *Ytemp, *RtlRtr;
    lapack_int info,n,lda,ldy,l,r,k,p;

    int i,j,w;
    int m;

    AIOfile.initialize(params);

    n = params.n; l=params.l; r=params.r; p = l+r;

    int y_amount = params.m;
    int y_block_size = params.mb;//kk

    int a_amount = params.t;
    int a_block_size = 1;

    int a_iters = a_amount/a_block_size;

    lda = n; ldy = n;
    k = l;



    type_precision* backupAL;

    AIOfile.load_AL(&backupAL);

    type_precision* Q = new type_precision[n*p];
    type_precision* AL = Q;
    copy_vec(backupAL,AL,n*l);


    type_precision tau[k];

    type_precision* QL;
    type_precision* QR = &Q[n*l];

    int y_iters = y_amount/y_block_size;


    type_precision* Y;

    type_precision Rtr[l*r];
    type_precision* Ar;
    type_precision Rbr[r*r];

    get_ticks(start_tick);
    //!Generate RTL
    info = LAPACKE_sgeqrf(STORAGE_TYPE,n,l,AL,lda,tau);
    assert(info == 0,"QR decomp");

    type_precision* R =prepare_R(AL,n,l,r);
    type_precision* Rtl = R;

    QL=AL;//same as Q
    AL=backupAL;
    //!generate QL
    info = LAPACKE_sorgqr(STORAGE_TYPE,n,l,k,QL,lda,tau);

    assert(info == 0,"Q form");

    //printf("\n\n%%Computations\n%%");


    for(i = 0; i < a_iters; i++)
    {

        if(a_iters < 10 || (i%(a_iters/10))==0)
        {
            cout << "%" << flush;
        }

        AIOfile.load_ARblock(&Ar, a_block_size);
        AIOfile.reset_Y();

        //!Rtr=Q1'*AR
        cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,
                        l,r*a_block_size,n,1.0,QL,n,Ar,n,0.0,Rtr,l);

        type_precision* Atemp = replicate_vec(Ar,n*r*a_block_size);

        //!Atemp=AR-Q1*Rtr
        cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                        n,r*a_block_size,l,-1.0,QL,n,Rtr,l,1.0,Atemp,n);

        type_precision* QrRbr = QR;
        copy_vec(Atemp,QrRbr,n*r);
        info = LAPACKE_sgeqrf(STORAGE_TYPE,n,r,QrRbr,lda,tau);
        assert(info == 0,"QR decomp of QrRbr");
        type_precision* Rbrtemp = prepare_R(QrRbr,n,r,0);//reuse of old function
        copy_vec(Rbrtemp,Rbr,r*r);
        info = LAPACKE_sorgqr(STORAGE_TYPE,n,r,r,QrRbr,lda,tau);
        assert(info == 0,"QR form");
        free(Rbrtemp);

        update_R(R, Rtr, Rbr,p, p,r);

//        matlab_print_matrix("RL",p,l,Rtl);
//        matlab_print_matrix("Rtr",l,r,Rtr);
//        matlab_print_matrix("Rbr",r,r,Rbr);
//        matlab_print_matrix("RR",p,p,R);
//        matlab_print_matrix("QQ",n,p,Q);
        //matlab_print_matrix("ALL",n,l,AL);
        //matlab_print_matrix("ARR",n,r,Ar);
        //matlab_print_matrix("R",p,p,R);
        int top_idx=0;

        for(j = 0; j < y_iters; j++)
        {
            AIOfile.load_Yblock(&Y,y_block_size);
            type_precision Qy[y_block_size*p];


            //!qy=Q'*Y
            cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,
                                p,y_block_size,n,1.0,Q,n,Y,n,0.0,Qy,p);

            type_precision* B = Qy;

            //!b=R\qy
            cblas_strsm(CblasColMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,
                                                p,y_block_size,1.0,R,p,Qy,p);

            //matlab_print_matrix("bcomp",p,y_block_size,B);


            if( ForceCheck)
            {
                check_result(AL,Ar,n,p,y_block_size,r,Y,B);
            }
        }

        free(Atemp);


    }


    get_ticks(end_tick);
    out.duration = ticks2sec(end_tick-start_tick,cpu_freq);
    out.gflops = 2.0*(n*l*l/1000.0 + a_amount*(2.0*n*l*r/1000.0+n*r*r/1000.0+y_amount/1000.0*(p*n+p*p*p)))/1000.0/1000.0;


    AIOfile.finalize();


    free(Rtl);
    free(AL);
    delete []Q;
}


void  Algorithm::partialQR_Blocked_Rtl(struct Settings params, struct Outputs &out)
{
    int max_threads = params.threads;

    //srand (time(NULL));

    int cpu_cores = max_threads;


    blas_set_num_threads(max_threads);
    omp_set_num_threads(max_threads);




    cputime_type start_tick,start_tick2, end_tick;
    double acc_pre=0;
    double acc_trsm=0;
    double acc_other = 0;
    double acc_atemp=0;
    double acc_rtr=0;
    double acc_rqr=0;
    double acc_RTL_QLY=0;


    type_precision *AL,*RL, *RQy_top, *Ytemp, *RtlRtr;
    lapack_int info,n,lda,ldy,l,r,k,p;

    int i,j,w;
    int m;


    AIOfile.initialize(params);

    n = params.n; l=params.l; r=params.r; p = l+r;



    double lvl_3_cache_size = 6;


    int y_amount = params.t;
    int y_block_size = params.tb;//kk

    int a_amount = params.m;
    int a_block_size = params.mb;

    int a_iters = (a_amount + a_block_size - 1) / a_block_size;

    int y_iters = (y_amount+ y_block_size - 1)/y_block_size;

    int qy_idx=0;

    out.gflops = n/1000.0*l/1000.0*l/1000.0 + gemm_flops(l,n,y_amount,0) + y_amount*l/1000.0*l/1000.0/1000.0 +
                    a_amount * (gemm_flops(l,n,r,0) + l/1000.0*l/1000.0/1000.0 + gemm_flops(n,l,r,1) +
                                    n/1000.0*r/1000.0*r/1000.0 + gemm_flops(r,n,y_amount,0) + y_amount/1000.0*r/1000.0*r/1000.0 +
                                            y_amount/1000.0*l/1000.0/1000.0);

    int sch_block_size = max(1.0,min((double)a_block_size/(double)max_threads, (1.0*1000*1000)/(double)((sizeof(type_precision)*n*r))));

    //cout << endl<< "taskChunk "<< sch_block_size <<endl;

    lda = n; ldy = n;
    k = l;

    AIOfile.load_AL(&AL);

    type_precision* backupAL = replicate_vec(AL,n*l);



    RQy_top = new type_precision[l*y_amount];

    type_precision* Bfinal = (type_precision*)malloc(max_threads*p*y_block_size*sizeof(type_precision));

    type_precision tau[k];

    type_precision* QL;

    type_precision* Rtr = (type_precision*)malloc(l*r*a_block_size*sizeof(type_precision));
    type_precision* Qy_bot = (type_precision*)malloc(a_block_size*y_block_size*r*sizeof(type_precision));

    type_precision* B_top = new type_precision[max_threads*y_block_size*l];
    type_precision* B_bot = new type_precision[max_threads*y_block_size*r];

    type_precision* Ar;
    type_precision Rbr[r*r*a_block_size];



    //printf("\n\n%%Computations\n%%");
    for(i = 0; i < a_iters; i++)
    {

        if(a_iters >= 10 && (i%(a_iters/10))==0)
        {
            //cout << "*" << flush;
        }
        for(j = 0; j < y_iters; j++)
        {

            if(a_iters < 10 && (y_iters < 10 || (j%(y_iters/10))==0))
            {
                //cout << "*" << flush;
            }
        }
    }
    cout << endl;


    //!Start of Computations

    get_ticks(start_tick);


    info = LAPACKE_sgeqrf(STORAGE_TYPE,n,l,AL,lda,tau);
    assert(info == 0,"QR decomp");

    type_precision* Rtl=prepare_R(AL,n,l,0);

    QL=AL;
    AL = backupAL;

    info = LAPACKE_sorgqr(STORAGE_TYPE,n,l,k,QL,lda,tau);

    assert(info == 0,"Q form");

    matlab_print_matrix("QL",n,l,QL);

    qy_idx=0;

    //printf("\n%%Preparing IO\n");
    //!Prepare File for the Y


    type_precision* Y;

    //#pragma omp parallel default(shared)
//    {
//        //#pragma omp for private(Y,y_block_size,qy_idx) nowait schedule(static,1)
//        for(int j = 0; j < y_iters; j++)
//        {
//            //qy_idx = j*y_block_size*l;
//            //#pragma omp critical
//            {AIOfile.load_Yblock(&Y,y_block_size);}
//
//            //!qy_top=QL'*Y
//            cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,
//                            l,y_block_size,n,1.0,QL,n,Y,n,0.0,&RQy_top[qy_idx],l);
//
//            //!K | RtlQlY =RTL\qy_top
//            cblas_strsm(CblasColMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,
//                                            l,y_block_size,1.0,Rtl,l,&RQy_top[qy_idx],l);
//            qy_idx += y_block_size*l;
//        }
//    }

    get_ticks(end_tick);

    out.acc_pre = ticks2sec(end_tick-start_tick,cpu_freq);

    //cout << "\npre " << out.acc_pre << endl;


    y_block_size = params.tb;


    for(int i = 0; i < a_iters; i++)
    {
        cout << i<<":i\n" << flush;
        if(a_iters >= 10 && (i%(a_iters/10))==0)
        {
            cout << AIOfile.io_overhead << flush;
            AIOfile.io_overhead = "*";
        }

        get_ticks(start_tick2);

        AIOfile.load_ARblock(&Ar, a_block_size);

        get_ticks(end_tick);
        out.acc_loadxr += ticks2sec(end_tick-start_tick2,cpu_freq);






        get_ticks(start_tick2);

        type_precision* Qr;
        type_precision* Atemp;

        #pragma omp parallel default(shared)
        {

             #pragma omp sections
             {
                   {
                       //!Rtr=Q1'*AR
                        cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,
                                            l,r*a_block_size,n,1.0,QL,n,Ar,n,0.0,Rtr,l);
                        //matlab_print_matrix("QL'*Ar",l,r*a_block_size,Rtr);

                   }
                   #pragma omp section
                   {
                        Qr = (type_precision*)malloc(n*r*a_block_size*sizeof(type_precision));
                   }
                   #pragma omp section
                   {
                        Atemp = replicate_vec(Ar,n*r*a_block_size);
                   }


             }
             #pragma omp barrier
     }



        //!Atemp=AR-Q1*Rtr
        cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                        n,r*a_block_size,l,-1.0,QL,n,Rtr,l,1.0,Atemp,n);

        //!H=RTL\Rtr
        cblas_strsm(CblasColMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,
                                    l,r*a_block_size,1.0,Rtl,l,Rtr,l);

        #pragma omp parallel default(shared)
        {
            #pragma omp for private(tau) nowait  schedule(static,sch_block_size)
            for(int ii = 0; ii<a_block_size; ii++)
            {
                type_precision* QrRbr = &Qr[n*r*ii];
                copy_vec(&Atemp[n*r*ii],QrRbr,n*r);


                info = LAPACKE_sgeqrf(STORAGE_TYPE,n,r,QrRbr,lda,tau);
                //assert(info == 0,"QR decomp of QrRbr");
                type_precision* Rbrtemp = prepare_R(QrRbr,n,r,0);//reuse of old function
                copy_vec(Rbrtemp,&Rbr[r*r*ii],r*r);


                info = LAPACKE_sorgqr(STORAGE_TYPE,n,r,r,QrRbr,lda,tau);
                //assert(info == 0,"QR form");
                free(Rbrtemp);
            }
        }


        //matlab_print_matrix("RTR",l,r,Rtr);
        //matlab_print_matrix("RTL",l,l,Rtl);




        get_ticks(end_tick);
        out.firstloop += ticks2sec(end_tick-start_tick2,cpu_freq);

        RtlRtr=Rtr;

        int top_idx=0;
        int qy_bot_idx=0;
        for(int j = 0; j < y_iters; j++)
        {

            if(a_iters < 10 && (y_iters < 10 || (j%(y_iters/10))==0))
            {
                cout << AIOfile.io_overhead << flush;
                AIOfile.io_overhead = "*";
            }

            get_ticks(start_tick2);

            AIOfile.load_Yblock(&Y,y_block_size);

            get_ticks(end_tick);
            out.acc_loady += ticks2sec(end_tick-start_tick2,cpu_freq);

            cout << j<<":j\n" << flush;
            get_ticks(start_tick2);
            if(i == 0)
            {
                matlab_print_matrix("QLpre",n,l,QL);

                matlab_print_matrix("Ypre",n,y_block_size,Y);

                 //!qy_top=QL'*Y
                cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,
                                l,y_block_size,n,1.0,QL,n,Y,n,0.0,&RQy_top[qy_idx],l);

               matlab_print_matrix("Ypost",n,y_block_size,Y);

                matlab_print_matrix("QL'*Y",l,y_block_size,&RQy_top[qy_idx]);

                //!K | RtlQlY =RTL\qy_top
                cblas_strsm(CblasColMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,
                                                l,y_block_size,1.0,Rtl,l,&RQy_top[qy_idx],l);
                qy_idx += y_block_size*l;


            }
            get_ticks(end_tick);
            out.acc_RTL_QLY += ticks2sec(end_tick-start_tick2,cpu_freq);


            get_ticks(start_tick2);

            matlab_print_matrix("Yprepre",n,y_block_size,Y);

            //!qy_bot=Qr'*Y
            cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,
                            r*a_block_size,y_block_size,n,1.0,Qr,n,Y,n,0.0,Qy_bot,r*a_block_size);
            //matlab_print_matrix("Qr'*Y",r*a_block_size,y_block_size,Qy_bot);


            matlab_print_matrix("Ypostpost",n,y_block_size,Y);


            get_ticks(end_tick);
            out.acc_gemm += ticks2sec(end_tick-start_tick2,cpu_freq);


            get_ticks(start_tick2);


            #pragma omp parallel default(shared)
            {
                #pragma omp for nowait  schedule(static,sch_block_size)
                for(int ii = 0; ii<a_block_size; ii++)
                {

                    int thread_id = omp_get_thread_num();
                    int b_bot_idx = thread_id*y_block_size*r;
                    int b_top_idx = thread_id*y_block_size*l;
                    int b_final_idx = thread_id*y_block_size*p;

                    extract_subMatrix(Qy_bot,&B_bot[b_bot_idx],r*a_block_size,y_block_size,ii*r,(ii+1)*r,0, y_block_size);
                    copy_vec(&RQy_top[top_idx],&B_top[b_top_idx],l*y_block_size);

                    //!bBot=Rbr\qy_bot
                    cblas_strsm(CblasColMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,
                                                        r,y_block_size,1.0,&Rbr[ii*r*r],r,&B_bot[b_bot_idx],r);

                    //!Btop=K-H*bBot     Btop=(Rtl\(Ql'*Y))-(Rtl\Rtr)*(Rbr\(Qr'*Y));
                    cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,
                                l,y_block_size,r,-1.0,&RtlRtr[ii*l*r],l,&B_bot[b_bot_idx],r,1.0,&B_top[b_top_idx],l);

//                    matlab_print_matrix("top",l,y_block_size,&B_top[b_top_idx]);
//                    matlab_print_matrix("bot",r,y_block_size,&B_bot[b_bot_idx]);


                    if(OUTPUT)
                    {
                        prepare_Bfinal(&Bfinal[b_final_idx], &B_top[b_top_idx],&B_bot[b_bot_idx],p,y_block_size,r);
                        AIOfile.write_B(&Bfinal[b_final_idx],p,y_block_size);
                    }

                    if( ForceCheck)
                    {
                        #pragma omp critical
                        {
                            if(!OUTPUT)
                                prepare_Bfinal(&Bfinal[b_final_idx], &B_top[b_top_idx],&B_bot[b_bot_idx],p,y_block_size,r);

                            prepare_Bfinal(&Bfinal[b_final_idx], &B_top[b_top_idx],&B_bot[b_bot_idx],p,y_block_size,r);
                            //matlab_print_matrix("outsol",p,y_block_size,&Bfinal[b_final_idx]);
                            check_result(AL,&Ar[ii*n*r],n,p,y_block_size,r,Y,&Bfinal[b_final_idx]);
                        }
                    }

                }
            }

            get_ticks(end_tick);
            out.acc_b += ticks2sec(end_tick-start_tick2,cpu_freq);



            top_idx += y_block_size*l;

        }

        if(i < (a_iters-1))
            AIOfile.reset_Y();

        free(Atemp);
        free(Qr);

    }
    get_ticks(end_tick);

    out.duration = ticks2sec(end_tick-start_tick,cpu_freq);

    AIOfile.finalize();

    free(Rtr);

    delete []B_top;
    delete []B_bot;

    free(Rtl);
    delete []RQy_top;
    free(Qy_bot);
    free(AL);
    free(Bfinal);
}

