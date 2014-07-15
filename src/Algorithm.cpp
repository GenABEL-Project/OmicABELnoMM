#include "Algorithm.h"



Algorithm::Algorithm()
{
    // ctor
}


Algorithm::~Algorithm()
{
    // dtor
}


void Algorithm::solve(struct Settings params, struct Outputs &out, int type)
{
    switch (type)
    {
        case P_NEQ_B_OPT_MD:
            partialNEQ_Blocked_STL_MD(params, out);
            break;

        default:
            break;
    }
}


void Algorithm::extract_subMatrix(type_precision* source, type_precision* dest,
                                  int dim1_source, int dim2_source,
                                  int dim1_ini, int dim1_end, int dim2_ini,
                                  int dim2_end)
{
    int i, j, idx = 0;
    int size, source_ini;
    for (i = dim2_ini; i < dim2_end; i++)
    {
        j = dim1_ini;
        source_ini = i * dim1_source+j;
        size = dim1_end - dim1_ini;
        memcpy( (type_precision*)&dest[idx],
                (type_precision*)&source[source_ini],
                size * sizeof(type_precision) );
//        for (j = dim1_ini; j<dim1_end; j++)
//        {
//            dest[idx] = source[i*dim1_source+j];
//            idx++;
//        }
        idx += size;
    }
}


void Algorithm::prepare_Bfinal(type_precision* bfinal, type_precision* bsource, int a_amount, int y_amount, int p)
{
//    // memcpy are faster version of the fors
//    int i, k, w, top_idx, bot_idx;
//    int size;
//    top_idx = 0;
//    bot_idx = 0;
//    for (k = 0; k < dim2_b; k++)
//    {
//        size = k * dim1_b + (dim1_b - dim1_b_bot) - (k * dim1_b);
//        memcpy( (type_precision*)&bfinal[k * dim1_b],
//                (type_precision*)&top[top_idx],
//                size * sizeof(type_precision) );
////        for (i = k*dim1_b; i < k*dim1_b+(dim1_b-dim1_b_bot); i++)
////        {
////            bfinal[i] = top[top_idx];
////            top_idx++;
////        }
//        top_idx += size;
//        i = k * dim1_b + size;
//        w = i;
//
//        size = w + dim1_b_bot - w;
//        memcpy( (type_precision*)&bfinal[w],
//                (type_precision*)&bot[bot_idx],
//                size * sizeof(type_precision) );
////        for (i = w; i < w+dim1_b_bot; i++)
////        {
////            bfinal[i] = bot[bot_idx];
////            bot_idx++;
////        }
//        bot_idx += size;
//    }
}


void Algorithm::prepare_QY(type_precision* qy, type_precision* top,
                           type_precision* bot, int dim1_QY,
                           int dim2_QY, int dim1_qy_bot, int bot_blocks )
{
    int i, k, w, top_idx, bot_idx;
    top_idx = 0;
    bot_idx = 0;
    for (k = 0; k < dim2_QY; k++)
    {
        for (i = k * dim1_QY; i < (k+1) * dim1_QY - dim1_qy_bot; i++)
        {
            qy[i] = top[top_idx];
            top_idx++;
        }
        w = i;

        for (i = w; i < w + dim1_qy_bot; i++)
        {
            qy[i] = bot[bot_idx];
            bot_idx++;
        }
        bot_idx += (bot_blocks-1) * dim1_qy_bot;
    }
}


type_precision* Algorithm::extract_R(type_precision* A, int dim1_A, int dim2_A)
{
    type_precision* R = (type_precision*)calloc(dim2_A * dim2_A,
                                                sizeof(type_precision));
    int i, j;

    int R_idx = 0;

    for (i = 0; i < dim2_A; i++)
    {
        for (j = 0; j <= i; j++)
        {
            R[R_idx] = A[j + i * dim1_A];
            R_idx++;
        }
        R_idx = dim2_A * (i+1);
    }
    return R;
}


type_precision* Algorithm::prepare_R(type_precision* RL, int dim1_A,
                                     int dim2_AL, int dim2_AR)
{
    int R_dim = (dim2_AR + dim2_AL);
    type_precision* R = new type_precision[R_dim * R_dim];

    int i, j;

    int RL_idx = 0;
    int R_idx = 0;

    for (i = 0; i < dim2_AL; i++)
    {
        for (j = 0; j <= i; j++)
        {
            RL_idx   = i * dim1_A + j;
            R_idx    = i * R_dim + j;
            R[R_idx] = RL[RL_idx];
        }
    }
    return R;
}


void Algorithm::update_R(type_precision* R, type_precision* topRr,
                         type_precision* botRr, int dim1, int dim2, int r)
{
    int i, j, w;
    int max = dim1 * dim2;
    int rtr_idx = 0;
    int rbr_idx = 0;
    for (j = r; j > 0; j--)
    {
        for (i = max - dim1 * j; i < max - dim1 * j + dim2 - r; i++)
        {
            R[i] = topRr[rtr_idx];
            rtr_idx++;
        }
        w = i;

        for (i = w; i < w + r; i++)
        {
            R[i] = botRr[rbr_idx];
            rbr_idx++;
        }
    }
}


void Algorithm::build_S(type_precision* S, type_precision* Stl,
                        type_precision* Str, type_precision* Sbr, int l, int r)
{
    int Sidx;
    int p = l + r;

    for (int i = 0; i < l; i++)
    {
        Sidx = i * p;
        for (int j = 0; j <= i; j++)
        {
            S[Sidx] = Stl[j+i*l];
            Sidx++;
        }
    }

    for (int i = 0; i < r; i++)
    {
        Sidx = l*p + p*i;
        for (int j = 0; j < l; j++)
        {
            S[Sidx] = Str[j + i * l];
            Sidx++;
        }
    }

    for (int i = 0; i < r; i++)
    {
        Sidx = l*p + l + p*i;
        for (int j= 0; j <= i; j++)
        {
            S[Sidx] = Sbr[j + i * r];
            Sidx++;
        }
    }
}


void Algorithm::check_result(type_precision* AL, type_precision* AR,
                             int rowsA, int colsA, int rhs, int colsAR,
                             type_precision* y, type_precision* res)
{
    type_precision* A = (type_precision*)malloc(rowsA * colsA *
                                                sizeof(type_precision));

    int i, ar_idx = 0;
    for (i = 0; i < rowsA * (colsA - colsAR); i++)
    {
        A[i] =  AL[i];
    }

    for (i = rowsA * (colsA - colsAR); i < rowsA * colsA; i++)
    {
        A[i] = AR[ar_idx];
        ar_idx++;
    }


    type_precision* ynew = replicate_vec(y, rowsA*rhs);
    type_precision* new_sol = (type_precision*)malloc(colsA * rhs *
                                                      sizeof(type_precision));

    lapack_int info = LAPACKE_sgels(STORAGE_TYPE, 'N', rowsA, colsA, rhs, A,
                                    rowsA, ynew, rowsA);
    myassert(info == 0, "Error Check");


    int index = 0;
    int index_new = 0;
    for (i = 0; i < rhs; i++)

    {
         copy_vec(&ynew[index], &new_sol[index_new], colsA);
         index += rowsA;
         index_new += colsA;
    }


//    if (PRINT)
//        printf("\n Btop=(Rtl\\(Ql'*Y))-(Rtl\\Rtr)*(Rbr\\(Qr'*Y)); \n [Btop ; Rbr\\Qr'*Y] - bcomputed \n");

    cblas_saxpy(rhs * colsA, -1.0, res, 1, new_sol, 1);
    type_precision u_norm = cblas_snrm2(rhs * colsA, new_sol, 1);
    //
    if (abs(u_norm) >= 0.0001 || isnan(u_norm))
    {
        fflush(stdout);
        matlab_print_matrix("AL", rowsA, colsA-colsAR, AL);
        matlab_print_matrix("AR", rowsA, colsAR, AR);
        matlab_print_matrix("Y", rowsA, rhs, y);
        printf("\nA = [AL AR]; [Q, R] = qr(A, 0); rr = R\\(Q'*Y)\n");
        matlab_print_matrix("bcomputed", colsA, rhs, res);
        matlab_print_matrix("newsol", colsA, rhs, ynew);
        printf("\n%%\tnrom: %0.2g", u_norm);
         exit(1);
    }
    else
    {
        matlab_print_matrix("bcomputed", colsA, rhs, res);
        matlab_print_matrix("newsol", colsA, rhs, ynew);
        //printf("\n%%\tnrom: %0.2g", u_norm);
    }


    // cout << "\t**************";
    free(ynew);
    free(new_sol);
    free(A);
}

///////////////////////////////


void Algorithm::partialNEQ_Blocked_STL_MD(struct Settings params,
                                          struct Outputs &out)
{
    int max_threads = params.threads;


    srand(time(NULL));

    blas_set_num_threads(max_threads);
    omp_set_num_threads(max_threads);


    //type_precision *Ytemp;
    lapack_int info, n, lda, l, r, p;

    cputime_type start_tick, start_tick2, start_tick3, end_tick;

    AIOfile.initialize(params);

    n = params.n; l = params.l; r = params.r; p = l+r;

    int y_amount = params.t;
    int y_block_size = params.tb;  // kk

    int a_amount = params.m;
    int a_block_size = params.mb;

    int a_iters = (a_amount + a_block_size - 1) / a_block_size;

    int y_iters = (y_amount + y_block_size - 1) / y_block_size;


    lda = n;


    for (int j = 0; j < y_iters && !params.ForceCheck; j++)
    {
        if (!params.ForceCheck &&  y_iters >= 10 && (j%(y_iters/10)) == 0 )
        {
            cout << "*" << flush;
        }

        for (int i = 0; i < a_iters; i++)
        {

                if ( !params.ForceCheck &&  y_iters < 10 &&
                (  (a_iters >= 10 && (i%(a_iters/(10/y_iters))) == 0) || (a_iters < (10/y_iters)) ))
                {
                    cout << "*" << flush;
                }

        }
    }

    if(!params.ForceCheck)
        cout << endl;

    //add memalign

    //type_precision Stl[l*l];
    //type_precision Str[l*r*a_block_size];
    type_precision* Stl = new type_precision[l*l*1];
    type_precision* Str = new type_precision[l*r*a_block_size*1];

    type_precision* Sbr = new type_precision[r *  r * a_block_size];
    type_precision* Ay = new type_precision[p * a_block_size];
    //type_precision* B_resorted = new type_precision[p * a_block_size*y_block_size];

    type_precision* S = new type_precision[p * p];


    type_precision* Ay_top = new type_precision[l * y_amount];
    type_precision* Ay_bot = new type_precision[y_block_size * a_block_size * r];

    type_precision* y_residual = new type_precision[n * y_block_size ];
    type_precision* y_res_norms = new type_precision[a_block_size];

    list<long int>* al_nan_idxs = new list<long int>[1];
    list<long int>* y_nan_idxs = new list<long int>[y_block_size];
    list<long int>* ar_nan_idxs = new list<long int>[a_block_size];


    type_precision* A = new type_precision[n * p * 1];
    type_precision* AR = new type_precision[n * r * a_block_size * 1];
//  type_precision* AL = new type_precision[n * l * 1];
    type_precision* AL = A;

    type_precision* B = Ay;

    type_precision* backupAR;  // = new type_precision[n*r*a_block_size];
    type_precision* backupAL;  // = new type_precision[n*l];


    AIOfile.load_AL(&backupAL);

    //pthread_barrier_wait(&(AIOfile.Fhandler->finalize_barrier));

    replace_nans(al_nan_idxs,1, backupAL, n, l);
    al_nan_idxs->push_back(1);

    //LAPACKE_dgesdd()

    copy_vec(backupAL, AL, n*l);


    type_precision* Y;

    // printf("\n\n%%Computations\n%%");


    get_ticks(start_tick);

    for (int j = 0; j < y_iters; j++)
    {
        if (!params.ForceCheck &&  y_iters >= 10 && (j%(y_iters/10)) == 0 && !params.ForceCheck)
        {
            cout << AIOfile.io_overhead << flush;
            AIOfile.io_overhead = "*";
        }

        get_ticks(start_tick2);

        AIOfile.load_Yblock(&Y, y_block_size);

        get_ticks(end_tick);
        out.acc_loady += ticks2sec(end_tick,start_tick2);

        get_ticks(start_tick2);
        replace_nans(&y_nan_idxs[0],y_block_size, Y, n,1);
        //out.acc_other += ticks2sec(end_tick,start_tick2);


        get_ticks(start_tick2);

        //! Ay_top = AL'*Y
        cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    l, y_block_size, n, 1.0, AL, n, Y, n, 0.0,
                    &Ay_top[j * l * y_block_size], l);

        get_ticks(end_tick);
        out.acc_gemm += ticks2sec(end_tick,start_tick2);


        for (int i = 0; i < a_iters; i++)
        {

            if (!params.ForceCheck && y_iters < 10 &&
                (  (a_iters >= 10 && (i%(a_iters/(10/y_iters))) == 0) || (a_iters < (10/y_iters)) ))
            {
                cout << "*" << flush;
            }

            get_ticks(start_tick2);

            AIOfile.load_ARblock(&backupAR, a_block_size);

            get_ticks(end_tick);
            out.acc_loadxr += ticks2sec(end_tick,start_tick2);

            get_ticks(start_tick2);

            replace_nans(&ar_nan_idxs[0],a_block_size, backupAR, n, r);
            replace_with_zeros(al_nan_idxs, backupAR,  n, r, a_block_size);

            copy_vec(backupAR, AR, n *  r * a_block_size);

            get_ticks(end_tick);
            //out.acc_other += ticks2sec(end_tick,start_tick2);



            get_ticks(start_tick2);

            //! Ay_bot = AR'*Y
            cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                        r * a_block_size, y_block_size, n, 1.0, AR, n, Y, n,
                        0.0, Ay_bot, r * a_block_size);

            get_ticks(end_tick);
            out.acc_gemm += ticks2sec(end_tick,start_tick2);


            get_ticks(start_tick3);
            for (int jj = 0; jj < y_block_size; jj++)
            {

                //int thread_id = 0 * omp_get_thread_num();//so far singel thread version only


                get_ticks(start_tick2);

                copy_vec(backupAL, AL, n * l);//try to remove!

                replace_with_zeros(&y_nan_idxs[jj], AL, n, l, 1);


                get_ticks(end_tick);//2%
                out.acc_other += ticks2sec(end_tick,start_tick2);

                 get_ticks(start_tick2);

                //! Generate Stl
                cblas_ssyrk(CblasColMajor, CblasUpper, CblasTrans,
                            l, n, 1.0, AL, lda, 0.0, Stl, l);

                get_ticks(end_tick);
                out.acc_stl += ticks2sec(end_tick,start_tick2);


//                get_ticks(start_tick2);

                //copy_vec(backupAR,AR, n*r*a_block_size);//!10%//try to remove!


                //replace_with_zeros(&y_nan_idxs[jj], AR,  n, r, a_block_size);



//                get_ticks(end_tick);
//                out.acc_other += ticks2sec(end_tick,start_tick2);

                get_ticks(start_tick2);

                //! Generate Str
                cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                            l, r * a_block_size, n, 1.0, AL,
                            n, AR, n, 0.0, Str, l);//!45

                get_ticks(end_tick);
                out.acc_str += ticks2sec(end_tick,start_tick2);



                //type_precision Sbr[r *  r * a_block_size*1];
                //type_precision Ay[p * a_block_size*1];

                blas_set_num_threads(1);
                omp_set_num_threads(max_threads);

                #pragma omp parallel default(shared)
                {

                #pragma omp for nowait
                for (int ii= 0; ii < a_block_size; ii++)
                {
                    //cout << omp_get_thread_num() << endl << flush;

                    get_ticks(start_tick2);

                    copy_vec(&backupAR[ii*r*n],&AR[ii*r*n], n*r);//!10%//try to remove!
                    replace_with_zeros(&y_nan_idxs[jj], &AR[ii*r*n],  n, r, 1);

                    get_ticks(end_tick);
                    out.acc_other += ticks2sec(end_tick,start_tick2);


                    get_ticks(start_tick2);

                    //! Generate Sbr
                    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                                r, r, n, 1.0, &AR[ii*r*n], n,
                                &AR[ii * r * n], n, 0.0,
                                &Sbr[ii * r * r], r);


                    get_ticks(end_tick);
                    out.acc_sbr += ticks2sec(end_tick,start_tick2 );

                    get_ticks(start_tick2);

                    copy_vec(&Ay_top[j*l*y_block_size+jj*l], &Ay[ii*p], l);
                    copy_vec(&Ay_bot[ii*r + jj*r*a_block_size],
                             &Ay[l+ii*p], r);


                    type_precision S2[p * p];

                    //! Rebuild S
                    build_S(S2, Stl, &Str[ii*r*l], &Sbr[ii*r*r], l, r);
                    // matlab_print_matrix("S", p, p, S);


                    get_ticks(end_tick);//5%
                    out.acc_other += ticks2sec(end_tick,start_tick2);

                    get_ticks(start_tick2);


                    //! b = S\Ay
                    info = LAPACKE_sposv(STORAGE_TYPE, 'U', p, 1, S2, p,
                                         &Ay[ii*p], p);

                    myassert(info == 0, "S\\Ay");


                    get_ticks(end_tick);
                    out.acc_solve += ticks2sec(end_tick,start_tick2 );


                    if (params.ForceCheck)
                    {
                        #pragma omp critical
                        {
                            check_result(AL, &AR[ii*r*n], n, p,
                                         1, r, &Y[jj*n], &Ay[ii*p]);
                        }
                    }
                }
                }

                /************statistics**************/
//                type_precision* T;
//                type_precision* R2;
//                type_precision* P;

                //hpc_statistics(n,A,a_block_size,Y,y_block_size,B,p,T,R2,P);

                /**************************/

                blas_set_num_threads(max_threads);

                get_ticks(start_tick2);

                AIOfile.write_B(B, p, a_block_size);

                get_ticks(end_tick);
                out.acc_storeb += ticks2sec(end_tick,start_tick2);





            }

            get_ticks(end_tick);
            out.acc_real_innerloops += ticks2sec(end_tick ,start_tick3);


        }
         AIOfile.reset_AR();
    }

    get_ticks(end_tick);

    out.duration = ticks2sec(end_tick ,start_tick);


    out.gflops = y_iters * (gemm_flops(l, n, params.tb*1, 0) +
                             a_iters * ( gemm_flops(params.mb*r, n, params.tb*1, 0) +
                             params.tb *(  gemm_flops(l, 1*n, l, 0) + gemm_flops(l, n, params.mb *r, 0) +
                                    params.mb * (
                                         gemm_flops(r, 1*n, r, 0) +
                                         (p * p * p / 3.0) /
                                         1000.0/1000.0/1000.0 )) ));

    AIOfile.finalize();

    delete []Ay_top;
    delete []Ay_bot;
    delete []AR;
    //delete []AL;
    delete []A;
    delete []Stl;
    delete []Str;
    delete []Sbr;
    delete []Ay;
    delete []S;
    delete []y_nan_idxs;
    delete []y_residual;
    delete []y_res_norms;
    // delete []backupAL;
    // delete []backupAR;
}

void Algorithm::hpc_statistics(list<long int>* indexs_AL,list<long int>* indexs_AR, list<long int>* indexs_Y, int n,
                type_precision* A, int a_amount, type_precision* y, int jj, type_precision* B, int p, int l,int r, type_precision* T, type_precision* R2,type_precision* P)
{
        type_precision tsx;
        type_precision tsxx;
        //type_precision tsx2;

        type_precision xb;
        type_precision sy;
        type_precision Syy;

        type_precision sym;
        type_precision st;
        type_precision SST;

        type_precision t;
        type_precision t1;

        type_precision* sab =  new type_precision[n];//used in res=y-y_fit


        for (int ii= 0; ii < a_amount; ii +=r)
        {

            type_precision* b = &B[jj*p*a_amount+ii*p];
            memset(sab,0,sizeof(type_precision)*n);



            for (int h=0; h< l; h++)
            {
            }

            list<long int> nans;//!decomment   = indexs_AR[ii].merge(indexs_Y[jj]);
            nans.unique();
            int n_not_nans = nans.size();

            type_precision* x = 0;//!remove old and decomment  = &AR[ii*n*p+h*n];


            for (int h=r; h< p; h++)
            {
                int b_idx = jj*p*a_amount+ii*p+h;
                tsx = 0;
                tsxx = 0;
                Syy = 0;

                for (int k=0; k < n; k++)
                {
                        tsx += x[k];
                        tsxx += x[k]*x[k];
                        xb = x[k]*b[h];
                        sab[k] -= xb;
                }

                for (int k=0; k < n; k++)
                {
                    sy=y[k]-sab[k];
                    Syy += sy*sy;
                }

                sym = 0;
                SST = 0;
                for (int k=0; k < n; k++)
                {
                    sym += y[k];
                }
                sym = sym/n_not_nans;

                for(int k=0; k < n; k++)
                {
                    st = y[k] - sym;
                    SST += st*st;
                }

                R2[b_idx] = (1-Syy)/SST;

                type_precision sig_res=sqrt(Syy/n);
                type_precision tsx2=tsx*tsx/n_not_nans;
                type_precision var_x=sqrt(tsxx-tsx2);

                t=b[h]*sig_res/var_x;
                //T[b_idx] =t;

                if(t < 1.28)
                {
                    t1=4.4-t;
                    P[b_idx]=0.5-0.1*t1*t;
                }
                else
                {
                    P[b_idx]=1-erf(t);
                }
            }


        }
        delete []sab;

        //t_students_cdf(y_amount,a_amount,p,T,P,n);

}

void Algorithm::t_students_cdf(int y_amount,int a_amount,int p, type_precision* T, type_precision* P, int deg_freedom)
{


}

//!*************************************************!//

//for reference
void Algorithm::partialNEQ_Blocked_STL(struct Settings params,
                                       struct Outputs &out)
{
//    int max_threads = params.threads;
//
//    srand(time(NULL));
//
//    blas_set_num_threads(max_threads);
//
//    type_precision *Ytemp;
//    lapack_int info, n, lda,ldy, l, r, p;
//
//    int i,k;
//
//    cputime_type start_tick, start_tick2, end_tick;
//
//    AIOfile.initialize(params);
//
//    n = params.n; l = params.l; r = params.r; p = l+r;
//
//    int y_amount = params.m;
//    int y_block_size = params.mb;  // kk
//
//    int a_amount = params.t;
//    int a_block_size = params.tb;
//
//    int a_iters = (a_amount + a_block_size - 1) / a_block_size;
//
//    lda = n; ldy = n;
//    k = p;
//
//
//
//    type_precision* backupAL;
//
//    AIOfile.load_AL(&backupAL);
//
//
//    type_precision Stl[l * l];
//    type_precision Str[l * r * a_block_size];
//
//
//    type_precision* Ay_top = new type_precision[l * y_amount];
//    type_precision* Ay_bot = new type_precision[y_block_size * a_block_size * r];
//    type_precision* A = new type_precision[n * p];
//    type_precision* AR = new type_precision[n * r * a_block_size];
//    type_precision* AL = new type_precision[n * l];
//
//    copy_vec(backupAL, AL, n * l);
//    copy_vec(backupAL, A, n * l);
//
//
//    int y_iters = y_amount/y_block_size;
//
//    type_precision* Y;
//    type_precision* backupAR;
//
//    // printf("\n\n%%Computations\n%%");
//
//
//    get_ticks(start_tick);
//
//
//    //! Generate S
//    cblas_ssyrk(CblasColMajor, CblasUpper, CblasTrans,
//                l, n, 1.0, backupAL, lda, 0.0, Stl, l);
//
//    for (int j = 0; j < y_iters; j++)
//    {
//        AIOfile.load_Yblock(&Y, y_block_size);
//
//        //! Ay_top = AL'*Y
//        cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
//                    l, y_block_size, n, 1.0, AL, n, Y, n, 0.0,
//                    &Ay_top[j * l * y_block_size], l);
//    }
//
//
//    for (int i = 0; i < a_iters; i++)
//    {
//        if (a_iters > 10 && (i%(a_iters/10)) == 0)
//        {
//            cout << "%" << flush;
//        }
//
//
//        AIOfile.load_ARblock(&backupAR, a_block_size);
//        AIOfile.reset_Y();
//        copy_vec(backupAR, AR, n*r*a_block_size);
//
//        // matlab_print_matrix("A", n, p, A);
//        //! Generate Str
//        cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
//                    l, r * a_block_size, n, 1.0, AL, n, AR, n, 0.0, Str, l);
//
//        type_precision Sbr[r*r*a_block_size];
//
//        #pragma omp parallel for default(shared) schedule(static)
//        for (int ii= 0; ii < a_block_size; ii++)
//        {
//            //! Generate Sbr
//            cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
//                        r, r, n, 1.0, &AR[ii*r*n], n, &AR[ii*r*n],
//                        n, 0.0, &Sbr[ii*r*r], r);
//        }
//
//        type_precision Ay[y_block_size*p];
//        type_precision S[p*p];
//
//
//        for (int j = 0; j < y_iters; j++)
//        {
//            if (a_iters < 10 && (y_iters < 10 || (j%(y_iters/10)) == 0))
//            {
//                cout << "%" << flush;
//            }
//
//
//            AIOfile.load_Yblock(&Y, y_block_size);
//
//            //! Ay_bot = AR'*Y
//            cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
//                        r * a_block_size, y_block_size, n, 1.0, AR, n, Y, n,
//                        0.0, Ay_bot, r * a_block_size);
//
//            #pragma omp parallel for private(S, Ay) default(shared) schedule(static)
//            for (int ii= 0; ii < a_block_size; ii++)
//            {
//                // matlab_print_matrix("Y", n, y_block_size, Y);
//
//                //! Rebuild AY
//
//                for (int jj = 0; jj < y_block_size; jj++)
//                {
//                    copy_vec(&Ay_top[j*l*y_block_size+jj*l], &Ay[jj*p], l);
//                    copy_vec(&Ay_bot[ii*r + jj*r*a_block_size], &Ay[l+jj*p], r);
//                }
//
//                type_precision* B = Ay;
//                // matlab_print_matrix("Ay_top", l, y_block_size,&Ay_top[j*l*y_block_size]);
//                // matlab_print_matrix("Ay_bot", r, y_block_size, Ay_bot);
//                // matlab_print_matrix("Ay", p, y_block_size, Ay);
//
//                //! Rebuild S
//                build_S(S, Stl, &Str[ii*r*l], &Sbr[ii*r*r], l, r);
//                // matlab_print_matrix("S", p, p, S);
//
//                //! b = S\Ay
//                info = LAPACKE_sposv(STORAGE_TYPE, 'U', p, y_block_size,
//                                     S, p, Ay, p);
//                // assert(info == 0,"POSV");
//
//
//                if (ForceCheck)
//                {
//                    check_result(backupAL, &AR[ii * r * n], n, p,
//                                 y_block_size, r, Y, B);
//                }
//            }
//        }
//    }
//
//    get_ticks(end_tick);
//    out.duration = ticks2sec(end_tick - start_tick, cpu_freq);
//    // out.gflops = a_amount/1000.0*(n*l*r+n*r*r+y_amount*(n*r+p*p*p)))/1000.0/1000.0;
//    out.gflops = gemm_flops(l, n, l, 0) + gemm_flops(l, n, y_amount, 0) +
//                 a_amount * (gemm_flops(l, n, r, 0) +
//                             gemm_flops(r, n, r, 0) +
//                             gemm_flops(r, n, y_amount, 0) +
//                             (y_amount/1000.0) *
//                             (p * p * p / 3.0) / 1000.0/1000.0);
//
//
//    AIOfile.finalize();
//
//    delete []Ay_top;
//    delete []Ay_bot;
//    delete []A;
//    delete []AL;
//    delete []AR;
}


