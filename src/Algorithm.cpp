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


void Algorithm::prepare_Bfinal(type_precision* bfinal, type_precision* top,
                               type_precision* bot, int dim1_b,
                               int dim2_b, int dim1_b_bot)
{
    // memcpy are faster version of the fors
    int i, k, w, top_idx, bot_idx;
    int size;
    top_idx = 0;
    bot_idx = 0;
    for (k = 0; k < dim2_b; k++)
    {
        size = k * dim1_b + (dim1_b - dim1_b_bot) - (k * dim1_b);
        memcpy( (type_precision*)&bfinal[k * dim1_b],
                (type_precision*)&top[top_idx],
                size * sizeof(type_precision) );
//        for (i = k*dim1_b; i < k*dim1_b+(dim1_b-dim1_b_bot); i++)
//        {
//            bfinal[i] = top[top_idx];
//            top_idx++;
//        }
        top_idx += size;
        i = k * dim1_b + size;
        w = i;

        size = w + dim1_b_bot - w;
        memcpy( (type_precision*)&bfinal[w],
                (type_precision*)&bot[bot_idx],
                size * sizeof(type_precision) );
//        for (i = w; i < w+dim1_b_bot; i++)
//        {
//            bfinal[i] = bot[bot_idx];
//            bot_idx++;
//        }
        bot_idx += size;
    }
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
    assert(info == 0, "Error Check");


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


    //type_precision *Ytemp;
    lapack_int info, n, lda, l, r, p;

    cputime_type start_tick, start_tick2, end_tick;

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
        if (y_iters >= 40 && (j%(y_iters/40)) == 0)
        {
            cout << "*" << flush;
        }

        for (int i = 0; i < a_iters; i++)
        {
            for (int jj = 0; jj < y_block_size; jj++)
            {
                if (y_iters < 40 &&
                    (y_block_size < 3 || (jj%(y_block_size/3)) == 0)
                    )
                {
                    cout << "*" << flush;
                }
            }
        }
    }

    if(!params.ForceCheck)
        cout << endl;



    //type_precision Stl[l*l];
    //type_precision Str[l*r*a_block_size];
    type_precision* Stl = new type_precision[l*l];
    type_precision* Str = new type_precision[l*r*a_block_size];

    type_precision* Sbr = new type_precision[r *  r * a_block_size];
    type_precision* Ay = new type_precision[p * a_block_size];

    type_precision* S = new type_precision[p * p];

    type_precision* Ay_top = new type_precision[l * y_amount];
    type_precision* Ay_bot = new type_precision[y_block_size * a_block_size * r];

    list<long int>* y_nan_idxs = new list<long int>[y_block_size];

    type_precision* AR = new type_precision[n * r * a_block_size * 1];
    type_precision* AL = new type_precision[n * l * 1];
    type_precision* backupAR;  // = new type_precision[n*r*a_block_size];
    type_precision* backupAL;  // = new type_precision[n*l];


    AIOfile.load_AL(&backupAL);
    //int total_al_nans = replace_nans(0, backupAL, n, l);
    replace_nans(0, backupAL, n, l);

    copy_vec(backupAL, AL, n*l);


    type_precision* Y;

    // printf("\n\n%%Computations\n%%");


    get_ticks(start_tick);

    for (int j = 0; j < y_iters; j++)
    {
        if (y_iters >= 40 && (j%(y_iters/40)) == 0 && !params.ForceCheck)
        {
            cout << AIOfile.io_overhead << flush;
            AIOfile.io_overhead = "*";
        }

        AIOfile.load_Yblock(&Y, y_block_size);



        //int total_y_nans = replace_nans(&y_nan_idxs[0], Y, n, y_block_size);
        replace_nans(&y_nan_idxs[0], Y, n, y_block_size);


        //! Ay_top = AL'*Y
        cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    l, y_block_size, n, 1.0, AL, n, Y, n, 0.0,
                    &Ay_top[j * l * y_block_size], l);


        for (int i = 0; i < a_iters; i++)
        {
            AIOfile.load_ARblock(&backupAR, a_block_size);
            //int total_ar_nans = replace_nans(0, backupAR, n, a_block_size * r);
            replace_nans(0, backupAR, n, a_block_size * r);
            copy_vec(backupAR, AR, n *  r * a_block_size);

            get_ticks(start_tick2);

            get_ticks(end_tick);
            out.acc_pre += ticks2sec(end_tick-start_tick2, cpu_freq);

            get_ticks(start_tick2);
            //! Ay_bot = AR'*Y
            cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                        r * a_block_size, y_block_size, n, 1.0, AR, n, Y, n,
                        0.0, Ay_bot, r * a_block_size);

            get_ticks(end_tick);
            out.firstloop += ticks2sec(end_tick - start_tick2, cpu_freq);

            //#pragma omp parallel default(shared)
            {
                for (int jj = 0; jj < y_block_size; jj++)
                {
                    int thread_id = 0 * omp_get_thread_num();
                    int aL_idx = thread_id * l * n;
                    int aR_idx = thread_id * r * n * a_block_size;

                    if (y_iters < 40 &&
                                  (y_block_size < 3 ||
                                                  (jj%(y_block_size/3)) == 0) && !params.ForceCheck)
                    {
                        cout << AIOfile.io_overhead << flush;
                        AIOfile.io_overhead = "*";
                    }

                    copy_vec(backupAL, &AL[aL_idx], n * l);

                    replace_with_zeros(&y_nan_idxs[jj], &AL[aL_idx], n, l, 1);
                    //! Generate Stl
                    cblas_ssyrk(CblasColMajor, CblasUpper, CblasTrans,
                                l, n, 1.0, &AL[aL_idx], lda, 0.0, Stl, l);

                    copy_vec(backupAR,&AR[aR_idx], n*r*a_block_size);
                    replace_with_zeros(&y_nan_idxs[jj], &AR[aR_idx],
                                       n, r, a_block_size);

                    get_ticks(start_tick2);
                    //! Generate Str
                    cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                                l, r * a_block_size, n, 1.0, &AL[aL_idx],
                                n, &AR[aR_idx], n, 0.0, Str, l);
                    get_ticks(end_tick);
                    out.acc_RTL_QLY += ticks2sec(end_tick - start_tick2,
                                                 cpu_freq);


                    //type_precision Sbr[r *  r * a_block_size];
                    //type_precision Ay[p * a_block_size];


                    //#pragma omp for nowait  schedule(dynamic)
                    for (int ii= 0; ii < a_block_size; ii++)
                    {
                        get_ticks(start_tick2);
                        //! Generate Sbr
                        cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                                    r, r, n, 1.0, &AR[aR_idx+ii*r*n], n,
                                    &AR[aR_idx + ii * r * n], n, 0.0,
                                    &Sbr[ii * r * r], r);
                        get_ticks(end_tick);
                        out.acc_gemm += ticks2sec(end_tick - start_tick2,
                                                  cpu_freq);

                        get_ticks(start_tick2);

                        copy_vec(&Ay_top[j*l*y_block_size+jj*l], &Ay[ii*p], l);
                        copy_vec(&Ay_bot[ii*r + jj*r*a_block_size],
                                 &Ay[l+ii*p], r);


                        //type_precision* B = Ay;
                        //type_precision S[p * p];


                        //! Rebuild S
                        build_S(S, Stl, &Str[ii*r*l], &Sbr[ii*r*r], l, r);
                        // matlab_print_matrix("S", p, p, S);

                        get_ticks(end_tick);
                        out.acc_b += ticks2sec(end_tick - start_tick2,
                                               cpu_freq);

                        get_ticks(start_tick2);

                        //! b = S\Ay
                        info = LAPACKE_sposv(STORAGE_TYPE, 'U', p, 1, S, p,
                                             &Ay[ii*p], p);

                        get_ticks(end_tick);
                        out.acc_loadxr += ticks2sec(end_tick - start_tick2,
                                                    cpu_freq);
                        assert(info == 0, "POSV");

                        if (params.ForceCheck)
                        {
                            #pragma omp critical
                            {
                                check_result(AL, &AR[aR_idx+ii*r*n], n, p,
                                             1, r, &Y[jj*n], &Ay[ii*p]);
                            }
                        }
                    }

                    AIOfile.write_B(Ay, p, a_block_size);
                }
            }
        }
         AIOfile.reset_AR();
    }

    get_ticks(end_tick);
    out.duration = ticks2sec(end_tick - start_tick, cpu_freq);
    // out.gflops = a_amount/1000.0*(n*l*r+n*r*r+y_amount*(n*r+p*p*p)))/1000.0/1000.0;
    out.gflops = y_amount * (gemm_flops(l, n, 1, 0) +
                             a_amount * (gemm_flops(r, n, 1, 0) +
                                         gemm_flops(l, n, l, 0) +
                                         gemm_flops(l, n, r, 0) +
                                         gemm_flops(r, n, r, 0) +
                                         (p * p * p / 3.0) /
                                         1000.0/1000.0/1000.0));

    AIOfile.finalize();

    delete []Ay_top;
    delete []Ay_bot;
    delete []AR;
    delete []AL;
    delete []Stl;
    delete []Str;
    delete []Sbr;
    delete []Ay;
    delete []S;
    delete []y_nan_idxs;
    // delete []backupAL;
    // delete []backupAR;
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


