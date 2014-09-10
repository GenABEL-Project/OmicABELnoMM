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
                             type_precision* y, type_precision* res,struct Settings params,int iX,int iiX, int jY, int jjY)
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

    //matlab_print_matrix("A", rowsA, colsA, A);


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



    if(!params.use_fake_files)
    {
        //cout << ".";
        FILE* fp_Bprecomputed = fopen("examples/bpre.fvd", "rb");
        if(fp_Bprecomputed == 0)
        {
            cout << "Error Reading File of precomputed B values " << endl;
        }
        else
        {
            type_precision precomp_betas[colsA];

            int file_pos = jY*params.tb*params.m*colsA + jjY*params.m*colsA+iX*params.mb*colsA+iiX*colsA;
            fseek ( fp_Bprecomputed , file_pos*sizeof(type_precision) , SEEK_SET );


            size_t result = fread (precomp_betas,sizeof(type_precision),colsA,fp_Bprecomputed); result++;

            //matlab_print_matrix("res", colsA, 1, res);
            //matlab_print_matrix("new_sol", colsA, 1, new_sol);
            //matlab_print_matrix("precomp_betas", colsA, 1, precomp_betas);

            cblas_saxpy(colsA, -1.0, res, 1, precomp_betas, 1);
            type_precision u_norm = cblas_snrm2(colsA, precomp_betas, 1);
            //cout << "pa:"<< u_norm ;


            if (fabs(u_norm) > 0.001 || isnan(u_norm))
            {
                fseek ( fp_Bprecomputed , file_pos*sizeof(type_precision) , SEEK_SET );
                result = fread (precomp_betas,sizeof(type_precision),colsA,fp_Bprecomputed); result++;
                cblas_saxpy(colsA, 1.0, res, 1, precomp_betas, 1);
                type_precision u_norm2 = cblas_snrm2(colsA, precomp_betas, 1);
                //cout << u_norm2 << " ";
                //cout << file_pos <<"c "<< u_norm2 << "\n";

                if(fabs(u_norm2) > 0.001 || isnan(u_norm))
                {
                    fseek ( fp_Bprecomputed , file_pos*sizeof(type_precision) , SEEK_SET );
                    result = fread (precomp_betas,sizeof(type_precision),colsA,fp_Bprecomputed); result++;
                    fflush(stdout);

//                    matlab_print_matrix("AL", rowsA, colsA-colsAR, AL);
//                    matlab_print_matrix("AR", rowsA, colsAR, AR);
//                    matlab_print_matrix("Y", rowsA, rhs, y);
                    matlab_print_matrix("res", colsA, 1, res);
                    matlab_print_matrix("precomp_betas", colsA, 1, precomp_betas);
                    cout << file_pos <<" "<< u_norm << "\n";
                    cout << file_pos <<" "<< u_norm2 << "\n";
                    //printf("\n%%\tBeta computed nrom: %0.2g\n", u_norm2);
                    exit(1);
                }

            }

            fseek ( fp_Bprecomputed , file_pos*sizeof(type_precision) , SEEK_SET );
            result = fread (precomp_betas,sizeof(type_precision),colsA,fp_Bprecomputed); result++;
            cblas_saxpy(colsA, -1.0, new_sol, 1, precomp_betas, 1);
            u_norm = cblas_snrm2(colsA, precomp_betas, 1);
            //cout << " pl:" << u_norm ;
            if (fabs(u_norm) > 0.001 || isnan(u_norm))
            {
                cout << "lapack betas do not match!" << endl;
                exit(1);
            }
        }


        fclose(fp_Bprecomputed);
    }
    //else
    {


    //    if (PRINT)
    //        printf("\n Btop=(Rtl\\(Ql'*Y))-(Rtl\\Rtr)*(Rbr\\(Qr'*Y)); \n [Btop ; Rbr\\Qr'*Y] - bcomputed \n");

        cblas_saxpy(rhs * colsA, -1.0, res, 1, new_sol, 1);
        type_precision u_norm = cblas_snrm2(rhs * colsA, new_sol, 1);
        //
        if (fabs(u_norm) >= 0.1 || isnan(u_norm))
        {
            fflush(stdout);
           // matlab_print_matrix("AL", rowsA, colsA-colsAR, AL);
            //matlab_print_matrix("AR", rowsA, colsAR, AR);
            //matlab_print_matrix("A", rowsA, colsA, A);
            //matlab_print_matrix("Y", rowsA, rhs, y);
//            printf("\nA = [AL AR]; [Q, R] = qr(A, 0); rr = R\\(Q'*Y)\n");
            matlab_print_matrix("bcomputed", colsA, rhs, res);
            matlab_print_matrix("newsol", colsA, rhs, ynew);
            printf("\n%%\tnrom: %0.5g\n\n", u_norm);
             exit(1);
        }
        else
        {
            //matlab_print_matrix("bcomputed", colsA, rhs, res);
            //matlab_print_matrix("newsol", colsA, rhs, ynew);
           //printf("%%%0.3g \n", u_norm);
           //matlab_print_matrix("AL", rowsA, colsA-colsAR, AL);
//            matlab_print_matrix("AR", rowsA, colsAR, AR);

            //matlab_print_matrix("Y", rowsA, rhs, y);

           //cout << " la:" << u_norm << "\n";
        }




    }


    // cout << "\t**************";
    free(ynew);
    free(new_sol);
    free(A);
}


void Algorithm::applyDefaultParams(struct Settings &params)
{

    params.fname_excludelist = "";
    params.ForceCheck = false;
    params.use_fake_files = false;
    params.disp_cov = false;
    params.storePInd = false;
    params.storeBin = false;
    params.dosages = false;
    params.threads = 1;
    params.r = 1;
    params.model = -1;

    params.use_interactions = false;
    params.keep_depVar = true;
    params.expand_depVar = false;
    params.fname_interactions = "";

    params.minR2store = 0.00001;
    params.minR2disp = 0.000001;

    params.minPstore = 0.1;
    params.minPdisp = 0.05;


}

///////////////////////////////


void Algorithm::partialNEQ_Blocked_STL_MD(struct Settings params,
                                          struct Outputs &out)
{



    srand(time(NULL));



    //type_precision *Ytemp;
    lapack_int info, n, lda, l, r, p;

    cputime_type start_tick, start_tick2, start_tick3, end_tick;


    if(params.minR2disp > params.minR2store || params.storeBin)
        params.minR2store = params.minR2disp;

    if(params.minPdisp > params.minPstore || params.storeBin)
        params.minPstore = params.minPdisp;

    AIOwrapper AIOfile;//leave here to avoid memory errors of reusing old threads
    AIOfile.initialize(params);//THIS HAS TO BE DONE FIRST! ALWAYS

    //cout << params.n <<  "\n";
    //parameters read on AIOfile.initialize are then copied localy
    n = params.n; l = params.l; r = params.r; p = l+r;
    disp_cov = params.disp_cov;
    storePInd = params.storePInd;

    minR2store = params.minR2store;
    //params.minR2disp = params.minR2disp;//passed to AIOwrapper

    minTstore = getTvalue(params.minPstore);


    max_threads = params.threads;
    blas_set_num_threads(max_threads);
    omp_set_num_threads(max_threads);


    int y_amount = params.t;
    int y_block_size = params.tb;  // kk
    //cout << "yt:"<< y_amount << " oybz:"<<y_block_size << flush;

    int a_amount = params.m;
    int a_block_size = params.mb;

    int a_iters = (a_amount + a_block_size - 1) / a_block_size;

    int y_iters = (y_amount + y_block_size - 1) / y_block_size;

    //cout << "yiters:" <<  y_iters << " aiters:" << a_iters << endl;


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

    list<long int>* al_nan_idxs = new list<long int>[l];
    list<long int>* y_nan_idxs = new list<long int>[y_block_size];
    list<long int>* ar_nan_idxs = new list<long int>[a_block_size*r];
    int* Ymiss=new int[y_block_size];


    type_precision* A = new type_precision[n * p * 1];
    type_precision* AR = new type_precision[n * r * a_block_size * 1];

    //Sum of squares for A parts
    type_precision* ssAR = new type_precision[p*r*a_block_size];
    type_precision* ssAL = new type_precision[p*l];
    type_precision* ssY = new type_precision[y_block_size];

    SYY = new type_precision[a_block_size];


    list < resultH >* sigResults;



//  type_precision* AL = new type_precision[n * l * 1];
    type_precision* AL = A;

    type_precision* B = Ay;

    type_precision* backupAR;  // = new type_precision[n*r*a_block_size];
    type_precision* backupAL;  // = new type_precision[n*l];


    AIOfile.load_AL(&backupAL);

    //pthread_barrier_wait(&(AIOfile.Fhandler->finalize_barrier));

    //matlab_print_matrix("AL", n,l,backupAL);

//    replace_nans(al_nan_idxs,1, backupAL, n, l);
//    for (int i = 1; i < l; i++)
//    {
//        al_nan_idxs[i]=al_nan_idxs[0];
//    }


    sumSquares(backupAL,l,n,ssAL,0);


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
        //cout << "ybz:"<< y_block_size << " " << flush;

        get_ticks(end_tick);
        out.acc_loady += ticks2sec(end_tick,start_tick2);

        get_ticks(start_tick2);

        replace_nans(&y_nan_idxs[0],y_block_size, Y, n,1);
        sumSquares(Y,y_block_size,n,ssY,y_nan_idxs);

        for (int jj = 0; jj < y_block_size; jj++)
        {
//            list<long int> nans = al_nan_idxs[0];
//            list<long int> nans2 = y_nan_idxs[jj];
//            nans.merge(nans2);
//            nans.sort();//!might not be needed
//            nans.unique();
            Ymiss[jj] = y_nan_idxs[jj].size();
        }


        //matlab_print_matrix("Y", n, 1, Y);
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

            replace_nans(ar_nan_idxs, a_block_size*r, backupAR, n , 1);

            sumSquares(backupAR,a_block_size,n,ssAR,ar_nan_idxs);

            replace_nans_avgs(a_block_size, backupAR, n, r, ar_nan_idxs);

            //replace_with_zeros(al_nan_idxs, backupAR,  n, r, a_block_size);



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


                get_ticks(start_tick2);
                copy_vec(backupAR,AR, n*r*a_block_size);//!10%//try to remove!
                replace_with_zeros(&y_nan_idxs[jj], AR,  n, r, a_block_size);
                get_ticks(end_tick);
                out.acc_other += ticks2sec(end_tick,start_tick2);

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
//                    //cout << omp_get_thread_num() << endl << flush;
//
//                    get_ticks(start_tick2);
//
//                    copy_vec(&backupAR[ii*r*n],&AR[ii*r*n], n*r);//!10%//try to remove!
//                    replace_with_zeros(&y_nan_idxs[jj], &AR[ii*r*n],  n, r, 1);
//
//                    get_ticks(end_tick);
//                    out.acc_other += ticks2sec(end_tick,start_tick2);


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
                            //replace_with_zeros(al_nan_idxs, &Y[jj*n],  n, 1, 1);
                            check_result(AL, &AR[ii*r*n], n, p,
                                         1, r, &Y[jj*n], &Ay[ii*p],params,i,ii,j,jj);
                        }
                    }
                }
                }

                get_ticks(start_tick2);

                AIOfile.getCurrentWriteBuffers(sigResults);

                get_ticks(end_tick);
                out.acc_storeb += ticks2sec(end_tick,start_tick2);


                get_ticks(start_tick2);
                hpc_statistics( Ymiss[jj],n,AL,AR,a_block_size,&Y[jj*n],jj,ssY[jj],B,p,l,r,ssAL,ssAR,j*params.tb, i*params.mb*r, sigResults);
                get_ticks(end_tick);
                out.acc_stats += ticks2sec(end_tick,start_tick2);
                /**************************/

                blas_set_num_threads(max_threads);

                get_ticks(start_tick2);

                AIOfile.write_OutFiles(sigResults);


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


    out.gflops = gemm_flops(l, n, params.t*1, 0) + gemm_flops(params.m*r, n, params.t*1, 0)+
                             params.t *(  gemm_flops(l, 1*n, l, 0) + gemm_flops(l, n, params.m *r, 0) +
                                    params.m * (
                                         gemm_flops(r, 1*n, r, 0) +
                                         (p * p * p / 3.0) /
                                         1000.0/1000.0/1000.0 ))   ;

     out.gflops += params.t/1000.0*params.m/1000.0*((params.n)*(params.p*2+2))/1000.0;//hpcstatistics

     //cout << (params.t/1000.0*params.m/1000.0*(params.n*(params.p*2+2))/1000.0)/out.acc_stats;


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

    delete []SYY;


}

void Algorithm::hpc_SSY(int n,int p, int l,int r, type_precision* __restrict  AL, type_precision* __restrict  AR, int a_amount, type_precision* __restrict y, type_precision* __restrict  B )
{


    float* __restrict xl;

    float Xl_n[l*n];
    float yred[n];
    int n_idx = 0;
    int xl_n_idx = 0;
    for (int h=0; h < l; h++)
    {
        for (int k=0; k < n; k++)
        {
            xl = &AL[n*h];
            if(AL[k]!=0 && y[k] !=0)
            {
                if(h==0)
                {
                    yred[n_idx] =  y[k];
                    n_idx++;
                }
                Xl_n[xl_n_idx] = xl[k];
                xl_n_idx++;
            }
        }
    }

//        matlab_print_matrix("yred",1,n_idx,yred);
//        matlab_print_matrix("Y",1,n,y);

//    matlab_print_matrix("Xl",1,n-n/2,&AL[n/2*(l-1)]);
//    matlab_print_matrix("Xl_n",1,n_idx-n_idx/2,&Xl_n[n_idx/2*(l-1)]);

    //cout << n_idx << " ";

    //#pragma omp parallel for schedule(static) default(shared)
    for (int ii= 0; ii < a_amount; ii++)
    {
        int ar_idx = ii*(n*r);

        float* __restrict b = &B[ii*p];

        float* __restrict xr;
        float syvec[n_idx];
        float* __restrict sy = syvec;



        float XR_n[r*n_idx];

        int xr_idx=0;
        for (int h=0; h < r; h++)
        {
            xr = &AR[ar_idx+n*h];
            for (int k=0; k < n; k++)
            {
                if(AL[k]!=0 && y[k] !=0)
                {
                    XR_n[xr_idx] = xr[k];
                    xr_idx++;
                }
            }
        }

//
//    matlab_print_matrix("XR",1,n,&AR[ar_idx+n*(r-1)]);
//
//    matlab_print_matrix("XR_n",1,n_idx,&XR_n[n_idx*(r-1)]);

//        for (int k=0; k < n_idx; k++)
//        {
//            sy[k] = yred[k];
//        }
        memcpy(sy,yred,n_idx*sizeof(float));
//        matlab_print_matrix("yred",1,n_idx,yred);
//        matlab_print_matrix("SY",1,n_idx,sy);


        SYY[ii] = 0;

        for (int h=0; h < l; h++)
        {
            xl = &Xl_n[n_idx*h];

            for (int k=0; k < n_idx; k++)
            {
                sy[k] -= xl[k]*b[h];
            }
        }

        for (int h=l; h < p; h++)
        {
            xr = &XR_n[n_idx*(h-l)];

            for (int k=0; k < n_idx; k++)
            {
                sy[k] -= xr[k]*b[h];
            }
        }


        for (int k=0; k < n_idx; k++)
        {
            SYY[ii] += sy[k]*sy[k];
        }

    }

}


void Algorithm::hpc_statistics(int Ymiss, int n,
                type_precision* __restrict  AL, type_precision* __restrict  AR, int a_amount, type_precision* __restrict y, int jj, float varY, type_precision* __restrict  B,
                int p, int l,int r,type_precision* __restrict var_xL,type_precision* __restrict var_xR, int y_blck_offset, int A_blck_offset, list < resultH >* __restrict sigResults)
{



        int start_h=l;
        if(disp_cov)
            start_h = 1;

        type_precision SST= varY;


        sigResults->clear();

        int n_not_nans = Ymiss;
        //cout << " " <<  nans.size() << "\t";
        int n_corrected = n-n_not_nans-p;
        float varFactor = (float)(n_corrected+p)/(float)n;


        resultH res[max_threads];
        for(int i= 0; i < max_threads; i++)
        {
            res[i].nUsed = n_corrected+p;
            res[i].Y_name_idx = y_blck_offset+jj;
            res[i].nUsedPct = varFactor;
            res[i].ARoffset = A_blck_offset;
        }


        //!****************************

        hpc_SSY(n,p,l,r,AL,AR,a_amount,y,B);

        //!***************************

        #pragma omp parallel for schedule(static) default(shared)
        for (int ii= 0; ii < a_amount; ii++)
        {

            int thread_id = omp_get_thread_num();

            type_precision Syy;

            long double ptemp;
            float SE;

            long double t;
            type_precision t1;

            type_precision* b = &B[ii*p];


            Syy=SYY[ii];
            //cout << Syy << " ";

            //R2[ay_idx] = 1-(Syy/(SST*varFactor));
            res[thread_id].R2 = 1-((Syy/(n_corrected-1))/(SST/(n-1)));
            //R2 = SReg/(SST*varFactor);

            //cout << n_corrected <<" "<< Syy/(n_corrected-p-1) << " " << SST/(n-1) << " " << R2[ay_idx] << endl;
            type_precision sig_res=sqrt(Syy/(n_corrected));

            //res.R2 = R2;

            bool hasTminStorage = false;

            res[thread_id].res_p.clear();

            //cout << minR2store << " " << res.R2 << endl;

            for (int h=start_h; h < p && ( res[thread_id].R2 > minR2store || storePInd ); h++)
            {
                //cout << " " << sig_res<< " " <<  var_xR[h] << endl;
                if(h < l)
                {
                    SE = sig_res/sqrt(var_xL[h]*varFactor);
                    //cout << " " << sig_res<< " " <<  var_xL[h] << " "<< varFactor <<  endl;
                }
                else
                {
                    SE = sig_res/sqrt(var_xR[h]*varFactor);
                    //cout << " " << sig_res<< " " <<  var_xR[h] << " "<< varFactor << endl;
                }


                t=fabs(b[h]/SE);
                //cout << t << " ";



                if( t > minTstore || storePInd)
                {

                    //cout << R2 << " " << Syy << " " << SST <<  "\t:\t";
                     //cout << (Syy/(n_corrected-1)) << " " << (SST/(n-1)) << endl;
//                     if(n>900)
//                        cout << jj<< ":" << ii << " " << b[h] << " " <<  R2 << " " <<  SReg <<  " " << SST*varFactor <<" " << 1.0-((Syy/(n_corrected-1))/(SST*varFactor/(n-1))) <<  " "<< Syy << " " << t << endl;

                    hasTminStorage = true;

                    if( t < 1.6)
                    {
                        t1=4.4-t;
                        ptemp=0.5-0.1*t1*t;
                    }
                    else
                    {
                        ptemp=erfc(t*one_oversqrt2);
                    }

                    result rdata;

                    rdata.B = b[h];
                    rdata.T = t;
                    rdata.P = ptemp;
                    rdata.SE = SE;


                    rdata.AR_name_idx = ii*r;
                    rdata.AL_name_idx = h;

                    if(h >= l)
                    {
                        rdata.AR_name_idx += h-l;
                        rdata.AL_name_idx = -1;
                    }


                    res[thread_id].res_p.push_back(rdata);

                }
            }

            if(hasTminStorage)
            {
                #pragma omp critical
                {
                    sigResults->push_back(res[thread_id]);
                }
            }

        }





}

void Algorithm::sumSquares(type_precision* Data, int cols, int rows, type_precision* ssData, list<long int>* indexs_Data)
{

    for (int h=0; h< cols; h++)
    {
        double sx = 0;
        double sxx = 0;
        type_precision* x = &Data[h*rows];

        int n = rows;
        if(indexs_Data)
            n =n-indexs_Data[h].size();

        for (int k=0; k < rows; k++)
        {
                sx += x[k];
                sxx += x[k]*x[k];
        }
        ssData[h] = (type_precision)(sxx - sx*sx/n);
    }

}

//void Algorithm::mean(type_precision* Data, int block, int cols, int rows, type_precision* ssData, list<long int>* indexs_Data)
//{
//    for (int b=0; b< block; b++)
//    {
//        for (int h=0; h< cols; h++)
//        {
//            double sx = 0;
//            double sxx = 0;
//            type_precision* x = &Data[h*rows];
//
//            int n = rows;
//            if(indexs_Data)
//                n =n-indexs_A[b].size();
//
//            for (int k=0; k < rows; k++)
//            {
//                    sx += x[k];
//            }
//            ssData[h] = (type_precision)(sxx - sx*sx/n);
//        }
//    }
//}




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


