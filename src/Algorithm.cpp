#include "Algorithm.h"



/*!
 * \brief All initialization happens in the resective solver method.


 *
 */
Algorithm::Algorithm()
{
    // ctor

}



/*!
 * \brief All cleaning happens in the respective solver method.


 *
 */
Algorithm::~Algorithm()
{
    // dtor
}


/*!
 * \brief In the presence of more than one algorothm for solving the LSQ, coose the resective one.


 *
 */
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





/*!
 * \brief Rebuilds the SPD matrix S=X^T*X from its parts.


 * \param the parts that were computed separately and its dimensions.
 */
void Algorithm::build_S(float* S, float* Stl,
                        float* Str, float* Sbr, int l, int r)
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

/*!
 * \brief Used to check validity of a single OLS.


 * \param Parts of the matrix A (X), Vector Y, results (Coefficients B), dimensions.
 * \warning iX,iiX,jY,jjY and related data to precomputed results are for testing puroses and are broken atm.
 */
void Algorithm::check_result(float* AL, float* AR,
                             int rowsA, int colsA, int rhs, int colsAR,
                             float* y, float* res,struct Settings params,int iX,int iiX, int jY, int jjY)
{
    float* A = (float*)malloc(rowsA * colsA *
                                                sizeof(float));

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


    float* ynew = replicate_vec(y, rowsA*rhs);
    float* new_sol = (float*)malloc(colsA * rhs *
                                                      sizeof(float));

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

    float* precomp_betas = new float[colsA];

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
            //float precomp_betas[colsA];

            int file_pos = jY*params.tb*params.m*colsA + jjY*params.m*colsA+iX*params.mb*colsA+iiX*colsA;
            fseek ( fp_Bprecomputed , file_pos*sizeof(float) , SEEK_SET );


            size_t result = fread (precomp_betas,sizeof(float),colsA,fp_Bprecomputed); result++;

            //matlab_print_matrix("res", colsA, 1, res);
            //matlab_print_matrix("new_sol", colsA, 1, new_sol);
            //matlab_print_matrix("precomp_betas", colsA, 1, precomp_betas);

            cblas_saxpy(colsA, -1.0, res, 1, precomp_betas, 1);
            float u_norm = cblas_snrm2(colsA, precomp_betas, 1);
            //cout << "pa:"<< u_norm ;


            if (fabs(u_norm) > 0.001 || isnan(u_norm))
            {
                fseek ( fp_Bprecomputed , file_pos*sizeof(float) , SEEK_SET );
                result = fread (precomp_betas,sizeof(float),colsA,fp_Bprecomputed); result++;
                cblas_saxpy(colsA, 1.0, res, 1, precomp_betas, 1);
                float u_norm2 = cblas_snrm2(colsA, precomp_betas, 1);
                //cout << u_norm2 << " ";
                //cout << file_pos <<"c "<< u_norm2 << "\n";

                if(fabs(u_norm2) > 0.001 || isnan(u_norm))
                {
                    fseek ( fp_Bprecomputed , file_pos*sizeof(float) , SEEK_SET );
                    result = fread (precomp_betas,sizeof(float),colsA,fp_Bprecomputed); result++;
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

            fseek ( fp_Bprecomputed , file_pos*sizeof(float) , SEEK_SET );
            result = fread (precomp_betas,sizeof(float),colsA,fp_Bprecomputed); result++;
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
        float u_norm = cblas_snrm2(rhs * colsA, new_sol, 1);
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
    delete []precomp_betas;
}


/*!
 * \brief Initializes parameters passed to the solvers.

    Should be called before the usage of any paramaters and
    must be consistent with the default behaivior in the absense of
    any other user specified changes.

 * \param parameters to initlialize.
 */
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
    params.mpi_id = 0;
    params.mpi_num_threads = 1;

    params.use_interactions = false;
    params.keep_depVar = false;
    params.use_multiple_interaction_sets = false;
    params.fname_interactions = "";
    params.ini_IDX_interactions = -1;
    params.end_IDX_interactions = -1;

    params.minR2store = 0.00001;
    params.minR2disp = 0.000001;

    params.minPstore = 0.1;
    params.minPdisp = 0.05;

    params.limit_m = INT_MAX;
    params.limit_t = INT_MAX;
    params.limit_n = INT_MAX;


}

///////////////////////////////


/*!
 * \brief Main OLS algorithm based on NEQ with Missing Data Suport.


 * \param Parameter, both default and/or user defined, Output structure to store non file results
 */
void Algorithm::partialNEQ_Blocked_STL_MD(struct Settings params,
                                          struct Outputs &out)
{



    srand(time(NULL));



    //float *Ytemp;
    lapack_int info, n, lda, l, r, p;

    cputime_type start_tick, start_tick2, start_tick3, end_tick;


    if(params.minR2disp < params.minR2store || params.storeBin)
        params.minR2store = params.minR2disp;

    if(params.minPdisp > params.minPstore || params.storeBin)
        params.minPstore = params.minPdisp;

    AIOwrapper AIOfile;//leave here to avoid memory errors of reusing old threads
    AIOfile.initialize(params);//THIS HAS TO BE DONE FIRST! ALWAYS


    total_results = 0;

    //cout << params.n <<  "\n";
    //parameters read on AIOfile.initialize are then copied localy
    n = params.n; l = params.l; r = params.r; p = l+r;
    disp_cov = params.disp_cov;
    storePInd = params.storePInd;

    minR2store = params.minR2store;
    //params.minR2disp = params.minR2disp;//passed to AIOwrapper

    minTstore = getTvalue(params.minPstore);
    minPdisp   = params.minPdisp;


    max_threads = params.threads;
    blas_set_num_threads(max_threads);
    omp_set_num_threads(max_threads);


    int y_amount = params.t;
    int y_block_size = params.tb;

    //cout << "yt:"<< y_amount << " oybz:"<<y_block_size << flush;

    int a_amount = params.m;
    int a_block_size = params.mb;


    //cout << r << endl;


    int a_iters = (a_amount + a_block_size - 1) / a_block_size;

    int y_iters = (y_amount + y_block_size - 1) / y_block_size;

    //cout << "yiters:" <<  y_iters << " aiters:" << a_iters << endl;


    lda = n;
    if(!params.ForceCheck && params.mpi_id == 0)
    {
        cout << endl;
    }

    for (int j = 0; j < y_iters && !params.ForceCheck && params.mpi_id == 0; j++)
    {
        if (!params.ForceCheck &&  ((y_iters < 10 && y_iters > 2 )|| (y_iters >= 10 && (j%(y_iters/10)) == 0 )) )
        {
            cout << "*" << flush;
        }

        for (int i = 0; i < a_iters; i++)
        {

                if ( !params.ForceCheck &&  y_iters <= 2 &&
                (  (a_iters >= 10 && (i%(a_iters/(10))) == 0) || (a_iters < (10)) ))
                {
                    cout << "*" << flush;
                }

        }
    }

    if(!params.ForceCheck && params.mpi_id == 0)
        cout << endl;



    float* Stl = new float[l*l*1];
    float* Str = new float[l*r*a_block_size*1];


    float* Stl_corr = new float[l*l*y_block_size];
    float* Sbr_corr = new float[r*r*a_block_size*y_block_size];
    float* Str_corr = new float[l*r*a_block_size*y_block_size];


    float* Sbr = new float[r *  r * a_block_size];
    float* Ay = new float[p * a_block_size * y_block_size];



    float* S2global = new float[p*p*max_threads];


    float* Ay_top = new float[l * y_amount];
    float* Ay_bot = new float[y_block_size * a_block_size * r];

    //float* y_residual = new float[n * y_block_size ];
    //float* y_res_norms = new float[a_block_size];

    //list<long int>* al_nan_idxs = new list<long int>[l];
    list<long int>* y_nan_idxs = new list<long int>[y_block_size];
    list<long int>* ar_nan_idxs = new list<long int>[a_block_size*r];
    int* Ymiss=new int[y_block_size];




    float* A = new float[n * p * 1];
    float* AR;// = new float[n * r * a_block_size * 1];

    //Sum of squares for A parts
    float* ssA = new float[p*a_block_size];

    float* ssY = new float[y_block_size];

    SYY = new float[a_block_size*y_block_size];


    //float* cholSAL = new float[l*l];


    list < resultH >* sigResults;



//  float* AL = new float[n * l * 1];
    float* AL = A;

    float* B = Ay;

    float* backupAR;  // = new float[n*r*a_block_size];
    float* backupAL;  // = new float[n*l];


    AIOfile.load_AL(&backupAL);



    //matlab_print_matrix("AL", n,l,backupAL);





    copy_vec(backupAL, AL, n*l);

    get_ticks(start_tick2);

    //! Generate Stl
    cblas_ssyrk(CblasColMajor, CblasUpper, CblasTrans,
                l, n, 1.0, AL, lda, 0.0, Stl, l);

    get_ticks(end_tick);
    out.acc_stl += ticks2sec(end_tick,start_tick2);


    //linear dep of columns check
    float* singular_vals = new float[l];
    float* u = new float[n*n];
    float* vt = new float[l*l];
    float* superb = new float[l-1];
//    for (int j = 0; j < n; j++)
//    {
//        AL[j] = AL[j+n];
//    }


    info = LAPACKE_sgesvd(STORAGE_TYPE, 'N','N',n,l,AL,n,singular_vals,u,n,vt,l,superb);info++;

    for (int j = 0; j < l; j++)
    {
        if(singular_vals[j] < 0.000001)
        {
            if(params.mpi_id == 0)
            {
                cout << endl;
                cout << "Error, the covariate matrix contains linearly dependent columns. Please review and change it accordingly." << singular_vals[j];
                cout << "Make sure there are no duplciates or closesly realted covariates (complement).";
            }
            exit(1);
        }
    }

    delete []singular_vals;
    delete []u;
    delete []vt;
    delete []superb;



    float* Y;

    // printf("\n\n%%Computations\n%%");

    copy_vec(backupAL, AL, n*l);


    get_ticks(start_tick);

    for (int j = 0; j < y_iters; j++)
    {
        if (!params.ForceCheck && params.mpi_id == 0 && ((y_iters < 10 && y_iters > 2 )|| (y_iters >= 10 && (j%(y_iters/10)) == 0 )) )
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
//            for (list< long int  >::iterator it=y_nan_idxs[jj].begin(); it != y_nan_idxs[jj].end(); ++it)
//            {
//                list< int > jjs = y_missings_jj[*it];
//                jjs.push_back(jj);
//                y_missings_jj[*it] = jjs;
//            }


            Ymiss[jj] = y_nan_idxs[jj].size();
        }


        //matlab_print_matrix("Y", n, y_block_size, Y);
        get_ticks(end_tick);
        out.acc_other += ticks2sec(end_tick,start_tick2);


        get_ticks(start_tick2);

        //! Ay_top = AL'*Y
        cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    l, y_block_size, n, 1.0, AL, n, Y, n, 0.0,
                    &Ay_top[j * l * y_block_size], l);

        get_ticks(end_tick);
        out.acc_gemm += ticks2sec(end_tick,start_tick2);


        get_ticks(start_tick2);
        generate_correction_missings_Stl(y_block_size,1,1,n,l,l,AL,AL,Stl_corr,y_nan_idxs);

        get_ticks(end_tick);
        out.acc_scorrect += ticks2sec(end_tick,start_tick2);


        for (int i = 0; i < a_iters; i++)
        {

            if (!params.ForceCheck && y_iters <= 2 && params.mpi_id == 0 &&
                (  (a_iters >= 10 && (i%(a_iters/(10))) == 0) || (a_iters < (10)) ))
            {
                cout << AIOfile.io_overhead << flush;
                AIOfile.io_overhead = "*";
            }

            get_ticks(start_tick2);
            //cout << "^"  << flush;
            AIOfile.load_ARblock(&backupAR, a_block_size);
           // cout << "^" << endl << flush;

            get_ticks(end_tick);
            out.acc_loadxr += ticks2sec(end_tick,start_tick2);

            get_ticks(start_tick2);

            replace_nans(ar_nan_idxs, a_block_size*r, backupAR, n , 1);

            //sumSquares(backupAR,a_block_size*r,n,ssAR,0);

            replace_nans_avgs(a_block_size, backupAR, n, r, ar_nan_idxs);

            get_ticks(end_tick);
            out.acc_other += ticks2sec(end_tick,start_tick2);

            AR = backupAR;


            get_ticks(start_tick2);

            //! Generate Str
            cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                                l, r * a_block_size, n, 1.0, AL,
                                            n, AR, n, 0.0, Str, l);//!

            get_ticks(end_tick);
            out.acc_str += ticks2sec(end_tick,start_tick2);


            get_ticks(start_tick2);

            blas_set_num_threads(1);
            omp_set_num_threads(max_threads);

            #pragma omp parallel default(shared)
            {

            #pragma omp for nowait
            for (int ii= 0; ii < a_block_size; ii++)
            {
                //! Generate Sbr
                cblas_sgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                            r, r, n, 1.0, &AR[ii*r*n], n,
                            &AR[ii * r * n], n, 0.0,
                            &Sbr[ii * r * r], r);

                float* S2 = &S2global[p*p*omp_get_thread_num()];

                //! build S for X'X-1 for fake "variance of x"
                build_S(S2, Stl, &Str[ii*r*l], &Sbr[ii*r*r], l, r);



                info = LAPACKE_spotrf(STORAGE_TYPE,'U',p,S2,p);info++;

                info = LAPACKE_spotri(STORAGE_TYPE,'U',p,S2,p);info++;

                for (int h= 0; h < p; h++)
                    ssA[ii*p+h] = S2[h*p+h];//diagonal elements ARE the variance


            }




            }

            blas_set_num_threads(max_threads);



            get_ticks(end_tick);
            out.acc_sbr += ticks2sec(end_tick,start_tick2 );

            get_ticks(start_tick2);

            generate_correction_missings_Str(y_block_size,a_block_size,n,l,r,AL,AR,Str_corr,Sbr_corr,y_nan_idxs);

            get_ticks(end_tick);
            out.acc_scorrect += ticks2sec(end_tick,start_tick2);



            //cout << "mb:" << a_block_size << " ";

            //matlab_print_matrix("ARb",n,a_block_size*r,backupAR);



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
                //cout << "y_block_size:" << y_block_size << " a_block_size:" << a_block_size <<" r:" << r << "\n";

                blas_set_num_threads(1);
                omp_set_num_threads(max_threads);

                #pragma omp parallel default(shared)
                {

                #pragma omp for nowait
                for (int ii= 0; ii < a_block_size; ii++)
                {
//                    //cout << omp_get_thread_num() << endl << flush;
                    cputime_type start_tick2;

                    get_ticks(start_tick2);

                    copy_vec(&Ay_top[j*l*y_block_size+jj*l], &Ay[jj*p*a_block_size+ ii*p], l);
                    copy_vec(&Ay_bot[ii*r + jj*r*a_block_size],
                             &Ay[jj*p*a_block_size+l+ii*p], r);


                    float* S2 = &S2global[p*p*omp_get_thread_num()];

                    //! Rebuild S
                    build_S(S2, Stl, &Str[ii*r*l], &Sbr[ii*r*r], l, r);

                    if(omp_get_thread_num()==0)
                    {
                    get_ticks(end_tick);//5%
                    out.acc_other += ticks2sec(end_tick,start_tick2);
                    }

                    get_ticks(start_tick2);
                    //!correct S from missings of Y
                    correct_missings_S(p,l,r, S2, &Stl_corr[jj*l*l], &Sbr_corr[jj*r*r*a_block_size+ii*r*r], &Str_corr[jj*l*r*a_block_size+ii*l*r]);

                    if(omp_get_thread_num()==0)
                    {
                    get_ticks(end_tick);
                    out.acc_inner_scorrect += ticks2sec(end_tick,start_tick2);
                    }


                    //matlab_print_matrix("S", p, p, S);

                    get_ticks(start_tick2);


                    //! b = S\Ay
                    info = LAPACKE_sposv(STORAGE_TYPE, 'U', p, 1, S2, p,
                                         &Ay[jj*p*a_block_size+ii*p], p);

                    //myassert(info == 0, "S\\Ay",info)
//                    if(info < 0)
//                        cout << info << " ";
                    info++;

                    if(omp_get_thread_num()==0)
                    {
                    get_ticks(end_tick);
                    out.acc_solve += ticks2sec(end_tick,start_tick2 );
                    }


                    if (params.ForceCheck)
                    {
                        #pragma omp critical
                        {
                            //replace_with_zeros(al_nan_idxs, &Y[jj*n],  n, 1, 1);
                            check_result(AL, &AR[ii*r*n], n, p,
                                         1, r, &Y[jj*n], &Ay[jj*p*a_block_size+ii*p],params,i,ii,j,jj);
                        }
                    }
                }
                }



                get_ticks(start_tick2);
                AIOfile.getCurrentWriteBuffers(sigResults);
                get_ticks(end_tick);
                out.acc_storeb += ticks2sec(end_tick,start_tick2);


                omp_set_num_threads(max_threads);
            }

            omp_set_num_threads(max_threads);
            //!****************************
            get_ticks(start_tick2);
            hpc_SSY(n,p,l,r,AL,AR,a_block_size,y_block_size,Y,B);

            for (int jj = 0; jj < y_block_size; jj++)
            {
                hpc_statistics( Ymiss[jj],n,AIOfile.realN,AL,AR,a_block_size,&Y[jj*n],jj,ssY[jj],&B[jj*a_block_size*p],p,l,r,ssA,j*params.tb, i*params.mb*r, sigResults);
            }
            get_ticks(end_tick);
            out.acc_stats += ticks2sec(end_tick,start_tick2);


            get_ticks(start_tick2);
            AIOfile.write_OutFiles(sigResults);
            get_ticks(end_tick);
            out.acc_storeb += ticks2sec(end_tick,start_tick2);


            blas_set_num_threads(max_threads);

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

     out.total_sig_results = total_results;

    {
    AIOfile.finalize();

    delete []Ay_top;
    delete []Ay_bot;
    delete []Ymiss;
    //delete []AR;
    //delete []AL;
    delete []A;
    delete []Stl;
    delete []Str;
    delete []Sbr;
    delete []Ay;
    //delete []S;
    delete []y_nan_idxs;
    delete []ar_nan_idxs;


    delete []Stl_corr;
    delete []Str_corr;
    delete []Sbr_corr;
    // delete []backupAL;
    // delete []backupAR;
    delete []ssY;
    delete []ssA;


    delete []SYY;
    delete []S2global;
    }


}

/*!
 * \brief Highperformance calculation of residuals between blocks of Y and Xb: Y-Xb for all n of all elements within a block


 * \param Matrix blocks, dimensions
 * \warning This function is 50% of the load of the entire runtime. It is memorybound at the moment.
 * Consider optimizing it with BLAS or any other options.
 */
void Algorithm::hpc_SSY(int n,int p, int l,int r, float*  AL, float*  ARb, int a_block_size,
                                    int y_block_size,
                                    float*  Yb, float*  Bb )
{
     float*  syvec = new float[n*max_threads];

    omp_set_num_threads(max_threads);

    int y_sub_block = min(100,y_block_size);
    int orig_y_sub_block = y_sub_block;
    int y_iters = (y_block_size+y_sub_block-1)/y_sub_block;


    int a_sub_block = min(25,a_block_size);
    int a_iters = (a_block_size+a_sub_block-1)/a_sub_block;
    int orig_a_sub_block = a_sub_block;

    //cout << "a_iters:" << a_iters << " y_iters:" << y_iters << " " << y_block_size<< endl;



    for (int j = 0; j < y_iters; j++)
    {
        float* Y = &Yb[j*orig_y_sub_block*n];
        float* By = &Bb[j*orig_y_sub_block*a_block_size*p];
        for (int i = 0; i < a_iters; i++)
        {
            float* AR = &ARb[i*orig_a_sub_block*r*n];
            float* Bxy = &By[i*orig_a_sub_block*p];

            y_sub_block = min(orig_y_sub_block, y_block_size-j*orig_y_sub_block);
            for (int jj = 0; jj < y_sub_block; jj++)
            {
                a_sub_block = min(orig_a_sub_block, a_block_size-i*orig_a_sub_block);

                #pragma omp parallel default(shared)
                {
                #pragma omp  for nowait
                for (int ii= 0; ii < a_sub_block; ii++)
                {
                    int thead_id = omp_get_thread_num();

                    float* sy = &syvec[thead_id*n];
                    float* B = &Bxy[jj*a_block_size*p+ii*p];

                    memset(sy,0,n*sizeof(float));

                    int ssy_idx = j*orig_y_sub_block*a_block_size   +   i*orig_a_sub_block    +    jj*a_block_size    +  ii;
//                    if(ssy_idx >= a_block_size*y_block_size)
//                    {
//                            cout << "err:"<< ssy_idx <<":"<< a_block_size*y_block_size << endl;
//                            cout << j*orig_y_sub_block*a_block_size << endl;
//                            cout << i*orig_a_sub_block << endl;
//                            cout << jj*a_block_size << endl;
//                            cout <<j <<" "  << i<<" "  << jj <<":" <<y_sub_block  <<" "  << ii<< ":" << a_sub_block << " "  <<endl;
//                            exit(1);
//                    }


                    SYY[ssy_idx]=singleSSY(B,&Y[jj*n], sy,AL,&AR[ii*n*r],l,r,p,n,thead_id );

                }
                }
            }
        }
    }


    delete []syvec;
}


/*!
 * \brief Highperformance calculation of residuals between single Y and Xb: Y-Xb for all n


 * \param matrix X, Vector Y and position
 * \warning This function is 50% of the load of the entire runtime. It is memorybound at the moment.
 * Consider optimizing it with BLAS or any other options.
 */
inline float Algorithm::singleSSY(float* __restrict b, float* __restrict y, float* __restrict sy,float* __restrict Xl_n,float* __restrict XR_n,
                     int l,int r, int p, int n, int thead_id  )
{

    long double sum = 0;

    float* __restrict xl;
    float* __restrict xr;




    for (int h=0; h < l; h++)
    {
        xl = &Xl_n[(n)*h];
        float bval = b[h];

        for (int k=0; k < n; k++)
        {
            sy[k] += xl[k]*bval;
        }
    }


    for (int h=l; h < p; h++)
    {
        xr = &XR_n[(n)*(h-l)];

        float bval = b[h];

        for (int k=0; k < n; k++)
        {
            sy[k] += xr[k]*bval;

        }
    }

    for (int k=0; k < n; k++)
    {
        if(y[k] != 0)
        {
            sy[k] -= y[k];
            sum += sy[k]*sy[k];
        }
    }

    return (float)sum;

}


/*!
 * \brief Highperformance calculation of statistics resulting from OLS

    Will calculate residuals,T statistics and Pvalues based on the available data of X, Y, variances and coefficients

 * \param TODO
 */
void Algorithm::hpc_statistics(int Ymiss, int n, int realN,
                float* __restrict  AL, float* __restrict  AR, int a_block_size, float* __restrict y, int jj, float varY, float* __restrict  B,
                int p, int l,int r,float* __restrict var_x, int y_blck_offset, int A_blck_offset, list < resultH >* __restrict sigResults)
{


        int start_h=l;
        if(disp_cov)
            start_h = 1;

        float SST= varY;




        int n_not_nans = Ymiss;
        //cout << " " <<  nans.size() << "\t";
        int n_corrected = n-n_not_nans-p;
        float varFactor = (float)(n_corrected+p)/(float)n;


        resultH* res = new resultH[max_threads];
        for(int i= 0; i < max_threads; i++)
        {
            res[i].nUsed = n_corrected+p;
            res[i].Y_name_idx = y_blck_offset+jj;
            res[i].nUsedPct = (float)(n_corrected+p)/(float)realN;
            res[i].ARoffset = A_blck_offset;
        }

//        cout << endl;
//        for (int ii= 0; ii < a_block_size*r; ii++)
//        {
//            cout << var_xR[ii] << " ";
//        }


        //#pragma omp parallel for schedule(static) default(shared)//fix
        for (int ii= 0; ii < a_block_size; ii++)
        {

            int thread_id = omp_get_thread_num();

            float Syy;

            long double ptemp;
            float SE;

            long double t;
            float t1;

            float* b = &B[ii*p];


            Syy=SYY[jj*a_block_size+ii];
            //cout << Syy << " ";

            //R2[ay_idx] = 1-(Syy/(SST*varFactor));
            res[thread_id].R2 = 1-((Syy/(n_corrected-1))/(SST/(n-1)));
            //R2 = SReg/(SST*varFactor);

            //cout << n_corrected <<" "<< Syy/(n_corrected-p-1) << " " << SST/(n-1) << " " << R2[ay_idx] << endl;
            if(Syy < 0.0)
            {
                cout << "Error with Syy < 0 \n";exit(1);
            }
            float sig_res=sqrt(Syy/(n_corrected));


            //cout << A_blck_offset << " " << ii;
//            for (int h=start_h; h < p ; h++)
//            {
//                //if(jj==(ii+A_blck_offset))
//                {
//                    //cout << "b:" << b[h] << " R2:" << res[thread_id].R2 << " sig:" << sig_res << " vxr:" <<  var_xR[ii*r+h-l]   << " SE" <<  SE <<  " t:" << t;
//                    if( b[h] < 0.0)
//                        cout << "b:" << b[h] << " ";
//                }
//            }






            bool hasTminStorage = false;

            res[thread_id].res_p.clear();

            //cout << res[thread_id].R2 << " ";

            for (int h=(p-1);  Syy > 0 && h >= start_h && ( res[thread_id].R2 > minR2store || storePInd ); h--)
            {
                //cout << " " << sig_res<< " " <<  var_xR[ii*r+h-l] << endl;


                if(var_x[ii*p+h] < 0.0)
                {
                    cout << "Error with var_x < 0 idx:" << ii*p+h << "\n";exit(1);
                }
                SE = sig_res*sqrt(var_x[ii*p+h]*varFactor);


                t=fabs(b[h]/SE);
                //cout << t << endl;


//                if(jj==(ii+A_blck_offset))
//                {
//                    cout << "b:" << b[h] << " R2:" << res[thread_id].R2 << " sig:" << sig_res << " vxr:" <<  var_xR[ii*r+h-l]   << " SE" <<  SE <<  " t:" << t;
//                    cout << endl;
//                }



                if( t > minTstore || storePInd)
                {

                    //cout << R2 << " " << Syy << " " << SST <<  "\t:\t";
                     //cout << (Syy/(n_corrected-1)) << " " << (SST/(n-1)) << endl;
//                     if(n>900)
//                        cout << jj<< ":" << ii << " " << b[h] << " " <<  R2 << " " <<  SReg <<  " " << SST*varFactor <<" " << 1.0-((Syy/(n_corrected-1))/(SST*varFactor/(n-1))) <<  " "<< Syy << " " << t << endl;
                    if(h >= l)
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

                    //cout << "found:"  << jj << "-"<< ii<< ":"<< A_blck_offset << "\n";
                    rdata.AR_name_idx = ii*r;
                    rdata.AL_name_idx = h;

                    if(h >= l)
                    {
                        rdata.AR_name_idx += h-l;
                        rdata.AL_name_idx = -1;
                    }


                    res[thread_id].res_p.push_back(rdata);

                    if(minPdisp > ptemp)
                    {
                        total_results++;
                    }

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


        delete []res;





}

/*!
 * \brief Simple sum of squares for blocks of vectors (Used for Y)
    Note that the method uses x as generic name, but it is meant for singel vector of uncorrelated variables.

 * \param block of Vectors, dimensions, destination block of results, real size of the data (list<long int>* indexs_Data)
 */
void Algorithm::sumSquares(float* Data, int cols, int rows, float* ssData, list<long int>* indexs_Data)
{
    //cout <<"\nvars" << endl;
    for (int h=0; h< cols; h++)
    {
        long double sx = 0;
        long double sxx = 0;
        float* x = &Data[h*rows];

        int n = rows;
        if(indexs_Data)
            n =n-indexs_Data[h].size();

        for (int k=0; k < rows; k++)
        {
                sx += x[k];
                sxx += x[k]*x[k];
        }
        ssData[h] = (float)(sxx - (sx*sx)/n);

        //cout << ssData[h] << " ";
        if(ssData[h] < 0.0)
        {
            cout << "SS failed SS:" << ssData[h] << " sxx:"<< sxx <<" sx2/n:" << sx*sx/n << endl;
        }

    }
    //cout << endl;

}


/*!
 * \brief Adds up missing values in a vector product form to later subtract them from precomputed matrix products.


 * \param Blocks of data, dimnsions, destination blocks of data, missing indexes.
 */
void Algorithm::generate_correction_missings_Str(int y_block_size, int xr_block_size, int n, int l, int r,
    float* __restrict XL_block, float* __restrict XR_block, float* __restrict Str_block, float* __restrict Sbr_block,
        list<long int>* index_nans)
{

    memset (Str_block,0,y_block_size*l*r*xr_block_size*sizeof(float));
    memset (Sbr_block,0,y_block_size*r*r*xr_block_size*sizeof(float));

    float* __restrict XL = XL_block;


    //int n_iters;



    for(int jj = 0; jj < y_block_size; jj++)
    {

        //n_iters = (index_nans[jj].size() + max_n -1)/max_n;



        #pragma omp parallel default(shared)
        {
        #pragma omp for nowait
        for(int b = 0; b < xr_block_size; b++)
        {

            float* __restrict Str = &Str_block[jj*xr_block_size*l*r+ b*l*r];
            float* __restrict Sbr = &Sbr_block[jj*xr_block_size*r*r+ b*r*r];
            float* __restrict XR = &XR_block[b*r*n];

            for(int i = 0; i < r; i++)
            {
                    for(int j = 0; j < l; j++)
                    {
                        for (list< long int  >::iterator it=index_nans[jj].begin(); it != index_nans[jj].end(); ++it)
                        {
                             int k = *it;
                            Str[i*l+j] += XL[j*n+k]*XR[i*n+k];
                            if(j<r)
                                Sbr[i*r+j] += XR[j*n+k]*XR[i*n+k];
                        }
                    }
            }


        }

        }
    }




}


/*!
 * \brief Adds up missing values in a vector product form to later subtract them from precomputed matrix products.


 * \param Blocks of data, dimnsions, destination blocks of data, missing indexes.
 */
void Algorithm::generate_correction_missings_Stl(int y_block_size, int x1_block_size, int x2_block_size, int n, int cols1, int cols2, float* X1_block, float* X2_block, float* S_block,
        list<long int>* index_nans)
{

    memset (S_block,0,y_block_size*cols1*cols2*x1_block_size*sizeof(float));
    //cout << "gen\n";

    for(int yb = 0; yb < y_block_size; yb++)
    {

        #pragma omp parallel default(shared)
        {
        #pragma omp for nowait
        for(int b = 0; b < x1_block_size; b++)
        {
            float* X1 = &X1_block[b*cols1*n];
            float* X2 = X1;

            float* S = &S_block[yb*x1_block_size*cols1*cols2+ b*cols1*cols2];


            for(int i = 0; i < cols1; i++)
            {
                for(int j = 0; j < cols2; j++)
                {
                    for (list< long int  >::iterator it=index_nans[yb].begin(); it != index_nans[yb].end(); ++it)
                    {
                        int k = *it;
                        S[i*cols1+j] += X1[i*n+k]*X2[j*n+k];
                    }
                    //cout << S[i*cols1+j] << " ";
                }
            }
            //cout << endl;
        }
        }

    }
}

/*!
 * \brief Will do the actual substraction of corrections on the SPD matrix S based on it parts.


 * \param parts of the corrected SPD containing the corrections, dimensions and estination (original S matrix).
 */
void Algorithm::correct_missings_S(int p,int l,int r, float* __restrict S, float*__restrict  Stl, float* __restrict Sbr, float* __restrict Str)
{
    for(int i = 0; i < l; i++)
    {
        for(int j = 0; j < l; j++)
        {
            S[i*p+j] -= Stl[i*l+j];
        }
    }
    //cout << "corr: ";
    for(int i = 0; i < r; i++)
    {
        for(int j = 0; j < r; j++)
        {
            S[(l*p+l)+i*p+j] -= Sbr[i*r+j];
            //cout << Sbr[i*r+j] <<" ";
        }
    }
    //cout << endl;

    for(int i = 0; i < r; i++)
    {
        for(int j = 0; j < l; j++)
        {
            S[l*p+j] -= Str[i*l+j];
        }
    }


}







