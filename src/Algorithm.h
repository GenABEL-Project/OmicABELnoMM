#ifndef ALGORITHM_H
#define ALGORITHM_H

#include "Definitions.h"
#include "Utility.h"
#include "AIOwrapper.h"


#define FULL_NEQ 1
#define FULL_NEQ_STR "FULL_NEQ"

#define P_NEQ 2
#define P_NEQ_STR "P_NEQ"

#define P_NEQ_B_OPT 3
#define P_NEQ_B_OPT_STR "P_NEQ_B_OPT"

#define FULL_QR 4
#define FULL_QR_STR "FULL_QR"

#define P_QR 5
#define P_QR_STR "P_QR"

#define P_QR_B_OPT 6
#define P_QR_B_OPT_STR "P_QR_B_OPT"

#define P_NEQ_B_OPT_MD 7
#define P_NEQ_B_OPT_MD_STR "P_NEQ_B_OPT_MD"



class Algorithm
{
    public:

        Algorithm();

        virtual ~Algorithm();

        void solve(struct Settings params, struct Outputs &out, int type);

        void partialNEQ_Blocked_STL_MD(struct Settings params, struct Outputs &out);

        void partialNEQ_Blocked_STL(struct Settings params, struct Outputs &out);

    protected:
    private:

        AIOwrapper AIOfile;

        void t_students_cdf(int y_amount,int a_amount,int p, type_precision* T, type_precision* P, int deg_freedom);

        void hpc_statistics(list<long int>* indexs_AL,list<long int>* indexs_AR, list<long int>* indexs_Y, int n,
                type_precision* A, int a_amount, type_precision* y, int jj, type_precision* B, int p, int l,int r, type_precision* T, type_precision* R2,type_precision* P);


        void check_result(type_precision* AL, type_precision* AR,int rowsA,int colsA, int rhs,int colsAR,
                                                    type_precision* y, type_precision* res);

        void update_R(type_precision* R, type_precision* topRr, type_precision* botRr,int dim1, int dim2, int r);

        type_precision* extract_R(type_precision* A,int dim1_A, int dim2_A);

        type_precision* prepare_R(type_precision* RL,int dim1_A, int dim2_AL,int dim2_AR);

        void prepare_QY(type_precision* qy, type_precision* top,type_precision* bot,int dim1_QY, int dim2_QY,int dim1_qy_bot,int bot_blocks );

        void prepare_Bfinal(type_precision* bfinal, type_precision* bsource, int a_amount, int y_amount, int p);

        void extract_subMatrix(type_precision* source, type_precision* dest,int dim1_source, int dim2_source,
                                                            int dim1_ini,int dim1_end,int dim2_ini,int dim2_end);
        void build_S(type_precision* S,type_precision* Stl,type_precision* Str,type_precision* Sbr,int l,int r);

};

#endif // ALGORITHM_H
