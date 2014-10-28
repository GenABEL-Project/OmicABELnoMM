#ifndef ALGORITHM_H
#define ALGORITHM_H

#include "Definitions.h"
#include "Utility.h"
#include "AIOwrapper.h"

#include <limits>
#include <iomanip>

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

    void applyDefaultParams(struct Settings &params);


 protected:


 private:
    list < resultH > sigResults;

    int max_threads;

    unsigned int total_results;

    float minPdisp;


    float minTstore;
    float minR2store;
    bool storePInd;
    bool disp_cov;


    float* SYY;


    void generate_correction_missings_Str(int y_block_size,
                                          int xr_block_size,
                                          int n, int l, int r,
                                          float* __restrict XL_block,
                                          float* __restrict XR_block,
                                          float* __restrict Str_block,
                                          float* __restrict Sbr_block,
                                          list<long int>* index_nans);

    void generate_correction_missings_Stl(int y_block_size,
                                          int x1_block_size,
                                          int x2_block_size,
                                          int n, int cols1, int cols2,
                                          float* X1_block,
                                          float* X2_block,
                                          float* S_block,
                                          list< long int>* index_nans);

    void correct_missings_S(int p, int l, int r, float* S, float* Stl,
                            float* Sbr, float* Str);


    void t_students_cdf(int y_amount, int a_amount, int p, float* T,
                        float* P, int deg_freedom);

    void sumSquares(float* Data, int cols, int rows, float* ssData,
                    list<long int>* indexs_Data);

    void hpc_statistics(int YxALmiss, int n, int realN, float* AL,
                        float* AR, int a_amount, float* y, int jj,
                        float varY, float* B, int p, int l, int r,
                        float* var_x, int y_blck_offset, int A_blck_offset,
                        list < resultH >* sigResults);


    void check_result(float* AL, float* AR, int rowsA, int colsA, int rhs,
                      int colsAR,float* y, float* res, struct Settings params,
                      int iX, int iiX, int jY, int jjY);





    void prepare_Bfinal(float* bfinal, float* bsource, int a_amount,
                        int y_amount, int p);


    void build_S(float* S, float* Stl, float* Str, float* Sbr, int l, int r);


    void hpc_SSY(int n, int p, int l, int r, float* __restrict  AL,
                 float* __restrict  AR, int a_amount, int y_block_size,
                 float* __restrict y, float* __restrict B);


    void singleSSY( float* __restrict sum, float* __restrict b,
                    float* __restrict b2, float* __restrict b3,
                    float* __restrict b4, float* __restrict y1,
                    float* __restrict y2, float* __restrict y3,
                    float* __restrict y4, float* __restrict sy,
                    float* __restrict Xl_n, float* __restrict XR_n,
                    int l, int r, int p, int n, int thead_id);


    inline float singleSSY(float* __restrict b, float* __restrict y,
                           float* __restrict sy, float* __restrict Xl_n,
                           float* __restrict XR_n, int l, int r, int p,
                           int n, int thead_id);
};

#endif // ALGORITHM_H
