/* The functions in this file are written by myself, aiming to 
 * solve the transport-equations */

/* $Id: other_functions.h, v 0.1 2015/09/04 10:11:00 zlb Exp $ */

#ifndef OTHER_FUNCTIONS_H
#define OTHER_FUNCTIONS_H
#endif

#include "phg.h"


/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/* beta_function.h
 * this function is to compute the value of \beta_{lm}.
 * the int l,int m are the subscript of \beta_{lm}, 
 * int c==1 <==> 1*(\beta_{lm}); c==-1 <==> (-1)*(\beta_{lm}) */
//FLOAT beta_function(INT l, INT m, INT c);


/* alpha_function.h
 * this function is to compute the value of \alpha_{lm} */
//FLOAT alpha_function(INT l, INT m, INT c);


/* gamma_function.h
 * this function is to compute the value of \gamma_{lm} */
//FLOAT gamma_function(INT l, INT m, INT c);


/* eta_function.h
 * this function is to compute the value of \eta_{lm} */
//FLOAT eta_function(INT l, INT m, INT c);


/* alpha_function.c 
 * this function is to compute the value of \alpha_{lm}*/
FLOAT
alpha_function(INT l, INT m, FLOAT c);


/* lambda_function.c 
 * this function is to compute the value of \lambda_{lm}*/
FLOAT
lambda_function(INT l, INT m, FLOAT c);
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/* build_D_x_matrix.h; build_D_y_matrix.h; build_D_z_matrix.h.
 * These functions is written for building the matrix D_x, D_y, D_z,
 * moreover the double *D_x is to store the matrix D_x which elements 
 * are arranged by row.
 * The int N stands for P_N approximate, that is l==N in the P_N 
 * approximate. */
void build_D_x_matrix(INT nY, INT PN, FLOAT **D_x);
void build_D_y_matrix(INT nY, INT PN, FLOAT **D_y);
void build_D_z_matrix(INT nY, INT PN, FLOAT **D_z);
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/*
 * Next considering the matrix A(or A^T) multiply matrix B
 * here, we just take the consideration that A and B are same 
 * order: A(n*n), B(n*n)
 * 
 * the BOOLEAN is defined in phg.h, the values are TRUE(1) and FLASE(0),
 * when tran=TRUE stands for C= A^T * B, tran=FALSE stands for C= A * B.
 * And more, here we just consider A and B are the same size, so the code
 * is the simplest.
 * if C= A * B, then C_{ij}=\sum_k A_{ik}*B_{kj}.
 * if C= A^T *B, then C_{ij}=\sum_k A_{ki}*B_{kj}.
 */

void
MatA_multiply_MatB(BOOLEAN tran, INT n, FLOAT **A, FLOAT **B, FLOAT **C);



/*
 * The following code is also to compute the matrixes A*B=C,
 * but this time the A and B may have different rows and colls.
 * Also:
 * tran=TRUE stands for C= A^T * B, tran=FALSE stands for C= A * B.
 */
void
MatA_x_MatB(BOOLEAN tran, INT rowA, INT colA, FLOAT **A, 
        INT rowB, INT colB, FLOAT **B, FLOAT **C);

/**********************************************************************************/
/*--------------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/*
 * Compute C=A*M*B, or C=A^T*M*B
 * when tran1=TRUE, tran2=FLASE stands for C= A^T * M * B, 
 * tran1=FALSE, tran2=FLASE stands for C= A * M * B.
 * And more, here we just consider A and B are the same size, 
 * so the code is simple.
 *
 * In the function name Compute_D_xx_x, 
 * xx may stands for xx, xy, xz, yx, yy, yz, zx, zy, zz.
 * x may stands for \delta, I.
 *
 */
/*
void
Compute_D_xx_x(FLOAT **A, BOOLEAN tran1, FLOAT **B, BOOLEAN tran2, 
        INT n, FLOAT **C);
*/
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
void
test_2_dim_pointer(int n, double **A);
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
void
arrangeMatrixInRows(int row, int col, double **A, double *C);
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/*
 * The following funcion is to compute the integration of the...
 *
 * 下面的函数是在单元e上计算第m，n个基函数的偏导数乘积的积分，
 * 其中ParGradi(ParGradj)取0时表示对x的偏导数，1表示对y的偏导数，
 * 2表示对z的偏导数。
 *
 * 此函数是在 quad.c 中 phgQuadGradBasAGradBas 函数的基础上改的，
 * 并且只用到了A=NULL 的情况，所以A！=NULL 的都删掉了。
 */
FLOAT
phgQuadBasParGradi_BasParGradj(ELEMENT *e, int ParGradi, DOF *u, int n, 
        int ParGradj, DOF *v, int m, int order);
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/*
 * The following funcion is to compute the integration of the...
 *
 * 下面的函数是在单元e上计算第n个基函数的偏导数乘积的积分，
 * 其中ParGradi取0时表示不求导，1表示对x的偏导数，
 * 2表示对y的偏导数，3表示对z的偏导数。
 *
 * 此函数是在 quad.c 中 phgQuadGradBasAGradBas 函数的基础上改的，
 * 并且只用到了A=NULL 的情况，所以A！=NULL 的都删掉了。
 */
FLOAT
phgQuadBasParGradi(ELEMENT *e, int ParGradi, DOF *u, int n, int order);
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/*
 * To assemble the matrixes: F_xx, F_xy ...
 */
void
assemble_Fxx_matrixes(MAT *F_xx, MAT *F_xy, MAT *F_xz, MAT *F_x0, MAT *F_0x, 
        MAT *F_yx, MAT *F_yy, MAT *F_yz, MAT *F_y0, MAT *F_0y,
        MAT *F_zx, MAT *F_zy, MAT *F_zz, MAT *F_z0, MAT *F_0z,
        MAT *F_00);
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/* 
 * 
 */ 
void
build_rhs(FLOAT **D_x, FLOAT **D_y, FLOAT **D_z, DOF *u_F, VEC *rhs, INT nY);
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
void
build_linear_system(MAT *F_xx, MAT *F_xy, MAT *F_xz, MAT *F_x0, MAT *F_0x, 
        MAT *F_yx, MAT *F_yy, MAT *F_yz, MAT *F_y0, MAT *F_0y,
        MAT *F_zx, MAT *F_zy, MAT *F_zz, MAT *F_z0, MAT *F_0z,
        MAT *F_00, MAT *XmFbd_00, MAT *XpFbd_00, MAT *YmFbd_00,
        MAT *YpFbd_00, MAT *ZmFbd_00, MAT *ZpFbd_00,
        FLOAT **D_x, FLOAT **D_y, FLOAT **D_z, VEC *rhs, INT nY);
/*
 * To build the linear system. 
 * this function is aimed to combine the build_rhs() and the assemble_Fxx_matrixes()
 * functions, more this function will add one function which will build boundary face
 * matrix.
 *
 * Input: All of the parameters are the input.
 *
 */
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/




/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/




/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/
