/* The functions in this file are written by myself, aiming to 
 * solve the transport-equations */

/* $Id: other_functions.c, v 0.1 2015/09/04 10:11:00 zlb Exp $ */
#include "phg.h"
#include "phg/quad-gauss.h"
#include "phg/quad-permu.h"
#include "other_functions.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

# define BasisOrder(u, e, i) (!DofIsHP(u) ? (u)->type->order : \
	(u)->hp->info->types[(u)->hp->elem_order[e->index]]->order)

/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* beta_function.c */
/*
FLOAT
beta_function(INT l, INT m, INT c){
    FLOAT r;
    r=c*sqrt(((FLOAT)l-m)*(l-m-1)/(4*(2*l-1)*(2*l+1)));
    return r;
}//endof-beta_function
*/


/* alpha_function.c */
/*
FLOAT
alpha_function(INT l, INT m, INT c){
    FLOAT r;
    r=c*sqrt(((FLOAT)l+m+2)*(l+m+1)/(4*(2*l+1)*(2*l+3)));
    return r;
}//endof-alpha_function
*/


/* gamma_function.c */
/*
FLOAT
gamma_function(INT l, INT m, INT c){
    FLOAT r;
    r=c*sqrt(((FLOAT)l-m)*(l+m)/((2*l-1)*(2*l+1)));
    return r;
}//endof-gamma_function
*/


/* eta_function.c */
/*
FLOAT
eta_function(INT l, INT m, INT c){
    FLOAT r;
    r=c*sqrt(((FLOAT)l+m+1)*(l-m+1)/((2*l+3)*(2*l+1)));
    return r;
}//endof-eta_function
*/

/* alpha_function.c */
FLOAT
alpha_function(INT l, INT m, FLOAT c){
    FLOAT r;
    FLOAT numerator;//分子
    FLOAT denominator;//分母
    numerator=(FLOAT)(l+m+1)*(l+m+2);
    denominator=(FLOAT)4*(2*l+1)*(2*l+3);
    r=c*(-1.0)*sqrt(numerator/denominator);
    return r;
}//endof-alpha_function


/* lambda_function.c */
FLOAT
lambda_function(INT l, INT m, FLOAT c){
    FLOAT r;
    FLOAT numerator;//分子
    FLOAT denominator;//分母
    numerator=(FLOAT)(l+m+1)*(l-m+1);
    denominator=(FLOAT)(2*l+1)*(2*l+3);
    r=c*sqrt(numerator/denominator);
    return r;
}//endof-lambda_function
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/* build_D_x_matrix.c */
void
build_D_x_matrix(INT nY, INT PN, FLOAT **D_x)
{
    /* it is easy to compute that in P_N(0,1,...N) approximate the number 
     * of how many Y_l^m is (N+1)*(N+1), so the D_x is a n*n matrix */
    ///*
    //INT n;
    //n=(N+1)*(N+1);
    INT i,j;
    for(i=0;i<nY;i++){
        for(j=0;j<nY;j++)
            *(*(D_x+i)+j)=0.0;
    }
    // Just only to give an initial value to matrix D_x 
    //*/

    INT l,m;
    INT row,col;//ptr;//ptr is to transform row and col to 1 dim in row storage
    for(l=0;l<=PN;l++){
        for(m=-l;m<=l;m++){
            col=l*l+l+m;
            if(m>=1){
                if(l-1>=0){
                    row=(l-1)*(l-1)+(l-1)+(m-1);
                    *(*(D_x+row)+col)=alpha_function(l-1,m-1,1.0);
                    if(m+1<=l-1){
                        row=(l-1)*(l-1)+(l-1)+(m+1);
                        *(*(D_x+row)+col)=alpha_function(l-1,-m+1,-1.0);
                    }
                }//endof-if(l-1>=0)

                if(l+1<=PN){
                    row=(l+1)*(l+1)+(l+1)+(m-1);
                    *(*(D_x+row)+col)=alpha_function(l,-m,-1.0);
                    row=(l+1)*(l+1)+(l+1)+(m+1);
                    *(*(D_x+row)+col)=alpha_function(l,m,1.0);
                }
            }//endof-if(m>=1)

            if(m==0){
                if(l-1<=0 && m+1<=l-1){
                    row=(l-1)*(l-1)+(l-1)+(m+1);
                    *(*(D_x+row)+col)=alpha_function(l-1,-m-1,-sqrt(2.0));
                }
                if(l+1<=PN){
                    row=(l+1)*(l+1)+(l+1)+(m+1);
                    *(*(D_x+row)+col)=alpha_function(l,m,sqrt(2.0));
                }
            }//endof-if(m==0)

            if(m==-1){
                if(l-1>=0 && m-1>=-(l-1)){
                    row=(l-1)*(l-1)+(l-1)+(m-1);
                    *(*(D_x+row)+col)=alpha_function(l-1,m-1,-1.0);
                }
                if(l+1<=PN){
                    row=(l+1)*(l+1)+(l+1)+(m-1);
                    *(*(D_x+row)+col)=alpha_function(l,-m,1.0);
                }
            }//endof-if(m==-1)

            if(m<=-2){
                if(l-1>=0){
                    row=(l-1)*(l-1)+(l-1)+(m+1);
                    *(*(D_x+row)+col)=alpha_function(l-1,-m-1,1.0);
                    if(m-1>=-(l-1)){
                        row=(l-1)*(l-1)+(l-1)+(m-1);
                        *(*(D_x+row)+col)=alpha_function(l-1,m-1,1.0);
                    }
                }

                if(l+1<=PN){
                    row=(l+1)*(l+1)+(l+1)+(m+1);
                    *(*(D_x+row)+col)=alpha_function(l,m,-1.0);
                    
                    row=(l+1)*(l+1)+(l+1)+(m-1);
                    *(*(D_x+row)+col)=alpha_function(l,-m,1.0);
                }
            }//endof-if(m<=-2)
            
        }//endof-for(m=-l;...)
    }//endof-for(l=0;...)
}//endof-build_D_matrix


/* build_D_y_matrix.c */
void
build_D_y_matrix(INT nY, INT PN, FLOAT **D_y)
{
    /* it is easy to compute that in P_N approximate the number 
     * of how many Y_l^m is (N+1)*(N+1), so the D_x is a n*n matrix */
    ///*
    //INT n;
    //n=(N+1)*(N+1);
    INT i,j;
    for(i=0;i<nY;i++){
        for(j=0;j<nY;j++)
            *(*(D_y+i)+j)=0.0;
    }
    // Just only to give an initial value to matrix D_y 
    //*/
    
    INT l,m;
    INT row,col;//ptr;//ptr is to transform row and col to 1 dim in row storage
    for(l=0;l<=PN;l++){
        for(m=-l;m<=l;m++){
           col=l*l+l+m;
           if(m>=1){
               if(l-1>=0){
                   row=(l-1)*(l-1)+(l-1)-(m-1);
                   *(*(D_y+row)+col)=alpha_function(l-1,m-1,-1.0);
                   if(m+1<=l-1){
                       row=(l-1)*(l-1)+(l-1)-(m+1);
                       *(*(D_y+row)+col)=alpha_function(l-1,-m-1,1.0);
                   }
               }
               if(l+1<=PN){
                   row=(l+1)*(l+1)+(l+1)-(m-1);
                   *(*(D_y+row)+col)=alpha_function(l,-m,-1.0);
                   row=(l+1)*(l+1)+(l+1)-(m+1);
                   *(*(D_y+row)+col)=alpha_function(l,m,1.0);
               }
           }//endof-if(m>=1)

           if(m==1){
               if(l-1>=0 && m+1<=l-1){
                   row=(l-1)*(l-1)+(l-1)-(m-1);
                   *(*(D_y+row)+col)=alpha_function(l-1,-m-1,-1.0);
               }
                if(l+1<=PN){
                    row=(l+1)*(l+1)+(l+1)-(m+1);
                    *(*(D_y+row)+col)=alpha_function(l,m,1.0);
                }
           }//endof-if(m==1)

           //m==0时，系数都为0，所以不用处理

           if(m<=-1){
               if(l-1>=0){
                   row=(l-1)*(l-1)+(l-1)-(m+1);
                   *(*(D_y+row)+col)=alpha_function(l-1,-m-1,1.0);
                   if(-(m-1)<=l-1){
                       row=(l-1)*(l-1)+(l-1)-(m-1);
                       *(*(D_y+row)+col)=alpha_function(l-1,m-1,1.0);
                   }
               }
               if(l+1<=PN){
                   row=(l+1)*(l+1)+(l+1)-(m+1);
                   *(*(D_y+row)+col)=alpha_function(l,m,-1.0);
                   row=(l+1)*(l+1)+(l+1)-(m-1);
                   *(*(D_y+row)+col)=alpha_function(l,-m,-1.0);
               }
           }//endof-if(m<=-1)
        }//endof-for(m=-l;...)
    }//endof-for(l=0;...)
}//endof-build_D_matrix


/* build_D_z_matrix.c */
void
build_D_z_matrix(INT nY, INT PN, FLOAT **D_z){
    ///*
    //INT n=(N+1)*(N+1);
    INT i,j;
    for(i=0;i<nY;i++){
        for(j=0;j<nY;j++)
	    //printf("i=%d, j=%d \n",i,j);
            *(*(D_z+i)+j)=0.0;
    }
    //printf("testing build_D_z_matrix\n");
    //*/
    INT l,m;
    INT col,row;//ptr;
    for(l=0;l<=PN;l++){
        for(m=-l;m<=l;m++){
            col=l*l+l+m;
            if(l-1>=0 && m!=-l && m!=l){
                row=(l-1)*(l-1)+(l-1)+m;
                *(*(D_z+row)+col)=lambda_function(l-1,m,1.0);
            }
            if(l+1<=PN){
                row=(l+1)*(l+1)+(l+1)+m;
                *(*(D_z+row)+col)=lambda_function(l,m,1.0);
                //printf("right2: l=%d,m=%d, row=%d, col=%d, D_Z[%d][%d]=%f \n",l,m,row,col,row,col,*(*(D_z+row)+col));
            }
        }//endof-for(m=-l)
    }//endof-for(l=0)
}//endof-build_D_z_matrix
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*
 * Next considering the matrix A(or A^T) multiply matrix B
 * here, we just take the consideration that A and B are same 
 * order: A(n*n) B(n*n)
 * */
void
MatA_multiply_MatB(BOOLEAN tran, INT n, FLOAT **A, FLOAT **B, FLOAT **C){
    /*
     * the BOOLEAN is defined in phg.h, the values are TRUE(1) and FLASE(0),
     * when tran=TRUE stands for C= A^T * B, tran=FALSE stands for C= A * B.
     * And more, here we just consider A and B are the same size, so the code
     * is the simplest.
     * if C= A * B, then C_{ij}=\sum_k A_{ik}*B_{kj}.
     * if C= A^T *B, then C_{ij}=\sum_k A_{ki}*B_{kj}.
     * */
    INT i,j,k;
    FLOAT val=0.0;
    if(tran==TRUE){
        for(i=0;i<n;i++){
            for(j=0;j<n;j++){
                val=0.0;
                for(k=0;k<n;k++)
                    val+=A[k][i]*B[k][j];
                *(*(C+i)+j)=val;
            }//endof-for(j=0;...)
        }//endof-for(i=0;...)
    }//endof-if(tran==TRUE)
    
    if(tran==FALSE){
        for(i=0;i<n;i++){
            for(j=0;j<n;j++){
                val=0.0;
                for(k=0;k<n;k++)
                    val+=A[i][k]*B[k][j];
                *(*(C+i)+j)=val;
            }//endof-for(j=0;...)
        }//endof-for(i=0;...)
    }//endof-if(tran==FALSE)

}//endof-MatA_multiply_MatB



/*
 * The following code is also to compute the matrixes A*B=C,
 * but this time the A and B may have different rows and colls.
 * Also:
 * tran=TRUE stands for C= A^T * B, tran=FALSE stands for C= A * B.
 */
void
MatA_x_MatB(BOOLEAN tran, INT rowA, INT colA, FLOAT **A, 
        INT rowB, INT colB, FLOAT **B, FLOAT **C){
    
    INT i,j,k;
    FLOAT r;

    if(tran==FALSE){
        if(colA!=rowB){
            printf("In function 'MatA_x_MatB()', tran==FLASE,but colA!=rowB, abort!\n");
            abort();
        }//endof_if(colA!=rowB)

        for(i=0;i<rowA;++i){
            for(j=0;j<colB;++j){
                r=0.0;
                for(k=0;k<rowB;++k)
                    r+=A[i][k]*B[k][j];
                *(*(C+i)+j)=r;
            }
        }//endof_for(i=0;i<rowA;++i)
    }//endof_if(tran==FLASE)

    if(tran==TRUE){
        if(rowA!=rowB){
            printf("In function 'MatA_x_MatB()', tran==TRUE, but rowA!=rowB, abort!");
            abort();
        }//endof_if(colA!=rowB)

        for(i=0;i<colA;++i){
            for(j=0;j<colB;++j){
                r=0.0;
                for(k=0;k<rowB;++k)
                    r+=A[k][i]*B[k][j];
                *(*(C+i)+j)=r;
            }
        }//endof_for(i=0;i<rowA;++i)
    }//endof_if(tran==TRUE)

}//endof
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*
 * This code is to compute C=A*M*B, or C=A^T * M * B
 * first using the function MatA_multiply_MatB(), compute tmp = A^T*M, 
 * then also using function MatA_multiply_MatB(), compute C = tmp*B,
 *
 */
/*
void
Compute_D_xx_x(FLOAT **A, BOOLEAN tran1, FLOAT **B, BOOLEAN tran2, 
        INT n, FLOAT **C)
{
    INT i;
    FLOAT **tmp;
    tmp=(FLOAT **)malloc(sizeof(FLOAT *)*n);
    for(i=0;i<n;i++){
        *(tmp+i)=(FLOAT *)malloc(n*sizeof(FLOAT));
    }
    
    MatA_multiply_MatB(tran1, n, A, M, tmp);
    MatA_multiply_MatB(tran2, n, tmp, B, C);

}//endof_Compute_D_xx_x
*/
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
void
test_2_dim_pointer(int n, double **A){
    int i;
    int j;
    for(i=0;i<n;i++)
        *(A+i)=(double*)malloc(n*sizeof(double));
    
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            A[i][j]=i+j*1.0;

    printf("entering test_2_dim_pointer\n");
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            printf("i=%d, j=%d",i,j);
            printf("A[%d][%d]=%f  ",i,j,A[i][j]);
        }
    }


}
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
void
arrangeMatrixInRows(int row, int col, double **A, double *C)
{
    int i, j;
    for(i=0; i<row; ++i){
        for(j=0; j<col; ++j){
            *(C+i*row+j)=*(*(A+i)+j);
        }
    }
}
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*
 * The following funcion is to compute the integration of the...
 *
 * 下面的函数是在单元e上计算第m，n个基函数的偏导数乘积的积分，
 * 其中ParGradi(ParGradj)取0时表示不求导，1表示对x求偏导数，
 * 2表示对y的偏导数，3表示对z的偏导数。
 *
 * 此函数是在 quad.c 中 phgQuadGradBasAGradBas 函数的基础上改的，
 * 并且只用到了A=NULL 的情况，所以A！=NULL 的都删掉了。
 */
FLOAT
phgQuadBasParGradi_BasParGradj(ELEMENT *e, int ParGradi, DOF *u, int n, 
        int ParGradj, DOF *v, int m, int order)
{
    int i, j, nvalues = DofTypeDim(u);//对于Lagrange元，nvalues=1.
    int Pin1, Pin2, Pout1, Pout2;
    FLOAT d, d0;
    const FLOAT *g1, *g2, *gNull1, *gNull2, *gPGrad1, *gPGrad2, *w;
    QUAD *quad;

    assert(!SpecialDofType(u->type) && !SpecialDofType(v->type));
    assert(nvalues == DofTypeDim(v));

    if (order < 0) {
	    order = BasisOrder(u, e, n) - 1 + BasisOrder(v, e, m) - 1;
	    if (order < 0)
	        order = 0;
    }
    quad = phgQuadGetQuad3D(order);

    d = 0.;
    gNull1 = phgQuadGetBasisValues(e, u, n, quad);
    gNull2 = phgQuadGetBasisValues(e, v, m, quad);
    gPGrad1 = phgQuadGetBasisGradient(e, u, n, quad);
    gPGrad2 = phgQuadGetBasisGradient(e, v, m, quad);
    w = quad->weights;

    /* 下面的代码中"Pin, Pout" 的赋值完全是为了后面的for循环的代码统一 */
    if(ParGradi==0){
        g1=gNull1;
        Pin1=0;//这里因为不求导，所以g1就是基函数本身的高斯积分点的值
        Pout1=1;
    }
    else{
        g1=gPGrad1;
        Pin1=ParGradi-1;
        //这里就是根据gPGrad1的数据结构来的，具体的数据结构可以参看笔记
        Pout1=Dim;
    }

    if(ParGradj==0){
        g2=gNull2;
        Pin2=0;
        Pout2=1;
    }
    else{
        g2=gPGrad2;
        Pin2=ParGradj-1;
        Pout2=Dim;
    }
    
	for (i = 0; i < quad->npoints; i++) {
	    d0 = 0.;
	    for (j = 0; j < nvalues; j++) {
		    d0 += *(g1+Pin1) * *(g2+Pin2);
		    g1 += Pout1;
		    g2 += Pout2;
	    }
	    d += d0 * (*(w++));
	}
	return d * phgGeomGetVolume(u->g, e);
    
}
/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------*/
/*
 * The following funcion is to compute the integration of the...
 *
 * 下面的函数是在单元e上计算第n个基函数的偏导数的积分，
 * 其中ParGradi取0时表示不求导，1表示对x求偏导数，
 * 2表示对y的偏导数，3表示对z的偏导数。
 *
 * 此函数是在 quad.c 中 phgQuadGradBasAGradBas 函数的基础上改的，
 * 并且只用到了A=NULL 的情况，所以A！=NULL 的都删掉了。
 */
FLOAT
phgQuadBasParGradi(ELEMENT *e, int ParGradi, DOF *u, int n, int order)
{
    int i, j, nvalues = DofTypeDim(u);//对于Lagrange元，nvalues=1.
    int Pin1, Pout1;
    FLOAT d, d0;
    const FLOAT *g1, *gNull1, *gPGrad1, *w;
    QUAD *quad;

    if(nvalues != 1){
        printf("\n In function 'phgQuadBasParGradi(...)' DofTypeDim(u)!=1, quit! \n");
    }

    assert(nvalues == 1);

    if (order < 0) {
	    order = BasisOrder(u, e, n)-1;
	    if (order < 0)
	        order = 0;
    }
    quad = phgQuadGetQuad3D(order);
    //printf("in phgQuadBasParGradi() 'order'=%d, num of GaussPoints:%d \n",order,quad->npoints);//myDebug

    d = 0.;
    gNull1 = phgQuadGetBasisValues(e, u, n, quad);
    gPGrad1 = phgQuadGetBasisGradient(e, u, n, quad);
    w = quad->weights;

    /* 下面的代码中"Pin, Pout" 的赋值完全是为了后面的for循环的代码统一 */
    if(ParGradi==0){
        g1=gNull1;
        Pin1=0;//这里因为不求导，所以g1就是基函数本身的高斯积分点的值
        Pout1=1;
    }
    else{
        g1=gPGrad1;
        Pin1=ParGradi-1;
        //这里就是根据gPGrad1的数据结构来的，具体的数据结构可以参看笔记
        Pout1=Dim;
    }
    
	for (i = 0; i < quad->npoints; i++) {
	    d0 = 0.;
	    for (j = 0; j < nvalues; j++) {
		    d0 += *(g1+Pin1);
		    g1 += Pout1;
	    }
	    d += d0 * (*(w++));
	}
	return d * phgGeomGetVolume(u->g, e); 
}//endof-phgQuadBasParGradi(...)
/*--------------------------------------------------------------------------------*/
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
        MAT *F_00){
    DOF *u_h = F_xx->rmap->dofs[0];
    MAP *mapF = F_xx->rmap;
    int N = u_h->type->nbas;	/* number of basis functions in an element */
    GRID *g = u_h->g;
    ELEMENT *e;
    int i, j;
    FLOAT matF_xx[N][N], matF_xy[N][N], matF_xz[N][N], matF_x0[N][N], matF_0x[N][N];
    FLOAT matF_yx[N][N], matF_yy[N][N], matF_yz[N][N], matF_y0[N][N], matF_0y[N][N];
    FLOAT matF_zx[N][N], matF_zy[N][N], matF_zz[N][N], matF_z0[N][N], matF_0z[N][N];
    FLOAT matF_00[N][N];
    INT I[N];

    assert(u_h->dim == 1);

    ForAllElements(g, e) {
        for (i = 0; i < N; ++i) {
            I[i] = phgMapE2L(mapF, 0, e, i);
            for(j = 0; j <= i; ++j){
                matF_xx[i][j]=matF_xx[j][i]=phgQuadBasParGradi_BasParGradj(e, 1, u_h, j, 1, u_h, i, 4);
                matF_xy[i][j]=matF_xy[j][i]=phgQuadBasParGradi_BasParGradj(e, 1, u_h, j, 2, u_h, i, 4);
                matF_xz[i][j]=matF_xz[j][i]=phgQuadBasParGradi_BasParGradj(e, 1, u_h, j, 3, u_h, i, 4);
                matF_x0[i][j]=matF_x0[j][i]=phgQuadBasParGradi_BasParGradj(e, 1, u_h, j, 0, u_h, i, 4);
                matF_0x[i][j]=matF_0x[j][i]=phgQuadBasParGradi_BasParGradj(e, 0, u_h, j, 1, u_h, i, 4);

                matF_yx[i][j]=matF_yx[j][i]=phgQuadBasParGradi_BasParGradj(e, 2, u_h, j, 1, u_h, i, 4);
                matF_yy[i][j]=matF_yy[j][i]=phgQuadBasParGradi_BasParGradj(e, 2, u_h, j, 2, u_h, i, 4);
                matF_yz[i][j]=matF_yz[j][i]=phgQuadBasParGradi_BasParGradj(e, 2, u_h, j, 3, u_h, i, 4);
                matF_y0[i][j]=matF_y0[j][i]=phgQuadBasParGradi_BasParGradj(e, 2, u_h, j, 0, u_h, i, 4);
                matF_0y[i][j]=matF_0y[j][i]=phgQuadBasParGradi_BasParGradj(e, 0, u_h, j, 2, u_h, i, 4);

                matF_zx[i][j]=matF_zx[j][i]=phgQuadBasParGradi_BasParGradj(e, 3, u_h, j, 1, u_h, i, 4);
                matF_zy[i][j]=matF_zy[j][i]=phgQuadBasParGradi_BasParGradj(e, 3, u_h, j, 2, u_h, i, 4);
                matF_zz[i][j]=matF_zz[j][i]=phgQuadBasParGradi_BasParGradj(e, 3, u_h, j, 3, u_h, i, 4);
                matF_z0[i][j]=matF_z0[j][i]=phgQuadBasParGradi_BasParGradj(e, 3, u_h, j, 0, u_h, i, 4);
                matF_0z[i][j]=matF_0z[j][i]=phgQuadBasParGradi_BasParGradj(e, 0, u_h, j, 3, u_h, i, 4);

                matF_00[i][j]=matF_00[j][i]=phgQuadBasParGradi_BasParGradj(e, 0, u_h, j, 0, u_h, i, 4);
            }//endof_for(j = 0; j <= i; ++j)
        }//endof_for(i = 0; i <= N; ++i)

        for(i = 0;i < N;++i){/* add entries to MAT *F_xx, *F_xy ... */
            phgMatAddEntries(F_xx, 1, I+i, N, I, &matF_xx[i][0]);
            phgMatAddEntries(F_xy, 1, I+i, N, I, &matF_xy[i][0]);
            phgMatAddEntries(F_xz, 1, I+i, N, I, &matF_xz[i][0]);
            phgMatAddEntries(F_x0, 1, I+i, N, I, &matF_x0[i][0]);
            phgMatAddEntries(F_0x, 1, I+i, N, I, &matF_0x[i][0]);

            phgMatAddEntries(F_yx, 1, I+i, N, I, &matF_yx[i][0]);
            phgMatAddEntries(F_yy, 1, I+i, N, I, &matF_yy[i][0]);
            phgMatAddEntries(F_yz, 1, I+i, N, I, &matF_yz[i][0]);
            phgMatAddEntries(F_y0, 1, I+i, N, I, &matF_y0[i][0]);
            phgMatAddEntries(F_0y, 1, I+i, N, I, &matF_0y[i][0]);

            phgMatAddEntries(F_zx, 1, I+i, N, I, &matF_zx[i][0]);
            phgMatAddEntries(F_zy, 1, I+i, N, I, &matF_zy[i][0]);
            phgMatAddEntries(F_zz, 1, I+i, N, I, &matF_zz[i][0]);
            phgMatAddEntries(F_z0, 1, I+i, N, I, &matF_z0[i][0]);
            phgMatAddEntries(F_0z, 1, I+i, N, I, &matF_0z[i][0]);

            phgMatAddEntries(F_00, 1, I+i, N, I, &matF_00[i][0]);
        }//endof_for(i = 0;i < N;++i){/* add entries to MAT *F_xx, *F_xy ... */
    }//endof_ForAllElements(g, e)
}//endof_assemble_Fxx_matrixes()
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/* 
 * 
 */ 
void
build_rhs(FLOAT **D_x, FLOAT **D_y, FLOAT **D_z, DOF *u_F, VEC *rhs, INT nY)
{
    extern FLOAT sigma_t;
    extern FLOAT sigma_s;
    extern FLOAT q_0;
    GRID *g = u_F->g;
    ELEMENT *e;
    MAP *mapu_F=phgMapCreate(u_F,NULL);
    FLOAT coefRHS_x[nY], coefRHS_y[nY], coefRHS_z[nY], coefRHS_0[nY];
    int N = u_F->type->nbas;
    int num = DofGetDataCountGlobal(u_F);
    int i, j, k;

    VEC *rhs_xf, *rhs_yf, *rhs_zf, *rhs_0f;
    FLOAT vec_xf[N], vec_yf[N], vec_zf[N], vec_0f[N];// vec_0f1[N];//myDebug
    INT I[N];

    //DOF *f_c=phgDofNew(g, DOF_CONSTANT, 1, "f_c", DofNoAction);//myDebug
    //phgDofSetDataByValuesV(f_c, (FLOAT)1.0);//myDebug

    for(i=0;i<nY;++i){
        coefRHS_x[i]=q_0*D_x[0][i];
        coefRHS_y[i]=q_0*D_y[0][i];
        coefRHS_z[i]=q_0*D_z[0][i];
        coefRHS_0[i]=0;
    }
    coefRHS_0[0]=q_0*(sigma_t - sigma_s);

    rhs_xf=phgMapCreateVec(mapu_F,1);
    phgVecDisassemble(rhs_xf);
    rhs_yf=phgMapCreateVec(mapu_F,1);
    phgVecDisassemble(rhs_yf);
    rhs_zf=phgMapCreateVec(mapu_F,1);
    phgVecDisassemble(rhs_zf);
    rhs_0f=phgMapCreateVec(mapu_F,1);
    phgVecDisassemble(rhs_0f);

    assert(u_F->dim == 1);
    
    ForAllElements(g, e){ 
        for(i=0;i<N;++i){
            I[i] = phgMapE2L(mapu_F, 0, e, i);
            vec_0f[i]=phgQuadBasParGradi(e, 0, u_F, i, 4);
            //phgQuadDofTimesBas(e,f_c,u_F,i,4,vec_0f1+i);//myDebug
            vec_xf[i]=phgQuadBasParGradi(e, 1, u_F, i, 4);
            vec_yf[i]=phgQuadBasParGradi(e, 2, u_F, i, 4);
            vec_zf[i]=phgQuadBasParGradi(e, 3, u_F, i, 4);
        }//endof-for(i=0;i<N;++i)

	/* myDebug
        printf("test in build_rhs: test1\n");

        for(i=0;i<N;++i){
            printf("vec_0f[%d]=%f  ",i,*(vec_0f+i));
            printf("vec_0f1[%d]=%f \n",i,*(vec_0f1+i));
        }
	*/

        phgVecAddEntries(rhs_0f, 0, N, I, vec_0f);
        phgVecAddEntries(rhs_xf, 0, N, I, vec_xf);
        phgVecAddEntries(rhs_yf, 0, N, I, vec_yf);
        phgVecAddEntries(rhs_zf, 0, N, I, vec_zf); 
    }//endof-ForAllElements(g,e)

    k=0;
    for(i=0; i<nY; ++i){
        for(j=0; j<num; ++j){
            *(rhs->data + k) = *(rhs_xf->data + j) * coefRHS_x[i] + *(rhs_yf->data + j) * coefRHS_y[i] +
                *(rhs_zf->data + j) * coefRHS_z[i] + *(rhs_0f->data + j) * coefRHS_0[i];
            ++k;
        }
    }

    phgMapDestroy(&mapu_F);

    phgVecDestroy(&rhs_xf);
    phgVecDestroy(&rhs_yf);
    phgVecDestroy(&rhs_zf);
    phgVecDestroy(&rhs_0f);


}//endof-build_rhs
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/*
 * To build the linear system. 
 * this function is aimed to combine the build_rhs() and the assemble_Fxx_matrixes()
 * functions, more this function will add one function which will build boundary face
 * matrix. 
 */
void
build_linear_system(MAT *F_xx, MAT *F_xy, MAT *F_xz, MAT *F_x0, MAT *F_0x, 
        MAT *F_yx, MAT *F_yy, MAT *F_yz, MAT *F_y0, MAT *F_0y,
        MAT *F_zx, MAT *F_zy, MAT *F_zz, MAT *F_z0, MAT *F_0z,
        MAT *F_00,
        FLOAT **D_x, FLOAT **D_y, FLOAT **D_z, VEC *rhs, INT nY)
{
    DOF *u_h = F_xx->rmap->dofs[0];
    MAP *mapF = F_xx->rmap;
    int N = u_h->type->nbas;	/* number of basis functions in an element */
    GRID *g = u_h->g;
    ELEMENT *e;
    //int i, j;
    FLOAT matF_xx[N][N], matF_xy[N][N], matF_xz[N][N], matF_x0[N][N], matF_0x[N][N];
    FLOAT matF_yx[N][N], matF_yy[N][N], matF_yz[N][N], matF_y0[N][N], matF_0y[N][N];
    FLOAT matF_zx[N][N], matF_zy[N][N], matF_zz[N][N], matF_z0[N][N], matF_0z[N][N];
    FLOAT matF_00[N][N];
    FLOAT matF_bd_00[N][N],
    INT I[N];


    extern FLOAT sigma_t;
    extern FLOAT sigma_s;
    extern FLOAT q_0;
    //MAP *mapu_F=phgMapCreate(u_F,NULL); //下面即用上面的 mapF 代替
    FLOAT coefRHS_x[nY], coefRHS_y[nY], coefRHS_z[nY], coefRHS_0[nY];
    int num = DofGetDataCountGlobal(u_h);
    int i, j, k;

    VEC *rhs_xf, *rhs_yf, *rhs_zf, *rhs_0f;
    FLOAT vec_xf[N], vec_yf[N], vec_zf[N], vec_0f[N];
    INT I[N];

    for(i=0;i<nY;++i){
        coefRHS_x[i]=q_0*D_x[0][i];
        coefRHS_y[i]=q_0*D_y[0][i];
        coefRHS_z[i]=q_0*D_z[0][i];
        coefRHS_0[i]=0;
    }
    coefRHS_0[0]=q_0*(sigma_t - sigma_s);

    rhs_xf=phgMapCreateVec(mapF,1);
    phgVecDisassemble(rhs_xf);
    rhs_yf=phgMapCreateVec(mapF,1);
    phgVecDisassemble(rhs_yf);
    rhs_zf=phgMapCreateVec(mapF,1);
    phgVecDisassemble(rhs_zf);
    rhs_0f=phgMapCreateVec(mapF,1);
    phgVecDisassemble(rhs_0f);

    assert(u_h->dim == 1);

    ForAllElements(g, e) {
        for (i = 0; i < N; ++i) {
            I[i] = phgMapE2L(mapF, 0, e, i);
	        /* build matF_xx ... */
            for(j = 0; j <= i; ++j){
                matF_xx[i][j]=matF_xx[j][i]=phgQuadBasParGradi_BasParGradj(e, 1, u_h, j, 1, u_h, i, 4);
                matF_xy[i][j]=matF_xy[j][i]=phgQuadBasParGradi_BasParGradj(e, 1, u_h, j, 2, u_h, i, 4);
                matF_xz[i][j]=matF_xz[j][i]=phgQuadBasParGradi_BasParGradj(e, 1, u_h, j, 3, u_h, i, 4);
                matF_x0[i][j]=matF_x0[j][i]=phgQuadBasParGradi_BasParGradj(e, 1, u_h, j, 0, u_h, i, 4);
                matF_0x[i][j]=matF_0x[j][i]=phgQuadBasParGradi_BasParGradj(e, 0, u_h, j, 1, u_h, i, 4);

                matF_yx[i][j]=matF_yx[j][i]=phgQuadBasParGradi_BasParGradj(e, 2, u_h, j, 1, u_h, i, 4);
                matF_yy[i][j]=matF_yy[j][i]=phgQuadBasParGradi_BasParGradj(e, 2, u_h, j, 2, u_h, i, 4);
                matF_yz[i][j]=matF_yz[j][i]=phgQuadBasParGradi_BasParGradj(e, 2, u_h, j, 3, u_h, i, 4);
                matF_y0[i][j]=matF_y0[j][i]=phgQuadBasParGradi_BasParGradj(e, 2, u_h, j, 0, u_h, i, 4);
                matF_0y[i][j]=matF_0y[j][i]=phgQuadBasParGradi_BasParGradj(e, 0, u_h, j, 2, u_h, i, 4);

                matF_zx[i][j]=matF_zx[j][i]=phgQuadBasParGradi_BasParGradj(e, 3, u_h, j, 1, u_h, i, 4);
                matF_zy[i][j]=matF_zy[j][i]=phgQuadBasParGradi_BasParGradj(e, 3, u_h, j, 2, u_h, i, 4);
                matF_zz[i][j]=matF_zz[j][i]=phgQuadBasParGradi_BasParGradj(e, 3, u_h, j, 3, u_h, i, 4);
                matF_z0[i][j]=matF_z0[j][i]=phgQuadBasParGradi_BasParGradj(e, 3, u_h, j, 0, u_h, i, 4);
                matF_0z[i][j]=matF_0z[j][i]=phgQuadBasParGradi_BasParGradj(e, 0, u_h, j, 3, u_h, i, 4);

                matF_00[i][j]=matF_00[j][i]=phgQuadBasParGradi_BasParGradj(e, 0, u_h, j, 0, u_h, i, 4);


                /* build the boundary matrixes, there are 4 faces in one element */
                for(k=1;k<=4;++k){
                    if(e->bound_type[i]==BDRY_USER0){
                        
                    }
                }

            }//endof_for(j = 0; j <= i; ++j)

	        /* build rhs vec_0f ... */
            vec_0f[i]=phgQuadBasParGradi(e, 0, u_h, i, 4);
            vec_xf[i]=phgQuadBasParGradi(e, 1, u_h, i, 4);
            vec_yf[i]=phgQuadBasParGradi(e, 2, u_h, i, 4);
            vec_zf[i]=phgQuadBasParGradi(e, 3, u_h, i, 4);


        }//endof_for(i = 0; i < N; ++i)

	/* add entries to MAT *F_xx, *F_xy ... */
        for(i = 0;i < N;++i){
            phgMatAddEntries(F_xx, 1, I+i, N, I, &matF_xx[i][0]);
            phgMatAddEntries(F_xy, 1, I+i, N, I, &matF_xy[i][0]);
            phgMatAddEntries(F_xz, 1, I+i, N, I, &matF_xz[i][0]);
            phgMatAddEntries(F_x0, 1, I+i, N, I, &matF_x0[i][0]);
            phgMatAddEntries(F_0x, 1, I+i, N, I, &matF_0x[i][0]);

            phgMatAddEntries(F_yx, 1, I+i, N, I, &matF_yx[i][0]);
            phgMatAddEntries(F_yy, 1, I+i, N, I, &matF_yy[i][0]);
            phgMatAddEntries(F_yz, 1, I+i, N, I, &matF_yz[i][0]);
            phgMatAddEntries(F_y0, 1, I+i, N, I, &matF_y0[i][0]);
            phgMatAddEntries(F_0y, 1, I+i, N, I, &matF_0y[i][0]);

            phgMatAddEntries(F_zx, 1, I+i, N, I, &matF_zx[i][0]);
            phgMatAddEntries(F_zy, 1, I+i, N, I, &matF_zy[i][0]);
            phgMatAddEntries(F_zz, 1, I+i, N, I, &matF_zz[i][0]);
            phgMatAddEntries(F_z0, 1, I+i, N, I, &matF_z0[i][0]);
            phgMatAddEntries(F_0z, 1, I+i, N, I, &matF_0z[i][0]);

            phgMatAddEntries(F_00, 1, I+i, N, I, &matF_00[i][0]);
        }//endof_for(i = 0;i < N;++i)/* add entries to MAT *F_xx, *F_xy ... */

	/* add entries to rhs VEC *rhs_0f ... */
        phgVecAddEntries(rhs_0f, 0, N, I, vec_0f);
        phgVecAddEntries(rhs_xf, 0, N, I, vec_xf);
        phgVecAddEntries(rhs_yf, 0, N, I, vec_yf);
        phgVecAddEntries(rhs_zf, 0, N, I, vec_zf);

	
    }//endof_ForAllElements(g, e)

    k=0;
    for(i=0; i<nY; ++i){
        for(j=0; j<num; ++j){
            *(rhs->data + k) = *(rhs_xf->data + j) * coefRHS_x[i] + *(rhs_yf->data + j) * coefRHS_y[i] +
                *(rhs_zf->data + j) * coefRHS_z[i] + *(rhs_0f->data + j) * coefRHS_0[i];
            ++k;
        }
    }


    phgMapDestroy(&mapu_h);

    phgVecDestroy(&rhs_xf);
    phgVecDestroy(&rhs_yf);
    phgVecDestroy(&rhs_zf);
    phgVecDestroy(&rhs_0f);
}//endof_build_linear_system()
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

