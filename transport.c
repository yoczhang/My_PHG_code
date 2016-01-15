/* Parallel Hierarchical Grid -- an adaptive finite element library.
 *
 * Copyright (C) 2005-2010 State Key Laboratory of Scientific and
 * Engineering Computing, Chinese Academy of Sciences. */

/* This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA */

/* $Id: transport.c,v 1.116 2014/06/16 03:12:52 zlb Exp $
 *
 */

/* 
 * Declare:
 * this transoprt.c is modified through simplest.c
 */

#include "phg.h"
#include "other_functions.h"
#include "int_associated_legendre_polyns.h"
#include "handle_boundary_conditions.h"
#include <string.h>
#include <math.h>  

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
# undef ELEMENT
typedef SIMPLEX ELEMENT;
#endif

#define PN 3  /*  P3 approximate */
#define nY ((PN+1)*(PN+1))
#define Space_basis DOF_P1
    /* 
     * it is easy to compute that in P_N(0,1,...N) approximate the number 
     * of how many Y_l^m is (N+1)*(N+1). 
     */

FLOAT sigma_t = 1.0;
FLOAT sigma_s = 0.5;
FLOAT q_0 = 1.0;
/* the lambda is coefficient of boundary matrixes */
FLOAT lambda=2.0;

/*
 * 1,2 standsfor the faces x-,x+
 * 3,4 standsfor the faces y-,y+
 * 5,6 standsfor the faces z-,z+
 *
 * if the following arry is changed, remember to 
 * modify the function:
 * void build_coefD_xx_bd_( ).
 */
INT fixed_bd_num[3]={1,3,5};
INT reflected_bd_num[3]={2,4,6};

static FLOAT a = 1.0;

static double
elapsed_time(GRID *g, BOOLEAN flag, double mflops)
/* returns and prints elapsed time since last call to this function */
{
    static double t0 = 0.0;
    double et, tt[3];
    size_t mem;
 
    phgGetTime(tt);
    et = tt[2] - t0;
    t0 = tt[2];

    mem = phgMemoryUsage(g, NULL);

    if (flag) {
	if (mflops <= 0)
	    phgPrintf("[%0.4lgMB %0.4lfs]\n", mem / (1024.0 * 1024.0), et);
	else
	    phgPrintf("[%0.4lgMB %0.4lfs %0.4lgGF]\n", mem / (1024.0 * 1024.0),
			 et, mflops*1e-3);
    }

    return et;
}


int
main(int argc, char *argv[])
{
    INT periodicity = 0 /* X_MASK | Y_MASK | Z_MASK */;
    INT mem_max = 900;
    FLOAT tol = 1e-3;
    INT pre_refines = 0;
    INT depth = 1;	/* <=0 => uniformly refine -depth or 3 times */
    INT submesh_threshold = 1000;
    FLOAT lif_threshold = 1.2;

    BOOLEAN tranT=TRUE;
    BOOLEAN tranF=FALSE;
    INT im, jm;
    /* here the Gauss_order standsfor there are 'Gauss_order' Gauss points */
    INT Gauss_order=5;
    FLOAT *Gauss_points_l, *Gauss_points_r;
    FLOAT *Gauss_weights_l, *Gauss_weights_r;

    char *fn = "./mytest_bmarker.mesh";//mytest.mesh is the simplest mesh: is a cube and has only 6 elements.
    GRID *g;
    ELEMENT *e;
    DOF *u_F; // to generate the F_xx matrixes map.
    DOF *u_solver; //to generate the SOLVER *solver. 
    DOF *grad_u, *error;
    SOLVER *solver;
    FLOAT L2error, indicator;
    size_t mem_peak;

    MAP *rmap, *cmap, *rhsmap;

    /*
     * 下面是声明块矩阵中的组成元素，pmatF_xx.. 的含义
     * 同 maxwell-complex.c 的 pmat.
     * 并且下面的 F_xx... 在不强调的情况下都是在四面体
     * 内的体积分.
     * XmFbd_00 则是在边界面X-上的积分，00 则表示不对基函数
     * 求偏导，直接对基函数的积分.
     * 同样，pmatXmFbd_00 则是对应Fbd_00 的.
     * Xm表示X-, Xp表示X+
     */
    MAT *F_xx, *pmatF_xx[nY*nY], *F_xy, *pmatF_xy[nY*nY], *F_xz, *pmatF_xz[nY*nY], *F_x0, *pmatF_x0[nY*nY], *F_0x, *pmatF_0x[nY*nY];
    MAT *F_yx, *pmatF_yx[nY*nY], *F_yy, *pmatF_yy[nY*nY], *F_yz, *pmatF_yz[nY*nY], *F_y0, *pmatF_y0[nY*nY], *F_0y, *pmatF_0y[nY*nY]; 
    MAT *F_zx, *pmatF_zx[nY*nY], *F_zy, *pmatF_zy[nY*nY], *F_zz, *pmatF_zz[nY*nY], *F_z0, *pmatF_z0[nY*nY], *F_0z, *pmatF_0z[nY*nY];
    MAT *F_00, *pmatF_00[nY*nY];

    MAT *XmFbd_00, *pmatXmFbd_00[nY*nY], *XpFbd_00, *pmatXpFbd_00[nY*nY];
    MAT *YmFbd_00, *pmatYmFbd_00[nY*nY], *YpFbd_00, *pmatYpFbd_00[nY*nY];
    MAT *ZmFbd_00, *pmatZmFbd_00[nY*nY], *ZpFbd_00, *pmatZpFbd_00[nY*nY];
    
    VEC *Q=NULL, *rhs=NULL;

    MAT *BlockMat_DxF[16]; //Because in the last there will be 16 blockmatrixes to add up together.
    MAT *matDF; // matDF=BlockMat_DxF[0]+BlockMat_DxF[1]+...+BlockMat_DxF[15];
    MAT *BlockMat_DxF_rhs[4];
    MAT *matDF_rhs;

    /*
     * 注意下面的BlockMat_DxFbd 的大小是 6，这是因为目前所用的网格是一个立方体，
     * 有6个面，每个面都是一个边界
     */
    MAT *BlockMat_DxF_bd[6];
    MAT *matDF_bd;

    /*
     * 根据maxwell-complex.c中分块矩阵的形式，系数矩阵和基矩阵都要排列成一维的向量形式,
     * 具体可以参考 maxwell-complex.c
     */
    FLOAT coefD_xx[nY*nY], coefD_xy[nY*nY], coefD_xz[nY*nY], coefD_x0[nY*nY], coefD_0x[nY*nY];
    FLOAT coefD_yx[nY*nY], coefD_yy[nY*nY], coefD_yz[nY*nY], coefD_y0[nY*nY], coefD_0y[nY*nY];
    FLOAT coefD_zx[nY*nY], coefD_zy[nY*nY], coefD_zz[nY*nY], coefD_z0[nY*nY], coefD_0z[nY*nY];
    FLOAT coefD_00[nY*nY];

    FLOAT coefD_x0_rhs[nY*nY], coefD_y0_rhs[nY*nY], coefD_z0_rhs[nY*nY], coefD_00_rhs[nY*nY];

    FLOAT coefD_Xm_bd[nY*nY], coefD_Xp_bd[nY*nY];
    FLOAT coefD_Ym_bd[nY*nY], coefD_Yp_bd[nY*nY];
    FLOAT coefD_Zm_bd[nY*nY], coefD_Zp_bd[nY*nY];

    /*
     * 接下来声明的变量FLOAT**D_x 是由下面的关系产生的
     * \Omega_x \cdot \vec{Y}^T = \vec{Y}^T \cdot D_x.
     */
    FLOAT **D_x, **D_y, **D_z;

    /* 
     * Compute D_xx=(D_x)^T * D_x
     * ...
     * Dd_x0 stands for D^{\Delta}_{x0} ...
     * I0 stands for the nY*nY identity matrix, 
     * Delta0 stands for the nY*nY \underline{\Delta}^0 matrix.
     */
    FLOAT **D_xx, **D_xy, **D_xz, **D_x0, **Dd_x0, **D_0x, **Dd_0x;
    FLOAT **D_yx, **D_yy, **D_yz, **D_y0, **Dd_y0, **D_0y, **Dd_0y;
    FLOAT **D_zx, **D_zy, **D_zz, **D_z0, **Dd_z0, **D_0z, **Dd_0z;
    FLOAT **Ie, **Delta0;
    FLOAT **temp;


    phgOptionsRegisterFloat("a", "Coefficient", &a);
    phgOptionsRegisterFloat("tol", "Tolerance", &tol);
    phgOptionsRegisterInt("-pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterInt("mem_max", "Maximum memory (MB)", &mem_max);
    phgOptionsRegisterInt("periodicity", "Set periodicity", &periodicity);
    phgOptionsRegisterFloat("-lif_threshold", "LIF threshold", &lif_threshold);
    phgOptionsRegisterInt("-submesh_threshold", "Submesh threshold",
							&submesh_threshold);

    phgInit(&argc, &argv);
    phgOptionsShowUsed();
    if (depth == 0)
	    depth = -3;

    g = phgNewGrid(-1);
    phgSetPeriodicity(g, periodicity);

    if (!phgImport(g, fn, FALSE))
	    phgError(1, "can't read file \"%s\".\n", fn);

    phgRefineAllElements(g, pre_refines);



    
    u_F = phgDofNew(g, Space_basis, 1, "u_F", DofInterpolation);
    u_F->DB_mask = BDRY_MASK;

    u_solver = phgDofNew(g, Space_basis, nY, "u_solver", DofInterpolation);
    // It's very important to comprehend the parameter "nY" in phgDofNew().

    error = phgDofNew(g, DOF_P0, 1, "error", DofNoAction);
    
    D_x=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    D_y=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    D_z=(FLOAT **)malloc(sizeof(FLOAT *)*nY);

    for(im=0;im<nY;++im){
        *(D_x+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(D_y+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(D_z+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
    }

    for(im=0;im<nY;++im){
        for(jm=0;jm<nY;++jm){
            *(*(D_x+im)+jm)=0.0;
            *(*(D_y+im)+jm)=0.0;
            *(*(D_z+im)+jm)=0.0;
        }
    }
 
    build_D_x_matrix(nY, PN, D_x);
    build_D_y_matrix(nY, PN, D_y);
    build_D_z_matrix(nY, PN, D_z);

    D_xx=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    D_xy=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    D_xz=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    D_xz=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    D_x0=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    Dd_x0=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    D_0x=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    Dd_0x=(FLOAT **)malloc(sizeof(FLOAT *)*nY);

    D_yx=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    D_yy=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    D_yz=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    D_y0=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    Dd_y0=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    D_0y=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    Dd_0y=(FLOAT **)malloc(sizeof(FLOAT *)*nY);

    D_zx=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    D_zy=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    D_zz=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    D_z0=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    Dd_z0=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    D_0z=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    Dd_0z=(FLOAT **)malloc(sizeof(FLOAT *)*nY);

    Ie=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    Delta0=(FLOAT **)malloc(sizeof(FLOAT *)*nY);
    
    temp=(FLOAT **)malloc(sizeof(FLOAT *)*nY);


    for(im=0;im<nY;im++){
        *(D_xx+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(D_xy+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(D_xz+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(D_x0+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(Dd_x0+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(D_0x+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(Dd_0x+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
    
        *(D_yx+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(D_yy+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(D_yz+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(D_y0+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(Dd_y0+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(D_0y+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(Dd_0y+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));

        *(D_zx+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(D_zy+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(D_zz+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(D_z0+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(Dd_z0+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(D_0z+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(Dd_0z+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));

        *(Ie+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
        *(Delta0+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));

        *(temp+im)=(FLOAT *)malloc(nY*sizeof(FLOAT));
    }//endof_for(im=0;im<nY;im++)

    for(im=0;im<nY;++im){
        for(jm=0;jm<nY;++jm){
            *(*(Ie+im)+jm)=0.0;
            *(*(Delta0+im)+jm)=0.0;
        }
        *(*(Ie+im)+im)=1.0;
    }
    *(*(Delta0+0)+0)=1.0;


    MatA_multiply_MatB(tranT, nY, D_x, D_x, D_xx);
    MatA_multiply_MatB(tranT, nY, D_x, D_y, D_xy);
    MatA_multiply_MatB(tranT, nY, D_x, D_z, D_xz);
    MatA_multiply_MatB(tranT, nY, D_x, Ie, D_x0);
    MatA_multiply_MatB(tranT, nY, D_x, Delta0, Dd_x0);
    MatA_multiply_MatB(tranF, nY, Ie, D_x, D_0x);
    MatA_multiply_MatB(tranF, nY, Delta0, D_x, Dd_0x);
    
    MatA_multiply_MatB(tranT, nY, D_y, D_x, D_yx);
    MatA_multiply_MatB(tranT, nY, D_y, D_y, D_yy);
    MatA_multiply_MatB(tranT, nY, D_y, D_z, D_yz);
    MatA_multiply_MatB(tranT, nY, D_y, Ie, D_y0);
    MatA_multiply_MatB(tranT, nY, D_y, Delta0, Dd_y0);
    MatA_multiply_MatB(tranF, nY, Ie, D_y, D_0y);
    MatA_multiply_MatB(tranF, nY, Delta0, D_y, Dd_0y);

    MatA_multiply_MatB(tranT, nY, D_z, D_x, D_zx);
    MatA_multiply_MatB(tranT, nY, D_z, D_y, D_zy);
    MatA_multiply_MatB(tranT, nY, D_z, D_z, D_zz);
    MatA_multiply_MatB(tranT, nY, D_z, Ie, D_z0);
    MatA_multiply_MatB(tranT, nY, D_z, Delta0, Dd_z0);
    MatA_multiply_MatB(tranF, nY, Ie, D_z, D_0z);
    MatA_multiply_MatB(tranF, nY, Delta0, D_z, Dd_0z);

    //myDebug
    printf("test10\n"); 
    
    arrangeMatrixInRows(nY,nY,D_xx,coefD_xx);
    /*
     * the coef_xx as the (grad_x f_i, grad_x f_j) 
     * block mattix's coefficients. In PHG's set, the
     * coefficients matrix should be arranged in rows.
     */
    arrangeMatrixInRows(nY,nY,D_xy,coefD_xy);
    arrangeMatrixInRows(nY,nY,D_xz,coefD_xz);
    for(im=0;im<nY;++im){
        for(jm=0;jm<nY;++jm){
            *(*(temp+im)+jm)=*(*(D_x0+im)+jm)*sigma_t - *(*(Dd_x0+im)+jm)*sigma_s;
        }
    }
    arrangeMatrixInRows(nY,nY,temp,coefD_x0);

    arrangeMatrixInRows(nY,nY,D_yx,coefD_yx);
    arrangeMatrixInRows(nY,nY,D_yy,coefD_yy);
    arrangeMatrixInRows(nY,nY,D_yz,coefD_yz);
    for(im=0;im<nY;++im){
        for(jm=0;jm<nY;++jm){
            *(*(temp+im)+jm)=*(*(D_y0+im)+jm)*sigma_t - *(*(Dd_y0+im)+jm)*sigma_s;
        }
    }
    arrangeMatrixInRows(nY,nY,temp,coefD_y0);

    arrangeMatrixInRows(nY,nY,D_zx,coefD_zx);
    arrangeMatrixInRows(nY,nY,D_zy,coefD_zy);
    arrangeMatrixInRows(nY,nY,D_zz,coefD_zz);
    for(im=0;im<nY;++im){
        for(jm=0;jm<nY;++jm){
            *(*(temp+im)+jm)=*(*(D_z0+im)+jm)*sigma_t - *(*(Dd_z0+im)+jm)*sigma_s;
        }
    }
    arrangeMatrixInRows(nY,nY,temp,coefD_z0);

    for(im=0;im<nY;++im){
        for(jm=0;jm<nY;++jm){
            *(*(temp+im)+jm)=*(*(D_0x+im)+jm)*sigma_t - *(*(Dd_0x+im)+jm)*sigma_s;
        }
    }
    arrangeMatrixInRows(nY,nY,temp,coefD_0x);

    for(im=0;im<nY;++im){
        for(jm=0;jm<nY;++jm){
            *(*(temp+im)+jm)=*(*(D_0y+im)+jm)*sigma_t - *(*(Dd_0y+im)+jm)*sigma_s;
        }
    }
    arrangeMatrixInRows(nY,nY,temp,coefD_0y);

    for(im=0;im<nY;++im){
        for(jm=0;jm<nY;++jm){
            *(*(temp+im)+jm)=*(*(D_0z+im)+jm)*sigma_t - *(*(Dd_0z+im)+jm)*sigma_s;
        }
    }
    arrangeMatrixInRows(nY,nY,temp,coefD_0z);

    for(im=0;im<nY;++im){
        for(jm=0;jm<nY;++jm){
            *(*(temp+im)+jm)=*(*(Ie+im)+jm)*sigma_t*sigma_t - 
                *(*(Delta0+im)+jm)*2*(sigma_t*sigma_s-sigma_s*sigma_s);
        }
    }
    arrangeMatrixInRows(nY,nY,temp,coefD_00);

    arrangeMatrixInRows(nY,nY,D_x0,coefD_x0_rhs);
    arrangeMatrixInRows(nY,nY,D_y0,coefD_y0_rhs);
    arrangeMatrixInRows(nY,nY,D_z0,coefD_z0_rhs);
    for(im=0;im<nY;++im){
        for(jm=0;jm<nY;++jm){
            *(*(temp+im)+jm)=*(*(Ie+im)+jm)*sigma_t - *(*(Delta0+im)+jm)*sigma_s;
        }
    }
    arrangeMatrixInRows(nY,nY,temp,coefD_00_rhs);

    Gauss_points_l=(FLOAT *)malloc(Gauss_order * sizeof(FLOAT));
    Gauss_weights_l=(FLOAT *)malloc(Gauss_order * sizeof(FLOAT));
    get_Gauss_points_weights(-1, 0, Gauss_order, Gauss_points_l, Gauss_weights_l);

    Gauss_points_r=(FLOAT *)malloc(Gauss_order * sizeof(FLOAT));
    Gauss_weights_r=(FLOAT *)malloc(Gauss_order * sizeof(FLOAT));
    get_Gauss_points_weights(0, 1, Gauss_order, Gauss_points_r, Gauss_weights_r); 


    /*---------------------------------------------------------------------------------------*/
    // following is the test 

    /*
    FLOAT len_v;
    len_v=0.0;

    len_v=int_associated_legendre_polyns(5,3,5,2,Gauss_order,Gauss_points_l,Gauss_weights_l);
    */

    /*
    for(im=0;im<Gauss_order;im++){
        printf("GaussPoints[%d]=%f ,GaussWeights[%d]=%f \n",im,*(Gauss_points_l+im),im,*(Gauss_weights_l+im));
    }

    printf("P_5^3=%f \n",len_v);
    printf("PI=%f \n",M_PI);
    printf("fllowing is return 0\n");
    return 0;
    */
    

    /*
    FLOAT **C_1;
    int ni=(PN+1)*(PN+2)/2;
    int nj=nY;

    C_1=(FLOAT **)malloc(sizeof(FLOAT *)*ni);
    for(im=0;im<ni;im++){
        *(C_1+im)=(FLOAT *)malloc(nj*sizeof(FLOAT));
    }

    buildMat_fixed_boundary_conditions(PN, 5, Gauss_order, Gauss_points_l, Gauss_weights_l, 
            Gauss_points_r, Gauss_weights_r, C_1);

    for(im=0;im<ni;im++){
        for(jm=0;jm<nj;jm++){
            if(*(*(C_1+im)+jm)!=0.0)
                printf("C_1[%d][%d]=%f   ",im,jm,*(*(C_1+im)+jm));
        }
        printf("\n");
    }
    */

    /*
    FLOAT **C_2X, **C_2Y, **C_2Z;
    int in_XY;
    int in_Z;
    int jn;
    
    in_XY=PN*(PN+1)/2;
    in_Z=PN*(PN+1)/2; // there the in_xy and in_z exactly equivalent.
    jn=nY;

    C_2X=(FLOAT **)malloc(sizeof(FLOAT *)*in_XY);
    C_2Y=(FLOAT **)malloc(sizeof(FLOAT *)*in_XY);
    C_2Z=(FLOAT **)malloc(sizeof(FLOAT *)*in_XY);
    for(im=0;im<in_XY;im++){
        *(C_2X+im)=(FLOAT *)malloc(jn*sizeof(FLOAT));
        *(C_2Y+im)=(FLOAT *)malloc(jn*sizeof(FLOAT));
        *(C_2Z+im)=(FLOAT *)malloc(jn*sizeof(FLOAT));
    }


    buildMat_reflective_boundary_conditions(PN, C_2X, C_2Y, C_2Z);
    */

    /*
    for(im=0;im<in_XY;im++){
        for(jm=0;jm<jn;jm++){
            if(*(*(C_2Z+im)+jm)!=0)
            printf("C_2Z[%d][%d]=%f   ",im,jm,*(*(C_2Z+im)+jm));
        }
        printf("\n\n");
    }
    */

    /*
    printf("DofTypeDim(u_solver)=%d \n",DofTypeDim(u_solver));
    printf("DofDim(u_solver)=%d \n",DofDim(u_solver));

    build_coefD_xx_bd(PN,Gauss_order,Gauss_points_l,Gauss_weights_l,
            Gauss_points_r,Gauss_weights_r,coefD_Xm_bd,coefD_Xp_bd,coefD_Ym_bd,
            coefD_Yp_bd,coefD_Zm_bd,coefD_Zp_bd);


    printf("fllowing is return 0\n");
    return 0;
    */

    // ending the test
    /*---------------------------------------------------------------------------------------*/

    build_coefD_xx_bd(PN,Gauss_order,Gauss_points_l,Gauss_weights_l,
            Gauss_points_r,Gauss_weights_r,coefD_Xm_bd,coefD_Xp_bd,coefD_Ym_bd,
            coefD_Yp_bd,coefD_Zm_bd,coefD_Zp_bd);


    /*
    printf("fllowing is return 0\n");
    return 0;
    */


    while (TRUE) 
    {
	    elapsed_time(g, FALSE, 0.);
        phgPrintf("\n------ %"dFMT" DOF, %"dFMT" elements, mesh LIF = %lg\n",
	        nY * DofGetDataCountGlobal(u_F), g->nleaf_global, (double)g->lif);
        
	    if (phgBalanceGrid(g, lif_threshold, submesh_threshold, NULL, 0.)) {
	        phgPrintf("------ Repartition mesh: nprocs = %d, LIF = %lg ",
			        g->nprocs, (double)g->lif);
	        elapsed_time(g, TRUE, 0.);
	    }//endof_if()
        
	    phgPrintf("Set up linear solver: ");
	    solver = phgSolverCreate(SOLVER_GMRES, u_solver, NULL);
	    solver->mat->handle_bdry_eqns = FALSE;
        
        if(!solver->mat->handle_bdry_eqns)
            printf("solver->mat->handle_bdry_eqns==FALSE \n");
        /* to check whether the solver->mat->handle_bdry_eqns is FLASE
        */
        
	    phgPrintf("solver LIF = %lg \n", (double)solver->rhs->map->lif);


        /*------- creating the Mat *F_xx, *F_xy ... VEC *Q ------*/
        /*-------------------------------------------------------*/
        rmap = phgMapCreate(u_F,NULL);
        cmap = phgMapCreate(u_F,NULL);
        rhsmap = phgMapCreate(u_solver,NULL);

        F_xx = phgMapCreateMat(rmap, cmap);
        F_xx->handle_bdry_eqns = solver->mat->handle_bdry_eqns;
        F_xy = phgMapCreateMat(rmap, cmap);
        F_xy->handle_bdry_eqns = solver->mat->handle_bdry_eqns;
        F_xz = phgMapCreateMat(rmap, cmap);
        F_xz->handle_bdry_eqns = solver->mat->handle_bdry_eqns;
        F_x0 = phgMapCreateMat(rmap, cmap);
        F_x0->handle_bdry_eqns = solver->mat->handle_bdry_eqns;
        F_0x = phgMapCreateMat(rmap, cmap);
        F_0x->handle_bdry_eqns = solver->mat->handle_bdry_eqns;

        F_yx = phgMapCreateMat(rmap, cmap);
        F_yx->handle_bdry_eqns = solver->mat->handle_bdry_eqns;
        F_yy = phgMapCreateMat(rmap, cmap);
        F_yy->handle_bdry_eqns = solver->mat->handle_bdry_eqns;
        F_yz = phgMapCreateMat(rmap, cmap);
        F_yz->handle_bdry_eqns = solver->mat->handle_bdry_eqns;
        F_y0 = phgMapCreateMat(rmap, cmap);
        F_y0->handle_bdry_eqns = solver->mat->handle_bdry_eqns;
        F_0y = phgMapCreateMat(rmap, cmap);
        F_0y->handle_bdry_eqns = solver->mat->handle_bdry_eqns;

        F_zx = phgMapCreateMat(rmap, cmap);
        F_zx->handle_bdry_eqns = solver->mat->handle_bdry_eqns;
        F_zy = phgMapCreateMat(rmap, cmap);
        F_zy->handle_bdry_eqns = solver->mat->handle_bdry_eqns;
        F_zz = phgMapCreateMat(rmap, cmap);
        F_zz->handle_bdry_eqns = solver->mat->handle_bdry_eqns;
        F_z0 = phgMapCreateMat(rmap, cmap);
        F_z0->handle_bdry_eqns = solver->mat->handle_bdry_eqns;
        F_0z = phgMapCreateMat(rmap, cmap);
        F_0z->handle_bdry_eqns = solver->mat->handle_bdry_eqns;

        F_00 = phgMapCreateMat(rmap, cmap);
        F_00->handle_bdry_eqns = solver->mat->handle_bdry_eqns;

        XmFbd_00 = phgMapCreateMat(rmap, cmap);
        XmFbd_00->handle_bdry_eqns = solver->mat->handle_bdry_eqns;
        XpFbd_00 = phgMapCreateMat(rmap, cmap);
        XpFbd_00->handle_bdry_eqns = solver->mat->handle_bdry_eqns;
        YmFbd_00 = phgMapCreateMat(rmap, cmap);
        YmFbd_00->handle_bdry_eqns = solver->mat->handle_bdry_eqns;
        YpFbd_00 = phgMapCreateMat(rmap, cmap);
        YpFbd_00->handle_bdry_eqns = solver->mat->handle_bdry_eqns;
        ZmFbd_00 = phgMapCreateMat(rmap, cmap);
        ZmFbd_00->handle_bdry_eqns = solver->mat->handle_bdry_eqns;
        ZpFbd_00 = phgMapCreateMat(rmap, cmap);
        ZpFbd_00->handle_bdry_eqns = solver->mat->handle_bdry_eqns;

        Q = phgMapCreateVec(rhsmap, 1);
        /*-------------------------------------------------------*/
        /*------- creating the Mat *F_xx, *F_xy ... VEC *Q ------*/
        //myDebug
        printf("test11\n");



        /*------------ assemble the block matrix ... ------------*/
        /*-------------------------------------------------------*/
        for(im=0; im<nY*nY; ++im){
            pmatF_xx[im]=F_xx;
            pmatF_xy[im]=F_xy;
            pmatF_xz[im]=F_xz;
            pmatF_x0[im]=F_x0;
            pmatF_0x[im]=F_0x;

            pmatF_yx[im]=F_yx;
            pmatF_yy[im]=F_yy;
            pmatF_yz[im]=F_yz;
            pmatF_y0[im]=F_y0;
            pmatF_0y[im]=F_0y;

            pmatF_zx[im]=F_zx;
            pmatF_zy[im]=F_zy;
            pmatF_zz[im]=F_zz;
            pmatF_z0[im]=F_z0;
            pmatF_0z[im]=F_0z;

            pmatF_00[im]=F_00;

            pmatXmFbd_00[im]=XmFbd_00;
            pmatXpFbd_00[im]=XpFbd_00;
            pmatYmFbd_00[im]=YmFbd_00;
            pmatYpFbd_00[im]=YpFbd_00;
            pmatZmFbd_00[im]=ZmFbd_00;
            pmatZpFbd_00[im]=ZpFbd_00;
        }//endof_for(im=0;...)
        
        BlockMat_DxF[0]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_xx, coefD_xx, NULL);
        BlockMat_DxF[1]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_xy, coefD_xy, NULL);
        BlockMat_DxF[2]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_xz, coefD_xz, NULL);
        BlockMat_DxF[3]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_x0, coefD_x0, NULL);
        BlockMat_DxF[4]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_0x, coefD_0x, NULL);

        BlockMat_DxF[5]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_yx, coefD_yx, NULL);
        BlockMat_DxF[6]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_yy, coefD_yy, NULL);
        BlockMat_DxF[7]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_yz, coefD_yz, NULL);
        BlockMat_DxF[8]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_y0, coefD_y0, NULL);
        BlockMat_DxF[9]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_0y, coefD_0y, NULL);

        BlockMat_DxF[10]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_zx, coefD_zx, NULL);
        BlockMat_DxF[11]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_zy, coefD_zy, NULL);
        BlockMat_DxF[12]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_zz, coefD_zz, NULL);
        BlockMat_DxF[13]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_z0, coefD_z0, NULL);
        BlockMat_DxF[14]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_0z, coefD_0z, NULL);       

        BlockMat_DxF[15]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_00, coefD_00, NULL);

        BlockMat_DxF_rhs[0]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_x0, coefD_x0_rhs, NULL);
        BlockMat_DxF_rhs[1]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_y0, coefD_y0_rhs, NULL);
        BlockMat_DxF_rhs[2]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_z0, coefD_z0_rhs, NULL);
        BlockMat_DxF_rhs[3]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatF_00, coefD_00_rhs, NULL);

        if(coefD_Xm_bd!=NULL)
            BlockMat_DxF_bd[0]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatXmFbd_00, coefD_Xm_bd, NULL);
        else
            BlockMat_DxF_bd[0]=NULL;

        if(coefD_Xp_bd!=NULL)
            BlockMat_DxF_bd[1]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatXpFbd_00, coefD_Xp_bd, NULL);
        else
            BlockMat_DxF_bd[1]=NULL;

        if(coefD_Ym_bd!=NULL)
            BlockMat_DxF_bd[2]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatYmFbd_00, coefD_Ym_bd, NULL);
        else
            BlockMat_DxF_bd[2]=NULL;

        if(coefD_Yp_bd!=NULL)
            BlockMat_DxF_bd[3]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatYpFbd_00, coefD_Yp_bd, NULL);
        else
            BlockMat_DxF_bd[3]=NULL;

        if(coefD_Zm_bd!=NULL)
            BlockMat_DxF_bd[4]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatZmFbd_00, coefD_Zm_bd, NULL);
        else
            BlockMat_DxF_bd[4]=NULL;

        if(coefD_Zp_bd!=NULL)
            BlockMat_DxF_bd[5]=phgMatCreateBlockMatrix(g->comm, nY, nY, pmatZpFbd_00, coefD_Zp_bd, NULL);
        else
            BlockMat_DxF_bd[5]=NULL;
        /*-------------------------------------------------------*/
        /*------------ assemble the block matrix ... ------------*/


        /*--------- Add the block matrixes to MAT *matDF -------*/
        /*------------------------------------------------------*/
        matDF=phgMatAXPBY(1.0,BlockMat_DxF[0],0.0,&matDF);
        for(im=1;im<16;++im){
            matDF=phgMatAXPBY(1.0,BlockMat_DxF[im],1.0,&matDF);
        }
        /*------------------------------------------------------*/
        /*--------- Add the block matrixes to MAT *matDF -------*/



        /*------- Add the block matrixes to MAT *matDF_rhs -----*/
        /*------------------------------------------------------*/
        matDF_rhs=phgMatAXPBY(1.0,BlockMat_DxF_rhs[0],0.0,&matDF_rhs);
        for(im=1;im<4;++im){
            matDF_rhs=phgMatAXPBY(1.0,BlockMat_DxF_rhs[im],1.0,&matDF_rhs);
        }
        /*------------------------------------------------------*/
        /*------- Add the block matrixes to MAT *matDF_rhs -----*/


        /*------- Add the block matrixes to MAT *matDF_bd  -----*/
        /*------------------------------------------------------*/
        /*
        for(im=0;im<6;++im){
            if(BlockMat_DxF_bd[im]!=NULL){
                printf("add to matDF_bd , BlockMat_DxF_bd[%d]!=NULL\n",im);
                matDF_bd=phgMatAXPBY(1.0,BlockMat_DxF_bd[im],0.0,&matDF_bd);
                im++;
            }
            matDF_bd=phgMatAXPBY(1.0,BlockMat_DxF_bd[im],1.0,&matDF_rhs);
        }
        */
        jm=0;
        for(im=0;im<6;++im){
            if(BlockMat_DxF_bd[im]!=NULL){
                jm=im;
                break;
            }
        }

        printf("jm=%d\n",jm);
        matDF_bd=phgMatAXPBY(1.0,BlockMat_DxF_bd[jm],0.0,&matDF_bd);
        for(im=jm+1;im<6;++im){
            matDF_bd=phgMatAXPBY(1.0,BlockMat_DxF_bd[im],1.0,&matDF_rhs);
        }
        /*------------------------------------------------------*/
        /*------- Add the block matrixes to MAT *matDF_bd  -----*/

 

        /*------------ Disassemble the MAT *F_xx ... -----------*/
        /*------------------------------------------------------*/
        /* Through the operation phgMatAXPBY(), here we must 
         * disassemble the MAT *F_xx ... before assemble them.
         * */
        phgMatDisassemble(F_xx);
        phgMatDisassemble(F_xy);
        phgMatDisassemble(F_xz);
        phgMatDisassemble(F_x0);
        phgMatDisassemble(F_0x);

        phgMatDisassemble(F_yx);
        phgMatDisassemble(F_yy);
        phgMatDisassemble(F_yz);
        phgMatDisassemble(F_y0);
        phgMatDisassemble(F_0y);

        phgMatDisassemble(F_zx);
        phgMatDisassemble(F_zy);
        phgMatDisassemble(F_zz);
        phgMatDisassemble(F_z0);
        phgMatDisassemble(F_0z);

        phgMatDisassemble(F_00); 

        phgMatDisassemble(XmFbd_00);
        phgMatDisassemble(XpFbd_00);
        phgMatDisassemble(YmFbd_00);
        phgMatDisassemble(YpFbd_00);
        phgMatDisassemble(ZmFbd_00);
        phgMatDisassemble(ZpFbd_00);
        /*------------------------------------------------------*/
        /*------------ Disassemble the MAT *F_xx ... -----------*/

        //printf("next is return 0\n");
        //return 0;

        phgMatDestroy(&solver->mat);
        /* 下面是将原矩阵 matDF 与边界矩阵 matDF_bd 相加 */
        matDF=phgMatAXPBY(lambda,matDF_bd,1.0,&matDF); 
        solver->mat=matDF;
        solver->rhs->mat=solver->mat;

        phgPrintf("Assemble the matrixes: F_xx, F_xy ... \n");
        /*
        assemble_Fxx_matrixes(F_xx, F_xy, F_xz, F_x0, F_0x, 
                F_yx, F_yy, F_yz, F_y0, F_0y, 
                F_zx, F_zy, F_zz, F_z0, F_0z, 
                F_00);
        */



        /*------------------ Build the RHS ---------------------*/
        /*------------------------------------------------------*/
        /*
         * According the equation's RHS expression, we can use the
         * matrixes which we have built above.
         */
        //printf("Build the RHS \n");

        /*--------------- test the length of VEC *Q -------------*/
        /*-------------------------------------------------------*/
        int lengthQ;
        lengthQ=Q->map->nlocal*Q->nvec;
        //lengthQ=Q->map->nlocal;
        printf("length of Q = %d \n",lengthQ);

        if(lengthQ != nY * DofGetDataCountGlobal(u_F)){
            printf("the length of Q is not matching the total number of DOFs(nY * DofGetDataCountGlobal(u_F)), quit!");
        }

        assert(lengthQ == nY * DofGetDataCountGlobal(u_F));

        printf("build the RHS \n");
        //build_rhs(D_x, D_y, D_z, u_F, rhs, nY);

        for(im=0;im<lengthQ;++im){
            *(Q->data+im)=1.0;
        }//assignment to VEC *Q.
        /*-------------------------------------------------------*/
        /*--------------- test the length of VEC *Q -------------*/

        rhs=solver->rhs;

        int lengthrhs;
        lengthrhs=rhs->map->nlocal*rhs->nvec;
        //lengthQ=Q->map->nlocal;
        printf("length of rhs = %d \n",lengthrhs);

        /*
        build_rhs(D_x, D_y, D_z, u_F, rhs, nY);
        */

        //rhs=phgMatVec(0, 1.0, matDF_rhs, Q, 0.0, &rhs);
        /* rhs= matDF_rhs*Q */
        /*------------------------------------------------------*/
        /*------------------ Build the RHS ---------------------*/
       /* 
        for(im=0;im<lengthQ;++im){
            printf("%d.rhs=%f   \n",im,*(solver->rhs->data+im));
        }

        /*
        printf("\n");
        for(im=0;im<lengthQ;++im){
            printf("rhs=%f   \n",*(Q->data+im));
        }
        */
        build_linear_system(F_xx, F_xy, F_xz, F_x0, F_0x, F_yx, F_yy, F_yz, 
                F_y0, F_0y, F_zx, F_zy, F_zz, F_z0, F_0z, F_00, XmFbd_00, 
                XpFbd_00, YmFbd_00, YpFbd_00, ZmFbd_00, ZpFbd_00, D_x, D_y, 
                D_z, rhs, nY);



        elapsed_time(g,TRUE,0.);

        phgPrintf("Solve linear system: ");
        phgSolverSolve(solver,TRUE,u_solver,NULL);
        phgPrintf("nits=%d, resid=%0.4lg ",solver->nits,
                (double)solver->residual);
        phgSolverDestroy(&solver);
        elapsed_time(g,TRUE,0.);



        printf("next is return 0\n");
        return 0;

        phgMatDestroy(&F_xx);
        phgMatDestroy(&F_xy);
        phgMatDestroy(&F_xz);
        phgMatDestroy(&F_x0);
        phgMatDestroy(&F_0x);

        phgMatDestroy(&F_yx);
        phgMatDestroy(&F_yy);
        phgMatDestroy(&F_yz);
        phgMatDestroy(&F_y0);
        phgMatDestroy(&F_0y);

        phgMatDestroy(&F_zx);
        phgMatDestroy(&F_zy);
        phgMatDestroy(&F_zz);
        phgMatDestroy(&F_z0);
        phgMatDestroy(&F_0z);

        phgMatDestroy(&F_00);

        phgVecDestroy(&Q);
    }//endof_while(TRUE)


#if 0
    phgPrintf("Final mesh written to \"%s\".\n",
	phgExportDX(g, "simplest.dx", u_h, error, NULL));
#elif 0
    phgPrintf("Final mesh written to \"%s\".\n",
	phgExportVTK(g, "simplest.vtk", u_h, error, NULL));
#endif

    free(coefD_xx); free(coefD_xy); free(coefD_xz); free(coefD_x0); free(coefD_0x);
    free(coefD_yx); free(coefD_yy); free(coefD_yz); free(coefD_y0); free(coefD_0y);
    free(coefD_zx); free(coefD_zy); free(coefD_zz); free(coefD_z0); free(coefD_0z);
    free(coefD_00); 
    free(coefD_x0_rhs); free(coefD_y0_rhs); free(coefD_z0_rhs); free(coefD_00_rhs);

    free(coefD_Xm_bd); free(coefD_Xp_bd); 
    free(coefD_Ym_bd); free(coefD_Yp_bd);
    free(coefD_Zm_bd); free(coefD_Zp_bd);

    phgDofFree(&u_F);
    phgDofFree(&u_solver);
    phgDofFree(&error);

    phgFreeGrid(&g);

    phgFinalize();

    return 0;
}
