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
    INT Gauss_order=5;
    FLOAT *Gauss_points_l, *Gauss_points_r;
    FLOAT *Gauss_weights_l, *Gauss_weights_r;

    char *fn = "./mytest.mesh";//mytest.mesh is the simplest mesh: is a cube and has only 6 elements.
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
     * 注意下面的BlockMat_DxFbd 的大小是 6，这是因为目前所用的网格是一个立方体，
     * 有6个面，每个面都是一个边界
     */
    MAT *BlockMat_DxF_bd[6];
    MAT *matDF_bd;



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
   

    for(im=0;im<nY;++im){
        for(jm=0;jm<nY;++jm){
            *(*(Ie+im)+jm)=0.0;
            *(*(Delta0+im)+jm)=0.0;
        }
        *(*(Ie+im)+im)=1.0;
    }
    *(*(Delta0+0)+0)=1.0;


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


    printf("fllowing is return 0\n");
    return 0;



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
