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
#include <string.h>
#include <math.h>  

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
# undef ELEMENT
typedef SIMPLEX ELEMENT;
#endif

#define PN 3  /*  P3 approximate */
#define nY ((PN+1)*(PN+1))
    /* 
     * it is easy to compute that in P_N(0,1,...N) approximate the number 
     * of how many Y_l^m is (N+1)*(N+1). 
     */
#define sigma_t 1.0
#define sigma_s 0.6

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

    //myDebug
    printf("test0\n");

    char *fn = "./mytest.mesh";//mytest.mesh is the simplest mesh: is a cube and has only 6 elements.
    GRID *g;
    ELEMENT *e;
    DOF *u_F; // to generate the F_xx matrixes map.
    DOF *u_solver; //to generate the SOLVER *solver. 
    DOF *rhs, *grad_u, *error;
    SOLVER *solver;
    FLOAT L2error, indicator;
    size_t mem_peak;
    MAT *F_xx, *pmatF_xx[nY*nY], *F_xy, *pmatF_xy[nY*nY], *F_xz, *pmatF_xz[nY*nY], *F_x0, *pmatF_x0[nY*nY], *F_0x, *pmatF_0x[nY*nY];
    MAT *F_yx, *pmatF_yx[nY*nY], *F_yy, *pmatF_yy[nY*nY], *F_yz, *pmatF_yz[nY*nY], *F_y0, *pmatF_y0[nY*nY], *F_0y, *pmatF_0y[nY*nY]; 
    MAT *F_zx, *pmatF_zx[nY*nY], *F_zy, *pmatF_zy[nY*nY], *F_zz, *pmatF_zz[nY*nY], *F_z0, *pmatF_z0[nY*nY], *F_0z, *pmatF_0z[nY*nY];
    MAT *F_00, *pmatF_00[nY*nY];

    //myDebug
    printf("test1\n");

    MAT *BlockMat_DxF[16]; //Because in the last there will be 16 blockmatrixes to add up together.
    MAT *matDF; // matDF=BlockMat_DxF[0]+BlockMat_DxF[1]+...+BlockMat_DxF[15];

    FLOAT coefD_xx[nY*nY], coefD_xy[nY*nY], coefD_xz[nY*nY], coefD_x0[nY*nY], coefD_0x[nY*nY];
    FLOAT coefD_yx[nY*nY], coefD_yy[nY*nY], coefD_yz[nY*nY], coefD_y0[nY*nY], coefD_0y[nY*nY];
    FLOAT coefD_zx[nY*nY], coefD_zy[nY*nY], coefD_zz[nY*nY], coefD_z0[nY*nY], coefD_0z[nY*nY];
    FLOAT coefD_00[nY*nY];

    MAP *rmap, *cmap;

    //myDebug
    printf("test2\n");

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

    //myDebug
    printf("test3\n");
    
    u_F = phgDofNew(g, DOF_DEFAULT, 1, "u_F", DofInterpolation);
    u_F->DB_mask = BDRY_MASK;

    //myDebug
    printf("test4\n");

    u_solver = phgDofNew(g, DOF_DEFAULT, nY, "u_solver", DofInterpolation);
    // It's very important to comprehend the parameter "nY" in phgDofNew().
    
    //myDebug
    printf("test5\n");

    rhs = phgDofNew(g, DOF_CONSTANT, 1, "rhs", DofNoAction);
    phgDofSetDataByValuesV(rhs, (FLOAT)1.0);

    error = phgDofNew(g, DOF_P0, 1, "error", DofNoAction);

    //myDebug
    printf("test6\n");

    rmap = phgMapCreate(u_F,NULL);
    cmap = phgMapCreate(u_F,NULL);

    //myDebug
    printf("test7\n"); 

    
    
    int im, jm;
    FLOAT **D_x, **D_y, **D_z;
    
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

    //myDebug
    printf("test8\n");
 
    build_D_x_matrix(nY, PN, D_x);
    build_D_y_matrix(nY, PN, D_y);
    build_D_z_matrix(nY, PN, D_z);

    //myDebug
    printf("test9\n");

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

    BOOLEAN tranT=TRUE;
    BOOLEAN tranF=FALSE;

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
    
    ArrangMatrixInRows(nY,nY,D_xx,coefD_xx);
    /*
     * the coef_xx as the (grad_x f_i, grad_x f_j) 
     * block mattix's coefficients. In PHG's set, the
     * coefficients matrix should be arranged in rows.
     */
    ArrangMatrixInRows(nY,nY,D_xy,coefD_xy);
    ArrangMatrixInRows(nY,nY,D_xz,coefD_xz);
    for(im=0;im<nY;++im){
        for(jm=0;jm<nY;++jm){
            *(*(temp+im)+jm)=*(*(D_x0+im)+jm)*sigma_t - *(*(Dd_x0+im)+jm)*sigma_s;
        }
    }
    ArrangMatrixInRows(nY,nY,temp,coefD_x0);

    ArrangMatrixInRows(nY,nY,D_yx,coefD_yx);
    ArrangMatrixInRows(nY,nY,D_yy,coefD_yy);
    ArrangMatrixInRows(nY,nY,D_yz,coefD_yz);
    for(im=0;im<nY;++im){
        for(jm=0;jm<nY;++jm){
            *(*(temp+im)+jm)=*(*(D_y0+im)+jm)*sigma_t - *(*(Dd_y0+im)+jm)*sigma_s;
        }
    }
    ArrangMatrixInRows(nY,nY,temp,coefD_y0);

    ArrangMatrixInRows(nY,nY,D_zx,coefD_zx);
    ArrangMatrixInRows(nY,nY,D_zy,coefD_zy);
    ArrangMatrixInRows(nY,nY,D_zz,coefD_zz);
    for(im=0;im<nY;++im){
        for(jm=0;jm<nY;++jm){
            *(*(temp+im)+jm)=*(*(D_z0+im)+jm)*sigma_t - *(*(Dd_z0+im)+jm)*sigma_s;
        }
    }
    ArrangMatrixInRows(nY,nY,temp,coefD_z0);

    for(im=0;im<nY;++im){
        for(jm=0;jm<nY;++jm){
            *(*(temp+im)+jm)=*(*(D_0x+im)+jm)*sigma_t - *(*(Dd_0x+im)+jm)*sigma_s;
        }
    }
    ArrangMatrixInRows(nY,nY,temp,coefD_x0);

    for(im=0;im<nY;++im){
        for(jm=0;jm<nY;++jm){
            *(*(temp+im)+jm)=*(*(D_0y+im)+jm)*sigma_t - *(*(Dd_0y+im)+jm)*sigma_s;
        }
    }
    ArrangMatrixInRows(nY,nY,temp,coefD_0y);

    for(im=0;im<nY;++im){
        for(jm=0;jm<nY;++jm){
            *(*(temp+im)+jm)=*(*(D_0z+im)+jm)*sigma_t - *(*(Dd_0z+im)+jm)*sigma_s;
        }
    }
    ArrangMatrixInRows(nY,nY,temp,coefD_0z);

    for(im=0;im<nY;++im){
        for(jm=0;jm<nY;++jm){
            *(*(temp+im)+jm)=*(*(Ie+im)+jm)*sigma_t*sigma_t - 
                *(*(Delta0+im)+jm)*2*(sigma_t*sigma_s-sigma_s*sigma_s);
        }
    }
    ArrangMatrixInRows(nY,nY,temp,coefD_00);

    //printf("next is return 0\n");
    //return 0;

    while (TRUE) 
    {
	    elapsed_time(g, FALSE, 0.);
        phgPrintf("\n------ %"dFMT" DOF, %"dFMT" elements, mesh LIF = %lg\n",
	        2 * DofGetDataCountGlobal(u_F), g->nleaf_global, (double)g->lif);
        
	    if (phgBalanceGrid(g, lif_threshold, submesh_threshold, NULL, 0.)) {
	        phgPrintf("------ Repartition mesh: nprocs = %d, LIF = %lg ",
			        g->nprocs, (double)g->lif);
	        elapsed_time(g, TRUE, 0.);
	    }//endof_if()
        
	    phgPrintf("Set up linear solver: ");
	    solver = phgSolverCreate(SOLVER_DEFAULT, u_solver, NULL);
	    solver->mat->handle_bdry_eqns = FALSE;
        
        if(!solver->mat->handle_bdry_eqns)
            printf("solver->mat->handle_bdry_eqns==FALSE \n");
        /* to check whether the solver->mat->handle_bdry_eqns is FLASE
        */
        
	    phgPrintf("solver LIF = %lg ", (double)solver->rhs->map->lif);


        /*---------- creating the Mat *F_xx, *F_xy ... ----------*/
        /*-------------------------------------------------------*/
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
        /*-------------------------------------------------------*/
        /*---------- creating the Mat *F_xx, *F_xy ... ----------*/
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
        /*------------ Disassemble the MAT *F_xx ... -----------*/
        /*------------------------------------------------------*/

        //printf("next is return 0\n");
        //return 0;

        phgMatDestroy(&solver->mat);
        solver->mat=matDF;
        solver->rhs->mat=solver->mat;

        phgPrintf("Assemble the matrixes: F_xx, F_xy ... ");
        Assemble_Fxx_matrixes(F_xx, F_xy, F_xz, F_x0, F_0x, 
                F_yx, F_yy, F_yz, F_y0, F_0y, 
                F_zx, F_zy, F_zz, F_z0, F_0z, 
                F_00);

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
    }//endof_while(TRUE)


#if 0
    phgPrintf("Final mesh written to \"%s\".\n",
	phgExportDX(g, "simplest.dx", u_h, error, NULL));
#elif 0
    phgPrintf("Final mesh written to \"%s\".\n",
	phgExportVTK(g, "simplest.vtk", u_h, error, NULL));
#endif

    phgDofFree(&u_F);
    phgDofFree(&u_solver);
    phgDofFree(&error);

    phgFreeGrid(&g);

    phgFinalize();

    return 0;
}
