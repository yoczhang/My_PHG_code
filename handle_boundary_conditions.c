/*
 * The fix_boundary_conditions.c is aimed to deal the 
 * boundary conditions 
 */

#include "phg.h"
#include "other_functions.h"
#include "int_associated_legendre_polyns.h"
#include "handle_boundary_conditions.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <malloc.h>

/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/*
 * To compute the coefficient D(m)
 */
FLOAT 
D(int m)
{
    int flag;
    int delta;
    FLOAT v;

    m=abs(m);

    if(m%2==0)
        flag=1;
    else
        flag=-1;

    if(m==0)
        delta=1;
    else
        delta=0;

    v=flag*sqrt((FLOAT)(2-delta)/(4*M_PI));

    return v;
}//endof_D()

/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/*
 * To compute Ip(l, m, l1, m1)
 */
FLOAT
Ip(int l, int m, int l1, int m1, int Gauss_order, FLOAT *Gauss_points, 
        FLOAT *Gauss_weights)
{
    FLOAT v;
    v=int_associated_legendre_polyns(l,m,l1,m1,Gauss_order,Gauss_points,Gauss_weights);
    return v;
}//endof_Ip()
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/




/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
FLOAT
I(int l, int m, int l1, int m1, int lable, int Gauss_order, FLOAT *Gauss_points_l, 
        FLOAT *Gauss_weights_l, FLOAT *Gauss_points_r, FLOAT *Gauss_weights_r)
/*
 * int lable standsfor the different faces.
 * lable=1,2 standsfor the faces x-,x+
 * lable=3,4 standsfor the faces y-,y+
 * lable=5,6 standsfor the faces z-,z+
 */
{
    FLOAT v;
    FLOAT *Gauss_points, *Gauss_weights;

    if(lable==1 || lable==2 || lable==3 || lable==4){
        if(l==l1 && m==m1 && m!=0)
            v=D(m)*D(m1)*(M_PI);
        else if(l==l1 && m==m1 && m==0)
            v=D(m)*D(m1)*2*(M_PI);
        else
            v=0.0;

        return v;
    }
    else if(lable==5){
        Gauss_points=Gauss_points_r;
        Gauss_weights=Gauss_weights_r;
        if(m==m1 && m!=0)
            v=D(m)*D(m1)*(M_PI)*Ip(l,m,l1,m1,Gauss_order,Gauss_points,Gauss_weights);
        else if(m==m1 && m==0)
            v=D(m)*D(m1)*(2*M_PI)*Ip(l,m,l1,m1,Gauss_order,Gauss_points,Gauss_weights);
        else
            v=0.0;

        return v;
    }
    else if(lable==6){
        Gauss_points=Gauss_points_l;
        Gauss_weights=Gauss_weights_l;
        if(m==m1 && m!=0)
            v=D(m)*D(m1)*(M_PI)*Ip(l,m,l1,m1,Gauss_order,Gauss_points,Gauss_weights);
        else if(m==m1 && m==0)
            v=D(m)*D(m1)*(2*M_PI)*Ip(l,m,l1,m1,Gauss_order,Gauss_points,Gauss_weights);
        else
            v=0.0;

        return v;
    }
    else{
        printf("In the function 'I()', the 'lable' is wrong! exit!\n");
        exit(-1);
    }
}
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
void
buildMat_fixed_boundary_conditions(int N, int lable, int Gauss_order, 
        FLOAT *Gauss_points_l, FLOAT *Gauss_weights_l, FLOAT *Gauss_points_r, 
        FLOAT *Gauss_weights_r, FLOAT **C_1)
/*
 * This function is aimed to build the fixed boundary conditions' matrixes.
 *
 * Input: N, the PN approximating number, is odd.
 *        lable, the lable of the boundary face.
 *        Gauss_order, the number of Gausspoints.
 *        Gauss_points_l, the Gausspoints on interval [-1,0].
 *        Gauss_weights_l, the Gaussweights according Gauss_points_l.
 *        Gauss_points_r, the Gausspoints on interval [0,1].
 *        Gauss_weights_r, the Gaussweights according Gauss_points_r.
 *        
 * Output: C_1, the matrix based on the fixed boundary face.
 */
{
    int i, j;
    int l, m, l1, m1;

    //l1 is odd, so l1 initialization is 1.
    for(l1=1;l1<=N;l1+=2){
        for(m1=-l1;m1<=l1;m1++){
            i=l1*(l1+1)/2+m1;
                
            for(l=0;l<=N;l++){
                for(m=-l;m<=l;m++){
                    j=l*l+l+m;
                    *(*(C_1+i)+j)=I(l, m, l1, m1, lable, Gauss_order, Gauss_points_l, 
                            Gauss_weights_l, Gauss_points_r, Gauss_weights_r);
                }//endof_for(m=-l;m<=l;m++)
            }//endof_for(l=0;l<=N;l++)

        }//endof_for(m=-l1;m<=l1;m++)
    }//endof_for(l1=1;l1<=N;l1+=2)

}//endof_fixed_boundary_conditins()

/**********************************************************************************/
/*--------------------------------------------------------------------------------*/




/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
void
buildMat_reflective_boundary_conditions(int N, FLOAT **C_2X, FLOAT **C_2Y, FLOAT **C_2Z)
/*
 * This function is aimed to build the reflective boundary conditions' matrixes.
 *
 * Input: N, the PN approximating number, is odd.
 *
 * Output: just notice that the matrix on face X- or X+ is the same, the same to Y and Z.
 *         C_2X, the matrix based on the reflective boundary faces x- or x+.
 *         C_2Y, the matrix based on the reflective boundary faces y- or y+.
 *         C_2Z, the matrix based on the reflective boundary faces z- or z+.
 */
{
    int i, j;
    int l, m;

    int leap_pre, leap_now;

    int in_XY;// the number of rows based on x or y boundary faces.
    int in_Z; // the number of rows based on z boundary faces.
    int jn;   // the number of columns.

    in_XY=N*(N+1)/2;
    in_Z=N*(N+1)/2; // there the in_xy and in_z exactly equivalent.
    jn=(N+1)*(N+1);

    //initialization
    for(i=0;i<in_XY;i++){
        for(j=0;j<jn;j++){
            *(*(C_2X+i)+j)=0.0;
            *(*(C_2Y+i)+j)=0.0;
            *(*(C_2Z+i)+j)=0.0;
        }
    }//endof_for(i=0;...)

    /* assignment the C_2X matrix */
    for(l=1;l<=N;l++){
        for(m=1;m<=l;m++){
            i=l*(l-1)/2+m-1;

            if(m%2==1)
                j=l*l+l+m;
            else
                j=l*l+l+(-m);

            *(*(C_2X+i)+j)=1.0;
        }
    }

    /* assignment the C_2Y matrix */
    for(l=1;l<=N;l++){
        for(m=1;m<=l;m++){
            i=l*(l-1)/2+m-1;
            j=l*l+l+(-m);
            *(*(C_2Y+i)+j)=1.0;
        }
    }

    /* assignment the C_2Z matrix */
    for(l=1;l<=N;l++){
        leap_pre=0;
        leap_now=0;

        for(m=-l;m<=l;m++){
            if((l+m)%2==1)
                leap_now++;

            if(leap_now-leap_pre==1){
                i=l*(l-1)/2+(leap_now-1);
                j=l*l+l+m;

                *(*(C_2Z+i)+j)=1.0;

                leap_pre=leap_now;
            }

        }//endof_for(m=-l;...)
    }//endof_for(l=1;...)

}//endof_buildMat_reflective_boundary_conditions()

/**********************************************************************************/
/*--------------------------------------------------------------------------------*/





/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
void
build_coefD_xx_bd(int PN, int Gauss_order, FLOAT *Gauss_points_l, FLOAT *Gauss_weights_l, 
        FLOAT *Gauss_points_r, FLOAT *Gauss_weights_r, 
        FLOAT *coefD_Xm_bd, FLOAT *coefD_Xp_bd, FLOAT *coefD_Ym_bd, 
        FLOAT *coefD_Yp_bd, FLOAT *coefD_Zm_bd, FLOAT *coefD_Zp_bd)
/*
 * This function is aimed to build the coefD_Xm_bd, ...
 * which is coefficient matrixes of boundary faces used in the 
 * block matrix, but coefD_Xm_bd is a vec arranged by rows of 
 * the matrix. More details see "程序过程整理1" page:39.
 *
 * Input: PN, the PN approximating number, is odd.
 *        lable, the lable of the boundary face.
 *        Gauss_order, the number of Gausspoints.
 *        Gauss_points_l, the Gausspoints on interval [-1,0].
 *        Gauss_weights_l, the Gaussweights according Gauss_points_l.
 *        Gauss_points_r, the Gausspoints on interval [0,1].
 *        Gauss_weights_r, the Gaussweights according Gauss_points_r.
 * 
 * Output: coefD_Xm_bd, the X- boundary face coefficient vec. 
 *         coefD_Xp_bd, the X+ boundary face coefficient vec. 
 *         ...
 *         ...
 */
{
    extern INT fixed_bd_num[3];
    extern INT reflected_bd_num[3];

    int len_fixed=sizeof(fixed_bd_num)/sizeof(fixed_bd_num[0]);
    int len_refl=sizeof(reflected_bd_num)/sizeof(reflected_bd_num[0]);
    int i,j,nY;
    int fixed_bd_row, fixed_bd_col, refl_bd_row, refl_bd_col;
    FLOAT **coefXm, **coefXp, **coefYm, **coefYp, **coefZm, **coefZp;
    FLOAT **refl_coefX_cache, **refl_coefY_cache, **refl_coefZ_cache;
    FLOAT **matAxmatB;// 用来存储矩阵 A*B 的结果.
    int fg_coefXm, fg_coefXp, fg_coefYm, fg_coefYp, fg_coefZm, fg_coefZp;
    /* fg_coefXm 用来标识 coefXm 有没有用到，可能并不是所有面都会用到. */
    
    nY=(PN+1)*(PN+1);

    fixed_bd_row=(PN+1)*(PN+2)/2;
    fixed_bd_col=nY;

    refl_bd_row=PN*(PN+1)/2;   /* 这里要注意，其实对于 refl_bd_row X 和 Y 是一样的，但  */
                               /* 对于 Z 的两个面，仅仅也恰好是这个值.                  */
    refl_bd_col=nY;    

    /* refl_coefX_cache ... 用来临时存储 reflective 边界面形成的
     * 矩阵，在接下来的程序中用来给 coefXm ...进行赋值。*/
    refl_coefX_cache=(FLOAT **)malloc(sizeof(FLOAT *)*refl_bd_row);
    refl_coefY_cache=(FLOAT **)malloc(sizeof(FLOAT *)*refl_bd_row);
    refl_coefZ_cache=(FLOAT **)malloc(sizeof(FLOAT *)*refl_bd_row);
    for(i=0;i<refl_bd_row;++i){
        *(refl_coefX_cache+i)=(FLOAT *)malloc(refl_bd_col*sizeof(FLOAT));
        *(refl_coefY_cache+i)=(FLOAT *)malloc(refl_bd_col*sizeof(FLOAT)); 
        *(refl_coefZ_cache+i)=(FLOAT *)malloc(refl_bd_col*sizeof(FLOAT));
    }

    buildMat_reflective_boundary_conditions(PN, refl_coefX_cache, 
            refl_coefY_cache, refl_coefZ_cache);


    /* 对 coefXm 申请内存，并建立 coefXm，并注意在这套程序中
     * 边界面 X- 是用 1 来做标识的. */
    fg_coefXm=0;

    if(fg_coefXm==0){
        for(i=0;i<len_fixed;++i){
            if(1==fixed_bd_num[i]){
                fg_coefXm++;

                coefXm=(FLOAT **)malloc(sizeof(FLOAT *)*fixed_bd_row);
                for(j=0;j<fixed_bd_row;++j)
                    *(coefXm+j)=(FLOAT *)malloc(fixed_bd_col*sizeof(FLOAT));

                buildMat_fixed_boundary_conditions(PN, 1, Gauss_order, Gauss_points_l, 
                        Gauss_weights_l, Gauss_points_r, Gauss_weights_r, coefXm);
            }//endof_if(1==fixed_bd_num[i])
        }//endof_for(i=0;i<len_fixed;++i)


        if(fg_coefXm!=0){
            /* 因为要重复使用 matAxmatB 这块内存，所以用 calloc 函数，
             * calloc 与 malloc 的区别就是，calloc 申请内存后会自动
             * 初始化，而 malloc 申请内存后要用 memset 函数初始化。
             * 其实一般 malloc 申请内存后是要初始化的，对于要重复使用
             * 的内存一定要初始化。*/
            matAxmatB=(FLOAT **)calloc(nY,sizeof(FLOAT *));
            for(j=0;j<nY;++j)
                *(matAxmatB+j)=(FLOAT *)calloc(nY,sizeof(FLOAT));

            MatA_x_MatB(TRUE,fixed_bd_row,fixed_bd_col,coefXm,
                    fixed_bd_row,fixed_bd_col,coefXm,matAxmatB);

            arrangeMatrixInRows(nY,nY,matAxmatB,coefD_Xm_bd);

            for(j=0;j<fixed_bd_row;++j){
                free(coefXm[j]);
                coefXm[j]=NULL;
            }
            free(coefXm);
            coefXm=NULL;

             for(j=0;j<nY;++j){
                free(matAxmatB[j]);
                matAxmatB[j]=NULL;
             }
             free(matAxmatB);
             matAxmatB=NULL;
        }//endof_if(fg_coefXm!=0)
    }//endof_if(fg_coefXm==0)

    if(fg_coefXm==0){
        for(i=0;i<len_refl;++i){
            if(1==reflected_bd_num[i]){
                fg_coefXm++;

                coefXm=refl_coefX_cache;
            }//endof_if(1==reflected_bd_num[i])
        }//endof_for(i=0;i<len_refl;++i)

        if(fg_coefXm!=0){
            matAxmatB=(FLOAT **)calloc(nY,sizeof(FLOAT *));
            for(j=0;j<nY;++j)
                *(matAxmatB+j)=(FLOAT *)calloc(nY,sizeof(FLOAT));

            MatA_x_MatB(TRUE,refl_bd_row,refl_bd_col,coefXm,
                    refl_bd_row,refl_bd_col,coefXm,matAxmatB);

            arrangeMatrixInRows(nY,nY,matAxmatB,coefD_Xm_bd);

             for(j=0;j<nY;++j){
                free(matAxmatB[j]);
                matAxmatB[j]=NULL;
             }
             free(matAxmatB);
             matAxmatB=NULL;
        }//endof_if(fg_coefXm!=0)
    }//endof_if(fg_coefXm==0)

    if(fg_coefXm==0)
        coefD_Xm_bd=NULL;   


    /* 对 coefXp 申请内存，并建立 coefXp，并注意在这套程序中
     * 边界面 X+ 是用 2 来做标识的. */
    fg_coefXp=0;

    if(fg_coefXp==0){
        for(i=0;i<len_fixed;++i){
            if(2==fixed_bd_num[i]){
                fg_coefXp++;

                coefXp=(FLOAT **)malloc(sizeof(FLOAT *)*fixed_bd_row);
                for(j=0;j<fixed_bd_row;++j)
                    *(coefXp+j)=(FLOAT *)malloc(fixed_bd_col*sizeof(FLOAT));

                buildMat_fixed_boundary_conditions(PN, 2, Gauss_order, Gauss_points_l, 
                        Gauss_weights_l, Gauss_points_r, Gauss_weights_r, coefXp);
            }//endof_if(2==fixed_bd_num[i])
        }//endof_for(i=0;i<len_fixed;++i)

        if(fg_coefXp!=0){
            matAxmatB=(FLOAT **)calloc(nY,sizeof(FLOAT *));
            for(j=0;j<nY;++j)
                *(matAxmatB+j)=(FLOAT *)calloc(nY,sizeof(FLOAT));

            MatA_x_MatB(TRUE,fixed_bd_row,fixed_bd_col,coefXp,
                    fixed_bd_row,fixed_bd_col,coefXp,matAxmatB);

            arrangeMatrixInRows(nY,nY,matAxmatB,coefD_Xp_bd);

            for(j=0;j<fixed_bd_row;++j){
                free(coefXp[j]);
                coefXp[j]=NULL;
            }
            free(coefXp);
            coefXp=NULL;

             for(j=0;j<nY;++j){
                free(matAxmatB[j]);
                matAxmatB[j]=NULL;
             }
             free(matAxmatB);
             matAxmatB=NULL;
        }//endof_if(fg_coefXp!=0)
    }//endof_if(fg_coefXp==0)


    if(fg_coefXp==0){
        for(i=0;i<len_refl;++i){
            if(2==reflected_bd_num[i]){
                fg_coefXp++;

                coefXp=refl_coefX_cache;
            }//endof_if(2==reflected_bd_num[i])
        }//endof_for(i=0;i<len_refl;++i)

        if(fg_coefXp!=0){
            matAxmatB=(FLOAT **)calloc(nY,sizeof(FLOAT *));
            for(j=0;j<nY;++j)
                *(matAxmatB+j)=(FLOAT *)calloc(nY,sizeof(FLOAT));

            MatA_x_MatB(TRUE,refl_bd_row,refl_bd_col,coefXp,
                    refl_bd_row,refl_bd_col,coefXp,matAxmatB);

            arrangeMatrixInRows(nY,nY,matAxmatB,coefD_Xp_bd);

            for(j=0;j<nY;++j){
                free(matAxmatB[j]);
                matAxmatB[j]=NULL;
            }
            free(matAxmatB);
            matAxmatB=NULL;
        }//endof_if(fg_coefXp!=0)
    }//endof_if(fg_coefXp==0)

    if(fg_coefXp==0)
        coefD_Xp_bd=NULL;


    /* 对 coefYm 申请内存，并建立 coefYm，并注意在这套程序中
     * 边界面 Y- 是用 3 来做标识的. */
    fg_coefYm=0;

    if(fg_coefYm==0){
        for(i=0;i<len_fixed;++i){
            if(3==fixed_bd_num[i]){
                fg_coefYm++;

                coefYm=(FLOAT **)malloc(sizeof(FLOAT *)*fixed_bd_row);
                for(j=0;j<fixed_bd_row;++j)
                    *(coefYm+j)=(FLOAT *)malloc(fixed_bd_col*sizeof(FLOAT));

                buildMat_fixed_boundary_conditions(PN, 3, Gauss_order, Gauss_points_l, 
                        Gauss_weights_l, Gauss_points_r, Gauss_weights_r, coefYm);
            }//endof_if(3==fixed_bd_num[i])
        }//endof_for(i=0;i<len_fixed;++i)

        if(fg_coefYm!=0){
            matAxmatB=(FLOAT **)calloc(nY,sizeof(FLOAT *));
            for(j=0;j<nY;++j)
                *(matAxmatB+j)=(FLOAT *)calloc(nY,sizeof(FLOAT));

            MatA_x_MatB(TRUE,fixed_bd_row,fixed_bd_col,coefYm,
                    fixed_bd_row,fixed_bd_col,coefYm,matAxmatB);

            arrangeMatrixInRows(nY,nY,matAxmatB,coefD_Ym_bd);

            for(j=0;j<fixed_bd_row;++j){
                free(coefYm[j]);
                coefYm[j]=NULL;
            }
            free(coefYm);
            coefYm=NULL;

             for(j=0;j<nY;++j){
                free(matAxmatB[j]);
                matAxmatB[j]=NULL;
             }
             free(matAxmatB);
             matAxmatB=NULL;
        }//endof_if(fg_coefYm!=0)
    }//endof_if(fg_coefYm==0)

    if(fg_coefYm==0){
        for(i=0;i<len_refl;++i){
            if(3==reflected_bd_num[i]){
                fg_coefYm++;

                coefYm=refl_coefY_cache;
            }//endof_if(3==reflected_bd_num[i])
        }//endof_for(i=0;i<len_refl;++i)

        if(fg_coefYm!=0){
            matAxmatB=(FLOAT **)calloc(nY,sizeof(FLOAT *));
            for(j=0;j<nY;++j)
                *(matAxmatB+j)=(FLOAT *)calloc(nY,sizeof(FLOAT));

            MatA_x_MatB(TRUE,refl_bd_row,refl_bd_col,coefYm,
                    refl_bd_row,refl_bd_col,coefYm,matAxmatB);

            arrangeMatrixInRows(nY,nY,matAxmatB,coefD_Ym_bd);

             for(j=0;j<nY;++j){
                free(matAxmatB[j]);
                matAxmatB[j]=NULL;
             }
             free(matAxmatB);
             matAxmatB=NULL;
        }//endof_if(fg_coefYm!=0)
    }//endof_if(fg_coefYm==0)

    if(fg_coefYm==0)
        coefD_Ym_bd=NULL;   


    /* 对 coefYp 申请内存，并建立 coefYp，并注意在这套程序中
     * 边界面 Y+ 是用 4 来做标识的. */
    fg_coefYp=0;

    if(fg_coefYp==0){
        for(i=0;i<len_fixed;++i){
            if(4==fixed_bd_num[i]){
                fg_coefYp++;

                coefYp=(FLOAT **)malloc(sizeof(FLOAT *)*fixed_bd_row);
                for(j=0;j<fixed_bd_row;++j)
                    *(coefYp+j)=(FLOAT *)malloc(fixed_bd_col*sizeof(FLOAT));

                buildMat_fixed_boundary_conditions(PN, 4, Gauss_order, Gauss_points_l, 
                        Gauss_weights_l, Gauss_points_r, Gauss_weights_r, coefYp);
            }//endof_if(4==fixed_bd_num[i])
        }//endof_for(i=0;i<len_fixed;++i)

        if(fg_coefYp!=0){
            matAxmatB=(FLOAT **)calloc(nY,sizeof(FLOAT *));
            for(j=0;j<nY;++j)
                *(matAxmatB+j)=(FLOAT *)calloc(nY,sizeof(FLOAT));

            MatA_x_MatB(TRUE,fixed_bd_row,fixed_bd_col,coefYp,
                    fixed_bd_row,fixed_bd_col,coefYp,matAxmatB);

            arrangeMatrixInRows(nY,nY,matAxmatB,coefD_Yp_bd);

            for(j=0;j<fixed_bd_row;++j){
                free(coefYp[j]);
                coefYp[j]=NULL;
            }
            free(coefYp);
            coefYp=NULL;

             for(j=0;j<nY;++j){
                free(matAxmatB[j]);
                matAxmatB[j]=NULL;
             }
             free(matAxmatB);
             matAxmatB=NULL;
        }//endof_if(fg_coefYp!=0)
    }//endof_if(fg_coefYp==0)

    if(fg_coefYp==0){
        for(i=0;i<len_refl;++i){
            if(4==reflected_bd_num[i]){
                fg_coefYp++;

                coefYp=refl_coefY_cache;
            }//endof_if(4==reflected_bd_num[i])
        }//endof_for(i=0;i<len_refl;++i)

        if(fg_coefYp!=0){
            matAxmatB=(FLOAT **)calloc(nY,sizeof(FLOAT *));
            for(j=0;j<nY;++j)
                *(matAxmatB+j)=(FLOAT *)calloc(nY,sizeof(FLOAT));

            MatA_x_MatB(TRUE,refl_bd_row,refl_bd_col,coefYp,
                    refl_bd_row,refl_bd_col,coefYp,matAxmatB);

            arrangeMatrixInRows(nY,nY,matAxmatB,coefD_Yp_bd);

             for(j=0;j<nY;++j){
                free(matAxmatB[j]);
                matAxmatB[j]=NULL;
             }
             free(matAxmatB);
             matAxmatB=NULL;
        }//endof_if(fg_coefYp!=0)
    }//endof_if(fg_coefYp==0)

    if(fg_coefYp==0)
        coefD_Yp_bd=NULL;



    /* 对 coefZm 申请内存，并建立 coefZm，并注意在这套程序中
     * 边界面 Z- 是用 5 来做标识的. */
    fg_coefZm=0;

    if(fg_coefZm==0){
        for(i=0;i<len_fixed;++i){
            if(5==fixed_bd_num[i]){
                fg_coefZm++;

                coefZm=(FLOAT **)malloc(sizeof(FLOAT *)*fixed_bd_row);
                for(j=0;j<fixed_bd_row;++j)
                    *(coefZm+j)=(FLOAT *)malloc(fixed_bd_col*sizeof(FLOAT));

                buildMat_fixed_boundary_conditions(PN, 5, Gauss_order, Gauss_points_l, 
                        Gauss_weights_l, Gauss_points_r, Gauss_weights_r, coefZm);
            }//endof_if(5==fixed_bd_num[i])
        }//endof_for(i=0;i<len_fixed;++i)

        if(fg_coefZm!=0){
            matAxmatB=(FLOAT **)calloc(nY,sizeof(FLOAT *));
            for(j=0;j<nY;++j)
                *(matAxmatB+j)=(FLOAT *)calloc(nY,sizeof(FLOAT));

            MatA_x_MatB(TRUE,fixed_bd_row,fixed_bd_col,coefZm,
                    fixed_bd_row,fixed_bd_col,coefZm,matAxmatB);

            arrangeMatrixInRows(nY,nY,matAxmatB,coefD_Zm_bd);

            for(j=0;j<fixed_bd_row;++j){
                free(coefZm[j]);
                coefZm[j]=NULL;
            }
            free(coefZm);
            coefZm=NULL;

             for(j=0;j<nY;++j){
                free(matAxmatB[j]);
                matAxmatB[j]=NULL;
             }
             free(matAxmatB);
             matAxmatB=NULL;
        }//endof_if(fg_coefZm!=0)
    }//endof_if(fg_coefZm==0)

    if(fg_coefZm==0){
        for(i=0;i<len_refl;++i){
            if(5==reflected_bd_num[i]){
                fg_coefZm++;

                coefZm=refl_coefZ_cache;
            }//endof_if(5==reflected_bd_num[i])
        }//endof_for(i=0;i<len_refl;++i)

        if(fg_coefZm!=0){
            matAxmatB=(FLOAT **)calloc(nY,sizeof(FLOAT *));
            for(j=0;j<nY;++j)
                *(matAxmatB+j)=(FLOAT *)calloc(nY,sizeof(FLOAT));

            MatA_x_MatB(TRUE,refl_bd_row,refl_bd_col,coefZm,
                    refl_bd_row,refl_bd_col,coefZm,matAxmatB);

            arrangeMatrixInRows(nY,nY,matAxmatB,coefD_Zm_bd);

             for(j=0;j<nY;++j){
                free(matAxmatB[j]);
                matAxmatB=NULL;
             }
             free(matAxmatB);
             matAxmatB[j]=NULL;
        }//endof_if(fg_coefZm!=0)
    }//endof_if(fg_coefZm==0)

    if(fg_coefZm==0)
        coefD_Zm_bd=NULL;   


    /* 对 coefZp 申请内存，并建立 coefZp，并注意在这套程序中
     * 边界面 Z+ 是用 6 来做标识的. */
    fg_coefZp=0;

    if(fg_coefZp==0){
        for(i=0;i<len_fixed;++i){
            if(6==fixed_bd_num[i]){
                fg_coefZp++;

                coefZp=(FLOAT **)malloc(sizeof(FLOAT *)*fixed_bd_row);
                for(j=0;j<fixed_bd_row;++j)
                    *(coefZp+j)=(FLOAT *)malloc(fixed_bd_col*sizeof(FLOAT));

                buildMat_fixed_boundary_conditions(PN, 6, Gauss_order, Gauss_points_l, 
                        Gauss_weights_l, Gauss_points_r, Gauss_weights_r, coefZp);
            }//endof_if(6==fixed_bd_num[i])
        }//endof_for(i=0;i<len_fixed;++i)

        if(fg_coefZp!=0){
            matAxmatB=(FLOAT **)calloc(nY,sizeof(FLOAT *));
            for(j=0;j<nY;++j)
                *(matAxmatB+j)=(FLOAT *)calloc(nY,sizeof(FLOAT));

            MatA_x_MatB(TRUE,fixed_bd_row,fixed_bd_col,coefZp,
                    fixed_bd_row,fixed_bd_col,coefZp,matAxmatB);

            arrangeMatrixInRows(nY,nY,matAxmatB,coefD_Zp_bd);

            for(j=0;j<fixed_bd_row;++j){
                free(coefZp[j]);
                coefZp[j]=NULL;
            }
            free(coefZp);
            coefZp=NULL;

             for(j=0;j<nY;++j){
                free(matAxmatB[j]);
                matAxmatB[j]=NULL;
             }
             free(matAxmatB);
             matAxmatB=NULL;
        }//endof_if(fg_coefZp!=0)
    }//endof_if(fg_coefZp==0)

    if(fg_coefZp==0){
        for(i=0;i<len_refl;++i){
            if(6==reflected_bd_num[i]){
                fg_coefZp++;

                coefZp=refl_coefZ_cache;
            }//endof_if(6==reflected_bd_num[i])
        }//endof_for(i=0;i<len_refl;++i)

        if(fg_coefZp!=0){
            matAxmatB=(FLOAT **)calloc(nY,sizeof(FLOAT *));
            for(j=0;j<nY;++j)
                *(matAxmatB+j)=(FLOAT *)calloc(nY,sizeof(FLOAT));

            MatA_x_MatB(TRUE,refl_bd_row,refl_bd_col,coefZp,
                    refl_bd_row,refl_bd_col,coefZp,matAxmatB);

            arrangeMatrixInRows(nY,nY,matAxmatB,coefD_Zp_bd);

             for(j=0;j<nY;++j){
                free(matAxmatB[j]);
                matAxmatB[j]=NULL;
             }
             free(matAxmatB);
             matAxmatB=NULL;
        }//endof_if(fg_coefZp!=0)
    }//endof_if(fg_coefZp==0)

    if(fg_coefZp==0)
        coefD_Zp_bd=NULL;

    for(j=0;j<refl_bd_row;++j){
        free(refl_coefX_cache[j]);
        free(refl_coefY_cache[j]);
        free(refl_coefZ_cache[j]);
        refl_coefX_cache[j]=NULL;
        refl_coefY_cache[j]=NULL;
        refl_coefZ_cache[j]=NULL;
    }
    free(refl_coefX_cache);
    free(refl_coefY_cache);
    free(refl_coefZ_cache);
    refl_coefX_cache=NULL;
    refl_coefY_cache=NULL;
    refl_coefZ_cache=NULL;

}//endof_build_coefD_xx_bd_()

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

