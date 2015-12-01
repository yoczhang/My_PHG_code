/*
 * This test2.c file is aimed to test the block matrix defined 
 * in PHG, the problems' specific description will be found in 
 * the email sended to LiuHui "关于分块矩阵".
 * 
 * To see more details, look at "草稿本1", page 8 and page 9.
 *
 * This code doesn't complete, LiuHui has denied my thought, 
 * he says PHG may not support the double level block matrix.
 *
 */

#include "phg.h"

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
# undef ELEMENT
typedef SIMPLEX ELEMENT;
typedef USE_MPI 1;
#endif
 
#include <string.h>
#include <math.h>

int 
main(int argc, char *argv[])
{
    
    char *fn = "../mytest.mesh";
    GRID *g;
    
    phgInit(&argc, &argv);

    g = phgNewGrid(-1);

    if (!phgImport(g, fn, TRUE))
	phgError(1, "can't read file \"%s\".\n", fn);

    /*-----------------------------------------------------------------*/
    /*--------------- the following is my test ------------------------*/
    ELEMENT *e;
    DOF *u;
    DOF *X;
    MAT *A11, *A12, *A21, *A22, *B11, *B12, *B21, *B22, *A, *B;
    MAT pmatA[4], pmatB[4], pmatM[4];
    FLOAT coefA[4]={1.0,1.0,1.0,1.0};
    FLOAT coefB[4]={1.0,1.0,1.0,1.0};
    FLOAT coefM[4]={1.0,-1.0,1.0,1.0};

    SOLVER *solver=NULL;

    u=phgDofNew(g, DOF_DEFAULT, 1, "u", DofInterpolation);
    X=phgDofNew(g, DOF_DEFAULT, 4, "X", DofInterpolation);



}//endof_main()
