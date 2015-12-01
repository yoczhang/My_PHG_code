/*
 * This test1.c file is aimed to make the difference of 
 * the PHG functions phgQuadDofTimesBas() and phgQuadDofDotBas().
 *
 * =_=! I think the phgQuadDofDotBas() is used the ND1(ND2...) basises 
 * in maxwell problems, grnerally, we just use the phgQuadDofTimesBas() 
 * to compute the quad. 
 * More deatils we can see the papery code "maxwell-comples.c" I have printed.
 *
 * Be very careful the important codes "phgOptionsPreset("-dof_type ND1");"
 * in the following main() code.
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
  
  
    phgOptionsPreset("-dof_type ND1");//This is a very important option, seting DOF_DEFAULT is ND1.

    phgInit(&argc, &argv);

    g = phgNewGrid(-1);

    if (!phgImport(g, fn, TRUE))
	phgError(1, "can't read file \"%s\".\n", fn);


    /*-----------------------------------------------------------------*/
    /*--------------- the following is my test ------------------------*/
    ELEMENT *e;
    DOF *u_re, *u1, *u_h;
    SOLVER *pc=NULL;
    int i;

    u_re = phgDofNew(g, DOF_DEFAULT, 1, "u_re", DofInterpolation);
    u_re->DB_mask = BDRY_MASK;
       
    u1=phgDofNew(g,DOF_CONSTANT,3,"u1",DofNoAction);
    phgDofSetDataByValuesV(u1,(FLOAT)1.0,(FLOAT)1.0,(FLOAT)1.0);

    pc = phgSolverCreate(SOLVER_PCG, u_re, NULL);
    u_h = pc->mat->rmap->dofs[0];
    int N=u_h->type->nbas;
    FLOAT C1[N];
    
    int nvalues1;
    nvalues1=DofDim(u1);
    printf("nvalues1(DofDim(u1))=%d \n",nvalues1);

    int nvalues2;
    nvalues2=DofTypeDim(u_h);
    printf("nvalues2(DofTypeDim(u_h))=%d \n",nvalues2);

    ForAllElements(g,e){
        for(i=0;i<N;++i){
            C1[i]=phgQuadDofDotBas(e,u1,u_h,i,QUAD_DEFAULT);
        }
    }//endof_ForAllElements()


    return 0;

}
