/*
 * This test4.c is aimed to test the mesh with 
 * boundary markers.
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
    
    char *fn = "../mytest_bmarker.mesh";
    GRID *g;
    
    phgInit(&argc, &argv);

    g = phgNewGrid(-1);

    if (!phgImport(g, fn, TRUE))
	phgError(1, "can't read file \"%s\".\n", fn);

    /*-----------------------------------------------------------------*/
    /*--------------- the following is my test ------------------------*/
    
    INT i;
    int n=0;

    ELEMENT *e;

    ForAllElements(g,e){
    
        printf("\n/*****************************************************/\n");
        printf("IN %dth ELEMENT(BEGIN FROM THE 0th ELEMENT)\n",n);

        // there has 4 facets in one element
        for(i=0;i<4;i++){
            /*
            printf("in the %dth element, the %dth face is the usr's boundary %d.\n",
                        n,i,e->bound_type[i]);
            */
            printf("the %dth face's valuve is faces[%d]=%d\n",i,i,e->faces[i]);

            if(e->bound_type[i]==BDRY_USER0)
                printf("in the %dth element, the %dth face is the usr's boundary %d.\n",
                        n,i,e->bound_type[i]);

            if(e->bound_type[i]==BDRY_USER1)
                printf("in the %dth element, the %dth face is the usr's boundary %d.\n",
                        n,i,e->bound_type[i]);

            if(e->bound_type[i]==BDRY_USER2)
                printf("in the %dth element, the %dth face is the usr's boundary %d.\n",
                        n,i,e->bound_type[i]);

            if(e->bound_type[i]==BDRY_USER3)
                printf("in the %dth element, the %dth face is the usr's boundary %d.\n",
                        n,i,e->bound_type[i]);
        }

        n++;
        printf("/*****************************************************/\n");
    }//endof_ForAllElements

    phgFreeGrid(&g);
    phgFinalize();

    return 0;

}//endof_main()
