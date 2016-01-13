/*
 * The fix_boundary_conditions.h is aimed to deal the 
 * boundary conditions 
 */

# include "phg.h"
#include "other_functions.h"
# include "int_associated_legendre_polyns.h"
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>
#include <malloc.h>

/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/*
 * To compute the coefficient D(m)
 */
FLOAT 
D(int m);
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/*
 * To compute Ip(l, m, l1, m1)
 */
FLOAT
Ip(int l, int m, int l1, int m1, int Gauss_order, FLOAT *Gauss_points, 
        FLOAT *Gauss_weights);
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/




/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
FLOAT
I(int l, int m, int l1, int m1, int lable, int Gauss_order, FLOAT *Gauss_points_l, 
        FLOAT *Gauss_weights_l, FLOAT *Gauss_points_r, FLOAT *Gauss_weights_r);
/*
 * int lable standsfor the different faces.
 * lable=1,2 standsfor the faces x-,x+
 * lable=3,4 standsfor the faces y-,y+
 * lable=5,6 standsfor the faces z-,z+
 */
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
void
fixed_boundary_conditions(int N, int lable, int Gauss_order, FLOAT *Gauss_points_l, 
        FLOAT *Gauss_weights_l, FLOAT *Gauss_points_r, FLOAT *Gauss_weights_r, FLOAT **C_1);
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
 * Output: C_1, the matirx.
 */
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/




/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
void
buildMat_reflective_boundary_conditions(int N, FLOAT **C_2X, FLOAT **C_2Y, 
        FLOAT **C_2Z);
/*
 * This function is aimed to build the reflective boundary conditions' matrixes.
 *
 * Input: N, the PN approximating number, is odd.
 *
 * Output: C_2X, the matrix based on the boundary faces x- or x+.
 *         C_2Y, the matrix based on the boundary faces y- or y+.
 *         C_2Z, the matrix based on the boundary faces z- or z+.
 */
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/




/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
void
build_coefD_xx_bd(int PN, int Gauss_order, FLOAT *Gauss_points_l, FLOAT *Gauss_weights_l, 
        FLOAT *Gauss_points_r, FLOAT *Gauss_weights_r, 
        FLOAT *coefD_Xm_bd, FLOAT *coefD_Xp_bd, FLOAT *coefD_Ym_bd, 
        FLOAT *coefD_Yp_bd, FLOAT *coefD_Zm_bd, FLOAT *coefD_Zp_bd);
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

/**********************************************************************************/
/*--------------------------------------------------------------------------------*/





/*--------------------------------------------------------------------------------*/
/**********************************************************************************/


/**********************************************************************************/
/*--------------------------------------------------------------------------------*/

