/*
 * The int_associated_legendre_ployns.c is aimed to computed integration 
 * of associated legendre ploynoimals using the Gauss-Legender quad rules.
 *
 */

# include "phg.h"
# include "int_associated_legendre_polyns.h"
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <string.h>

/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
INT 
factorial(n)
{
    INT v;
    v=1;
    while(n>0){
        v*=n--;
    }
    return v;

}
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
void
get_Gauss_points_weights(INT left, INT right, INT Gauss_order, FLOAT *Gauss_points, 
        FLOAT *Gauss_weights)
{
    INT kind;
    FLOAT alpha;
    FLOAT beta;

    /*
     * 下面这3个具体的赋值，其实是要根据 legendre_rule.c 中
     * 的 main() 函数来的，具体信息可以查看 legendre_rule.c
     */
    kind=1;
    alpha=0.0;
    beta=0.0;

    cgqf(Gauss_order, kind, alpha, beta, left, right, Gauss_points, Gauss_weights);

}
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/*
 * To compute \int_0^1{P_l^m(x) * P_l^m(x)}
 */
FLOAT
int_associated_legendre_polyns(int l, int m, int l1, int m1, int mm, FLOAT *x, 
        FLOAT *wts)
/* mm 其实为所计算的点的个数，其他函数调用 int_associated_legendre_polyns() 
 * 这个函数时的 mm 取值时为高斯积分点个数
 */

{
    int i;
    //int n_fac=factorial(l+m);
    FLOAT v;
    FLOAT *val_pm_polyn;
    FLOAT *val_pm_polyn1;
    FLOAT C_lm;
    FLOAT C_l1m1;

    val_pm_polyn=(FLOAT *)malloc(mm * sizeof(FLOAT));
    val_pm_polyn1=(FLOAT *)malloc(mm * sizeof(FLOAT));

    /*
     * 因为自己需要求 P_l^{|m|}(x) 的值,
     * 细节部分"程序过程整理 1"Page:40
     * 并且这里 pm_polynomial_value(mm, l, m, x) 函数返回的值
     * 是 P_0^m(x<mm>), P_2^m(x{mm}), P_3^m(x{mm})...P_l^m(x{mm}).
     * 要注意到每个 P_l^m(x{mm}) 含义都是包含 mm 个值，即
     * P_l^m(x[0]), P_l^m(x[1])... P_l^m(x[mm-1]).
     * 所以为了求得 P_l^m 的 mm 个(高斯积分点数) 的值，在后面要 +l*mm.
     */    
    m=abs(m);
    m1=abs(m1);
    val_pm_polyn=pm_polynomial_value(mm, l, m, x)+l*mm;
    val_pm_polyn1=pm_polynomial_value(mm, l1, m1, x)+l1*mm;

    //for(i=0;i<mm;i++)
        //printf("val_pm_polyn[%d]=%f \n",i,*(val_pm_polyn+i));

    C_lm=sqrt( ((FLOAT)(2*l+1)*factorial(l-m))/factorial(l+m) );
    C_l1m1=sqrt( ((FLOAT)(2*l1+1)*factorial(l1-m1))/factorial(l1+m1) );
    //printf("C_lm=%f\n",C_lm);

    v=0.0;
    for(i=0;i<mm;i++)
        v+=*(wts+i) * *(val_pm_polyn+i) * *(val_pm_polyn1+i);

    /*
     * 因为 pm_polynomial_value() 所计算的 P_l^m 里面包含了
     * 乘以(-1)^m(这一点已经确认),而自己在离散时所定义的 P_l^m 里面
     * 是没有乘以(-1)^m 的，所以，下面在(|m|+|m1|)%2!=0时会有一个转换。
     */
    if((m+m1)%2!=0)
        v=-v;

    v*=C_lm*C_l1m1;

    //printf("%d!=%d \n",l+m,n_fac);
    return v;
}//endof_int_associated_legendre_ployns
/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/
/*
 * Generating the associated_legendre_ployns values in the points x[], the details 
 * in the following interpretations.
 */
FLOAT *
pm_polynomial_value(int mm, int n, int m, FLOAT *x)//fllowing has my Debug.

/******************************************************************************/
/*
  Purpose:

    PM_POLYNOMIAL_VALUE evaluates the Legendre polynomials Pm(n,m,x).

  Differential equation:

    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0

  First terms:

    M = 0  ( = Legendre polynomials of first kind P(N,X) )

    Pm(0,0,x) =    1
    Pm(1,0,x) =    1 X
    Pm(2,0,x) = (  3 X^2 -   1)/2
    Pm(3,0,x) = (  5 X^3 -   3 X)/2
    Pm(4,0,x) = ( 35 X^4 -  30 X^2 +   3)/8
    Pm(5,0,x) = ( 63 X^5 -  70 X^3 +  15 X)/8
    Pm(6,0,x) = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
    Pm(7,0,x) = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16

    M = 1

    Pm(0,1,x) =   0
    Pm(1,1,x) =   1 * SQRT(1-X^2)
    Pm(2,1,x) =   3 * SQRT(1-X^2) * X
    Pm(3,1,x) = 1.5 * SQRT(1-X^2) * (5*X^2-1)
    Pm(4,1,x) = 2.5 * SQRT(1-X^2) * (7*X^3-3*X)

    M = 2

    Pm(0,2,x) =   0
    Pm(1,2,x) =   0
    Pm(2,2,x) =   3 * (1-X^2)
    Pm(3,2,x) =  15 * (1-X^2) * X
    Pm(4,2,x) = 7.5 * (1-X^2) * (7*X^2-1)

    M = 3

    Pm(0,3,x) =   0
    Pm(1,3,x) =   0
    Pm(2,3,x) =   0
    Pm(3,3,x) =  15 * (1-X^2)^1.5
    Pm(4,3,x) = 105 * (1-X^2)^1.5 * X

    M = 4

    Pm(0,4,x) =   0
    Pm(1,4,x) =   0
    Pm(2,4,x) =   0
    Pm(3,4,x) =   0
    Pm(4,4,x) = 105 * (1-X^2)^2

  Recursion:

    if N < M:
      Pm(N,M,x) = 0
    if N = M:
      Pm(N,M,x) = (2*M-1)!! * (1-X*X)^(M/2) where N!! means the product of
      all the odd integers less than or equal to N.
    if N = M+1:
      Pm(N,M,x) = X*(2*M+1)*Pm(M,M,x)
    if M+1 < N:
      Pm(N,M,x) = ( X*(2*N-1)*Pm(N-1,M,x) - (N+M-1)*Pm(N-2,M,x) )/(N-M)

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Reference:

    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.

  Parameters:

    Input, int MM, the number of evaluation points.

    Input, int N, the maximum first index of the Legendre
    function, which must be at least 0.

    Input, int M, the second index of the Legendre function,
    which must be at least 0, and no greater than N.

    Input, double X[MM], the point at which the function is to be
    evaluated.

    Output, double PM_POLYNOMIAL_VALUE[MM*(N+1)], the function values.
    output, P_0^m(x[mm]), P_1^m(x[mm]), P_2^m(x[mm]), ..., P_n^m(x[mm])
*/
{
  double fact;
  int i;
  int j;
  int k;
  double *v;

  v = ( double * ) malloc ( mm*(n+1) * sizeof ( double ) );

  // myDebug
  // the 'j' standsfor the Associated Legender ploynomial "P_n^m" subindex 'n'
  // I have confirmed.  
  for ( j = 0; j < n + 1; j++ )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = 0.0;
    }
  }
/*
  J = M is the first nonzero function.
*/
  if ( m <= n )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+m*mm] = 1.0;
    }

    fact = 1.0;
    for ( k = 0; k < m; k++ )
    {
      for ( i = 0; i < mm; i++ )
      {
        v[i+m*mm] = - v[i+m*mm] * fact * sqrt ( 1.0 - x[i] * x[i] );
      }
      fact = fact + 2.0;
    }
  }
/*
  J = M + 1 is the second nonzero function.
*/
  if ( m + 1 <= n )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+(m+1)*mm] = x[i] * ( double ) ( 2 * m + 1 ) * v[i+m*mm];
    }
  }
/*
  Now we use a three term recurrence.
*/
  for ( j = m + 2; j <= n; j++ )
  {
    for ( i = 0; i < mm; i++ )
    {
      v[i+j*mm] = ( ( double ) ( 2 * j     - 1 ) * x[i] * v[i+(j-1)*mm]
                  + ( double ) (   - j - m + 1 ) *        v[i+(j-2)*mm] )
                  / ( double ) (     j - m     );
    }
  }

  //return v;
  /*
  printf("----------------- in pm_polynomial_value() ------------------------\n");
    for(i=0;i<mm*(n+1);i++){
        printf("v[%d]=%f\n",i,*(v+i));
    }
    printf("----------------- in pm_polynomial_value() ------------------------\n");
    */
    return v;

}


/**********************************************************************************/
/*--------------------------------------------------------------------------------*/
void cdgqf ( int nt, int kind, double alpha, double beta, double t[], 
  double wts[] )

/******************************************************************************/
/*
  Purpose:

    CDGQF computes a Gauss quadrature formula with default A, B and simple knots.

  Discussion:

    This routine computes all the knots and weights of a Gauss quadrature
    formula with a classical weight function with default values for A and B,
    and only simple knots.

    There are no moments checks and no printing is done.

    Use routine EIQFS to evaluate a quadrature computed by CGQFS.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

  Parameters:

    Input, int NT, the number of knots.

    Input, int KIND, the rule.
    1, Legendre,             (a,b)       1.0
    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta

    Input, double ALPHA, the value of Alpha, if needed.

    Input, double BETA, the value of Beta, if needed.

    Output, double T[NT], the knots.

    Output, double WTS[NT], the weights.
*/
{
  double *aj;
  double *bj;
  double zemu;

  parchk ( kind, 2 * nt, alpha, beta );
/*
  Get the Jacobi matrix and zero-th moment.
*/
  aj = ( double * ) malloc ( nt * sizeof ( double ) );
  bj = ( double * ) malloc ( nt * sizeof ( double ) );

  zemu = class_matrix ( kind, nt, alpha, beta, aj, bj );
/*
  Compute the knots and weights.
*/
  sgqf ( nt, aj, bj, zemu, t, wts );

  free ( aj );
  free ( bj );

  return;
}
/******************************************************************************/

void cgqf ( int nt, int kind, double alpha, double beta, double a, double b, 
  double t[], double wts[] )

/******************************************************************************/
/*
  Purpose:

    CGQF computes knots and weights of a Gauss quadrature formula.

  Discussion:

    The user may specify the interval (A,B).

    Only simple knots are produced.

    Use routine EIQFS to evaluate this quadrature formula.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 September 2013

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

  Parameters:

    Input, int NT, the number of knots.

    Input, int KIND, the rule.
    1, Legendre,             (a,b)       1.0
    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^-0.5)
    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta

    Input, double ALPHA, the value of Alpha, if needed.

    Input, double BETA, the value of Beta, if needed.

    Input, double A, B, the interval endpoints.

    Output, double T[NT], the knots.

    Output, double WTS[NT], the weights.
*/
{
  int i;
  int *mlt;
  int *ndx;
/*
  Compute the Gauss quadrature formula for default values of A and B.
*/
  cdgqf ( nt, kind, alpha, beta, t, wts );
/*
  Prepare to scale the quadrature formula to other weight function with 
  valid A and B.
*/
  mlt = ( int * ) malloc ( nt * sizeof ( int ) );
  for ( i = 0; i < nt; i++ )
  {
    mlt[i] = 1;
  }
  ndx = ( int * ) malloc ( nt * sizeof ( int ) );
  for ( i = 0; i < nt; i++ )
  {
    ndx[i] = i + 1;
  }
  scqf ( nt, t, mlt, wts, nt, ndx, wts, t, kind, alpha, beta, a, b );

  free ( mlt );
  free ( ndx );

  return;
}
/******************************************************************************/

double class_matrix ( int kind, int m, double alpha, double beta, double aj[], 
  double bj[] )

/******************************************************************************/
/*
  Purpose:

    CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.

  Discussion:

    This routine computes the diagonal AJ and sub-diagonal BJ
    elements of the order M tridiagonal symmetric Jacobi matrix
    associated with the polynomials orthogonal with respect to
    the weight function specified by KIND.

    For weight functions 1-7, M elements are defined in BJ even
    though only M-1 are needed.  For weight function 8, BJ(M) is
    set to zero.

    The zero-th moment of the weight function is returned in ZEMU.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

  Parameters:

    Input, int KIND, the rule.
    1, Legendre,             (a,b)       1.0
    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta

    Input, int M, the order of the Jacobi matrix.

    Input, double ALPHA, the value of Alpha, if needed.

    Input, double BETA, the value of Beta, if needed.

    Output, double AJ[M], BJ[M], the diagonal and subdiagonal
    of the Jacobi matrix.

    Output, double CLASS_MATRIX, the zero-th moment.
*/
{
  double a2b2;
  double ab;
  double aba;
  double abi;
  double abj;
  double abti;
  double apone;
  int i;
  const double pi = 3.14159265358979323846264338327950;
  double temp;
  double temp2;
  double zemu;

  temp = r8_epsilon ( );

  parchk ( kind, 2 * m - 1, alpha, beta );

  temp2 = 0.5;

  if ( 500.0 * temp < fabs ( pow ( r8_gamma ( temp2 ), 2 ) - pi ) )
  {
    printf ( "\n" );
    printf ( "CLASS_MATRIX - Fatal error!\n" );
    printf ( "  Gamma function does not match machine parameters.\n" );
    exit ( 1 );
  }

  if ( kind == 1 )
  {
    ab = 0.0;

    zemu = 2.0 / ( ab + 1.0 );

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    for ( i = 1; i <= m; i++ )
    {
      abi = i + ab * ( i % 2 );
      abj = 2 * i + ab;
      bj[i-1] = sqrt ( abi * abi / ( abj * abj - 1.0 ) );
    }
  }
  else if ( kind == 2 )
  {
    zemu = pi;

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    bj[0] =  sqrt ( 0.5 );
    for ( i = 1; i < m; i++ )
    {
      bj[i] = 0.5;
    }
  }
  else if ( kind == 3 )
  {
    ab = alpha * 2.0;
    zemu = pow ( 2.0, ab + 1.0 ) * pow ( r8_gamma ( alpha + 1.0 ), 2 )
      / r8_gamma ( ab + 2.0 );

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    bj[0] = sqrt ( 1.0 / ( 2.0 * alpha + 3.0 ) );
    for ( i = 2; i <= m; i++ )
    {
      bj[i-1] = sqrt ( i * ( i + ab ) / ( 4.0 * pow ( i + alpha, 2 ) - 1.0 ) );
    }
  }
  else if ( kind == 4 )
  {
    ab = alpha + beta;
    abi = 2.0 + ab;
    zemu = pow ( 2.0, ab + 1.0 ) * r8_gamma ( alpha + 1.0 ) 
      * r8_gamma ( beta + 1.0 ) / r8_gamma ( abi );
    aj[0] = ( beta - alpha ) / abi;
    bj[0] = sqrt ( 4.0 * ( 1.0 + alpha ) * ( 1.0 + beta ) 
      / ( ( abi + 1.0 ) * abi * abi ) );
    a2b2 = beta * beta - alpha * alpha;

    for ( i = 2; i <= m; i++ )
    {
      abi = 2.0 * i + ab;
      aj[i-1] = a2b2 / ( ( abi - 2.0 ) * abi );
      abi = abi * abi;
      bj[i-1] = sqrt ( 4.0 * i * ( i + alpha ) * ( i + beta ) * ( i + ab ) 
        / ( ( abi - 1.0 ) * abi ) );
    }
  }
  else if ( kind == 5 )
  {
    zemu = r8_gamma ( alpha + 1.0 );

    for ( i = 1; i <= m; i++ )
    {
      aj[i-1] = 2.0 * i - 1.0 + alpha;
      bj[i-1] = sqrt ( i * ( i + alpha ) );
    }
  }
  else if ( kind == 6 )
  {
    zemu = r8_gamma ( ( alpha + 1.0 ) / 2.0 );

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    for ( i = 1; i <= m; i++ )
    {
      bj[i-1] = sqrt ( ( i + alpha * ( i % 2 ) ) / 2.0 );
    }
  }
  else if ( kind == 7 )
  {
    ab = alpha;
    zemu = 2.0 / ( ab + 1.0 );

    for ( i = 0; i < m; i++ )
    {
      aj[i] = 0.0;
    }

    for ( i = 1; i <= m; i++ )
    {
      abi = i + ab * ( i % 2 );
      abj = 2 * i + ab;
      bj[i-1] = sqrt ( abi * abi / ( abj * abj - 1.0 ) );
    }
  }
  else if ( kind == 8 )
  {
    ab = alpha + beta;
    zemu = r8_gamma ( alpha + 1.0 ) * r8_gamma ( - ( ab + 1.0 ) ) 
      / r8_gamma ( - beta );
    apone = alpha + 1.0;
    aba = ab * apone;
    aj[0] = - apone / ( ab + 2.0 );
    bj[0] = - aj[0] * ( beta + 1.0 ) / ( ab + 2.0 ) / ( ab + 3.0 );
    for ( i = 2; i <= m; i++ )
    {
      abti = ab + 2.0 * i;
      aj[i-1] = aba + 2.0 * ( ab + i ) * ( i - 1 );
      aj[i-1] = - aj[i-1] / abti / ( abti - 2.0 );
    }

    for ( i = 2; i <= m - 1; i++ )
    {
      abti = ab + 2.0 * i;
      bj[i-1] = i * ( alpha + i ) / ( abti - 1.0 ) * ( beta + i ) 
        / ( abti * abti ) * ( ab + i ) / ( abti + 1.0 );
    }
    bj[m-1] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      bj[i] =  sqrt ( bj[i] );
    }
  }
  else
  {
    printf ( "\n" );
    printf ( "CLASS_MATRIX - Fatal error!\n" );
    printf ( "  Illegal value of KIND = %d.\n", kind );
    exit ( 1 );
  }

  return zemu;
}
/******************************************************************************/

void imtqlx ( int n, double d[], double e[], double z[] )

/******************************************************************************/
/*
  Purpose:

    IMTQLX diagonalizes a symmetric tridiagonal matrix.

  Discussion:

    This routine is a slightly modified version of the EISPACK routine to 
    perform the implicit QL algorithm on a symmetric tridiagonal matrix. 

    The authors thank the authors of EISPACK for permission to use this
    routine. 

    It has been modified to produce the product Q' * Z, where Z is an input 
    vector and Q is the orthogonal matrix diagonalizing the input matrix.  
    The changes consist (essentially) of applying the orthogonal transformations
    directly to Z as they are generated.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

    Roger Martin, James Wilkinson,
    The Implicit QL Algorithm,
    Numerische Mathematik,
    Volume 12, Number 5, December 1968, pages 377-383.

  Parameters:

    Input, int N, the order of the matrix.

    Input/output, double D(N), the diagonal entries of the matrix.
    On output, the information in D has been overwritten.

    Input/output, double E(N), the subdiagonal entries of the 
    matrix, in entries E(1) through E(N-1).  On output, the information in
    E has been overwritten.

    Input/output, double Z(N).  On input, a vector.  On output,
    the value of Q' * Z, where Q is the matrix that diagonalizes the
    input symmetric tridiagonal matrix.
*/
{
  double b;
  double c;
  double f;
  double g;
  int i;
  int ii;
  int itn = 30;
  int j;
  int k;
  int l;
  int m=0;
  int mml;
  double p;
  double prec;
  double r;
  double s;

  prec = r8_epsilon ( );

  if ( n == 1 )
  {
    return;
  }

  e[n-1] = 0.0;

  for ( l = 1; l <= n; l++ )
  {
    j = 0;
    for ( ; ; )
    {
      for ( m = l; m <= n; m++ )
      {
        if ( m == n )
        {
          break;
        }

        if ( fabs ( e[m-1] ) <= prec * ( fabs ( d[m-1] ) + fabs ( d[m] ) ) )
        {
          break;
        }
      }
      p = d[l-1];
      if ( m == l )
      {
        break;
      }
      if ( itn <= j )
      {
        printf ( "\n" );
        printf ( "IMTQLX - Fatal error!\n" );
        printf ( "  Iteration limit exceeded\n" );
        exit ( 1 );
      }
      j = j + 1;
      g = ( d[l] - p ) / ( 2.0 * e[l-1] );
      r =  sqrt ( g * g + 1.0 );
      g = d[m-1] - p + e[l-1] / ( g + fabs ( r ) * r8_sign ( g ) );
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - l;

      for ( ii = 1; ii <= mml; ii++ )
      {
        i = m - ii;
        f = s * e[i-1];
        b = c * e[i-1];

        if ( fabs ( g ) <= fabs ( f ) )
        {
          c = g / f;
          r =  sqrt ( c * c + 1.0 );
          e[i] = f * r;
          s = 1.0 / r;
          c = c * s;
        }
        else
        {
          s = f / g;
          r =  sqrt ( s * s + 1.0 );
          e[i] = g * r;
          c = 1.0 / r;
          s = s * c;
        }
        g = d[i] - p;
        r = ( d[i-1] - g ) * s + 2.0 * c * b;
        p = s * r;
        d[i] = g + p;
        g = c * r - b;
        f = z[i];
        z[i] = s * z[i-1] + c * f;
        z[i-1] = c * z[i-1] - s * f;
      }
      d[l-1] = d[l-1] - p;
      e[l-1] = g;
      e[m-1] = 0.0;
    }
  }
/*
  Sorting.
*/
  for ( ii = 2; ii <= m; ii++ )
  {
    i = ii - 1;
    k = i;
    p = d[i-1];

    for ( j = ii; j <= n; j++ )
    {
      if ( d[j-1] < p )
      {
         k = j;
         p = d[j-1];
      }
    }

    if ( k != i )
    {
      d[k-1] = d[i-1];
      d[i-1] = p;
      p = z[i-1];
      z[i-1] = z[k-1];
      z[k-1] = p;
    }
  }
  return;
}
/******************************************************************************/

void parchk ( int kind, int m, double alpha, double beta )

/******************************************************************************/
/*
  Purpose:

    PARCHK checks parameters ALPHA and BETA for classical weight functions. 

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

  Parameters:

    Input, int KIND, the rule.
    1, Legendre,             (a,b)       1.0
    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta

    Input, int M, the order of the highest moment to
    be calculated.  This value is only needed when KIND = 8.

    Input, double ALPHA, BETA, the parameters, if required
    by the value of KIND.
*/
{
  double tmp;

  if ( kind <= 0 )
  {
    printf ( "\n" );
    printf ( "PARCHK - Fatal error!\n" );
    printf ( "  KIND <= 0.\n" );
    exit ( 1 );
  }
/*
  Check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
*/
  if ( 3 <= kind && alpha <= -1.0 )
  {
    printf ( "\n" );
    printf ( "PARCHK - Fatal error!\n" );
    printf ( "  3 <= KIND and ALPHA <= -1.\n" );
    exit ( 1 );
  }
/*
  Check BETA for Jacobi.
*/
  if ( kind == 4 && beta <= -1.0 )
  {
    printf ( "\n" );
    printf ( "PARCHK - Fatal error!\n" );
    printf ( "  KIND == 4 and BETA <= -1.0.\n" );
    exit ( 1 );
  }
/*
  Check ALPHA and BETA for rational.
*/
  if ( kind == 8 )
  {
    tmp = alpha + beta + m + 1.0;
    if ( 0.0 <= tmp || tmp <= beta )
    {
      printf ( "\n" );
      printf ( "PARCHK - Fatal error!\n" );
      printf ( "  KIND == 8 but condition on ALPHA and BETA fails.\n" );
      exit ( 1 );
    }
  }
  return;
}
/******************************************************************************/

double r8_epsilon ( )

/******************************************************************************/
/*
  Purpose:

    R8_EPSILON returns the R8 round off unit.

  Discussion:

    R8_EPSILON is a number R which is a power of 2 with the property that,
    to the precision of the computer's arithmetic,
      1 < 1 + R
    but
      1 = ( 1 + R / 2 )

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 September 2012

  Author:

    John Burkardt

  Parameters:

    Output, double R8_EPSILON, the R8 round-off unit.
*/
{
  const double value = 2.220446049250313E-016;

  return value;
}
/******************************************************************************/

double r8_gamma ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_GAMMA evaluates Gamma(X) for a real argument.

  Discussion:

    The C math library includes the GAMMA ( X ) function which should generally
    be used instead of this function.

    This routine calculates the gamma function for a real argument X.

    Computation is based on an algorithm outlined in reference 1.
    The program uses rational functions that approximate the gamma
    function to at least 20 significant decimal digits.  Coefficients
    for the approximation over the interval (1,2) are unpublished.
    Those for the approximation for 12 <= X are from reference 2.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by William Cody, Laura Stoltz.
    C version by John Burkardt.

  Reference:

    William Cody,
    An Overview of Software Development for Special Functions,
    in Numerical Analysis Dundee, 1975,
    edited by GA Watson,
    Lecture Notes in Mathematics 506,
    Springer, 1976.

    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
    Charles Mesztenyi, John Rice, Henry Thatcher,
    Christoph Witzgall,
    Computer Approximations,
    Wiley, 1968,
    LC: QA297.C64.

  Parameters:

    Input, double X, the argument of the function.

    Output, double R8_GAMMA, the value of the function.
*/
{
  double c[7] = {
   -1.910444077728E-03, 
    8.4171387781295E-04, 
   -5.952379913043012E-04, 
    7.93650793500350248E-04, 
   -2.777777777777681622553E-03, 
    8.333333333333333331554247E-02, 
    5.7083835261E-03 };
  double eps = 2.22E-16;
  double fact;
  int i;
  int n;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01, 
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02, 
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04, 
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  int parity;
  const double pi = 3.1415926535897932384626434;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02, 
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03, 
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03, 
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  double res;
  const double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;

  parity = 0;
  fact = 1.0;
  n = 0;
  y = x;
/*
  Argument is negative.
*/
  if ( y <= 0.0 )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != 0.0 )
    {
      if ( y1 != ( double ) ( int ) ( y1 * 0.5 ) * 2.0 )
      {
        parity = 1;
      }

      fact = - pi / sin ( pi * res );
      y = y + 1.0;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
/*
  Argument is positive.
*/
  if ( y < eps )
  {
/*
  Argument < EPS.
*/
    if ( xminin <= y )
    {
      res = 1.0 / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < 12.0 )
  {
    y1 = y;
/*
  0.0 < argument < 1.0.
*/
    if ( y < 1.0 )
    {
      z = y;
      y = y + 1.0;
    }
/*
  1.0 < argument < 12.0.
  Reduce argument if necessary.
*/
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - 1.0;
    }
/*
  Evaluate approximation for 1.0 < argument < 2.0.
*/
    xnum = 0.0;
    xden = 1.0;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + 1.0;
/*
  Adjust result for case  0.0 < argument < 1.0.
*/
    if ( y1 < y )
    {
      res = res / y1;
    }
/*
  Adjust result for case 2.0 < argument < 12.0.
*/
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + 1.0;
      }
    }
  }
  else
  {
/*
  Evaluate for 12.0 <= argument.
*/
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - 0.5 ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
/*
  Final adjustments and return.
*/
  if ( parity )
  {
    res = - res;
  }

  if ( fact != 1.0 )
  {
    res = fact / res;
  }

  value = res;

  return value;
}
/******************************************************************************/

double r8_sign ( double x )

/******************************************************************************/
/*
  Purpose:

    R8_SIGN returns the sign of an R8.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 October 2004

  Author:

    John Burkardt

  Parameters:

    Input, double X, the number whose sign is desired.

    Output, double R8_SIGN, the sign of X.
*/
{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  } 
  else
  {
    value = 1.0;
  }
  return value;
}
/******************************************************************************/

void r8mat_write ( char *output_filename, int m, int n, double table[] )

/******************************************************************************/
/*
  Purpose:

    R8MAT_WRITE writes an R8MAT file.

  Discussion:

    An R8MAT is an array of R8's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    01 June 2009

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output filename.

    Input, int M, the spatial dimension.

    Input, int N, the number of points.

    Input, double TABLE[M*N], the data.
*/
{
  int i;
  int j;
  FILE *output;
/*
  Open the file.
*/
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8MAT_WRITE - Fatal error!\n" );
    fprintf ( stderr, "  Could not open the output file.\n" );
    exit ( 1 );
  }
/*
  Write the data.
*/
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %24.16g", table[i+j*m] );
    }
    fprintf ( output, "\n" );
  }
/*
  Close the file.
*/
  fclose ( output );

  return;
}
/******************************************************************************/

void rule_write ( int order, char *filename, double x[], double w[],
  double r[] )

/******************************************************************************/
/*
  Purpose:

    RULE_WRITE writes a quadrature rule to three files.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 October 2012

  Author:

    John Burkardt

  Parameters:

    Input, int ORDER, the order of the rule.

    Input, double A, the left endpoint.

    Input, double B, the right endpoint.

    Input, char *FILENAME, specifies the output filenames.
    "filename_w.txt", "filename_x.txt", "filename_r.txt"
    defining weights, abscissas, and region.
*/
{
  char filename_r[80];
  char filename_w[80];
  char filename_x[80];

  strcpy ( filename_r, filename );
  strcat ( filename_r, "_r.txt" );
  strcpy ( filename_w, filename );
  strcat ( filename_w, "_w.txt" );
  strcpy ( filename_x, filename );
  strcat ( filename_x, "_x.txt" );

  printf ( "\n" );
  printf ( "  Creating quadrature files.\n" );
  printf ( "\n" );
  printf ( "  Root file name is     \"%s\".\n", filename );
  printf ( "\n" );
  printf ( "  Weight file will be   \"%s\".\n", filename_w );
  printf ( "  Abscissa file will be \"%s\".\n", filename_x );
  printf ( "  Region file will be   \"%s\".\n", filename_r );

  r8mat_write ( filename_w, 1, order, w );
  r8mat_write ( filename_x, 1, order, x );
  r8mat_write ( filename_r, 1, 2,     r );

  return;
}
/******************************************************************************/

void scqf ( int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], 
  double swts[], double st[], int kind, double alpha, double beta, double a, 
  double b )

/******************************************************************************/
/*
  Purpose:

    SCQF scales a quadrature formula to a nonstandard interval.

  Discussion:

    The arrays WTS and SWTS may coincide.

    The arrays T and ST may coincide.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

  Parameters:

    Input, int NT, the number of knots.

    Input, double T[NT], the original knots.

    Input, int MLT[NT], the multiplicity of the knots.

    Input, double WTS[NWTS], the weights.

    Input, int NWTS, the number of weights.

    Input, int NDX[NT], used to index the array WTS.  
    For more details see the comments in CAWIQ.

    Output, double SWTS[NWTS], the scaled weights.

    Output, double ST[NT], the scaled knots.

    Input, int KIND, the rule.
    1, Legendre,             (a,b)       1.0
    2, Chebyshev Type 1,     (a,b)       ((b-x)*(x-a))^(-0.5)
    3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
    4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
    5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
    6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
    7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
    8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta

    Input, double ALPHA, the value of Alpha, if needed.

    Input, double BETA, the value of Beta, if needed.

    Input, double A, B, the interval endpoints.
*/
{
  double al;
  double be;
  int i;
  int k;
  int l;
  double p;
  double shft;
  double slp;
  double temp;
  double tmp;

  temp = r8_epsilon ( );

  parchk ( kind, 1, alpha, beta );

  if ( kind == 1 )
  {
    al = 0.0;
    be = 0.0;
    if ( ( b - a ) <= temp )
    {
      printf ( "\n" );
      printf ( "SCQF - Fatal error!\n" );
      printf ( "  B - A too small.\n" );
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 2 )
  {
    al = -0.5;
    be = -0.5;
    if ( ( b - a ) <= temp )
    {
      printf ( "\n" );
      printf ( "SCQF - Fatal error!\n" );
      printf ( "  B - A too small.\n" );
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 3 )
  {
    al = alpha;
    be = alpha;
    if ( ( b - a ) <= temp )
    {
      printf ( "\n" );
      printf ( "SCQF - Fatal error!\n" );
      printf ( "  B - A too small.\n" );
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 4 )
  {
    al = alpha;
    be = beta;

    if ( ( b - a ) <= temp )
    {
      printf ( "\n" );
      printf ( "SCQF - Fatal error!\n" );
      printf ( "  B - A too small.\n" );
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 5 )
  {
    if ( b <= 0.0 )
    {
      printf ( "\n" );
      printf ( "SCQF - Fatal error!\n" );
      printf ( "  B <= 0\n" );
      exit ( 1 );
    }
    shft = a;
    slp = 1.0 / b;
    al = alpha;
    be = 0.0;
  }
  else if ( kind == 6 )
  {
    if ( b <= 0.0 )
    {
      printf ( "\n" );
      printf ( "SCQF - Fatal error!\n" );
      printf ( "  B <= 0.\n" );
      exit ( 1 );
    }
    shft = a;
    slp = 1.0 / sqrt ( b );
    al = alpha;
    be = 0.0;
  }
  else if ( kind == 7 )
  {
    al = alpha;
    be = 0.0;
    if ( ( b - a ) <= temp )
    {
      printf ( "\n" );
      printf ( "SCQF - Fatal error!\n" );
      printf ( "  B - A too small.\n" );
      exit ( 1 );
    }
    shft = ( a + b ) / 2.0;
    slp = ( b - a ) / 2.0;
  }
  else if ( kind == 8 )
  {
    if ( a + b <= 0.0 )
    {
      printf ( "\n" );
      printf ( "SCQF - Fatal error!\n" );
      printf ( "  A + B <= 0.\n" );
      exit ( 1 );
    }
    shft = a;
    slp = a + b;
    al = alpha;
    be = beta;
  }

  p = pow ( slp, al + be + 1.0 );

  for ( k = 0; k < nt; k++ )
  {
    st[k] = shft + slp * t[k];
    l = abs ( ndx[k] );

    if ( l != 0 )
    {
      tmp = p;
      for ( i = l - 1; i <= l - 1 + mlt[k] - 1; i++ )
      {
        swts[i] = wts[i] * tmp;
        tmp = tmp * slp;
      }
    }
  }
  return;
}
/******************************************************************************/

void sgqf ( int nt, double aj[], double bj[], double zemu, double t[], 
  double wts[] )

/******************************************************************************/
/*
  Purpose:

    SGQF computes knots and weights of a Gauss Quadrature formula.

  Discussion:

    This routine computes all the knots and weights of a Gauss quadrature
    formula with simple knots from the Jacobi matrix and the zero-th
    moment of the weight function, using the Golub-Welsch technique.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    11 January 2010

  Author:

    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
    C version by John Burkardt.

  Reference:

    Sylvan Elhay, Jaroslav Kautsky,
    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
    Interpolatory Quadrature,
    ACM Transactions on Mathematical Software,
    Volume 13, Number 4, December 1987, pages 399-415.

  Parameters:

    Input, int NT, the number of knots.

    Input, double AJ[NT], the diagonal of the Jacobi matrix.

    Input/output, double BJ[NT], the subdiagonal of the Jacobi 
    matrix, in entries 1 through NT-1.  On output, BJ has been overwritten.

    Input, double ZEMU, the zero-th moment of the weight function.

    Output, double T[NT], the knots.

    Output, double WTS[NT], the weights.
*/
{
  int i;
/*
  Exit if the zero-th moment is not positive.
*/
  if ( zemu <= 0.0 )
  {
    printf ( "\n" );
    printf ( "SGQF - Fatal error!\n" );
    printf ( "  ZEMU <= 0.\n" );
    exit ( 1 );
  }
/*
  Set up vectors for IMTQLX.
*/
  for ( i = 0; i < nt; i++ )
  {
    t[i] = aj[i];
  }
  wts[0] = sqrt ( zemu );
  for ( i = 1; i < nt; i++ )
  {
    wts[i] = 0.0;
  }
/*
  Diagonalize the Jacobi matrix.
*/
  imtqlx ( nt, t, bj, wts );

  for ( i = 0; i < nt; i++ )
  {
    wts[i] = wts[i] * wts[i];
  }

  return;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len=0;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len=strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}


/*--------------------------------------------------------------------------------*/
/**********************************************************************************/


/**********************************************************************************/
/*--------------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------------*/
/**********************************************************************************/


/**********************************************************************************/
/*--------------------------------------------------------------------------------*/
