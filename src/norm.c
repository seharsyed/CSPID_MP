#include <stdio.h>
#include<stdlib.h>
#include<math.h>

/*file for small(pen-paper) test vectors*/

double vecnorm( int n, double a1[], double a2[]);
double **dmatrix ( int nrl, int nrh, int ncl, int nch );

int  main(){
double *x;
double sum;
int n = 6;


sum = 0.0;
x = (double *)malloc(n * sizeof(double));
for (int i = 0;i<n; i++){
x[i] = i+1.0;
}

sum = vecnorm(6, x, x);
printf("The vector norm is = %e\n", sum);
return 0;
}


/******************************************************************************/ 
 double vecnorm( int n, double a1[], double a2[])
 /******************************************************************************/
 /*Parameters:
 
     Input, int N, the number of entries in the vectors.
 
     Input, double A1[N], A2[N], the two vectors to be considered.
 
     Output, SQRT of the dot product of the vectors.
 */
 
 {
   int i;
   double value;
 
   value = 0.0;
   for ( i = 0; i < n; i++ ){
     value += a1[i] * a2[i];
   }
   return sqrt(value);
 } 

/******************************************************************************/

double **dmatrix ( int nrl, int nrh, int ncl, int nch )

/******************************************************************************/
/*
  Purpose:

    DMATRIX allocates a double matrix.

  Discussion:

    The matrix will have a subscript range m[nrl...nrh][ncl...nch] .

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 August 2006

  Author:

    Lili Ju

  Parameters:

    Input, int NRL, NRH, the low and high row indices.

    Input, int NCL, NCH, the low and high column indices.

    Output, double **DMATRIX, a doubly-dimensioned array with
    the requested row and column ranges.
*/
{
  int i;
  double **m;
  int nrow = nrh - nrl + 1;
  int ncol = nch - ncl + 1;
/*
  Allocate pointers to the rows.
*/
  m = ( double ** ) malloc ( (size_t) ( ( nrow + 1 ) * sizeof ( double* ) ) );

  if ( ! m )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DMATRIX - Fatal error!\n" );
    fprintf ( stderr, "  Failure allocating pointers to rows.\n");
    exit ( 1 );
  }
  m = m + 1;
  m = m - nrl;
/*
  Allocate each row and set pointers to them.
*/
  m[nrl] = ( double * ) malloc ( (size_t) ( ( nrow * ncol + 1 ) * sizeof ( double ) ) );

  if ( ! m[nrl] )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "DMATRIX - Fatal error!\n" );
    fprintf ( stderr, "  Failure allocating rows.\n");
    exit ( 1 );
  }
  m[nrl] = m[nrl] + 1;
  m[nrl] = m[nrl] - ncl;

  for ( i = nrl + 1; i <= nrh; i++ )
  {
    m[i] = m[i-1] + ncol;
  }
/*
  Return the pointer to the array of pointers to the rows;
*/
  return m;
}
