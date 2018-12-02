/************************************
  This is going to be the main file
 *************************************
 Author:  Sehar Naveed
 Institute:  Unibz
 Date: November 12, 2018
***************************************/
/*Implementation of Block Arnoldi*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>

#include "mmio.h"
#include "clock.h"
#include "coo.h"
#include "csr.h"
#include "f2c.h"

/* Reading the Matrix from MM. The matrix is in CSR Format*/

unsigned int iterations=1000;

void parse_args(int argc, char *argv[]);
void print_vector(char* pre, double *v, unsigned int size);
double randf(double low,double high);
double **dmatrix ( int nrl, int nrh, int ncl, int nch );
void free_dmatrix ( double **m, int nrl, int nrh, int ncl, int nch );
double  vecnorm(int n, double *a1, double *a2); 
void print_matrix(double **arr, int rows, int cols);

int main(int argc, char *argv[])
{
    unsigned int i;
    int sym;
    int ret_code;
    CSR_Matrix csr;
    //COO_Matrix coo;//
    //int sym1;//
    Clock clock;
    MM_typecode matcode;
    char* filename;
    double *w;
    double **B, *y;
    FILE *f;       //This file is used for reading RHS//
    int M, N, nz, nrhs;
    int bloc;
    int *I, *J, k, j;
    double *val;
    double a;

parse_args(argc, argv);

 filename = argv[1];

 /***************************************
  * If the matrix Instance is COO
  * ***********************************/

      //Initialize COO Matrix first//
     // coo_init_matrix(&coo);

      //Load COO Matrix
     // printf("Loading matrix \"%s\"\n", filename);
      
     // sym1 = coo_load_matrix(filename, &coo);

     // Print Matrix data
     // printf("COO matrix data:\n");
     // coo_print_matrix(&coo);

/********************************************
 * COO to CSR Matrix - Conversion and loading
 * ********************************************/ 
      // Initialize matrices.
      csr_init_matrix(&csr);
      
     // Load matrix from file into COO.
      printf("Loading matrix \"%s\"\n", filename);
      sym = csr_load_matrix(filename, &csr);
    
     if(sym) printf("Matrix is symmetric\n");
     else    printf("Matrix is general (non-symmetric)\n");
 
      // Print Matrix data
      printf("CSR matrix data:\n");
      csr_print_matrix(&csr);


/*****************************************
 *Random Right Hand Side Generation 
 ****************************************/

    //Random Matrix Generation for RHS//

 printf("Enter the desired  no of right hand sides for matrix B\n");
 scanf("%d",&nrhs);
   
   printf(" The value of rows = %d\n",csr.rows);
   printf(" The value of cols = %d\n",nrhs);
 
  B = dmatrix(0,csr.rows, 0, nrhs);

  for(i = 0; i <csr.rows; i++){
	  for (j=0;j<nrhs; j++){
              B[i][j] = randf(0,1);
      }
    }

print_matrix(B, csr.rows, nrhs);


/*********************************
 * * TODO 
 *********************************/
 /*1. Norm  (status done! but values are not precise something to do with pointer perhaps) 
   2. Initialize V, H and E
   3. QR factors
   4. Modified Gram-Schmit Loop
   5. Least Square
   6. LU factors (Preconditioning) */

//First initialize a vector for computing norm//
w = (double *)malloc(csr.rows * sizeof(double));
  
 printf("\nThe norm of each rhs is \n");

  for (k = 0; k<nrhs;k++){
   w[k] = vecnorm(csr.rows, B[k], B[k]);
  } 

 print_vector("\nnorm =\n ", w, csr.rows);





//Free resources  
//free(w);
free_dmatrix ( B, 0, csr.rows, 0, nrhs );
exit(EXIT_SUCCESS); //Exit the main function 
}


/*********************************
 * Functions 
*********************************/

void print_vector(char* pre, double *v, unsigned int size){
       unsigned int i;
      printf("%s", pre);
         for(i = 0; i < size; i++){
         //printf("%.1f ", v[i]);
         printf("%e \n", v[i]);
     }
     printf("\n");
 }

     void parse_args(int argc, char *argv[]){
     int i;
     if (argc < 2) {
         printf("Usage: %s input_matrix [iterations]\n", argv[0]);
         exit(EXIT_FAILURE);
     }
       if(argc >= 3){
         i = atoi(argv[2]);
         if (i <= 0){
             printf("Invalid number of iterations.\n");
             exit(EXIT_FAILURE);
         }
             iterations = i;
     }
 }

double randf(double low,double high){
 return (rand()/(double)(RAND_MAX))*fabs(low-high)+low;
}


void print_matrix(double **arr, int rows, int cols){
   
    for(int i = 0; i <rows; i++){
        for (int j=0;j<cols; j++){
              
               printf("\t%e\t",arr[i][j]);
       }
      printf("\n");
    }

}


/*************************************
 *Allocating 2D-Array
Parameters:

    Input, int NRL, NRH, the low and high row indices.

    Input, int NCL, NCH, the low and high column indices.

    Output, double **DMATRIX, a doubly-dimensioned array with
    the requested row and column ranges.

 ************************************/
double **dmatrix ( int nrl, int nrh, int ncl, int nch )

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

/******************************************************************************/

void free_dmatrix ( double **m, int nrl, int nrh, int ncl, int nch )

/******************************************************************************/
/*
  Purpose:

    FREE_DMATRIX frees a double matrix allocated by DMATRIX 

  Parameters:

    Input, int NRL, NRH, the low and high row indices.

    Input, int NCL, NCH, the low and high column indices.

    Input, double **M, the pointer to the doubly-dimensioned array,
    previously created by a call to DMATRIX.
*/
{
  free ( ( char * ) ( m[nrl] + ncl - 1 ) );
  free ( ( char * ) ( m + nrl - 1 ) );

  return;
}

/******************************************************************************/

double vecnorm( int n, double *a1, double *a2 )

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

