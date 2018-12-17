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
//#include <cblas.h>
//#include <lapacke.h>

#include "mmio.h"
#include "clock.h"
#include "coo.h"
#include "csr.h"

/* Reading the Matrix from MM. The matrix is in CSR Format*/

unsigned int maxit=1000;

void parse_args(int argc, char *argv[]);
void print_vector(char* pre, double *v, unsigned int size);
double randf(double low,double high);
double **dmatrix ( int nrl, int nrh, int ncl, int nch );
void free_dmatrix ( double **m, int nrl, int nrh, int ncl, int nch );
double  vecnorm(int n, double a1[], double a2[]); 
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
    double *w,*e, *relres;
    double **B, **R0, *V, **H, **E, **scal ;
    double **T;
    FILE *f;       //This file is used for reading RHS//
    int M, N, nz, nrhs;
    int work, lwork, info;
    int bloc;
    int restart, iter, maxit, m; //m = inner+p from matlab code//
    int *I, *J, k, j;
    double *val, *tau;
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
     printf("\n\nReading Matrix Market file Data\n");
      // Initialize matrices.
      csr_init_matrix(&csr);
      
     // Load matrix from file into COO.
      printf("Loading matrix \"%s\"\n", filename);
      sym = csr_load_matrix(filename, &csr);
    
     if(sym) printf("\n\nMatrix is symmetric\n");
     else    printf("Matrix is general (non-symmetric)\n");
 
      // Print Matrix data
//      printf("CSR matrix data:\n");
  //    csr_print_matrix(&csr);

       printf("\n\nReading : Successful\n");
/*****************************************
 *Random Right Hand Side Generation 
 ****************************************/

    //Random Matrix Generation for RHS//
// printf("Enter the restart value = \n ");
// scanf("%d", &restart);
restart = 10;

printf("\nEnter the desired  no of right hand sides for matrix B\n");
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

//Temperoray Transpose of B 
  T = dmatrix(0,nrhs,0, csr.rows);

   for(i=0; i<csr.rows; ++i)
        for(j=0; j<nrhs; ++j)
        {
            T[j][i] = B[i][j];
        }

//printf("The transpose of B is\n");
//print_matrix(T, nrhs, csr.rows);


/*********************************
 * * TODO 
 *********************************/
 /*1. Norm  (Done!) 
   2. Initialize V, H and E (Done) 
   3. QR factors  (Doing) 
   4. Modified Gram-Schmit Loop (Doing) 
   5. Least Square
   6. LU factors (Preconditioning) */

//Initialization of vectors for computing norm
w = (double *)malloc(nrhs* sizeof(double));    //Allocation of vector Norm//
//e = (double *)malloc(nrhs * sizeof(double));  
//relres = (double *)malloc(nrhs * sizeof(double));  //Allocation of Relative Residual Vector//
tau = (double *)malloc(nrhs* sizeof(double));  
/********************************
*Norm and Reidual Norm
********************************/

 printf("\nThe norm of each rhs is \n");
  for (k = 0; k<nrhs;k++){
   w[k] = vecnorm(csr.rows, T[k], T[k]);
     if (w[k]==0.0){
        w[k] = 1.0;
        }
  } 

print_vector("\nnorm =\n ", w, nrhs);
printf("\n");

//Check the routine of Sparse Matrix Vector Multiplication//


/*
//for (k=0;k<csr.rows; k++){
//e[k] = vecnorm(nrhs, R0[k], R0[k]);
//}
//Residual Norm 
 //for (k=0; k<csr.rows; k++){
 //relres[k] = e[k]/w[k];
// }
//print_vector("\n Relative Residual Norm =\n ", relres, csr.rows);
/****************************************
//Initialization of V-space, H and E//
   //Start of while block// 
****************************************/
//R0 = B-A*X, where X = 0.0//

R0 = B;               //Here I will add the Multiplication of co-efficient matrix with initial guess// 
m = restart+nrhs;
iter = 0;

//for(iter <= maxit) {
//	V =  dmatrix(0, csr.rows, 0, m);    //Orthogonal Subspace V
        V = calloc(csr.rows*m, sizeof(double));
	H = dmatrix(0, m, 0, restart);      //Hessenberg Matrix
	E = dmatrix(0, m, 0, nrhs);
	scal = dmatrix(0,nrhs,0,nrhs);      //R factor from QR factorization of R0//

	for ( i = 0; i < m; i++ ) {
      		for ( j = 0; j < restart; j++ ) {
       			 H[i][j] = 0.0;
     		 }
    	}	
/*
	for ( j = 0; j <m; j++ ) {
      		for (i = 0; i < csr.rows; i++ ){
       		 V[i][j] = 0.0;
      		}
    */

	for ( i = 0; i < m; i++ ) {
       		for ( j = 0; j < nrhs; j++ ){
           	E[i][j] = 0.0;
       		}
	}

/******************************
 * Construction of V space
 * *****************************/

/*LAPACKE routine for QR factors
LAPACKE_ dgeqrf(m, n, a, lda, tau, work, lwork, info);
 Input Arguments: m = number of rows in matrix A used in the computation, 
		  n = number of columns in matrix A used in the computation.
                  a = m by n general matrix A whose QR factorization is to be computed.
                lda = is the leading dimension of the array specified for a.
               lwork = 0 for best performance
 /
lwork =0;
dgeqrf(csr.rows,nrhs,R0, csr.rows,tau,work, lwork, info);


//}  //End of for loop before QR-factors 
// Printing matrices for Debugging 

printf("\n\n The V space is \n");
print_matrix(V, csr.rows, m);
printf("\n\n The Hessenberg matrix is \n");
print_matrix(H, m, restart);
printf("\n\n The Matrix E is \n");
print_matrix(E, m, nrhs);*/

//print_matrix(V, csr.rows, m);
//Free resources  
free(w);
//free(relres);
//free(e);
//free_dmatrix(V,0, csr.rows, 0, m);
free_dmatrix(H, 0, m, 0, restart);
free_dmatrix(E, 0, m, 0, nrhs);
free_dmatrix ( B, 0, csr.rows, 0, nrhs );
free_dmatrix(T,0,nrhs,0, csr.rows);
free_dmatrix(scal,0,nrhs,0,nrhs);
free(V);

exit(EXIT_SUCCESS); //Exit the main function 
}// End of main function


/*********************************
 * Functions 
*********************************/

void print_vector(char* pre, double *v, unsigned int size){
       unsigned int i;
      printf("%s", pre);
         for(i = 0; i < size; i++){
         //printf("%.1f ", v[i]);
         printf("%e \t", v[i]);
     }
     printf("\t");
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
             maxit = i;
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
