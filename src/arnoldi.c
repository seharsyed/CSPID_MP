#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <lapacke.h>
#include <string.h>

#include "clock.h"
#include "coo.h"
#include "csr.h"

unsigned int maxit = 1000;

void parse_args(int argc, char *argv[]);
double randf(double low,double high);
void print_matrix(double *arr, int rows, int cols);
double vecnorm( int n, double a1[], double a2[]);
void print_vector(char* pre, double *v, unsigned int size);
double **dmatrix ( int nrl, int nrh, int ncl, int nch );
void matriscopy (double * destmat, double * srcmat, int rowcount, int columncount);

int main (int argc, char * argv[])
{

double *B,*w,*tau, *scal, *V, *H, *E;
double *relres, *e, *nrm;
int rows, rhs;
int M, N, nz, work, lwork;
int initer, iter, i, j, k;
int info,ldb,lda, k_in;
int sym;
int ret_code;
CSR_Matrix csr;
int restart, m;
Clock clock;
MM_typecode matcode;
char* filename;
double **T1;
double *T;
/* compute sparse matrix matrix multiplication 
*/

parse_args(argc, argv);
 
filename = argv[1];  //Passing on the file 

/***********************************************************
* Reading the Matrix from MM. The matrix is in CSR Format*
************************************************************/
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

/********************************
*Initialization 
********************************/
iter = 0;
rows = csr.rows;
rhs = 10;
restart = 10;
m = restart+rhs; //In matlab m = inner+p

//Initialize and allocate B

B = calloc(rows*rhs, sizeof(double));
nrm = calloc(rhs, sizeof(double)); 
w = calloc(rows, sizeof(double));  //Allocation of Vector Norm//
relres = (double *)malloc(rhs* sizeof(double)); // Allocation of Relative Residual//
tau = calloc(rhs,sizeof(double));
scal = calloc(rhs*rhs, sizeof(double));
T = calloc(rhs*rows,sizeof(double));
//V = calloc(rhs*rows, sizeof(double));
V = calloc(m*rows,sizeof(double));
H = calloc(m*restart, sizeof(double));
E = calloc(m*rhs,sizeof(double));

/*******************************
*Generate Random RHS Matrix 
*******************************/

for(i =0; i<rows;i++){
   for (j = 0; j<rhs; j++){
               B[i*rhs+j]= randf(0,1);
    }
}

printf("\n\nThe randomly generated RHS is\n"); 
print_matrix(B,rows,rhs);

/***********************************************************
*Transpose of B/ Calculation Norm and Relative Residual Norm 
************************************************************/

for(i=0; i<rows; ++i){
        for(j=0; j<rhs; ++j){
             T[j*rows+i] = B[i*rhs+j];
         }
}


/*T = dmatrix(0,rhs,0, csr.rows);
k = 0; 
    for(i=0; i<csr.rows; i++){
         for(j=0; j<rhs; j++){
            T[j][i] = B[k];
            k++;
        }
}
/*
printf("The transpose of B is\n");
print_matrix(T, rhs, rows);
*/


printf("\n\nThe norm of each RHS is \n");

ldb = rows;
  for (k = 0; k<rhs;k++){
        // for (i = 0; i <rows;i++){
    nrm[k] = vecnorm(rows,&T[k*ldb], &T[k*ldb]);
      if (nrm[k]==0.0){
         nrm[k] = 1.0;
         }
}
 
print_vector("\nnorm =\n ", nrm, rhs);
 printf("\n");


/*************************
*Relative Residual Norm 
**************************/

/*
for (k=0;k<csr.rows; k++){
 e[k] = vecnorm(nrhs, R0[k], R0[k]);
 }
 //Residual Norm 
  for (k=0; k<csr.rows; k++){
  relres[k] = e[k]/w[k];
  }
 print_vector("\n Relative Residual Norm =\n ", relres, csr.rows);
*/


/********************
QR FACTORS 
*********************/


printf("\nQR factorization started\n");

lda = rhs;
info = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, rows, rhs, B, lda, tau );

print_matrix(B,rows, rhs);


printf("The R factor is\n\n");
for (i=0;i<rhs;i++){
  for(j=0;j<rhs;j++){
             if(i<=j){
                scal[i*rhs+j] = B[i*rhs+j];
     }
  }
}

print_matrix(scal,rhs,rhs);

/* Extracting V as Q factor of B */

printf("\n\nThe Q factor is\n");
info = LAPACKE_dorgqr(LAPACK_ROW_MAJOR, rows, rhs, rhs, B, lda, tau);

/*
if (info /= 0){ 
    printf("DQRSL returns info = %d", info);
    exit;
  } 
*/

print_matrix(B, rows, rhs);
//printf("\n\n The transpose of B--V is \n");

// Allocating the transpose of Q to V//
for (i =0;i<rows; i++){
   for (j=0;j<rhs;j++){
        V[j*rows+i] = B[i*rhs+j];
}
}
//print_matrix(V, rhs, rows);


/********************************************************************
/*Pointer Artithmetic for allocating values to V
for(i=0;i<rows;i++){
  for(j=0;j<rhs;j++){
      V[i*rhs+j]=B[i*rhs+j];
          // if(j==9){
         printf("\n\n I am here and my value is %d\n", j);         
          V[i*rhs+restart+j]=B[i*rhs+j];
              }
} 
}
*/
/*
for (j = 0;j<rhs;j++){
   for (i = 0;i<rows;i++){
       printf("%.2f\t", B[i*rhs+j]);
}
printf("\n");
}  I will come to it later to understand  why this technique didn't work 
************************************************************************

/*****************************
Modified Gram-Schmidt Portion
******************************/
/*
csr_mvp_sym2(&csr,&V[9],w);
printf("\nThe first row of V orthogonal is \nn");
print_vector("\nV[3]\n", &V[7], rows);
 
//Sparse Matrix Vector Multiplication 
/*for (int initer = rhs;initer<m-2;initer++){
         k_in = initer - rhs+1;
            csr_mvp_sym2(&csr,&V[k_in],w);

}*/

//print_vector("\nw =\n ", w, rows);
/*
for(j = 0;j <rhs;j++){
   for(i = 0;i <rows;i++){
      w[j] = vecnorm(i, &B[j*rows], &B[j*rows]);
   }  
}
*/

//print_vector("\nnorm =\n ", w, rhs);

/*************************************
Printing Matrix for Debugging
*************************************/

printf("\n\nThe Orthogonal basis V is:\n");
print_matrix(V,m,rows);

/*
printf("\n\nThe Hessenberg H is:\n");
print_matrix(H,m,restart);

printf("\n\nThe Identity matrix E is: \n");
print_matrix(E,m,rhs);
*/
printf("\n\n");

/******************************
Free Resources
******************************/
free(B);
free(w);
free(V);
free(H);
free(E);
free(tau);
free(scal);

exit(EXIT_SUCCESS);

} //End of main program 
/****************************************
*Functions 
***************************************/

void print_matrix(double *arr, int rows, int cols){
 
     for(int i = 0; i <rows; i++){
         for (int j=0;j<cols; j++){
 
                printf("\t%.2f\t",arr[i*cols+j]);

        }
       printf("\n");
     } 
 }

double randf(double low,double high){
return (rand()/(double)(RAND_MAX))*fabs(low-high)+low;
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


void print_vector(char* pre, double *v, unsigned int size){
        unsigned int i;
       printf("%s", pre);
         for(i = 0; i < size; i++){
//          printf("%.2f\t ", v[i]);
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
 /* 
   Return the pointer to the array of pointers to the rows;
 */
   return m;
 }


void matriscopy (double * destmat, double * srcmat, int rowcount, int columncount)
{
  int i, j;
  double (*dst)[columncount];
  double (*src)[columncount];
  dst = (double (*)[columncount])destmat;
  src = (double (*)[columncount])srcmat;
  for (j=0; j<columncount; j++) /* rad-nr */
    for (i=0; i<rowcount; i++) /* kolumn-nr */
      dst[j][i] = src[j][i];
}
