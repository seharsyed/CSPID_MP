
/*********************************
Implementation of Block Arnoldi
*********************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <string.h>
#include <stddef.h>

#include <mmio.h>
#include "clock.h"
#include "coo.h"
#include "csr.h"

unsigned int iterations=1000;

double randf(double low,double high);
void print_matrix(double *arr, int rows, int cols);
double vecnorm( int n, double *a1, double *a2);
void print_vector(char* pre, double *v, unsigned int size);
void matriscopy (double * destmat, double * srcmat, int rowcount, int columncount);
void parse_args(int argc, char *argv[]);

int main (int argc, const char * argv[])
{

int sym;
int ret_code;
CSR_Matrix csr;
Clock clock;
MM_typecode matcode;
char* filename;


double *B,*w,*tau, *scal, *V, *H, *E;
int rows, rhs;
int i, j, k, k_in;
int info, lda;
int restart, m;

/* Reading the Matrix from MM. The matrix is in CSR Format*/

/* compute sparse matrix matrix multiplication 
*/
parse_args(argc, argv);

filename = argv[1];

/********************************************
  * COO to CSR Matrix - Conversion and loading
 **********************************************/
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
   

rows = csr.rows;
//rows = 48;
rhs = 10;
restart = 10;
m = restart+rhs; //In matlab m = inner+p


/*******************************
//Initialization and allocation
********************************/

B = calloc(rows*rhs, sizeof(double));
w = (double *)malloc(rows* sizeof(double)); 
tau = calloc(rhs,sizeof(double));
scal = calloc(rhs*rhs, sizeof(double));

//V = calloc(rhs*rows, sizeof(double));
V = calloc(m*rows,sizeof(double));
H = calloc(m*restart, sizeof(double));
E = calloc(m*rhs,sizeof(double));

/*********************************
Random Matrix Generation for RHS 
*********************************/

for(i =0; i<rows;i++){
   for (j = 0; j<rhs; j++){
          // for(j =0; j<rhs;j++){
               B[i*rhs+j]= randf(0,1);
    // printf("\t%e\t",B[i*rhs+j]);
    }
 // printf("\n");
} 

print_matrix(B,rows,rhs);


/***********************************
QR Factorization 
***********************************/

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
printf("\n\n The transpose of B--V is \n");

// Allocating the transpose of Q to V//
for (i =0;i<rows; i++){
   for (j=0;j<rhs;j++){
        V[j*rows+i] = B[i*rhs+j];
}
}
print_matrix(V, rhs, rows);


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

        //Sparse Matrix Vector Multiplication 
        printf("\n\nCalculating Matrix Vector product of A and each column of V\n");

        clock_start(&clock);

     for(i=rhs; i < m-1; i++){
        k_in = i-rhs+1;
        csr_mvp(&csr,&V[k_in],w);
	}
       clock_stop(&clock);
       printf("CSR mvp2 time: %.2fus \n\n", clock.sec);

//print_vector("\nnorm =\n ", w, rhs);

/*************************************
Printing Matrix for Debugging
*************************************/

printf("\n\nThe Orthogonal basis V is:\n");
print_matrix(V,m,rows);

printf("\n\nThe Hessenberg H is:\n");
print_matrix(H,m,restart);

//printf("\n\nThe Identity matrix E is: \n");
//print_matrix(E,m,rhs);

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

return 0;

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
 double vecnorm( int n, double *a1, double *a2)
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
              iterations = i;
      }
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
