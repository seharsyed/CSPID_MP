/****************************************
 Simple Block GMRES Method

Author: Sehar Naveed
        sehar.naveed@stud-inf.unibz.it
        January 2019

Supervisor: Dr. Bruno Carpentieri
            Dr. Flavio Vella

Institute: Free University of Bolzano
**************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <lapacke.h>
#include <string.h>
#include <cblas.h>
#include <blas_sparse.h>

#include "clock.h"
#include "coo.h"
#include "csr.h"

unsigned int maxit = 1000;

void get_trans(double *x, double *y, int rows, int cols);
void parse_args(int argc, char *argv[]);
double randf(double low,double high);
void print_matrix(double *arr, int rows, int cols);
double vecnorm( int n, double a1[], double a2[]);
void print_vector(char* pre, double *v, unsigned int size);
double **dmatrix ( int nrl, int nrh, int ncl, int nch );
void matriscopy (double * destmat, double * srcmat, int rowcount, int columncount);
double dot_product(double v[], double u[],  int n);
void subtract(double xx[], double yy[], double result[], int num);
void scalvec(int n, double sa, double *sx, double *sy, int incx);
void GRot(double dx, double dy, double cs, double sn);
void matrixadd(double xx[], double yy[], double result[], int num);

int main (int argc, char * argv[])
{

double *B,*X, *R0, *tmp,*tau, *scal, *Q;
double *R, *H, *V, *w, *C, *S, *tmp1;
double *nrm, *relres, *vtol, *e;
double tol, eps; 

int rows, rhs, p;
int M, N, nz, work, lwork;
int initer, iter, i, j, k;
int info,ldb,lda, k_in;
int sym, sym1;


int ret_code;
CSR_Matrix A;
//CSC_Matrix csc;

int restart, m;
FILE *fp;
char * line = NULL;
size_t len = 0;
ssize_t read;
Clock clock;
MM_typecode matcode;
char* filename;
double *temp;
int istat, stat; 
#define max(a,b) ((a)<(b)?(b):(a))

/* compute sparse matrix matrix multiplication 
*/

parse_args(argc, argv);
 
filename = argv[1];  //Passing on the file 

/***********************************************************
* Reading the Matrix from MM. The matrix is in CSR Format*
************************************************************/

/*TODO:Small tasks: File Management

       1 Change the name from csr to A and the corresponding variables  : DONE
       2 Create an array with vector tolerance : Done
       3 Shift transpose and print functions to utility and make header file 
       4 Create Residual Matrix with residual block*/

  /********************************************
   * COO to CSR Matrix - Conversion and loading
   * ********************************************/
       printf("\n\nReading Matrix Market file Data\n");
        // Initialize matrices.
       csr_init_matrix(&A);
  
       // Load matrix from file into COO.
        printf("Loading matrix \"%s\"\n", filename);
       sym = csr_load_matrix(filename, &A);
  
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
rows = A.rows;
rhs = 10;                   //Change it to get rhs from user
restart = 10;              //Change it later to get it from user
eps = 0.1;
tol = 1e-6;
p = rhs;
// m = restart+rhs; //In matlab m = inner+p

//Initialize and allocate B

B = calloc(rows*rhs, sizeof(double)); //RHS
X = calloc(rows*rhs, sizeof(double)); //Solution Array
R0 = calloc(rows*rhs, sizeof(double)); //Initial Residual Array
tmp = calloc(rhs*rows, sizeof(double)); //Temporary Array for keeping transpose of Matrices
R = calloc(rows*rhs, sizeof(double)); //Residual Array
tmp1 =  calloc(rhs*rows, sizeof(double));

nrm = calloc(rhs, sizeof(double)); 
w = calloc(rows, sizeof(double));  //Allocation of Vector Norm//
vtol = calloc(rhs, sizeof(double)); //Vector of tolerance
e = calloc(rhs, sizeof(double));

//y = calloc(rows, sizeof(double));  //Allocation of temperoray vector//

relres = (double *)malloc(rhs* sizeof(double)); // Allocation of Relative Residual//

tau = calloc(rhs,sizeof(double)); //array for LAPACK calculation

scal = calloc(rhs*rhs, sizeof(double));  //R factor of QR factorization of R 

//E = calloc(m*rhs,sizeof(double));

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
printf("\n\nTranspose of B \n");
get_trans(B,tmp,rows,rhs);
print_matrix(tmp,rhs,rows);

printf("\n\nThe norm of each RHS is \n");

ldb = rows;
  for (k = 0; k<rhs;k++){
        // for (i = 0; i <rows;i++){
    nrm[k] = vecnorm(rows,&tmp[k*ldb], &tmp[k*ldb]);
      if (nrm[k]==0.0){
         nrm[k] = 1.0;
         }
}
 
print_vector("\nnorm =\n ", nrm, rhs);
 printf("\n");

//R0 = B-A*X , since initial value is zero so here we are using R0 = B

for (k=0;k<rhs; k++){
 e[k] = vecnorm(rows, &tmp[k*ldb], &tmp[k*ldb]);
 }
 //Residual Norm 
  for (k=0; k<rhs; k++){
  relres[k] = e[k]/nrm[k];
  }
 //print_vector("\n\n Relative Residual Norm =\n ", relres, rhs);

//Initialtization of vector of tolerance
for (k =0;k<rhs;k++){
        vtol[k] = tol;
    }

// print_vector("\n\n Vector of Tolerance =\n ", vtol, rhs);

for(k = 0;k<rhs;k++){
    if(relres[k]<vtol[k])
     exit(0);
}

/*****************************
Orthogonal Space Calculation 
*****************************/

 m = restart+rhs;
 
 V = (double*)calloc(m*rows, sizeof(double)); //Orthogonal Basis
 H = (double*)calloc(m*restart, sizeof(double)); //Hessenberg Matrix
 S = (double*)calloc(restart*rhs, sizeof(double));
//Resultant Matrix after Multiplication V and S
 C = (double*)calloc(rows*rhs, sizeof(double));
/********************
 QR FACTORiZATION
*********************/


printf("\nQR factorization started\n");

//lda = rhs;
info = LAPACKE_dgeqrf(LAPACK_ROW_MAJOR, rows, rhs, B, rhs, tau ); //QR Factors 

//print_matrix(B,rows, rhs);


printf("The R factor from QR factorization is\n\n");
for (i=0;i<rhs;i++){
  for(j=0;j<rhs;j++){
             if(i<=j){
                scal[i*rhs+j] = B[i*rhs+j];
     }
  }
}

printf("\nThe R factor is\n");
print_matrix(scal,rhs,rhs);

/* Extracting V as Q factor of B */

printf("\n\nThe Q factor is\n");
info = LAPACKE_dorgqr(LAPACK_ROW_MAJOR, rows, rhs, rhs, B, rhs, tau);

    /*
   if (info /= 0){ 
     printf("DQRSL returns info = %d", info);
     exit;
    } 
 */

      //printf("\n\n The transpose of B--V is \n");

       // Extracting the Q factor of B//
         for (i =0;i<rows; i++){
           for (j=0;j<rhs;j++){
             V[j*rows+i] = B[i*rhs+j];
          }  
      }
     print_matrix(V,m,rows);
     printf("\n\nQR Factorization, extraction of V and R completed successfully\n");

/*****************************
Block Arnoldi Variant
******************************/

for (int initer = p;initer<m;initer++){
         k_in = initer - p;
            csr_mvp_sym2(&A,&V[k_in*rows],w); //Sparse-Matrix Vector Multiplication */ 
        /**********************
         Modified Gram-Schmidt 
         **********************/        
      for (i = 0;i <initer; i++){
                  H[i*restart+k_in]= dot_product(&V[i*rows], w, rows);
                     cblas_daxpy (rows,-H[i*restart+k_in], &V[i*rows],1, w,1);
            }
       /*************************/

           H[initer*restart+k_in] = vecnorm(rows,w, w);
           cblas_dscal(rows, 1.0/H[initer*restart+k_in],w, 1);
           cblas_dcopy(rows, w, 1, &V[initer*rows], 1); 
  
       /**************
       End of Arnoldi
       ***************/

       /**************
        Givens Rotation
       **************/

       //Reading the S matrix from MATLB Source, I need to port this step in C
       // E=[eye(p,p);zeros(k_in,p)]*scal;
       // S = H(1:k_in + p,1:k_in)\E(1:k_in + p,:); 

     //Matrix Read      
   
        fp = fopen("S.txt", "r");//Right Now I am reading S obtained from MATLAB results
        if (fp == NULL)
        exit(0);

        while(!feof(fp)){
              for(i=0;i<restart;i++){
                  for(j=0;j<p;j++){
                fscanf(fp,"%lf",&S[i*restart+j]);
         }
        }
     }  //End of while loop for reading matrix


  /*
     for (i =0; i<pd;i++){   
      scalvec(pd, Sigma_title[i], &VT[i*pd], &VT1[i*pd], 1);
      }
*/
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, rows, p, p, 1.0, V, rows, S, p, 0.0, C, p);  //V(:,1:k_in)*S = C
      //X = X0+(V*S = C)  
      //num = rows*cols

      // matrixadd(double xx[], double yy[], double result[], int num);
        matrixadd(X, C, X,rows*rhs);  
       //R = B-A*X;

       //Assigning transpose of X to tmp so that its columns gets called for multiplication

       get_trans(X, tmp1, rows, rhs);
       printf("\n\nThe transpose of X is\n\n");
      
      //print_matrix(tmp,rhs, rows);
        
//       printf("\n\nThe matrix R before multiplication is\n");
  //     print_matrix(R, rows, rhs);

        for (i = 0;i<rhs;i++){
       csr_mvp_sym(&A,&tmp1[i*rows],&R[i*rows]);     
       }
   
   subtract(tmp, R, R, rows*rhs);       

//    get_trans(R,R,rows,rhs); 
    //printf("\n\nThe matrix R is\n");
     //print_matrix(R, rows, rhs);
       //subtract( 

} //End of outer for loop

printf("\n\nThe matrix R is\n");
print_matrix(R, rhs, rows);

/************
Debugging
*************/
//printf("\n\nThe Orthogonal basis V is\n");
//print_matrix(V,m, rows);

print_vector("w is\n",w,rows);
 
printf("\n\nThe Hessenberg Matrix is\n\n");
print_matrix(H,m, restart);

printf("\n\nThe matrix S is\n");
print_matrix(S, p, p);

printf("\n\nV and S gives\n");
print_matrix(C, rows, p);

printf("\n\n The final solution is\n");
print_matrix(X, rows, p);

printf("\n\nThe matrix R is\n");
print_matrix(R, rows, rhs);
/******************************
Free Resources
******************************/
free(B);free(X);free(R0);
free(scal);free(tmp);
free(V); free(H);free(w); 
free(R);
free(w);free(nrm); free(relres);
//free(E);
free(tau);
exit(EXIT_SUCCESS);

} //End of main program 

/****************************************
*Functions 
***************************************/

void print_matrix(double *arr, int rows, int cols){
 
     for(int i = 0; i <rows; i++){
         for (int j=0;j<cols; j++){
 
                printf("\t%e\t",arr[i*cols+j]);

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
//          printf("%.3f\t ", v[i]);
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

double dot_product(double v[], double u[], int n)
{
    double result = 0.0;
    for (int i = 0; i < n; i++)
        result += v[i]*u[i];
    return result;
}

void subtract(double xx[], double yy[], double result[], int num) {
    for (int ii = 0; ii < num; ii++) {
        result[ii] = xx[ii]-yy[ii];
    }
}

void matrixadd(double xx[], double yy[], double result[], int num) {
     for (int ii = 0; ii < num; ii++) {
         result[ii] = xx[ii]+yy[ii];
     }
 }


void get_trans(double *x, double *y, int rows, int cols)
{
      unsigned int i,j;
         for(i=0; i<rows; ++i){
            for(j=0; j<cols; ++j){
              y[j*rows+i] = x[i*cols+j];
          }
 } 
}


void scalvec(int n, double sa, double *sx, double *sy, int incx)
{
  long int i, m, nincx, nn, iincx;
  double ssa;

    // scales a vector by a constant.   
    // uses unrolled loops for increment equal to 1.   
    // jack dongarra, linpack, 3/11/78.   
    // modified 3/93 to return if incx .le. 0.   
    // modified 12/3/93, array(1) declarations changed to array(*) 

  // Dereference inputs 
  nn = n;
  iincx = incx;
  ssa = sa;

  if (nn > 0 && iincx > 0)
  {
    if (iincx == 1) // code for increment equal to 1 
    {
      m = nn-4;
      for (i = 0; i < m; i += 5)
      {
        sy[i] = ssa * sx[i];
        sy[i+1] = ssa * sx[i+1];
        sy[i+2] = ssa * sx[i+2];
        sy[i+3] = ssa * sx[i+3];
        sy[i+4] = ssa * sx[i+4];
      }
      for ( ; i < nn; ++i) // clean-up loop 
        sy[i] = ssa * sx[i];
    }
    else // code for increment not equal to 1 
    {
      nincx = nn * iincx;
      for (i = 0; i < nincx; i += iincx)
        sy[i] = ssa * sx[i];
    }
  }

} 

void GRot(double dx, double dy, double cs, double sn)
{
  double temp;

  if (dy == 0.0) {
    cs = 1.0;
    sn = 0.0;
  } 
   else if (abs(dy) > abs(dx)) {
    temp = dx / dy;
    sn = 1.0 / sqrt( 1.0 + temp*temp );
    cs = temp * sn;
  } 
   else {
    temp = dy / dx;
    cs = 1.0 / sqrt( 1.0 + temp*temp );
    sn = temp * cs;
  }
}
