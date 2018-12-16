#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include<string.h>

double randf(double low,double high);
void print_matrix(double *arr, int rows, int cols);
double vecnorm( int n, double *a1, double *a2);
void print_vector(char* pre, double *v, unsigned int size);
void matriscopy (double * destmat, double * srcmat, int rowcount, int columncount);

int main (int argc, const char * argv[])
{

double *B,*w,*tau, *scal, *V, *H, *E;
int rows, rhs;
int i, j, k;
int info, lda;
int restart, m;

/* compute sparse matrix matrix multiplication 
*/

rows = 48;
rhs = 10;
restart = 10;
m = restart+rhs; //In matlab m = inner+p

//Initialize and allocate B

B = calloc(rows*rhs, sizeof(double));
w = (double *)malloc(rhs* sizeof(double)); 
tau = calloc(rhs,sizeof(double));
scal = calloc(rhs*rhs, sizeof(double));
V = calloc(rhs*rows, sizeof(double));
H = calloc(m*restart, sizeof(double));
E = calloc(m*rhs,sizeof(double));

for(i =0; i<rows;i++){
   for (j = 0; j<rhs; j++){
          // for(j =0; j<rhs;j++){
               B[i*rhs+j]= randf(0,1);
    // printf("\t%e\t",B[i*rhs+j]);
    }
 // printf("\n");
} 

print_matrix(B,rows,rhs);

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

printf("\n\nAfter rellocation\n");
V = (double*)realloc(V,sizeof(double)*m*rows);

for (i =rhs;i<m;i++){
  for(j=0;j<rows;j++){
    V[i*rows+j] = 0.0;
}
}


print_matrix(V, m, rows);

/*
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
 for (i = 0; i<rows; i++){
        int width = (B[i]) + 1;
*/


//V = realloc(V, rows*m*sizeof(double));
//memcpy (V, B, sizeof(double)*rows*rhs);
/*
for (j = 0;j<rhs;j++){
   for (i = 0;i<rows;i++){
       printf("%.2f\t", B[i*rhs+j]);
}
printf("\n");
}

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
/*
printf("\n\nThe Orthogonal basis V is:\n");
print_matrix(V,rows,m);
*/
printf("\n\nThe Hessenberg H is:\n");
print_matrix(H,m,restart);

printf("\n\nThe Identity matrix E is: \n");
print_matrix(E,m,rhs);

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
