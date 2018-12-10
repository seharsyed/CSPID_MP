#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

double randf(double low,double high);
void print_matrix(double *arr, int rows, int cols);
double vecnorm( int n, double *a1, double *a2);
void print_vector(char* pre, double *v, unsigned int size);

int main (int argc, const char * argv[])
{

double *B,*w,*tau, *scal;
int rows, rhs;
int i, j, k;
int info, lda;
/* I am initializing B as a single pointer matrix as LAPACKE only accepts single pointer. 
Generate B as random matrix in column major 
After computing QR factors check the Fortran code for extracting Q 
Find a way to save csr matrix as single pointer instead of double! 
then compute sparse matrix matrix multiplication 

*/

rows = 48;
rhs = 10;

//Initialize and allocate B

B = calloc(rows*rhs, sizeof(double));
w = (double *)malloc(rhs* sizeof(double)); 
tau = calloc(rows,sizeof(double));
scal = calloc(rhs*rhs, sizeof(double));

for(j =0; j<rhs;j++){
   for (i = 0; i<rows; i++){
          // for(j =0; j<rhs;j++){
               B[i*rhs+j]= randf(0,1);
    // printf("\t%e\t",B[i*rhs+j]);
    }
 // printf("\n");
} 

print_matrix(B,rows,rhs);

printf("\nQR factorization started\n");

lda = rows;
info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, rows, rhs, B, lda, tau );

print_matrix(B,rows, rhs);

/*
printf("The R factor is\n\n");
for (j=0;j<rhs;j++){
  for(i=0;i<rhs;i++){
             if(i>j){
                scal[i+j*rhs] = B[i+j*rhs];
     }
  }
}

print_matrix(scal,rhs,rhs);

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
printf("\n\n");

free(B);
free(w);
free(scal);

return 0;

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


