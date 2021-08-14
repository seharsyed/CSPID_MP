//This file contains the routines of addition of matrices

/*************************************************************
TO DO I need to add all these routines to a single header file
*************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void add_matrix(double xx[], double yy[], double result[], int rows, int cols);
void print_matrix(double *arr, int rows, int cols);

int main(){
double *A, *B, *C;
int i, j;
int rhs;

rhs = 3;
A = calloc(rhs*rhs, sizeof(double));
B = calloc(rhs*rhs, sizeof(double));

printf("Enter the elements of first matrix\n");
 
   for (i = 0; i< rhs; i++){
      for (j = 0 ; j < rhs; j++){
         scanf("%e", &A[i*rhs+j]);
 } 
}

printf("Enter the elements of B matrix\n");
  
    for (i = 0; i< rhs; i++){
       for (j = 0 ; j < rhs; j++){
  scanf("%e", &B[i*rhs+j]);
  }   
}


add_matrix(A, B, A, rhs, rhs);
print_matrix (A, rhs, rhs);


}//End of main function 

void add_matrix(double xx[], double yy[], double result[], int rows, int cols)
{
for (int i = 0; i< rows; i++){
  for (int j = 0 ; j < cols; j++){

    result[i*cols+j] = xx[i*cols+j] +yy[i*cols+j];
  // printf("%.2e\t", A[i*cols+j]);
  }  
//printf("\n");
}
}

void print_matrix(double *arr, int rows, int cols){ 
      for(int i = 0; i <rows; i++){
          for (int j=0;j<cols; j++){

                 printf("\t%e\t",arr[i*cols+j]);
 
         }
        printf("\n");
      }
  }









//}//end of main function
