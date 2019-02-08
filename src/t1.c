#include <stdio.h>
#include <stdlib.h>


void print_matrix(int *arr, int rows, int cols);

int main(){


int *E, *Eritz;
int i, j, ii,m;
int k_in, pd;

k_in = 10;
pd =3;

E = (int*)calloc(k_in*k_in, sizeof(int));
Eritz = (int*)malloc(pd*k_in*sizeof(int));

for(i = 0;i<k_in;i++){
    for (j =0; j<k_in;j++){
      if (i==j)
        E[i*k_in+j] = 1.0;
      else
        E[i*k_in+j] = 0.0;
  }  
}

printf("\n\nThe Identity Matrix is\n\n");
print_matrix(E, k_in, k_in);

for(i = k_in-pd,ii=0; i<k_in, ii<pd; i++, ii++){
             for(j =0;j<k_in;j++){
      
    // printf("\n\n The value assigned to ERitz is %d\n", E[i*k_in+j]);
    Eritz[ii*k_in+j] = E[i*k_in+j];
  }  
}

printf("\n\nE Ritz is\n\n");
print_matrix(Eritz, pd, k_in);


}//End of main function


/**************
Print Matrix
**************/

void print_matrix(int *arr, int rows, int cols){
 
      for(int i = 0; i <rows; i++){
          for (int j=0;j<cols; j++){
                  printf("\t%d\t",arr[i*cols+j]);
         }
        printf("\n");
      }
  }

