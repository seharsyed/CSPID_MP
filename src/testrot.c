//This is the file made to test the implementation of givens rotation on a small hessenberg matrix
//I will read a small H matrix then apply the givens rotation on it. Check the same result on Matlab. 

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <string.h>

int main ()
{

double *c;
double *s;
double *h;
double *y; // starting as a vector and will implement it later for matrix
double *g;

int m; //maximum no of inner iterations, mr < size of matrix 
//mr = m,
int restart;

int i;
int j;
int k;

    printf("Enter number of Rows :");
    scanf("%d",&m);
    printf("Enter number of Cols :");
    scanf("%d",&restart);

c = ( double * ) malloc ( m * sizeof ( double ) );
s = ( double * ) malloc ( m * sizeof ( double ) );
g = ( double * ) malloc ( ( m + 1 ) * sizeof ( double ) );
h = (double * ) malloc ((m*restart)* sizeof (double));

printf("\nEnter matrix elements :\n");
    for(i=0;i< m;i++)
    {
        for(j=0;j< restart;j++)
        {
            printf("Enter element [%d,%d] : ",i+1,j+1);
            scanf("%le",&h[i*restart+j]);
        }
    }

 printf("\nMatrix is :\n");
    for(i=0;i< m;i++)
    {
        for(j=0;j< restart;j++)


        {
            printf("%e\t",h[i*restart+j]);
        }
        printf("\n");   /*new line after row elements*/
    }







}// End of main function 
