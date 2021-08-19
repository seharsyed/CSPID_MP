//File for checking backward substitution and updated solution 
//

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <string.h>



void print_vector(char* pre, double *v, unsigned int size);
void print_matrix(double *arr, int rows, int cols);
void GeneratePlaneRotation(double dx, double dy, double cs, double sn);
void mult_givens(double cs, double sn, int k, double *g);

//void ApplyPlaneRotation(double dx, double dy, double cs, double sn);

int main ()
{

double *cs;
double *sn;
double *h;
double *y; // starting as a vector and will implement it later for matrix
double *g;
double *e;

int m; //maximum no of inner iterations, mr < size of matrix 
//mr = m,
int restart;

int i;
int j;
int k;
double temp;
double h1;
double h2;

    printf("Enter number of Rows :");
    scanf("%d",&m);
    printf("Enter number of Cols :");
    scanf("%d",&restart);

cs = ( double * ) malloc ( m * sizeof ( double ) );
sn = ( double * ) malloc ( m * sizeof ( double ) );
//g = ( double * ) malloc ( ( m + 1 ) * sizeof ( double ) );//rotation matrix
y = ( double * ) malloc ( ( m + 1 ) * sizeof ( double ) );
h = (double * ) malloc ((m*restart)* sizeof (double));
g = (double *) malloc (m*m*sizeof(double)); //Identity matrix

//Declaration of matrix from user 
   printf("\nEnter matrix elements :\n");
    for(i=0;i< m;i++)
    {
        for(j=0;j< restart;j++)
        {
            printf("Enter element [%d,%d] : ",i+1,j+1);
            scanf("%le",&h[i*restart+j]);
        }
    }

//Print the entered matrix
printf("\nThe entered matrix is: \n");
print_matrix(h,m, restart);



    //Declaration of rotation matrix I am starting with identity and then discover the real 



//declaration of identity matrix 

   for (k = 0; k < m; k++){
    g[k*m+k] = 1.0;
    }

//Reading the identity matrix
printf("\n\nThe matrix of identity is: \n\n");
print_matrix(g,m,m);


//Starting the rotation here 
/*********************************************************************************
//This routine is working but somehow picking up the wrong values in second row
//One way is to double check and print the values that it picks at each step
**********************************************************************/
for (k =0;k<m;k++){

//for (k = 0;k<m;k++){

 if ( 0 < k )
      {
        for ( i = 0; i < k + 2; i++ )
       { 
         y[i] = h[i*restart+k];
        }
printf("\n\n**** I am assigning values of h to y*****\n\n");
print_vector("\n\n vector Y is\n\n",y ,m);
    
       for ( j = 0; j < k; j++ ) 
       {
          mult_givens ( cs[j], sn[j], j, y );
       }

 printf("\n\n**** Checking the multiplcation of multiple givens routine\n\n");
 print_vector("\n\n vector Y after multiplication is\n\n",y ,m);


         for ( i = 0; i < k + 2; i++ ) 
       {
          h[i*restart+k] = y[i];
       }
     }

/**********************************************************************
 TO DO:1 Replace the following line of codes with routine
       2 Understand the working of givens matrix in updating the results
       3 Implement backward substitution 
************************************************************************/

     temp = sqrt ( h[k*restart+k] * h[k*restart+k] + h[(k+1)*restart+k] * h[(k+1)*restart+k] );
     h1 = h[k*restart+k];
     h2 = h[(k+1)*restart+k];

      printf("\n\nvalue of temp/r at loop %d is %f\n",k+1, temp);

      printf("\n\nvalue of x at loop %d is %f\n", k+1,h1);
       printf("\n\nvalue of y at loop %d is %f\n", k+1,h2);


      cs[k] = h[k*restart+k] / temp;
      printf("\n\nvalue of c at loop %d is %f\n", k+1,cs[k]);
      
      sn[k] = -h[(k+1)*restart+k] / temp;
      printf("\n\nvalue of s at loop %d is %f\n", k+1,sn[k]);


  //    printf("\n\n\nMatrix after applying  Rotation at loop  %d is :\n", k+1);
    //  print_matrix(h,m,restart);
  
      h[k*restart+k] = cs[k] * h[k*restart+k] - sn[k] * h[(k+1)*restart+k];
//updted value of h[k*restart+k]
      h[(k+1)*restart+k] = 0.0;
      mult_givens ( cs[k], sn[k], k, g);      // this e should be the givens matrix 

      //printf("\n\n\nMatrix after update at last step is :\n");
    //  print_matrix(h,m,restart);

       printf("\n****************************************\n");
      printf("\nEnd of loop\n");
      printf("\n****************************************\n");
      
   }////End of outer loop 


//****************************
//Printing values for DEBUGGING
//***************************//

printf("\n\n\nMatrix after Rotation is :\n");
print_matrix(h,m,restart);

print_vector("\n\nC vector is\n\n",cs,m);

print_vector("\n\nS vector is\n\n",sn,m);


printf("\n\n\n Givens Matrix is :\n");
print_matrix(g,m,m);
}// End of main function 


//______________________________________________________________________________________
//User Defined Functions
//______________________________________________________________________________________

void GeneratePlaneRotation(double dx, double dy, double cs, double sn)
{
  double temp;	
  if (dy == 0.0) {
    cs = 1.0;
    sn = 0.0;
  } else if (fabs(dy) > fabs(dx)) {
    temp = dx / dy;
    sn = 1.0 / sqrt( 1.0 + temp*temp );
    cs = temp * sn;
  } else {
    temp = dy / dx;
    cs = 1.0 / sqrt( 1.0 + temp*temp );
    sn = temp * cs;
  }
}


//template<class Real>
void mult_givens(double cs, double sn, int k, double *g)
{
	/*Input, double C, S, the cosine and sine of a Givens
    rotation.

    Input, int K, indicates the location of the first vector entry.

    Input/output, double G[K+2], the vector to be modified.  On output,
    the Givens rotation has been applied to entries G(K) and G(K+1).
*/

        double g1;
	double g2;
         g1  =  cs * g[k] - sn * g[k+1];
         g2  =  sn * g[k]+ cs * g[k+1];
          
	 g[k]=g1;
	 g[k+1]= g2;
	 return;
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

void print_matrix(double *arr, int rows, int cols){

     for(int i = 0; i <rows; i++){
         for (int j=0;j<cols; j++){
                printf("\t%e\t",arr[i*cols+j]);

        }
       printf("\n");
     } }
