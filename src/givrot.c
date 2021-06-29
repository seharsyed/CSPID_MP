//This is the file made ito test the implementation of givens rotation on a small hessenberg matrix
//I will read a small H matrix then apply the givens rotation on it. Check the same result on Matlab. 

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <string.h>




void GeneratePlaneRotation(double dx, double dy, double cs, double sn);
//void ApplyPlaneRotation(double dx, double dy, double cs, double sn);

int main ()
{

double *cs;
double *sn;
double *h;
double *y; // starting as a vector and will implement it later for matrix
double *g;

int m; //maximum no of inner iterations, mr < size of matrix 
//mr = m,
int restart;

int i;
int j;
int k;
double temp;

    printf("Enter number of Rows :");
    scanf("%d",&m);
    printf("Enter number of Cols :");
    scanf("%d",&restart);

cs = ( double * ) malloc ( m * sizeof ( double ) );
sn = ( double * ) malloc ( m * sizeof ( double ) );
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

//Declaration of rotation matrix I am starting with identity and then discover the real 




    for (k = 0; k < m; k++)
    g[k*restart+k] = 1.0;

//Starting the rotation here! 
for (i = 0; i < m; i++){   
        for (k = 0;k<i; k++){	
  
	      	temp = cs[k]*h[k*restart+i]+sn[k]*h[(k+1)*restart+i];
   h[(k+1)*restart+i]= -sn[k]*h[k*restart+i]+cs[k]*h[(k+1)*restart+i];
   h[k*restart+i]= temp;
	}
    
 
	    
	//	ApplyPlaneRotation(h[j*restart+i], h[(j+1)*restart+i], cs[j], sn[j]);
      
       GeneratePlaneRotation(h[i*restart+i], h[(i+1)*restart+i], cs[i], sn[i]);
      // ApplyPlaneRotation(h[i*restart+i], h[(i+1)*restart+i],cs[i], sn[i]);
    //   ApplyPlaneRotation(g[i], g[i+1], cs[i], sn[i]);
      
   }//End of outer loop 




printf("\n\n\nMatrix after Rotation is :\n");
    for(i=0;i< m;i++)
    {
        for(j=0;j< restart;j++)


        {
            printf("%e\t",h[i*restart+j]);
        }
        printf("\n");   /*new line after row elements*/
    }


printf("\n\n\nMatrix g after Rotation is :\n");
    for(i=0;i< m;i++)
    {
        for(j=0;j< restart;j++)


        {
            printf("%e\t",g[i*restart+j]);
        }
        printf("\n");   /*new line after row elements*/
    }




}// End of main function 

//template<class Real>
void GeneratePlaneRotation(double dx, double dy, double cs, double sn)
{
  double temp;	
  if (dy == 0.0) {
    cs = 1.0;
    sn = 0.0;
  } else if (abs(dy) > abs(dx)) {
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
void ApplyPlaneRotation(double dx, double dy, double cs, double sn)
{
        double temp;
	temp  =  cs * dx + sn * dy;
         dy = -sn * dx + cs * dy;
          
	 dx = temp;
}

