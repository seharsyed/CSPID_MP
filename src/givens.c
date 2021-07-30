//This is the file made ito test the ementation of givens rotation on a small hessenberg matrix
 //I will read a small H matrix then apply the givens rotation on it. Check the same result on Matlab. 
 //This ROT1.C is implementation of C code into givens 
 //
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
 double r;
double newrow;
double newcol;

     printf("Enter number of Rows :");
     scanf("%d",&m);
     printf("Enter number of Cols :");
     scanf("%d",&restart);

 cs = ( double * ) malloc ( m * sizeof ( double ) );
 sn = ( double * ) malloc ( m * sizeof ( double ) );
// g = ( double * ) malloc ( ( m + 1 ) * sizeof ( double ) );//rotation matrix
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
 for (i = 0;i<m-1;i++){

       r  = sqrt ( h[i*restart+i] * h[i*restart+i] + h[(i+1)*restart+i] * h[(i+1)*restart+i] );
       cs[i] = h[i*restart+i] /r;
       sn[i] = -h[(i+1)*restart+k]/r;

 for (j = 0;j<m;j++){ 
       newrow  = cs[j] * h[i*restart+j] - sn[j] * h[(i+1)*restart+j];
       h[(i+1)*restart+j] = sn[j] * h[i*restart+j] + cs[j] * h[(i+1)*restart+j];
       h[i*restart+j] = newrow;

       newcol = cs[j] * g[j*restart+i] - sn[j] * g[j*restart+(i+1)];
       g[j*restart+(i+1)]= sn[j] * g[j*restart+i] + cs[j] * g[j*restart+(i+1)];
       g[j*restart+i]= newcol;
       }

    }//End of outer loop 


 //****************************
 //Printing values for DEBUGGING
 //***************************//

 printf("\n\n\nMatrix after Rotation is :\n");
 print_matrix(h,m,restart);

 print_vector("\n\nC vector is\n\n",cs,m);

 print_vector("\n\nS vector is\n\n",sn,m);






 }// End of main function 


 //______________________________________________________________________________________
 //User Defined Functions
 //________________________________iiiii______________________________________________________

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
 	 g[k+2]= g2;
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