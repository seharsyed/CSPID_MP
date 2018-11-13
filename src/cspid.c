/************************************
  This is going to be the main file
 *************************************
 Author:  Sehar Naveed
 Institute:  Unibz
 Date: November 12, 2018
***************************************/
/*Implementation of Block Arnoldi*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include "mmio.h"



/* Reading the Matrix from MM. The matrix is in COO Format*/


int main(int argc, char *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;
    int i, *I, *J, k, j;
    double *val;

    if (argc < 2)
	{
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	}
    else
    {
        if ((f = fopen(argv[1], "r")) == NULL)
            exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);


    /* reseve memory for matrices */

 if (mm_is_symmetric(matcode)){

    I = (int *) malloc(nz *2* sizeof(int));
    J = (int *) malloc(nz *2* sizeof(int));
    val = (double *) malloc(nz *2* sizeof(double));
 }

 else {
    I = (int *) malloc(nz * sizeof(int));
    J = (int *) malloc(nz * sizeof(int));
    val = (double *) malloc(nz * sizeof(double));
  }

    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */


 k=0;
  for (i=0; i<nz; i++)  {
    if (mm_is_pattern(matcode)){
      fscanf(f, "%d %d", &I[i], &J[i]);
      I[i] --;  /* adjust from 1-based to 0-based */
      J[i] --;

      val[i] = random_double(-1, 1);
    }
    else if (mm_is_real(matcode)){
      fscanf(f, "%d %d %le\n", &I[i], &J[i], &val[i]);
      (I)[i] --;  /* adjust from 1-based to 0-based */
      (J)[i] --;
    }

    if (mm_is_symmetric(matcode)){
      if ( I[i] != J[i] ){
	I[nz+k] = J[i];
	J[nz+k] = I[i];
	val[nz+k] = val[i];
	k++;
      }
    }
  }
  nz += k;

   if (f !=stdin) fclose(f);

    coo2csr_in (M, nz, val, I, J);

    /************************/
    /* now write out matrix */
    /************************/

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    for (i=0; i<M; i++){
    
    
         fprintf(stdout, "%d %d\n",i, I[i+1]-I[i]);
            for (j=I[i]; j<I[i+1]; j++)
                fprintf(f, "%d %2.8le\n", J[j], val[j]);
    
    }

	return 0;
}

/*********************************
 * TODO 
*********************************/

/*1. Matrix vector multiplication of sparse matrix in COO format with a vector
  2. Creat a Random matrix B in Matlab/C and allocate the values to array. 
  3. Compute Norm */ 
