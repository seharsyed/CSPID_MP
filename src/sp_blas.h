/**********************************************************
File for creating Blas SParse Handle and Data from MM   
 
 **********************************************************/
#ifndef SP_BLAS
#define SP_BLAS

#include "coo.h"
#include <blas_sparse.h>

blas_sparse_matrix AA;


int spblas_load_matrix(char* filename, blas_sparse_matrix *m);
void coo2spblas(COO_Matrix*in, blas_sparse_matrix *out);


int spblas_load_matrix(char* filename, blas_sparse_matrix *m){
         int sym;
         COO_Matrix temp;
         coo_init_matrix(&temp);
         sym = coo_load_matrix(filename, &temp);

        // coo222csr(&temp,m); 
          coo2spblas(&temp,&m);

         coo_free_matrix(&temp);
         return sym;
 }
 

void coo2spblas(COO_Matrix *in, blas_sparse_matrix AA){
    unsigned int i;
   //unsigned int tot = 0;

    COO_Matrix coo;

   // Copy in matrix
   coo_copy(in,&coo);

   coo_reorder_by_rows(&coo);

   ///Creating Handle
   AA = BLAS_duscr_begin(coo.rows, coo.rows);

   //Insert entries 
     for (i = 0;i<coo.nz;i++){
     BLAS_duscr_insert_entry(AA, coo.val[i],coo.row[i],coo.col[i]);
     }
    //Comple Construction of sparse Matrix
    BLAS_duscr_end(AA);

    coo_free_matrix(&coo);
}//End of coo2spblas function



#endif // End of Header file 
