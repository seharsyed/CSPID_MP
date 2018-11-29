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
#include "clock.h"
#include "coo.h"
#include "csr.h"

/* Reading the Matrix from MM. The matrix is in CSR Format*/

unsigned int iterations=1000;

void parse_args(int argc, char *argv[]);
void print_vector(char* pre, double *v, unsigned int size);

int main(int argc, char *argv[])
{
    unsigned int i;
    int sym;
    int ret_code;
    CSR_Matrix csr;
    //COO_Matrix coo;//
    //int sym1;//
    Clock clock;
    MM_typecode matcode;
    char* filename;
    double **B;
    double *x, *y;
    FILE *f;       //This file is used for reading RHS//
    int M, N, nz, nrhs;
    int *I, *J, k, j;
    double *val;

parse_args(argc, argv);

 filename = argv[1];

 /***************************************
  * If the matrix Instance is COO
  * ***********************************/

      //Initialize COO Matrix first//
     // coo_init_matrix(&coo);

      //Load COO Matrix
     // printf("Loading matrix \"%s\"\n", filename);
      
     // sym1 = coo_load_matrix(filename, &coo);

     // Print Matrix data
     // printf("COO matrix data:\n");
     // coo_print_matrix(&coo);

/********************************************
 * COO to CSR Matrix - Conversion and loading
 * ********************************************/ 
      // Initialize matrices.
      csr_init_matrix(&csr);
      
     // Load matrix from file into COO.
      printf("Loading matrix \"%s\"\n", filename);
      sym = csr_load_matrix(filename, &csr);
    
     if(sym) printf("Matrix is symmetric\n");
     else    printf("Matrix is general (non-symmetric)\n");
 
      // Print Matrix data
      printf("CSR matrix data:\n");
      csr_print_matrix(&csr);


/*****************************************
 *Random Right Hand Side Generation 
 ****************************************/

    //Random Matrix Generation for RHS//
    /* Here I am reading the matrix B generated in Matlab and allocating the values to C.
    *  I will change it later to create the random matrix in C*/


/*    printf("Enter the no of desired right hand sides\n");
	scanf("%d",&nrhs);*/
/*
    B = (double **)malloc(csr.rows * sizeof(double));
      
     for(i = 0; i < csr.rows; i++){
         for(j=0;j<nrhs;j++){
             B[i][j]=rand();
             printf("%f\t",B[i][j]);
			}
			printf("\n");
          }*/
    nrhs = 10;
    x = (double *)malloc(csr.cols *nrhs* sizeof(double));
	for(i = 0; i < csr.cols; i++){
	  for (j=0;j<nrhs; j++){
              x[i][j] = rand();
	          printf("\t%f\t",x[i][j]);
      }
      printf("\n");
    }

    //print_vector("\t%d\t", x[i][j]);//


exit(EXIT_SUCCESS);
}

/*********************************
 * TODO 
*********************************/
/*
  1. Creat a Random matrix B in Matlab/C and allocate the values to array. 
  2. Compute Norm */

void print_vector(char* pre, double *v, unsigned int size){
       unsigned int i;
      printf("%s", pre);
         for(i = 0; i < size; i++){
         //printf("%.1f ", v[i]);
         printf("%e \n", v[i]);
     }
     printf("\n");
 }

     void parse_args(int argc, char *argv[]){
     int i;
     if (argc < 2) {
         printf("Usage: %s input_matrix [iterations]\n", argv[0]);
         exit(EXIT_FAILURE);
     }
       if(argc >= 3){
         i = atoi(argv[2]);
         if (i <= 0){
             printf("Invalid number of iterations.\n");
             exit(EXIT_FAILURE);
         }
             iterations = i;
     }
 }

