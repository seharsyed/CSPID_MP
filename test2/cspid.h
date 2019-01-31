#ifndef C_SPID
#define C_SPID

#include <sys/time.h>


/*********************
*Clock
*********************/

/**
 * This structure is used to measure and retrieve
 * execution time.
 */
typedef struct {
	struct timeval begin; // Stores initial time.
	struct timeval end;   // Stores ending time.
	double sec;           // Stores measured time in seconds.
} Clock;

void clock_start(Clock *clock);
void clock_stop(Clock *clock);

/**
 * Begins measuring time.
 */
void clock_start(Clock *clock){
	gettimeofday(&(clock->begin),NULL);
}

/**
 * Stops measuring time.
 * Time measure is store in clock->sec,
 * in seconds.
 */
void clock_stop(Clock *clock){
	gettimeofday(&(clock->end),NULL);
	
	clock->sec = (clock->end).tv_sec;
	clock->sec -= (clock->begin).tv_sec;
	clock->sec *= 1000000.0;
	clock->sec += (clock->end).tv_usec;
	clock->sec -= (clock->begin).tv_usec;
	//clock->sec /= 1000000.0;
} //END of CLOCK 

/************************
*Matrix Market Definitions
*************************/

#define MM_MAX_LINE_LENGTH 1025
#define MatrixMarketBanner "%%MatrixMarket"
#define MM_MAX_TOKEN_LENGTH 64

typedef char MM_typecode[4];

char *mm_typecode_to_str(MM_typecode matcode);

int mm_read_banner(FILE *f, MM_typecode *matcode);
int mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz);
int mm_read_mtx_array_size(FILE *f, int *M, int *N);

int mm_write_banner(FILE *f, MM_typecode matcode);
int mm_write_mtx_crd_size(FILE *f, int M, int N, int nz);
int mm_write_mtx_array_size(FILE *f, int M, int N);


/********************* MM_typecode query fucntions ***************************/

#define mm_is_matrix(typecode)	((typecode)[0]=='M')

#define mm_is_sparse(typecode)	((typecode)[1]=='C')
#define mm_is_coordinate(typecode)((typecode)[1]=='C')
#define mm_is_dense(typecode)	((typecode)[1]=='A')
#define mm_is_array(typecode)	((typecode)[1]=='A')

#define mm_is_complex(typecode)	((typecode)[2]=='C')
#define mm_is_real(typecode)		((typecode)[2]=='R')
#define mm_is_pattern(typecode)	((typecode)[2]=='P')
#define mm_is_integer(typecode) ((typecode)[2]=='I')

#define mm_is_symmetric(typecode)((typecode)[3]=='S')
#define mm_is_general(typecode)	((typecode)[3]=='G')
#define mm_is_skew(typecode)	((typecode)[3]=='K')
#define mm_is_hermitian(typecode)((typecode)[3]=='H')

int mm_is_valid(MM_typecode matcode);		/* too complex for a macro */


/********************* MM_typecode modify fucntions ***************************/

#define mm_set_matrix(typecode)	((*typecode)[0]='M')
#define mm_set_coordinate(typecode)	((*typecode)[1]='C')
#define mm_set_array(typecode)	((*typecode)[1]='A')
#define mm_set_dense(typecode)	mm_set_array(typecode)
#define mm_set_sparse(typecode)	mm_set_coordinate(typecode)

#define mm_set_complex(typecode)((*typecode)[2]='C')
#define mm_set_real(typecode)	((*typecode)[2]='R')
#define mm_set_pattern(typecode)((*typecode)[2]='P')
#define mm_set_integer(typecode)((*typecode)[2]='I')


#define mm_set_symmetric(typecode)((*typecode)[3]='S')
#define mm_set_general(typecode)((*typecode)[3]='G')
#define mm_set_skew(typecode)	((*typecode)[3]='K')
#define mm_set_hermitian(typecode)((*typecode)[3]='H')

#define mm_clear_typecode(typecode) ((*typecode)[0]=(*typecode)[1]= \
									(*typecode)[2]=' ',(*typecode)[3]='G')

#define mm_initialize_typecode(typecode) mm_clear_typecode(typecode)


/********************* Matrix Market error codes ***************************/


#define MM_COULD_NOT_READ_FILE	11
#define MM_PREMATURE_EOF		12
#define MM_NOT_MTX				13
#define MM_NO_HEADER			14
#define MM_UNSUPPORTED_TYPE		15
#define MM_LINE_TOO_LONG		16
#define MM_COULD_NOT_WRITE_FILE	17


/******************** Matrix Market internal definitions ********************
   MM_matrix_typecode: 4-character sequence
				    ojbect 		sparse/   	data        storage 
						  		dense     	type        scheme
   string position:	 [0]        [1]			[2]         [3]
   Matrix typecode:  M(atrix)  C(oord)		R(eal)   	G(eneral)
						        A(array)	C(omplex)   H(ermitian)
											P(attern)   S(ymmetric)
								    		I(nteger)	K(kew)
 ***********************************************************************/

#define MM_MTX_STR		"matrix"
#define MM_ARRAY_STR	"array"
#define MM_DENSE_STR	"array"
#define MM_COORDINATE_STR "coordinate" 
#define MM_SPARSE_STR	"coordinate"
#define MM_COMPLEX_STR	"complex"
#define MM_REAL_STR		"real"
#define MM_INT_STR		"integer"
#define MM_GENERAL_STR  "general"
#define MM_SYMM_STR		"symmetric"
#define MM_HERM_STR		"hermitian"
#define MM_SKEW_STR		"skew-symmetric"
#define MM_PATTERN_STR  "pattern"


/*  high level routines */

int mm_write_mtx_crd(char fname[], int M, int N, int nz, int II[], int JJ[],
		 double val[], MM_typecode matcode);
int mm_read_mtx_crd_data(FILE *f, int M, int N, int nz, int II[], int JJ[],
		double val[], MM_typecode matcode);
int mm_read_mtx_crd_entry(FILE *f, int *II, int *JJ, double *real, double *img,
			MM_typecode matcode);

int mm_read_unsymmetric_sparse(const char *fname, int *M_, int *N_, int *nz_,
                double **val_, int **I_, int **J_);


/********************
*COO ROUTINES 
********************/


/**
 * Sparse matrix structure stored
 * in coordinate list (COO) format.
 */
typedef struct {
	double *val;       // Value of each non-zero element.
	unsigned int *row; // Corresponding row coordinate.
	unsigned int *col; // Corresponding column coordinate.
	unsigned int nz;   // Total non-zero elements.
	unsigned int rows; // Total rows of matrix.
	unsigned int cols; // Total columns of matrix.
} COO_Matrix;

void coo_init_matrix(COO_Matrix *m);
void coo_free_matrix(COO_Matrix *m);
void coo_print_matrix(COO_Matrix *m);
int coo_load_matrix(char* filename, COO_Matrix *coo);
void coo_copy(COO_Matrix *in, COO_Matrix *out);
void coo_reorder_by_rows(COO_Matrix *m);

void cooBottomUpMerge(
	unsigned int *A, unsigned int *A2, double *A3, 
	unsigned int iLeft, unsigned int iRight, unsigned int iEnd, 
	unsigned int *B, unsigned int *B2, double *B3 );
void cooCopyArray(unsigned int *B, unsigned int *A, unsigned int n);
void cooCopyArrayDouble(double *B, double *A, unsigned int n);
unsigned int cooMergeMin(unsigned int a, unsigned int b);
void cooMergeSort(unsigned int *A, unsigned int *A2, double *A3, unsigned int n);

/**
 * Initialize matrix m.
 */
void coo_init_matrix(COO_Matrix *m){
	m->val = NULL;
	m->row = NULL;
	m->col = NULL;
	m->nz = m->rows = m->cols = 0;
}

/**
 * Free resources associated with matrix m.
 */
void coo_free_matrix(COO_Matrix *m){
	if(m->val != NULL) free(m->val);
	if(m->row != NULL) free(m->row);
	if(m->col != NULL) free(m->col);
	m->val = NULL;
	m->row = NULL;
	m->col = NULL;
	m->nz = m->rows = m->cols = 0;
}

/**
 * Prints matrix m data on console.
 */
void coo_print_matrix(COO_Matrix *m){
	unsigned int i;
	printf("val= ");
	for(i=0; i < m->nz; i++){
		printf("%.2g ", m->val[i]);
	}
	printf("\nrow= ");
	for(i=0; i < m->nz; i++){
		printf("%d ", m->row[i]);
	}
	printf("\ncol= ");
	for(i=0; i < m->nz; i++){
		printf("%d ", m->col[i]);
	}
	printf("\nnz= %d\n", m->nz);
	printf("rows= %d\n", m->rows);
	printf("cols= %d\n", m->cols);
}

/**
 * Load COO matrix from file filename.
 * File must be in Matrix Market format.
 */
int coo_load_matrix(char* filename, COO_Matrix *coo){
	FILE *file;
	MM_typecode code;
	int m,n,nz,i,ival;
	unsigned int dc = 0;
	
	// Open file
	file = fopen(filename, "r");
	if(file == NULL){
		fprintf(stderr, "ERROR: can't open file \"%s\"\n", filename);
		exit(EXIT_FAILURE);
	}
	
	if (mm_read_banner(file, &code) != 0){
		fprintf(stderr, "ERROR: Could not process Matrix Market banner.\n");
		exit(EXIT_FAILURE);
	}
	
	// Check if input matrix is supported
	if(!mm_is_matrix(code) || !mm_is_sparse(code) || !mm_is_coordinate(code) ){
		fprintf(stderr, "Sorry, this application does not support ");
		fprintf(stderr, "Market Market type: [%s]\n", mm_typecode_to_str(code));
		exit(EXIT_FAILURE);
	}
	
	// Read Matrix size
	if ((mm_read_mtx_crd_size(file, &m, &n, &nz)) != 0){
		fprintf(stderr, "ERROR: Could not read matrix size.\n");
		exit(EXIT_FAILURE);
	}
	coo->rows = m;
	coo->cols = n;
	coo->nz = nz;
	
	// Reserve space for matrix data
	coo->val = (double *) malloc(nz * sizeof(double));
	coo->row = (unsigned int *) malloc(nz * sizeof(unsigned int));
	coo->col = (unsigned int *) malloc(nz * sizeof(unsigned int));
	
	// Read matrix
	for(i = 0; i < nz; i++){
		if(mm_is_pattern(code)){
			if(fscanf(file, "%d %d\n", &(coo->row[i]), &(coo->col[i]) ) < 2){
				fprintf(stderr, "ERROR: reading matrix data.\n");
				exit(EXIT_FAILURE);
			}
			coo->val[i] = 1.0;
		}else if (mm_is_real(code)){
			if(fscanf(file, "%d %d %lg\n", &(coo->row[i]), &(coo->col[i]), &(coo->val[i])) < 3){
				fprintf(stderr, "ERROR: reading matrix data.\n");
				exit(EXIT_FAILURE);
			}
		}else if (mm_is_integer(code)){
			if(fscanf(file, "%d %d %d\n", &(coo->row[i]), &(coo->col[i]), &ival) < 3){
				fprintf(stderr, "ERROR: reading matrix data.\n");
				exit(EXIT_FAILURE);
			}
			coo->val[i] = (double)ival;
		}
		coo->row[i]--;
		coo->col[i]--;
		if(coo->row[i] == coo->col[i]) { ++dc; }
	}
	
	fclose(file);
	
	return mm_is_symmetric(code);
}

void coo_copy(COO_Matrix *in, COO_Matrix *out){
	unsigned int i;
	out->val = (double *) malloc(in->nz * sizeof(double));
	out->row = (unsigned int *) malloc(in->nz * sizeof(unsigned int));
	out->col = (unsigned int *) malloc(in->nz * sizeof(unsigned int));
	out->nz = in->nz;
	out->rows = in->rows;
	out->cols = in->cols;
	
	for(i = 0; i < in->nz; i++){
		out->val[i] = in->val[i]; 
	}
	for(i = 0; i < in->nz; i++){
		out->row[i] = in->row[i]; 
	}
	for(i = 0; i < in->nz; i++){
		out->col[i] = in->col[i]; 
	}
}

void coo_reorder_by_rows(COO_Matrix *m){
	cooMergeSort(m->row, m->col, m->val, m->nz);
}

/**************************************
 * Merge sort utility functions.
 * Code adapted from:
 * http://en.wikipedia.org/wiki/Merge_sort
 *************************************/
void cooBottomUpMerge(
	unsigned int *A, unsigned int *A2, double *A3, 
	unsigned int iLeft, unsigned int iRight, unsigned int iEnd, 
	unsigned int *B, unsigned int *B2, double *B3 ) {
	
	unsigned int i0 = iLeft;
	unsigned int i1 = iRight;
	unsigned int j;
 
	/* While there are elements in the left or right runs */
	for (j = iLeft; j < iEnd; j++){
		/* If left run head exists and is <= existing right run head */
		if (i0 < iRight && (i1 >= iEnd || A[i0] <= A[i1])){
			B[j] = A[i0];
			B2[j] = A2[i0];
			B3[j] = A3[i0];
			i0++;
		} else {
			B[j] = A[i1];
			B2[j] = A2[i1];
			B3[j] = A3[i1];
			i1++;
		}
	}
}
 
void cooCopyArray(unsigned int *B, unsigned int *A, unsigned int n) {
	unsigned int i;
   for(i = 0; i < n; i++){
        A[i] = B[i];
	}
}

void cooCopyArrayDouble(double *B, double *A, unsigned int n) {
	unsigned int i;
   for(i = 0; i < n; i++){
        A[i] = B[i];
	}
}

inline unsigned int cooMergeMin(unsigned int a, unsigned int b){
	return (a<b)?a:b;
}

/**
 * Use MergeSort algorithm to order the elements of A, from high to low.
 * Elements of A2 and A3 are moved correspondingly to elements of A.
 */
void cooMergeSort(unsigned int *A, unsigned int *A2, double *A3, unsigned int n) {
	/* array A[] has the items to sort; array B[] is a work array */
	unsigned int width,i;
	unsigned int *t;
	double *td;
	unsigned int swaps = 0;
	unsigned int *B = (unsigned int *)malloc(n * sizeof(unsigned int));
	unsigned int *B2 = (unsigned int *)malloc(n * sizeof(unsigned int));
	double *B3 = (double *)malloc(n * sizeof(double));
	/* Each 1-element run in A is already "sorted". */
	/* Make successively longer sorted runs of length 2, 4, 8, 16... until whole array is sorted. */
	for (width = 1; width < n; width = 2 * width){
      /* Array A is full of runs of length width. */
      for (i = 0; i < n; i = i + 2 * width){
         /* Merge two runs: A[i:i+width-1] and A[i+width:i+2*width-1] to B[] */
         /* or copy A[i:n-1] to B[] ( if(i+width >= n) ) */
         cooBottomUpMerge(A,A2,A3, i, cooMergeMin(i+width, n), cooMergeMin(i+2*width, n), B,B2,B3);
      }
      /* Now work array B is full of runs of length 2*width. */
      /* Copy array B to array A for next iteration. */
      /* A more efficient implementation would swap the roles of A and B */
      //cooCopyArray(B, A, n);
		//cooCopyArray(B2, A2, n);
		//cooCopyArrayDouble(B3, A3, n);
		t = B; B = A; A = t;
		t = B2; B2 = A2; A2 = t;
		td = B3; B3 = A3; A3 = td;
		swaps++;
      /* Now array A is full of runs of length 2*width. */
   }
	if(swaps%2 == 1){
		//printf("Swap A and B\n");
		t = B; B = A; A = t;
		t = B2; B2 = A2; A2 = t;
		td = B3; B3 = A3; A3 = td;
		cooCopyArray(B, A, n);
		cooCopyArray(B2, A2, n);
		cooCopyArrayDouble(B3, A3, n);
	}
	
	free(B);
	free(B2);
	free(B3);
} //End of COO

/***************************
CSR Routines
***************************/

typedef struct {
	double *val;       // Value of each non-zero element.
	unsigned int *col; // Corresponding column coordinate.
	unsigned int *ptr; // Initial position of each row.
	unsigned int nz;   // Total non-zero elements.
	unsigned int rows; // Total rows of matrix.
	unsigned int cols; // Total columns of matrix.
} CSR_Matrix;


void csr_init_matrix(CSR_Matrix *m);
void csr_free_matrix(CSR_Matrix *m);
void csr_print_matrix(CSR_Matrix *m);
void coo2csr(COO_Matrix *in, CSR_Matrix *out);
void coo22csr(COO_Matrix *in, CSR_Matrix *out);
void coo222csr(COO_Matrix *in, CSR_Matrix *out);
int csr_load_matrix(char* filename, CSR_Matrix *m);


/**
 * Initialize matrix m.
 */
void csr_init_matrix(CSR_Matrix *m){
	m->val = NULL;
	m->col = NULL;
	m->ptr = NULL;
	m->nz = m->rows = m->cols = 0;
}

/**
 * Free resources associated with matrix m.
 */
void csr_free_matrix(CSR_Matrix *m){
	if(m->val != NULL) free(m->val);
	if(m->col != NULL) free(m->col);
	if(m->ptr != NULL) free(m->ptr);
	m->val = NULL;
	m->col = NULL;
	m->ptr = NULL;
	m->nz = m->rows = m->cols = 0;
}

/**
 * Prints matrix m data on console.
 */
void csr_print_matrix(CSR_Matrix *m){
	// unsigned int i,j;
	// printf("val= ");
	// j = 1;
	// for(i=0; i < m->nz; i++){
		// if(i == m->ptr[j]){
			// printf("|");
			// j++;
		// }
		// printf("%.2g ", m->val[i]);
	// }
	// j = 1;
	// printf("\ncol= ");
	// for(i=0; i < m->nz; i++){
		// if(i == m->ptr[j]){
			// printf("|");
			// j++;
		// }
		// printf("%d ", m->col[i]);
	// }
	// printf("\nptr= ");
	// for(i=0; i < m->rows+1; i++){
		// printf("%d ", m->ptr[i]);
	// }
	printf("\tnz= %d\n", m->nz);
	printf("\trows= %d\n", m->rows);
	printf("\tcols= %d\n", m->cols);
}

/**
 * Copies the COO matrix 'in' into CSR matrix 'out'.
 */
void coo2csr(COO_Matrix *in, CSR_Matrix *out){
	unsigned int i,j;
	unsigned int tot = 0;
	
	// Init mem for CSR matrix
	out->val = (double *) malloc(in->nz * sizeof(double));
	out->col = (unsigned int *) malloc(in->nz * sizeof(unsigned int));
	out->ptr = (unsigned int *) malloc(((in->rows)+1) * sizeof(unsigned int));
	out->nz = in->nz;
	out->rows = in->rows;
	out->cols = in->cols;
	
	// Copy elements row by row
	out->ptr[0] = tot;
	for(i = 0; i < out->rows; i++){
		for(j = 0; j < in->nz; j++){
			if(in->row[j] == i){
				out->val[tot] = in->val[j];
				out->col[tot] = in->col[j];
				tot++;
			}
		}
		out->ptr[i+1] = tot;
	}
}

/**
 * Copies the COO matrix 'in' into CSR matrix 'out'.
 */
void coo222csr(COO_Matrix *in, CSR_Matrix *out){
	unsigned int i;
	unsigned int tot = 0;
	COO_Matrix coo;
	
	// Init mem for CSR matrix
	out->val = (double *) malloc(in->nz * sizeof(double));
	out->col = (unsigned int *) malloc(in->nz * sizeof(unsigned int));
	out->ptr = (unsigned int *) malloc(((in->rows)+1) * sizeof(unsigned int));
	out->nz = in->nz;
	out->rows = in->rows;
	out->cols = in->cols;
	
	// Copy in matrix
	coo_copy(in,&coo);
	//coo_print_matrix(&coo);
	coo_reorder_by_rows(&coo);
	//coo_print_matrix(&coo);
	
	// Copy elements row by row
	out->ptr[0] = tot;
	for(i = 0; i < coo.rows; i++){
		while(tot < coo.nz && coo.row[tot] == i){
			out->val[tot] = coo.val[tot];
			out->col[tot] = coo.col[tot];
			tot++;
		}
		out->ptr[i+1] = tot;
	}
	
	coo_free_matrix(&coo);
}

/**
 * Copies the COO matrix 'in' into CSR matrix 'out'.
 * Faster than previous function. (I hope...)
 * TODO: FIX this function doesn't work
 */
void coo22csr(COO_Matrix *in, CSR_Matrix *out){
	unsigned int i,j;
	unsigned int tot;
	unsigned int *col_start = (unsigned int *)malloc(in->cols * sizeof(unsigned int));
	
	// Init mem for CSR matrix
	out->val = (double *) malloc(in->nz * sizeof(double));
	out->col = (unsigned int *) malloc(in->nz * sizeof(unsigned int));
	out->ptr = (unsigned int *) malloc(((in->rows)+1) * sizeof(unsigned int));
	out->nz = in->nz;
	out->rows = in->rows;
	out->cols = in->cols;
	
	// Find column starts.
	//printf("col_start= ");
	j = 0;
	for(i = 0; i < in->cols; i++){
		col_start[i] = j;
		//printf("%d ", col_start[i]);
		while(j < in->nz && in->col[j] == i){
			j++;
		}		
	}
	//printf("\n");
	
	// Copy elements row by row
	tot = 0;
	out->ptr[0] = tot;
	for(i = 0; i < in->rows; i++){
		//printf("For row %d\n",i);
		for(j = 0; j < in->cols; j++){
			//printf("\trow: %d\n",in->row[col_start[j]] );
			if(in->row[col_start[j]] == i){
				// Copy element
				out->val[tot] = in->val[col_start[j]];
				out->col[tot] = in->col[col_start[j]];
				tot++;
				// Advance column start
				col_start[j]++;
			}
		}
		out->ptr[i+1] = tot;
	}
	
	free(col_start);
}

/**
 * Load CSR matrix from file filename.
 * File must be in Matrix Market format.
 */
int csr_load_matrix(char* filename, CSR_Matrix *m){
	int sym;
	COO_Matrix temp;
	coo_init_matrix(&temp);
	sym = coo_load_matrix(filename, &temp);
	coo222csr(&temp,m);
	coo_free_matrix(&temp);
	return sym;

} //END OF CSR 

#endif //END OF HEADER CSPID_H 
