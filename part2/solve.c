#include <stdlib.h>
#include <stdio.h>
#include "matrix_io.h"
#define NULL_INPUT -9
#define NON_SQUARE_MATRIX -10
#define INCOMPATIBLE_DIMENSIONS -11
#define MEMORORY_ALLOCATION_ERROR -12

/* C prototype for LAPACK routine DGESV */
void dgesv_(const int *n,    /* columns/rows in A          */
            const int *nrhs, /* number of right-hand sides */
            double *A,       /* array A                    */
            const int *lda,  /* leading dimension of A     */
            int *ipiv,       /* pivoting array             */
            double *B,       /* array B                    */
            const int *ldb,  /* leading dimension of B     */
            int *info        /* status code                */
            );


/* local_call_dgesv : wrapper for LAPACK's DGESV routine


Purpose:
Solves system of equations A*x=b using LAPACK's DGESV routine
Upon exit, the input vector b is overwriten by the solution x.

Return value:
The function returns the output `info` from DGESV with the
following exceptions: the return value is

   -9 if the input A is NULL and/or the input B is NULL
   -10 if A is not a square matrix
   -11 if the dimensions of A and b are incompatible
   -12 in case of memory allocation errors.
*/
int local_call_dgesv(matrix_t * A, vector_t * b) {
  // Testing if input is NULL
  if (A == NULL || b == NULL)
  {
    return NULL_INPUT;
  }
  // Testing if A or v are NULL
  if (A->A == NULL || b->v == NULL)
  {
    return NULL_INPUT;
  }
  // Testing if A[0] is NULL
  if (A->A[0] == NULL)
  {
    return NULL_INPUT;
  }
  // Testing for illegal dimensions
  if (A->m <= 0 || A->n <= 0 || b->n <= 0)
  {
    return MATRIX_IO_FAILURE;
  }
  // Testing for non-square matrix
  if (A->m != A->n)
  {
    return NON_SQUARE_MATRIX;
  }
  // testing for incompatible dimensions.
  if (A->n != b->n)
  {
    return INCOMPATIBLE_DIMENSIONS;
  }

  // Creating local variables
  int n = A->n, LDB = n, nrhs = 1, info = 0;
  int * IPIV = malloc(n*sizeof(int));
  if(IPIV == NULL)
  {
    return MEMORORY_ALLOCATION_ERROR;
  }

  // Creating empty matrix for argument in dgesv_
  double ** temp = malloc(n*sizeof(double*));
  if (temp == NULL)
  {
    return MEMORORY_ALLOCATION_ERROR;
  }
  temp[0] = malloc(n*n*sizeof(double));
  if (temp[0] == NULL)
  {
    free(temp);
    return MEMORORY_ALLOCATION_ERROR;
  }
  for (size_t i = 1; i < n; i++) {
    temp[i] = temp[0] + i*n;
  }

  // Transposing A->A into temp
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      temp[i][j] = A->A[j][i];
    }
  }

  // Calling function
  dgesv_(&n, &nrhs, *temp, &n, IPIV, b->v, &LDB, &info);

  // Trash collection
  free(IPIV);
  free(temp[0]);
  free(temp);

  return info;
}


int main(int argc, char *argv[]) {

  // 4 arguments needed, if not:
  if (argc != 4) {
    fprintf(stderr,"Usage: %s A b x\n", argv[0]);
    return EXIT_FAILURE;
  }

  // Opening matrix A
  matrix_t * A = read_matrix(argv[1]);
  // Error handling
  if(A == NULL)
  {
    return EXIT_FAILURE;
  }

  // Opening vector b
  vector_t * b = read_vector(argv[2]);
  // Error handling
  if(b == NULL)
  {
    free_matrix(A);
    return EXIT_FAILURE;
  }

  // Solve the system
  int returnval;
  returnval = local_call_dgesv(A,b);

  // Error handling
  if(returnval != 0)
  {
    if(returnval == NULL_INPUT)
    {
      free_matrix(A);
      free_vector(b);
      fprintf(stderr, "Testing for NULL before sending into function failed!\n");
      return EXIT_FAILURE;
    }
    if(returnval == NON_SQUARE_MATRIX)
    {
      free_matrix(A);
      free_vector(b);
      fprintf(stderr, "The matrix found in %s was non-square!\n", argv[1]);
      return EXIT_FAILURE;
    }
    if(returnval == INCOMPATIBLE_DIMENSIONS)
    {
      free_matrix(A);
      free_vector(b);
      fprintf(stderr, "The matrix from %s, and the vector from %s had incompatible dimensions!\n", argv[1], argv[2]);
      return EXIT_FAILURE;
    }
    if(returnval == MEMORORY_ALLOCATION_ERROR)
    {
      free_matrix(A);
      free_vector(b);
      fprintf(stderr, "Memory allocation failed trying to solve the system Ax=b\n");
      return EXIT_FAILURE;
    }
    // singular system (other failure)
    free_matrix(A);
    free_vector(b);
    fprintf(stderr, "Failed to solve the system Ax=b\n");
    return EXIT_FAILURE;
  }
  // Computing solution succeeded

  // Writing solution
  returnval = write_vector(argv[3],b);
  // Error handling
  if(returnval != MATRIX_IO_SUCCESS)
  {
    free_matrix(A);
    free_vector(b);
    return EXIT_FAILURE;
  }

  // Trash collection
  free_matrix(A);
  free_vector(b);

  return EXIT_SUCCESS;
}
