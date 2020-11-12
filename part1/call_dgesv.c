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

/* call_dgesv : wrapper for LAPACK's DGESV routine

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
int call_dgesv(matrix_t * A, vector_t * b) {

  if (A == NULL || b == NULL)
  {
    return NULL_INPUT;
  }
  
  if (A->A == NULL || b->v == NULL)
  {
    return NULL_INPUT;
  }
  
  if (A->A[0] == NULL)
  {
    return NULL_INPUT;
  }
  
  if (A->m <= 0 || A->n <= 0 || b->n <= 0)
  {
    return MATRIX_IO_FAILURE;
  }
  
  if (A->m != A->n)
  {
    return NON_SQUARE_MATRIX;
  }
  
  if (A->n != b->n)
  {
    return INCOMPATIBLE_DIMENSIONS;
  }
  
  int n = A->n, LDB = 1, nrhs = 1, info = 0;
  int * IPIV = malloc(n*sizeof(int));

  if (IPIV == NULL) 
  {
    return MEMORORY_ALLOCATION_ERROR;
  }
  
  dgesv_(&n, &nrhs, *(A->A), &n, IPIV, b->v, &LDB, &info);  
  
  free_matrix(IPIV);

  return info;
}
