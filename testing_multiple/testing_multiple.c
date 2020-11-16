#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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

int vectorcompare(vector_t * ans, vector_t * b, char * output)
{
  if(ans->n != b->n)
  {
    output = "Dimension mismatch";
    return -1;
  }
  double epsilon = 1e-15;
  for(int i=0;i<ans->n;i++)
  {
    if(ans->v[i]-b->v[i] >= epsilon)
    {
      sprintf(output, "Absolute discreptency detected at line %d\n", i);
      return -2;
    }
    /*if((ans->v[i]-b->v[i])/ >= epsilon)
    {
      fprintf(output, "Absolute discreptency detected at line %d\n", i);
      return -2;
    }*/
  }
  return 0;
}


int main() {

  char * A_test = "Ai.txt";
  char * b_test = "bi.txt";
  char * x_test = "xi.txt";

  // Update this for every new testcase
  int testcases = 1;

  if(testcases >= 10)
  {
    fprintf(stderr, "Please don't have more than 9 testcases :)\n");
    return 0;
  }

  // init
  matrix_t * A;
  vector_t * b;
  vector_t * x_res;     // For later checking of correctness
  int returnval;
  FILE * x_print;

  // Running all tesetcases
  for(int i=0;i<testcases;i++)
  {
    // Changing filenames
    itoa(i, &A_test[1], 1);
    itoa(i, &b_test[1], 1);
    itoa(i, &x_test[1], 1);

    // Opening matrix A
    A = read_matrix(A_test);
    // Error handling
    if(A == NULL)
    {
      x_print = fopen(x_test, "w");
      fprintf(x_print, "ERROR A_read");
      fclose(x_print);
      continue;
    }

    // Opening vector b
    b = read_vector(b_test);
    // Error handling
    if(b == NULL)
    {
      free_matrix(A);
      x_print = fopen(x_test, "w");
      fprintf(x_print, "ERROR b_read");
      fclose(x_print);
      continue;
    }

    // Solve the system
    returnval = local_call_dgesv(A,b);

    // Error handling
    if(returnval != 0)
    {
      free_matrix(A);
      free_vector(b);
      if(returnval == NULL_INPUT)
      {
        x_print = fopen(x_test, "w");
        fprintf(x_print, "ERROR send_null_%d", returnval);
        fclose(x_print);
        continue;
      }
      if(returnval == NON_SQUARE_MATRIX)
      {
        x_print = fopen(x_test, "w");
        fprintf(x_print, "ERROR matrix_nonsquare_%d", returnval);
        fclose(x_print);
        continue;
      }
      if(returnval == INCOMPATIBLE_DIMENSIONS)
      {
        x_print = fopen(x_test, "w");
        fprintf(x_print, "ERROR incomp_dim_%d", returnval);
        fclose(x_print);
        continue;
      }
      if(returnval == MEMORORY_ALLOCATION_ERROR)
      {
        x_print = fopen(x_test, "w");
        fprintf(x_print, "ERROR malloc_failed_%d", returnval);
        fclose(x_print);
        continue;
      }
      // singular system (other failure)
      x_print = fopen(x_test, "w");
      fprintf(x_print, "ERROR solving_failed_%d", returnval);
      fclose(x_print);
      continue;
    }
    // Computing solution succeeded

    // Writing solution
    returnval = write_vector(x_test,b);
    // Error handling
    if(returnval != MATRIX_IO_SUCCESS)
    {
      free_matrix(A);
      free_vector(b);
      x_print = fopen(x_test, "w");
      fprintf(x_print, "ERROR writing_solution");
      fclose(x_print);
      continue;
    }

    // Trash collection
    free_matrix(A);
    free_vector(b);
  }

  // Testing that it works:
  char * ans = "xi_.txt";
  char errormessage[100] = "NONE";
  char erroroutput[100] = "NONE";
  for(int i=0;i<testcases;i++)
  {
    itoa(i, &ans[1], 1);
    itoa(i, &b_test[1], 1);

    // testing for errors:
    x_print = fopen(ans, "r");
    fgets(errormessage, 101, x_print);
    fclose(x_print);
    x_print = fopen(b_test, "r");
    fgets(erroroutput, 101, x_print);
    fclose(x_print);

    if(errormessage[0] == 'E' || erroroutput[0] == 'E')
    {
      fprintf(stderr, "Error was detected in testcase %d, this is the output from both files\nAnswer:\n%s\n\nProgram output:\n%s\n", i, errormessage, erroroutput);
    }
    else
    {
      b = read_vector(b_test);
      x_res = read_vector(ans);

      returnval = vectorcompare(x_res, b, errormessage);
      if(returnval != 0)
      {
        continue;
      }
      else
      {
        fprintf(stderr, "Testcase %d succeeded!\n", i);
      }
    }
  }

  return EXIT_SUCCESS;
}
