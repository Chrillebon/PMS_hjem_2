#include <stdio.h>
#include "matrix_io.h"


int main(void) {    

    /* Read matrix and vector from txt file */
    matrix_t * A = read_matrix("A.txt");
    vector_t * b = read_vector("b.txt");

    print_matrix(A);
    printf("*x = \n");
    print_vector(b);

    int info = call_dgesv(A,b);

    /* Present solution */
    if (info == 0)
    {
        printf("Solution is: x = \n");
        print_vector(b);
    } else {
        printf("No solution was found!\n");
    }
    
    return 0;
}



