#include <stdio.h>
#include "matrix_io.h"


int main(void) {    

    matrix_t * A = read_matrix("A.txt");
    vector_t * b = read_vector("b.txt");

    vector_t * x = b;

    int info = call_dgesv(A,b);

    if (info == 0)
    {
        print_matrix(A);
        printf("*");
        print_vector(x);
        printf("=");
        print_vector(b);
    } else {
        printf("No solution was found!\n");
    }

    return 0;
}



