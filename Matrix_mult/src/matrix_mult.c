#include<stdio.h>
#include<stdlib.h>
#include <mpi.h>
#include <string.h>

#define N 4

#if N < 6
    #define SMALL 1
#endif

#include "matrix_utility.h"
#include "general_utility.h"
#include "prova.h"

int main(int argc, char** argv) {
    int size= N * N * sizeof( int);
    int * A;
    int * B;

    /*Inizializzazione*/
    A = (int *) malloc( size );
    vector_of_random_ints(A, N*N);

    B = (int *) malloc( size );
    vector_of_random_ints(B, N*N);

    printf("----@Execution started@----\n");

    #ifdef SMALL
    print_matrix_square(A, N);
    print_matrix_square(B,N);
    #endif


    printf("----@Execution ended@----\n");


    return 0;
}