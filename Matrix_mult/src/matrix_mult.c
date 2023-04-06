#include<stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>

#include "linearized_matrix_utility.h"
#include "general_utility.h"
#include "prova.h"

#define N 6

#if N < 6
    #define SMALL 1
#endif

int main(int argc, char** argv) {
    
    int size= N * N * sizeof( int);
    double *A, *B, *C;

    int n_proc, irank;

    MPI_Init ( & argc , & argv ) ;
    MPI_Comm_rank ( MPI_COMM_WORLD , & irank ) ;
    MPI_Comm_size ( MPI_COMM_WORLD , & n_proc ) ;

    //A = (double *) malloc( size );
    //array_of_random_doubles(A, N*N);

    //B = (double *) malloc( size );
    //create_identity_matrix(B, N);

    //C = (double *) malloc( size );
    //create_null_matrix(C, N*N);

    /*printf("----@I am %d @----\n", irank );
    #ifdef SMALL
    print_matrix_square(A,N);
    print_matrix_square(B,N);
    #endif



    #ifdef SMALL
    //print_matrix_square(C,N);
    #endif

    printf("----@Execution ended@----\n");*/
    MPI_Finalize();

    return 0;
}