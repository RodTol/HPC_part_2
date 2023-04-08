#include <stdio.h>
#include <mpi.h>

#include "headers/linearized_matrix_utility.h"
#include "headers/general_utility.h"
#include "headers/prova.h"

#define MASTER 0
#define COMM MPI_COMM_WORLD
#define N 6

int main(int argc, char** argv) {
    
    double *A, *B, *C;

    int n_proc_tot = 1, irank = 0;
    int n_loc = N, offset = 0, rest=0;
    int i_local = 0, j_global = 0;

    /*MPI setup*/
    MPI_Init ( & argc , & argv ) ;
    MPI_Comm_rank ( COMM , & irank ) ;
    MPI_Comm_size ( COMM , & n_proc_tot ) ;

    /*Parameters for the matrix distribution*/
    n_loc = N/n_proc_tot;
    rest = N % n_proc_tot;
    offset = 0;
    if (irank == MASTER) { printf("-------------------------------------------\n"
                             "N: %d, n_proc_tot: %d, n_loc: %d, rest: %d \n"
                             "-------------------------------------------\n", N, n_proc_tot, n_loc, rest); }
    
    /*Conditional exit if there's not an even distribution of the matrix*/
    if (rest != 0)
    {
        if (irank == MASTER){
            printf_red();
            printf("WARNING: ");
            printf_reset();
            printf("The matrix cannot be evenly distributed among the processors\n");
        }
        MPI_Barrier(COMM);
        MPI_Finalize();
        _Exit(0);
    }
    
    /*Flags for debugging purpouse*/
    #if n_loc < 4
        #if N < 9
            #define SMALL 1
        #endif
    #endif

    int size= N * n_loc * sizeof( double);
    A = (double *) malloc( size );
    array_of_random_doubles(A, N*n_loc);

    B = (double *) malloc( size );
    array_of_random_doubles(B, N*n_loc);

    C = (double *) malloc( size );
    create_null_array(C, N*n_loc);

    #ifdef SMALL
    MPI_Barrier(COMM);
    printf("----@I am %d @----\n", irank );
    print_matrix(A, n_loc, N);
    printf("\n");
    print_matrix(B, n_loc, N);
    printf("\n");
    #endif

    #ifdef SMALL
    print_matrix(C, n_loc, N);
    printf("\n");
    #endif



    /*If the matrices are small enough (N < 6 and n_proc < 4)
    the program will print the complete form*/
    #ifdef SMALL
        MPI_Barrier(COMM);
        if ( irank == 0 ) {
            printf("\nMatrix A\n");
            print_matrix ( A , n_loc, N) ;
            for (int count = 1; count < n_proc_tot ; count ++ ) {
                MPI_Recv ( A , n_loc * N , MPI_DOUBLE , count ,
                 count , COMM , MPI_STATUS_IGNORE ) ;
                print_matrix ( A , n_loc, N ) ;
            }
        }
        else {
            MPI_Send ( A , n_loc * N , MPI_DOUBLE , 0 ,
             irank , COMM );
        } 
    #endif

    /*Final output and deallocation of the memory*/
    free(A); free(B); free(C);
    MPI_Barrier(COMM);
    if (irank == MASTER) {
        printf("\n----@");
        printf_green();
        printf(" Execution ended with a succes! ");
        printf_reset();
        printf("@----\n");
    }
    MPI_Finalize();

    return 0;
}