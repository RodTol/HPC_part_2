#include <stdio.h>
#include <mpi.h>

#include "headers/linearized_matrix_utility.h"
#include "headers/general_utility.h"
#include "headers/prova.h"

#define MASTER 0
#define COMM MPI_COMM_WORLD

int main(int argc, char** argv) {
    
    double *A, *B, *C, *B_col;

    int N = atoi(argv[1]);
    int n_proc_tot = 1, irank = 0;
    int n_loc = N, rest=0;
    double start_compute, end_compute, start_comm,
     end_comm, compute_total, comm_total;
    //int i_local = 0, j_global = 0, offset = 0;

    /*MPI setup*/
    MPI_Init ( & argc , & argv ) ;
    MPI_Comm_rank ( COMM , & irank ) ;
    MPI_Comm_size ( COMM , & n_proc_tot ) ;

    MPI_Datatype Even_block;

    /*Parameters for the matrix distribution*/
    n_loc = N/n_proc_tot;
    rest = N % n_proc_tot;
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

    
    /*Allocation and initialisation*/
    int size= N * n_loc * sizeof( double);
    A = (double *) malloc( size );
    array_of_random_doubles(A, n_loc*N);

    B = (double *) malloc( size );
    //create_identity_matrix_distributed(B, irank, n_loc, N, offset, n_proc_tot);
    array_of_random_doubles(B, n_loc*N);

    C = (double *) malloc( size );
    create_null_array(C, n_loc*N);

    /*The datatype that will represent the columns of B to be
    * multiplied.
    (#of tot elements, #taken from each proc, #distance from, Type, name)*/
    MPI_Type_vector(n_loc, n_loc, N, MPI_DOUBLE, &Even_block);
    MPI_Type_commit(&Even_block);

    /*Local buffer for the multiplication*/
    B_col = (double *) malloc( n_loc * N * sizeof(double) );
    create_null_array(B_col, n_loc*N);

    /*The multiplication*/
    for (int count = 0; count < n_proc_tot; count++) {

        start_comm = MPI_Wtime();
        
        MPI_Allgather(&B[n_loc * count], 1, Even_block, B_col, n_loc*n_loc, MPI_DOUBLE, COMM);
        
        end_comm = MPI_Wtime();
        start_compute = MPI_Wtime();

        matrix_multiplication(A, B_col, C, N, n_loc, count);

        end_compute = MPI_Wtime();

        compute_total += end_compute-start_compute;
        comm_total += end_comm-start_comm;

#ifdef DEBUG
            MPI_Barrier(COMM);
            if (irank == MASTER) {
                printf("\nMatrix C at count = %d \n", count);
            }
            print_matrix_distributed(C, irank, n_loc, N, n_proc_tot, COMM);
            MPI_Barrier(COMM);
#endif
    }

    /*If the matrices are DEBUG enough (N < 6 and n_proc < 4)
    the program will print the complete form*/
#ifdef DEBUG
        MPI_Barrier(COMM);
        if (irank == MASTER) {
            printf("\nMatrix A \n");
        }
        print_matrix_distributed(A, irank, n_loc, N, n_proc_tot, COMM);
        
        if (irank == MASTER) {
            printf("\nMatrix B \n");
        }
        print_matrix_distributed(B, irank, n_loc, N, n_proc_tot, COMM);

        if (irank == MASTER) {
            printf("\nMatrix C \n");
        }
        print_matrix_distributed(C, irank, n_loc, N, n_proc_tot, COMM);
#endif

    /*Final output and deallocation of the memory*/
    free(A); free(B); free(C), free(B_col);
    MPI_Barrier(COMM);

    FILE* file;
    char* title = "times.dat";
    if (irank == MASTER) {
        if (!file_exists(title)) {
            file = fopen(title, "w");
            fprintf(file, "N, n_proc_tot, comm_total, compute_total\n");
            fclose(file);
        }

        file = fopen(title, "a");
        fprintf(file, "%d %d %15.12f %15.12f\n", N, n_proc_tot, comm_total, compute_total);
        fclose(file);
    }

    if (irank == MASTER) {
        printf("\n----@");
        printf_green();
        printf("Execution ended with a succes!");
        printf_reset();
        printf("@----\n");
        printf_yellow();
        printf("Communication time: %15.12f \n", comm_total);
        printf("Computational time: %15.12f \n", compute_total);
        printf_red();
        printf("Total time:         %15.12f \n", compute_total+comm_total);
        printf_reset();
    }

    MPI_Finalize();



    return 0;
}
