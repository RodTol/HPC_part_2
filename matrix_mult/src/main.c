#include <stdio.h>
#include <mpi.h>

#ifdef DGEMM
	#include <cblas.h>
#endif

#ifdef GPU
	#include"gpu_computation.cu" 
#endif

#include "headers/linearized_matrix_utility.h"
#include "headers/general_utility.h"
#include "headers/data_distr_utility.h"

#define MASTER 0
#define COMM MPI_COMM_WORLD

int main(int argc, char** argv) {
    
    double *A, *B, *C, *B_col;

    int N = atoi(argv[1]);
    int n_proc_tot = 1, irank = 0;
    int n_loc, rest;
    double start_compute, end_compute, start_comm,
     end_comm, compute_total=0.0, comm_total=0.0;
    int *n_rows_local, *displacement, *n_elements_local, *displacement_col;
    
    /*MPI setup*/
    MPI_Init ( & argc , & argv ) ;
    MPI_Comm_rank ( COMM , & irank ) ;
    MPI_Comm_size ( COMM , & n_proc_tot ) ;

    /*Parameters for the matrix distribution*/
    n_loc = N/n_proc_tot;
    rest = N % n_proc_tot;
    if (irank == MASTER) { printf("-------------------------------------------\n"
                             "N: %d, n_proc_tot: %d, n_loc: %d, rest: %d \n"
                             "-------------------------------------------\n", N, n_proc_tot, n_loc, rest);
#ifdef DGEMM
        printf_red();
        printf("Using OpenBLAS DGEMM for computation\n");
        printf_reset();
#elif GPU
        printf_red();
        printf("Using cuBLAS DGEMM for computation\n");
        printf_reset();
#else 
        printf_red();
        printf("Using basic MPI for computation\n");
        printf_reset();
#endif
    }
    
    /*Now I can have a rest so, I need to calculate the
    correct sizes before the allocation. Every process needs
    this array*/
    n_rows_local = (int *) malloc( n_proc_tot * sizeof(int) );
    calculate_n_rows(n_rows_local, n_loc, rest, n_proc_tot);

#ifdef DEBUG
    if (irank == MASTER) {
        printf("\n # of rows for each processor:\n");
        for (int i = 0; i < n_proc_tot; i++) {
            printf("(rank: %d rows: %d)\n", i, n_rows_local[i]);
        }    
    }
    MPI_Barrier(COMM);
#endif

    /*This the offset. It will be necessary to create the column
    block for the multiplication. All the process needs this array*/
    displacement = (int *) malloc(n_proc_tot*sizeof(int));
    calculate_displ(displacement, n_rows_local, n_proc_tot);

#ifdef DEBUG
    if (irank == MASTER) {
        printf("\n # displacement for each processor:\n");
        for (int i = 0; i < n_proc_tot; i++) {
            printf("(rank: %d rows: %d)\n", i, displacement[i]);
        }    
    }
    MPI_Barrier(COMM);
#endif

    /*Allocation and initialisation*/
    int size= N * n_rows_local[irank] * sizeof( double);
    A = (double *) malloc( size );
    array_of_random_doubles(A, n_rows_local[irank]*N);

    B = (double *) malloc( size );
    create_identity_matrix_distributed(B, irank, n_rows_local[irank], N, displacement[irank]);
    //array_of_random_doubles(B, n_rows_local[irank]*N);

    C = (double *) malloc( size );
    create_null_array(C, n_rows_local[irank]*N);
    
    /*Column buffer for the multiplication*/
    B_col = (double *) malloc( (n_loc+1) * N * sizeof(double) );
    create_null_array(B_col, (n_loc+1)*N);

    /*How many elements each processor should give to the 
    column buffer. This array will be updated at each step, since
    the # is not constant*/
    n_elements_local = (int *) malloc ( n_proc_tot * sizeof(int) );

    /*This is the array of the offset for each processor. It
    will be updated at each step of the multiplication with the
    correct size, since the # of rows is not constant.*/
    displacement_col = (int *) malloc ( n_proc_tot * sizeof(int) );

#ifdef GPU
    double *dev_A, *dev_B_col, *dev_C;

    cudaEvent_t start, stop;
    float TotalTime, max_TotalTime, computation_Time, max_computation_Time;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    initialise_cuda(A, &dev_A, &dev_B_col, &dev_C, n_rows_local, N, n_loc, irank);
 #ifdef DEBUG   
    printf("Initialize cuda completed");
 #endif
#endif

    /*The multiplication*/
    for (int count = 0; count < n_proc_tot; count++) {

        start_comm = MPI_Wtime();
        
        /*I update the columnd displacement*/
        calculate_n_elements(n_elements_local, n_rows_local, count, n_proc_tot);  
        calculate_displ_col(displacement_col, n_elements_local, n_proc_tot);

#ifdef DEBUG
        if (irank == MASTER) {
            printf("\n displacement_col for each procs, and #of elements:\n");
            for (int i = 0; i < n_proc_tot; i++) {
                printf("(rank: %d displacement_col: %d elements: %d)\n", i, displacement_col[i], n_elements_local[i]);
            }    
        }
        MPI_Barrier(COMM);
#endif

        MPI_build_column(n_rows_local, displacement, n_elements_local, displacement_col,
           B, B_col, irank, count, N);

    
        end_comm = MPI_Wtime();
        start_compute = MPI_Wtime();

#ifdef DEBUG
        //printf("\n I am: %d and this is my B_col at count: %d \n", irank, count);
        //print_matrix(B_col, N, n_rows_local[count]);
        //printf("\n");
        //MPI_Barrier(COMM);
#endif

#ifdef DGEMM
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
            n_rows_local[irank], n_rows_local[count], N, // m, n, k
            1.0, A, N, B_col, n_rows_local[count], 0.0, C + displacement[count], N);
#elif GPU
        cublasHandle_t handle;
        computation(count, B_col, dev_A, dev_B_col, dev_C,
         n_rows_local, displacement, N, n_loc, irank, &computation_Time, handle);
#else
        matrix_multiplication(A, B_col, C, N, n_rows_local, displacement, irank, count);
#endif

#ifdef GPU
        cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&TotalTime, start, stop);
        MPI_Reduce(&TotalTime, &max_TotalTime, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&computation_Time, &max_computation_Time, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);     

        comm_total = TotalTime - max_computation_Time;
        compute_total = max_computation_Time;
#endif

        end_compute = MPI_Wtime();
        compute_total += end_compute-start_compute;
        comm_total += end_comm-start_comm;


#ifdef DEBUG
        if (irank == MASTER) {
            printf("\nMatrix C at count = %d \n", count);
        }
        print_matrix_distributed(C, irank, n_rows_local, N, n_proc_tot, COMM);
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
        print_matrix_distributed(A, irank, n_rows_local, N, n_proc_tot, COMM);
        
        if (irank == MASTER) {
            printf("\nMatrix B \n");
        }
        print_matrix_distributed(B, irank, n_rows_local, N, n_proc_tot, COMM);

        if (irank == MASTER) {
            printf("\nMatrix C \n");
        }
        print_matrix_distributed(C, irank, n_rows_local, N, n_proc_tot, COMM);
#endif

    /*Final output and deallocation of the memory*/
    free(A); free(B); free(C), free(B_col);

#ifdef GPU
    cudaMemcpy(C, dev_C, n_rows_local[irank] * N * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(dev_A);
    cudaFree(dev_B_col);
    cudaFree(dev_C);
#endif

    MPI_Barrier(COMM);

    FILE* file;
    char* title = "times.dat";
    if (irank == MASTER) {
        if (!file_exists(title)) {
            file = fopen(title, "w");
            fprintf(file, "         N, n_proc_tot, comm_total, compute_total\n");
            fclose(file);
        }

        file = fopen(title, "a");
#ifdef DGEMM
        fprintf(file, "OpenBLAS ");

#elif GPU
        fprintf(file, "cuBLAS   ");
#else 
        fprintf(file, "MPI      ");
#endif
        fprintf(file, "%d  %d   %15.12f %15.12f\n", N, n_proc_tot, comm_total, compute_total);
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