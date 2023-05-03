#include "headers/gpu_computation.h"
#include <stdlib.h>
#include <unistd.h>

void initialise_cuda(double *A, double **dev_A, double **dev_B_col, double **dev_C,
 int *n_rows_local, int N, int n_loc, int irank, cublasHandle_t *handle) {
    int n_gpus;
    cudaGetDeviceCount(&n_gpus);
    cudaSetDevice(irank % n_gpus);
    cublasCreate(handle);
    // Allocate memory on the device
    cudaMalloc( (void **) dev_A, n_rows_local[irank] * N * sizeof(double) );
    cudaMalloc( (void **) dev_B_col, N * (n_loc + 1) * sizeof(double));
    cudaMalloc( (void **) dev_C, n_rows_local[irank] * N * sizeof(double) );

    cudaMemcpy(*dev_A, A, n_rows_local[irank] * N * sizeof(double), cudaMemcpyHostToDevice);
}

/*
void print_matrix_2(double * A, int dim_1, int dim_2 ) {
  //fprintf( stdout, "\n");
  for(int i = 0; i < dim_1; i++ ){
    for(int j = 0; j < dim_2; j++ ){
      fprintf( stdout, "%.3g ", A[ j + ( i * dim_2 ) ] );
    }
    fprintf( stdout, "\n");
  }
}
*/

void computation(int count, double *B_col, double *dev_A, double *dev_B_col, double *dev_C,
 int *n_rows_local, int* displacement, int N, int n_loc, int irank, float *computation_Time, cublasHandle_t handle) {

/*
#ifdef DEBUG
    double *A_tmp;
    int size= N * n_rows_local[irank] * sizeof( double);
    A_tmp = (double *) malloc( size );
    sleep(2*irank);

    printf("\n I am: %d and this is my dev_A at count: %d before the copy\n", irank, count);
    cudaMemcpy(A_tmp, dev_A, size, cudaMemcpyDeviceToHost);
    print_matrix_2(A_tmp,n_rows_local[irank],N);
    printf("\n");
    free(A_tmp);

    printf("\n I am: %d and this is my B_col at count: %d before the copy\n", irank, count);
    print_matrix_2(B_col, N, n_rows_local[count]);
    printf("\n");
#endif
*/
    cudaMemcpy(dev_B_col, B_col, N * (n_loc + 1) * sizeof(double), cudaMemcpyHostToDevice);

    const double alpha = 1.0, beta = 0.0;
    cudaEvent_t start_c, stop_c;
    float time;

    cudaEventCreate(&start_c);
    cudaEventCreate(&stop_c);
    cudaEventRecord(start_c, 0);

/*
#ifdef DEBUG
    double *B_col_tmp;
    int size_B_col= N*n_rows_local[count] * sizeof(double);
    B_col_tmp = (double *) malloc(size_B_col);

    printf("\n I am: %d and this is my dev_B_col at count: %d after the copy\n", irank, count);
    cudaMemcpy(B_col_tmp, dev_B_col, size_B_col, cudaMemcpyDeviceToHost);
    print_matrix_2(B_col_tmp, N, n_rows_local[count]);
    printf("\n");
#endif
*/
    cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, n_rows_local[count], n_rows_local[irank], N, &alpha,
     dev_B_col, n_rows_local[count], dev_A, N, &beta, dev_C + displacement[count], N);

/*
#ifdef DEBUG
    double *C_tmp;
    int size_C = N*n_rows_local[irank] * sizeof(double);
    C_tmp = (double *) malloc(size_C);

    printf("\n I am: %d and this is my dev_C at count: %d after the copy\n", irank, count);
    cudaMemcpy(C_tmp, dev_C, size_C, cudaMemcpyDeviceToHost);
    print_matrix_2(C_tmp, n_rows_local[irank], N);
    printf("\n");
#endif
*/
   cudaEventRecord(stop_c, 0);
   cudaEventSynchronize(stop_c);
   cudaEventElapsedTime(&time, start_c, stop_c);
   *computation_Time += time;
}
