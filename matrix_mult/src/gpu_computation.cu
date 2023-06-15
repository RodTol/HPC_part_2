#include "headers/gpu_computation.h"

void initialise_cuda(double *A, double **dev_A, double **dev_B_col, double **dev_C,
 int *n_rows_local, int N, int n_loc, int irank, cublasHandle_t *handle) {
    int n_gpus;
    cudaGetDeviceCount(&n_gpus);
    cudaSetDevice(irank % n_gpus);
    cublasCreate(handle);

#ifdef DEBUG
    if (irank==0)
    {
      printf("I found %i devices\n", n_gpus);
    }
    
#endif

    
    // Allocate memory on the device
    cudaError_t errA = cudaMalloc( (void **) dev_A, n_rows_local[irank] * N * sizeof(double) );
    if (errA != cudaSuccess) printf("Error allocating memory of A on the device: %s\n", cudaGetErrorString(errA));
    cudaError_t errB = cudaMalloc( (void **) dev_B_col, N * (n_loc + 1) * sizeof(double));
    if (errB != cudaSuccess) printf("Error allocating memory of B_col on the device: %s\n", cudaGetErrorString(errB));
    cudaError_t errC = cudaMalloc( (void **) dev_C, n_rows_local[irank] * N * sizeof(double) );
    if (errC != cudaSuccess) printf("Error allocating memory of C on the device: %s\n", cudaGetErrorString(errC));


    cudaError_t err = cudaMemcpy(*dev_A, A, n_rows_local[irank] * N * sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) printf("Error copying A on the device: %s\n", cudaGetErrorString(err));

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
    float time;
    cudaEvent_t start,stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start,0);
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
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);  
    *computation_Time += time;

#ifdef DEBUG
    if (irank==0) {
      printf("Computational time for cuda at count %d : %15.17f", count, *computation_Time);
    }
    
#endif
}
