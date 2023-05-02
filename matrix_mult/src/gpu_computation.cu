#include "headers/gpu_computation.h"

void initialise_cuda(double *A, double **dev_A, double **dev_B_col, double **dev_C,
 int *n_rows_local, int N, int n_loc, int irank) {
    int n_gpus;
    cudaGetDeviceCount(&n_gpus);
    cudaSetDevice(irank % n_gpus);
    // Allocate memory on the device
    cudaMalloc( (void **) dev_A, n_rows_local[irank] * N * sizeof(double) );
    cudaMalloc( (void **) dev_B_col, N * (n_loc + 1) * sizeof(double));
    cudaMalloc( (void **) dev_C, n_rows_local[irank] * N * sizeof(double) );

    cudaMemcpy(*dev_A, A, n_rows_local[irank] * N * sizeof(double), cudaMemcpyHostToDevice);
}

void computation(int count, double *B_col, double *dev_A, double *dev_B_col, double *dev_C,
 int *n_row_local, int* displacement, int N, int n_loc, int irank, float *computation_Time, cublasHandle_t handle) {
    cudaMemcpy(dev_B_col, B_col, N * (n_loc + 1) * sizeof(double), cudaMemcpyHostToDevice);
    
    const double alpha = 1.0, beta = 0.0;
    cudaEvent_t start, stop;
    float time;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    
    cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, n_row_local[count], n_row_local[irank], N, &alpha,
     dev_B_col, n_row_local[count], dev_A, N, &beta, dev_C + displacement[count], N);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    *computation_Time += time;
}