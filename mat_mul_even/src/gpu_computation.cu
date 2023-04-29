#include "headers/gpu_computation.h"

void alloc(const int* row_counts, const int n_loc, const int rank, const int N, const double *A,
                        double **d_A, double **d_B_col, double **d_C, cublasHandle_t *handle) {

    int n_of_gpus;
    cudaGetDeviceCount(&n_of_gpus);
    cudaSetDevice(rank % n_of_gpus);
    cublasCreate(handle);
    // Allocate memory on the device
    cudaMalloc( (void **)d_A, row_counts[rank] * N * sizeof(double) );
    cudaError_t err = cudaMalloc((void **)d_B_col, N * (n_loc + 1) * sizeof(double));
    if (err != cudaSuccess) printf("Error allocating memory on the device: %s\n", cudaGetErrorString(err));
    cudaMalloc( (void **)d_C, row_counts[rank] * N * sizeof(double) );
    // copy A to the device (C is already allocated on the gpu)
    cudaError_t err1 = cudaMemcpy(*d_A, A, row_counts[rank] * N * sizeof(double), cudaMemcpyHostToDevice);
    if (err1 != cudaSuccess) { printf("Error on copying A to d_A: %s\n", cudaGetErrorString(err1)); }
}