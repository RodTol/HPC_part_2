#include <stdio.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>

void initialise_cuda(double *A, double **dev_A, double **dev_B_col, double **dev_C,
 int *n_rows_local, int N, int n_loc, int irank, cublasHandle_t *handle);

void computation(int count, double *B_col, double *dev_A, double *dev_B_col, double *dev_C,
 int *n_row_local, int* displacement, int N, int n_loc, int irank, float *computation_Time, cublasHandle_t handle);
