#include <stdlib.h>
#include <stdio.h>

#define THREADS_PER_BLOCK 25
#define N_cols 5
#define N (N_cols * N_cols)

void random_ints(int * a, int dim)
{
   int i;
   for (i = 0; i < dim; ++i) {
    a[i] = rand()%10;
   }
}

void print_matrix_square( int * A, int dim ){

  int i , j;
  
  for( i = 0; i < dim; i++ ){
    for( j = 0; j < dim; j++ ){
      fprintf( stdout, "%d ", A[ j + ( i * dim ) ] );
    }
    fprintf( stdout, "\n");
  }
}

__global__ void transpose( int *dev_A, int dim) {
    int temp;
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    
    int j = idx%dim;  
    int i = idx/dim;

    int idx_trasp = j*dim +i;

    /*Check indexing is correct*/
    printf("I am Thread id: %d Block id: %d\n", threadIdx.x, blockIdx.x );
    printf("i: %d j: %d value: %d\n", i, j , dev_A[idx] );

    if (j > i)
    {
      temp = dev_A[idx];
      dev_A[idx] = dev_A[idx_trasp];
      dev_A[idx_trasp] = temp;
    }

}


int main(int argc, char * argv[]) {
    int size=N * sizeof( int);
    
    int * A;
    int * dev_A;

    /*Inizializzazione*/
    A = (int *) malloc( size );
    //memset( A, 0, N * sizeof(int) );
    random_ints(A, N);

    cudaMalloc( (void**)&dev_A, size );
    cudaMemcpy( dev_A, A, size, cudaMemcpyHostToDevice );

    print_matrix_square(A, N_cols);

    transpose <<<N/THREADS_PER_BLOCK, THREADS_PER_BLOCK>>> (dev_A, N_cols); 

    cudaMemcpy( A, dev_A, size, cudaMemcpyDeviceToHost );
    
    printf("\n");
    print_matrix_square(A, N_cols);

    cudaFree(dev_A);
    free(A);

    return 0;
}