#include <stdlib.h>
#include <stdio.h>

#define THREADS_PER_BLOCK 256
#define N_rows 4
#define N (N_rows * N_rows)

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

void transpose( int *A, int dim) {
    int temp;
    int * ptr1;
    int * ptr2;
    int offset = 0;
    for (int j = 0; j < dim; j++)
    {
        for (int i = 1; i < dim; i++)
        {

            temp = A[i];
            *ptr1 = *ptr2;
            *ptr2 = temp;    
        }
        offset++  
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
    cudaMemcpy( dev_a, a, size, cudaMemcpyHostToDevice );

    print_matrix_square(A, N_rows);

    return 0;
}