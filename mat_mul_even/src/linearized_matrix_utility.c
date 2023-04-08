#include "headers/linearized_matrix_utility.h"

void create_identity_matrix (double * A, int dim) {
  memset( A, 0, dim * dim * sizeof(double) );
  for (int i = 0; i < dim; i++){
    A[i + i*dim] = 1.0;
  }
  
}

void print_matrix_square(double * A, int dim ){
  fprintf( stdout, "\n");
  for(int i = 0; i < dim; i++ ){
    for(int j = 0; j < dim; j++ ){
      fprintf( stdout, "%.f ", A[ j + ( i * dim ) ] );
    }
    fprintf( stdout, "\n");
  }
}

/**
 * @brief This function print a matrix of dimension 
 * dim1 X dim2
 * 
 * @param A matrix
 * @param dim_1 rows dimension
 * @param dim_2 cols dimension
 */
void print_matrix(double * A, int dim_1, int dim_2 ) {
  fprintf( stdout, "\n");
  for(int i = 0; i < dim_1; i++ ){
    for(int j = 0; j < dim_2; j++ ){
      fprintf( stdout, "%.f ", A[ j + ( i * dim_2 ) ] );
    }
    fprintf( stdout, "\n");
  }
}

