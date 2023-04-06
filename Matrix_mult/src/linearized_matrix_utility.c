#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linearized_matrix_utility.h"


void create_null_matrix (double * A, int dim) {
  memset( A, 0, dim * dim * sizeof(double) );
}

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

