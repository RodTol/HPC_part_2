#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include "linearized_matrix_utility.h"

void create_identity_matrix (int * A, int dim) {
  for (int i = 0; i < dim; i++){
    A[i + i*dim] = 1;
  }
  
}

void print_matrix_square(int * A, int dim ){

  int i , j;
  
  fprintf( stdout, "\n");
  for( i = 0; i < dim; i++ ){
    for( j = 0; j < dim; j++ ){
      fprintf( stdout, "%d ", A[ j + ( i * dim ) ] );
    }
    fprintf( stdout, "\n");
  }
}

