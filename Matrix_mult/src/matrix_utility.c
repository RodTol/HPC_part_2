#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>

#include "matrix_utility.h"


void print_matrix_square( int * A, int dim ){

  int i , j;
  
  fprintf( stdout, "\n");
  for( i = 0; i < dim; i++ ){
    for( j = 0; j < dim; j++ ){
      fprintf( stdout, "%d ", A[ j + ( i * dim ) ] );
    }
    fprintf( stdout, "\n");
  }
}

