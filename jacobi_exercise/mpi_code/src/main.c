#include "headers/general_utility.h"
#include "headers/linearized_matrix_utility.h"
#include "headers/data_distr_utility.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

#define COMM MPI_COMM_WORLD
#define MASTER 0

/*** function declarations ***/

// save matrix to file
void save_gnuplot( double *M, size_t dim );

// evolve Jacobi
void evolve( double * matrix, double *matrix_new, size_t dimension );

// return the elapsed time
double seconds( void );

/*** end function declaration ***/

int main(int argc, char* argv[]){

  // timing variables
  double t_start, t_end, increment;

  // indexes for loops
  size_t i, j, it;
  
  // initialize matrix
  double *matrix, *matrix_new, *tmp_matrix;

  size_t dimension = 0, iterations = 0, row_peek = 0, col_peek = 0;
  size_t byte_dimension = 0;

  int n_proc_tot, irank;
  int n_loc, rest;
  int* n_rows_local;

  MPI_Init ( & argc , & argv ) ;
  MPI_Comm_rank ( COMM , & irank ) ;
  MPI_Comm_size ( COMM , & n_proc_tot ) ;

  // check on input parameters
  if(irank == MASTER && argc != 5) {
    printf_red();
    fprintf(stderr,"\nwrong number of arguments. Usage: ./a.out dim it n m\n");
    printf_reset();
    return 1;
  }

  dimension = atoi(argv[1]);
  iterations = atoi(argv[2]);
  row_peek = atoi(argv[3]);
  col_peek = atoi(argv[4]);

  if (irank == MASTER) {
    printf_yellow();
    printf("--------------------------------------------------------------------------\n"
           "Matrix size: %zu, # of iteration: %zu, Elements for checking =  Mat[%zu,%zu]\n"
           "--------------------------------------------------------------------------\n",
          dimension, iterations, row_peek, col_peek);
    printf_reset();
  }

  if((row_peek > dimension) || (col_peek > dimension)){
    if (irank == MASTER) {
      printf_red();
      fprintf(stderr, "Cannot Peek a matrix element outside of the matrix dimension\n");
      fprintf(stderr, "Arguments n and m must be smaller than %zu\n", dimension);
      printf_reset();
      return 1;
    }
  }

  /*Calculate matrix distribution sizes*/
  n_loc = dimension/n_proc_tot;
  rest = dimension % n_proc_tot;
  if (irank == MASTER) {
    printf_yellow();
    printf("-------------------------------------------\n"
           "N: %zu, n_proc_tot: %d, n_loc: %d, rest: %d \n"
           "-------------------------------------------\n", dimension, n_proc_tot, n_loc, rest);
    printf_reset();
  }

  /*Now I can have a rest so, I need to calculate the
  correct sizes before the allocation. Every process needs
  this array*/
  n_rows_local = (int *) malloc( n_proc_tot * sizeof(int) );
  calculate_n_rows(n_rows_local, n_loc, rest, n_proc_tot);

#ifdef DEBUG
    if (irank == MASTER) {
        printf("\n # of rows for each processor:\n");
        for (int i = 0; i < n_proc_tot; i++) {
            printf("(rank: %d rows: %d)\n", i, n_rows_local[i]);
        }    
    }
    MPI_Barrier(COMM);
#endif

  /*Allocation of the spaces for each processor. Remember that each
  matrix need extra 2 rows and cols for the ghots layer*/
  byte_dimension = sizeof(double) * ( dimension + 2 ) * ( n_loc + 2 );
  matrix = ( double* )malloc( byte_dimension );
  matrix_new = ( double* )malloc( byte_dimension );

  /*Both are set to zero*/
  memset( matrix, 0, byte_dimension );
  memset( matrix_new, 0, byte_dimension );

  //fill initial values  
  
  /*for( i = 1; i <= dimension; ++i )
    for( j = 1; j <= dimension; ++j )
      matrix[ ( i * ( dimension + 2 ) ) + j ] = 0.5;
	*/
      
  // set up borders 
  increment = 100.0 / ( dimension + 1 );
  
  for( i=1; i <= dimension+1; ++i ){
    matrix[ i * ( dimension + 2 ) ] = i * increment;
    matrix[ ( ( dimension + 1 ) * ( dimension + 2 ) ) + ( dimension + 1 - i ) ] = i * increment;
    matrix_new[ i * ( dimension + 2 ) ] = i * increment;
    matrix_new[ ( ( dimension + 1 ) * ( dimension + 2 ) ) + ( dimension + 1 - i ) ] = i * increment;
  }
  
  // start algorithm
  t_start = seconds();
  for( it = 0; it < iterations; ++it ){
    
    evolve( matrix, matrix_new, dimension );

    // swap the pointers
    tmp_matrix = matrix;
    matrix = matrix_new;
    matrix_new = tmp_matrix;

  }
  t_end = seconds();
  
  printf( "\nelapsed time = %f seconds\n", t_end - t_start );
  printf( "\nmatrix[%zu,%zu] = %f\n", row_peek, col_peek, matrix[ ( row_peek + 1 ) * ( dimension + 2 ) + ( col_peek + 1 ) ] );

  save_gnuplot( matrix, dimension );
  
  free( matrix );
  free( matrix_new );

  return 0;
}

void evolve( double * matrix, double *matrix_new, size_t dimension ){
  
  size_t i , j;

  //This will be a row dominant program.
  for( i = 1 ; i <= dimension; ++i )
    for( j = 1; j <= dimension; ++j )
      matrix_new[ ( i * ( dimension + 2 ) ) + j ] = ( 0.25 ) * 
	( matrix[ ( ( i - 1 ) * ( dimension + 2 ) ) + j ] + 
	  matrix[ ( i * ( dimension + 2 ) ) + ( j + 1 ) ] + 	  
	  matrix[ ( ( i + 1 ) * ( dimension + 2 ) ) + j ] + 
	  matrix[ ( i * ( dimension + 2 ) ) + ( j - 1 ) ] ); 
}

void save_gnuplot( double *M, size_t dimension ){
  
  size_t i , j;
  const double h = 0.1;
  FILE *file;

  file = fopen( "solution.dat", "w" );

  for( i = 0; i < dimension + 2; ++i )
    for( j = 0; j < dimension + 2; ++j )
      fprintf(file, "%f\t%f\t%f\n", h * j, -h * i, M[ ( i * ( dimension + 2 ) ) + j ] );

  fclose( file );

}

// A Simple timer for measuring the walltime
double seconds(){

    struct timeval tmp;
    double sec;
    gettimeofday( &tmp, (struct timezone *)0 );
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;
}

