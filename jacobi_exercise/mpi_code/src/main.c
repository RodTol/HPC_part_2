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
// evolve Jacobi
void evolve( double * matrix, double *matrix_new, size_t dimension );

// return the elapsed time
double seconds( void );

/*** end function declaration ***/

int main(int argc, char* argv[]){

  // timing variables
  double t_start, t_end, increment;

  // indexes for loops
  size_t it;
  
  // initialize matrix
  double *matrix, *matrix_new, *tmp_matrix;

  size_t dimension = 0, iterations = 0, row_peek = 0, col_peek = 0;
  size_t matrix_local_dimension = 0;

  int n_proc_tot, irank;
  int n_loc, rest, dim_2_local;
  int *n_rows_local, *dim_1_local, *displacement;

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

  displacement = (int *) malloc(n_proc_tot*sizeof(int));
  calculate_displ(displacement, n_rows_local, n_proc_tot);

#ifdef DEBUG
    if (irank == MASTER) {
        printf("\n # displacement for each processor:\n");
        for (int i = 0; i < n_proc_tot; i++) {
            printf("(rank: %d rows: %d)\n", i, displacement[i]);
        }
        printf("\n");    
    }
    MPI_Barrier(COMM);
#endif

  /*Allocation of the spaces for each processor. Remember that each
  matrix need extra 2 rows for the ghots layer, but for the first and the last
  process these are extra since they have a boundary*/
  dim_1_local = (int *) malloc(n_proc_tot*sizeof(int));
  for (int count = 0; count < n_proc_tot; count++) {
    dim_1_local[count] = n_rows_local[count] + 2;
  }
  dim_2_local = dimension;
  
#ifdef DEBUG
    if (irank == MASTER) {
        printf("\n # of rows for each processor (including ghosts cells):\n");
        for (int i = 0; i < n_proc_tot; i++) {
            printf("(rank: %d rows: %d)\n", i, dim_1_local[i]);
        }
        printf("\n");        
    }
    MPI_Barrier(COMM);
#endif

  /*Each process knows only its size*/
  matrix_local_dimension = sizeof(double) * ( dim_2_local ) * ( dim_1_local[irank] );
  matrix = ( double* )malloc( matrix_local_dimension );
  matrix_new = ( double* )malloc( matrix_local_dimension );

  /*Both are set to zero*/
  memset( matrix, 0, matrix_local_dimension );
  memset( matrix_new, 0, matrix_local_dimension );

t_start = seconds();
  create_jacobi_start_distributed(matrix, irank, dim_1_local, dim_2_local,
   displacement, n_proc_tot);
  create_jacobi_start_distributed(matrix_new, irank, dim_1_local, dim_2_local,
   displacement, n_proc_tot);
  MPI_Barrier(COMM);
t_end = seconds();

#ifdef DEBUG
  print_matrix_distributed(matrix, irank, dim_1_local, dim_2_local,
    n_proc_tot, COMM, true);
#endif
  print_matrix_distributed_gnuplot(matrix, irank, dim_1_local, dim_2_local,
    displacement, n_proc_tot, COMM, "initial.dat");

  if (irank==MASTER) {
    printf_red();
    printf( "\nElapsed initialisation time = %f seconds\n", t_end - t_start );
    printf_reset();
  }

  // start algorithm
  t_start = seconds();
  for( it = 0; it < iterations; ++it ){

    evolve_mpi( matrix, matrix_new, dimension);

    // swap the pointers
    tmp_matrix = matrix;
    matrix = matrix_new;
    matrix_new = tmp_matrix;
  }
  t_end = seconds();

  if (irank==MASTER) {
    printf_red();
    printf( "\nElapsed computation time = %f seconds\n", t_end - t_start );
    printf_reset();
  }

  print_matrix_distributed_gnuplot(matrix, irank, dim_1_local, dim_2_local,
    displacement, n_proc_tot, COMM, "solution.dat");

  free( matrix );
  free( matrix_new );

  MPI_Finalize();

  return 0;
}


void evolve_mpi( double * matrix, double *matrix_new, size_t dimension ){
  
  size_t i , j;

  //This will be a row dominant program.
  for( i = 1 ; i <= dimension; ++i ) {
    for( j = 1; j <= dimension; ++j ) {
      matrix_new[ ( i * ( dimension + 2 ) ) + j ] = ( 0.25 ) * 
        ( matrix[ ( ( i - 1 ) * ( dimension + 2 ) ) + j ] + 
          matrix[ ( i * ( dimension + 2 ) ) + ( j + 1 ) ] + 	  
          matrix[ ( ( i + 1 ) * ( dimension + 2 ) ) + j ] + 
          matrix[ ( i * ( dimension + 2 ) ) + ( j - 1 ) ] ); 
    }
  }

}

// A Simple timer for measuring the walltime
double seconds(){

    struct timeval tmp;
    double sec;
    gettimeofday( &tmp, (struct timezone *)0 );
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;
}

