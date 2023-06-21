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
//Send and receive ghosts layers
void ghost_layer_transfer(double * matrix, int irank, int n_proc_tot, int * dim_1_local, int dim_2_local);

// evolve Jacobi
void evolve_mpi( double * matrix, double *matrix_new, int * dim_1_local, int dim_2_local, int irank );

// return the elapsed time
double seconds( void );

/*** end function declaration ***/

int main(int argc, char* argv[]){

  // timing variables
  double t_start, t_end, time, max_time;

  // indexes for loops
  size_t it;
  
  // initialize matrix
  double *matrix, *matrix_new, *tmp_matrix; 

  size_t dimension = 0, iterations = 0, row_peek = 0, col_peek = 0;
  size_t matrix_local_dimension = 0, matrix_with_borders_dim=0;

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

  // print some info
  if (irank == MASTER) {
    printf_yellow();
    printf("--------------------------------------------------------------------------\n"
           "Matrix size: %zu, # of iteration: %zu, Elements for checking =  Mat[%zu,%zu]\n"
           "--------------------------------------------------------------------------\n",
          dimension, iterations, row_peek, col_peek);
    printf_reset();
  }
   
  matrix_with_borders_dim =  dimension+2;

  if((row_peek > matrix_with_borders_dim) || (col_peek > matrix_with_borders_dim+2)){
    if (irank == MASTER) {
      printf_red();
      fprintf(stderr, "Cannot Peek a matrix element outside of the matrix dimension\n");
      fprintf(stderr, "Arguments n and m must be smaller than %zu\n", matrix_with_borders_dim);
      printf_reset();
      return 1;
    }
  }

  /*Calculate matrix distribution sizes*/
  n_loc = matrix_with_borders_dim/n_proc_tot;
  rest = matrix_with_borders_dim % n_proc_tot;
  if (irank == MASTER) {
    printf_yellow();
    printf("------------------------------------------------------------\n"
           "Total matrix size: %zu, n_proc_tot: %d, n_loc: %d, rest: %d \n"
           "------------------------------------------------------------\n",
            matrix_with_borders_dim, n_proc_tot, n_loc, rest);
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
  matrix need extra 2 rows for the ghots layer. The last and first
  process already have an extra row since they have the boundary*/
  dim_1_local = (int *) malloc(n_proc_tot*sizeof(int));
  
  dim_1_local[0] = n_rows_local[0]+1;
  dim_1_local[n_proc_tot-1] = n_rows_local[n_proc_tot-1]+1;
  for (int count = 1; count < n_proc_tot-1; count++) {
    dim_1_local[count] = n_rows_local[count] + 2;
  }
  dim_2_local = matrix_with_borders_dim;
  
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
  tmp_matrix = ( double* )malloc( matrix_local_dimension );


  /*Both are set to zero*/
  memset( matrix, 0, matrix_local_dimension );
  memset( matrix_new, 0, matrix_local_dimension );
  memset( tmp_matrix, 0, matrix_local_dimension );


  t_start = seconds();
  create_jacobi_start_distributed(matrix, irank, dim_1_local, dim_2_local,
   displacement, n_proc_tot);
  create_jacobi_start_distributed(matrix_new, irank, dim_1_local, dim_2_local,
   displacement, n_proc_tot);
  MPI_Barrier(COMM);
  t_end = seconds();
  time = t_end-t_start;
  
#ifdef DEBUG
  if (dimension <=11) {  print_matrix_distributed(matrix, irank, dim_1_local, dim_2_local,
    n_proc_tot, COMM, true); }
#endif

  MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (irank==MASTER) {
    printf_red();
    printf( "\nElapsed initialisation time = %f seconds\n", max_time );
    printf_reset();
  }


  if (irank==MASTER) {
    printf_yellow();
    printf("\n");
    printf("-------Printing the initial matrix-------\n");
    printf_reset();
  }

  print_matrix_distributed_file(matrix, irank, dim_1_local, dim_2_local,
    displacement, n_proc_tot, COMM, "initial.dat");

  MPI_Barrier(MPI_COMM_WORLD);

  // start algorithm
  t_start = seconds();
  for( it = 0; it < iterations; ++it ){
    ghost_layer_transfer(matrix, irank, n_proc_tot, dim_1_local, dim_2_local);
    evolve_mpi(matrix, matrix_new, dim_1_local, dim_2_local, irank);

    // swap the pointers
    tmp_matrix = matrix;
    matrix = matrix_new;
    matrix_new = tmp_matrix;
  }
  t_end = seconds();
  time = t_end-t_start;
  MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (irank==MASTER) {
    printf_red();
    printf( "\nElapsed computation time = %f seconds\n", max_time );
    printf_reset();
  }

#ifdef DEBUG
  MPI_Barrier(COMM);
  if (irank==MASTER) {
    printf_yellow();
    printf("\n");
    printf("-------Evolution executed-------\n");
    printf_reset();
  }
  if (dimension <=11) {print_matrix_distributed(matrix, irank, dim_1_local, dim_2_local,
    n_proc_tot, COMM, true);}
#endif

  MPI_Barrier(MPI_COMM_WORLD);

  if (irank==MASTER) {
    printf_yellow();
    printf("\n");
    printf("-------Printing the final matrix-------\n");
    printf_reset();
  }

  print_matrix_distributed_file(matrix, irank, dim_1_local, dim_2_local,
    displacement, n_proc_tot, COMM, "solution.dat");


 /*I save the result in a times.dat file*/
    FILE* file;
    char* title = "times.dat";
    if (irank == MASTER) {
        if (!file_exists(title)) {
            file = fopen(title, "w");
            fprintf(file, "size, n_proc_tot, it, time\n");
            fclose(file);
        }

        file = fopen(title, "a");
        fprintf(file, "%ld  %d   %ld %15.12f\n", dimension, n_proc_tot, it, max_time);
        fclose(file);
    }

  free(matrix);
  free(matrix_new);

  /*I think it's optional because one pointer is set
  equal to another*/
  //free(tmp_matrix);

  free(n_rows_local);
  free(dim_1_local);
  free(displacement);

  MPI_Finalize();

  return 0;
}

void ghost_layer_transfer(double * matrix, int irank, int n_proc_tot, int * dim_1_local, int dim_2_local) {
  
  int next = irank+1;
  int previous = irank-1;

  if (irank==0) {
    /*Send the last actual row to next and receive
    from next in its bottom ghost layer*/
    MPI_Sendrecv(matrix+(dim_1_local[irank]-2)*dim_2_local, dim_2_local,
     MPI_DOUBLE, next, 0,
     matrix+(dim_1_local[irank]-1)*dim_2_local, dim_2_local,
     MPI_DOUBLE, next, 0, COMM, MPI_STATUS_IGNORE);
  } else if (irank==n_proc_tot-1) {
    /*Send the first actual row to previous and receive
    from previous in its top ghost layer*/
    MPI_Sendrecv(matrix+dim_2_local, dim_2_local,
     MPI_DOUBLE, previous, 0,
     matrix, dim_2_local,
     MPI_DOUBLE, previous, 0, COMM, MPI_STATUS_IGNORE);
  } else {
    /*Send the first actual row to previous and receive
    from previous in its top ghost layer*/
    MPI_Sendrecv(matrix+dim_2_local, dim_2_local,
     MPI_DOUBLE, previous, 0,
     matrix, dim_2_local,
     MPI_DOUBLE, previous, 0, COMM, MPI_STATUS_IGNORE);
    /*Send the last actual row to next and receive
    from next in its bottom ghost layer*/
    MPI_Sendrecv(matrix+(dim_1_local[irank]-2)*dim_2_local, dim_2_local,
     MPI_DOUBLE, next, 0,
     matrix+(dim_1_local[irank]-1)*dim_2_local, dim_2_local,
     MPI_DOUBLE, next, 0, COMM, MPI_STATUS_IGNORE);
  }

}

void evolve_mpi( double * matrix, double *matrix_new, int * dim_1_local, int dim_2_local, int irank ) {
  size_t i , j;

  //This will be a row dominant program.
  for( i = 1 ; i <= dim_1_local[irank]-2; ++i ) {
    for( j = 1; j <= dim_2_local-2; ++j ) {
      matrix_new[ linear_index(i,j,dim_1_local[irank],dim_2_local) ] = ( 0.25 ) * 
        ( matrix[ linear_index(i-1,j,dim_1_local[irank],dim_2_local) ] + 
          matrix[ linear_index(i,j+1,dim_1_local[irank],dim_2_local) ] + 	  
          matrix[ linear_index(i+1,j,dim_1_local[irank],dim_2_local) ] + 
          matrix[ linear_index(i,j-1,dim_1_local[irank],dim_2_local) ] ); 
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

