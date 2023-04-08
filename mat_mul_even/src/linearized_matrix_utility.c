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
  //fprintf( stdout, "\n");
  for(int i = 0; i < dim_1; i++ ){
    for(int j = 0; j < dim_2; j++ ){
      fprintf( stdout, "%.f ", A[ j + ( i * dim_2 ) ] );
    }
    fprintf( stdout, "\n");
  }
}

/**
 * @brief This function is used to print the full matrix, 
 * scattered among n_proc_tot processor
 * 
 * @param A the matrix
 * @param irank rank of each processor
 * @param dim_1 rows of each submatrix
 * @param dim_2 cols of each submatrix
 * @param n_proc_tot number of total processors
 * @param COMM MPI communicator
*/
void print_matrix_distributed (double * A, int irank,
 int dim_1 , int dim_2, int n_proc_tot, MPI_Comm COMM) {
  if ( irank == 0 ) {
            print_matrix ( A , dim_1, dim_2) ;
            for (int count = 1; count < n_proc_tot ; count ++ ) {
                MPI_Recv ( A , dim_1 * dim_2 , MPI_DOUBLE , count ,
                 count , COMM , MPI_STATUS_IGNORE ) ;
                print_matrix ( A , dim_1, dim_2 ) ;
            }
        }
        else {
            MPI_Send ( A , dim_1 * dim_2 , MPI_DOUBLE , 0 ,
             irank , COMM );
        }
}

