#include "headers/linearized_matrix_utility.h"

/**
 * @brief This function creates a local identity
 * matrix of size dim X dim, made of doubles
 * 
 * @param A the matrix
 * @param dim the size of the matrix
*/
void create_identity_matrix (double * A, int dim) {
  memset( A, 0, dim * dim * sizeof(double) );
  for (int i = 0; i < dim; i++){
    A[i + i*dim] = 1.0;
  }
  
}

/**
 * @brief This function prints a local matrix
 * of size dim X dim, made of doubles
 * 
 * @param A the matrix
 * @param dim the size of the matrix
*/
void print_matrix_square(double * A, int dim ){
  fprintf( stdout, "\n");
  for(int i = 0; i < dim; i++ ){
    for(int j = 0; j < dim; j++ ){
      fprintf( stdout, "%.3g ", A[ j + ( i * dim ) ] );
    }
    fprintf( stdout, "\n");
  }
}

/**
 * @brief This function return the linear index for 
 * a matrix of dimension dim1xdim2 given the row index i
 * and the column index j
 * 
 * @param i row index
 * @param j column index
 * @param dim_1 rows dimension
 * @param dim_2 cols dimension
*/
int linear_index ( int i, int j, int dim1, int dim2)
{
  return dim2*i + j; 
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
      fprintf( stdout, "%.3g ", A[ j + ( i * dim_2 ) ] );
    }
    fprintf( stdout, "\n");
  }
}

/**
 * @brief This function is used to print the full matrix, 
 * scattered among n_proc_tot processor.
 * @note this is a NON SCALABLE function because
 * the copy on a buffer requires that the matrix doesn't
 * occuy more than half on a single node.
 * 
 * @param A the matrix
 * @param irank rank of each processor
 * @param dim_1 vector of #of rows for each submatrix
 * @param dim_2 cols of each submatrix
 * @param n_proc_tot number of total processors
 * @param COMM MPI communicator
*/
void print_matrix_distributed (double * A, int irank,
 int* dim_1 , int dim_2, int n_proc_tot, MPI_Comm COMM) {
  
  if ( irank == 0 ) {
      /*This variable are only for irank==0
      so they can be moved inside the if*/
      double * A_tmp;
      print_matrix ( A , dim_1[irank], dim_2) ;

      for (int count = 1; count < n_proc_tot ; count ++ ) {
          int size= dim_1[count] * dim_2 * sizeof( double);
          A_tmp = (double *) malloc( size );

          MPI_Recv ( A_tmp , dim_1[count] * dim_2 , MPI_DOUBLE , count ,
            count , COMM , MPI_STATUS_IGNORE ) ;
          print_matrix ( A_tmp , dim_1[count], dim_2 ) ;
      }
  }
  else {
      MPI_Send ( A , dim_1[irank] * dim_2 , MPI_DOUBLE , 0 ,
        irank , COMM );
  }
}

/**
 * @brief This function initialize a given NxN matrix,
 * distributed among n_processors, to the identity
 * 
 * @param A the matrix
 * @param irank rank of each processor
 * @param dim_1 rows of each submatrix
 * @param dim_2 cols of each submatrix
 * @param offset offset if rest!=0
*/
void create_identity_matrix_distributed (double * A, int irank,
 int dim_1 , int dim_2,  int offset) {
  int j_glob = 0;
  memset ( A , 0 , dim_1 * dim_2 * sizeof ( double ) ) ;
  for (int i_loc = 0; i_loc < dim_1 ; i_loc ++ ) {
    j_glob = i_loc + offset ;
    A [ j_glob + ( i_loc * dim_2 ) ] = 1.0;
  }
}


/**
 * @brief This function implements the naive computation of the
 * matrix multiplication
 * 
 * @param A The first matrix
 * @param B_col The transposed second matrix
 * @param C The result matrix 
 * @param N Size of the matrix
 * @param n_rows_local array that stores the # of local rows in each process
 * @param displacement array that stores the offsets of each process
 * @param irank rank of each process
 * @param count index of the multiplication
 */
void matrix_multiplication(double* A, double* B_col, double* C, 
  int N, int* n_rows_local, int* displacement, int irank, int count) {
  for (int i = 0; i < n_rows_local[irank]; i++) {       // A is n_loc x N (i,j)
      for (int k = 0; k < n_rows_local[count]; k++) {   // B_col is N x n_loc (j,k)
        for (int j = 0; j < N; j++) {   // C slice is n_loc x n_loc -> C is n_loc x N (same as A)
          C[i * N + k + displacement[count]] += A[i * N + j] * B_col[j* n_rows_local[count] + k];
        }
      }
  }
}