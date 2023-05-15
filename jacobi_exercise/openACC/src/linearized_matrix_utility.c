#include "headers/linearized_matrix_utility.h"
#include "headers/general_utility.h"
/**
 * @brief This function creates a local identity matrix
 * of size dim X dim 
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
 * @brief This function print a local matrix of dimension 
 * dim1 X dim2
 * 
 * @param A matrix
 * @param dim_1 rows dimension
 * @param dim_2 cols dimension
 */
void print_matrix(double * A, int dim_1, int dim_2, bool ghost) {
  
  /*Top ghost*/
  if (ghost) printf_red();
  for (int j = 0; j < dim_2; j++) {
    printf("%.3f ", A[linear_index(0,j,dim_1,dim_2)]);
  }
  printf("\n");
  if (ghost) printf_reset();

  /*Middle*/
  for(int i = 1; i < dim_1-1; i++ ){
    for(int j = 0; j < dim_2; j++ ){
      if (j==0 || j==dim_2-1) {
        if (ghost) printf_red();
      }
      printf("%.3f ", A[linear_index(i,j,dim_1,dim_2)]);
      if (j==0 || j==dim_2-1) {
        if (ghost) printf_reset();
      }
    }
    printf("\n");
  }

  /*Bottom ghost*/
  if (ghost) printf_red();
  for (int j = 0; j < dim_2; j++) {
    printf("%.3f ", A[linear_index(dim_1-1,j,dim_1,dim_2)]);
  }
  printf("\n");
  if (ghost) printf_reset();


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
 int* dim_1 , int dim_2, int n_proc_tot, MPI_Comm COMM, bool divisor) {
  
  if ( irank == 0 ) {
      /*This variable are only for irank==0
      so they can be moved inside the if*/
      double * A_tmp;
      print_matrix ( A , dim_1[irank], dim_2, true);
      if (divisor)     printf("-----------------------------------\n");
      for (int count = 1; count < n_proc_tot ; count ++ ) {
          int size= dim_1[count] * dim_2 * sizeof( double);
          A_tmp = (double *) malloc( size );

          MPI_Recv ( A_tmp , dim_1[count] * dim_2 , MPI_DOUBLE , count ,
            count , COMM , MPI_STATUS_IGNORE ) ;
          print_matrix ( A_tmp , dim_1[count], dim_2, true) ;
          if (divisor) printf("-----------------------------------\n");
      }
  }
  else {
      MPI_Send ( A , dim_1[irank] * dim_2 , MPI_DOUBLE , 0 ,
        irank , COMM );
  }
}

/**
 * @brief This function is used to print the full matrix, 
 * scattered among n_proc_tot processor, on a file on a serial way.
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
void print_matrix_distributed_file (double * A, int irank,
 int* dim_1 , int dim_2, int* displacement, int n_proc_tot, MPI_Comm COMM, char filename[]) {
  
  FILE *file;
  const double h = 0.1;

  /*NOTE: I want only the actual matrix, without the ghost layer*/
  if ( irank == 0 ) {
      file = fopen( filename, "w" );
      /*First proc skip last row(ghost layer)*/
      for(int i = 0; i < dim_1[0]-1; ++i ) {
         for(int j = 0; j < dim_2 ; ++j ) {
          fprintf(file, "%f\t%f\t%f\n", h * (j), -h * (i), A[linear_index(i,j,dim_1[0],dim_2)] );
        }
      }

      double * A_tmp;

      for (int count = 1; count < n_proc_tot-1 ; count ++ ) {
        size_t size= dim_1[count] * dim_2 * sizeof( double);
        A_tmp = (double *) malloc( size );
        MPI_Recv ( A_tmp , dim_1[count] * dim_2 , MPI_DOUBLE , count ,
          count , COMM , MPI_STATUS_IGNORE );
        /*All other proc skips last and first rows*/
        for(int i = 1; i < dim_1[count]-1; ++i ) {
          for(int j = 0; j < dim_2; ++j ) {
            fprintf(file, "%f\t%f\t%f\n", h * (j), -h * (i-1+displacement[count]), A_tmp[linear_index(i,j,dim_1[count],dim_2)] );
          }
        }
      }
      size_t size= dim_1[n_proc_tot-1] * dim_2 * sizeof( double);
      A_tmp = (double *) malloc( size );
      MPI_Recv ( A_tmp , dim_1[n_proc_tot-1] * dim_2 , MPI_DOUBLE , n_proc_tot-1 ,
        n_proc_tot-1 , COMM , MPI_STATUS_IGNORE );
      /*Last proc skips first row*/
      for(int i = 1; i < dim_1[n_proc_tot-1]; ++i ) {
        for(int j = 0; j < dim_2; ++j ) {
          fprintf(file, "%f\t%f\t%f\n", h * (j), -h * (i-1+displacement[n_proc_tot-1]), A_tmp[linear_index(i,j,dim_1[n_proc_tot-1],dim_2)] );
        }
      }
      fclose( file );
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
 int dim_1 , int dim_2, int offset) {
  int j_glob = 0;
  memset ( A , 0 , dim_1 * dim_2 * sizeof ( double ) ) ;
  for (int i_loc = 0; i_loc < dim_1 ; i_loc ++ ) {
    j_glob = i_loc + offset ;
    A [ j_glob + ( i_loc * dim_2 ) ] = 1.0;
  }
}

void create_jacobi_start_distributed (double * A, int irank,
 int *dim_1 , int dim_2, int* offset, int n_proc_tot) {
  
  /*Firstly I set everything to 0.5*/
  for (int i = 0; i < dim_1[irank]; i++) {
    for (int j = 0; j < dim_2; j++) {
      A[linear_index(i,j,dim_1[irank], dim_2)] = 0.5;
    }
  }

  /*Top boundary is set to 0*/
  for (int j = 0; j < dim_2; j++) {
    if (irank==0) {
      A[linear_index(0,j,dim_1[irank],dim_2)] = 0;
    }
  }

  /*Right boundary is set to o*/
  for (int i = 0; i < dim_1[irank]; i++) {
    A[linear_index(i,dim_2-1,dim_1[irank],dim_2)] = 0;
  }
  
  double increment = 100.0 / (dim_2-1);
  /*Last row*/
  if (irank==n_proc_tot-1) {
    for (int j = 0; j < dim_2; j++) {
      A[linear_index(dim_1[irank]-1,j,dim_1[irank],dim_2)] = (dim_2-j-1)*increment;
    }
  }
  
  /*Left column*/
  if (irank==0) {
    A[linear_index(0,0,dim_1[irank],dim_2)] = 0;
  }
  /*I don't need to touch the bottom left corner because
  is already setted up by the row initialisation*/

  for (int i = 1; i < dim_1[irank]-1; i++) {
    int displacement = offset[irank];
    if (irank!=0) {
      displacement = displacement -1; 
    }
    A[linear_index(i,0,dim_1[irank],dim_2)] = (i+displacement)*increment;
  }

}

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
