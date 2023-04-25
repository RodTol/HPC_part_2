#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

void create_identity_matrix (double * A, int dim);

void print_matrix_square(double * A, int dim );

int linear_index ( int i, int j, int dim1, int dim2);

void print_matrix(double * A, int dim_1, int dim_2);

void print_matrix_distributed (double * A, int irank,
 int dim_1 , int dim_2, int n_proc_tot, MPI_Comm COMM);

void create_identity_matrix_distributed (double * A, int irank,
 int dim_1 , int dim_2,  int offset, int n_proc_tot);

void matrix_multiplication();