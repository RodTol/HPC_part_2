#include <stdio.h>
#include <mpi.h>

void calculate_n_rows(int *n_rows_local, const int n_loc, const int rest, const int n_proc_tot);

void calculate_displ(int* displacement, int* n_rows_local, int n_proc_tot);

void calculate_displ_col(int* displacement, int* n_rows_local, int n_proc_tot);

void calculate_n_elements(int *n_elements_local, const int *n_rows_local, const int count, const int n_proc_tot);

void MPI_build_column(int* n_rows_local, int* displacement, int* n_elements_local, int* displacement_col,
    double* B, double* B_col, int irank, int count, int N);
