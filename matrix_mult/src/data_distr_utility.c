#include "headers/data_distr_utility.h"

void calculate_n_rows(int *n_rows_local, const int n_loc, const int rest, const int n_proc_tot) {    
    for (int i = 0; i < n_proc_tot; i++) {
        if (i < rest) {
            n_rows_local[i] = n_loc + 1;
        } else {
            n_rows_local[i] = n_loc;
        }
    }
}

void calculate_displ(int* displacement, int* n_rows_local, int n_proc_tot) {
    displacement[0] = 0;
    for (int i = 0; i < n_proc_tot-1; i++) {
        displacement[i+1] = displacement[i] + n_rows_local[i];
    }
}

void calculate_displ_col(int* displacement, int* n_rows_local, int n_proc_tot) {
    displacement[0] = 0;
    for (int i = 0; i < n_proc_tot-1; i++) {
        displacement[i+1] = displacement[i] + n_rows_local[i];
    }
}

void calculate_n_elements(int *n_elements_local, const int *n_rows_local, const int count, const int n_proc_tot) {
    for (int i = 0; i < n_proc_tot; i++) {
        n_elements_local[i]  = n_rows_local[i] * n_rows_local[count];
    }
    
}

void MPI_build_column(int* n_rows_local, int* displacement, int* n_elements_local, int* displacement_col,
    double* B, double* B_col, int irank, int count, int N) {
    
    MPI_Datatype my_column_block;
    MPI_Type_vector(n_rows_local[irank], n_rows_local[count], N, MPI_DOUBLE, &my_column_block);
    MPI_Type_commit(&my_column_block);

    MPI_Allgatherv(B + displacement[count], 1, my_column_block, B_col, n_elements_local, displacement_col, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Type_free(&my_column_block);

}