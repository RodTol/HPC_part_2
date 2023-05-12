#include "headers/data_distr_utility.h"

/**
 * @brief This function computes the # of rows for each processor and
 * store it in array. Instead of the dimension of the matrix, this 
 * function needs already dimension/n_proc_tot and the  rest, 
 * dimension % n_proc_tot
 * 
 * @param n_rows_local The array made to store the values 
 * @param n_loc dimension/n_proc_tot
 * @param rest dimension % n_proc_tot 
 * @param n_proc_tot The # of processor
 */
void calculate_n_rows(int *n_rows_local, const int n_loc, const int rest, const int n_proc_tot) {    
    for (int i = 0; i < n_proc_tot; i++) {
        if (i < rest) {
            n_rows_local[i] = n_loc + 1;
        } else {
            n_rows_local[i] = n_loc;
        }
    }
}

/**
 * @brief This function given the # of processors,
 * compute the offset for each of them 
 * 
 * @param displacement The array where to store all the offsets
 * @param n_rows_local Array with the # of rows for each process
 * @param n_proc_tot # of total process
 */
void calculate_displ(int* displacement, int* n_rows_local, int n_proc_tot) {
    displacement[0] = 0;
    for (int i = 0; i < n_proc_tot-1; i++) {
        displacement[i+1] = displacement[i] + n_rows_local[i];
    }
}

/**
 * @brief This function computes the total number of elements
 * that process[count] should send to process[i]
 * 
 * @param n_elements_local Array designed to store the # of elements
 * @param n_rows_local Array with the # of rows for each process
 * @param count Array to select for which processor we are computing
 * @param n_proc_tot # of total process
 */
void calculate_n_elements(int *n_elements_local, const int *n_rows_local, const int count, const int n_proc_tot) {
    for (int i = 0; i < n_proc_tot; i++) {
        n_elements_local[i]  = n_rows_local[i] * n_rows_local[count];
    }
    
}

/**
 * @brief This function computes the offsets needed to access
 * the columns matrices in the correct way and store them in
 * an array
 * 
 * @param displacement_col Array of the offsets 
 * @param n_elements_count # of total elements
 * @param n_proc_tot # of total process
 */
void calculate_displ_col(int* displacement_col, int* n_elements_count, int n_proc_tot) {
    displacement_col[0] = 0;
    for (int i = 0; i < n_proc_tot-1; i++) {
        displacement_col[i+1] = displacement_col[i] + n_elements_count[i];
    }
}

/**
 * @brief This function encapsulates all the MPI_calls to build the 
 * column matrix needed for the multiplication
 * 
 * @param n_rows_local array with all the # of rows for each process
 * @param displacement array of the offsets
 * @param n_elements_local array with the number of elements for this
 * column matrix
 * @param displacement_col array with the offsets for the column matrix 
 * @param B The multiplier matrix
 * @param B_col The column matrix
 * @param irank Rank of the process
 * @param count Rank of the process whose storing the multiplier
 * @param N Size of B
 */
void MPI_build_column(int* n_rows_local, int* displacement, int* n_elements_local, int* displacement_col,
    double* B, double* B_col, int irank, int count, int N) {
    
    MPI_Datatype my_column_block;
    MPI_Type_vector(n_rows_local[irank], n_rows_local[count], N, MPI_DOUBLE, &my_column_block);
    MPI_Type_commit(&my_column_block);

    MPI_Allgatherv(B + displacement[count], 1, my_column_block, B_col, n_elements_local, displacement_col, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Type_free(&my_column_block);

}