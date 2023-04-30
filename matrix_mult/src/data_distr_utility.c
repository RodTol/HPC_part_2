#include "headers/data_distr_utility.h"

void calculate_n_rows(int *n_rows_local,const int n_loc,const int rest, const int irank,const int n_proc_tot) {    
    for (int i = 0; i < n_proc_tot; i++) {
        if (i < rest) {
            n_rows_local[i] = n_loc + 1;
        } else {
            n_rows_local[i] = n_loc;
        }
    }
}