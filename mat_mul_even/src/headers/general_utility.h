
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

void printf_red ();
void printf_yellow ();
void printf_green();
void printf_reset ();

void create_null_array (double * A, int dim);
void array_of_random_doubles(double * A, int dim);