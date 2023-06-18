
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <unistd.h>
#include <stdbool.h>


void printf_red ();
void printf_yellow ();
void printf_green();
void printf_reset ();

bool file_exists(const char *filename);

void create_null_array (double * A, int dim);
void array_of_random_doubles(double * A, int dim);
double randfrom(double min, double max);