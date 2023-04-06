#include <stdlib.h>
#include <stdio.h>

#include "general_utility.h"

void vector_of_random_ints(int * a, int dim)
{
   int i;
   for (i = 0; i < dim; ++i) {
    a[i] = rand()%10;
   }
}