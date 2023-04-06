#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "general_utility.h"

void array_of_random_doubles(double * A, int dim)
{
   double max = 4.0;
   double min = 0.0;
   double div = RAND_MAX/(max-min);

   for (int i = 0; i < dim; ++i) {
      srand(i + time(NULL));
      A[i] =  min + (rand() / div);
   }
}
