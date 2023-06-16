#include "headers/general_utility.h"

/**
 * @brief Color Red for the output of C program
 */
void printf_red () {
  printf("\033[1;31m");
}

/**
 * @brief Color Yellow for the output of C program
 */
void printf_yellow (){
  printf("\033[1;33m");
}

/**
 * @brief Color Green for the output of C program
 */
void printf_green (){
  printf("\033[1;32m");
}

/**
 * @brief Color reset for the output of C program
 */
void printf_reset () {
  printf("\033[0m");
}

/**
 * @brief  return true if the file specified by the
 *  filename exists
 * 
 * @param filename name of the file
*/

bool file_exists(const char *filename) {
    return access(filename, F_OK) == 0;
}

/**
 * @brief Create a null array object
 * 
 * @param A the array
 * @param dim length of the array
 */
void create_null_array (double * A, int dim) {
  memset( A, 0, dim * sizeof(double) );
}

/**
 * @brief Create an array of random doubles
 * 
 * @param A the array
 * @param dim length of the array
 */
void array_of_random_doubles(double * A, int dim)
{
   int id;
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   double max = 5.0;
   double min = 0.0;
   double div = RAND_MAX/(max-min);

   for (int i = 0; i < dim; ++i) {
      srand(i + id + time(NULL));
      A[i] =  min + (rand() / div);
#ifdef DEBUG      
      if (i % 5000==0 && i >=5000)
      {
        printf("\n Generated 5000 number \n");
      }
#endif      
   }
}
