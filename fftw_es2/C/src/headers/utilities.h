
#ifndef _FFTW_UTLITIES_
#define _FFTW_UTLITIES_
#include <complex.h>
#include <sys/time.h>
#include <fftw3.h>
#include <mpi.h>
#include <stdbool.h>
#define pi 3.14159265358979323846
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

/**
 * @brief This struct is created to store and save all the parameters
 * and information that will be used during the execution
 * 
 */
typedef struct {

  //fftw plans
  fftw_plan fw_plan_2d; 
  fftw_plan fw_plan_1d; 
  fftw_plan bw_plan_2d;
  fftw_plan bw_plan_1d;

  fftw_complex *data;
  fftw_complex *data_redistributed;
  
  ptrdiff_t global_size_grid;
  ptrdiff_t local_size_grid;
  ptrdiff_t all_to_all_block_size;
  //First data distribution parameters
  ptrdiff_t local_n1;
  ptrdiff_t local_n1_offset;
  //Second data distribution parameters
  ptrdiff_t local_n2;
  ptrdiff_t local_n2_offset;
  // Dimension of the grid
  ptrdiff_t n1;
  ptrdiff_t n2;
  ptrdiff_t n3;
  /*Parameters for all to all communication*/
  MPI_Datatype *send_type;
  MPI_Datatype *recv_type;
  int *send_n_of_blocks;
  int *recv_n_of_blocks;
  int *send_displacement;
  int *recv_displacement;

  MPI_Comm mpi_comm;  
  
} fftw_dist_handler;



double seconds();
int index_f ( int i1, int i2, int i3, int n1, int n2, int n3 );

void plot_data_1d( char* name, int n1, int n2, int n3, int dir, double* data );
void plot_data_2d( char* name, int n1, int n2, int n3, int n1_local,
 int n1_local_offset, int dir, double* data);
void init_fftw( fftw_dist_handler* fft, int n1, int n2, int n3, MPI_Comm comm );
void close_fftw( fftw_dist_handler* fft );

void derivative( fftw_dist_handler* fft,int n1, int n2, int n3, double L1, double L2, double L3, int ipol, double* data, double* deriv );
void fft_3d( fftw_dist_handler* fft, double *data_direct, fftw_complex* data_rec, bool direct_to_reciprocal );

#endif
