#ifndef _FFTW_UTLITIES_
#define _FFTW_UTLITIES_
#include <complex.h>
#include <sys/time.h>
#include <fftw3.h>
#include <mpi.h>
#include <fftw3-mpi.h>
#include <stdbool.h>
#include <unistd.h>

#define pi 3.14159265358979323846

/**
 * @brief This struct is created to store and save all the parameters
 * and information that will be used during the execution
 * 
 */
typedef struct {

  fftw_plan fw_plan; 
  fftw_plan bw_plan;
  fftw_complex *fftw_data;
  ptrdiff_t global_size_grid; /* Global size of the FFT grid */
  ptrdiff_t local_size_grid;  /* Local size of the FFT grid */
  ptrdiff_t local_n1;         /* Local dimension of n1 */
  ptrdiff_t local_n1_offset;  /* Offset due to possible rests */
  MPI_Comm mpi_comm;
  
} fftw_mpi_handler;



double seconds();
int index_f ( int i1, int i2, int i3, int n1, int n2, int n3 );


void plot_data_1d( char* name, int n1, int n2, int n3, int n1_local, int n1_local_offset, int dir, double* data );
void plot_data_2d( char* name, int n1, int n2, int n3, int n1_local, int n1_local_offset, int dir, double* data );
void init_fftw( fftw_mpi_handler* fft, int n1, int n2, int n3, MPI_Comm mpi_comm );
void close_fftw( fftw_mpi_handler* fft );

void derivative( fftw_mpi_handler* fft,int n1, int n2, int n3, double L1, double L2, double L3, int ipol, double* data, double* deriv );
void fft_3d( fftw_mpi_handler* fft, double *data_direct, fftw_complex* data_rec, bool direct_to_reciprocal );

#endif
