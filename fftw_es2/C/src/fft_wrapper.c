#include <string.h>
#include <stdlib.h>
#include "headers/utilities.h"


double seconds(){
/* Return the second elapsed since Epoch (00:00:00 UTC, January 1, 1970) */
  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
  return sec;
}

/* 
 *  Index linearization is computed following row-major order.
 */
int index_f ( int i1, int i2, int i3, int n1, int n2, int n3){
  return n3*n2*i1 + n3*i2 + i3; 
}


/**
 * @brief This function initialize the fft handler with all the parameters and
 * variables that will be used during the execution
 * 
 * @param fft the handler
 * @param n* grid size
 * @param comm MPI communicator
 */
void init_fftw(fftw_dist_handler *fft, int n1, int n2, int n3, MPI_Comm comm){
  
  // I initialize again the MPI environment
  int n_proc_tot, irank;
  MPI_Comm_size( comm, &n_proc_tot );
  MPI_Comm_rank( comm, &irank );
  fft->mpi_comm = comm;

  // I want a symmetric distribution of the data among the n_proc_tot processes
  if( ( ( n1 % n_proc_tot ) || ( n2 % n_proc_tot ) ) && !irank ){    
    fprintf( stdout, "\nN1 dimension must be multiple of the number of processes. The program will be aborted...\n\n" );
    MPI_Abort( comm, 1 );
  }

  fft->n1 = n1;
  fft->n2 = n2;
  fft->n3 = n3;

  // I compute the local dimension of the grid when the matrix is distributed
  // along the n1 dimension. Usally there's a fftw_mpi routine for this but here
  // we will do it by hand
  fft->local_n1 = n1/n_proc_tot; 
  fft->local_n1_offset = fft->local_n1*irank;

  // Same but along the n2 dimension
  fft->local_n2 = n2/n_proc_tot; 
  fft->local_n2_offset = fft->local_n2*irank;

  fft->global_size_grid = n1*n2*n3;
  fft->local_size_grid = fft->local_n1*n2*n3;
  // This is the dimensio of the volume that the all_to_all routine will send
  fft->all_to_all_block_size = fft->local_n1*fft->local_n2*n3;

  fft->data = ( fftw_complex* ) fftw_malloc( fft->local_size_grid * sizeof( fftw_complex ) );
  fft->data_redistributed = ( fftw_complex* ) fftw_malloc( fft->local_size_grid * sizeof( fftw_complex ) );

  /**
   * We create the fft plans using the fftw advanced interface. Both in the distribution
   * along the n1 and n2 dimension, the data are contiguos. The first plan will
   * work on the n2 and n3 dimension while the last one only on n1.
   * The advanced interface need some extra parameters: howmany, inembed, istridem onemebed,
   * ostride, odist
   * 
   */

  // The advanced interface requires an array of dimension on which perform the fft
  int dims[] = {n2,n3};
  int dim_1[] = {n1};
  // First the 2 dimensional plan
  fft->fw_plan_2d = fftw_plan_many_dft(2, dims, fft->local_n1,
    fft->data, dims, 1, fft->n2*fft->n3, fft->data, 
    dims, 1, fft->n2*fft->n3, FFTW_FORWARD, FFTW_ESTIMATE);

  // Then the 1 dimensional
  fft->fw_plan_1d = fftw_plan_many_dft(1, dim_1, fft->local_n2*fft->n3,
    fft->data_redistributed, dim_1, fft->local_n2*fft->n3, 1, fft->data_redistributed,
    dim_1, fft->local_n2*fft->n3, 1, FFTW_FORWARD, FFTW_ESTIMATE);

  //For both plan we also need the inverse transformation, so we change the sign
  fft->bw_plan_2d = fftw_plan_many_dft(2, dims, fft->local_n1,
    fft->data, dims, 1, fft->n2*fft->n3, fft->data, 
    dims, 1, fft->n2*fft->n3, FFTW_BACKWARD, FFTW_ESTIMATE);
  fft->bw_plan_1d = fftw_plan_many_dft(1, dim_1, fft->local_n2*fft->n3,
    fft->data_redistributed, dim_1, fft->local_n2*fft->n3, 1, fft->data_redistributed,
    dim_1, fft->local_n2*fft->n3, 1, FFTW_BACKWARD, FFTW_ESTIMATE);

  // We create an MPI_Datatype to select the columns of the local array
  MPI_Datatype column_block;
  MPI_Type_vector(fft->local_n1, fft->local_n2*n3, n2*n3, MPI_C_DOUBLE_COMPLEX, &column_block);
  MPI_Type_commit(&column_block);

  // Parameters for the  MPI_Alltoallw
  fft->send_n_of_blocks = (int *)malloc(n_proc_tot * sizeof(int));
  fft->recv_n_of_blocks = (int *)malloc(n_proc_tot * sizeof(int));
  fft->send_displacement = (int *)malloc(n_proc_tot * sizeof(int));
  fft->recv_displacement = (int *)malloc(n_proc_tot * sizeof(int));
  fft->send_type = (MPI_Datatype *)malloc(n_proc_tot * sizeof(MPI_Datatype));
  fft->recv_type = (MPI_Datatype *)malloc(n_proc_tot * sizeof(MPI_Datatype));
  
  for (int i = 0; i < n_proc_tot; i++) {
    fft->send_type[i] = column_block;
    fft->recv_type[i] = MPI_C_DOUBLE_COMPLEX; //equivalent of fftw_complex
    fft->send_n_of_blocks[i] = 1;
    fft->recv_n_of_blocks[i] = fft->all_to_all_block_size;
  }

  // We select what data each array has to send and where to put the data that are received
  fft->send_displacement[0] = 0*sizeof(fftw_complex);
  fft->recv_displacement[0] = 0*sizeof(fftw_complex);
  for (int i = 0; i < n_proc_tot-1; i++)
  {
    fft->send_displacement[i+1] = fft->send_displacement[i] + fft->local_n2*n3*sizeof(fftw_complex);
    fft->recv_displacement[i+1] = fft->recv_displacement[i] + fft->all_to_all_block_size*sizeof(fftw_complex); 
  }

}

/**
 * @brief This function is meant to deallocate and close the fftw process
 * 
 * @param fft 
 */
void close_fftw( fftw_dist_handler *fft ){

    fftw_destroy_plan( fft->bw_plan_2d );
    fftw_destroy_plan( fft->bw_plan_1d );

    fftw_destroy_plan( fft->fw_plan_2d );
    fftw_destroy_plan( fft->fw_plan_1d );

    fftw_free( fft->data );
    fftw_free( fft->data_redistributed);

    free(fft->send_n_of_blocks);
    free(fft->recv_n_of_blocks);
    free(fft->send_displacement); 
    free(fft->recv_displacement);
    free(fft->send_type);
    free(fft->recv_type);

    fftw_cleanup();
}

/**
 * @brief This subroutine uses fftw to calculate 3-dimensional discrete FFTs.
 * The data in direct space is assumed to be real-valued
 * The data in reciprocal space is complex. 
 * direct_to_reciprocal indicates in which direction the FFT is to be calculated
 * 
 * Note that for real data in direct space (like here), we have
 * F(N-j) = conj(F(j)) where F is the array in reciprocal space.
 * Here, we do not make use of this property.
 * Also, we do not use the special (time-saving) routines of FFTW which
 * allow one to save time and memory for such real-to-complex transforms.
 *
 * f: array in direct space
 * F: array in reciprocal space
 * 
 * F(k) = \sum_{l=0}^{N-1} exp(- 2 \pi I k*l/N) f(l)
 * f(l) = 1/N \sum_{k=0}^{N-1} exp(+ 2 \pi I k*l/N) F(k)
 * 
 * @param fft the handler with all the information
 * @param data_direct the real data
 * @param data_rec the complex data
 * @param direct_to_reciprocal the direction of the fft
 */
void fft_3d( fftw_dist_handler* fft, double *data_direct, fftw_complex* data_rec, bool direct_to_reciprocal ){
  double fac;
  int local_size_grid = fft->local_size_grid;
  fftw_complex * data = fft->data;
  fftw_complex * data_redistributed = fft->data_redistributed;
  
  //I need the number of process
  int n_proc_tot;
  MPI_Comm_size( fft->mpi_comm, &n_proc_tot );
    
  // Now distinguish in which direction the FFT is performed
  if( direct_to_reciprocal ) {
    // Since fft->data is fftw_complex, we need to make data_direct complex
    for(int i = 0; i < local_size_grid; i++) {
      data[i]  = data_direct[i] + 0.0 * I;
    } 

    //I perform the first fft on the n2-n3 plan locally
    fftw_execute(fft->fw_plan_2d);
     
    //Perform all_to_all to have the data distributed along n2 direction
    MPI_Alltoallw(data, fft->send_n_of_blocks, fft->send_displacement, fft->send_type,
     data_redistributed, fft->recv_n_of_blocks, fft->recv_displacement, fft->recv_type, MPI_COMM_WORLD);
    
    //Perform fft on n1 direction
    fftw_execute(fft->fw_plan_1d);

    // Perform an Alltoall communication to get the data to the original ordering
    MPI_Alltoallw(data_redistributed, fft->recv_n_of_blocks, fft->recv_displacement, fft->recv_type,
     data, fft->send_n_of_blocks, fft->send_displacement, fft->send_type, MPI_COMM_WORLD);    
    
    //Copy the data into the data_rec array
    memcpy(data_rec, data, fft->local_size_grid*sizeof(fftw_complex));
  } else {
    //Copy the complex data_rec into the data array
    memcpy(data, data_rec, fft->local_size_grid*sizeof(fftw_complex));

    //Perform the reverse transform on n2 and n3
    fftw_execute(fft->bw_plan_2d);
    
    //Perform all_to_all to have the data distributed along n2 direction
    MPI_Alltoallw(data, fft->send_n_of_blocks, fft->send_displacement, fft->send_type,
     data_redistributed, fft->recv_n_of_blocks, fft->recv_displacement, fft->recv_type, MPI_COMM_WORLD);
    
    //Perform the reverse transform on n1
    fftw_execute(fft->bw_plan_1d);

    // Perform an Alltoall communication to get the data to the original ordering
    MPI_Alltoallw(data_redistributed, fft->recv_n_of_blocks, fft->recv_displacement, fft->recv_type,
     data, fft->send_n_of_blocks, fft->send_displacement, fft->send_type, MPI_COMM_WORLD);
    
    //I normalize the data
    fac = 1.0 / ( fft->global_size_grid );
    for(int i = 0; i < fft->local_size_grid; ++i ) {
      data_direct[i] = creal(data[i])*fac;
    }

  }
  
  
}

