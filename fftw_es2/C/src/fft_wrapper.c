/* Assignement:
 * Here you have to modify the includes, the array sizes and the fftw calls, to use the fftw-mpi
 *
 * Regarding the fftw calls. here is the substitution 
 * fftw_plan_dft_3d -> fftw_mpi_plan_dft_3d
 * ftw_execute_dft  > fftw_mpi_execute_dft 
 * use fftw_mpi_local_size_3d for local size of the arrays
 * 
 * Created by G.P. Brandino, I. Girotto, R. Gebauer
 * Last revision: March 2016
 *
 */ 

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
 *  For more informtion see FFTW documentation:
 *  http://www.fftw.org/doc/Row_002dmajor-Format.html#Row_002dmajor-Format
 *
 */
int index_f ( int i1, int i2, int i3, int n1, int n2, int n3){
  return n3*n2*i1 + n3*i2 + i3; 
}

void init_fftw(fftw_dist_handler *fft, int n1, int n2, int n3, MPI_Comm comm){
  
  int n_proc_tot, irank;
  MPI_Comm_size( comm, &n_proc_tot );
  MPI_Comm_rank( comm, &irank );
  fft->mpi_comm = comm;

  /*
   *  Allocate a distributed grid for complex FFT using aligned memory allocation
   *  See details here:
   *  http://www.fftw.org/fftw3_doc/Allocating-aligned-memory-in-Fortran.html#Allocating-aligned-memory-in-Fortran
   *  HINT: initialize all global and local dimensions. Consider the first dimension being multiple of the number of processes
   *
   */
  /*Già nella consegna, non vuole i resti!!*/
  if( ( ( n1 % n_proc_tot ) || ( n2 % n_proc_tot ) ) && !irank ){    
    fprintf( stdout, "\nN1 dimension must be multiple of the number of processes. The program will be aborted...\n\n" );
    MPI_Abort( comm, 1 );
  }

  /* Fill the missing parts */
  fft->n1 = n1;
  fft->n2 = n2;
  fft->n3 = n3;

  /*Come prima direzione di divisione, scegliamo 1.
  Quindi ogni proc ha n1_loc*n2*n3 dati.
  Ci sarebbe la funzione che calcola le cose, ma 
  facciamo a mano per garantire portabilità */
  fft->local_n1 = n1/n_proc_tot; 
  fft->local_n1_offset = fft->local_n1*irank;

  /*Seconda divisione dati*/
  fft->local_n2 = n2/n_proc_tot; 
  fft->local_n2_offset = fft->local_n2*irank;

  fft->global_size_grid = n1*n2*n3;
  /*Dimensione della slab. NOTA: questa è unica perchè
  il numero di dati è sempre lo stesso, sia quando taglio
  in direzione 1 che in 2. Cambia ovviamente come li leggo*/
  fft->local_size_grid = fft->local_n1*n2*n3;
  /*Dimensione del volume da mandare nella all to all*/
  fft->all_to_all_block_size = fft->local_n1*fft->local_n2*n3;

  fft->data = ( fftw_complex* ) fftw_malloc( fft->local_size_grid * sizeof( fftw_complex ) );
  fft->data_redistributed = ( fftw_complex* ) fftw_malloc( fft->local_size_grid * sizeof( fftw_complex ) );
  /*
   * Allocate fft->fftw_data and create an FFTW plan for each 1D FFT among all dimensions
   *
   * Questo sarebbe il punto 1, noi andiamo direttamente al punto 2
   * quindi ho solo una piano 2D e un piano 1D.
   * La fft sarà fatta prima sul piano 2-3 e poi sulla singola
   * direzione 1 
   */
  /*Per usare l'interfaccia avanzata, devo passare un vettore
  con le dimensioni sui cui fare la fft*/
  int dims[] = {n2,n3};
  int dim_1[] = {n1};
  /*Usiamo l'interfaccia avanzata ma assumendo i dati contigui*/
  fft->fw_plan_i1 = fftw_plan_many_dft(2, dims, fft->local_n1,
    fft->data, dims, 1, fft->n2*fft->n3, fft->data, 
    dims, 1, fft->n2*fft->n3, FFTW_FORWARD, FFTW_ESTIMATE);

  /*Qua posso farlo come se avessi i dati contigui, ma dovr
  implementare il riordine. Oppure, introdurre già lo strife*/
  fft->fw_plan_i2 = fftw_plan_many_dft(1, dim_1, fft->local_n2*fft->n3,
    fft->data_redistributed, dim_1, fft->local_n2*fft->n3, 1, fft->data_redistributed,
    dim_1, fft->local_n2*fft->n3, 1, FFTW_FORWARD, FFTW_ESTIMATE);
  //fft->fw_plan_i3 = NULL;

  fft->bw_plan_i1 = fftw_plan_many_dft(2, dims, fft->local_n1,
    fft->data, dims, 1, fft->n2*fft->n3, fft->data, 
    dims, 1, fft->n2*fft->n3, FFTW_BACKWARD, FFTW_ESTIMATE);
  fft->bw_plan_i2 = fftw_plan_many_dft(1, dim_1, fft->local_n2*fft->n3,
    fft->data_redistributed, dim_1, fft->local_n2*fft->n3, 1, fft->data_redistributed,
    dim_1, fft->local_n2*fft->n3, 1, FFTW_BACKWARD, FFTW_ESTIMATE);
  //fft->bw_plan_i3 = NULL;

  MPI_Datatype column_block;
  MPI_Type_vector(fft->local_n1, fft->local_n2*n3, n2*n3, MPI_C_DOUBLE_COMPLEX, &column_block);
  MPI_Type_commit(&column_block);

  fft->send_n_of_blocks = (int *)malloc(n_proc_tot * sizeof(int));
  fft->recv_n_of_blocks = (int *)malloc(n_proc_tot * sizeof(int));
  fft->send_displacement = (int *)malloc(n_proc_tot * sizeof(int));
  fft->recv_displacement = (int *)malloc(n_proc_tot * sizeof(int));
  fft->send_type = (MPI_Datatype *)malloc(n_proc_tot * sizeof(MPI_Datatype));
  fft->recv_type = (MPI_Datatype *)malloc(n_proc_tot * sizeof(MPI_Datatype));
  
  for (int i = 0; i < n_proc_tot; i++) {
    fft->send_type[i] = column_block;
    fft->recv_type[i] = MPI_C_DOUBLE_COMPLEX;
    fft->send_n_of_blocks[i] = 1;
    fft->recv_n_of_blocks[i] = fft->all_to_all_block_size;
  }

  fft->send_displacement[0] = 0*sizeof(fftw_complex);
  fft->recv_displacement[0] = 0*sizeof(fftw_complex);
  for (int i = 0; i < n_proc_tot-1; i++)
  {
    fft->send_displacement[i+1] = fft->send_displacement[i] + fft->local_n2*n3*sizeof(fftw_complex);
    fft->recv_displacement[i+1] = fft->recv_displacement[i] + fft->all_to_all_block_size*sizeof(fftw_complex); 
  }

}

void close_fftw( fftw_dist_handler *fft ){

    fftw_destroy_plan( fft->bw_plan_i1 );
    fftw_destroy_plan( fft->bw_plan_i2 );
    //fftw_destroy_plan( fft->bw_plan_i3 );

    fftw_destroy_plan( fft->fw_plan_i1 );
    fftw_destroy_plan( fft->fw_plan_i2 );
    //fftw_destroy_plan( fft->fw_plan_i3 );

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

/* This subroutine uses fftw to calculate 3-dimensional discrete FFTs.
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
 */
void fft_3d( fftw_dist_handler* fft, double *data_direct, fftw_complex* data_rec,
 bool direct_to_reciprocal ){

  double fac;
  int n_proc_tot;
  int i1, i2, i3;
  int n1 = fft->n1, n2 = fft->n2, n3 = fft->n3;
  int local_size_grid = fft->local_size_grid;
  fftw_complex * data = fft->data;
  fftw_complex * data_redistributed = fft->data_redistributed;

  MPI_Comm_size( fft->mpi_comm, &n_proc_tot );
    
  // Now distinguish in which direction the FFT is performed
  if( direct_to_reciprocal ) {

    /*Rendo i valori complessi (per la memoria ?)*/
    for(int i = 0; i < local_size_grid; i++) {
      data[i]  = data_direct[i] + 0.0 * I;
    } 

    /*I perform the first fft on the n2-n3 plan locally*/
    fftw_execute(fft->fw_plan_i1);
     
    //Perform all_to_all
    MPI_Alltoallw(data, fft->send_n_of_blocks, fft->send_displacement, fft->send_type,
     data_redistributed, fft->recv_n_of_blocks, fft->recv_displacement, fft->recv_type, MPI_COMM_WORLD);
    
    //Perform fft on n1 direction
    fftw_execute(fft->fw_plan_i2);

    // Perform an Alltoall communication 
    MPI_Alltoallw(data_redistributed, fft->recv_n_of_blocks, fft->recv_displacement, fft->recv_type,
     data, fft->send_n_of_blocks, fft->send_displacement, fft->send_type, MPI_COMM_WORLD);    
    
    /*Copy the data into the data_rec array*/
    memcpy(data_rec, data, fft->local_size_grid*sizeof(fftw_complex));
  } else {
    /*Copio i dati da cui farò l'inversa in fft->data*/
    memcpy(data, data_rec, fft->local_size_grid*sizeof(fftw_complex));

    /* Implement the reverse transform */
    fftw_execute(fft->bw_plan_i1);
     
    MPI_Alltoallw(data, fft->send_n_of_blocks, fft->send_displacement, fft->send_type,
     data_redistributed, fft->recv_n_of_blocks, fft->recv_displacement, fft->recv_type, MPI_COMM_WORLD);
    
    fftw_execute(fft->bw_plan_i2);

    MPI_Alltoallw(data_redistributed, fft->recv_n_of_blocks, fft->recv_displacement, fft->recv_type,
     data, fft->send_n_of_blocks, fft->send_displacement, fft->send_type, MPI_COMM_WORLD);
    
    /*Normalizzo sulla globale*/
    fac = 1.0 / ( fft->global_size_grid );

    for(int i = 0; i < fft->local_size_grid; ++i ) {
      data_direct[i] = creal(data[i])*fac;
    }

  }
  
  
}

