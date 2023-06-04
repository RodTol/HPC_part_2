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

#include <complex.h>
#include <mpi.h>
#include <fftw3-mpi.h>
#include <fftw3.h>
#include <stdbool.h>
#include <string.h>
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

  //int buffer_size = 0;
  fft->mpi_comm = comm;

  /*
   *  Allocate a distributed grid for complex FFT using aligned memory allocation
   *  See details here:
   *  http://www.fftw.org/fftw3_doc/Allocating-aligned-memory-in-Fortran.html#Allocating-aligned-memory-in-Fortran
   *  HINT: initialize all global and local dimensions. Consider the first dimension being multiple of the number of processes
   *
   */
  /*Già nella consegna, non vuole i resti!!*/
  if( ( ( n1 % n_proc_tot ) || ( n2 % n_proc_tot ) ) && !irank ){    fprintf( stdout, "\nN1 dimension must be multiple of the number of processes. The program will be aborted...\n\n" );
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
  /*Usiamo l'interfaccia avanzata ma assumendo i dati contigui*/
  fft->fw_plan_i1 = fftw_plan_many_dft(2, dims, fft->local_n1,
    fft->data, dims, 1, fft->n2*fft->n1, fft->data, 
    dims, 1, fft->n2*fft->n3, FFTW_FORWARD, FFTW_ESTIMATE);

  /*Qua posso farlo come se avessi i dati contigui, ma dovr
  implementare il riordine. Oppure, introdurre già lo strife*/
  fft->fw_plan_i2 = NULL;
  //fft->fw_plan_i3 = NULL;

  fft->bw_plan_i1 = fftw_plan_many_dft(2, dims, fft->local_n1,
    fft->data, dims, 1, fft->n2*fft->n1, fft->data, 
    dims, 1, fft->n2*fft->n3, FFTW_BACKWARD, FFTW_ESTIMATE);;
  fft->bw_plan_i2 = NULL;
  //fft->bw_plan_i3 = NULL;

}

void close_fftw( fftw_dist_handler *fft ){

    fftw_destroy_plan( fft->bw_plan_i1 );
    fftw_destroy_plan( fft->bw_plan_i2 );
    //fftw_destroy_plan( fft->bw_plan_i3 );

    fftw_destroy_plan( fft->fw_plan_i1 );
    fftw_destroy_plan( fft->fw_plan_i2 );
    //fftw_destroy_plan( fft->fw_plan_i3 );

    fftw_free( fft->data );
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
// void fft_3d( fftw_dist_handler* fft, double *data_direct, fftw_complex* data_rec, bool direct_to_reciprocal ){

//   double fac;
//   int i1, i2, i3, index, start_index, end_index, index_buf, i2_loc;
//   int n2 = fft->n2, n3 = fft->n3, n1 = fft->n1, n_proc_tot, block_dim, nblock;

//   /* Allocate buffers to send and receive data */

//   MPI_Comm_size( fft->mpi_comm, &n_proc_tot );
    
//   // Now distinguish in which direction the FFT is performed
//   if( direct_to_reciprocal ){
//     /* among i3 dimension */
//     for( i1 = 0; i1 < fft->local_n1; i1++ ){
//       for( i2 = 0; i2 < n2; i2++ ){
// 	for( i3 = 0; i3 < n3; i3++ ){
	  
// 	  // Fill the missing part 
//     }
      
//     /* among i2 dimension */
//     for( i1 = 0; i1 < fft->local_n1; i1++ ){
//       for( i3 = 0; i3 < n3; i3++ ){
// 	for( i2 = 0; i2 < n2; i2++ ){
	    
// 	  // Fill the missing part 
// 	}
//       }
//     }

//     /*
//      * Reorder the different data blocks to be contigous in memory.
//      * The new distribution will allow to use the Alltoall function
//      *
//      */

//     // Perform an Alltoall communication 

//     /*  among i1 dimension */
//     for( i3 = 0; i3 < n3; i3++ ){
//       for( i2 = 0; i2 < /*block_dim*/; i2++ ){
// 	for( i1 = 0; i1 < n1; i1++ ){
	  
// 	  // Fill the missing part 
// 	}
//       }
//     }

//     // Perform an Alltoall communication 

//     /*
//      * Reoder the different data blocks to be consistent with the initial distribution.
//      *
//      */      

//   }
//   else{
    
//     /* Implement the reverse transform */

//   }
  
// }

