/* Assignment: 
 * Parallelize the code, using fftw-mpi
 * This amount to 
 *   - distribute the data contained in diffusivity, conc, dconc, aux1, aux2 in the way fftw-mpi expects
 *   - modify the fftw calls in fftw-mpi in p_fftw_wrapper
 * You will need to modify the files
 *   - diffusion.c   
 *   - derivative.c
 *   - fftw_wrapper.c
 *   - plot_data.c
 * In these files you will find some HINTs, make good use of them :)
 *
 * Created by G.P. Brandino, I. Girotto, R. Gebauer
 * Last revision: March 2016
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "headers/utilities.h"
// HINT: include mpi and fftw3-mpi 
//       http://www.fftw.org/doc/MPI-Files-and-Data-Types.html#MPI-Files-and-Data-Types

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

/**
 * @brief Color Red for the output of C program
 */
void printf_red () {
  printf("\033[1;31m");
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


int main( int argc, char* argv[] ){

  // Dimensions of the system
  double L1 = 10., L2 = 10., L3 = 20.;
  // Grid size  
  int n1 = 48, n2 = 48, n3 = 96;
  // time step for time integration
  double dt = 2.e-3; 
  // number of time steps
  int nstep = 300; 
  // Radius of diffusion channel
  double rad_diff = 0.7;
  // Radius of starting concentration
  double rad_conc = 0.6;
  double start, end;

  double *diffusivity, *conc, *dconc, *aux1, *aux2;
  
  int i1, i2, i3, ipol, istep, index;

  double f1conc, f2conc, f3conc, f1diff, f2diff, f3diff, fac;
  double x1, x2 , x3, rr;
  double ss, r2mean, global_ss, global_r2mean;

  fftw_dist_handler fft_h;
  int irank, n_proc_tot;
  int n1_local, n1_local_offset;
  int n2_local, n2_local_offset;
  int local_size_grid, global_size_grid;

  /* 
    * Initializzation of the MPI environment 
    *
    */
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &irank);
  MPI_Comm_size( MPI_COMM_WORLD, &n_proc_tot);

  /*
    * Allocate distribute memory arrays
    * HINT: the arrays need to be distributed, so you have to set the correct sizes 
    *       Use fftw_mpi_local_size_3d to calculate sizes
    *       http://www.fftw.org/doc/MPI-Data-Distribution-Functions.html
    *
    * More hints are reported into the fft_wrapper.c file where the fft_init_mpi is defined
    *
    */

  /*
    * initialize the fftw system and local dimension
    * as the value returned from the parallel FFT grid initializzation 
    *
    */
  init_fftw( &fft_h, n1, n2, n3, MPI_COMM_WORLD );

  local_size_grid = (fft_h.local_n1) * n2 * n3;
  global_size_grid = n1 * n2 * n3;
  n1_local = fft_h.local_n1;
  n1_local_offset = fft_h.local_n1_offset;

  n1_local = fft_h.local_n1;
  n1_local_offset = fft_h.local_n1_offset;  

  n2_local = fft_h.local_n2;
  n2_local_offset = fft_h.local_n2_offset;  
  
  diffusivity = ( double* ) malloc( local_size_grid * sizeof(double) );
  conc = ( double* ) malloc( local_size_grid * sizeof(double) );
  dconc = ( double* ) malloc( local_size_grid * sizeof(double) );
  aux1 = ( double* ) malloc( local_size_grid * sizeof(double) );
  aux2 = ( double* ) malloc( local_size_grid * sizeof(double) );

  /*
   * Define the diffusivity inside the system and
   * the starting concentration
   *
   * ss is to integrate (and normalize) the concentration
   *
   */

  ss = 0.0;

  for (i3 = 0; i3 < n3; ++i3) {
    x3 = L3 * ((double)i3) / n3;
    f3diff = exp(-pow((x3 - 0.5 * L3) / rad_diff, 2));
    f3conc = exp(-pow((x3 - 0.5 * L3) / rad_conc, 2));
    for (i2 = 0; i2 < n2; ++i2) {
      x2 = L2 * ((double)i2) / n2;
      f2diff = exp(-pow((x2 - 0.5 * L2) / rad_diff, 2));
      f2conc = exp(-pow((x2 - 0.5 * L2) / rad_conc, 2));

      if (irank == 1 && i3 == 24 && i2 == 22) {
        printf("I am %d, these are my X1: ", irank);
      }    
      /*Qua prova a far andare i1 da 0 a local_n1, ma sposta
        il fatto che metà va da uno e metà all'altro direttamente
        dentro x1!*/
      for (i1 = 0; i1 < n1_local; ++i1) {
        x1 = L1 * ((double)i1 + n1_local_offset) / n1;

        if (irank == 1 && i3 == 24 && i2 == 22) {
          printf("%.3f ", x1);
          if (i1 == n1_local-1) printf("\n");
        }

        f1diff = exp(-pow((x1 - 0.5 * L1) / rad_diff, 2));
        f1conc = exp(-pow((x1 - 0.5 * L1) / rad_conc, 2));

        index = index_f(i1, i2, i3, n1_local, n2, n3);
        diffusivity[index] = MAX(f1diff * f2diff, f2diff * f3diff);
        conc[index] = f1conc * f2conc * f3conc;
        ss += conc[index];
      }
    }
  }

  plot_data_2d("data/diffusivity1", n1, n2, n3, fft_h.local_n1, fft_h.local_n1_offset,
               1, diffusivity);

  plot_data_2d("data/diffusivity2", n1, n2, n3, fft_h.local_n1, fft_h.local_n1_offset,
               2, diffusivity);

  plot_data_2d("data/diffusivity3", n1, n2, n3, fft_h.local_n1, fft_h.local_n1_offset,
               3, diffusivity);

  /*Debuggin for r2mean and ss*/
  if (irank == 0) {
    double tmp_ss;
    printf("I am %d, this is my ss: %17.15f\n", irank, ss);
    for (int i = 1; i < n_proc_tot; i++) {
      MPI_Recv(&tmp_ss, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("I am %d, this is my ss: %17.15f\n", i, ss);
    }
  } else {
    MPI_Send(&ss, 1, MPI_DOUBLE, 0, irank, MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  /*
  * Now normalize the concentration and print
  *
  */
  fac = L1 * L2 * L3 / (global_size_grid);
  MPI_Allreduce(&ss, &global_ss, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  global_ss = 1.0 / (global_ss * fac);

  for (int count = 0; count < local_size_grid; ++count) {
    conc[count] *= global_ss;
  }

  plot_data_2d("data/concentration_init",n1,n2,n3,n1_local,n1_local_offset,
                2, conc);

  /*
    * Now everything is defined: system size, diffusivity inside the system, and
    * the starting concentration
    *
    * Start the dynamics
    *
    */

  if (irank == 0) {
    printf_red();
    printf("Dynamic is starting\n Local_size_grid is %d\n fac is: %17.15f \n",
   local_size_grid, fac);
    printf_reset();
  }
  
  /*Questa variabile indica ogni quanti frame voglio fare una stampa*/
  int interval = 10;
  start = seconds();

  /*Inizializzo dconc (si può mettere fuori perche
  è ridondante)*/
  for (i1=0; i1< local_size_grid; ++i1) dconc[i1] = 0.0;

  for (istep = 1; istep <= nstep; ++istep) {
    /*Devo fare la trasf sulle 3 direzioni e sommare i contributi*/
    
    for (ipol =1; ipol<=3; ++ipol ) {

        derivative(&fft_h, n1, n2, n3, L1, L2, L3, ipol, conc, aux1);

        for (i1=0; i1< local_size_grid; ++i1) {
          aux1[i1] *= diffusivity[i1];
        }

        derivative(&fft_h, n1, n2, n3, L1, L2, L3, ipol, aux1, aux2);

        for (i1=0; i1< local_size_grid; ++i1) {
          dconc[i1] += aux2[i1];
        } 
    }

    /*Aggiorno l'array della concentrazione*/
    for (i1=0; i1< local_size_grid; ++i1) {
      conc[i1] += dt*dconc[i1];
      dconc[i1] = 0.0;
    } 
    
    /*Voglio un frame ogni 10*/
    if (istep%interval == 1) {
      // Check the normalization of conc
      ss = 0.;
      r2mean = 0.;

      // HINT: the conc array is distributed, so only a part of it is on
      // each processor

      for (i3 = 0; i3 < n3; ++i3) {
        x3=L3*((double)i3)/n3 - 0.5*L3;
        for (i2 = 0; i2 < n2; ++i2) {
          x2=L2*((double)i2)/n2 - 0.5*L2;
          for (i1 = 0; i1 < n1_local; ++i1) {
            x1 = L1 * ( (double) (i1 + n1_local_offset) ) / n1 - 0.5 * L1;
            rr = pow( x1, 2)  + pow( x2, 2) + pow( x3, 2);
            index = index_f(i1, i2, i3, n1_local, n2, n3);
            ss += conc[index];
            r2mean += conc[index]*rr;
          }
        }
      }

      /*Debuggin for r2mean and ss*/
      /*if (irank == 0) {
        double tmp_ss;
        double tmp_r2mean;
        printf("I am %d, this is my ss: %17.15f and r2mean: %17.15f\n", irank, ss, r2mean);
        for (int i = 1; i < n_proc_tot; i++) {
          MPI_Recv(&tmp_ss, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Recv(&tmp_r2mean, 1, MPI_DOUBLE, i, i+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          printf("I am %d, this is my ss: %17.15f and r2mean: %17.15f\n", i, tmp_ss, tmp_r2mean);
        }
      } else {
        MPI_Send(&ss, 1, MPI_DOUBLE, 0, irank, MPI_COMM_WORLD);
        MPI_Send(&ss, 1, MPI_DOUBLE, 0, irank+1, MPI_COMM_WORLD);
      }
      MPI_Barrier(MPI_COMM_WORLD);*/

      /*
      * HINT: global values of ss and r2mean must be globally computed and
      distributed to all processes
      *
      */
      MPI_Allreduce( &ss, &global_ss, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      MPI_Allreduce( &r2mean, &global_r2mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
      
      global_ss = global_ss * fac;
      global_r2mean = global_r2mean * fac;

      end = seconds();

      if (irank == 0) {
        printf_green();
        printf("Time: %d , r2mean: %17.15f, ss: %17.15f, Elapsed time per iteration %f, rank: %d \n",
          istep, global_r2mean, global_ss, (end-start)/istep, irank);
        printf_reset();
      }

      // HINT: Use parallel version of output routines
      char title[80];
      sprintf(title, "data/concentration_%d", 1 + (istep - 1) / interval);
      plot_data_2d(title, n1, n2, n3, n1_local,
        n1_local_offset, 2, conc);
      //plot_data_1d("1d_conc", n1, n2, n3, fft_h.local_n1, fft_h.local_n1_offset,
        //3, conc);
    }
  }

  close_fftw(&fft_h);
  free(diffusivity);
  free(conc);
  free(dconc);
  free(aux1);
  free(aux2);

  MPI_Finalize();
  return 0;
} 
