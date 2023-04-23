/* Assignment:
 * Parallelize the code, using fftw-mpi
 * This amount to
 *   - distribute the data contained in diffusivity, conc, dconc, aux1, aux2 in
 * the way fftw-mpi expects
 *   - modify the fftw calls in fftw-mpi in p_fftw_wrapper
 * You will need to modify the files
 *   - diffusion.c
 *   - derivative.c
 *   - fftw_wrapper.c
 * In these files you will find some HINTs, make good use of them :)
 *
 * Created by G.P. Brandino, I. Girotto, R. Gebauer
 * Last revision: March 2016
 */

#include "utilities.h"
#include <complex.h>
#include <fftw3-mpi.h>
#include <fftw3.h>
#include <math.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
// HINT: include mpi and fftw3-mpi
//       http://www.fftw.org/doc/MPI-Files-and-Data-Types.html#MPI-Files-and-Data-Types

#define MAX(x, y) (((x) > (y)) ? (x) : (y))

int main(int argc, char **argv) {

  // Dimensions of the system
  double L1 = 10., L2 = 10., L3 = 20.;
  // Grid size
  int n1 = 48, n2 = 48, n3 = 96;
  // time step for time integration
  double dt = 2.e-3;
  // number of time steps
  int nstep = 101;
  // Radius of diffusion channel
  double rad_diff = 0.7;
  // Radius of starting concentration
  double rad_conc = 0.6;
  double start, end;

  double *diffusivity, *conc, *dconc, *aux1, *aux2;

  int i1, i2, i3, ipol, istep, index;

  double f1conc, f2conc, f3conc, f1diff, f2diff, f3diff, fac;
  double x1, x2, x3, rr, r2mean, r2mean_global;
  double ss, ss_global;

  int irank, n_proc_tot;

  fftw_mpi_handler fft_h;
  int local_size_grid, global_size_grid;

  /*
   * Initializzation of the MPI environment
   *
   */

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &irank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_proc_tot);
  /*
   * initialize the fftw system and local dimension
   * as the value returned from the parallel FFT grid initializzation
   *
   */
  init_fftw(&fft_h, n1, n2, n3, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

  if (irank == 0) {
    int tmp;
    printf("I am %d, this is my n1: %ld \n", irank, fft_h.local_n1);
    for (int i = 1; i < n_proc_tot; i++) {
      MPI_Recv(&tmp, 1, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("I am %d, this is my n1: %d \n", i, tmp);
    }
  } else {
    MPI_Send(&(fft_h.local_n1), 1, MPI_INT, 0, irank, MPI_COMM_WORLD);
  }
  /*
   * Allocate distribute memory arrays
   * HINT: the arrays need to be distributed, so you have to set the correct
   * sizes Use fftw_mpi_local_size_3d to calculate sizes
   *       http://www.fftw.org/doc/MPI-Data-Distribution-Functions.html
   *
   * More hints are reported into the fft_wrapper.c file where the fft_init_mpi
   * is defined
   *
   */

  /*Lui qua aggiunge una variabile local_grid_size per raccogliere tutto*/
  local_size_grid = (fft_h.local_n1) * n2 * n3;
  global_size_grid = n1 * n2 * n3;

  diffusivity = (double *)malloc(local_size_grid * sizeof(double));
  conc = (double *)malloc(local_size_grid * sizeof(double));
  dconc = (double *)malloc(local_size_grid * sizeof(double));
  aux1 = (double *)malloc(local_size_grid * sizeof(double));
  aux2 = (double *)malloc(local_size_grid * sizeof(double));

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
      /*Qua prova a far andare i1 da 0 a local_n1, ma sposta
        il fatto che metà va da uno e metà all'altro direttamente
        dentro x1!*/
      if (irank == 1 && i3 == 24 && i2 == 22) {
        printf("I am %d, these are my X1: ", irank);
      }
      
      for (i1 = 0; i1 < fft_h.local_n1; ++i1) {
        x1 = L1 * ((double)i1 + fft_h.local_n1_offset) / n1;
        if (irank == 1 && i3 == 24 && i2 == 22) {
          printf("%.3f ", x1);
          if (i1 == fft_h.local_n1-1) printf("\n");
        }

        f1diff = exp(-pow((x1 - 0.5 * L1) / rad_diff, 2));
        f1conc = exp(-pow((x1 - 0.5 * L1) / rad_conc, 2));

        index = index_f(i1, i2, i3, fft_h.local_n1, n2, n3);
        diffusivity[index] = MAX(f1diff * f2diff, f2diff * f3diff);
        conc[index] = f1conc * f2conc * f3conc;
        ss += conc[index];
      }
    }
  }
  /*
   * HINT: The parallel version of  the output routines is provided in the
   * mpi_output_routines folder
   *
   */
  plot_data_2d("diffusivity", n1, n2, n3, fft_h.local_n1, fft_h.local_n1_offset,
               2, diffusivity);

  /*Questa non funziona ma forse è giusto così ?*/               
  //plot_data_2d("diffusivity", n1, n2, n3, fft_h.local_n1, fft_h.local_n1_offset,
  //             3, diffusivity);

  /*Giusto che sia con la global size*/
  fac = L1 * L2 * L3 / (global_size_grid);

  /*
   * Now normalize the concentration
   *
   * HINT:the global ss must be computed and propagated to all processes
   *      ss = 1.0/(ss*fac);
   *
   */
  MPI_Allreduce(&ss, &ss_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  ss_global = 1.0 / (ss_global * fac);
  for (int count = 0; count < local_size_grid; ++count) {
    conc[count] *= ss_global;
  }
  plot_data_2d("concentration_init",n1,n2,n3,fft_h.local_n1,fft_h.local_n1_offset,
                2, conc);

  
  /*
   * Now everything is defined: system size, diffusivity inside the system, and
   * the starting concentration
   *
   * Start the dynamics
   *
   */
  // start = seconds();
  // for (istep = 1; istep <= nstep; ++istep)
  //   {
  //     for (i1=0; i1< n1*n2*n3; ++i1)
  // dconc[i1] = 0.0;
  //     for (ipol =1; ipol<=3; ++ipol )
  // {
  //   derivative(&fft_h, n1, n2, n3, L1, L2, L3, ipol, conc, aux1);
  //   for (i1=0; i1< n1*n2*n3; ++i1)
  //     {
  //             aux1[i1] *= diffusivity[i1];
  //     }
  //   derivative(&fft_h, n1, n2, n3, L1, L2, L3, ipol, aux1, aux2);
  //         // summing up contributions from the three spatial directions
  //         for (i1=0; i1< n1*n2*n3; ++i1)
  //     dconc[i1] += aux2[i1];
  // }
  //     for (i1=0; i1< n1*n2*n3; ++i1)
  // conc[i1] += dt*dconc[i1];

  //     if (istep%30 == 1)
  // {
  //         // Check the normalization of conc
  //         ss = 0.;
  //         r2mean = 0.;
  //         // HINT: the conc array is distributed, so only a part of it is on
  //         each processor for (i3 = 0; i3 < n3; ++i3)
  //     {
  //             x3=L3*((double)i3)/n3 - 0.5*L3;
  //             for (i2 = 0; i2 < n2; ++i2)
  //   {
  //                 x2=L2*((double)i2)/n2 - 0.5*L2;
  //                 for (i1 = 0; i1 < n1; ++i1)
  //       {
  // 	x1=L1*((double)i1)/n1 - 0.5*L1;
  // 	rr = pow( x1, 2)  + pow( x2, 2) + pow( x3, 2);
  // 	index = index_f(i1, i2, i3, n1, n2, n3);
  // 	ss += conc[index];
  // 	r2mean += conc[index]*rr;
  //       }
  //   }
  //     }

  //         /*
  //    * HINT: global values of ss and r2mean must be globally computed and
  //    distributed to all processes
  //    *
  //    */
  //         ss *= fac;
  //         r2mean *= fac;
  //         end = seconds();
  //         printf(" %d %17.15f %17.15f Elapsed time per iteration %f \n ",
  //         istep, r2mean, ss, (end-start)/istep);
  //         // HINT: Use parallel version of output routines
  //         plot_data_2d("concentration", n1, n2, n3, fft_h.local_n1,
  //         fft_h.local_n1_offset, 2, conc); plot_data_1d("1d_conc", n1, n2,
  //         n3, fft_h.local_n1, fft_h.local_n1_offset, 3, conc);
  // }

  //   }

  close_fftw(&fft_h);
  free(diffusivity);
  free(conc);
  free(dconc);
  free(aux1);
  free(aux2);

  /*
   * Finalize the MPI environment
   */
  MPI_Finalize();
  return 0;
}
