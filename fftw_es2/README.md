# Second exercise on fftw

## The Assignment
The objective is to implement a code that performs the evolution of a system using the fft, but without the use of the fftw-mpi library. In order to so, I have used the advanced interface of fftw and the MPI_Alltoall routine.

## The fft
The code works by splitting the data on the `n1` direction among all the process, and perform the fft on the `n2-n3` plan. After, it will reorder the data inside a buffer using the MPI_Alltoallw routine (see next paragraph) in such a way that the orginal matrix is now distributed among the process along the `n2` direction. It's now possible to perform the fft on the `n1` direction and then go back to the original data order. In this way there's no need to change the plot routines.  
To perform the actual fft I used the `fftw_plan_many_dft`, which allows to use only one call to perform the fft on a 2D or 1D plan. For both the 2D plan and the 1D plan, the data are considered contiguos. See [here](http://www.fftw.org/fftw3_doc/Advanced-Complex-DFTs.html) for information and examples on how to use it.

## Some information about the data transfer
I tested the usage of MPI_Alltoall on the file [all_to_all_test.c](all_to_all_test.c) and found out that it's easier to use the MPI_Alltoallw routine, instead of the standard one. This choice requires to insert all displacements in bytes. This is crucial since the sender and receiver use different datatypes and can lead to errors.  
I introduced a new MPI_Datatype, called column_block, which represents the matrix block that each process send and receive during the alltoall routine. The block will contain `n1_loc*n2_loc*n3` elements and, since we work in a row-major access, the next line is after `n2_loc*n3` elements.
The advantages of this approach is that both before and after the reordering, the data are contiguos through a row-major access.  

## Software stack
I used the following M100 modules:
1. profile/base
2. autoload
3. spectrum_mpi/10.3.1--binary
4. fftw/3.3.8--spectrum_mpi--10.3.1--binary

# How to use the code
The user just need to launch:
```
make run CORES=<number of process> flags=<info or debug for extra information>
```