# First exercise on FFTW

## The Assignment 
The objective is to implement a code that performs the evolution of a system using the fftw-mpi library. 
The library is quite straightforward for using. [here](http://www.fftw.org/fftw3_doc/2d-MPI-example.html) there's a very explanatory example.

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
NOTE: after a certain size, is necessary to decrease the time step, otherwise the algorithm won't converge