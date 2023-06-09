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


## How to use the code
The user just need to launch:
```
make run CORES=<number of process> flags=<info or debug for extra information>
```
NOTE: after a certain size, is necessary to decrease the time step, otherwise the algorithm won't converge

## Some results
Here are some result from the code
|    | grid_size    |   n_proc_tot |     time |
|---:|:-------------|-------------:|---------:|
|  0 | 512x512x1024 |           32 | 1051.17  |
|  1 | 512x512x1024 |           64 |  659.826 |
|  2 | 512x512x1024 |          128 |  454.705 |
|  3 | 512x512x1024 |          256 |  276.522 | 

<figure>
  <img
  src="../images/fft1.png"
  alt="1 Node"
  width="500" 
  height="400" >
</figure>