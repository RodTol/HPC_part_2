# The assignement
The objective of this assignment is to implement the matrix multiplication in 3 different ways:
1. Naive multiplication: just some for loops on the data
2. Use the of openblas dgemm
3. Use of the cublas dgemm
In all these cases, the matrices are distributed among one direction. There's no need for a specific size of the matrices because the code has been written with the capa

# Software Stack
For this solution I used the following modules on Marconi100:
1. profile/base
2. autoload   
3. gnu/8.4.0 
4. openblas/0.3.9--gnu--8.4.0   
5. cuda/11.0   
6. spectrum_mpi/10.4.0--binary  

Secondo me la max matrices è 44000*44000: questo perchè sono 3 matrici (A, B e C)

# How to run the code
I created a [launcher.sh](launcher.sh) file but it's useful only in some occasion. For a more general command use:
```
make run N=<size> CORES=<number of total process> n_socket=<# of process for each socket> n_node=<#of process for each node> compilation=<dgemm or gpu> flags=<debug>
```
The are some default values: n_socket = 16, n_node = 32. If the compilation variable is not specified, the naive versione will be compiled.