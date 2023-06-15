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