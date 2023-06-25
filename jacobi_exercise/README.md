## Software Stack
For the openACC solution I used:
```
module load hpc-sdk
```
and for the MPI one:
```
module load spectrum_mpi/
```

## Some result
Code execution is severly slowed by the writing on file of the initial and final matrix. The function is not parallelized and for large size of the problem, the printing can take several minutes.
