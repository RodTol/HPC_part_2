#!/bin/bash -x

C_FLAGS="-O3"
mpicc ${C_FLAGS} -c benchmark_LBE3D_IO.c -I/cineca/prod/opt/libraries/hdf5/1.12.0--spectrum_mpi--10.3.1/pgi--19.10--binary/include
mpicc ${C_FLAGS} -o benchmark_LBE3D_IO.x benchmark_LBE3D_IO.o -L/cineca/prod/opt/libraries/hdf5/1.12.0--spectrum_mpi--10.3.1/pgi--19.10--binary/lib -lhdf5 -lm

