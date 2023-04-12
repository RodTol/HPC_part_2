#!/bin/bash

#SBATCH --job-name="fftw_bench"

#SBATCH --partition=regular1

#SBATCH --nodes=8

#SBATCH --ntasks=320

#SBATCH --ntasks-per-node=40

#SBATCH --time=00:05:00


module load gnu8
module load openmpi3
module load fftw

mpirun -np $SLURM_NTASKS /home/mcastell/FFTW/diffusion.x >> /home/mcastell/FFTW/benchmark/bench.dat




