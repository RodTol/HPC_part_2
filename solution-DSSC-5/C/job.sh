#!/bin/bash
#SBATCH -A tra23_units
#SBATCH -p m100_usr_prod
#SBATCH --time 01:00:00       # format: HH:MM:SS
#SBATCH -N 1                  # nodes
#SBATCH --ntasks-per-node=6   # tasks out of 128
#SBATCH --gres=gpu:1          # gpus per node out of 4
#SBATCH --mem=246000          # memory per node out of 246000MB
#SBATCH --ntasks-per-core=1
#SBATCH --job-name=fftw_tolloi

module load autoload fftw/3.3.8--spectrum_mpi--10.3.1--binary
module load autoload gnuplot/5.2.6

cd /m100/home/usertrain/a08trb39/HPC_part_2/solution-DSSC-5/C
make clean
make flush

make

mpirun -np 3 diffusion.x

make plot
echo "JOB FINISHED"


