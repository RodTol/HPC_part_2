#!/bin/bash
#SBATCH -A tra23_units
#SBATCH -p m100_usr_prod
#SBATCH --time 03:00:00
#SBATCH -N 8                  # nodes
#SBATCH --ntasks-per-node=32
#SBATCH --gres=gpu:4          # gpus per node out of 4
#SBATCH --mem=246000          # memory per node out of 246000MB
#SBATCH --ntasks-per-core=1
#SBATCH --job-name=tolloi_mat_mul
#SBATCH -o ./output/run.out
#SBATCH -e ./output/err.out

module load autoload spectrum_mpi
module load autoload openblas
module load cuda

export OMP_NUM_THREADS=1

cd /m100/home/usertrain/a08trb39/HPC_part_2/matrix_mult

n_proc=32
((n_proc*=8))

for N in 1000 5000 10000 15000 20000 25000 30000
#for N in 30000
do
	make run N=$N CORES=$n_proc
	make run N=$N CORES=$n_proc compilation=dgemm
	make run N=$N CORES=$n_proc compilation=gpu
done

echo "JOB FINISHED"


