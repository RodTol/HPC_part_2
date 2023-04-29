#!/bin/bash
#SBATCH -A tra23_units
#SBATCH -p m100_usr_prod
#SBATCH --time 01:00:00
#SBATCH -N 1                  # nodes
#SBATCH --ntasks-per-node=25
#SBATCH --gres=gpu:1          # gpus per node out of 4
#SBATCH --mem=246000          # memory per node out of 246000MB
#SBATCH --ntasks-per-core=1
#SBATCH --job-name=tolloi_mat_mul
#SBATCH -o ./slurm_output/guppy_server_1_%j.out

module load autoload spectrum_mpi/10.3.1--binary

cd /m100/home/usertrain/a08trb39/HPC_part_2/mat_mul_even
n_proc=25

for N in 1000 5000 10000 15000 20000
do
	make run N=$N CORES=$n_proc 
done

echo "JOB FINISHED"


