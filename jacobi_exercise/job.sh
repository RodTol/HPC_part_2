#!/bin/bash
#SBATCH -A tra23_units_bis
#SBATCH -p m100_usr_prod
#SBATCH --time 03:00:00
#SBATCH -N 8                   # nodes
#SBATCH --ntasks-per-node=32
#SBATCH --gres=gpu:4          # gpus per node out of 4
#SBATCH --mem=246000          # memory per node out of 246000MB
#SBATCH --ntasks-per-core=1
#SBATCH --job-name=tolloi_JACOBI
#SBATCH -o ./output/run.out
#SBATCH -e ./output/err.out

module load spectrum_mpi

export OMP_NUM_THREADS=1

cd /m100/home/usertrain/a08tra76/HPC_part_2/jacobi_exercise/mpi_code


for n_nodes in 1 2 4 8
do
	n_proc=32
	((n_proc*=$n_nodes))
	make run dim=28000 time=1000 CORES=$n_proc flags=noprint
done


module purge
module load hpc-sdk
#We have open-mpi 3.1.5

cd /m100/home/usertrain/a08tra76/HPC_part_2/jacobi_exercise/openACC_code


for n_nodes in 1 2 4 8
do
	n_proc=4
	((n_proc*=$n_nodes))
	make run dim=28000 time=1000  CORES=$n_proc flags=noprint
done

echo "JOB FINISHED"


