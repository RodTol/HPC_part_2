#!/bin/bash
#SBATCH -A tra23_units_bis
#SBATCH -p m100_usr_prod
#SBATCH --time 03:00:00
#SBATCH -N 8                   # nodes
#SBATCH --ntasks-per-node=32
#SBATCH --gres=gpu:4          # gpus per node out of 4
#SBATCH --mem=246000          # memory per node out of 246000MB
#SBATCH --ntasks-per-core=1
#SBATCH --job-name=tolloi_fftw2
#SBATCH -o ./output/run.out
#SBATCH -e ./output/err.out

module purge
module load autoload fftw/

export OMP_NUM_THREADS=1

cd /m100/home/usertrain/a08tra76/HPC_part_2/fftw_es2/C

for n_nodes in 1 2 4 8
do
	n_proc=32
	((n_proc*=$n_nodes))
	make run CORES=$n_proc 
done

echo "JOB FINISHED"


