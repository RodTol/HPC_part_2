
#!/bin/bash
#SBATCH -A tra23_units_bis
#SBATCH -p m100_usr_prod
#SBATCH --time 03:00:00
#SBATCH -N 8                   # nodes
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

cd /m100/home/usertrain/a08tra76/HPC_part_2/jacobi_exercise/mpi_code

n_proc=32
for n_nodes in 1 2 4 8
do
	((n_proc*=$n_nodes))
	make run dim=64000 time=4000 n_socket=16 n_node=32 CORES=$n_proc
done

cd /m100/home/usertrain/a08tra76/HPC_part_2/jacobi_exercise/openACC_code

n_proc=4
for n_nodes in 1 2 4 8
do
	((n_proc*=$n_nodes))
	make run dim=64000 time=4000 n_socket=2 n_node=4 CORES=$n_proc
done

echo "JOB FINISHED"


