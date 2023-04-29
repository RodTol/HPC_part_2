#!/bin/bash
#SBATCH -A tra23_units
#SBATCH -p m100_usr_prod
#SBATCH --time 00:10:00     # format: HH:MM:SS
#SBATCH -N 1                # 1 node
#SBATCH --ntasks-per-node=8 # 8 tasks out of 128
#SBATHC --exclusive
#SBATCH --gres=gpu:1        # 1 gpus per node out of 4
#SBATCH --mem=8000          # memory per node out of 246000MB
#SBATCH --job-name=matrix_transpose
#SBATCH --constraint=gpureport

srun ./a.out         #in all the other cases
