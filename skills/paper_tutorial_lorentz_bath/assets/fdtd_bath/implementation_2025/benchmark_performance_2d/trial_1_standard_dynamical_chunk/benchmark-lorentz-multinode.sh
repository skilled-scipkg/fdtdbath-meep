#!/bin/bash
#SBATCH --job-name=lorentz-multinode
#SBATCH --partition=standard
#SBATCH --time=0:05:00
#SBATCH --nodes=2                  
#SBATCH --ntasks-per-node=48  
#SBATCH --cpus-per-task=1
##SBATCH --mem-per-cpu=4G

srun --mpi=pmix python ../run_sc_fabry_perot_2d_spectrum.py

