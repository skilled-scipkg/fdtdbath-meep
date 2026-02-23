#!/bin/bash

#SBATCH --job-name=outcav
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=taoeli
#SBATCH --ntasks-per-node=48

angle=0

for sigma in 0.002 0.02
do

DIR=incav_lb_lorentzian_sigma_$sigma\_angle_$angle
mkdir $DIR
cd $DIR
mpirun -np 48 python ../run_sc_fabry_perot_2d_spectrum.py -n 100 -c $angle -s $sigma -d 0.99 -k 0.01
cd ..

done
