#!/bin/bash

#SBATCH --job-name=incav_lb
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=1
#SBATCH --partition=taoeli

for sigma in 0.002 0.020
do

# in cavity, Lorentz  
DIR=incav_lb_lorentzian_sigma_$sigma
mkdir $DIR
cd $DIR
python ../run_sc_fabry_perot_1d.py -n 100 -s $sigma -d 0.99 -k 0.01
cd ..

done
