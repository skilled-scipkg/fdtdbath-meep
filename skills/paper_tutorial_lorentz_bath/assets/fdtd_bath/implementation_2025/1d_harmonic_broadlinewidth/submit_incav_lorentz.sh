#!/bin/bash

#SBATCH --job-name=incav_lorentz
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=1
#SBATCH --partition=taoeli

for sigma in 0.002 0.020
do

# in cavity, Lorentz  
DIR=incav_lorentz_sigma_$sigma
mkdir $DIR
cd $DIR
python ../run_sc_fabry_perot_1d.py -s $sigma
cd ..

done
