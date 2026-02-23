#!/bin/bash

#SBATCH --job-name=outcav
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=1
#SBATCH --partition=taoeli

# out cavity, Lorentz  
DIR=outcav_lorentz
mkdir $DIR
cd $DIR
python ../run_sc_fabry_perot_1d.py -r 1.0
cd ..

# out cavity, Lorentz-bath with uniform bath  
DIR=outcav_lb_uniform
mkdir $DIR
cd $DIR
python ../run_sc_fabry_perot_1d.py -r 1.0 -n 100 -l "uniform" -d 0.99 -k 0.01
cd ..

# out cavity, Lorentz-bath with lorentzian bath (default)  
DIR=outcav_lb_lorentzian
mkdir $DIR
cd $DIR
python ../run_sc_fabry_perot_1d.py -r 1.0 -n 100 -d 0.99 -k 0.01
cd ..
