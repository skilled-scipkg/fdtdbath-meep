#!/bin/bash

for TOTAL in 1 2 4 8 16 24; do

# decide on memory request:
if (( TOTAL > 16 )); then
  # more than 16 CPUs -> 4 GB per CPU
  MEM_ARG="--mem-per-cpu=4G"
else
  # 16 or fewer CPUs -> at least 64 GB total
  MEM_ARG="--mem=64G"
fi

DIR=lorentz_ncpu_$TOTAL
echo "working under $DIR"
mkdir $DIR
cp benchmark-lorentz-multinode.sh $DIR
cd $DIR

  sbatch --job-name=bench-${TOTAL}cpus \
         --nodes=1 \
         --ntasks-per-node=$TOTAL \
         ${MEM_ARG} \
         benchmark-lorentz-multinode.sh

cd ..

done
