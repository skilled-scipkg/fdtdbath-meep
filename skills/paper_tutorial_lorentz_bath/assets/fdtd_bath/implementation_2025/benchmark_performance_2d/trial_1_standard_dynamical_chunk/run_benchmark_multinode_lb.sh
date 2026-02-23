#!/bin/bash

for NODES in 2 3 4 6 8 10; do

MEM_ARG="--mem-per-cpu=4G"

TOTAL=$(( NODES * 24 ))
DIR=lb_multinode_ncpu_$TOTAL
echo "working under $DIR"
mkdir $DIR
cp benchmark-lb-multinode.sh $DIR
cd $DIR

  sbatch --job-name=bench-${NODES}nodes \
         --nodes=${NODES} \
         --ntasks-per-node=24 \
         ${MEM_ARG} \
         benchmark-lb-multinode.sh

cd ..

done
