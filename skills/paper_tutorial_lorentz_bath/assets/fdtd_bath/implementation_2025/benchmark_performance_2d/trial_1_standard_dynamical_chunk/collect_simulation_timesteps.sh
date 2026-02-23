#!/bin/bash

for METHOD in "lb" "lorentz"
do

FILENAME=info_$METHOD.txt
rm $FILENAME
touch $FILENAME

for NCPU in 1 2 4 8 16 24 48 72 96 192 240
do

DIR=$METHOD\_ncpu\_$NCPU

if (($NCPU > 24)); then
  DIR=$METHOD\_multinode_ncpu_$NCPU
fi

timestep=$(grep -i "on time step" $DIR/slurm*.out  | tail -n 1 | awk '{print $6}')

if [ -z "$timestep" ]; then
timestep=nan
fi

echo $DIR
echo "$NCPU $timestep"
echo "$NCPU $timestep" >> $FILENAME

done

done
