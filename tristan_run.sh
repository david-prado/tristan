#!/bin/bash

for steps in 1000 10000 100000 1000000
do
echo "STEPS=$steps" > steps.sh
[ -z "$JOBID" ] && JOBID=$(qsub -N tristan_${steps}s_110000p tristan.pbs) || JOBID=$(qsub -W depend=afterok:$JOBID -N tristan_${steps}s_110000p tristan.pbs)
done
