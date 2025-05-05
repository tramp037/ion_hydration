#!/bin/bash

# *** MODIFY INPUTS HERE ***

# Replicates to run
REP=1
# Starting replicate
START=1

# Slurm settings
SUBMIT=sbatch
TIME=24:00:00
PART=msismall
GPU=no
ACC=sarupria

# *** END INPUTS ***

for i in $(seq $START $(($START + $REP - 1))); do
    cd run_${i}
    if [ $SUBMIT == "sbatch" ]; then
	if [ $GPU == "no" ]; then
	    sbatch --time=${TIME} --partition=${PART} -A ${ACC} --job-name=!COMPONENTS!_${i}_!EFIELD! run_md.sh !NPT!
	else
            sbatch --time=${TIME} --partition=${PART} --gres=gpu:${GPU}:1 -A ${ACC} --job-name=!COMPONENTS!_${i}_!EFIELD! run_md.sh !NPT!
    	fi
    else
        bash run_md.sh !NPT!
    fi
    cd ..
done
