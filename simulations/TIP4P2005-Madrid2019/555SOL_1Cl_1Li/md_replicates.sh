#!/bin/bash

# *** MODIFY INPUTS HERE ***

# NPT or efield (y/n)
NPT=n

# Simulation time (ns)
SIMTIME=100
# Pressure (bar)
PRESSURE=1.0
# Temperature (K)
TEMPERATURE=300
# Electric Field (V/nm)
EFIELD=1.0
# Replicates to run
REP=16
# Starting replicate
START=1

# Slurm settings
SUBMIT=sbatch
TIME=24:00:00
PART=msismall
GPU=no
ACC=siepmann

# *** END INPUTS ***

if [ $NPT == "y" ]; then
    PATH_NAME=run_${SIMTIME}ns_${PRESSURE}bar_${TEMPERATURE}K_
else
    PATH_NAME=run_${SIMTIME}ns_${PRESSURE}bar_${TEMPERATURE}K_${EFIELD}V_
fi

for i in $(seq $START $(($START + $REP - 1))); do
    cd ${PATH_NAME}${i}
    if [ $SUBMIT == "sbatch" ]; then
	if [ $GPU == "no" ]; then
	    sbatch --time=${TIME} --partition=${PART} -A ${ACC} --job-name=${PATH_NAME} run_md.sh $NPT
	else
            sbatch --time=${TIME} --partition=${PART} --gres=gpu:${GPU}:1 -A ${ACC} --job-name=${PATH_NAME} run_md.sh $NPT
    	fi
    else
        bash run_md.sh $NPT
    fi
    cd ..
done
