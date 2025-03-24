#!/bin/bash

# *** MODIFY INPUTS HERE ***

# Simulations Folder
MAIN_DIR="$(git rev-parse --show-toplevel)"
echo "You are in ${MAIN_DIR} repo"
MPATH="${MAIN_DIR}/simulations/"

# List of cations and anion used
CATION=($1)
ANIONS=("Cl")
N_CAT=(1)
N_AN=(1)

# Water model
# WATERMODEL=TIP4P2005
WATERMODEL=SPCE
# IONMODEL=Madrid2019
IONMODEL=JCSPCE
N_SOL=555

# Starting index for naming
k_0=1
# Number of indexes to run
N=16
# Simulation time (ns)
SIMTIME=100

TEMPERATURE=300
PRESSURE=1.0

# step in the simulation process (em, eq, nvt, npt, md)
PROCESS=md

# *** END OF INPUTS ***


#Run analysis for each cation combination based on above parameters
NION=0
for i in "${N_CAT[@]}"; do
    NION=$((NION + i))
done
for i in "${N_AN[@]}"; do
    NION=$((NION + i))
done

WATERSTRING=${N_SOL}SOL
if [ "$NION" -eq 0 ]; then
    HEADPATH=${MPATH}${WATERMODEL}/"$WATERSTRING"/
    SIMPATH=${HEADPATH}run_"$SIMTIME"ns_"$PRESSURE"bar_"$TEMPERATURE"K
    WATER=y
else
    ANSTRING=""
    for i in "${!N_AN[@]}"; do
        if [ "${N_AN[$i]}" -ne 0 ]; then
            ANSTRING="${ANSTRING}${N_AN[$i]}${ANIONS[$i]}"
        fi
    done
    CATSTRING=""
    for i in "${!N_CAT[@]}"; do
        if [ "${N_CAT[$i]}" -ne 0 ]; then
            CATSTRING="${CATSTRING}${N_CAT[$i]}${CATION[$i]}"
        fi
    done
    HEADPATH=${MPATH}${WATERMODEL}-${IONMODEL}/"$WATERSTRING"_"$ANSTRING"_"$CATSTRING"/
    SIMPATH=${HEADPATH}run_"$SIMTIME"ns_"$PRESSURE"bar_"$TEMPERATURE"K
    WATER=n
fi
STDNAME=ions_${PROCESS}

# run job
source run_boxsize.sh
