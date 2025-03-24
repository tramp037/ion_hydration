#!/bin/bash

# *** MODIFY INPUTS HERE ***

# Simulations Folder
MAIN_DIR="$(git rev-parse --show-toplevel)"
echo "You are in ${MAIN_DIR} repo"
MPATH="${MAIN_DIR}/simulations/"

# Starting frame, final frame, and frames to skip
START=1
END=-1
SKIP=1

# Limit of RDF (in Angstroms)
LIMIT=10
# Number of bins
BINS_R=1000
BINS_T=100

# List of cations and anion used
CATION=($1)
ANIONS=("Cl")
N_CAT=(1)
N_AN=(1)

# Water model
WATERMODEL=TIP4P2005
IONMODEL=Madrid2019
N_SOL=555


# Starting index for naming
k_0=1
# Number of indexes to run
N=16
# Simulation time (ns)
SIMTIME=100

TEMPERATURE=300
PRESSURE=1.0
EFIELD=0.0

NPT=n

# step in the simulation process (em, eq, nvt, npt, md)
PROCESS=md

# output filename
OUTNAME=2d_radialdensity

# submit to slurm
SLURM=y
# slurm time limit
SLURMTIME=02:00:00

# *** END OF INPUTS ***


#Run analysis for each cation combination based on above parameters
NION=0
for i in "${N_CAT[@]}"; do
    NION=$((NION + i))
done
for i in "${N_AN[@]}"; do
    NION=$((NION + i))
done
for ((k = k_0; k < (k_0 + N); k++)); do
    WATERSTRING=${N_SOL}SOL
    if [ "$NION" -eq 0 ]; then
        if [ "$NPT" == "y" ]; then
            SIMPATH=${MPATH}${WATERMODEL}/"$WATERSTRING"/run_"$SIMTIME"ns_"$PRESSURE"bar_"$TEMPERATURE"K_"$k"/
        else
            SIMPATH=${MPATH}${WATERMODEL}/"$WATERSTRING"/run_"$SIMTIME"ns_"$PRESSURE"bar_"$TEMPERATURE"K_"$EFIELD"V_"$k"/
        fi
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
        if [ "$NPT" == "y" ]; then
            SIMPATH=${MPATH}${WATERMODEL}-${IONMODEL}/"$WATERSTRING"_"$ANSTRING"_"$CATSTRING"/run_"$SIMTIME"ns_"$PRESSURE"bar_"$TEMPERATURE"K_"$k"/
        else
            SIMPATH=${MPATH}${WATERMODEL}-${IONMODEL}/"$WATERSTRING"_"$ANSTRING"_"$CATSTRING"/run_"$SIMTIME"ns_"$PRESSURE"bar_"$TEMPERATURE"K_"$EFIELD"V_"$k"/
        fi
        WATER=n
    fi
    STDNAME=ions_${PROCESS}

    # define job name
    if [ "$WATER" == y ]; then
        JOB_NAME="2D_${WATERMODEL}_${k}"
    else
        JOB_NAME="2D_${WATERMODEL}-${IONMODEL}_${k}"
    fi

    # echo the first value f cation
    echo ${CATION[0]}
    # run job
    echo "Submitting "${JOB_NAME}
    if [ "$SLURM" == "y" ]; then
        sbatch --job-name="$JOB_NAME" -t "$SLURMTIME" run_2d.sh "${SIMPATH}" "${STDNAME}" "${OUTNAME}" "${CATION[0]}" "${ANIONS[0]}" "${LIMIT}" "${BINS_R}" "${BINS_T}" "${START}" "${END}" "${SKIP}"
    else
        bash run_2d.sh "${SIMPATH}" "${STDNAME}" "${OUTNAME}" "${CATION[0]}" "${ANIONS[0]}" "${LIMIT}" "${BINS_R}" "${BINS_T}" "${START}" "${END}" "${SKIP}"
    fi
done
