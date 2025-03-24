#!/bin/bash

# *** MODIFY INPUTS HERE ***

# Simulations Folder
MAIN_DIR="$(git rev-parse --show-toplevel)"
echo "You are in ${MAIN_DIR} repo"
MPATH="${MAIN_DIR}/simulations/"

# MD/MC?
MODE=md

# Executable
EXEC="/home/sarupria/shared/software/load_scripts/load_gromacs-2022.sh"

# List of cations and anions to simulate
# CATION1=("Li" "Li" "Li" "Li" "Na" "Na" "Na" "K" "K" "Rb")
# CATION2=("Na" "K" "Rb" "Cs" "K" "Rb" "Cs" "Rb" "Cs" "Cs")
CATION=("Na")
ANIONS=($1)

# Water model
# Currently supported: TIP4P (TIP4P), TIP4P2005 (TIP4P/2005), SPCE (SPC/E)
WATERMODEL=SPCE
# Currently supported: Madrid2019 (Madrid-2019), Dang (Dang), JCSPCE (JC for SPC/E)
# JC is the Joung Cheatham set of parameters for 3 different water models
IONMODEL=JCSPCE

# Box legnth (nm)
BOX_LEN=2.6

# Number of solvent molecules, cations, and anions
N_SOL=555
# One ion pair per 55 water molecules is 36 ion pairs per 1980 water molecules
N_CAT=(1)
N_AN=(1)

# Pressure (bar)
PRESSURE=1.0
# Temperature (K)
TEMPERATURE=300
# Electric Field (V/nm)
EFIELD=0.0
# Simulation time (ns)
SIM_TIME=100

# NpT?
NPT=y

# Starting index for naming
START=1
# Number of indexes to run
N=16

# *** END INPUTS ***

source ${EXEC}

source make_${MODE}.sh

# echo "bash ion_gen.sh "$MPATH" "$BOX_LEN" "$CATION" "$ANION" "$WATERMODEL" "$PRESSURE" "$TEPMERATURE" "$EFIELD" "$N_SOL" "$N_CAT" "$N_AN" "$START" "$N" "$SIM_TIME" "$RUN" "$SLURMTIME""

# source ion_gen.sh "$MPATH" "$BOX_LEN" "$CATION" "$ANIONS" "$WATERMODEL" "$PRESSURE" "$TEPMERATURE" "$EFIELD" "$N_SOL" "$N_CAT" "$N_AN" "$START" "$N" "$SIM_TIME" "$RUN" "$SLURMTIME"
