#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -e stderr-rdf
#SBATCH -o stdout-rdf
#SBATCH -p msismall
#SBATCH --mem=80g
#SBATCH --mail-type=FAIL  
#SBATCH --mail-user=tramp037@umn.edu 
module load conda
eval "$(conda shell.bash hook)"              
conda activate mdanalysis

# check if the number of arguments is correct
if [ "$#" -ne 11 ]; then
    echo "Usage: $0 <main simulation path> <simulation name> <output name> <cation> <anion> <limit> <bins radius> <bins theta> <start> <end> <skip>"
    exit 1
fi

MPATH=${1}
STDNAME=${2}
OUTNAME=${3}
CATION=${4}
ANION=${5}
LIMIT=${6}
BINS_R=${7}
BINS_T=${8}
START=${9}
END=${10}
SKIP=${11}

echo "$CATION"
echo "$ANION"

# python3 radial_density.py
python3 radial_density.py -p "$MPATH" -c "$STDNAME".gro -t "$STDNAME".xtc -o "$OUTNAME"-"$CATION"-OW -a OW "$CATION" -l "${LIMIT}" -b "${BINS_R}" "${BINS_T}" --start "$START" --end "$END" --skip "$SKIP" 
python3 radial_density.py -p "$MPATH" -c "$STDNAME".gro -t "$STDNAME".xtc -o "$OUTNAME"-"$CATION"-HW -a HW "$CATION" -l "${LIMIT}" -b "${BINS_R}" "${BINS_T}" --start "$START" --end "$END" --skip "$SKIP"
python3 radial_density.py -p "$MPATH" -c "$STDNAME".gro -t "$STDNAME".xtc -o "$OUTNAME"-"$ANION"-OW -a OW "$ANION" -l "${LIMIT}" -b "${BINS_R}" "${BINS_T}" --start "$START" --end "$END" --skip "$SKIP" 
python3 radial_density.py -p "$MPATH" -c "$STDNAME".gro -t "$STDNAME".xtc -o "$OUTNAME"-"$ANION"-HW -a HW "$ANION" -l "${LIMIT}" -b "${BINS_R}" "${BINS_T}" --start "$START" --end "$END" --skip "$SKIP"