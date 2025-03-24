#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -e stderr-so
#SBATCH -o stdout-so
#SBATCH -p msismall
#SBATCH --mem=80g
#SBATCH --mail-type=FAIL  
#SBATCH --mail-user=tramp037@umn.edu 
module load conda
eval "$(conda shell.bash hook)"              
conda activate mdanalysis

# check if the number of arguments is correct
if [ "$#" -lt 9 ]; then
    echo "Usage: $0 <main simulation path> <simulation name> <output name> <start> <end> <skip> <ions...>"
    exit 1
fi

MPATH=${1}
STDNAME=${2}
OUTNAME=${3}
XBINS=${4}
YBINS=${5}
START=${6}
END=${7}
SKIP=${8}
# the remaining arguments should be held in variable IONS
IONS=("${@:9}")

echo "Running analysis"

python3 solvation_field.py -p "$MPATH" -c "$STDNAME".gro -t "$STDNAME".xtc -o "$OUTNAME" -b "$XBINS" "$YBINS" -a "${IONS[@]}" --start "$START" --end "$END" --skip "$SKIP"
