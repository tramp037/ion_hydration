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
if [ "$#" -lt 7 ]; then
    echo "Usage: $0 <main simulation path> <simulation name> <output name> <start> <end> <skip> <ions...>"
    exit 1
fi

MPATH=${1}
STDNAME=${2}
OUTNAME=${3}
START=${4}
END=${5}
SKIP=${6}
# the remaining arguments should be held in variable IONS
IONS=("${@:7}")

echo "Running analysis"

for i in "${IONS[@]}"; do
    for j in "${IONS[@]}"; do
        if [ "$i" != "$j" ]; then
            echo "Calculating Distance for $i-$j"
            python3 ion_ion_distance.py -p "$MPATH" -c "$STDNAME".gro -t "$STDNAME".xtc -o "$i-$j-distance" -a $i $j --start "$START" --end "$END" --skip "$SKIP"
        fi
    done
done