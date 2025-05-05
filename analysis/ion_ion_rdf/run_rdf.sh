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
if [ "$#" -lt 10 ]; then
    echo "Usage: $0 <main simulation path> <simulation name> <output name> <limit> <bins> <blocks> <start> <end> <skip> <ions...>"
    exit 1
fi

MPATH=${1}
STDNAME=${2}
OUTNAME=${3}
LIMIT=${4}
BINS=${5}
BLOCKS=${6}
START=${7}
END=${8}
SKIP=${9}
# the remaining arguments should be held in variable IONS
IONS=("${@:10}")

#python3 rdf.py
for i in "${IONS[@]}"; do
    for j in "${IONS[@]}"; do
        if [ "$i" != "$j" ]; then
            echo "Calculating RDF for $i-$j"
            python3 ../scripts/structure.py rdf -p "$MPATH" -c "$STDNAME".gro -t "$STDNAME".xtc -o "$OUTNAME"-"$i"-"$j" -a "$j" "$i" -l "${LIMIT}" -b "${BINS}" -k "${BLOCKS}" -x 0 0 -y 0 0 -z 0 0 --start "$START" --end "$END" --skip "$SKIP"
        fi
    done
done
