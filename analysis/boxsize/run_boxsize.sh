#!/bin/bash
AVGBOXSIZE=0
N_BOXSIZE=0
for ((k = k_0; k < (k_0 + N); k++)); do
    INPATH=${SIMPATH}_$k/
    if [ ! -f ${INPATH}ions_${PROCESS}.log ]; then
        echo "Path ${INPATH}ions_${PROCESS}.log does not exist"
        continue
    fi
    # read the line after the line containing "Box-x" in the log file
    read BOXSIZE Y Z <<< $(awk '/Box-X/ {getline; print; exit}' ${INPATH}ions_${PROCESS}.log)

    # if $BOXSIZE ends in e+00, remove it
    if [[ $BOXSIZE == *e+00 ]]; then
        BOXSIZE=${BOXSIZE%e+00}
    else
        echo "Box size does not end in e+00"
        continue
    fi
    echo $BOXSIZE
    AVGBOXSIZE=$(echo $AVGBOXSIZE + $BOXSIZE | bc)
    N_BOXSIZE=$((N_BOXSIZE + 1))
done
AVGBOXSIZE=$(echo $AVGBOXSIZE / $N_BOXSIZE | bc -l)
echo "Average box size: $AVGBOXSIZE"
echo "Number of box sizes: $N_BOXSIZE"

# Save the average box size to a file
echo $AVGBOXSIZE > ${HEADPATH}boxlen.dat