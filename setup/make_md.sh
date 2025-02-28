#!/bin/bash

NION=0
for i in "${N_CAT[@]}"; do
    NION=$((NION + i))
done
for i in "${N_AN[@]}"; do
    NION=$((NION + i))
done

for ((k = START; k < (START + N); k++)); do
    # name the simulation path
    WATERSTRING=${N_SOL}SOL
    if [ "$NION" -eq 0 ]; then
        if [ "$NPT" == "y" ]; then
            SIMPATH=${MPATH}${WATERMODEL}/"$WATERSTRING"/run_"$SIM_TIME"ns_"$PRESSURE"bar_"$TEMPERATURE"K_"$k"/
        else
            SIMPATH=${MPATH}${WATERMODEL}/"$WATERSTRING"/run_"$SIM_TIME"ns_"$PRESSURE"bar_"$TEMPERATURE"K_"$EFIELD"V_"$k"/
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
            SIMPATH=${MPATH}${WATERMODEL}-${IONMODEL}/"$WATERSTRING"_"$ANSTRING"_"$CATSTRING"/run_"$SIM_TIME"ns_"$PRESSURE"bar_"$TEMPERATURE"K_"$k"/
        else
            SIMPATH=${MPATH}${WATERMODEL}-${IONMODEL}/"$WATERSTRING"_"$ANSTRING"_"$CATSTRING"/run_"$SIM_TIME"ns_"$PRESSURE"bar_"$TEMPERATURE"K_"$EFIELD"V_"$k"/
        fi
        WATER=n
    fi
    echo "Beginning simulation setup: ${SIMPATH}"

    # check if SIMPATH exists and if not, create the SIMPATH directory and all parent directories if it doesn't exist
    if [ -d "${SIMPATH}" ]; then
        echo "The folder "${SIMPATH}" already exists."
    else
        mkdir -p "${SIMPATH}"
        echo "The folder "${SIMPATH}" has been created."
    fi

    # copy run script if it doesn't exist
    if [ -f "${SIMPATH}/../md_replicates.sh" ]; then
        echo "The file "${SIMPATH}"/../md_replicates.sh already exists."
    else
        cp "${MPATH}"runscripts/md_replicates.sh "${SIMPATH}"/..
        echo "The file "${SIMPATH}"/../md_replicates.sh has been created."
    fi

    # name important files
    NAME="$SIMPATH"ion_solv.gro
    TOP="$SIMPATH"topol.top
    INDEX="$SIMPATH"index.ndx

    # create the ion box and solvate
    CAT_LIST="["
    N_CAT_LIST="["
    for i in "${!N_CAT[@]}"; do
        if [ "$i" -ne 0 ]; then
            CAT_LIST="${CAT_LIST},"
            N_CAT_LIST="${N_CAT_LIST},"
        fi
        CAT_LIST="${CAT_LIST}'${CATION[$i]}'"
        N_CAT_LIST="${N_CAT_LIST}${N_CAT[$i]}"
    done
    CAT_LIST="${CAT_LIST}]"
    N_CAT_LIST="${N_CAT_LIST}]"
    AN_LIST="["
    N_AN_LIST="["
    for i in "${!N_AN[@]}"; do
        if [ "$i" -ne 0 ]; then
            AN_LIST="${AN_LIST},"
            N_AN_LIST="${N_AN_LIST},"
        fi
        AN_LIST="${AN_LIST}'${ANIONS[$i]}'"
        N_AN_LIST="${N_AN_LIST}${N_AN[$i]}"
    done
    AN_LIST="${AN_LIST}]"
    N_AN_LIST="${N_AN_LIST}]"

    python3 -c 'import gro_io; gro_io.ion_box("'${SIMPATH}'",['${BOX_LEN}','${BOX_LEN}','${BOX_LEN}'],'${CAT_LIST}','${AN_LIST}','${N_CAT_LIST}','${N_AN_LIST}',test_n='${k}')'
    gmx solvate -cs "${MPATH}"forcefield/"${WATERMODEL}".gro -cp ${SIMPATH}ions.gro -o ${NAME} -maxsol "${N_SOL}"

    # renumber the atoms and create the index file
    ORDER="["
    for i in "${!CATION[@]}"; do
        ORDER="${ORDER}'${CATION[$i]}',"
    done
    for i in "${!ANIONS[@]}"; do
        ORDER="${ORDER}'${ANIONS[$i]}',"
    done
    ORDER="${ORDER}'SOL']"
    python3 -c 'import gro_io; gro_io.renumber("'${NAME}'","'${NAME}'",'${ORDER}')'
    python3 -c 'import gro_io; gro_io.index("'$NAME'","'$INDEX'","NONE")'

    # clean backup files
    rm "${SIMPATH}"\#*

    # write the number of ions and solvent molecules to the topology file
    echo "Writing the number of ions and solvent molecules to the topology file: ${TOP}"
    if [ -f "${TOP}" ]; then
        rm "${TOP}"
    fi
    if [ "$WATER" == "y" ]; then
        echo "#include \"${MPATH}forcefield/ff${WATERMODEL}.itp\"" >> ${TOP}
    else
        echo "#include \"${MPATH}forcefield/ff${WATERMODEL}-${IONMODEL}.itp\"" >> ${TOP}
    fi
    echo "" >> ${TOP}
    echo "" >> ${TOP}
    echo "[ system ]" >> ${TOP}
    echo "; Name" >> ${TOP}
    echo "Ions in water" >> ${TOP}
    echo "" >> ${TOP}
    echo "[ molecules ]" >> ${TOP}
    echo "; Compound        nmols" >> ${TOP}
    for i in "${!N_CAT[@]}"; do
        if [ "${N_CAT[$i]}" -ne 0 ]; then
            echo ${CATION[$i]}"                 "${N_CAT[$i]} >> ${TOP}
        fi
    done
    for i in "${!N_AN[@]}"; do
        if [ "${N_AN[$i]}" -ne 0 ]; then
            echo ${ANIONS[$i]}"                 "${N_AN[$i]} >> ${TOP}
        fi
    done
    echo "SOL                "${N_SOL} >> ${SIMPATH}topol.top

    # copy the mdp files and replace the placeholders
    cp "${MPATH}"mdp/minim.mdp "${SIMPATH}"
    cp "${MPATH}"mdp/nvt.mdp "${SIMPATH}"
    sed -i "s/!EFIELD!/"${EFIELD}"/g" "${SIMPATH}"nvt.mdp
    sed -i "s/!TEMPERATURE!/"${TEMPERATURE}"/g" "${SIMPATH}"nvt.mdp
    if [ ${NPT} == "y" ]; then
        cp "${MPATH}"mdp/npt.mdp "${SIMPATH}"
        sed -i '/!EFIELD!/d' "${SIMPATH}"npt.mdp
        sed -i "s/!TEMPERATURE!/"${TEMPERATURE}"/g" "${SIMPATH}"npt.mdp
        sed -i "s/!PRESSURE!/"${PRESSURE}"/g" "${SIMPATH}"npt.mdp
        cp "${MPATH}"mdp/md.mdp "${SIMPATH}"
        # calculate the number of steps, which is sim_time * 1000 / 0.002
        NSTEPS=$(echo "scale=0; ${SIM_TIME} * 500000" | bc)
        sed -i "s/!NSTEPS!/"${NSTEPS}"/g" "${SIMPATH}"md.mdp
        sed -i "/!EFIELD!/d" "${SIMPATH}"md.mdp
        sed -i "s/!TEMPERATURE!/"${TEMPERATURE}"/g" "${SIMPATH}"md.mdp
        sed -i "s/!PRESSURE!/"${PRESSURE}"/g" "${SIMPATH}"md.mdp
    else
        cp "${MPATH}"mdp/md_efield.mdp "${SIMPATH}"
        # calculate the number of steps, which is sim_time * 1000 / 0.002
        NSTEPS=$(echo "scale=0; ${SIM_TIME} * 500000" | bc)
        sed -i "s/!NSTEPS!/"${NSTEPS}"/g" "${SIMPATH}"md_efield.mdp
        sed -i "s/!EFIELD!/"${EFIELD}"/g" "${SIMPATH}"md_efield.mdp
        sed -i "s/!TEMPERATURE!/"${TEMPERATURE}"/g" "${SIMPATH}"md_efield.mdp
    fi

    cp "${MPATH}"runscripts/run_md.sh "${SIMPATH}"

    # # determine the job name
    # if [ "$WATER" == "y" ]; then
    #     JOB_NAME="${WATERSTRING}_${k}"
    # else
    #     JOB_NAME="${CATSTRING}${ANSTRING}_${k}"
    # fi

    # # submit the job if RUN is y
    # if [ "${RUN}" == "y" ]; then
    #     echo "submitting job: ${JOB_NAME}"
    #     sbatch --job-name="$JOB_NAME" -t "$SLURMTIME" ion_run.sh "$MPATH" "$SIMPATH" "$INDEX" "$TOP"
    # else
    #     echo "not submitting job: ${JOB_NAME}"
    # fi
done
