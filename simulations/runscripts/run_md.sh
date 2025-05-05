#!/bin/bash
#SBATCH --ntasks=4
#SBATCH --mem=40g
#SBATCH --tmp=40g
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tramp037@umn.edu
source /home/sarupria/shared/software/load_scripts/load_gromacs-2022.sh
export OMP_NUM_THREADS=4

NPT=$1

if [ ! -f ions_em.gro ]; then
    # running energy minimization
    gmx grompp -f minim.mdp -c ion_solv.gro -p topol.top -n index.ndx -o ions_em.tpr
    gmx mdrun -deffnm ions_em -ntmpi 1 -ntomp 4

    rm *step*
fi

# check to see if the energy minimization was successful
if [ ! -f ions_em.gro ]; then
    echo "ions_em.gro does not exist"
    exit 1
fi

# running NVT equilibration
gmx grompp -f nvt.mdp -c ions_em.gro -p topol.top -n index.ndx -o ions_nvt.tpr
gmx mdrun -deffnm ions_nvt -ntmpi 1 -ntomp 4

# check to see if the NVT equilibration was successful
if [ ! -f ions_nvt.gro ]; then
    echo "ions_nvt.gro does not exist"
    exit 1
fi

if [ $NPT == "y" ]; then

    # running NPT equilibration
    gmx grompp -f npt.mdp -c ions_nvt.gro -p topol.top -n index.ndx -o ions_npt.tpr -maxwarn 1
    gmx mdrun -deffnm ions_npt -ntmpi 1 -ntomp 4

    # check to see if the NPT equilibration was successful
    if [ ! -f ions_npt.gro ]; then
        echo "ions_npt.gro does not exist"
        exit 1
    fi

    # running production MD
    gmx grompp -f md.mdp -c ions_npt.gro -p topol.top -n index.ndx -o ions_md.tpr
    gmx mdrun -deffnm ions_md -ntmpi 1 -ntomp 4

    # check to see if the production MD was successful
    if [ ! -f ions_md.gro ]; then
        echo "ions_md.gro does not exist"
        exit 1
    else
        rm \#*
        echo "Simulation complete"
    fi

else

    # running production MD
    gmx grompp -f md_efield.mdp -c ions_nvt.gro -p topol.top -n index.ndx -o ions_md.tpr
    gmx mdrun -deffnm ions_md -ntmpi 1 -ntomp 4

    # check to see if the production MD was successful
    if [ ! -f ions_md.gro ]; then
        echo "ions_md.gro does not exist"
        exit 1
    else
        rm \#*
        echo "Simulation complete"
    fi
fi

