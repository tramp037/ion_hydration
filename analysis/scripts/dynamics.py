import numpy as np
import argparse
import sys
import os
import time
import math
import MDAnalysis as mda

def main():
    if sys.argv[1] == "msd":
        parser = argparse.ArgumentParser(description="Calculate mean square displacement of atoms along specified axes")
        parser.add_argument("density", type=str, help="Function to run")
        parser.add_argument("-p", "--simpath", type=str, default="./", help="Path to simulation directory")
        parser.add_argument("-c", "--confname", type=str, default="conf.gro", help="Name of configuration file")
        parser.add_argument("-t", "--trajname", type=str, default="traj.xtc", help="Name of trajectory file")
        parser.add_argument("-o", "--output", type=str, default="density", help="Output file name")
        parser.add_argument("-a", "--atoms", type=str, nargs="+", help="Name of atoms to calculate msd of")
        parser.add_argument("-x", "--axis", type=str, nargs="+", help="Axes of msd (x, y, z, or any combination)")
        parser.add_argument("--start", type=int, default=0, help="Starting frame")
        parser.add_argument("--end", type=int, default=-1, help="Ending frame")
        parser.add_argument("--skip", type=int, default=1, help="Frame interval")
    
    else:
        parser = argparse.ArgumentParser(description="Dynamical analysis of molecular dynamics simulations. Current functions: msd")
        parser.add_argument("function", type=str, help="Function to run")

    args = parser.parse_args()

    if sys.argv[1] == "msd":
        msd(args.simpath, args.confname, args.trajname, args.output, args.atoms, args.axis, args.start, args.end, args.skip)

def msd(simpath, confname, trajname, output, atoms, axis, start, end, skip):

    if not simpath.endswith("/"):
        simpath += "/"

    gropath = simpath + confname
    xtcpath = simpath + trajname
    stdout = simpath + output + ".txt"
    outpath = simpath + output + ".xvg"
        
    print("****************************************************************************")
    print("*                                                                          *")
    print("*                              dynamics.msd()                              *")
    print("*             Calculate the mean squared displacement of atoms             *")
    print("*                                                                          *")
    print("*                         Written By: Naomi Trampe                         *")
    print("*                                SAMPEL Lab                                *")
    print("*                         Last updated: 10/10/2023                         *")
    print("*                                                                          *")
    print("****************************************************************************\n")       

    # print command line input
    command = sys.argv[0]
    for arg in sys.argv[1:]:
        command += f" {arg}"
    print("Command line input:     " + command + "\n\n")

    
    with open(stdout, "w") as f:
        f.write("****************************************************************************\n")
        f.write("*                                                                          *\n")
        f.write("*                              dynamics.msd()                              *\n")
        f.write("*             Calculate the mean squared displacement of atoms             *\n")
        f.write("*                                                                          *\n")
        f.write("*                         Written By: Naomi Trampe                         *\n")
        f.write("*                                SAMPEL Lab                                *\n")
        f.write("*                         Last updated: 10/10/2023                         *\n")
        f.write("*                                                                          *\n")
        f.write("****************************************************************************\n\n")

        # print command line input
        command = sys.argv[0]
        for arg in sys.argv[1:]:
            command += f" {arg}"
        f.write("Command line input:     " + command + "\n\n")

    # exit if the simulation path does not exist
    if not os.path.isdir(simpath):
        print(f"Simulation path {simpath} does not exist\n")
        sys.exit(1)

    # exit if file gropath or xtcpath do not exist
    try:
        with open(gropath, "r") as f:
            pass
    except FileNotFoundError:
        print(f"File {gropath} not found\n")
        with open(stdout, "a") as f:
            f.write(f"File {gropath} not found\n")
        sys.exit(1)
    try:
        with open(xtcpath, "r") as f:
            pass
    except FileNotFoundError:
        print(f"File {xtcpath} not found\n")
        with open(stdout, "a") as f:
            f.write(f"File {xtcpath} not found\n")
        sys.exit(1)
    print(f"Using configuration:   {gropath}")
    print(f"Using trajectory:      {xtcpath}\n")
    with open(stdout, "a") as f:
        f.write(f"Using configuration:   {gropath}\n")
        f.write(f"Using trajectory:      {xtcpath}\n\n")
    
    # exit if atoms is not specified
    if atoms is None:
        print("Atoms not specified\n")
        with open(stdout, "a") as f:
            f.write("Atoms not specified\n")
        sys.exit(1)
    
    # create the mdanalysis universe
    system = mda.Universe(gropath, xtcpath)
    
    print("Loaded configuration and trajectory\n")
    with open(stdout, "a") as f:
        f.write("Loaded configuration and trajectory\n\n")

    # start the timer
    starttime = time.monotonic()
    
    # select the atoms of interest
    system.trajectory[start]
    atom = [system.select_atoms("name " + i) for i in atoms]
    # collect the initial positions of the atoms
    
    if len(axis) == 1:
        if axis[0] == "x":
            ax = 0
        elif axis[0] == "y":
            ax = 1
        elif axis[0] == "z":
            ax = 2
    else:
        print("Multiple axis not currently supported. Exiting...\n")
        with open(stdout, "a") as f:
            f.write("Multiple axis not currently supported. Exiting...\n\n")
        sys.exit(1) 

    dim = [round(system.dimensions[i],2) for i in range(3)]
    
    # calculate the msd
    msd = [[] for i in range(len(atoms))]
    time_values = []
    count = 0
    total = len(system.trajectory[start:end:skip])
    for ts in system.trajectory[start:end:skip]:
        time_values.append(ts.time/1000)
        if count != 0:
            last_pos = atom_pos
            atom_pos = [i.positions[:, ax] for i in atom]
            # for i in range(len(atoms)):
            #     for j in range(len(atom_pos[i])):
            #         if atom_pos[i][j] - last_pos[i][j] > 2 * dim[ax] / 3:
            #             #print(f"Moving atom {j} of {atoms[i]} at step {count}:")
            #             #print(last_pos[i][j],atom_pos[i][j])
            #             atom_0[i][j] += dim[ax]
            #             #print(atom_0[i][j] - dim[ax],atom_0[i][j],"\n")
            #         elif atom_pos[i][j] - last_pos[i][j] < -2 * dim[ax] / 3:
            #             #print(f"Moving atom {j} of {atoms[i]} at step {count}:")
            #             #print(last_pos[i][j],atom_pos[i][j])
            #             atom_0[i][j] -= dim[ax]
            #             #print(atom_0[i][j] + dim[ax],atom_0[i][j],"\n")
            dis_sum = [(i - j) ** 2 for i, j in zip(atom_pos, atom_0)]
            for i in range(len(atoms)):
                msd[i].append(np.average(dis_sum[i]))
                #msd[i].append(max(dis_sum[i]))
        else:
            atom_0 = [i.positions[:, ax] for i in atom]
            atom_pos = atom_0
            for i in range(len(atoms)):
                msd[i].append(0)
        
        # print progress
        if total < 100 or count % int(total / 100) == 0:
            print(f"Completing configuration {count}/{total}...")
            with open(stdout, "a") as f:
                f.write(f"Completing configuration {count}/{total}...\n")
            out = ""
            for i, a in enumerate(atoms):
                out += f"{a}: {msd[i][-1]:.2e}     "
            print(out + "\n")
            with open(stdout, "a") as f:
                f.write(out + "\n\n")
        count += 1

    # stop the timer and print the time elapsed
    endtime = time.monotonic()
    totaltime = (endtime - starttime)
    print(f"\nCalculation completed in {totaltime / 60:.2f} min\n")
    print(f"Processing time: {totaltime / total:.2e} s/frame")
    print(f"                 {total / totaltime * 60:.2f} frames/min\n")
    with open(stdout, "a") as f:
        f.write(f"\nCalculation completed in {totaltime / 60:.2f} min\n\n")
        f.write(f"Processing time: {totaltime / total:.2e} s/frame\n")
        f.write(f"                 {total / totaltime * 60:.2f} frames/min\n\n")
    
    print(f"Writing result to {outpath}\n")
    with open(outpath, "w") as f:
        for i in range(len(time_values)):
            f.write(f"{time_values[i]}")
            for j in range(len(atoms)):
                f.write(f" {msd[j][i]}")
            f.write("\n")

if __name__ == "__main__":
    main()
