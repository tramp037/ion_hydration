import numpy as np
import argparse
import sys
import os
import time
import datetime
import math
import MDAnalysis as mda

def main():
    if sys.argv[1] == "density":
        args = parse_density()
        density(args.simpath, args.confname, args.trajname, args.output, args.atoms, args.axis, start=args.start, end=args.end, skip=args.skip)

    elif sys.argv[1] == "rdf":
        args = parse_rdf()  
        rdf(args.simpath, args.confname, args.trajname, args.output, args.atoms, args.limit, args.bins, args.blocks, xwindow=args.xwindows, ywindow=args.ywindows, zwindow=args.zwindows, start=args.start, end=args.end, skip=args.skip)

    else:
        parser = argparse.ArgumentParser(description="Structural analysis of molecular dynamics simulations. Current functions: density, rdf")
        parser.add_argument("function", type=str, help="Function to run")
        args = parser.parse_args()

def density(simpath, confname, trajname, output, atoms, axis, start=0, end=-1, skip=1):
    """
    Calculate the density profile of atoms along a given axis

    Parameters
    ----------
    simpath : str
        Path to simulation directory
    confname : str
        Name of configuration file
    trajname : str
        Name of trajectory file
    output : str
        Output file name
    atoms : list
        Name of atoms to calculate density of
    axis : int
        Axis of density profile (0: x, 1: y, 2: z)
    start : int, optional
        Starting frame
    end : int, optional
        Ending frame
    skip : int, optional
        Frame interval
    """

    if not simpath.endswith("/"):
        simpath += "/"

    gropath = simpath + confname
    xtcpath = simpath + trajname
    stdout = simpath + output + ".txt"
    outpath = simpath + output + ".xvg"
        
    print("****************************************************************************")
    print("*                                                                          *")
    print("*                           structure.density()                            *")
    print("*                  Calculate the density profile of atoms                  *")
    print("*                                                                          *")
    print("*                         Written By: Naomi Trampe                         *")
    print("*                                SAMPEL Lab                                *")
    print("*                         Last updated: 05/06/2024                         *")
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
        f.write("*                           structure.density()                            *\n")
        f.write("*                  Calculate the density profile of atoms                  *\n")
        f.write("*                                                                          *\n")
        f.write("*                         Written By: Naomi Trampe                         *\n")
        f.write("*                                SAMPEL Lab                                *\n")
        f.write("*                         Last updated: 05/06/2024                         *\n")
        f.write("*                                                                          *\n")
        f.write("****************************************************************************\n\n")

        # print the date and time
        current_time = str(datetime.datetime.now())
        f.write("Running at: " + current_time + "\n\n")

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

    # find the minimum box size of the system along the axis of interest through the trajectory
    coordinates = [i.dimensions[axis] for i in system.trajectory[start:end:skip]]
    coord = round(min(coordinates),2)
    
    # select the atoms of interest
    system.trajectory[0]
    atom = [system.select_atoms("name " + i) for i in atoms]

    # divide the box into bins
    bins = math.floor(coord * 2) - 1
    binw = 0.5
    center_0 = round((1 - bins) * binw / 2,2)
    centers = np.array([center_0 + binw * i for i in range(bins)])
    # round the centers array
    centers = np.round(centers,2)
    edge_0 = round(center_0 - binw / 2,2)
    edges = np.array([edge_0 + binw * i for i in range(bins+1)])
    # round the edges array
    edges = np.round(edges,2)
    # initialize the count of atoms in each bin
    cnt = [np.zeros(bins) for i in range(len(atoms))]
    # calculate the total number of frames
    total = len(system.trajectory[start:end:skip])
    interval = 1 / binw

    # calculate the density profile
    count = 0
    for ts in system.trajectory[start:end:skip]:
        atom_pos = [i.positions[:, axis].round(2) - coord/2 for i in atom]
        # count all of the atoms in each bin, dividing by the volume of the bin
        for i in range(len(atoms)):      
            # find the bin the distance is in
            idx = edges.searchsorted(atom_pos[i],'right')
            for k in idx:
                if 0 < k < bins+1:
                    cnt[i][k-1] += interval / ts.dimensions[(axis + 1) % 3] / ts.dimensions[(axis + 2) % 3]
        
        # print progress
        if total < 100 or count % int(total / 100) == 0:
            print(f"Completing configuration {count}/{total}...")
            with open(stdout, "a") as f:
                f.write(f"Completing configuration {count}/{total}...\n")
            avg = [np.average(i)/(count+1) for i in cnt]
            out = ""
            for i, a in enumerate(atoms):
                out += f"{a}: {avg[i]:.2e}     "
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
    
    # write the result to the output file
    print(f"Writing result to {outpath}\n")
    with open(outpath, "w") as f:
        for binn in range(bins):
            f.write(f"{edges[binn] + binw / 2} ")
            for i in range(len(atoms)):
                f.write(f"{cnt[i][binn]/count} ")
            f.write("\n")

def rdf(simpath, confname, trajname, output, atoms, limit, bins, blocks, xwindow=[0,0], ywindow=[0,0], zwindow=[0,0], start=0, end=-1, skip=1):
    """
    Calculate the radial distribution function

    Parameters
    ----------
    simpath : str
        Path to simulation directory
    confname : str
        Name of configuration file
    trajname : str
        Name of trajectory file
    output : str
        Output file name
    atoms : list
        Name of atoms to calculate RDF of (order: atom2 atom1_1 atom1_2 ...)
    limit : float
        Limit of RDF
    bins : int
        Number of bins
    window : list, optional
        Windows to calculate RDF within
    start : int, optional
        Starting frame
    end : int, optional
        Ending frame
    skip : int, optional
        Frame interval
    """
    windows = [xwindow, ywindow, zwindow]

    if not simpath.endswith("/"):
        simpath += "/"

    gropath = simpath + confname
    xtcpath = simpath + trajname
    stdout = simpath + output + ".txt"
    outpath = simpath + output + ".xvg"
        
    print("****************************************************************************")
    print("*                                                                          *")
    print("*                             structure.rdf()                              *")
    print("*                Calculate the radial distribution function                *")
    print("*                                                                          *")
    print("*                         Written By: Naomi Trampe                         *")
    print("*                                SAMPEL Lab                                *")
    print("*                         Last updated: 09/11/2024                         *")
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
        f.write("*                             structure.rdf()                              *\n")
        f.write("*                Calculate the radial distribution function                *\n")
        f.write("*                                                                          *\n")
        f.write("*                         Written By: Naomi Trampe                         *\n")
        f.write("*                                SAMPEL Lab                                *\n")
        f.write("*                         Last updated: 09/11/2024                         *\n")
        f.write("*                                                                          *\n")
        f.write("****************************************************************************\n\n")

        # print the date and time
        current_time = str(datetime.datetime.now())
        f.write("Running at: " + current_time + "\n\n")

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
    system.trajectory[0]
    atom = [system.select_atoms("name " + i + "*") for i in atoms]
    n_atoms = [len(i) for i in atom]

    # divide the bins
    binw = limit / bins
    edges = np.array([binw * i for i in range(bins+1)])
    cnt = [[np.zeros(bins) for i in range(len(atoms)-1)] for k in range(blocks)]
    total = len(system.trajectory[start:end:skip])
    block_total = int(total / blocks)

    # lambda function to calculate the volume of a spherical shell
    vol = lambda r1, r2: (4 / 3) * np.pi * (r2 ** 3 - r1 ** 3)

    counts = [[0 for _ in range(len(atoms)-1)] for k in range(blocks)]
    # calculate the radial distribution function
    for block in range(blocks):
        block_start = block * block_total * skip
        block_end = (block + 1) * block_total * skip
        count = 0
        for ts in system.trajectory[block_start:block_end:skip]:
            atom_pos = [i.positions.T.round(2) for i in atom]
            atom_coord = [atom_pos[i+1] for i in range(len(atoms)-1)]
            #atom_coord = [0 for _ in range(len(atoms)-1)]
            box_dims = ts.dimensions[:3]


            # determine the mask of each first atom position
            # m = [0,0,0]
            # for a in range(len(atoms)-1):
            #     for i in range(3):
            #         for j in range(int(len(windows[i])/2)):
            #             # if the window has is a single value, set the mask to all True
            #             # print(windows[i][2*j], windows[i][2*j+1])
            #             if windows[i][2*j]-windows[i][2*j+1] == 0:
            #                 arr = np.ones_like(atom_pos[a+1][i], dtype=bool)
            #             # otherwise, set the mask to True if the atom is inside
            #             elif windows[i][0] > windows[i][1]:
            #                 arr = (atom_pos[a+1][i] >= windows[i][2*j]) | (atom_pos[a+1][i] < windows[i][2*j+1])
            #             else:
            #                 arr = (atom_pos[a+1][i] >= windows[i][2*j]) & (atom_pos[a+1][i] < windows[i][2*j+1])
            #             # if it's the first window, set the mask to the array
            #             if j == 0:
            #                 m[i] = arr
            #             # otherwise, combine the mask with the previous mask
            #             else:
            #                 m[i] = np.logical_or(m[i],arr)

            #     mask = np.array(m)
            #     # make sure the atom is inside the x, y, and z windows
            #     onemask = np.all(mask, axis=0)
            #     # number of atoms inside the window
            #     n_atoms[a+1] = np.sum(onemask)
            #     atom_coord[a] = atom_pos[a+1][:, onemask]

            # count all of the atoms in each bin, dividing by the volume of the system
            for i in range(len(atoms)-1):
                if n_atoms[i+1] != 0:
                    counts[block][i] += 1
                min = 200
                for j in range(n_atoms[i+1]):
                    # calculate the periodic distance
                    coordinate = atom_coord[i][:,j]
                    # print(coordinate)
                    # if coordinate[2] > 50 and coordinate[2] < min:
                    #     min = coordinate[2]
                    dx = np.min([np.abs(coordinate[0]-atom_pos[0][0]),np.abs(box_dims[0]-np.abs(coordinate[0]-atom_pos[0][0]))], axis=0)
                    dy = np.min([np.abs(coordinate[1]-atom_pos[0][1]),np.abs(box_dims[1]-np.abs(coordinate[1]-atom_pos[0][1]))], axis=0)
                    dz = np.min([np.abs(coordinate[2]-atom_pos[0][2]),np.abs(box_dims[2]-np.abs(coordinate[2]-atom_pos[0][2]))], axis=0)
                    dis = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
                    
                    # find the bin the distance is in
                    idx = edges.searchsorted(dis,'right')-1
                    for k in idx:
                        if k < bins:
                            cnt[block][i][k] += 1 * (box_dims[0] * box_dims[1] * box_dims[2]) / n_atoms[i+1]
                # print(min)
            # print progress
            if total < 100 or count % int(total / 100) == 0:
                print(f"Completing configuration {count}/{block_total} of block {block+1}/{blocks}...")
                with open(stdout, "a") as f:
                    f.write(f"Completing configuration {count}/{block_total} of block {block+1}/{blocks}...\n")
                last = [cnt[block][i][-1] / vol(edges[-2],edges[-1]) / n_atoms[0] / (counts[block][i]) for i in range(len(atoms)-1)]
                out = ""
                for i, a in enumerate(atoms[1:]):
                    out += f"{a}-{atoms[0]}: {last[i]:.2e}     "
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
    
    # write the result to the output file
    print(f"Writing result to {outpath}\n")
    with open(outpath, "w") as f:
        for binn in range(bins):
            f.write(f"{edges[binn] + binw / 2} ")
            for block in range(blocks):
                for i in range(len(atoms)-1):
                    # write the count of the bin, dividing by the volume of the bin, the number of atoms of each type, and the number of configurations analyzed
                    f.write(f"{cnt[block][i][binn] / vol(edges[binn],edges[binn + 1]) / n_atoms[0] / (counts[block][i] + 1)} ")
            f.write("\n")
    
def parse_density():
    parser = argparse.ArgumentParser(description="Calculate density of atoms along a given axis")
    parser.add_argument("density", type=str, help="Function to run")
    parser.add_argument("-p", "--simpath", type=str, default="./", help="Path to simulation directory")
    parser.add_argument("-c", "--confname", type=str, default="conf.gro", help="Name of configuration file")
    parser.add_argument("-t", "--trajname", type=str, default="traj.xtc", help="Name of trajectory file")
    parser.add_argument("-o", "--output", type=str, default="density", help="Output file name")
    parser.add_argument("-a", "--atoms", type=str, nargs="+", help="Name of atoms to calculate density of")
    parser.add_argument("-x", "--axis", type=int, default=0, help="Axis of density profile")
    parser.add_argument("--start", type=int, default=0, help="Starting frame")
    parser.add_argument("--end", type=int, default=-1, help="Ending frame")
    parser.add_argument("--skip", type=int, default=1, help="Frame interval")

    args = parser.parse_args()
    return args

def parse_rdf():
    parser = argparse.ArgumentParser(description="Calculate radial distribution function")
    parser.add_argument("rdf", type=str, help="Function to run")
    parser.add_argument("-p", "--simpath", type=str, default="./", help="Path to simulation directory")
    parser.add_argument("-c", "--confname", type=str, default="conf.gro", help="Name of configuration file")
    parser.add_argument("-t", "--trajname", type=str, default="traj.xtc", help="Name of trajectory file")
    parser.add_argument("-o", "--output", type=str, default="rdf", help="Output file name")
    parser.add_argument("-a", "--atoms", type=str, nargs="+", help="Name of atoms to calculate RDF of (order: atom2 atom1_1 atom1_2 ...)")
    parser.add_argument("-l", "--limit", type=float, default=10, help="Limit of RDF")
    parser.add_argument("-b", "--bins", type=int, default=100, help="Number of bins")
    parser.add_argument("-k", "--blocks", type=int, default=1, help="Number of blocks")
    parser.add_argument("-x", "--xwindows", type=float, nargs="+", default=[0,0], help="X-dim Windows to calculate RDF within")
    parser.add_argument("-y", "--ywindows", type=float, nargs="+", default=[0,0], help="Y-dim Windows to calculate RDF within")
    parser.add_argument("-z", "--zwindows", type=float, nargs="+", default=[0,0], help="Z-dim Windows to calculate RDF within")
    parser.add_argument("--start", type=int, default=0, help="Starting frame")
    parser.add_argument("--end", type=int, default=-1, help="Ending frame")
    parser.add_argument("--skip", type=int, default=1, help="Frame interval")  

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    main()
