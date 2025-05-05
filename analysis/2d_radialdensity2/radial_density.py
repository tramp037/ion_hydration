import numpy as np
import argparse
import sys
import os
import time
import datetime
import math
import MDAnalysis as mda

def main():
    args = parse_arguments()  
    radial_density(args.simpath, args.confname, args.trajname, args.output, args.atoms, args.limit, args.bins, window=args.windows, start=args.start, end=args.end, skip=args.skip)

def radial_density(simpath, confname, trajname, output, atoms, limit, bins, window=[[0,0],[0,0],[0,0]], start=0, end=-1, skip=1):
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


    if not simpath.endswith("/"):
        simpath += "/"

    # gropath = simpath + confname
    gropath = simpath + confname
    # xtcpath = simpath + trajname
    xtcpath = simpath + trajname
    stdout = simpath + output + ".txt"
    outpath = simpath + output + ".xvg"
        
    print("****************************************************************************")
    print("*                                                                          *")
    print("*                             radial_density2()                            *")
    print("*                Calculate the radial density in r and theta               *")
    print("*                                                                          *")
    print("*                         Written By: Naomi Trampe                         *")
    print("*                                SAMPEL Lab                                *")
    print("*                         Last updated: 03/18/2025                         *")
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
        f.write("*                             radial_density2()                            *\n")
        f.write("*                Calculate the radial density in r and theta               *\n")
        f.write("*                                                                          *\n")
        f.write("*                         Written By: Naomi Trampe                         *\n")
        f.write("*                                SAMPEL Lab                                *\n")
        f.write("*                         Last updated: 03/18/2025                         *\n")
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

    dim = system.dimensions[2]


    # binw_x = limit / bins[0]
    # edges_x = np.array([binw_x * i - limit/2 for i in range(bins[0]+1)])
    # centers_x = np.array([binw_x * (i + 0.5) - limit/2 for i in range(bins[0])])
    # binw_y = limit / bins[1]
    # edges_y = np.array([binw_y * i - limit/2 for i in range(bins[1]+1)])
    # centers_y = np.array([binw_y * (i + 0.5) - limit/2 for i in range(bins[1])])
    # cnt = np.zeros((bins[0],bins[1]))
    # total = len(system.trajectory[start:end:skip])

    # # divide the bins
    binw_r = limit / bins[0]
    edges_r = np.array([binw_r * i for i in range(bins[0]+1)])
    centers_r = np.array([binw_r * (i + 0.5) for i in range(bins[0])])
    binw_t = 2 / bins[1]
    edges_t = np.array([binw_t * i - 1 for i in range(bins[1]+1)])
    centers_t = np.array([binw_t * (i + 0.5) - 1 for i in range(bins[1])])
    cnt = np.zeros((bins[0],bins[1]))
    # cnt = [np.zeros(bins) for i in range(len(atoms)-1)]
    total = len(system.trajectory[start:end:skip])
    tot = total
    # num = [total for _ in range(len(atoms)-1)]

    # lambda function to calculate the volume of a spherical shell
    vol = lambda r1, r2: ((2*np.pi) / 3) * binw_t * (r2 ** 3 - r1 ** 3)

    ion_to_analyze=0

    # calculate the radial distribution function
    count = 0
    for ts in system.trajectory[start:end:skip]:
        atom_pos = [i.positions.T.round(2) for i in atom]
        atom_coord = [atom_pos[i+1] for i in range(len(atoms)-1)]
        #atom_coord = [0 for _ in range(len(atoms)-1)]
        box_dims = ts.dimensions[:3]
        # print(box_dims[2]*binw_x*binw_y)

        # print(atom_pos[0][:,0])
        # print(atom_coord[0][:,0])
        num_atoms = 0
        temp_counts = np.zeros((bins[0],bins[1]))
        for j in range(len(atom_coord[ion_to_analyze][0])):
            x2, y2, z2 = atom_coord[ion_to_analyze][:,j]
            num_atoms += 1
            # print("\n")
            # print(x2, z2)
            for i in range(len(atom_pos[0][0])):
                # if mem[0] < atom_coord[0][2][j] < mem[1]:
                x1, y1, z1 = atom_pos[0][:,i]
                dx = min(abs(x2-x1),abs(box_dims[0]-abs(x2-x1)))
                dy = min(abs(y2-y1),abs(box_dims[1]-abs(y2-y1)))
                dz = min(abs(z2-z1),abs(box_dims[2]-abs(z2-z1)))
                dis = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
                # dx = x2 - x1
                # if dx > box_dims[0] / 2:
                #     dx -= box_dims[0]
                # elif dx < -box_dims[0] / 2:
                #     dx += box_dims[0]
                # dy = y2 - y1
                # if dy > box_dims[1] / 2:
                #     dy -= box_dims[1]
                # elif dy < -box_dims[1] / 2:
                #     dy += box_dims[1]
                # idx_x = edges_x.searchsorted(dx,'right')-1
                # if 0 <= idx_x < bins[0]:
                #     idx_y = edges_y.searchsorted(dy,'right')-1
                #     if 0 <= idx_y < bins[1]:
                #         cnt[idx_x][idx_y] += 1/(box_dims[2]*binw_x*binw_y)/total/n_atoms[1]
                idx_r = edges_r.searchsorted(dis,'right')-1
                if idx_r < bins[0]:
                    if z1-z2 > box_dims[2]/2:
                        z1 -= box_dims[2]
                    elif z1-z2 < -box_dims[2]/2:
                        z1 += box_dims[2]
                    idx_t = edges_t.searchsorted((z1-z2)/dis, 'right')-1
                    if idx_t == bins[1]:
                        idx_t -= 1
                    temp_counts[idx_r][idx_t] += 1
        if num_atoms != 0:
            temp_counts /= num_atoms
            cnt += temp_counts
        else:
            tot -= 1
        

        count += 1        
        if total < 100 or count % int(total / 100) == 0:
            print(f"Completing configuration {count}/{total}...")
            with open(stdout, "a") as f:
                f.write(f"Completing configuration {count}/{total}...\n")
            # print(cnt)
            currenttime = time.monotonic()
            elapsedtime = currenttime - starttime
            print(f"Elapsed time: {elapsedtime:.2f} s")
            with open(stdout, "a") as f:
                f.write(f"Elapsed time: {elapsedtime:.2f} s\n")
            remainingtime = (total - count) * (elapsedtime / count)
            print(f"Estimated remaining time: {np.floor(remainingtime/60):.0f} min {remainingtime%60:.0f} s\n")
            with open(stdout, "a") as f:
                f.write(f"Estimated remaining time: {np.floor(remainingtime/60):.0f} min {remainingtime%60:.0f} s\n\n")




        # # determine the mask of each first atom position
        # m = [0,0,0]
        # for a in range(len(atoms)-1):
        #     for i in range(3):
        #         for j in range(int(len(window[i])/2)):
        #             # if the window has is a single value, set the mask to all True
        #             if window[i][2*j]-window[i][2*j+1] == 0:
        #                 arr = np.ones_like(atom_pos[a+1][i], dtype=bool)
        #             # otherwise, set the mask to True if the atom is inside
        #             elif window[i][0] > window[i][1]:
        #                 arr = (atom_pos[a+1][i] >= window[i][2*j]) | (atom_pos[a+1][i] < window[i][2*j+1])
        #             else:
        #                 arr = (atom_pos[a+1][i] >= window[i][2*j]) & (atom_pos[a+1][i] < window[i][2*j+1])
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

        # # count all of the atoms in each bin, dividing by the volume of the system
        # for i in range(len(atoms)-1):
        #     if n_atoms[i+1] != 0:
        #         counts[i] += 1
        #     for j in range(n_atoms[i+1]):
        #         # calculate the periodic distance
        #         coordinate = atom_coord[i][:,j]
        #         dx = np.min([np.abs(coordinate[0]-atom_pos[0][0]),np.abs(box_dims[0]-np.abs(coordinate[0]-atom_pos[0][0]))], axis=0)
        #         dy = np.min([np.abs(coordinate[1]-atom_pos[0][1]),np.abs(box_dims[1]-np.abs(coordinate[1]-atom_pos[0][1]))], axis=0)
        #         dz = np.min([np.abs(coordinate[2]-atom_pos[0][2]),np.abs(box_dims[2]-np.abs(coordinate[2]-atom_pos[0][2]))], axis=0)
        #         dis = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
                
        #         # find the bin the distance is in
        #         idx = edges.searchsorted(dis,'right')-1
        #         for k in idx:
        #             if k < bins:
        #                 cnt[i][k] += 1 * (box_dims[0] * box_dims[1] * box_dims[2]) / n_atoms[i+1]
        
        # # print progress
        # if total < 100 or count % int(total / 100) == 0:
        #     print(f"Completing configuration {count}/{total}...")
        #     with open(stdout, "a") as f:
        #         f.write(f"Completing configuration {count}/{total}...\n")
        #     last = [cnt[i][-1] / vol(edges[-2],edges[-1]) / n_atoms[0] / (counts[i]) for i in range(len(atoms)-1)]
        #     out = ""
        #     for i, a in enumerate(atoms[1:]):
        #         out += f"{a}-{atoms[0]}: {last[i]:.2e}     "
        #     print(out + "\n")
        #     with open(stdout, "a") as f:
        #         f.write(out + "\n\n")
        # count += 1

    # preprocess cnt
    for i in range(bins[0]):
        for j in range(bins[1]):
            if tot == 0:
                sys.exit("No atoms found in the window")
            cnt[i][j] = cnt[i][j]/vol(edges_r[i],edges_r[i+1])/tot
            # cnt[i][j] /= vol

    # stop the timer and print the time elapsed
    endtime = time.monotonic()
    totaltime = (endtime - starttime)
    print(tot)
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
        f.write(f"{bins[0]} {bins[1]}\n")
        for i in centers_r:
        # for i in centers_x:
            f.write(f"{i} ")
        f.write("\n")
        for i in range(len(centers_t)):
            f.write(f"{centers_t[i]} ")
            for j in range(len(centers_r)):
        # for i in range(len(centers_y)):
            # f.write(f"{centers_y[i]} ")
            # for j in range(len(centers_x)):
                f.write(f"{cnt[j][i]} ")
            f.write("\n")
        
    #     for binn in range(bins):
    #         f.write(f"{edges[binn] + binw / 2} ")
    #         for i in range(len(atoms)-1):
    #             # write the count of the bin, dividing by the volume of the bin, the number of atoms of each type, and the number of configurations analyzed
    #             f.write(f"{cnt[i][binn] / vol(edges[binn],edges[binn + 1]) / n_atoms[0] / (counts[i] + 1)} ")
    #         f.write("\n")


def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate radial distribution function")
    parser.add_argument("-p", "--simpath", type=str, default="./", help="Path to simulation directory")
    parser.add_argument("-c", "--confname", type=str, default="conf.gro", help="Name of configuration file")
    parser.add_argument("-t", "--trajname", type=str, default="traj.xtc", help="Name of trajectory file")
    parser.add_argument("-o", "--output", type=str, default="rdf", help="Output file name")
    parser.add_argument("-a", "--atoms", type=str, nargs="+", help="Name of atoms to calculate RDF of (order: atom2 atom1_1 atom1_2 ...)")
    parser.add_argument("-l", "--limit", type=float, default=10, help="Limit of RDF")
    parser.add_argument("-b", "--bins", type=int, nargs=2, default=[100, 10], help="Number of bins (radial, angular)")
    parser.add_argument("-w", "--windows", type=float, nargs=6, default=[[0,0],[0,0],[0,0]], help="Windows to calculate RDF within")
    parser.add_argument("--start", type=int, default=0, help="Starting frame")
    parser.add_argument("--end", type=int, default=-1, help="Ending frame")
    parser.add_argument("--skip", type=int, default=1, help="Frame interval")  

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    main()

# cation = "Cs"

# # Define paths here
# #MPATH = "/home/sarupria/tramp037/zeo-ion/LTAI/uc1/Cl/Li/K/test_1/"
# MPATH = "/scratch.global/tramp037/zeolite_ion_selectivity/simulations/ions/Cl"+cation+"/run_0/"
# #GROPATH = MPATH + "zeo_md.gro"
# GROPATH = MPATH + "ions_md.gro"
# #XTCPATH = MPATH + "zeo_md.xtc"
# XTCPATH = MPATH + "ions_md.xtc"
# STDOUT = MPATH + "rdf.txt"
# OUTPATH = MPATH + "RDF_"+cation+"_OW.xvg"

# bulk = [61.6,9.92]
# surface = [11.42,12.42,23.34,24.34,47.18,48.18,59.1,60.1]
# pore = [13.92,21.84,49.68,57.6]


# #COORDINATES = [35.76, 35.76, 71.52]
# #WINDOW = [[0,0],[0,0],surface]
# WINDOW=[[0,0],[0,0],[0,0]]
# WALL = [[0,0],[0,0],[0,0]]

# # Select atoms here
# # Options:
# #       "name <ion>"              - ion
# #       "name OW"                 - water oxygen
# #       "name O* and not name OW" - zeolite oxygen
# ATOM1 = "name "+cation
# ATOM2 = "name OW"

# # Define limit of the RDF
# ATOM1 = "name "+cation
# ATOM2 = "name OW"
# LIM = 10
# BINS = int(LIM) * 10

# START = 0
# END = -1
# SKIP = 10


# def main():
#     rdf()


# def rdf():
#     binw = LIM / BINS
#     with open(STDOUT, "w") as f:
#         f.write("Beginning calculation...\n")
#     univ = mda.Universe(GROPATH, XTCPATH)
#     with open(STDOUT, "a") as f:
#         f.write("Loaded configuration and trajectory...\n")

#     one = univ.select_atoms(ATOM1)
#     two = univ.select_atoms(ATOM2)
#     total = len(univ.trajectory[START:END:SKIP])
    
#     one_pos = np.zeros((3, total, len(one)))
#     two_pos = np.zeros((3, total, len(two)))
#     box_dims = np.zeros((3, total))
    
#     starttime = time.monotonic()
    
#     for t, ts in enumerate(univ.trajectory[START:END:SKIP]):
#         one_pos[:, t] = one.positions.T
#         two_pos[:, t] = two.positions.T
#         box_dims[:, t] = ts.dimensions[:3]

#     b_list = np.linspace(0, LIM, BINS)
#     hist = np.zeros(len(b_list))
#     count = 0
#     num = total

#     vol = lambda r1, r2: (4 / 3) * np.pi * (r2 ** 3 - r1 ** 3)
    
#     # for each timestep
#     for t in range(total):
#         # determine mask for the first atom position
#         m = [0,0,0]
#         for i in range(3):
#             for j in range(int(len(WINDOW[i])/2)):
#                 # if the window has is a single value, set the mask to all True
#                 if WINDOW[i][2*j]-WINDOW[i][2*j+1] == 0:
#                     arr = np.ones_like(one_pos[:, t][i], dtype=bool)
#                 # otherwise, set the mask to True if the atom is inside
#                 elif WINDOW[i][0] > WINDOW[i][1]:
#                     arr = (one_pos[:, t][i] > WINDOW[i][2*j]) | (one_pos[:, t][i] < WINDOW[i][2*j+1])
#                 else:
#                     arr = (one_pos[:, t][i] > WINDOW[i][2*j]) & (one_pos[:, t][i] < WINDOW[i][2*j+1])    
#                 # if it's the first window, set the mask to the array
#                 if j == 0:
#                     m[i] = arr
#                 # if it's not the first window, combine the mask with the previous mask
#                 else:
#                     m[i] = np.logical_or(m[i],arr)
#         mask = np.array(m)
#         # make sure the atom is inside the x, y, and z windows
#         onemask = np.all(mask, axis=0)
#         # number of atoms inside the window
#         len1 = np.sum(onemask)
#         # if there are no atoms inside the window, subtract from the total number of configurations
#         if len1 == 0:
#             num -= 1
#         # coordinates of the first atom inside the window
#         one_coords = one_pos[:, t][:, onemask]
        
#         # determine mask for the second atom position
#         m = [0,0,0]
#         for i in range(3):
#             for j in range(int(len(WALL[i])/2)):
#                 if WALL[i][2*j]-WALL[i][2*j+1] == 0:
#                     arr = np.ones_like(two_pos[:, t][i], dtype=bool)
#                 elif WALL[i][0] > WALL[i][1]:
#                     arr = (two_pos[:, t][i] > WALL[i][2*j]) | (two_pos[:, t][i] < WALL[i][2*j+1])
#                 else:
#                     arr = (two_pos[:, t][i] > WALL[i][2*j]) & (two_pos[:, t][i] < WALL[i][2*j+1])    
#                 if j == 0:
#                     m[i] = arr
#                 else:
#                     m[i] = np.logical_or(m[i],arr)
#         mask = np.array(m)
#         twomask = np.all(mask, axis=0)
#         len2 = np.sum(twomask)
#         if len2 == 0:
#             num -= 1
#         two_coords = two_pos[:, t][:, twomask]

#         for i in range(len(one_coords[0])):
#             x1, y1, z1 = one_coords[:, i]

#             for j in range(len(two_coords[0])):
#                 x2, y2, z2 = two_coords[:, j]
#                 dx = min(abs(x2-x1),abs(box_dims[0, t]-abs(x2-x1)))
#                 dy = min(abs(y2-y1),abs(box_dims[1, t]-abs(y2-y1)))
#                 dz = min(abs(z2-z1),abs(box_dims[2, t]-abs(z2-z1)))

#                 dis = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

#                 for b in range(len(b_list)):
#                     if b_list[b-1] <= dis < b_list[b]:
#                         hist[b] += 1 / (len1) * (box_dims[0, 0] * box_dims[1, 0] * box_dims[2, 0])


#         if total < 100:
#             with open(STDOUT, "a") as f:
#                 f.write(f"Completing configuration {count}/{total}...\n")
#         elif count % int(total / 100) == 0:
#             with open(STDOUT, "a") as f:
#                 f.write(f"Completing configuration {count}/{total}...\n")
#         count += 1

#     for i in range(1, len(b_list)):
#         hist[i] /= num
#         hist[i] /= len(two)
#         hist[i] /= vol(b_list[i - 1], b_list[i])

#     with open(STDOUT, "a") as f:
#         f.write("Calculation complete. Writing result...\n")

#     with open(OUTPATH, "w") as f:
#         for i in range(len(b_list)):
#             f.write(f"{b_list[i]} {hist[i]}\n")

#     endtime = time.monotonic()
#     totaltime = (endtime - starttime) / 60
#     with open(STDOUT, "a") as f:
#         f.write(f"Calculation completed in {totaltime:.2f} min\n")


# if __name__ == "__main__":
#     main()
