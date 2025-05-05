import numpy as np 
import pandas as pd
import argparse
import sys
import time
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase

class Diffusivity(AnalysisBase):
    """
    Calculate the MSD of the provided molecules in the system

    Parameters
    ----------
    gropath : str
        Path to GRO file
    xtcpath : str
        Path to XTC file
    atomss : list
        List of different atoms to calculate the MSD for
    bins : int
        Number of bins to divide the angle range into

    """
    def __init__(self, gropath, xtcpath, atoms, blocks, start=0, end=-1, skip=1, **kwargs):
        # Initialize the universe
        self.u = mda.Universe(gropath, xtcpath)
        
        # Set the number of bins
        self.nblocks = blocks

        # Create a list of atoms for each atom
        self.atoms = []
        for atom in atoms:
            self.atoms.append(self.u.select_atoms(f"name {atom}"))

        # Set the start, end, and skip values
        self.start = start
        self.end = end
        self.skip = skip

        # Initialize the AnalysisBase class
        super().__init__(self.u.trajectory, verbose=False, **kwargs)

    def _prepare(self):

        # Initialize the starting position matrix
        self._startpos = []
        self._startsimtime = 0
        for atom in self.atoms:
            self._startpos.append(np.zeros((len(atom), 3)))

        # Calculate the total number of frames and the block size
        self._totalframes = len(self.u.trajectory[self.start:self.end:self.skip])
        self._blocksize = int(self._totalframes/self.nblocks)

        # Initialize the diffusivity matrix
        self.results.msdx = np.zeros((len(self.atoms), self.nblocks, self._blocksize))
        self.results.msdy = np.zeros((len(self.atoms), self.nblocks, self._blocksize))
        self.results.msdz = np.zeros((len(self.atoms), self.nblocks, self._blocksize))

        # initilaze the time list
        self.results.simtime = np.zeros(self._blocksize)

        # start the clock
        self._starttime = time.monotonic()


    def _single_frame(self):

        # Check that the frame index is within one of the blocks
        if int(self._frame_index/self._blocksize) < self.nblocks:
            # Check if the frame index is the start of a block
            if self._frame_index % self._blocksize == 0:      
                # Set the start time for the block    
                self._startsimtime = self._ts.time/1000
            # add the time to the list if it is the first block
            if int(self._frame_index/self._blocksize) == 1:
                self.results.simtime[int(self._frame_index % self._blocksize)] = self._ts.time/1000 - self._startsimtime

            # Calculate the MSD for each atom
            for atom in self.atoms:
                for i in range(len(atom)):
                    # Set the starting position if it is the start of a block
                    if self._frame_index % self._blocksize == 0:
                        self._startpos[self.atoms.index(atom)][i] = atom[i].position
                    
                    # Calculate the MSD
                    self.results.msdx[self.atoms.index(atom), int(self._frame_index/self._blocksize), self._frame_index % self._blocksize] += (atom[i].position[0] - self._startpos[self.atoms.index(atom)][i][0])**2 / len(atom)
                    self.results.msdy[self.atoms.index(atom), int(self._frame_index/self._blocksize), self._frame_index % self._blocksize] += (atom[i].position[1] - self._startpos[self.atoms.index(atom)][i][1])**2 / len(atom)
                    self.results.msdz[self.atoms.index(atom), int(self._frame_index/self._blocksize), self._frame_index % self._blocksize] += (atom[i].position[2] - self._startpos[self.atoms.index(atom)][i][2])**2 / len(atom)
                    #np.linalg.norm(atom[i].position - self._startpos[self.atoms.index(atom)][i])**2 / len(atom)



    def _conclude(self):
        # stop the clock
        self._endtime = time.monotonic()
        self.results.time = self._endtime - self._starttime

        # print(self.results.msd[0][0])
        # print(self.results.simtime)

    

def main():
    args = parse_arguments()

    if not args.simpath.endswith("/"):
        args.simpath += "/"


    gropath = args.simpath + args.confname
    xtcpath = args.simpath + args.trajname
    stdout = args.simpath + args.output + ".txt"
    outpathx = args.simpath + args.output + "x"
    outpathy = args.simpath + args.output + "y"
    outpathz = args.simpath + args.output + "z"

    msd = Diffusivity(gropath, xtcpath, args.atoms, args.blocks, start=args.start, stop=args.end, skip=args.skip)
    msd.run(verbose=True, start=args.start, stop=args.end, step=args.skip)

    # Save the results to a text file
    for atom in range(len(args.atoms)):
        with open(outpathx + f"_{args.atoms[atom]}.xvg", "w") as f:
            for i in range(len(msd.results.simtime)):
                f.write(f"{msd.results.simtime[i]} ")
                for j in range(args.blocks):
                    f.write(f"{msd.results.msdx[atom][j][i]} ")
                f.write("\n")
        with open(outpathy + f"_{args.atoms[atom]}.xvg", "w") as f:
            for i in range(len(msd.results.simtime)):
                f.write(f"{msd.results.simtime[i]} ")
                for j in range(args.blocks):
                    f.write(f"{msd.results.msdy[atom][j][i]} ")
                f.write("\n")
        with open(outpathz + f"_{args.atoms[atom]}.xvg", "w") as f:
            for i in range(len(msd.results.simtime)):
                f.write(f"{msd.results.simtime[i]} ")
                for j in range(args.blocks):
                    f.write(f"{msd.results.msdz[atom][j][i]} ")
                f.write("\n")
    

def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate radial distribution function")
    parser.add_argument("-p", "--simpath", type=str, default="./", help="Path to simulation directory")
    parser.add_argument("-c", "--confname", type=str, default="conf.gro", help="Name of configuration file")
    parser.add_argument("-t", "--trajname", type=str, default="traj.xtc", help="Name of trajectory file")
    parser.add_argument("-o", "--output", type=str, default="msd", help="Output file name")
    parser.add_argument("-a", "--atoms", type=str, nargs="+", help="Name of the atoms ...)")
    parser.add_argument("-b", "--blocks", type=int, default=1, help="Number of blocks for block averaging")
    parser.add_argument("--start", type=int, default=0, help="Starting frame")
    parser.add_argument("--end", type=int, default=-1, help="Ending frame")
    parser.add_argument("--skip", type=int, default=1, help="Frame interval")  

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    main()