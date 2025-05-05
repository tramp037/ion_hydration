import numpy as np 
import pandas as pd
import argparse
import sys
import time
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase

class IonIonDistance(AnalysisBase):
    """
    Calculate the orientation of solvent molecules in the first solvation shell of ions

    Parameters
    ----------
    gropath : str
        Path to GRO file
    xtcpath : str
        Path to XTC file
    ions : list
        List of different ions in the system
    bins : int
        Number of bins to divide the angle range into

    """
    def __init__(self, gropath, xtcpath, ions, **kwargs):
        # Initialize the universe
        self.u = mda.Universe(gropath, xtcpath)

        # Create a list of atoms for each ion
        self.ions = []
        for ion in ions:
            self.ions.append(self.u.select_atoms(f"name {ion}"))

        # Initialize the AnalysisBase class
        super().__init__(self.u.trajectory, verbose=False, **kwargs)

    def _prepare(self):
        
        # Initialize the distance array
        self.results.distance = np.zeros(len(self.u.trajectory))

        # Initialize the time array
        self.results.time = np.zeros(len(self.u.trajectory))

        # start the clock
        self._starttime = time.monotonic()


    def _single_frame(self):
        # Get the dimensions of the box
        self.dimensions = self.u.dimensions
        self.results.time[self.u.trajectory.frame] = self.u.trajectory.time

        # get the position of each ion
        ion_pos = [ion.positions for ion in self.ions]
        # adjust for pbc
        for j in range(3):
            if ion_pos[0][0][j] - ion_pos[1][0][j] > self.dimensions[j] / 2:
                ion_pos[1][0][j] += self.dimensions[j]
            elif ion_pos[0][0][j] - ion_pos[1][0][j] < -self.dimensions[j] / 2:
                ion_pos[1][0][j] -= self.dimensions[j]
        distance = np.linalg.norm(ion_pos[0] - ion_pos[1])
        self.results.distance[self.u.trajectory.frame] = distance

    def _conclude(self):
        
        # stop the clock
        self._endtime = time.monotonic()
        self.results.simtime = self._endtime - self._starttime

    

def main():
    args = parse_arguments()

    if not args.simpath.endswith("/"):
        args.simpath += "/"

    gropath = args.simpath + args.confname
    xtcpath = args.simpath + args.trajname
    stdout = args.simpath + args.output + ".txt"
    outpath = args.simpath + args.output + f".xvg"

    iid = IonIonDistance(gropath, xtcpath, args.atoms)
    iid.run(verbose=True, start=args.start, stop=args.end, step=args.skip)
    with open(outpath, "w") as f:
        f.write(f"# Time (ps) Distance (nm)\n")
        for i in range(1,len(iid.results.time)-1):
            f.write(f"{iid.results.time[i]} {iid.results.distance[i]}\n")


    
    # # Save the results to a text file
    # for ion in range(len(args.atoms)):
    #     with open(outpath[ion], "w") as f:
    #         f.write(f"{args.bins[0]} {args.bins[1]}\n")
    #         for i in range(len(sf.bins[1])):
    #             f.write(f"{sf.bins[1][i]} ")
    #         f.write("\n")
    #         for j in range(len(sf.bins[0])):
    #             f.write(f"{sf.bins[0][j]} ")
    #             for i in range(len(sf.bins[1])):
    #                 f.write(f"{sf.results.orientation[ion, j, i]} ")
    #             f.write("\n")
    

def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate radial distribution function")
    parser.add_argument("-p", "--simpath", type=str, default="./", help="Path to simulation directory")
    parser.add_argument("-c", "--confname", type=str, default="conf.gro", help="Name of configuration file")
    parser.add_argument("-t", "--trajname", type=str, default="traj.xtc", help="Name of trajectory file")
    parser.add_argument("-o", "--output", type=str, default="rdf", help="Output file name")
    parser.add_argument("-a", "--atoms", type=str, nargs="+", help="Name of the ions ...)")
    parser.add_argument("--start", type=int, default=0, help="Starting frame")
    parser.add_argument("--end", type=int, default=-1, help="Ending frame")
    parser.add_argument("--skip", type=int, default=1, help="Frame interval")  

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    main()
