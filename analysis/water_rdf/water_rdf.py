import numpy as np 
import pandas as pd
import argparse
import sys
import time
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase

class WaterRDF(AnalysisBase):
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
    def __init__(self, gropath, xtcpath, atoms, limit, bins, **kwargs):
        # Initialize the universe
        self.u = mda.Universe(gropath, xtcpath)
        
        # Set the number of bins
        self.nbins = bins
        # Set the limit of the RDF
        self.limit = limit

        # Create a list of atoms for each ion
        self.atoms = []
        for a in range(len(atoms)):
            if a == 0:
                self.atom = atoms[a]
                self.n_atoma = self.u.select_atoms(f"name {atoms[a]}*").n_atoms
            else:  
                self.atoms.append(self.u.select_atoms(f"name {atoms[a]}*"))

        # Initialize the AnalysisBase class
        super().__init__(self.u.trajectory, verbose=False, **kwargs)

    def _prepare(self):
        
        # Initialize the orientation matrix
        self.results.rdf = np.zeros((len(self.atoms), self.nbins))

        # Create the edges and bins for the histogram
        self._edges = np.linspace(0, self.limit, (self.nbins+1))
        self._bins = self._edges[:-1] + np.diff(self._edges)/2
        self._count = [0 for i in self.atoms]

        # Function to calculate the volume of a shell
        self.vol = lambda r1, r2: (4 / 3) * np.pi * (r2 ** 3 - r1 ** 3)

        # start the clock
        self._starttime = time.monotonic()


    def _single_frame(self):
        # Get the dimensions of the box
        dimensions = np.array(self.u.dimensions)

        # Loop over each ion
        for atom in self.atoms:
            self._count[self.atoms.index(atom)] += 1
            for i in range(len(atom)):
                shell = self.u.select_atoms(f"(around {self.limit} resid {atom[i].resid}) and (name {self.atom}*) and not (resid {atom[i].resid})")
                if len(shell) == 0:
                    print("No atoms found around atom")
                    continue

                shell_pos = shell.positions

                dx = np.min([np.abs(atom[i].position[0] - shell_pos[:,0]), np.abs(dimensions[0] - np.abs(atom[i].position[0] - shell_pos[:,0]))], axis=0)
                dy = np.min([np.abs(atom[i].position[1] - shell_pos[:,1]), np.abs(dimensions[1] - np.abs(atom[i].position[1] - shell_pos[:,1]))], axis=0)
                dz = np.min([np.abs(atom[i].position[2] - shell_pos[:,2]), np.abs(dimensions[2] - np.abs(atom[i].position[2] - shell_pos[:,2]))], axis=0)
                r = np.sqrt(dx**2 + dy**2 + dz**2)

                idx = self._edges.searchsorted(r, "right")-1
                for j in idx:
                    if j < self.nbins:
                        self.results.rdf[self.atoms.index(atom), j] += 1 * dimensions[0] * dimensions[1] * dimensions[2] / len(atom)

                # for j in range(len(shell)):
                #     # Calculate the distance between the ion and the shell atom
                #     vec = shell[j].position - atom[i].position
                #     for k in range(3):
                #         if vec[k] > dimensions[k]/2:
                #             vec[k] -= dimensions[k]
                #         elif vec[k] < -dimensions[k]/2:
                #             vec[k] += dimensions[k]
                #     r = np.linalg.norm(vec)
                #     print(r)
                #     sys.exit(0)


                #     # Bin the distance
                #     bin_index = np.digitize(r, self._edges)
                #     if bin_index == self.nbins+1:
                #         continue
                #     self.results.rdf[self.atoms.index(atom), bin_index-1] += 1 * dimensions[0] * dimensions[1] * dimensions[2] / len(atom)

    def _conclude(self):
        # Normalize the orientation matrix and save as a DataFrame
        for i in range(len(self.atoms)):
            if self._count[i] == 0:
                continue
            for b in range(self.nbins):
                self.results.rdf[i, b] = self.results.rdf[i, b]/self._count[i]/self.vol(self._edges[b], self._edges[b+1])/self.n_atoma
        self.results.rdf = pd.DataFrame(self.results.rdf.T, index=self._bins, columns=[atom.names[0] for atom in self.atoms])
        
        # stop the clock
        self._endtime = time.monotonic()
        self.results.time = self._endtime - self._starttime

    

def main():
    args = parse_arguments()

    if not args.simpath.endswith("/"):
        args.simpath += "/"

    gropath = args.simpath + args.confname
    xtcpath = args.simpath + args.trajname
    stdout = args.simpath + args.output + ".txt"
    outpath = args.simpath + args.output + ".xvg"

    wr = WaterRDF(gropath, xtcpath, args.atoms, args.limit, args.bins)
    wr.run(verbose=True, start=args.start, stop=args.end, step=args.skip)
    
    # # Save the results to a text file
    with open(outpath, "w") as f:
        for i in range(len(wr.results.rdf.index)):
            f.write(f"{wr.results.rdf.index[i]}")
            for j in range(len(wr.results.rdf.columns)):
                f.write(f" {wr.results.rdf.iloc[i, j]}")
            f.write("\n")
    

def parse_arguments():
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