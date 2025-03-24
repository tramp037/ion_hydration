import numpy as np 
import pandas as pd
import argparse
import sys
import time
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase

class SolvationOrientation(AnalysisBase):
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
    def __init__(self, gropath, xtcpath, ions, solv, bins, **kwargs):
        # Initialize the universe
        self.u = mda.Universe(gropath, xtcpath)
        
        # Set the number of bins
        self.nbins = bins

        # Set the solvation shell
        self.shell = solv

        # Create a list of atoms for each ion
        self.ions = []
        for ion in ions:
            self.ions.append(self.u.select_atoms(f"name {ion}"))

        # Initialize the AnalysisBase class
        super().__init__(self.u.trajectory, verbose=False, **kwargs)

    def _prepare(self):
        # Create a dictionary of ion solvation shell limits (ion-OW distance in angstroms)
        self._dict = {}
        self._dict['Li'] = [2.63, 4.93]
        self._dict['Na'] = [3.14, 5.43]
        self._dict['K'] = [3.57, 5.77]
        self._dict['Rb'] = [3.58, 5.75]
        self._dict['Cs'] = [3.67, 5.84]
        self._dict['Cl'] = [3.66, 5.97]
        
        # Initialize the orientation matrix
        self.results.orientation = np.zeros((len(self.ions), self.nbins))

        # Create the edges and bins for the histogram
        self._edges = np.linspace(-1, 1, (self.nbins+1))
        self._bins = self._edges[:-1] + np.diff(self._edges)/2
        self._count = [0 for i in self.ions]

        # start the clock
        self._starttime = time.monotonic()


    def _single_frame(self):
        # Get the dimensions of the box
        dimensions = np.array(self.u.dimensions)

        # Loop over each ion
        for ion in self.ions:
            # shells = [0 for i in ion]
            # hydrogens = [0 for i in ion]
            for i in range(len(ion)):
            #     shells[i] = self.u.select_atoms(f"(around {self._dict[ion[i].name]} resid {ion[i].resid}) and (name OW)")
            #     selection_string = ' or '.join(f'resid {atom.resid} and name HW*' for atom in shells[i])
            #     hydrogens[i] = self.u.select_atoms(selection_string)

            # ow_pos = [shell.positions.T for shell in shells]
            # hw_pos = [np.array([hydrogen[2*i:2*i+2].center_of_mass() for i in range(int(len(hydrogen)/2))]).T for hydrogen in hydrogens]
            # ion_ow = [ion[i].position[:, np.newaxis] - ow_pos[i] for i in range(len(ion))]
            # hw_ow = [hw_pos[i] - ow_pos[i] for i in range(len(ion))]

            # for j in range(3):
            #     for i in ion_ow:
            #         i[j] = np.where(i[j] > dimensions[j]/2, i[j] - dimensions[j], i[j])
            #         i[j] = np.where(i[j] < -dimensions[j]/2, i[j] + dimensions[j], i[j])
            #     for i in hw_ow:
            #         i[j] = np.where(i[j] > dimensions[j]/2, i[j] - dimensions[j], i[j])
            #         i[j] = np.where(i[j] < -dimensions[j]/2, i[j] + dimensions[j], i[j])
            
            # angles = [np.arccos(np.sum(ion_ow[i]*hw_ow[i], axis=0)/(np.linalg.norm(ion_ow[i], axis=0)*np.linalg.norm(hw_ow[i], axis=0))) for i in range(len(ion))]
            # for i in angles:
            #     bin_index = np.digitize(i*180/np.pi, self._edges)
            #     if any(bin_index) == self.nbins:
            #         print(i*180/np.pi)
            #     self.results.orientation[self.ions.index(ion), bin_index-1] += 1
            #     self._count += 1
        
                # Select the solvent molecules in the first solvation shell of the ion
                if self.shell == 1:
                    shell = self.u.select_atoms(f"(around {self._dict[ion[i].name][0]} resid {ion[i].resid}) and (name OW)")
                else:
                    shell = self.u.select_atoms(f"((around {self._dict[ion[i].name][self.shell-1]} resid {ion[i].resid}) and (name OW)) and not (around {self._dict[ion[i].name][self.shell-2]} resid {ion[i].resid})")
                # select all full water molecules from shell
                # waters = self.u.atoms[[]]
                # for j in range(len(shell)):
                #     waters = waters + self.u.select_atoms(f"resid {shell[j].resid} and name HW*")
                if len(shell) == 0:
                    print("No solvent molecules found around ion")
                    continue
                for j in range(len(shell)):

                    # Select the hydrogen atoms in the solvent molecule
                    hw = self.u.select_atoms(f"resid {shell[j].resid} and name HW*")
                    if len(hw) == 0:
                        sys.exit("Error, no hydrogen atoms found in solvent molecule. Exiting...")

                    h1_h2 = hw[0].position - hw[1].position
                    for k in range(3):
                        if h1_h2[k] > dimensions[k]/2:
                            h1_h2[k] -= dimensions[k]
                        elif h1_h2[k] < -dimensions[k]/2:
                            h1_h2[k] += dimensions[k]
                    hw_center = hw[1].position + h1_h2 / 2

                    # Calculate the OW-ion and OW-dipole vectors
                    ion_ow = shell[j].position - ion[i].position
                    ow_hw = shell[j].position - hw_center
                    for k in range(3):
                        if ion_ow[k] > dimensions[k]/2:
                            ion_ow[k] -= dimensions[k]
                        elif ion_ow[k] < -dimensions[k]/2:
                            ion_ow[k] += dimensions[k]
                        if ow_hw[k] > dimensions[k]/2:
                            ow_hw[k] -= dimensions[k]
                        elif ow_hw[k] < -dimensions[k]/2:
                            ow_hw[k] += dimensions[k]

                    # Normalize the vectors
                    ion_ow = ion_ow/np.linalg.norm(ion_ow)
                    ow_hw = ow_hw/np.linalg.norm(ow_hw)

                    # Calculate the angle between the vectors
                    dot = np.dot(ion_ow, ow_hw)
                    if dot > 1:
                        dot = 1
                    elif dot < -1:
                        dot = -1
                    # angle = np.arccos(dot)
                        

                    # Bin the angle
                    # bin_index = np.digitize(angle*180/np.pi, self._edges)
                    bin_index = np.digitize((dot), self._edges)
                    if bin_index == self.nbins+1:
                        bin_index -= 1
                    self.results.orientation[self.ions.index(ion), bin_index-1] += 1
                    self._count[self.ions.index(ion)] += 1

    def _conclude(self):
        # Normalize the orientation matrix and save as a DataFrame
        for i in range(len(self.ions)):
            if self._count[i] == 0:
                continue
            self.results.orientation[i] = self.results.orientation[i]/self._count[i]/np.diff(self._edges)
        self.results.orientation = pd.DataFrame(self.results.orientation.T, index=self._bins, columns=[ion.names[0] for ion in self.ions])
        
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

    so = SolvationOrientation(gropath, xtcpath, args.atoms, args.shell, args.bins)
    so.run(verbose=True, start=args.start, stop=args.end, step=args.skip)
    
    # Save the results to a text file
    with open(outpath, "w") as f:
        for i in range(len(so.results.orientation.index)):
            f.write(f"{so.results.orientation.index[i]}")
            for j in range(len(so.results.orientation.columns)):
                f.write(f" {so.results.orientation.iloc[i, j]}")
            f.write("\n")
    

def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate radial distribution function")
    parser.add_argument("-p", "--simpath", type=str, default="./", help="Path to simulation directory")
    parser.add_argument("-c", "--confname", type=str, default="conf.gro", help="Name of configuration file")
    parser.add_argument("-t", "--trajname", type=str, default="traj.xtc", help="Name of trajectory file")
    parser.add_argument("-o", "--output", type=str, default="rdf", help="Output file name")
    parser.add_argument("-a", "--atoms", type=str, nargs="+", help="Name of the ions ...)")
    parser.add_argument("-b", "--bins", type=int, default=36, help="Number of bins")
    parser.add_argument("-s", "--shell", type=int, default=1, help="Solvation shell")
    parser.add_argument("--start", type=int, default=0, help="Starting frame")
    parser.add_argument("--end", type=int, default=-1, help="Ending frame")
    parser.add_argument("--skip", type=int, default=1, help="Frame interval")  

    args = parser.parse_args()
    return args

if __name__ == "__main__":
    main()