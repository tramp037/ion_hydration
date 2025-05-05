import pandas as pd
import random
import sys

def readgro(input_file,footer=1,verbose=True):
    """
    Reads a gro file and returns the dimensions of the box and the gro file entries in a pandas DataFrame format

    Parameters
    ----------
    input_file : str
        Filename of the input gro file
    footer : int, optional
        Set the number of lines to skip from the footer, by default 1 (for unit cell dimensions)
    verbose : bool, optional
        Verbose output, by default True

    Returns
    -------
    dim : list
        List of the dimensions of the box
    df : pandas.DataFrame
        gro file entries in a pandas DataFrame format
    
    """    
    if verbose:
        print("\n****************************************************************************")
        print("*                                                                          *")
        print("*                             gro_io.readgro()                             *")
        print("*               Read in atoms and coordinates of a .gro file               *")
        print("*                                                                          *")
        print("*                         Written By: Naomi Trampe                         *")
        print("*                                SAMPEL Lab                                *")
        print("*                         Last updated: 09/19/2023                         *")
        print("*                                                                          *")
        print("****************************************************************************\n\n") 
        print("Reading gro file: ",input_file)
    
    # read in the gro file
    df = pd.read_fwf(input_file, skiprows=2, skipfooter=footer, header=None, widths=[5, 4, 6, 5, 8, 8, 8], names=["Resid", "Resname", "Atom", "Index", "x", "y", "z"])
    # add an index column
    for i in range(len(df)):
        df.loc[i,'Index'] = i+1
    
    # read in the dimensions of the box
    with open(input_file, "r") as file:
        dim = file.readlines()[-1].split()
        dim = [float(dim[0]),float(dim[1]),float(dim[2])]

    # return the dimensions and the gro file entries
    if verbose:
        print("Finished Reading gro file: ",input_file)
    return dim, df

def writegro(title,df,dim,output_file,verbose=True):
    """
    Writes a gro file from a pandas DataFrame
    
    Parameters
    ----------
    title : str
        header of the gro file
    df : pandas.DataFrame
        gro file coordinates in a pandas DataFrame format
    dim : list
        box dimensions for the gro file
    output_file : str
        name of the output gro file
    verbose : bool, optional
        Verbose output, by default True
    """    
    if verbose:
        print("\n****************************************************************************")
        print("*                                                                          *")
        print("*                            gro_io.writegro()                             *")
        print("*              Write out atoms and coordinates to a .gro file              *")
        print("*                                                                          *")
        print("*                         Written By: Naomi Trampe                         *")
        print("*                                SAMPEL Lab                                *")
        print("*                         Last updated: 09/19/2023                         *")
        print("*                                                                          *")
        print("****************************************************************************\n\n") 
        print("Writing gro file: ",output_file)
    
    # write out the gro file
    n = len(df)
    # current function does not work with > 99999 atoms, the index line will need to be modified to start again at 00001 
    if n > 99999:
        print("gro file has > 99999 atoms. The current function does not support this number of atoms.")
        sys.exit(1)
    else:
        with open(output_file,'w') as f:
            f.write(title+"\n")
            f.write(f'{n}'.rjust(5)+'\n')
            resid = 0
            for i in range(n):
                # check if the residue number has changed
                if resid == 0 or df.loc[i,'Resid'] != df.loc[i-1,'Resid']:
                    resid += 1
                # write out the residue number, residue name, atom name, atom number, and coordinates
                f.write(str(resid).rjust(5))
                f.write(df.loc[i,'Resname'].ljust(4))
                f.write(str(df.loc[i,'Atom']).rjust(6))
                f.write(str(i+1).rjust(5))
                f.write(f'{float(df.loc[i,"x"]):.3f}'.rjust(8))
                f.write(f'{float(df.loc[i,"y"]):.3f}'.rjust(8))
                f.write(f'{float(df.loc[i,"z"]):.3f}'.rjust(8))
                f.write('\n')
            # write out the dimensions of the box
            f.write(f"{float(dim[0]):.5f}".rjust(10)+f"{float(dim[1]):.5f}".rjust(10)+f"{float(dim[2]):.5f}".rjust(10)+"\n")

        if verbose:
            print("Finished Writing gro file: ",output_file,"\n")
        
def combine(file1,file2,outfile,dim,verbose=True):
    """
    Combines the two gro files into one gro file

    Parameters
    ----------
    file1 : str
        first gro file
    file2 : str
        second gro file
    outfile : str
        name of the output gro file
    a : float
        Box dimension in x  
    b : float
        Box dimension in y   
    c : float
        Box dimension in z
    verbose : bool, optional   
        Verbose output, by default True 
    """    
    if verbose:
        print("\n****************************************************************************")
        print("*                                                                          *")
        print("*                             gro_io.combine()                             *")
        print("*                        Concatenate two .gro files                        *")
        print("*                                                                          *")
        print("*                         Written By: Naomi Trampe                         *")
        print("*                                SAMPEL Lab                                *")
        print("*                         Last updated: 09/19/2023                         *")
        print("*                                                                          *")
        print("****************************************************************************\n\n") 
        print("Reading gro file 1:",file1)
    
    # read in the gro files
    d, df1 = readgro(file1,verbose=False)
    if verbose:
        print("Reading gro file 2:",file2)
    d, df2 = readgro(file2,verbose=False)

    # combine the gro files
    if verbose:
        print("Combining files...")
    frames = [df1,df2]
    gro = pd.concat(frames)
    n = len(gro)

    # renumber the atoms
    for i in range(n):
        gro.iloc[i,3]=i+1
        if verbose and (i+1) % int((n)/10) - (n % 10) == 0:
            print(f"{i+1}/{n}")
    gro = gro.reset_index(drop=True)

    # write out the gro file
    if verbose:
        print("Writing gro file: ",outfile)
    writegro("Membrane",gro,dim,outfile,verbose=False)
    if verbose:
        print("Finished Combining files:", file1, "and", file2, "\n")
    
def renumber(input_file,output_file,resname_order,verbose=True):
    """
    Renumber and sort the gro file entries

    Parameters
    ----------
    input_file : str
        gro file to be renumbered
    output_file : str
        name of the output gro file
    order : list
        order of the atoms
    verbose : bool, optional
        Verbose output, by default True
    """    
    if verbose:
        print("\n****************************************************************************")
        print("*                                                                          *")
        print("*                            gro_io.renumber()                             *")
        print("*                      Renumber atoms of a .gro file                       *")
        print("*                                                                          *")
        print("*                         Written By: Naomi Trampe                         *")
        print("*                                SAMPEL Lab                                *")
        print("*                         Last updated: 02/21/2025                         *")
        print("*                                                                          *")
        print("****************************************************************************\n\n") 
        print("Reading gro file: ",input_file)
    
    # read in the gro file
    dim, df = readgro(input_file,verbose=False)
    
    # renumber atoms
    if verbose:
        print("Renumbering atoms...")
    df['Resname'] = pd.Categorical(df['Resname'],categories=resname_order,ordered=True)
    df = df.sort_values(by=['Resname','Index']).reset_index(drop=True)
    
    # write out the gro file
    if verbose:
        print("Writing gro file: ",output_file)
    writegro("Ions in water",df,dim,output_file,verbose=False)
    if verbose:
        print("Finished Renumbering file:", input_file, "\n")

def write_indices(counts,label,f,verbose=True):
    if len(counts[label]) > 0:
            if verbose:
                print("Writing [ " + label + " ]")
            f.write("\n[ " + label + " ]\n")
            for cnt, i in enumerate(counts[label]):
                f.write(str(i).rjust(4))
                if cnt % 15 == 14 and cnt != len(counts[label])-1:
                    f.write("\n")
                else:
                    f.write(" ")

def index(input_file,output_file,verbose=True):
    """
    Write out the index file for the system

    Parameters
    ----------
    input_file : str
        gro file to be indexed
    output_file : str
        name of the output index file
    verbose : bool, optional
        Verbose output, by default True
    """

    if verbose:
        print("\n****************************************************************************")
        print("*                                                                          *")
        print("*                              gro_io.index()                              *")
        print("*                           Write an index file                            *")
        print("*                                                                          *")
        print("*                         Written By: Naomi Trampe                         *")
        print("*                                SAMPEL Lab                                *")
        print("*                         Last updated: 02/21/2025                         *")
        print("*                                                                          *")
        print("****************************************************************************\n\n") 
        print("Reading gro file: ",input_file)

    dim, df = readgro(input_file,verbose=False)
    
    # initialize the counts
    counts = {}
    # count for atoms in the solution
    counts["Solution"] = []
    # count for ions in the solution
    counts["Ions"] = []
    # count for solvent in the solution
    counts["Solvent"] = []
    counts["Pull"] = []

    if verbose:
        print("Writing index file: ",output_file)
    with open(output_file, 'w') as f:
        # write indices of the system
        if verbose:
            print("Writing [ System ]")
        f.write("[ System ]\n")
        cnt = 0
        for index, row in df.iterrows():
            f.write(str(row["Index"]).rjust(4))
            if cnt % 15 == 14 and cnt != len(df)-1:
                f.write("\n")
            else:
                f.write(" ")
            cnt += 1
           
            counts["Solution"].append(row["Index"])
            # count the number of ions and solvent
            if "SOL" not in row["Resname"] and "ETO" not in row["Resname"]:
                counts["Ions"].append(row["Index"])
            else:
                counts["Solvent"].append(row["Index"])

        # write indices of the solution
        write_indices(counts,"Solution",f,verbose=verbose)
        
        # write indices of the ions
        write_indices(counts,"Ions",f,verbose=verbose)
        
        # write indices of the solvent
        write_indices(counts,"Solvent",f,verbose=verbose)
        
        # # write indices of the pulled ion
        # if len(counts["Ions"]) > 0:
        #     if verbose:
        #         print("Writing [ Pull ]")
        #     f.write("\n[ Pull ]\n")
        #     f.write(str(counts["Ions"][0]).rjust(4))
        # f.write("\n")

def ion_box(outpath,dim,cation,anion,cat,an,test_n,buff=[0,0,0,0,0,0],loc=[0,0,0],verbose=True):
    
    """
    Creates a box of ions with the specified dimensions and number of each ion

    Parameters
    ----------
    outpath : str
        path to the output directory - output will be named outpath+'ions.gro'
    dim : list
        dimensions of the box
    cation : list
        list of cations
    anion : list
        list of anions
    cat : list
        list of the number of each cation
    an : list
        list of the number of each anion
    test_n : int
        test number, to be used as the random seed
    buff : list, optional
        buffer around the edges, by default [0,0,0,0,0,0] (-x,+x,-y,+y,-z,+z)
    loc : list, optional
        location of the first cation, by default [0,0,0]
    verbose : bool, optional
        Verbose output, by default True
    """
    
    if verbose:
        print("\n****************************************************************************")
        print("*                                                                          *")
        print("*                             gro_io.ion_box()                             *")
        print("*                           Create a box of ions                           *")
        print("*                                                                          *")
        print("*                         Written By: Naomi Trampe                         *")
        print("*                                SAMPEL Lab                                *")
        print("*                         Last updated: 10/26/2023                         *")
        print("*                                                                          *")
        print("****************************************************************************\n\n") 
        print("Creating ion box...")

    if test_n >= 0:
        random.seed(test_n)

    df = pd.DataFrame(columns=["Resid", "Resname", "Atom", "Index", "x", "y", "z"])

    # add all of the ions
    sum_ions = sum(cat) + sum(an)
 
    num = 0
    ion_cumulative = []
    for i in range(len(cat)):
        num += cat[i]
        ion_cumulative.append(num)
    for i in range(len(an)):
        num += an[i]
        ion_cumulative.append(num)
    all_ions = cation + anion
    n_all_ions = cat + an
    for i in range(sum_ions):
        # get the index of the first value of ion_cumulative that is greater than i
        for j in range(len(ion_cumulative)):
            if ion_cumulative[j] > i:
                idx = j
                break
        # add the ion
        ion = all_ions[idx]
        if verbose and i in ion_cumulative:
            print(f"Adding {n_all_ions[idx]} {ion}")
        too_close = True
        if loc != [0,0,0] and i == 0:
            x = loc[0]
            y = loc[1]
            z = loc[2]
            too_close = False
        cnt = 0
        while too_close == True:
            too_close = False
            x = round(random.uniform(buff[0], dim[0]-buff[1]), 3)
            y = round(random.uniform(buff[2], dim[1]-buff[3]), 3)
            z = round(random.uniform(buff[4], dim[2]-buff[5]), 3)
            for index, row in df.iterrows():
                if (row["x"] - x)**2 + (row["y"] - y)**2 + (row["z"] - z)**2 < 0.1**2:
                    too_close = True
                    break
            cnt += 1
            if cnt > 100:
                print("Error: too many iterations")
                break
        df.loc[len(df.index)] = [i + 1, ion, ion, len(df) + 1, x, y, z]
    if verbose:
        print("Writing gro file: ",outpath+'ions.gro')
    writegro("Ion Box",df,dim,outpath+'ions.gro',verbose=False)

def wrap(input_file,output_file,verbose=True):
    """
    Wraps the coordinates of the input gro file

    Parameters
    ----------
    input_file : str
        gro file to be wrapped
    output_file : str
        name of the output gro file
    verbose : bool, optional
        Verbose output, by default True
    """

    if verbose:
        print("\n****************************************************************************")
        print("*                                                                          *")
        print("*                              gro_io.wrap()                               *")
        print("*                      Wrap molecules outside the box                      *")
        print("*                                                                          *")
        print("*                         Written By: Naomi Trampe                         *")
        print("*                                SAMPEL Lab                                *")
        print("*                         Last updated: 09/19/2023                         *")
        print("*                                                                          *")
        print("****************************************************************************\n\n") 
        print("Reading gro file: ",input_file)
    
    dim, df = readgro(input_file,verbose=False)

    resids = []
    insides = [[],[],[]]
    if verbose:
        print("Determining molecules to wrap...")
    for index, row in df.iterrows():
        if row["Resid"] not in resids:
            resids.append(row["Resid"])
            insides[0].append(False)
            insides[1].append(False)
            insides[2].append(False)
        x = float(row['x'])
        y = float(row['y'])
        z = float(row['z'])
        if x <= float(dim[0]):
            insides[0][resids.index(row["Resid"])] = True
        if y <= float(dim[1]):
            insides[1][resids.index(row["Resid"])] = True
        if z <= float(dim[2]):
            insides[2][resids.index(row["Resid"])] = True
    if verbose:
        print("Wrapping molecules...")
    for index, row in df.iterrows():
        if not insides[0][resids.index(row["Resid"])]:
            df.loc[index,'x'] = row['x'] - float(dim[0])
        if not insides[1][resids.index(row["Resid"])]:
            df.loc[index,'y'] = row['y'] - float(dim[1])
        if not insides[2][resids.index(row["Resid"])]:
            df.loc[index,'z'] = row['z'] - float(dim[2])

    writegro("Wrapped",df,dim,output_file,verbose=False)


def cut(input_file, low_dim, up_dim, output_file,verbose=True):
    """
    cut a gro file down to the specified dimensions

    Parameters
    ----------
    input_file : str
        gro file to be cut
    low_dim : list
        lower dimensions of the box
    up_dim : list
        upper dimensions of the box
    output_file : str
        name of the output gro file
    verbose : bool, optional
        Verbose output, by default True
    """

    if verbose:
        print("\n****************************************************************************")
        print("*                                                                          *")
        print("*                               gro_io.cut()                               *")
        print("*                 Cut the box to the specified dimensions                  *")
        print("*                                                                          *")
        print("*                         Written By: Naomi Trampe                         *")
        print("*                                SAMPEL Lab                                *")
        print("*                         Last updated: 07/28/2024                         *")
        print("*                                                                          *")
        print("****************************************************************************\n\n") 
        print("Reading gro file: ",input_file)
    
    dim, df = readgro(input_file,verbose=False)

    resids = []
    keep = []
    if verbose:
        print("Determining molecules to cut...")
    for index, row in df.iterrows():
        if row["Resid"] not in resids:
            resids.append(row["Resid"])
            keep.append(True)
        x = float(row['x'])
        y = float(row['y'])
        z = float(row['z'])
        if x > float(up_dim[0]) or y > float(up_dim[1]) or z > float(up_dim[2]) or x < float(low_dim[0]) or y < float(low_dim[1]) or z < float(low_dim[2]):
            keep[resids.index(row["Resid"])] = False
    if verbose:
        print("Cutting box...")
    for index, row in df.iterrows():
        if not keep[resids.index(row["Resid"])]:
            df.drop(index, inplace=True)

    df = df.reset_index(drop=True)
            

    writegro("Cut",df,[up_dim[i]-low_dim[i] for i in range(3)],output_file,verbose=False)

def rotate(input_file, orient, output_file,verbose=True):
    """
    Rotate the gro file to the specified orientation

    Parameters
    ----------
    input_file : str
        gro file to be cut
    orient : list
        orientation of the box
    output_file : str
        name of the output gro file
    verbose : bool, optional
        Verbose output, by default True
    """

    if verbose:
        print("\n****************************************************************************")
        print("*                                                                          *")
        print("*                              gro_io.rotate()                             *")
        print("*                            Rotates a gro file                            *")
        print("*                                                                          *")
        print("*                         Written By: Naomi Trampe                         *")
        print("*                                SAMPEL Lab                                *")
        print("*                         Last updated: 07/29/2024                         *")
        print("*                                                                          *")
        print("****************************************************************************\n\n") 
        print("Reading gro file: ",input_file)
    
    dim, df = readgro(input_file,verbose=False)

    if verbose:
        print("Rotating box...")
    for index, row in df.iterrows():
        coords = [float(row['x']),float(row['y']),float(row['z'])]
        df.loc[index,'x'] = coords[orient[0]]
        df.loc[index,'y'] = coords[orient[1]]
        df.loc[index,'z'] = coords[orient[2]]

    df = df.reset_index(drop=True)
            
    writegro("Rotate",df,[dim[orient[0]], dim[orient[1]], dim[orient[2]]],output_file,verbose=False)