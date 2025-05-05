# setup

## This directory contains all the files needed for setting up simulations.

Contents:
------------
`setup.sh`
- This file is the master setup file, which can be edited to perform straightforward MD (other functions to be added later)

`make_$$.sh`
- These files are run when running `setup.sh` depending on the mode in the input. These set up the system for a particular type of simulation.

`gro_io.py`
- Auxiliary functions for modifying .gro configurations
