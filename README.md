# GROMACS-simulation-analysis
Shell and python scripts I use for analysis and presentation of GROMACS simulation results 


**hbond_occupancy.py**
This python script calculates hydrogen bonds occupancy from a protein-ligand simulation performed in GROMACS using the processed trajectory file.
Get dependencies : pip3 install mdtraj pandas tqdm

Usage: hbond_occupancy.py [-h] [--traj TRAJ] [--start START]

Protein-ligand H-bond occupancy calculator

optional arguments:
  -h, --help     show this help message and exit
  --traj TRAJ    Trajectory file name (default: md_0_10_fit.xtc)
  --start START  Starting frame index (default: 0)
