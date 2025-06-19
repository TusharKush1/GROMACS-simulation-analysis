#This script allows calculation of hydrogen bonds between protein and ligand from a gromacs simulation after processing the trajectory
#Author : Tushar Kushwaha

import mdtraj as md
from collections import Counter
import pandas as pd
from tqdm import tqdm
import argparse

# === Parse arguments ===
parser = argparse.ArgumentParser(description="Protein-ligand H-bond occupancy calculator")
parser.add_argument("--traj", type=str, default="md_0_10_fit.xtc",
                    help="Trajectory file name (default: md_0_10_fit.xtc)")
parser.add_argument("--start", type=int, default=0,
                    help="Starting frame index (default: 0)")
args = parser.parse_args()

traj_file = args.traj
start_frame = args.start

print(f"✅ Using trajectory: {traj_file}")
print(f"✅ Starting from frame: {start_frame}")

# === Load trajectory and topology ===
traj = md.load(traj_file, top='start.gro')

pair_counts = Counter()
n_frames = traj.n_frames

# === Loop over frames with progress bar ===
for i in tqdm(range(start_frame, n_frames), desc="Processing frames"):
    hbonds = md.baker_hubbard(traj[i], freq=0.1, periodic=False)
    for triplet in hbonds:
        donor = triplet[0]
        acceptor = triplet[2]
        d_atom = traj.topology.atom(donor)
        a_atom = traj.topology.atom(acceptor)
        if (
            (d_atom.residue.is_protein and a_atom.residue.name == 'UNL') or
            (a_atom.residue.is_protein and d_atom.residue.name == 'UNL')
        ):
            pair = tuple(sorted([donor, acceptor]))
            pair_counts[pair] += 1

# === Collect results ===
rows = []
for (d, a), count in pair_counts.items():
    d_atom = traj.topology.atom(d)
    a_atom = traj.topology.atom(a)
    if d_atom.residue.is_protein:
        prot, lig = d_atom, a_atom
    else:
        prot, lig = a_atom, d_atom
    rows.append([
        prot.residue.chain.index, prot.residue.name, prot.residue.index, prot.name,
        lig.residue.name, lig.residue.index, lig.name,
        100 * count / (n_frames - start_frame)
    ])

df = pd.DataFrame(rows, columns=[
    'Protein_Chain', 'Protein_ResName', 'Protein_ResID', 'Protein_Atom',
    'Ligand_ResName', 'Ligand_ResID', 'Ligand_Atom',
    'Occupancy_Percent'
])

df.to_csv('protein_ligand_hbonds.csv', index=False)
print('✅ Saved: protein_ligand_hbonds.csv')

