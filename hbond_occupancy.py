#This script allows calculation of hydrogen bonds between protein and ligand from a gromacs simulation after processing the trajectory
#Author : Tushar Kushwaha

import mdtraj as md
from collections import defaultdict
from tqdm import tqdm
import argparse
import pandas as pd

# === Parse arguments ===
parser = argparse.ArgumentParser(description="Protein-ligand residue-wise H-bond occupancy (text + CSV)")
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

n_frames = traj.n_frames
frame_range = range(start_frame, n_frames)

# === Track: protein residue -> frames where it H-bonded to ligand ===
residue_frames = defaultdict(set)

for i in tqdm(frame_range, desc="Processing frames"):
    hbonds = md.baker_hubbard(traj[i], freq=0.1, periodic=False)
    for donor, hydrogen, acceptor in hbonds:
        d_atom = traj.topology.atom(donor)
        a_atom = traj.topology.atom(acceptor)

        if (
            (d_atom.residue.is_protein and a_atom.residue.name == 'UNL') or
            (a_atom.residue.is_protein and d_atom.residue.name == 'UNL')
        ):
            if d_atom.residue.is_protein:
                prot_res = d_atom.residue
            else:
                prot_res = a_atom.residue

            residue_frames[prot_res].add(i)

# === Prepare results ===
results = []
rows = []

for res, frames in residue_frames.items():
    occupancy = 100 * len(frames) / (n_frames - start_frame)
    res_name = res.name.capitalize()
    res_id = res.index
    res_chain = res.chain.index
    results.append( (occupancy, f"{res_name}{res_id} - {occupancy:.1f}%") )
    rows.append([res_chain, res_name, res_id, occupancy])

# === Sort by occupancy descending ===
results.sort(reverse=True)
rows.sort(reverse=True, key=lambda x: x[3])

# === Write plain text output ===
with open('protein_ligand_residuewise_hbonds.txt', 'w') as f:
    for occ, line in results:
        f.write(line + "\n")

print("✅ Saved: protein_ligand_residuewise_hbonds.txt")

# === Write CSV output ===
df = pd.DataFrame(rows, columns=[
    'Protein_Chain', 'Protein_ResName', 'Protein_ResID', 'Occupancy_Percent'
])
df.to_csv('protein_ligand_residuewise_hbonds.csv', index=False)
print("✅ Saved: protein_ligand_residuewise_hbonds.csv")
