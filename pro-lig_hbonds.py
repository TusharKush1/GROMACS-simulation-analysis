"""
Hbond_Analysis-Ligand.py

Description:
------------
This script analyzes a molecular dynamics trajectory to identify and quantify 
residue-wise hydrogen bond occupancy between a protein and a ligand (resname 'UNL').

It performs the following steps:
1. Loads the trajectory (.xtc) and topology (.gro or .pdb).
2. Identifies protein-ligand hydrogen bonds using MDTraj's baker_hubbard algorithm.
3. Calculates the percentage of frames each protein residue forms an H-bond with the ligand.
4. Outputs:
   - A plain text summary of interacting residues and their occupancy.
   - A CSV file with structured data.
   - A PyMOL `.pml` script to visualize the interacting residues and ligand.

The generated PyMOL script:
- Loads the specified PDB structure (typically the first frame).
- Shows ligand and interacting residues as sticks.
- Labels residues with their number and % occupancy.
- Sets a white background and high-quality rendering.
- Zooms to the ligand and nearby hydrogen bonding residues.

Usage:
------
python analyze_and_visualize_hbonds.py \
    --traj md_0_10_fit.xtc \
    --top start.gro \
    --pdb 0.pdb \
    --start 0 \
    --outprefix pro-lig

Dependencies:
-------------
- mdtraj
- pandas
- tqdm

Author: [Your Name]
"""

import mdtraj as md
from collections import defaultdict
from tqdm import tqdm
import argparse
import pandas as pd

# === Parse arguments ===
parser = argparse.ArgumentParser(description="Protein-ligand H-bond occupancy + PyMOL script")
parser.add_argument("--traj", type=str, default="md_0_10_fit.xtc", help="Trajectory file name")
parser.add_argument("--top", type=str, default="start.gro", help="Topology file name")
parser.add_argument("--start", type=int, default=0, help="Starting frame index")
parser.add_argument("--outprefix", type=str, default="pro_lig", help="Output prefix")
parser.add_argument("--pdb", type=str, default="0.pdb", help="Starting structure PDB for PyMOL")
args = parser.parse_args()

traj_file = args.traj
top_file = args.top
start_frame = args.start
outprefix = args.outprefix
pdb_file = args.pdb

# === Load trajectory ===
print(f"✅ Loading trajectory: {traj_file}")
traj = md.load(traj_file, top=top_file)

n_frames = traj.n_frames
frame_range = range(start_frame, n_frames)

# === Compute residue-wise H-bond occupancy ===
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
            prot_res = d_atom.residue if d_atom.residue.is_protein else a_atom.residue
            residue_frames[prot_res].add(i)

# === Prepare output ===
results = []
rows = []

for res, frames in residue_frames.items():
    occupancy = 100 * len(frames) / (n_frames - start_frame)
    res_name = res.name.capitalize()
    res_id = res.resSeq
    res_chain = res.chain.index
    results.append((occupancy, f"{res_name}{res_id} - {occupancy:.1f}%"))
    rows.append([res_chain, res_name, res_id, occupancy])

results.sort(reverse=True)
rows.sort(key=lambda x: -x[3])

txt_file = f"{outprefix}_residuewise_hbonds.txt"
csv_file = f"{outprefix}_residuewise_hbonds.csv"
pml_file = f"{outprefix}_hbonds.pml"

# === Write TXT output ===
with open(txt_file, 'w') as f:
    for _, line in results:
        f.write(line + "\n")
print(f"✅ Saved: {txt_file}")

# === Write CSV output ===
df = pd.DataFrame(rows, columns=[
    'Protein_Chain', 'Protein_ResName', 'Protein_ResID', 'Occupancy_Percent'
])
df.to_csv(csv_file, index=False)
print(f"✅ Saved: {csv_file}")

# === Generate PML ===
hb_residues = []

with open(pml_file, 'w') as f:
    f.write(f"load {pdb_file}, complex\n")
    f.write("bg_color white\n")
    f.write("set antialias, 2\n")
    f.write("set ray_trace_mode, 1\n")
    f.write("set ray_trace_gain, 0.1\n")
    f.write("hide everything, all\n")
    f.write("show cartoon, complex\n")
    f.write("show sticks, resn UNL\n")
    f.write("color red, resn UNL\n")
    f.write("set stick_radius, 0.2\n")
    f.write("set label_color, black\n")
    f.write("set label_size, 14\n")

    for _, line in results:
        try:
            parts = line.split(' - ')
            res = parts[0].strip()
            occ = parts[1].strip()
            resn_letter = res[0].upper()
            resn = ''.join([c for c in res if c.isalpha()])
            resid = ''.join([c for c in res if c.isdigit()])
            sele = f"resn {resn} and resi {resid}"
            label_text = f"{resn_letter}{resid} - {occ}"
            f.write(f"select hb_{resn}{resid}, {sele}\n")
            f.write(f"show sticks, hb_{resn}{resid}\n")
            f.write(f"label hb_{resn}{resid} and name CA, \"{label_text}\"\n")
            hb_residues.append(f"hb_{resn}{resid}")
        except Exception as e:
            print(f"❌ Could not parse line: {line.strip()} ({e})")

    if hb_residues:
        hb_sele = " or ".join(hb_residues)
        f.write(f"zoom resn UNL or ({hb_sele})\n")
    else:
        f.write("zoom resn UNL\n")
    f.write("orient\n")

print(f"✅ Generated PyMOL script: {pml_file}")

