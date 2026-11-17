import mdtraj as md
import argparse
from collections import defaultdict
from tqdm import tqdm
import numpy as np


def load_ndx_groups(ndx_file):
    with open(ndx_file) as f:
        lines = f.readlines()

    groups = {}
    current_group = None
    for line in lines:
        if line.startswith('['):
            current_group = line.strip()[1:-1].strip()
            groups[current_group] = []
        else:
            if current_group:
                groups[current_group].extend([int(x)-1 for x in line.strip().split()])
    return groups

# === Parse arguments ===
parser = argparse.ArgumentParser(description="Print atom counts for selected groups")
parser.add_argument("-f", dest="traj", required=True, help="Trajectory file")
parser.add_argument("-s", dest="top", required=True, help="Topology file")
parser.add_argument("-n", dest="ndx", required=True, help="GROMACS index file")
parser.add_argument("-g1", dest="group1", required=True, help="First group name")
parser.add_argument("-g2", dest="group2", required=True, help="Second group name")
parser.add_argument("-o", dest="outprefix", type=str, default="pp", help="Output prefix")
args = parser.parse_args()

# === Load trajectory and index groups ===
traj = md.load(args.traj, top=args.top)
groups = load_ndx_groups(args.ndx)



g1_atoms = groups.get(args.group1, [])
g2_atoms = groups.get(args.group2, [])

print(f"Atoms in {args.group1}: {len(g1_atoms)}")
print(f"Atoms in {args.group2}: {len(g2_atoms)}")

import numpy as np
from collections import defaultdict

# 1) Map of ring atoms for each aromatic residue
ring_atoms_map = {
    'PHE': ['CG','CD1','CD2','CE1','CE2','CZ'],
    'TYR': ['CG','CD1','CD2','CE1','CE2','CZ'],
    'TRP': ['CD2','CE2','CE3','CZ2','CZ3','CH2'],
    'HIS': ['CG','ND1','CD2','CE1','NE2'],
}

# 2) Atom names to treat as “cation” 
cation_names = {'NZ',         # Lysine
               'NH1','NH2',  # Arginine
               'NE'}         # Arginine or protonated His

def compute_ring_centroids(traj, group):
    """Return (residues, centroids) for aromatic rings in `group`."""
    arom_res, atom_inds = [], []
    for idx in group:
        atom = traj.topology.atom(idx)
        res  = atom.residue
        if res.name in ring_atoms_map:
            inds = [a.index for a in res.atoms if a.name in ring_atoms_map[res.name]]
            if inds:
                arom_res.append(res)
                atom_inds.append(inds)
    n_frames = traj.n_frames
    cents    = np.zeros((n_frames, len(arom_res), 3))
    for j, inds in enumerate(atom_inds):
        # mean over ring atoms
        cents[:, j, :] = traj.xyz[:, inds, :].mean(axis=1)
    return arom_res, cents

def find_pi_cation(traj, group1, group2,
                   cutoff=0.6, min_occ=0.05):
    """
    Identify cation–π contacts between group1 and group2.
    Returns dict {(res_ring, res_cat): occupancy_percent}.
    """
    # 3) compute ring centroids for group1
    ring_res, centroids = compute_ring_centroids(traj, group1)
    # 4) collect cationic atoms in group2
    cation_list = []
    for idx in group2:
        atom = traj.topology.atom(idx)
        if atom.name in cation_names:
            cation_list.append((atom.residue, idx))

    result = {}
    # 5) for each ring and each cation atom, compute occupancy
    for j, res_r in enumerate(ring_res):
        cent_vals = centroids[:, j, :]  # shape (n_frames,3)
        for res_c, cat_idx in cation_list:
            # skip intra‑chain or same residue
            if res_r.chain.index == res_c.chain.index or res_r.index == res_c.index:
                continue
            cat_vals = traj.xyz[:, cat_idx, :]  # (n_frames,3)
            dists    = np.linalg.norm(cent_vals - cat_vals, axis=1)
            occ      = np.mean(dists < cutoff)
            if occ >= min_occ:
                key = (res_r, res_c) if res_r.index < res_c.index else (res_c, res_r)
                # if multiple atoms map to same residue pair, take max occupancy
                result[key] = max(result.get(key, 0), occ * 100)

    return result

# ——— After loading traj, groups, etc. ————————————————
pi_cation = find_pi_cation(traj, g1_atoms, g2_atoms,
                           cutoff=0.6, min_occ=0.05)

# ——— Write summary (.txt) —————————————————————————
txt = f"{args.outprefix}_pi_cation.txt"
with open(txt, 'w') as f:
    f.write("Aromatic_Res - Cationic_Res - Occupancy (%)\n")
    for (r1, r2), occ in sorted(pi_cation.items(), key=lambda x: -x[1]):
        c1 = chr(65 + r1.chain.index)
        c2 = chr(65 + r2.chain.index)
        f.write(f"{r1.name}{r1.resSeq}({c1}) - "
                f"{r2.name}{r2.resSeq}({c2}) - {occ:.1f}\n")
print(f"✅ Saved π–cation summary: {txt}")


# ——— Generate PyMOL script (.pml) —————————————————————
pml_file = f"{args.outprefix}_pi_cation.pml"
sel_set = set()

with open(pml_file, 'w') as f:
    f.write(f"load {args.top}, complex\n")
    f.write("bg_color white\n")
    f.write("set antialias, 2\n")
    f.write("set ray_trace_mode, 1\n")
    f.write("set ray_trace_gain, 0.1\n")
    f.write("hide everything, all\n")
    f.write("show cartoon, complex\n")
    f.write("set dash_color, orange\n")
    f.write("set dash_width, 2\n")
    f.write("set dash_gap, 0.3\n")
    f.write("set label_color, black\n")
    f.write("set label_size, 14\n\n")

    for idx, ((r1, r2), occ) in enumerate(sorted(pi_cation.items(), key=lambda x: -x[1])):
        if occ < 5.0 or r1.chain.index == r2.chain.index:
            continue
        c1 = chr(65 + r1.chain.index)
        c2 = chr(65 + r2.chain.index)
        sel1 = f"pc_{idx}_a"
        sel2 = f"pc_{idx}_b"

        # highlight both residues
        for sel, res, chain in ((sel1, r1, c1), (sel2, r2, c2)):
            if sel not in sel_set:
                f.write(f"select {sel}, chain {chain} and resi {res.resSeq} and resn {res.name}\n")
                f.write(f"show sticks, {sel}\n")
                sel_set.add(sel)

        # draw and label dashed line
        line = f"pi_cat_{idx}"
        f.write(f"distance {line}, {sel1}/CG, {sel2}/NZ\n")
        f.write(f"label {line}, \"{occ:.1f}%\"\n\n")

    if sel_set:
        f.write("zoom " + " or ".join(sel_set) + "\n")
    f.write("orient\n")

print(f"✅ Generated PyMOL script: {pml_file}")

