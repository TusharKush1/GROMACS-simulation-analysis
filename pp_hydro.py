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

from collections import defaultdict
import mdtraj as md

from collections import defaultdict
import numpy as np

def find_hydrophobic_contacts_fast(traj, group1, group2,
                                   cutoff=0.5, min_occ=0.05):
    # 1) Hydrophobic residue set
    hydro = {'ALA','VAL','LEU','ILE','MET','PHE','PRO','TRP','TYR'}
    natoms = traj.topology.n_atoms

    # 2) Prefilter 0‑based heavy atoms
    a1 = [i for i in group1 if i<natoms
          and traj.topology.atom(i).residue.name in hydro
          and traj.topology.atom(i).element.symbol!='H']
    a2 = [j for j in group2 if j<natoms
          and traj.topology.atom(j).residue.name in hydro
          and traj.topology.atom(j).element.symbol!='H']

    # 3) Build cross‑group pairs ↔ residue keys
    pairs, keymap = [], {}
    for i in a1:
        r1 = traj.topology.atom(i).residue
        for j in a2:
            r2 = traj.topology.atom(j).residue
            if r1.index==r2.index: continue
            k = tuple(sorted((r1, r2), key=lambda R: R.index))
            pairs.append([i, j])
            keymap[(i, j)] = k

    if not pairs:
        return {}

    # 4) Batch compute distances
    d = md.compute_distances(traj, pairs)  # shape=(n_frames,n_pairs)
    hits = (d < cutoff).sum(axis=0)
    occs = hits / traj.n_frames

    # 5) Build result: only keep ≥min_occ
    result = {
        keymap[tuple(pairs[p])]: occs[p] * 100
        for p in range(len(pairs))
        if occs[p] >= min_occ
    }
    return result


# === Call & report ===
hydro = find_hydrophobic_contacts_fast(traj, g1_atoms, g2_atoms,
                                       cutoff=0.5, min_occ=0.05)
txt = f"{args.outprefix}_hydrophobic.txt"
with open(txt, "w") as f:
    f.write("Res1 - Res2 - Occupancy (%)\n")
    for (r1, r2), occ in sorted(hydro.items(), key=lambda x:-x[1]):
        c1 = chr(65 + r1.chain.index)
        c2 = chr(65 + r2.chain.index)
        f.write(f"{r1.name}{r1.resSeq}({c1}) - "
                f"{r2.name}{r2.resSeq}({c2}) - {occ:.1f}\n")
print(f"✅ Saved hydrophobic contacts summary: {txt}")

# === Generate PyMOL script for hydrophobic contacts (same style as hb) ===
pml_file = f"{args.outprefix}_hydrophobic.pml"
hp_residues = set()
with open(pml_file, 'w') as f:
    f.write(f"load {args.top}, complex\n")
    f.write("bg_color white\n")
    f.write("set antialias, 2\n")
    f.write("set ray_trace_mode, 1\n")
    f.write("set ray_trace_gain, 0.1\n")
    f.write("hide everything, all\n")
    f.write("show cartoon, complex\n")
    f.write("set stick_radius, 0.2\n")
    f.write("set label_color, black\n")
    f.write("set label_size, 14\n\n")

    for (r1, r2), occ in sorted(hydro.items(), key=lambda x: -x[1]):
        if occ < 5.0 or r1.chain.index == r2.chain.index:
            continue
        # select and show both residues as sticks
        for res in (r1, r2):
            chain = chr(65 + res.chain.index)
            sel  = f"hp_{res.name}{res.resSeq}_{chain}"
            if sel not in hp_residues:
                f.write(f"select {sel}, chain {chain} and resi {res.resSeq} and resn {res.name}\n")
                f.write(f"show sticks, {sel}\n")
                f.write(f"label {sel} and name CA, \"{occ:.1f}%\"\n")
                hp_residues.add(sel)

    if hp_residues:
        f.write("\nzoom " + " or ".join(hp_residues) + "\n")
    f.write("orient\n")
print(f"✅ Generated PyMOL script: {pml_file}")


