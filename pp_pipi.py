#!/usr/bin/env python3
import mdtraj as md
import argparse
import numpy as np
from collections import defaultdict
from tqdm import tqdm

def load_ndx_groups(ndx_file):
    groups = {}
    current_group = None
    with open(ndx_file) as f:
        for line in f:
            if line.startswith('['):
                current_group = line.strip()[1:-1].strip()
                groups[current_group] = []
            elif current_group:
                groups[current_group].extend([int(x) - 1 for x in line.strip().split()])
    return groups

parser = argparse.ArgumentParser(description="Calculate π–π interactions between two groups frame-by-frame")
parser.add_argument("-f", dest="traj", required=True, help="Trajectory file")
parser.add_argument("-s", dest="top", required=True, help="Topology file")
parser.add_argument("-n", dest="ndx", required=True, help="GROMACS index file")
parser.add_argument("-g1", dest="group1", required=True, help="First group name")
parser.add_argument("-g2", dest="group2", required=True, help="Second group name")
parser.add_argument("-o", dest="outprefix", type=str, default="pp", help="Output prefix")
parser.add_argument("--cutoff", type=float, default=0.6, help="Distance cutoff (nm)")
args = parser.parse_args()

traj = md.load(args.traj, top=args.top)
groups = load_ndx_groups(args.ndx)
g1_atoms = groups.get(args.group1, [])
g2_atoms = groups.get(args.group2, [])

ring_atoms_map = {
    'PHE': ['CG','CD1','CD2','CE1','CE2','CZ'],
    'TYR': ['CG','CD1','CD2','CE1','CE2','CZ'],
    'TRP': ['CD2','CE2','CE3','CZ2','CZ3','CH2'],
    'HIS': ['CG','ND1','CD2','CE1','NE2'],
}

def compute_ring_centroids(traj, group):
    arom_res, atom_inds = [], []
    for idx in group:
        atom = traj.topology.atom(idx)
        res = atom.residue
        if res.name in ring_atoms_map:
            inds = [a.index for a in res.atoms if a.name in ring_atoms_map[res.name]]
            if inds and res not in arom_res:
                arom_res.append(res)
                atom_inds.append(inds)
    cents = np.zeros((traj.n_frames, len(arom_res), 3))
    for j, inds in enumerate(atom_inds):
        cents[:, j, :] = traj.xyz[:, inds, :].mean(axis=1)
    return arom_res, cents

# Compute centroids for both groups
res1, cent1 = compute_ring_centroids(traj, g1_atoms)
res2, cent2 = compute_ring_centroids(traj, g2_atoms)

inter_counts = defaultdict(int)
total_frames = traj.n_frames

# Calculate π–π interactions frame by frame
for i in tqdm(range(total_frames), desc="Analyzing frames"):
    for j, r1 in enumerate(res1):
        for k, r2 in enumerate(res2):
            if r1.index == r2.index or r1.chain.index == r2.chain.index:
                continue
            dist = np.linalg.norm(cent1[i, j, :] - cent2[i, k, :])
            if dist < args.cutoff:
                key = (r1, r2)
                inter_counts[key] += 1

# Compute occupancy
out_file = f"{args.outprefix}_pipi.txt"
with open(out_file, 'w') as f:
    f.write("Residue (Chain A) | Residue (Chain B) | Occupancy (%)\n")
    for (r1, r2), count in sorted(inter_counts.items(), key=lambda x: -x[1]):
        occ = (count / total_frames) * 100
        c1 = chr(65 + r1.chain.index)
        c2 = chr(65 + r2.chain.index)
        f.write(f"{r1.name}{r1.resSeq} ({c1}) | {r2.name}{r2.resSeq} ({c2}) | {occ:.1f}\n")

print(f"✅ π–π interaction occupancy saved to: {out_file}")
