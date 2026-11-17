import mdtraj as md
from collections import defaultdict
import argparse
import pandas as pd
from tqdm import tqdm

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
                groups[current_group].extend(map(int, line.strip().split()))
    return groups

# === Parse arguments ===
parser = argparse.ArgumentParser(description="Protein-protein hydrogen bond occupancy")
parser.add_argument("-f", dest="traj", required=True, help="Trajectory file")
parser.add_argument("-s", dest="top", required=True, help="Topology file")
parser.add_argument("-n", dest="ndx", required=True, help="GROMACS index file")
parser.add_argument("-g1", dest="group1", required=True, help="First group name")
parser.add_argument("-g2", dest="group2", required=True, help="Second group name")
parser.add_argument("-b", dest="start", type=int, default=0, help="Starting frame index")
parser.add_argument("-o", dest="outprefix", type=str, default="prot_prot", help="Output prefix")
args = parser.parse_args()

# === Load trajectory and index groups ===
traj = md.load(args.traj, top=args.top)
top = traj.topology
groups = load_ndx_groups(args.ndx)

g1_atoms = set(groups[args.group1])
g2_atoms = set(groups[args.group2])
frame_range = range(args.start, traj.n_frames)

# === Track H-bonding residue-residue pairs per frame ===
pair_frames = defaultdict(set)
detailed_pairs = defaultdict(set)

for i in tqdm(frame_range, desc="Analyzing frames"):
    hbonds = md.baker_hubbard(traj[i], freq=0.1, periodic=False)
    for d_idx, h_idx, a_idx in hbonds:
        d_atom = traj.topology.atom(d_idx)
        a_atom = traj.topology.atom(a_idx)
        d_res, a_res = d_atom.residue, a_atom.residue

        # Check if atoms belong to group1 and group2 (in any order)
        if ((d_idx in g1_atoms and a_idx in g2_atoms) or (a_idx in g1_atoms and d_idx in g2_atoms)):
            res_pair = tuple(sorted([(d_res, d_atom.name), (a_res, a_atom.name)], key=lambda x: x[0].index))
            resnames = tuple(sorted([d_res, a_res], key=lambda x: x.index))
            pair_frames[resnames].add(i)
            detailed_pairs[res_pair].add(i)

# === Write summary text file ===
txt_file = f"{args.outprefix}_hbonds.txt"
with open(txt_file, "w") as f:
    f.write("Chain1 - Chain2 - Occupancy %\n")
    for (r1, r2), frames in sorted(pair_frames.items(), key=lambda x: -len(x[1])):
        if r1.chain.index == r2.chain.index:
            continue
        occ = 100 * len(frames) / (traj.n_frames - args.start)
        chain1 = r1.chain.id if hasattr(r1.chain, 'id') and r1.chain.id else chr(65 + r1.chain.index)
        chain2 = r2.chain.id if hasattr(r2.chain, 'id') and r2.chain.id else chr(65 + r2.chain.index)
        f.write(f"{r1.name}{r1.resSeq}({chain1}) - {r2.name}{r2.resSeq}({chain2}) - {occ:.1f}%\n")
print(f"✅ Saved: {txt_file}")

# === Write detailed CSV file ===
csv_file = f"{args.outprefix}_hbonds_detailed.csv"
rows = []
for (r1, a1), (r2, a2) in detailed_pairs.keys():
    if r1.chain.index == r2.chain.index:
        continue
    frames = detailed_pairs[((r1, a1), (r2, a2))]
    occ = 100 * len(frames) / (traj.n_frames - args.start)
    chain1 = r1.chain.id if hasattr(r1.chain, 'id') and r1.chain.id else str(r1.chain.index)
    chain2 = r2.chain.id if hasattr(r2.chain, 'id') and r2.chain.id else str(r2.chain.index)
    rows.append([
        f"{r1.name}{r1.resSeq}({chain1})", a1,
        f"{r2.name}{r2.resSeq}({chain2})", a2,
        f"{occ:.1f}"
    ])
df = pd.DataFrame(rows, columns=["Donor_Res", "Donor_Atom", "Acceptor_Res", "Acceptor_Atom", "Occupancy"])
df.to_csv(csv_file, index=False)
print(f"✅ Saved: {csv_file}")

# === Generate PML visualization script ===
pml_file = f"{args.outprefix}_hbonds.pml"
hb_residues = set()

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
    f.write("set label_size, 14\n")

    for (r1, r2), frames in pair_frames.items():
        if r1.chain.index == r2.chain.index:
            continue
        occ = 100 * len(frames) / (traj.n_frames - args.start)
        for res in [r1, r2]:
            resn = res.name
            resid = res.resSeq
            chain = res.chain.id if hasattr(res.chain, 'id') and res.chain.id else chr(65 + res.chain.index)
            label = f"{resn}{resid}({chain}) - {occ:.1f}%"
            sele_name = f"hb_{resn}{resid}_{chain}"
            if sele_name not in hb_residues:
                f.write(f"select {sele_name}, chain {chain} and resi {resid} and resn {resn}\n")
                f.write(f"show sticks, {sele_name}\n")
                f.write(f"label {sele_name} and name CA, \"{label}\"\n")
                hb_residues.add(sele_name)

    f.write("zoom or ".join(hb_residues) + "\n")
    f.write("orient\n")

print(f"✅ Generated PyMOL script: {pml_file}")
