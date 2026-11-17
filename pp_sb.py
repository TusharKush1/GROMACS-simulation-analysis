import mdtraj as md
import argparse
from collections import defaultdict
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


def find_salt_bridges(traj, group1, group2, distance_cutoff=0.4):
    from collections import defaultdict
    from tqdm import tqdm

    acidic = {'ASP', 'GLU'}
    basic = {'ARG', 'LYS', 'HIS'}

    # Filter acidic atoms in group1
    acidic_atoms = []
    for idx in group1:
        atom = traj.topology.atom(idx)
        if atom.residue.name in acidic and (atom.name.startswith("OD") or atom.name.startswith("OE")):
            acidic_atoms.append(idx)

    # Filter basic atoms in group2
    basic_atoms = []
    for idx in group2:
        atom = traj.topology.atom(idx)
        if atom.residue.name in basic and (atom.name.startswith("NZ") or atom.name.startswith("NH") or atom.name.startswith("NE")):
            basic_atoms.append(idx)

    # Create all acidic-basic atom pairs
    pairs = []
    pair_info = {}
    for a1 in acidic_atoms:
        atom1 = traj.topology.atom(a1)
        res1 = atom1.residue
        for a2 in basic_atoms:
            atom2 = traj.topology.atom(a2)
            res2 = atom2.residue
            if res1.index == res2.index:
                continue
            key = tuple(sorted([res1, res2], key=lambda x: x.index))
            pairs.append([a1, a2])
            pair_info[(a1, a2)] = key

    print(f"Checking {len(pairs)} possible salt bridge pairs")

    # Compute distances for all frames in batch
    distances = md.compute_distances(traj, pairs)  # shape: (n_frames, n_pairs)

    salt_bridges = defaultdict(set)
    for frame_idx in tqdm(range(traj.n_frames), desc="Analyzing salt bridges"):
        for pair_idx, dist in enumerate(distances[frame_idx]):
            if dist < distance_cutoff:
                atom_pair = tuple(pairs[pair_idx])
                key = pair_info[atom_pair]
                salt_bridges[key].add(frame_idx)

    return salt_bridges


    
# === Call function ===
salt_bridges = find_salt_bridges(traj, g1_atoms, g2_atoms)

# === Write salt bridge summary (.txt) ===
txt_file = f"{args.outprefix}_sb.txt"
with open(txt_file, "w") as f:
    f.write("Acidic_Res - Basic_Res - Occupancy (%)\n")
    for res1, res2 in sorted(salt_bridges, key=lambda x: -len(salt_bridges[x])):
        occ = 100 * len(salt_bridges[(res1, res2)]) / traj.n_frames
        # Skip low‑occupancy or same‑chain pairs
        if occ < 5.0 or res1.chain.index == res2.chain.index:
            continue
        c1 = chr(65 + res1.chain.index)
        c2 = chr(65 + res2.chain.index)
        f.write(f"{res1.name}{res1.resSeq}({c1}) - "
                f"{res2.name}{res2.resSeq}({c2}) - {occ:.1f}\n")
print(f"✅ Saved salt bridge summary: {txt_file}")

# === Generate PyMOL visualization script (.pml) ===
pml_file = f"{args.outprefix}_sb.pml"
sb_residues = set()
with open(pml_file, 'w') as f:
    f.write(f"load {args.top}, complex\n")
    f.write("bg_color white\nhide everything, all\nshow cartoon, complex\n")
    f.write("set stick_radius, 0.2\nset label_color, black\nset label_size, 14\n")
    f.write("set dash_color, yellow\nset dash_width, 2\nset dash_gap, 0.3\n")

    for res1, res2 in sorted(salt_bridges, key=lambda x: -len(salt_bridges[x])):
        occ = 100 * len(salt_bridges[(res1, res2)]) / traj.n_frames
        if occ < 5.0 or res1.chain.index == res2.chain.index:
            continue
        c1 = chr(65 + res1.chain.index)
        c2 = chr(65 + res2.chain.index)

        sel1 = f"sb_{res1.name}{res1.resSeq}_{c1}"
        sel2 = f"sb_{res2.name}{res2.resSeq}_{c2}"

        # highlight the residues
        if sel1 not in sb_residues:
            f.write(f"select {sel1}, chain {c1} and resi {res1.resSeq}\n")
            f.write(f"show sticks, {sel1}\n")
            sb_residues.add(sel1)
        if sel2 not in sb_residues:
            f.write(f"select {sel2}, chain {c2} and resi {res2.resSeq}\n")
            f.write(f"show sticks, {sel2}\n")
            sb_residues.add(sel2)

        # draw dashed line
        line_name = f"sb_line_{res1.resSeq}_{res2.resSeq}"
        f.write(f"distance {line_name}, {sel1}/CA, {sel2}/CA\n")
        f.write(f"set dash_color, yellow, {line_name}\n")
        f.write(f"label {line_name}, \"{occ:.1f}%\"\n")

    if sb_residues:
        f.write("zoom " + " or ".join(sb_residues) + "\n")
    f.write("orient\n")

print(f"✅ Generated PyMOL script: {pml_file}")

