"""Shared helpers for the protein-protein interaction scripts.

Every pp_*.py script takes the same inputs (trajectory, structure, GROMACS
index file and two index groups) and writes the same kinds of output (a text
summary and a PyMOL script), so the common parts live here.

Author: Tushar Kushwaha
"""

import argparse

import mdtraj as md
import numpy as np

HIS_NAMES = {"HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP"}

# Atoms that define the aromatic ring used for centroid calculations.
RING_ATOMS = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
}
RING_ATOMS.update({name: ["CG", "ND1", "CD2", "CE1", "NE2"] for name in HIS_NAMES})


# --------------------------------------------------------------------------
# Input
# --------------------------------------------------------------------------

def build_parser(description, min_occ, cutoff=None):
    """Argument parser with the options shared by all pp_*.py scripts."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-f", dest="traj", required=True,
                        help="Trajectory file (e.g. .xtc), already processed "
                             "so that the complex is whole")
    parser.add_argument("-s", dest="top", required=True,
                        help="Structure file matching the trajectory; use a "
                             ".pdb with chain IDs so the two partners are "
                             "recognised as separate chains")
    parser.add_argument("-n", dest="ndx", required=True,
                        help="GROMACS index file")
    parser.add_argument("-g1", dest="group1", required=True,
                        help="Name of the first index group")
    parser.add_argument("-g2", dest="group2", required=True,
                        help="Name of the second index group")
    parser.add_argument("-b", dest="start", type=int, default=0,
                        help="First frame to analyse (default: 0)")
    parser.add_argument("-o", dest="outprefix", type=str, default="pp",
                        help="Output prefix (default: pp)")
    parser.add_argument("--min-occ", dest="min_occ", type=float, default=min_occ,
                        help="Report only interactions present in at least "
                             f"this percentage of frames (default: {min_occ})")
    parser.add_argument("--same-chain", dest="same_chain", action="store_true",
                        help="Also report pairs within the same chain "
                             "(by default only inter-chain pairs are kept)")
    if cutoff is not None:
        parser.add_argument("--cutoff", type=float, default=cutoff,
                            help=f"Distance cutoff in nm (default: {cutoff})")
    return parser


def load_ndx_groups(ndx_file):
    """Read a GROMACS index file into {group name: [0-based atom indices]}.

    GROMACS numbers atoms from 1 and MDTraj from 0, so 1 is subtracted here.
    """
    groups = {}
    current = None
    with open(ndx_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("["):
                current = line[1:-1].strip()
                groups[current] = []
            elif current is not None and line:
                groups[current].extend(int(x) - 1 for x in line.split())
    return groups


def load_system(args):
    """Load the trajectory and the two index groups named in `args`.

    Returns (traj, group1, group2). Frames before `args.start` are dropped,
    so occupancies are always relative to the analysed frames.
    """
    traj = md.load(args.traj, top=args.top)
    if not 0 <= args.start < traj.n_frames:
        raise SystemExit(
            f"-b {args.start} leaves no frames to analyse "
            f"(trajectory has {traj.n_frames} frames)."
        )
    if args.start:
        traj = traj[args.start:]

    groups = load_ndx_groups(args.ndx)
    selected = []
    for name in (args.group1, args.group2):
        if name not in groups:
            raise SystemExit(
                f"Group '{name}' not found in {args.ndx}. "
                f"Available groups: {', '.join(groups)}"
            )
        atoms = groups[name]
        if not atoms:
            raise SystemExit(f"Group '{name}' in {args.ndx} is empty.")
        if max(atoms) >= traj.topology.n_atoms or min(atoms) < 0:
            raise SystemExit(
                f"Group '{name}' refers to atoms outside the structure "
                f"({traj.topology.n_atoms} atoms). Was the index file made "
                "from a different system?"
            )
        print(f"Atoms in {name}: {len(atoms)}")
        selected.append(atoms)

    print(f"Frames analysed: {traj.n_frames}")
    return traj, selected[0], selected[1]


# --------------------------------------------------------------------------
# Residues and chains
# --------------------------------------------------------------------------

def chain_id(chain):
    """Chain ID from the structure file, or None if it has none."""
    for attr in ("chain_id", "id"):
        cid = getattr(chain, attr, None)
        if isinstance(cid, str) and cid.strip():
            return cid.strip()
    return None


def chain_label(chain):
    """Chain ID if present, otherwise A, B, C... by order in the file."""
    return chain_id(chain) or chr(65 + chain.index)


def res_label(res):
    return f"{res.name}{res.resSeq}({chain_label(res.chain)})"


def keep_pair(res1, res2, same_chain=False):
    """True if this residue pair should be counted."""
    if res1.index == res2.index:
        return False
    return same_chain or res1.chain.index != res2.chain.index


def sorted_pair(res1, res2):
    """Residue pair in a fixed order, so A-B and B-A are the same key."""
    return (res1, res2) if res1.index <= res2.index else (res2, res1)


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------

def contact_occupancy(traj, pairs, keys, cutoff, max_distances=20_000_000):
    """Occupancy (%) of each key, from distances between atom pairs.

    `pairs` is a list of [atom_i, atom_j] and `keys` gives, for each pair, the
    residue pair (or any hashable) it belongs to. A key is in contact in a
    frame if any of its atom pairs is closer than `cutoff` (nm). Frames are
    processed in chunks so that memory use stays bounded.
    """
    if not pairs:
        return {}

    key_ids = {}
    ids = np.array([key_ids.setdefault(k, len(key_ids)) for k in keys])
    order = np.argsort(ids, kind="stable")
    pairs_sorted = np.asarray(pairs)[order]
    starts = np.unique(ids[order], return_index=True)[1]

    counts = np.zeros(len(key_ids), dtype=int)
    chunk = max(1, max_distances // len(pairs_sorted))
    for begin in range(0, traj.n_frames, chunk):
        dist = md.compute_distances(traj[begin:begin + chunk], pairs_sorted)
        hit = np.logical_or.reduceat(dist < cutoff, starts, axis=1)
        counts += hit.sum(axis=0)

    return {key: 100.0 * counts[i] / traj.n_frames for key, i in key_ids.items()}


def ring_centroids(traj, group):
    """Aromatic residues with an atom in `group`, and their ring centroids.

    Returns (residues, centroids) where centroids has shape
    (n_frames, n_residues, 3).
    """
    residues, ring_indices = [], []
    seen = set()
    for idx in group:
        res = traj.topology.atom(idx).residue
        if res.index in seen or res.name not in RING_ATOMS:
            continue
        seen.add(res.index)
        inds = [a.index for a in res.atoms if a.name in RING_ATOMS[res.name]]
        if inds:
            residues.append(res)
            ring_indices.append(inds)

    centroids = np.zeros((traj.n_frames, len(residues), 3))
    for j, inds in enumerate(ring_indices):
        centroids[:, j, :] = traj.xyz[:, inds, :].mean(axis=1)
    return residues, centroids


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------

def filter_and_sort(occupancy, min_occ):
    """[(key, occupancy)] at or above `min_occ`, highest occupancy first."""
    kept = [(k, occ) for k, occ in occupancy.items() if occ >= min_occ and occ > 0]
    return sorted(kept, key=lambda item: -item[1])


def write_summary(path, header, contacts):
    """Write `Res1 - Res2 - occupancy` lines for [((res1, res2), occ)]."""
    with open(path, "w") as f:
        f.write(header + "\n")
        for (res1, res2), occ in contacts:
            f.write(f"{res_label(res1)} - {res_label(res2)} - {occ:.1f}\n")
    print(f"Saved: {path}")


def write_pml(path, structure, tag, contacts, atom_names=None, dash_color="yellow"):
    """Write a PyMOL script that shows the interacting residues.

    `contacts` is [((res1, res2), occ)], highest occupancy first. Each residue
    is shown as sticks and labelled with its highest occupancy. If
    `atom_names` is given as a function (res1, res2) -> (name1, name2), a
    dashed line is drawn between those two atoms for every pair.
    """
    selections = []

    def sel_name(res):
        return f"{tag}_{res.name}{res.resSeq}_{chain_label(res.chain)}"

    def sel_expr(res):
        parts = [f"resi {res.resSeq}", f"resn {res.name}"]
        if chain_id(res.chain):
            parts.insert(0, f"chain {chain_id(res.chain)}")
        return " and ".join(parts)

    lines = [
        f"load {structure}, complex",
        "bg_color white",
        "set antialias, 2",
        "set ray_trace_mode, 1",
        "set ray_trace_gain, 0.1",
        "hide everything, all",
        "show cartoon, complex",
        "set stick_radius, 0.2",
        "set label_color, black",
        "set label_size, 14",
        f"set dash_color, {dash_color}",
        "set dash_width, 2",
        "set dash_gap, 0.3",
        "",
    ]

    for n, ((res1, res2), occ) in enumerate(contacts, start=1):
        for res in (res1, res2):
            sel = sel_name(res)
            if sel not in selections:
                selections.append(sel)
                lines.append(f"select {sel}, {sel_expr(res)}")
                lines.append(f"show sticks, {sel}")
                lines.append(f'label {sel} and name CA, "{res_label(res)} - {occ:.1f}%"')
        if atom_names is not None:
            name1, name2 = atom_names(res1, res2)
            dash = f"{tag}_dash_{n}"
            lines.append(f"distance {dash}, {sel_name(res1)} and name {name1}, "
                         f"{sel_name(res2)} and name {name2}")
            lines.append(f"hide labels, {dash}")

    if selections:
        lines.append("")
        lines.append("orient " + " or ".join(selections))
    lines.append("deselect")

    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Saved: {path}")
