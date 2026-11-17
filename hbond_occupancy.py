#!/usr/bin/env python3
"""Residue-wise protein-ligand hydrogen bond occupancy from a GROMACS trajectory.

For every protein residue that hydrogen bonds to the ligand, report the
percentage of analysed frames in which at least one such bond is present.
Hydrogen bonds are detected with the Baker-Hubbard criteria as implemented
in MDTraj.

The trajectory should already be processed (molecules made whole and the
complex centred, e.g. with `gmx trjconv -pbc mol -center`), because periodic
images are not considered.

Author: Tushar Kushwaha
"""

import argparse
from collections import defaultdict

import mdtraj as md
import pandas as pd
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Protein-ligand residue-wise H-bond occupancy (text + CSV)"
    )
    parser.add_argument("--traj", type=str, default="md_0_10_fit.xtc",
                        help="Trajectory file name (default: md_0_10_fit.xtc)")
    parser.add_argument("--top", type=str, default="start.gro",
                        help="Topology/structure file matching the trajectory "
                             "(default: start.gro)")
    parser.add_argument("--ligand", type=str, default="UNL",
                        help="Residue name of the ligand (default: UNL)")
    parser.add_argument("--start", type=int, default=0,
                        help="Starting frame index (default: 0)")
    parser.add_argument("--out", type=str,
                        default="protein_ligand_residuewise_hbonds",
                        help="Output file prefix "
                             "(default: protein_ligand_residuewise_hbonds)")
    return parser.parse_args()


def residue_hbond_frames(traj, ligand, start_frame=0):
    """Map each protein residue to the set of frames where it H-bonds the ligand."""
    topology = traj.topology
    residue_frames = defaultdict(set)

    for i in tqdm(range(start_frame, traj.n_frames), desc="Processing frames"):
        hbonds = md.baker_hubbard(traj[i], freq=0.1, periodic=False)
        for donor, _hydrogen, acceptor in hbonds:
            d_res = topology.atom(donor).residue
            a_res = topology.atom(acceptor).residue

            if d_res.is_protein and a_res.name == ligand:
                residue_frames[d_res].add(i)
            elif a_res.is_protein and d_res.name == ligand:
                residue_frames[a_res].add(i)

    return residue_frames


def occupancy_table(residue_frames, n_analysed):
    """Build a table of occupancy per residue, highest occupancy first."""
    rows = [
        [res.chain.index, res.name.capitalize(), res.resSeq,
         100 * len(frames) / n_analysed]
        for res, frames in residue_frames.items()
    ]
    df = pd.DataFrame(rows, columns=[
        "Protein_Chain", "Protein_ResName", "Protein_ResID", "Occupancy_Percent"
    ])
    return df.sort_values("Occupancy_Percent", ascending=False).reset_index(drop=True)


def main():
    args = parse_args()

    print(f"Using trajectory: {args.traj}")
    print(f"Using topology: {args.top}")
    print(f"Ligand residue name: {args.ligand}")
    print(f"Starting from frame: {args.start}")

    traj = md.load(args.traj, top=args.top)

    n_analysed = traj.n_frames - args.start
    if n_analysed <= 0:
        raise SystemExit(
            f"--start {args.start} leaves no frames to analyse "
            f"(trajectory has {traj.n_frames} frames)."
        )
    if not any(res.name == args.ligand for res in traj.topology.residues):
        raise SystemExit(
            f"No residue named '{args.ligand}' in {args.top}. "
            "Set the ligand residue name with --ligand."
        )

    residue_frames = residue_hbond_frames(traj, args.ligand, args.start)
    df = occupancy_table(residue_frames, n_analysed)

    txt_file = f"{args.out}.txt"
    with open(txt_file, "w") as f:
        for row in df.itertuples(index=False):
            f.write(f"{row.Protein_ResName}{row.Protein_ResID} - "
                    f"{row.Occupancy_Percent:.1f}%\n")
    print(f"Saved: {txt_file}")

    csv_file = f"{args.out}.csv"
    df.to_csv(csv_file, index=False)
    print(f"Saved: {csv_file}")


if __name__ == "__main__":
    main()
