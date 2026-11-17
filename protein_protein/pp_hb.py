#!/usr/bin/env python3
"""Hydrogen bond occupancy between two protein groups.

Hydrogen bonds are detected in every frame with the Baker-Hubbard criteria
(MDTraj). A residue pair is counted in a frame if at least one hydrogen bond
links the two residues.

Output:
  <prefix>_hbonds.txt            residue pairs and occupancy
  <prefix>_hbonds_detailed.csv   donor and acceptor atoms and occupancy
  <prefix>_hbonds.pml            PyMOL script showing the residues

Author: Tushar Kushwaha
"""

from collections import defaultdict

import mdtraj as md
import pandas as pd
from tqdm import tqdm

from pp_common import (build_parser, filter_and_sort, keep_pair, load_system,
                       res_label, sorted_pair, write_pml, write_summary)


def main():
    parser = build_parser("Protein-protein hydrogen bond occupancy", min_occ=0.0)
    args = parser.parse_args()

    traj, g1, g2 = load_system(args)
    g1, g2 = set(g1), set(g2)
    top = traj.topology

    pair_frames = defaultdict(set)   # (res, res) -> frames
    atom_frames = defaultdict(set)   # (donor atom, acceptor atom) -> frames

    for i in tqdm(range(traj.n_frames), desc="Analyzing frames"):
        hbonds = md.baker_hubbard(traj[i], freq=0.1, periodic=False)
        for donor, _hydrogen, acceptor in hbonds:
            donor, acceptor = int(donor), int(acceptor)
            across = ((donor in g1 and acceptor in g2) or
                      (donor in g2 and acceptor in g1))
            if not across:
                continue
            d_atom, a_atom = top.atom(donor), top.atom(acceptor)
            if not keep_pair(d_atom.residue, a_atom.residue, args.same_chain):
                continue
            pair_frames[sorted_pair(d_atom.residue, a_atom.residue)].add(i)
            atom_frames[(d_atom, a_atom)].add(i)

    def percent(frames):
        return 100.0 * len(frames) / traj.n_frames

    contacts = filter_and_sort(
        {pair: percent(frames) for pair, frames in pair_frames.items()},
        args.min_occ)
    write_summary(f"{args.outprefix}_hbonds.txt",
                  "Res1 - Res2 - Occupancy (%)", contacts)

    detailed = filter_and_sort(
        {atoms: percent(frames) for atoms, frames in atom_frames.items()},
        args.min_occ)
    rows = [[res_label(d.residue), d.name, res_label(a.residue), a.name, f"{occ:.1f}"]
            for (d, a), occ in detailed]
    csv_file = f"{args.outprefix}_hbonds_detailed.csv"
    pd.DataFrame(rows, columns=["Donor_Res", "Donor_Atom", "Acceptor_Res",
                                "Acceptor_Atom", "Occupancy"]
                 ).to_csv(csv_file, index=False)
    print(f"Saved: {csv_file}")

    write_pml(f"{args.outprefix}_hbonds.pml", args.top, "hb", contacts)


if __name__ == "__main__":
    main()
