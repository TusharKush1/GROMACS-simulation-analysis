#!/usr/bin/env python3
"""Hydrophobic contact occupancy between two protein groups.

Two hydrophobic residues (Ala, Val, Leu, Ile, Met, Phe, Pro, Trp, Tyr) are in
contact in a frame when any pair of their side-chain carbon or sulfur atoms
is within the cutoff. Backbone atoms and polar side-chain atoms are not used.

Output:
  <prefix>_hydrophobic.txt   residue pairs and occupancy
  <prefix>_hydrophobic.pml   PyMOL script showing the residues

Author: Tushar Kushwaha
"""

from pp_common import (build_parser, contact_occupancy, filter_and_sort,
                       keep_pair, load_system, sorted_pair, write_pml,
                       write_summary)

HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "PRO", "TRP", "TYR"}
BACKBONE = {"N", "CA", "C", "O", "OXT"}


def hydrophobic_atoms(top, group):
    """Side-chain carbon and sulfur atoms of hydrophobic residues in `group`."""
    atoms = []
    for idx in group:
        atom = top.atom(idx)
        if (atom.residue.name in HYDROPHOBIC
                and atom.name not in BACKBONE
                and atom.element is not None
                and atom.element.symbol in ("C", "S")):
            atoms.append(atom)
    return atoms


def main():
    parser = build_parser("Protein-protein hydrophobic contact occupancy",
                          min_occ=5.0, cutoff=0.5)
    args = parser.parse_args()

    traj, g1, g2 = load_system(args)
    top = traj.topology

    pairs, keys = [], []
    atoms2 = hydrophobic_atoms(top, g2)
    for a1 in hydrophobic_atoms(top, g1):
        for a2 in atoms2:
            if keep_pair(a1.residue, a2.residue, args.same_chain):
                pairs.append([a1.index, a2.index])
                keys.append(sorted_pair(a1.residue, a2.residue))

    print(f"Checking {len(pairs)} possible hydrophobic atom pairs")
    occupancy = contact_occupancy(traj, pairs, keys, args.cutoff)
    contacts = filter_and_sort(occupancy, args.min_occ)

    write_summary(f"{args.outprefix}_hydrophobic.txt",
                  "Res1 - Res2 - Occupancy (%)", contacts)
    write_pml(f"{args.outprefix}_hydrophobic.pml", args.top, "hp", contacts)


if __name__ == "__main__":
    main()
