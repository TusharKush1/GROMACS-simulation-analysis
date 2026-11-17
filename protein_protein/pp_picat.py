#!/usr/bin/env python3
"""Cation-pi interaction occupancy between two protein groups.

A cation-pi contact is counted in a frame when a cationic nitrogen of Lys
(NZ) or Arg (NE, NH1, NH2) is within the cutoff of the ring centroid of Phe,
Tyr, Trp (six-membered ring) or His. Rings and cations are searched in both
groups, so the order of -g1 and -g2 does not matter.

Output:
  <prefix>_pi_cation.txt   residue pairs and occupancy
  <prefix>_pi_cation.pml   PyMOL script showing the residues

Author: Tushar Kushwaha
"""

from collections import defaultdict

import numpy as np

from pp_common import (build_parser, filter_and_sort, keep_pair, load_system,
                       ring_centroids, write_pml, write_summary)

CATION_ATOMS = {"LYS": {"NZ"}, "ARG": {"NE", "NH1", "NH2"}}
# Atoms joined by the dashed line in the PyMOL script.
RING_DASH_ATOM = {"TRP": "CD2"}          # CG for the other rings
CATION_DASH_ATOM = {"LYS": "NZ", "ARG": "CZ"}


def cation_atoms(top, group):
    """{residue: [atom indices]} for the cationic nitrogens in `group`."""
    cations = defaultdict(list)
    for idx in group:
        atom = top.atom(idx)
        if atom.name in CATION_ATOMS.get(atom.residue.name, ()):
            cations[atom.residue].append(idx)
    return cations


def main():
    parser = build_parser("Protein-protein cation-pi interaction occupancy",
                          min_occ=5.0, cutoff=0.6)
    args = parser.parse_args()

    traj, g1, g2 = load_system(args)
    top = traj.topology

    # Ring in one group and cation in the other, in both directions.
    occupancy = {}
    for ring_group, cation_group in ((g1, g2), (g2, g1)):
        rings, centroids = ring_centroids(traj, ring_group)
        cations = cation_atoms(top, cation_group)
        for j, ring_res in enumerate(rings):
            centroid = centroids[:, j, :]
            for cat_res, atom_inds in cations.items():
                if not keep_pair(ring_res, cat_res, args.same_chain):
                    continue
                # distance from every cationic N to the centroid, per frame
                diff = traj.xyz[:, atom_inds, :] - centroid[:, None, :]
                close = np.linalg.norm(diff, axis=2) < args.cutoff
                occ = 100.0 * np.mean(close.any(axis=1))
                key = (ring_res, cat_res)
                occupancy[key] = max(occupancy.get(key, 0.0), occ)

    contacts = filter_and_sort(occupancy, args.min_occ)

    write_summary(f"{args.outprefix}_pi_cation.txt",
                  "Aromatic_Res - Cationic_Res - Occupancy (%)", contacts)
    write_pml(f"{args.outprefix}_pi_cation.pml", args.top, "pc", contacts,
              atom_names=lambda ring, cat: (RING_DASH_ATOM.get(ring.name, "CG"),
                                            CATION_DASH_ATOM[cat.name]),
              dash_color="orange")


if __name__ == "__main__":
    main()
