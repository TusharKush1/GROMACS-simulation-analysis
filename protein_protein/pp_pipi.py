#!/usr/bin/env python3
"""Pi-pi interaction occupancy between two protein groups.

Aromatic rings of Phe, Tyr, Trp (six-membered ring) and His are compared
between the two groups. A pair is counted in a frame when the distance
between the two ring centroids is below the cutoff. Ring orientation is not
checked, so stacked and T-shaped arrangements are both included.

Output:
  <prefix>_pipi.txt   residue pairs and occupancy
  <prefix>_pipi.pml   PyMOL script showing the residues

Author: Tushar Kushwaha
"""

import numpy as np

from pp_common import (build_parser, filter_and_sort, keep_pair, load_system,
                       ring_centroids, sorted_pair, write_pml, write_summary)


def main():
    parser = build_parser("Protein-protein pi-pi interaction occupancy",
                          min_occ=0.0, cutoff=0.6)
    args = parser.parse_args()

    traj, g1, g2 = load_system(args)

    res1, cent1 = ring_centroids(traj, g1)
    res2, cent2 = ring_centroids(traj, g2)
    print(f"Aromatic rings: {len(res1)} in {args.group1}, "
          f"{len(res2)} in {args.group2}")

    occupancy = {}
    for j, r1 in enumerate(res1):
        for k, r2 in enumerate(res2):
            if not keep_pair(r1, r2, args.same_chain):
                continue
            dist = np.linalg.norm(cent1[:, j, :] - cent2[:, k, :], axis=1)
            occ = 100.0 * np.mean(dist < args.cutoff)
            key = sorted_pair(r1, r2)
            occupancy[key] = max(occupancy.get(key, 0.0), occ)

    contacts = filter_and_sort(occupancy, args.min_occ)

    write_summary(f"{args.outprefix}_pipi.txt",
                  "Res1 - Res2 - Occupancy (%)", contacts)
    write_pml(f"{args.outprefix}_pipi.pml", args.top, "pipi", contacts)


if __name__ == "__main__":
    main()
