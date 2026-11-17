#!/usr/bin/env python3
"""Salt bridge occupancy between two protein groups.

A salt bridge is counted in a frame when a side-chain oxygen of Asp or Glu
(OD*, OE*) is within the cutoff of a side-chain nitrogen of Lys (NZ),
Arg (NE, NH*) or His (ND1, NE2). Acidic and basic residues are searched in
both groups, so the order of -g1 and -g2 does not matter.

Histidine is counted as basic whatever its protonation state; use --no-his
to leave it out.

Output:
  <prefix>_sb.txt   residue pairs and occupancy
  <prefix>_sb.pml   PyMOL script showing the residues

Author: Tushar Kushwaha
"""

from pp_common import (HIS_NAMES, build_parser, contact_occupancy,
                       filter_and_sort, keep_pair, load_system, write_pml,
                       write_summary)

ACIDIC = {"ASP": ("OD",), "GLU": ("OE",)}
BASIC = {"LYS": ("NZ",), "ARG": ("NE", "NH")}
HIS_ATOMS = ("ND1", "NE2")


def charged_atoms(top, group, table):
    """Atoms in `group` whose residue and atom-name prefix are in `table`."""
    atoms = []
    for idx in group:
        atom = top.atom(idx)
        prefixes = table.get(atom.residue.name)
        if prefixes and atom.name.startswith(prefixes):
            atoms.append(atom)
    return atoms


def main():
    parser = build_parser("Protein-protein salt bridge occupancy",
                          min_occ=5.0, cutoff=0.4)
    parser.add_argument("--no-his", dest="no_his", action="store_true",
                        help="Do not count histidine as a basic residue")
    args = parser.parse_args()

    traj, g1, g2 = load_system(args)
    top = traj.topology

    basic = dict(BASIC)
    if not args.no_his:
        basic.update({name: HIS_ATOMS for name in HIS_NAMES})

    # Acid in one group and base in the other, in both directions.
    pairs, keys = [], []
    for acid_group, base_group in ((g1, g2), (g2, g1)):
        acids = charged_atoms(top, acid_group, ACIDIC)
        bases = charged_atoms(top, base_group, basic)
        for acid in acids:
            for base in bases:
                if keep_pair(acid.residue, base.residue, args.same_chain):
                    pairs.append([acid.index, base.index])
                    keys.append((acid.residue, base.residue))

    print(f"Checking {len(pairs)} possible salt bridge atom pairs")
    occupancy = contact_occupancy(traj, pairs, keys, args.cutoff)
    contacts = filter_and_sort(occupancy, args.min_occ)

    write_summary(f"{args.outprefix}_sb.txt",
                  "Acidic_Res - Basic_Res - Occupancy (%)", contacts)
    write_pml(f"{args.outprefix}_sb.pml", args.top, "sb", contacts,
              atom_names=lambda acid, base: ("CA", "CA"), dash_color="yellow")


if __name__ == "__main__":
    main()
