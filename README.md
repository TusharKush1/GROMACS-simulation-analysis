# GROMACS-simulation-analysis

Python scripts I use to analyse GROMACS molecular dynamics trajectories and present the results. Each script reports how often an interaction is present over the trajectory (its occupancy) and writes a PyMOL script that shows the interacting residues on the structure.

## Scripts

| Script | Interaction | Criterion | Default cutoff |
| --- | --- | --- | --- |
| [`protein_ligand/hbond_occupancy.py`](protein_ligand/hbond_occupancy.py) | Protein–ligand hydrogen bonds | Baker–Hubbard | – |
| [`protein_protein/pp_hb.py`](protein_protein/pp_hb.py) | Protein–protein hydrogen bonds | Baker–Hubbard | – |
| [`protein_protein/pp_sb.py`](protein_protein/pp_sb.py) | Salt bridges | Asp/Glu side-chain O to Lys/Arg/His side-chain N | 0.4 nm |
| [`protein_protein/pp_hydro.py`](protein_protein/pp_hydro.py) | Hydrophobic contacts | Side-chain C/S atoms of hydrophobic residues | 0.5 nm |
| [`protein_protein/pp_pipi.py`](protein_protein/pp_pipi.py) | π–π interactions | Ring centroid to ring centroid | 0.6 nm |
| [`protein_protein/pp_picat.py`](protein_protein/pp_picat.py) | Cation–π interactions | Lys/Arg cationic N to ring centroid | 0.6 nm |

`protein_protein/pp_common.py` holds the code shared by the `pp_*.py` scripts and must stay in the same folder.

## Setup

```bash
git clone https://github.com/TusharKush1/GROMACS-simulation-analysis.git
cd GROMACS-simulation-analysis
pip3 install -r requirements.txt
```

The scripts use [MDTraj](https://www.mdtraj.org), NumPy, pandas and tqdm. [PyMOL](https://pymol.org) is needed only to open the generated `.pml` files.

## Preparing the input

1. **Process the trajectory** so that the complex is whole and centred. Periodic images are not reconstructed by the scripts.

   ```bash
   gmx trjconv -s md.tpr -f md.xtc -o md_fit.xtc -pbc mol -center
   ```

2. **Provide a structure file** with the same atoms as the trajectory. For the protein–protein scripts use a `.pdb` that has chain IDs: by default only pairs from different chains are reported, and the PyMOL scripts select residues by chain.

3. **Make an index file** with one group per binding partner (protein–protein scripts only).

   ```bash
   gmx make_ndx -f complex.pdb -o index.ndx
   # at the prompt, for example:  chain A   then   chain B   then   q
   ```

## Protein–ligand hydrogen bonds

For every protein residue that hydrogen bonds to the ligand, `hbond_occupancy.py` reports the percentage of analysed frames in which at least one such bond is present.

```bash
python3 protein_ligand/hbond_occupancy.py --traj md_fit.xtc --top start.gro --ligand UNL
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--traj` | `md_0_10_fit.xtc` | Trajectory file |
| `--top` | `start.gro` | Structure file matching the trajectory |
| `--ligand` | `UNL` | Residue name of the ligand |
| `--start` | `0` | First frame to analyse, for skipping equilibration |
| `--outprefix` | `pro_lig` | Output prefix |
| `--pdb` | the `--top` file | Structure that the PyMOL script loads |

Output:

- `<prefix>_residuewise_hbonds.txt`: one line per residue, highest occupancy first
- `<prefix>_residuewise_hbonds.csv`: columns `Protein_Chain`, `Protein_ResName`, `Protein_ResID`, `Occupancy_Percent`
- `<prefix>_hbonds.pml`: PyMOL script showing the ligand and the hydrogen-bonding residues, labelled with occupancy

## Protein–protein interactions

The five `pp_*.py` scripts share one command line. Run them from any folder:

```bash
python3 protein_protein/pp_hb.py    -f md_fit.xtc -s complex.pdb -n index.ndx -g1 chA -g2 chB
python3 protein_protein/pp_sb.py    -f md_fit.xtc -s complex.pdb -n index.ndx -g1 chA -g2 chB
python3 protein_protein/pp_hydro.py -f md_fit.xtc -s complex.pdb -n index.ndx -g1 chA -g2 chB
python3 protein_protein/pp_pipi.py  -f md_fit.xtc -s complex.pdb -n index.ndx -g1 chA -g2 chB
python3 protein_protein/pp_picat.py -f md_fit.xtc -s complex.pdb -n index.ndx -g1 chA -g2 chB
```

| Option | Default | Meaning |
| --- | --- | --- |
| `-f` | required | Trajectory file |
| `-s` | required | Structure file matching the trajectory |
| `-n` | required | GROMACS index file |
| `-g1`, `-g2` | required | Names of the two index groups; their order does not matter |
| `-b` | `0` | First frame to analyse |
| `-o` | `pp` | Output prefix |
| `--min-occ` | `5` (`0` for `pp_hb.py` and `pp_pipi.py`) | Lowest occupancy (%) to report |
| `--cutoff` | see the table above | Distance cutoff in nm (not used by `pp_hb.py`) |
| `--same-chain` | off | Also report pairs within the same chain |
| `--no-his` | off | `pp_sb.py` only: do not count histidine as basic |

Each script writes a text summary with one line per residue pair, highest occupancy first:

```
Res1 - Res2 - Occupancy (%)
<RES><number>(<chain>) - <RES><number>(<chain>) - <occupancy>
```

and a PyMOL script that shows the residues as sticks with their occupancy as labels:

| Script | Summary | PyMOL script | Extra |
| --- | --- | --- | --- |
| `pp_hb.py` | `<prefix>_hbonds.txt` | `<prefix>_hbonds.pml` | `<prefix>_hbonds_detailed.csv` with donor and acceptor atoms |
| `pp_sb.py` | `<prefix>_sb.txt` | `<prefix>_sb.pml` | |
| `pp_hydro.py` | `<prefix>_hydrophobic.txt` | `<prefix>_hydrophobic.pml` | |
| `pp_pipi.py` | `<prefix>_pipi.txt` | `<prefix>_pipi.pml` | |
| `pp_picat.py` | `<prefix>_pi_cation.txt` | `<prefix>_pi_cation.pml` | |

Open a result in PyMOL from the folder that holds the structure file:

```bash
pymol pp_hbonds.pml
```

### How the interactions are defined

- **Occupancy** is the percentage of analysed frames in which a residue pair meets the criterion at least once, whichever atoms are involved.
- **Hydrogen bonds** use the Baker–Hubbard criteria as implemented in MDTraj, evaluated frame by frame.
- **Salt bridges**: an Asp `OD*` or Glu `OE*` atom within the cutoff of a Lys `NZ`, Arg `NE`/`NH*` or His `ND1`/`NE2` atom. Histidine is counted whatever its protonation state unless `--no-his` is given.
- **Hydrophobic contacts**: side-chain carbon or sulfur atoms of Ala, Val, Leu, Ile, Met, Phe, Pro, Trp and Tyr within the cutoff. Backbone atoms are not used.
- **π–π**: centroids of the rings of Phe, Tyr, Trp (six-membered ring) and His within the cutoff. Ring orientation is not checked, so stacked and T-shaped arrangements are both included.
- **Cation–π**: a Lys `NZ` or Arg `NE`/`NH1`/`NH2` atom within the cutoff of a ring centroid.

## License

MIT, see [LICENSE](LICENSE).
