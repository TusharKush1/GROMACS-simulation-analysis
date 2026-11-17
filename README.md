# GROMACS-simulation-analysis

Shell and Python scripts I use for the analysis and presentation of GROMACS molecular dynamics simulation results.

## Scripts

| Script | What it does |
| --- | --- |
| [`hbond_occupancy.py`](hbond_occupancy.py) | Residue-wise protein–ligand hydrogen bond occupancy from a processed trajectory |

## Setup

```bash
git clone https://github.com/TusharKush1/GROMACS-simulation-analysis.git
cd GROMACS-simulation-analysis
pip3 install -r requirements.txt
```

## hbond_occupancy.py

For every protein residue that hydrogen bonds to the ligand, the script reports the percentage of analysed frames in which at least one such bond is present. Hydrogen bonds are detected with the Baker–Hubbard criteria as implemented in [MDTraj](https://www.mdtraj.org).

### Input

- A processed trajectory (`.xtc`). Make molecules whole and centre the complex first, because periodic images are not considered:

  ```bash
  gmx trjconv -s md_0_10.tpr -f md_0_10.xtc -o md_0_10_fit.xtc -pbc mol -center
  ```

- A structure file with the same atoms as the trajectory (`.gro`).

### Usage

```bash
python3 hbond_occupancy.py --traj md_0_10_fit.xtc --top start.gro --ligand UNL
```

| Option | Default | Meaning |
| --- | --- | --- |
| `--traj` | `md_0_10_fit.xtc` | Trajectory file |
| `--top` | `start.gro` | Structure file matching the trajectory |
| `--ligand` | `UNL` | Residue name of the ligand |
| `--start` | `0` | First frame to analyse, for skipping equilibration |
| `--out` | `protein_ligand_residuewise_hbonds` | Output file prefix |

### Output

- `<prefix>.txt`: one line per residue, highest occupancy first, for example `Asp855 - 78.4%`
- `<prefix>.csv`: columns `Protein_Chain`, `Protein_ResName`, `Protein_ResID`, `Occupancy_Percent`

Residue numbers are those in the structure file.

## License

MIT, see [LICENSE](LICENSE).
