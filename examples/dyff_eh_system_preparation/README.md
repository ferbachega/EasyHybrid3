# From XYZ / MOL2 / PDB to a DYFF-1.0-eh system (.pkl)

`prepare_dyff_eh.py` builds a pDynamo system with EasyHybrid's DYFF force field from a structure file, step by step, and saves it as a `.pkl`. It runs the same steps as the GUI, *System > Prepare DYFF System* followed by *DYFF + residue library*.

```
python3 prepare_dyff_eh.py                              # the three files in inputs/
python3 prepare_dyff_eh.py inputs/theophylline.xyz
python3 prepare_dyff_eh.py inputs/valacyclovir.mol2
python3 prepare_dyff_eh.py inputs/4oas_raw.pdb --minimize
python3 prepare_dyff_eh.py my_molecule.xyz --ph 5.0 --out my_output
```

Options:

| option | effect |
|---|---|
| `--ph X` | pH for protonation (default 7.4) |
| `--no-solvent` | CM5 charges in the gas phase (default: implicit water, ALPB) |
| `--charges-only` | for atoms outside the library: CM5 charges only, DYFF bonded rules (faster) |
| `--minimize` | LBFGS minimization: hydrogens first, then all atoms |
| `--threads N` | threads for xTB (default 4) |
| `--out DIR` | output folder (default `output/`) |

## The steps

1. **Import** (`pBabel.ImportSystem`). Each format carries different information:
   - XYZ: coordinates only;
   - MOL2: atoms, bonds and bond orders;
   - PDB: residues and atom names.
2. **Prepare** (`pdynamo/dyff_prep.py`):
   - standard residues are recognized;
   - anything else (ligands, or the molecule of an XYZ/MOL2 file) goes through EasyHybrid's native chemistry perception: connectivity, bond orders and aromaticity, plus hydrogen counts from the geometry when the file has no hydrogens;
   - the structure is protonated at the chosen pH and hydrogens are added;
   - truncated residues are completed and metal sites detected.
   
   A structure that already has hydrogens keeps them, with its protonation as given.
3. **Residue library** (`pdynamo/residue_library.plan`):
   - residues that match the library by bond graph get precomputed xTB charges and bonded terms (amino acids, nucleotides, ...);
   - waters get TIP3P charges and ions their formal charge.
4. **xTB parameters for the rest** (`residue_library.rest_parameters`), computed on a capped subsystem:
   - CM5 charges from GFN2-xTB;
   - bonds and angles from the xTB Hessian (modified Seminario, calibrated to DFT);
   - torsions from the DFT torsion library.
5. **MM model** (`dyff_edits.dyff_model`):
   - DYFF types from the `dyff-1.0-eh` set (pDynamo's dyff-1.0 plus EasyHybrid's added types);
   - bonded terms replaced by the library/xTB terms;
   - charges, with the total made an integer;
   - small LJ on polar hydrogens, TIP3P oxygen LJ on waters, guanidinium correction.
6. **Non-bonded model** (`NBModelCutOff`, no periodic box) and a single-point energy. With `--minimize`, an LBFGS minimization follows.
7. **Save** to `output/<name>_dyff_eh.*`:
   - `.pkl`: the system;
   - `.pdb`: the prepared structure;
   - `.json`: a short report.

## The example inputs

| file | what it shows |
|---|---|
| `inputs/theophylline.xyz` | coordinates only (PDB CCD TEP, ideal geometry, all H): bonds and bond orders perceived from the geometry |
| `inputs/valacyclovir.mol2` | bonds and orders from the file; flexible molecule (11 of 12 rotatable bonds found in the torsion library) |
| `inputs/4oas_raw.pdb` | raw crystal structure (PLINDER 4oas, chain A + ligand 2SW, no H): protonation, ~850 H added, 3 truncated residues completed, 93 of 94 residues from the library, the ligand by xTB |

Run on 2026-10-15 (output/ has the results). For the raw PDB, three residues truncated in the crystal (Gln18, Glu69, Lys70) are completed, and their side-chain rotamers are chosen away from the other atoms. Energy: −3 163 kcal/mol as prepared, −6 929 kcal/mol after `--minimize`. Without the rotamer choice the completed Gln18 landed on Arg97, at +65 000 kcal/mol. Timings with xTB 6.7.1 and 4 threads: 3 s, 15 s and 81 s.

## Using the .pkl

- **EasyHybrid:** *File > Import System*, file type "pDynamo files (*.pkl, *.yaml)".
- **A pDynamo script:** the MM model class lives in EasyHybrid (`src/pdynamo/dyff_edits.py`), so EasyHybrid's `src/` must be in `sys.path` before the file is read:

```python
import sys
sys.path.insert(0, "/path/to/EasyHybrid3/src")
from pBabel import ImportSystem
system = ImportSystem("output/4oas_raw_dyff_eh.pkl")
system.Summary()
```

## Running it in NAMD

`../dyff_eh_namd_export/convert_to_namd.py` turns these `.pkl` files into NAMD input in CHARMM format (`.psf`, `.prm`, `.pdb`, `.crd`). It checks the conversion against pDynamo by single-point energies and a NAMD minimization.

## Requirements

- pDynamo3 with EasyHybrid's modifications, and numpy.
- xTB, for atoms outside the library. Put it in the `PATH` or under `~/programs/xtb-*/bin/xtb`. The CM5 step also accepts `$XTBHOME/bin/xtb`, but the Hessian step (`seminario.py`) does not.
- RDKit and Open Babel are not needed.
