# EasyHybrid DYFF system (.pkl) → NAMD (.psf / .prm / .pdb / .crd)

`convert_to_namd.py` takes a system built with EasyHybrid's DYFF (for example by `../dyff_eh_system_preparation/prepare_dyff_eh.py`) and writes NAMD input in CHARMM format. It then checks the conversion against pDynamo.

```
python3 convert_to_namd.py                    # every .pkl in ../dyff_eh_system_preparation/output
python3 convert_to_namd.py --namd-minimize    # ... plus a NAMD minimization compared with pDynamo
python3 convert_to_namd.py my_system.pkl --namd /path/to/namd3
```

The conversion itself is `src/util/dyff_eh_to_namd.py`, which also works as a command-line tool:

```
python3 src/util/dyff_eh_to_namd.py system.pkl output/prefix [--namd /path/namd3] [--no-validate]
```

## Output (`output/<name>/`)

| file | content |
|---|---|
| `<name>.psf` | topology: atoms, types, charges, masses, bonds, angles, dihedrals, impropers (XPLOR/EXT format) |
| `<name>.prm` | CHARMM parameters (the atom types are listed in the comments at the top) |
| `<name>.pdb` | coordinates, as NAMD reads them |
| `<name>.crd` | CHARMM extended coordinates (for CHARMM; NAMD does not read CRD) |
| `<name>_minimize.namd` | a NAMD input to start from (1000 minimization steps) |
| `<name>_singlepoint.namd/.log` | the single point used for the check |
| `<name>_report.json` | the numbers printed by the script |

Every NAMD input made for these files needs two DYFF-specific lines:

```
vdwGeometricSigma  yes      # DYFF combines sigma (and epsilon) geometrically
1-4scaling         0.5      # DYFF 1-4 electrostatics (NAMD 3 name: oneFourScaling)
```

They are in `<name>_minimize.namd`. The `.prm` title also declares them, so EasyHybrid's *Run NAMD* window (CHARMM input) adds both lines by itself.

## What is exact and what is not

| DYFF term | in CHARMM/NAMD | |
|---|---|---|
| bonds, harmonic | harmonic bonds | exact |
| torsions (cosine series, incl. the DFT torsion library) | one CHARMM line K(1 + cos(nφ − δ)) per periodicity, several lines per dihedral | exact, up to a constant per dihedral |
| charges, electrostatics, 1-4 scaling 0.5 | PSF charges, `1-4scaling` | exact |
| Lennard-Jones, OPLS form, geometric σ and ε, own 1-4 parameters | Rmin/2, ε and 1-4 columns + `vdwGeometricSigma` | exact (no NBFIX needed) |
| angles (cosine series) | harmonic angle with the same minimum and curvature | approximate |
| out-of-plane (cosine series in ω) | harmonic improper in ψ: same minimum and curvature, via dω/dψ | approximate |

**Atom types.** EasyHybrid's DYFF has parameters per bond, angle and torsion (residue library, xTB), but CHARMM looks parameters up by atom-type tuple. So each atom gets a type defined by its LJ type and the parameters of its own terms. Equal residues share types, and any tuple left with two parameter sets is split. For the protein example this gives 713 types for 1 639 atoms. The backbone types differ between neighbours because the terms across a peptide bond depend on the neighbouring residue. Waters would share their types.

**Pyramidal centres.** For amines (N:Tet), the minimum of the DYFF out-of-plane term is not the planar geometry. The improper's ψ0 is taken from the DYFF minimum nearest to the structure, so it keeps the chirality. A harmonic improper has one well only, so N inversion is not described.

## Check: the three example systems (2026-10-15, NAMD 3.0.1; 4oas re-made after the planar-group fixes)

Single point without cutoff, kcal/mol. "Expected" is what NAMD must give for the written files: pDynamo for the exact terms, pDynamo minus the series constants for the dihedrals, and the harmonic terms evaluated in Python for angles and impropers. Every expected value is matched to ≤ 0.0004 kcal/mol, except electrostatics in the protein (−0.030 of −7 810), which comes from the Coulomb constants of the two programs. The pDynamo energies are taken at the PDB-rounded coordinates.

| system | angles: pDynamo → harmonic | out-of-plane: pDynamo → harmonic |
|---|---|---|
| theophylline (21 atoms) | 1.52 → 1.51 | 0.0001 → 0.0001 |
| valacyclovir (43 atoms) | 7.51 → 7.49 | 0.31 → 0.33 |
| 4oas, protein + ligand (1 639 atoms) | 258.0 → 256.8 | 54.4 → 58.6 |

NAMD minimization (1000 steps) compared with the structure in the `.pkl`:

| system | RMSD all / heavy atoms (Å) | largest displacement (Å) |
|---|---|---|
| theophylline | 0.001 / 0.001 | 0.00 |
| valacyclovir | 0.005 / 0.003 | 0.02 |
| 4oas | 0.21 / 0.18 | 1.9 |

For 4oas, the DYFF energy computed by pDynamo at the NAMD minimum (−8 835 kcal/mol) is lower than at the `.pkl` structure (−8 782 kcal/mol). That structure was minimized with a cutoff and few steps, so it is not a full minimum either.

## Limits

- Periodic systems: the PSF/PRM conversion is the same, and the minimization input gets the cell, PME and a 12 Å cutoff. The single-point check is only done without a box.
- Only the four DYFF term types are written: harmonic bonds, cosine angles, cosine dihedrals and cosine out-of-plane. A system with other MM terms (restraints, QC region) stops with a message.
- `dyff_eh_to_namd.py` needs EasyHybrid's `src/` in `sys.path` to read the `.pkl`, because the MM model class lives in `pdynamo/dyff_edits.py`. The script sets the path itself.
