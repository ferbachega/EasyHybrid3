# Solvent library (Builder > Solvate)

One `.mol2` file per solvent. Copy a file here (or use **Add...** in the
Solvate window) and it appears in the solvent list.

Two kinds of file are accepted:

1. **A single molecule**: EasyHybrid builds a periodic template box at the
   solvent density the first time the solvent is used (rigid-body Monte
   Carlo overlap removal, Packmol-style) and caches it in `.cache/`. The box
   is not equilibrated: minimise and equilibrate (NPT) the solvated system.
2. **A pre-equilibrated periodic box**: many copies of the molecule, one
   SUBSTRUCTURE id per molecule (column 7 of the ATOM lines), every molecule
   with the same atoms in the same order, plus a `@<TRIPOS>CRYSIN` section
   with an orthorhombic cell (90 degree angles). Used as is (best quality).
   `water_tip3p.mol2` is such a file.

Metadata go in comment lines at the top of the file, before `@<TRIPOS>MOLECULE`:

```
# name: Methanol        label shown in the Solvate window
# residue: MOH          residue name of every solvent molecule (max. 4 characters)
# density: 0.792        g/cm3, required for single-molecule files
```

Atom names come from the file (made unique if they repeat); bond orders
come from the BOND section (`ar` bonds are Kekulized).
