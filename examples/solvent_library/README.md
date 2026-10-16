# EasyHybrid Solvent Library 0.1

Initial library of 30 common solvents.

## Contents

Each solvent directory contains:

- `<name>.xyz` — Cartesian coordinates in Å
- `<name>.pdb` — PDB/HETATM representation
- `solvent.json` — per-solvent metadata

The root `solvents.json` contains the complete library index.

## Important

The coordinates in this first package were generated with RDKit using ETKDGv3 followed by MMFF94 when parameters were available, otherwise UFF. They are intended as chemically reasonable starting geometries.

They are **not force-field parameter files**. In particular, the PDB/XYZ files do not define atomic partial charges, Lennard-Jones parameters, bonded parameters, or a complete NAMD/AMBER/CHARMM/OPLS/DYFF solvent model.

PubChem PUG REST can provide 3D compound records and SDF output. The individual `solvent.json` files contain the corresponding PubChem query URL for cross-checking or replacing the local geometries with PubChem conformers.

## Solvents

water, methanol, ethanol, 1-propanol, 2-propanol, 1-butanol, tert-butanol, acetone, methyl-ethyl-ketone, acetonitrile, propionitrile, dmso, dmf, nmp, thf, 2-methyltetrahydrofuran, diethyl-ether, diisopropyl-ether, dioxane, ethyl-acetate, butyl-acetate, dichloromethane, chloroform, carbon-tetrachloride, benzene, toluene, ethylbenzene, hexane, cyclohexane, heptane
