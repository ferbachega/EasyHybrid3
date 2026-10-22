#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  convert_to_namd.py
#
#  Step by step: an EasyHybrid DYFF system (.pkl, e.g. from
#  ../dyff_eh_system_preparation) -> NAMD input in CHARMM format:
#
#      <name>.psf   topology (atoms, charges, bonds, angles, dihedrals, impropers)
#      <name>.prm   parameters (DYFF -> CHARMM functional forms)
#      <name>.pdb   coordinates (what NAMD reads)
#      <name>.crd   CHARMM coordinates (for CHARMM itself; NAMD does not read it)
#      <name>_minimize.namd    a NAMD input to start from
#
#  then checks the conversion: single-point energy in NAMD vs pDynamo, component
#  by component, and (--namd-minimize) a NAMD minimization compared with pDynamo.
#
#  Usage (from the EasyHybrid3 root or from this folder):
#
#      python3 convert_to_namd.py                                   # every .pkl of ../dyff_eh_system_preparation/output
#      python3 convert_to_namd.py ../dyff_eh_system_preparation/output/4oas_raw_dyff_eh.pkl --namd-minimize
#      python3 convert_to_namd.py my_system.pkl --namd /path/to/namd3
#
#  The conversion itself is src/util/dyff_eh_to_namd.py (also a command-line tool):
#  bonds, dihedrals, charges and Lennard-Jones are exact; DYFF angle and
#  out-of-plane terms (cosine series) become harmonic angles and impropers with the
#  same minimum and curvature -- see the README.
#
#  Any NAMD input written for these files needs two DYFF-specific lines (they are in
#  <name>_minimize.namd, and EasyHybrid's "Run NAMD" window adds them by itself):
#
#      vdwGeometricSigma  yes
#      1-4scaling         0.5
#

import argparse
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EASYHYBRID = os.path.dirname(os.path.dirname(HERE))
for path in (os.path.join(EASYHYBRID, "src"), os.path.join(EASYHYBRID, "src", "graphics_engine", "src")):
    if path not in sys.path:
        sys.path.insert(0, path)

import numpy as np                                              # noqa: E402
from util import dyff_eh_to_namd as D                           # noqa: E402

DEFAULT_INPUTS = os.path.join(HERE, "..", "dyff_eh_system_preparation", "output", "*.pkl")


def step(n, text):
    print("\n[{}] {}".format(n, text))


def read_namd_coordinates(path):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in open(path) if l.startswith(("ATOM", "HETATM"))])


def process(pkl, outdir, namd=None, minimize=False):
    name = os.path.basename(pkl).replace("_dyff_eh.pkl", "").replace(".pkl", "")
    folder = os.path.join(outdir, name)
    prefix = os.path.join(folder, name)
    report = {"input": os.path.relpath(pkl)}
    print("=" * 78)
    print(" {}".format(os.path.relpath(pkl)))
    print("=" * 78)

    # -------------------------------------------------------------------------
    # STEP 1. Read the pDynamo system (DYFF MM model, charges, coordinates).
    # -------------------------------------------------------------------------
    step(1, "Read the system")
    system = D.load_system(pkl)
    print("    {} atoms, MM model {}".format(len(system.atoms), type(system.mmModel).__name__))

    # -------------------------------------------------------------------------
    # STEP 2. Convert. CHARMM looks parameters up by atom-type tuple, while an
    #   EasyHybrid DYFF system has parameters per bond / angle / torsion (residue
    #   library, xTB): atoms get the smallest set of types for which every type
    #   tuple has a single parameter set (equal residues share their types).
    # -------------------------------------------------------------------------
    step(2, "Convert to CHARMM format (PSF / PRM / PDB / CRD) + NAMD input")
    system, data = D.convert(system, prefix)
    s = D.summary(data)
    print("    {atoms} atoms, {atom_types} atom types; {bonds} bonds, {angles} angles, {dihedrals} dihedrals, "
          "{impropers} impropers".format(**s))
    print("    parameter lines: {} bonds, {} angles, {} dihedral tuples, {} impropers".format(
        s["bond_types"], s["angle_types"], s["dihedral_types"], s["improper_types"]))
    print("    files: {}.psf .prm .pdb .crd, {}_minimize.namd".format(os.path.relpath(prefix), os.path.basename(prefix)))
    report["conversion"] = {k: v for k, v in s.items() if k != "box"}

    if D.find_namd() is None and namd is None:
        print("\n    NAMD not found ($NAMD, namd3/namd2 in the PATH or ~/programs/NAMD*): no check.")
        json.dump(report, open(prefix + "_report.json", "w"), indent=1)
        return report

    # -------------------------------------------------------------------------
    # STEP 3. Single point, no cutoff: NAMD vs pDynamo per component. "expected"
    #   is what NAMD must give for the written files (exact terms: pDynamo; the
    #   dihedrals lose a constant per term; angles/impropers: the harmonic terms
    #   evaluated here) -- "diff" checks the files, pDynamo vs expected shows the
    #   size of the harmonic approximation.
    # -------------------------------------------------------------------------
    step(3, "Single-point energy: NAMD vs pDynamo (kcal/mol)")
    rows = D.validate(system, data, prefix, namd)
    print("    {:<6s} {:>13s} {:>13s} {:>13s} {:>9s}  {}".format("term", "pDynamo", "NAMD", "expected", "diff", ""))
    for comp, pd, nm, exp, note in rows:
        print("    {:<6s} {:13.4f} {:13.4f} {:13.4f} {:9.4f}  {}".format(comp, pd, nm, exp, nm - exp, note))
    report["single_point"] = [{"term": c, "pdynamo": p, "namd": n, "expected": e, "note": t} for c, p, n, e, t in rows]

    # -------------------------------------------------------------------------
    # STEP 4 (optional). Minimize in NAMD and compare the structure with the one
    #   in the .pkl, and the DYFF energy (pDynamo) at both.
    # -------------------------------------------------------------------------
    if minimize:
        step(4, "NAMD minimization (1000 steps) vs the pDynamo structure")
        D.run_namd(prefix + "_minimize.namd", namd, cpus=4)
        x1 = read_namd_coordinates(prefix + "_minimize.coor")
        x0 = np.array([[system.coordinates3[i, k] for k in range(3)] for i in range(len(system.atoms))])
        heavy = np.array([a.atomicNumber > 1 for a in system.atoms])
        d = np.sqrt(((x1 - x0) ** 2).sum(1))
        e0 = sum(D.pdynamo_energies(system).values())
        for i in range(len(x1)):
            for k in range(3):
                system.coordinates3[i, k] = x1[i, k]
        e1 = sum(D.pdynamo_energies(system).values())
        print("    RMSD to the .pkl structure: {:.3f} A (heavy atoms {:.3f} A, largest {:.2f} A)".format(
            float(np.sqrt((d ** 2).mean())), float(np.sqrt((d[heavy] ** 2).mean())), float(d.max())))
        print("    DYFF energy (pDynamo, no cutoff): {:.2f} at the .pkl structure, {:.2f} at the NAMD minimum".format(e0, e1))
        report["namd_minimization"] = {"rmsd": float(np.sqrt((d ** 2).mean())), "rmsd_heavy": float(np.sqrt((d[heavy] ** 2).mean())),
                                       "max_displacement": float(d.max()), "E_pdynamo_pkl": e0, "E_pdynamo_at_namd_minimum": e1}
    json.dump(report, open(prefix + "_report.json", "w"), indent=1)
    return report


def main():
    parser = argparse.ArgumentParser(description="EasyHybrid DYFF system (.pkl) -> NAMD (CHARMM-format PSF/PRM/PDB/CRD)")
    parser.add_argument("inputs", nargs="*", help=".pkl files (default: ../dyff_eh_system_preparation/output/*.pkl)")
    parser.add_argument("--namd", default=None, help="NAMD executable (default: $NAMD, PATH, ~/programs/NAMD*)")
    parser.add_argument("--namd-minimize", action="store_true", help="also minimize in NAMD and compare")
    parser.add_argument("--out", default=os.path.join(HERE, "output"))
    args = parser.parse_args()
    inputs = args.inputs or sorted(glob.glob(DEFAULT_INPUTS))
    if not inputs:
        sys.exit("no .pkl given and none in ../dyff_eh_system_preparation/output (run prepare_dyff_eh.py first)")
    for pkl in inputs:
        process(pkl, args.out, args.namd, args.namd_minimize)


if __name__ == "__main__":
    main()
