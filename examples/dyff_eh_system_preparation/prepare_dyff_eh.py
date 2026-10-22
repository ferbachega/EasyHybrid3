#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  prepare_dyff_eh.py
#
#  Step by step: a structure file (XYZ, MOL2 or PDB) -> a pDynamo system with the
#  EasyHybrid DYFF force field (DYFF-1.0-eh + residue library + xTB-derived terms)
#  -> a .pkl file that EasyHybrid and pDynamo scripts can open.
#
#  This is what EasyHybrid does in the GUI with  System > Prepare DYFF System
#  (followed by "DYFF + residue library"), written out as a plain script.
#
#  Usage (from the EasyHybrid3 root or from this folder):
#
#      python3 prepare_dyff_eh.py inputs/theophylline.xyz
#      python3 prepare_dyff_eh.py inputs/valacyclovir.mol2
#      python3 prepare_dyff_eh.py inputs/4oas_raw.pdb --minimize
#      python3 prepare_dyff_eh.py                       # the three inputs above
#
#  Options:  --ph 7.4         pH used to protonate the structure
#            --no-solvent     CM5 charges in the gas phase (default: implicit water, ALPB)
#            --charges-only   outside the library: only CM5 charges (no xTB Hessian/torsions)
#            --minimize       short LBFGS minimization (hydrogens first, then everything)
#            --threads N      threads for xTB (default 4)
#            --out DIR        output folder (default: ./output)
#
#  Output (DIR/<name>_dyff_eh.*):
#      .pkl   the pDynamo System with the DYFF-1.0-eh MM model, charges and NB model
#      .pdb   the prepared structure (hydrogens, protonation) to look at
#      .json  a short report (atoms, charge, library matches, energies)
#
#  Requirements: pDynamo3 (with EasyHybrid's modifications), numpy, and xTB (in the
#  PATH or as ~/programs/xtb-*/bin/xtb) for the atoms that are not in the residue
#  library. RDKit and Open Babel are NOT needed.
#
#  Opening the .pkl outside EasyHybrid: the MM model class lives in EasyHybrid
#  (pdynamo/dyff_edits.py), so put EasyHybrid's src/ folder in sys.path before
#  ImportSystem("....pkl") -- see the end of this file.
#

import argparse
import contextlib
import io
import json
import os
import sys
import time

# -----------------------------------------------------------------------------
# . EasyHybrid modules: src/ and the graphics engine (only for imports)
# -----------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
EASYHYBRID = os.path.dirname(os.path.dirname(HERE))
for path in (os.path.join(EASYHYBRID, "src"), os.path.join(EASYHYBRID, "src", "graphics_engine", "src")):
    if path not in sys.path:
        sys.path.insert(0, path)

import pBabel                                                   # noqa: F401,E402 (registers the file formats)
from pBabel import ImportSystem, ExportSystem                   # noqa: E402
from pCore import Pickle, Selection                             # noqa: E402
from pMolecule.NBModel import NBModelCutOff                     # noqa: E402
from pSimulation import LBFGSMinimize_SystemGeometry            # noqa: E402

from pdynamo import dyff_prep as DP                             # noqa: E402  system preparation
from pdynamo import residue_library as RL                       # noqa: E402  residue library + xTB rest
from pdynamo.dyff_edits import dyff_model, soften_polar_hydrogens, water_lj, EH_DYFF_SET   # noqa: E402
from pdynamo.pDynamo2EasyHybrid.session import _apply_dyff_guanidinium_correction           # noqa: E402

KJ_TO_KCAL = 1.0 / 4.184
EXAMPLES = ("inputs/theophylline.xyz", "inputs/valacyclovir.mol2", "inputs/4oas_raw.pdb")
EXAMPLES = ("/home/fernando/programs/EasyHybrid3_dev/examples/dyff_prep_bachega/1BX4.pdb")

@contextlib.contextmanager
def quiet():
    """ pDynamo prints a lot (some of it straight to the terminal's file
        descriptor, e.g. the PDB reader summary); the steps below print their
        own summary. """
    sys.stdout.flush()
    saved = os.dup(1)
    with open(os.devnull, "w") as devnull:
        os.dup2(devnull.fileno(), 1)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                yield
        finally:
            sys.stdout.flush()
            os.dup2(saved, 1)
            os.close(saved)


def step(n, text):
    print("\n[{}] {}".format(n, text))


def prepare(path, ph=7.4, solvent="water", bonded=True, minimize=False, threads=4, outdir="output"):
    name = os.path.splitext(os.path.basename(path))[0]
    report = {"input": path, "force_field": "DYFF ({}) + EasyHybrid residue library".format(EH_DYFF_SET)}
    t0 = time.time()
    print("=" * 78)
    print(" {}".format(os.path.relpath(path)))
    print("=" * 78)

    # -------------------------------------------------------------------------
    # STEP 1. Read the file. XYZ: coordinates only; MOL2: atoms, bonds and bond
    #         orders; PDB: residues and atom names (bonds are perceived later).
    # -------------------------------------------------------------------------
    step(1, "Import the structure")
    with quiet():
        raw = ImportSystem(path, log=None)
    n_bonds = len(list(raw.connectivity.bonds)) if raw.connectivity is not None else 0
    n_h = sum(1 for a in raw.atoms if a.atomicNumber == 1)
    print("    {} atoms ({} hydrogens), {} bonds in the file".format(len(raw.atoms), n_h, n_bonds))

    # -------------------------------------------------------------------------
    # STEP 2. Chemistry of the structure (EasyHybrid "Prepare DYFF System"):
    #   - residues are recognized (pDynamo PDB component library); anything else
    #     (ligands, a molecule read from XYZ/MOL2) gets EasyHybrid's native
    #     chemistry perception: connectivity, bond orders, aromaticity and, when
    #     there are no hydrogens, hydrogen counts from the geometry;
    #   - protonation at the chosen pH (residues and ligand groups), hydrogens
    #     added, truncated residues completed, metal sites detected.
    # -------------------------------------------------------------------------
    step(2, "Prepare the system: chemistry, protonation at pH {:.1f}, hydrogens".format(ph))
    info = DP.analyze(raw)
    if info["covalent"]:
        # . covalently bound residues (e.g. a ligand bonded to a Cys) are not handled
        #   by the preparation: dropped here, as the GUI does by default
        print("    removing {} covalently bound residue(s)".format(len(info["covalent"])))
        DP.remove_residues(info, list(info["covalent"]))
    DP.apply_site_protonation(info)                     # . donors of bonded metal sites
    with quiet():
        system, prep = DP.build(info, label=name, complete_truncated=True, ph=ph)
    if prep["undefined"]:
        raise ValueError("{} atoms without coordinates after the preparation".format(len(prep["undefined"])))
    total_formal = sum(int(getattr(a, "formalCharge", 0) or 0) for a in system.atoms)
    print("    {} atoms after preparation, total formal charge {:+d}".format(len(system.atoms), total_formal))
    for w in prep["warnings"][:5]:
        print("    note:", w)
    for site in prep.get("metal_sites", []):
        print("    metal site: {} ({}, charge {:+d}, multiplicity {})".format(site.get("label"), site.get("model"),
                                                                         site.get("charge", 0), site.get("multiplicity", 1)))
    report.update(atoms=len(system.atoms), formal_charge=total_formal)

    # -------------------------------------------------------------------------
    # STEP 3. Residue library: residues matched by their bond graph get charges
    #   and bonded terms precomputed with xTB (amino acids, nucleotides, ...);
    #   waters get TIP3P charges and ions their formal charge. The atoms left
    #   ("rest": ligands, unusual residues) are parametrized in step 4.
    # -------------------------------------------------------------------------
    step(3, "Residue library")
    plan = RL.plan(system)
    s = plan["summary"]
    print("    {} of {} residues matched ({})".format(s["matched"], s["residues"],
          ", ".join("{} {}".format(n, t) for t, n in sorted(s["templates"].items())) or "none"))
    if s.get("waters") or s.get("ions"):
        print("    {} waters, {} ions".format(s.get("waters", 0), s.get("ions", 0)))
    print("    {} atoms outside the library{}".format(s["rest_atoms"],
          ": " + ", ".join(s["unmatched"][:8]) if s["unmatched"] else ""))
    report["library"] = {"matched": s["matched"], "residues": s["residues"], "rest_atoms": s["rest_atoms"]}

    # -------------------------------------------------------------------------
    # STEP 4. Parameters for the atoms outside the library, on capped subsystems:
    #   - CM5 charges from GFN2-xTB (implicit water, ALPB, by default);
    #   - bonds and angles from the xTB Hessian (modified Seminario, calibrated
    #     to DFT by chemical environment);
    #   - torsions from the DFT torsion library (by chemical environment; the
    #     DYFF rules where none is found).
    # -------------------------------------------------------------------------
    charges = dict(plan["charges"])
    edits = {k: dict(v) for k, v in plan["edits"].items()}
    if plan["rest_atoms"]:
        step(4, "xTB parameters for {} atoms outside the library ({})".format(
            len(plan["rest_atoms"]), "CM5 + Hessian + torsions" if bonded else "CM5 charges only"))
        t4 = time.time()
        rest = RL.rest_parameters(system, plan["rest_atoms"], solvent=solvent, hessian=bonded, torsions=bonded,
                                  method="gfn2", threads=threads, site_states=DP.site_groups(system, prep),
                                  say=lambda text: print("    ...", text))
        charges.update(rest["charges"])
        for key in ("bond_terms", "angle_terms", "dihedral_bonds"):
            edits.setdefault(key, {}).update(rest["edits"][key])
        if bonded:
            r = rest["report"]
            print("    {} bonds and {} angles from the Hessian; {} of {} rotatable bonds from the torsion library ({:.0f} s)".format(
                r["bonds"], r["angles"], r["torsions"], r["rotatable"], time.time() - t4))
    else:
        step(4, "xTB parameters: nothing outside the library")

    # -------------------------------------------------------------------------
    # STEP 5. The DYFF-1.0-eh MM model:
    #   - atom types: DYFF-1.0 patterns + the EasyHybrid extensions (dyff-1.0-eh);
    #     metal types from the preparation (type_overrides);
    #   - bonded terms: DYFF rules, replaced by the library / xTB terms (edits);
    #   - charges: library + CM5, the total made integral;
    #   - polar hydrogens get a small LJ type, waters the TIP3P oxygen LJ;
    #   - guanidinium groups (Arg) corrected.
    # -------------------------------------------------------------------------
    step(5, "MM model: DYFF ({}) with library/xTB terms and charges".format(EH_DYFF_SET))
    high_spin = {i for i, st in DP.site_groups(system, prep).items()
                 if system.atoms[i].atomicNumber in (25, 26, 27) and int(st.get("multiplicity", 1)) >= 4}
    with quiet():
        n_guanidinium = _apply_dyff_guanidinium_correction(system)
        system.DefineMMModel(dyff_model("dyff-1.0", edits, prep.get("type_overrides"), high_spin_atoms=high_spin), log=None)
    n = len(system.atoms)
    q = [float(charges.get(i, 0.0)) for i in range(n)]
    residual = round(sum(q)) - sum(q)                   # . spread the rounding error over all atoms
    q = [v + residual / n for v in q]
    for i, v in enumerate(q):
        system.mmState.charges[i] = v
    n_polar = soften_polar_hydrogens(system)
    n_water = water_lj(system)
    print("    total charge {:+.3f}; {} polar H with the small LJ type; {} waters; {} guanidinium group(s)".format(
        sum(q), n_polar, n_water, n_guanidinium))
    report.update(total_charge=round(sum(q), 6))

    # -------------------------------------------------------------------------
    # STEP 6. Non-bonded model and a single-point energy (no periodic box here:
    #   add one, solvate etc. in EasyHybrid or in your own script if needed).
    # -------------------------------------------------------------------------
    step(6, "Non-bonded model and energy")
    system.DefineNBModel(NBModelCutOff.WithDefaults())
    with quiet():
        energy = system.Energy(log=None)
    print("    E = {:.2f} kcal/mol".format(energy * KJ_TO_KCAL))
    report["energy_kcal"] = round(energy * KJ_TO_KCAL, 3)
    if minimize:
        # . hydrogens first (heavy atoms fixed: relieves the clashes of added H),
        #   then every atom
        with quiet():
            system.freeAtoms = Selection.FromIterable([a.index for a in system.atoms if a.atomicNumber == 1])
            LBFGSMinimize_SystemGeometry(system, log=None, maximumIterations=200, rmsGradientTolerance=1.0)
            system.freeAtoms = None
            LBFGSMinimize_SystemGeometry(system, log=None, maximumIterations=500, rmsGradientTolerance=0.4)
            energy = system.Energy(log=None)
        print("    after minimization: E = {:.2f} kcal/mol".format(energy * KJ_TO_KCAL))
        report["energy_minimized_kcal"] = round(energy * KJ_TO_KCAL, 3)

    # -------------------------------------------------------------------------
    # STEP 7. Save: .pkl (the whole system, ready for EasyHybrid / pDynamo),
    #   .pdb (to look at) and .json (this report).
    # -------------------------------------------------------------------------
    step(7, "Save")
    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, name + "_dyff_eh")
    Pickle(base + ".pkl", system)
    with quiet():
        ExportSystem(base + ".pdb", system, log=None)
    report["seconds"] = round(time.time() - t0, 1)
    json.dump(report, open(base + ".json", "w"), indent=1)
    shown = os.path.relpath(base)
    print("    {0}.pkl  {0}.pdb  {0}.json  ({1:.0f} s in total)".format(shown, time.time() - t0))
    return base + ".pkl"


def main():
    parser = argparse.ArgumentParser(description="XYZ / MOL2 / PDB -> pDynamo system with EasyHybrid DYFF (dyff-1.0-eh) -> .pkl")
    parser.add_argument("inputs", nargs="*", help="structure files (default: the three files in inputs/)")
    parser.add_argument("--ph", type=float, default=7.4)
    parser.add_argument("--no-solvent", action="store_true", help="CM5 charges in the gas phase")
    parser.add_argument("--charges-only", action="store_true", help="outside the library: CM5 charges, DYFF bonded rules")
    parser.add_argument("--minimize", action="store_true")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--out", default=os.path.join(HERE, "output"))
    args = parser.parse_args()
    inputs = args.inputs or [os.path.join(HERE, f) for f in EXAMPLES]
    for path in inputs:
        pkl = prepare(path, ph=args.ph, solvent=None if args.no_solvent else "water", bonded=not args.charges_only,
                      minimize=args.minimize, threads=args.threads, outdir=args.out)
        # . the .pkl can be read back as any pDynamo system (EasyHybrid: File > Import
        #   System, pDynamo files) -- a quick check:
        with quiet():
            again = ImportSystem(pkl, log=None)
            e = again.Energy(log=None)
        print("    check: {} read back, {} atoms, E = {:.2f} kcal/mol".format(os.path.basename(pkl), len(again.atoms), e * KJ_TO_KCAL))


if __name__ == "__main__":
    main()
