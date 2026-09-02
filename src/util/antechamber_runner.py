#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: antechamber_runner
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
"""
antechamber_runner
===================

Pure-Python (no GTK) backend for parametrizing a ligand/small molecule
with AmberTools' `antechamber` (atom typing + partial charges) and
`parmchk2` (missing bonded parameters -> .frcmod). Used by
gui/windows/setup/windows_and_dialogs/system_windows/
prepare_ligand_antechamber.py (the "Prepare Ligand (Antechamber)"
window). Kept free of any GTK imports so it can be tested/run
standalone -- mirrors util/tleap_runner.py's shape.

The .mol2 antechamber produces (atom types + charges filled in) and the
.frcmod parmchk2 produces are exactly what util/tleap_runner.py's
build_tleap_script() already knows how to load (loadmol2 / loadamberparams
-- see its "Additional files" handling), so a ligand prepared here can be
added straight into a tleap run.
"""

import os
import shutil
import subprocess


def write_ligand_pdb(vobject, residue_name, output_path, frame=-1):
    """ Writes a single-residue PDB from a vismol object, for use as
        antechamber input -- NOT pdynamo.pDynamo2EasyHybrid.helpers'
        export_special_PDB(), even though this window also has access
        to a loaded system's vobject the same way the tLeap window does.

        Reason, found by actually running this against real antechamber
        on a real loaded object: export_special_PDB() writes each atom's
        OWN atom.name into the PDB's fixed-width atom-name column
        (4 characters). Real vismol atom names can be longer than that
        -- e.g. duplicate-name disambiguation can produce "H_alt10" (7
        chars) -- which silently shifts every fixed-column field after
        it on that line. Antechamber's own PDB reader is strict about
        column position and rejected the resulting file outright:
        "Coordinates must be in Columns 31-38, 39-46 and 47-54".

        Sidesteps this by not reusing the object's atom/residue naming
        at all: antechamber only needs element identity, coordinates,
        and connectivity (which it (re)computes from geometry for a
        small isolated molecule like this) -- not the original atom/
        residue names. Names are regenerated here as a short, guaranteed-
        4-characters-or-less, per-element sequence (C1, C2, H1, H2, ...),
        and the residue name is the one the user typed in this window
        (truncated to 3 characters, the field's conventional width),
        not whatever the vobject's own residue happens to be called.
    """
    element_counts = {}
    lines = []
    for index, atom in enumerate(vobject.atoms.values()):
        element = (atom.symbol or "X").strip() or "X"
        element_counts[element] = element_counts.get(element, 0) + 1
        atom_name = (element + str(element_counts[element]))[:4]
        x, y, z = atom.coords(frame)
        lines.append(
            "{:<6s}{:5d} {:<4s} {:3s} {:1s}{:>4s}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:<4s}\n".format(
                "ATOM  ", index + 1, atom_name, residue_name[:3], "A", "1",
                x, y, z, 1.0, 1.0, element))

    with open(output_path, "w") as pdb_file:
        pdb_file.writelines(lines)


def find_antechamber_executable():
    """ Locates `antechamber`: PATH first (shutil.which), then
        $AMBERHOME/bin/antechamber. Returns the path, or None. """
    return _find_amber_executable("antechamber")


def find_parmchk2_executable():
    """ Locates `parmchk2`: PATH first (shutil.which), then
        $AMBERHOME/bin/parmchk2. Returns the path, or None. """
    return _find_amber_executable("parmchk2")


def _find_amber_executable(name):
    command = shutil.which(name)
    if command:
        return command
    amberhome = os.environ.get("AMBERHOME")
    if amberhome:
        candidate = os.path.join(amberhome, "bin", name)
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return None


# Charge methods exposed in the UI. Deliberately NOT exposing 'resp' or
# the Gaussian/MOPAC-backed options antechamber also supports -- those
# need external QM packages configured outside this tool's scope. All
# three of these work out of the box with a stock AmberTools install:
# 'bcc' (AM1-BCC, via the bundled `sqm`) is the standard choice for
# MD-ready ligand charges; 'gas' (Gasteiger) is instant/empirical, no QM
# at all; 'abcg2' is a newer bcc-family method, also via sqm.
CHARGE_METHODS = ["bcc", "gas", "abcg2"]

# Atom type sets -- same two options already offered for GAFF in the
# tLeap window (util/tleap_runner.py's list_leaprc_files('gaff')).
ATOM_TYPES = ["gaff", "gaff2"]


def run_antechamber(input_path, output_mol2_path, charge_method, net_charge,
                     multiplicity, residue_name, atom_type,
                     antechamber_command, workdir, input_format=None):
    """ Runs antechamber to assign atom types and partial charges.

        input_format -- antechamber's -fi format code (e.g. "mol2",
            "pdb", "sdf", "ac"). Guessed from input_path's extension if
            not given.

        Returns a dict:
            {'ok': bool, 'stdout': str, 'stderr': str, 'returncode': int,
             'mol2': path or None}

        'ok' requires a zero exit code AND the output .mol2 existing
        with non-zero size -- unlike tleap, antechamber has no single
        reliable "N errors" summary line to key off of; a failed sqm
        convergence or a valence antechamber can't assign typically
        shows up as a non-zero exit and/or simply no (or an empty)
        output file, so checking for a real, non-empty output file is
        the most robust signal available.
    """
    os.makedirs(workdir, exist_ok=True)

    if input_format is None:
        input_format = os.path.splitext(input_path)[1].lstrip(".").lower() or "pdb"

    output_name = os.path.basename(output_mol2_path)

    command = [
        antechamber_command,
        "-i", input_path,
        "-fi", input_format,
        "-o", output_name,
        "-fo", "mol2",
        "-c", charge_method,
        "-nc", str(net_charge),
        "-m", str(multiplicity),
        "-rn", residue_name,
        "-at", atom_type,
        "-pf", "y",
    ]

    result = subprocess.run(command, cwd=workdir, capture_output=True, text=True)

    output_path = os.path.join(workdir, output_name)
    ok = (result.returncode == 0) and os.path.isfile(output_path) and os.path.getsize(output_path) > 0

    return {
        "ok": ok,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "returncode": result.returncode,
        "mol2": output_path if ok else None,
    }


def run_parmchk2(mol2_path, output_frcmod_path, atom_type, parmchk2_command, workdir):
    """ Runs parmchk2 on an antechamber-typed .mol2 to generate the
        .frcmod covering any bonded parameters (bonds/angles/dihedrals)
        gaff/gaff2 don't already have a value for.

        Returns a dict:
            {'ok': bool, 'stdout': str, 'stderr': str, 'returncode': int,
             'frcmod': path or None}

        'ok' requires a zero exit code AND the output .frcmod existing
        (parmchk2 always writes one, even if empty -- see run_parmchk2's
        docstring note in the caller UI: check the file's own ATTENTION/
        "0.0" penalty-score comments to judge parameter quality, this
        function only confirms the run itself succeeded).
    """
    os.makedirs(workdir, exist_ok=True)
    output_name = os.path.basename(output_frcmod_path)

    command = [
        parmchk2_command,
        "-i", os.path.basename(mol2_path),
        "-f", "mol2",
        "-o", output_name,
        "-s", atom_type,
    ]

    result = subprocess.run(command, cwd=workdir, capture_output=True, text=True)

    output_path = os.path.join(workdir, output_name)
    ok = (result.returncode == 0) and os.path.isfile(output_path)

    return {
        "ok": ok,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "returncode": result.returncode,
        "frcmod": output_path if ok else None,
    }
