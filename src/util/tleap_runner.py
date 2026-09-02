#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: tleap_runner
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
"""
tleap_runner
============

Pure-Python (no GTK) backend for preparing an AMBER system with tLeap
(AmberTools): locating the `tleap` executable, discovering which force
fields/water models are actually available in the local AmberTools
install, building a tleap input script, and running it as a subprocess.

Used by gui/windows/setup/windows_and_dialogs/system_windows/
prepare_amber_system.py (the "Prepare AMBER System" window). Kept free
of any GTK/pDynamo-session imports so it can be tested/run standalone.
"""

import os
import glob
import shutil
import subprocess
import re


def find_tleap_executable():
    """ Locates the `tleap` executable: PATH first (shutil.which), then
        $AMBERHOME/bin/tleap. Returns the path, or None if not found.
    """
    command = shutil.which("tleap")
    if command:
        return command

    amberhome = os.environ.get("AMBERHOME")
    if amberhome:
        candidate = os.path.join(amberhome, "bin", "tleap")
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate

    return None


def _leap_cmd_dir():
    """ $AMBERHOME/dat/leap/cmd -- where the leaprc.* files that
        list_leaprc_files() scans actually live. None if $AMBERHOME
        isn't set or doesn't look like a real AmberTools install.
    """
    amberhome = os.environ.get("AMBERHOME")
    if not amberhome:
        return None
    path = os.path.join(amberhome, "dat", "leap", "cmd")
    return path if os.path.isdir(path) else None


def list_leaprc_files(prefix):
    """ Scans $AMBERHOME/dat/leap/cmd/leaprc.<prefix>* and returns the
        part of each filename after "leaprc." (e.g. "protein.ff14SB",
        "water.tip3p", "gaff2"), sorted. Scanned rather than hardcoded
        so the list always matches whatever AmberTools version is
        actually installed, instead of going stale. Empty list if
        $AMBERHOME isn't set/found.
    """
    cmd_dir = _leap_cmd_dir()
    if cmd_dir is None:
        return []

    pattern = os.path.join(cmd_dir, "leaprc." + prefix + "*")
    names = []
    for path in sorted(glob.glob(pattern)):
        names.append(os.path.basename(path)[len("leaprc."):])
    return names


# Standard AMBER atomic-ion residue names (atomic_ions.lib, loaded by
# every leaprc.water.* file) offered in the "Neutralize/Ionize" UI.
# Not scanned -- this short, well-known set covers the common case;
# the "Additional files" list covers anything more exotic.
COMMON_CATIONS = ["Na+", "K+", "Mg2+", "Ca2+"]
COMMON_ANIONS  = ["Cl-"]


def build_tleap_script(pdb_path, protein_ff=None, glycam_ff=None, gaff_ff=None,
                        extra_files=None, bonds=None, solvate=None, ions=None,
                        output_basename="prepared_system"):
    """ Builds the text of a tleap input script.

        pdb_path -- structure to load (loadpdb).
        bonds -- list of ((resnum1, atomname1), (resnum2, atomname2))
            pairs, emitted as `bond mol.<resnum1>.<atomname1>
            mol.<resnum2>.<atomname2>` right after loadpdb (tleap's way
            of forcing a covalent bond it wouldn't otherwise detect by
            distance -- disulfide bridges, metal coordination, linking
            a ligand to the protein, ...). resnum/atomname must match
            exactly what ends up in pdb_path's residue-number/atom-name
            columns -- see prepare_amber_system.py's
            on_button_add_bond_from_picking_clicked, which reads them
            straight off a picked vismol Atom's `.residue.index`/`.name`
            (verified to be the very same values export_special_PDB()
            writes for that atom, since both come from the same
            underlying object).
        protein_ff, glycam_ff, gaff_ff -- leaprc names (as returned by
            list_leaprc_files(), e.g. "protein.ff14SB", "gaff2") to
            `source`, or None to skip.
        extra_files -- list of extra parameter file paths (.frcmod
            loaded via loadamberparams, .lib/.off via loadoff, .mol2
            via loadmol2 -- picked by file extension). A .mol2 is
            loaded into a tleap variable named after the file's own
            basename (e.g. "NUU.mol2" -> `NUU = loadmol2"..."`)
            instead of a bare `loadmol2 "..."` call, because tleap's
            residue-matching step (the "Matching PDB residue names to
            LEaP variables" phase of loadpdb) looks up a global
            variable with the SAME NAME as the PDB residue -- with no
            assignment, the unit is loaded but never becomes reachable
            under that name, so loadpdb can't substitute it for the
            matching residue in the structure. The basename is used
            (not the residue name embedded inside the mol2 file
            itself) because that's what the residue name in the PDB
            actually is, for both sources of .mol2 files this window
            accepts: prepare_ligand_antechamber.py always names its
            output "<residue_name>.mol2" (the same residue_name it
            also writes into the ligand PDB via write_ligand_pdb()),
            and a user manually adding a .mol2 via "Additional Files"
            is expected to follow the same AMBER convention (name the
            file after the 3-4 letter residue code used in the PDB).
        solvate -- None, or a dict {'water_model': <leaprc.water.*
            name>, 'box_type': 'box'|'oct', 'buffer': float}. The
            solvent BOX TEMPLATE is always TIP3PBOX regardless of the
            chosen water_model -- this is standard AMBER/tleap
            practice (verified against a real tleap run on this
            machine): the box template only supplies the initial
            packing lattice, the actual water residue/parameters come
            from whichever leaprc.water.* was sourced above, so there
            is no need to match box-variable names per water model
            (which, in solvents.lib, aren't even 1:1 with model names
            -- e.g. spc and spce both reuse SPCBOX).
        ions -- None, or a dict {'mode': 'neutralize'|'add',
            'cation': str, 'n_cation': int, 'anion': str, 'n_anion':
            int}. 'neutralize' emits two SEPARATE `addions unit ion 0`
            commands, one per ion -- verified on this machine that
            tleap's addIons rejects a second ion when the count is 0
            ("'0' is not allowed as the value for the second ion");
            each call independently neutralizes with whichever ion is
            actually needed for the unit's current net charge sign
            (a no-op if that ion isn't needed). 'add' issues a single
            combined `addions unit cation n_cation anion n_anion` call
            (valid there, since counts are non-zero).
        output_basename -- saveamberparm writes "<output_basename>.top"/
            ".crd" in the tleap run's working directory. ".top"/".crd"
            rather than AMBER's own usual "*.prmtop"/"*.inpcrd"
            convention -- verified on this machine that pDynamo3's own
            format auto-detection (pBabel.ExportImport, used by
            session.py's load_a_new_pDynamo_system_from_dict() to import
            the result back into EasyHybrid) only registers the
            extensions "top"/"TOP" for Amber topologies and "crd"/"CRD"
            for Amber coordinates (see AmberTopologyFileReader.py /
            AmberCrdFileReader.py's own _Importer.AddHandler(...) calls)
            -- "prmtop"/"inpcrd" raise "Unrecognized format" there, even
            though the file *contents* are the standard, valid AMBER
            format either way.

        Returns the script text (str). Does not touch disk -- pass the
        result to run_tleap().
    """
    lines = []

    if protein_ff:
        lines.append("source leaprc.{}".format(protein_ff))
    if gaff_ff:
        lines.append("source leaprc.{}".format(gaff_ff))
    if glycam_ff:
        lines.append("source leaprc.{}".format(glycam_ff))
    if solvate:
        lines.append("source leaprc.water.{}".format(solvate["water_model"]))

    for extra_file in (extra_files or []):
        ext = os.path.splitext(extra_file)[1].lower()
        if ext == ".frcmod":
            lines.append('loadamberparams "{}"'.format(extra_file))
        elif ext in (".lib", ".off"):
            lines.append('loadoff "{}"'.format(extra_file))
        elif ext == ".mol2":
            varname = os.path.splitext(os.path.basename(extra_file))[0]
            lines.append('{} = loadmol2 "{}"'.format(varname, extra_file))
        # Unrecognized extensions are silently skipped -- the "Additional
        # files" UI only lets the user add these three types (see
        # prepare_amber_system.py's on_button_add_file_clicked).

    lines.append('mol = loadpdb "{}"'.format(pdb_path))

    for (resnum1, atomname1), (resnum2, atomname2) in (bonds or []):
        lines.append("bond mol.{}.{} mol.{}.{}".format(resnum1, atomname1, resnum2, atomname2))

    if solvate:
        command = "solvateoct" if solvate.get("box_type") == "oct" else "solvatebox"
        lines.append("{} mol TIP3PBOX {}".format(command, solvate.get("buffer", 10.0)))

    if ions:
        mode = ions.get("mode")
        cation = ions.get("cation") or "Na+"
        anion  = ions.get("anion")  or "Cl-"
        if mode == "neutralize":
            # Two SEPARATE calls -- see the docstring above for why a
            # single "addions mol cation 0 anion 0" call is invalid.
            lines.append("addions mol {} 0".format(cation))
            lines.append("addions mol {} 0".format(anion))
        elif mode == "add":
            n_cation = ions.get("n_cation", 0)
            n_anion  = ions.get("n_anion", 0)
            lines.append("addions mol {} {} {} {}".format(cation, n_cation, anion, n_anion))

    lines.append("saveamberparm mol {0}.top {0}.crd".format(output_basename))
    lines.append("quit")

    return "\n".join(lines) + "\n"


# Residue names tleap's own leaprc.water.* files treat as water (see e.g.
# "HOH = TP3" / "WAT = TP3" in leaprc.water.tip3p) -- used by
# insert_ter_records() below to spot the solute/solvent boundary.
_SOLVENT_RESIDUE_NAMES = {"HOH", "WAT", "TIP3", "TIP", "T3P", "T4P", "SPC"}


def insert_ter_records(pdb_path):
    """ Inserts PDB "TER" records at chain-polymer boundaries in-place.

        export_special_PDB() (used to write the structure tleap loads)
        writes every atom as a plain ATOM record with no TER records at
        all, even between separate chains/molecules -- fine for
        visualization, but tleap's loadpdb relies on TER to know where
        one polymer unit ends, and *needs* that boundary to apply a
        residue's C-terminal variant (the OXT atom). Verified on this
        machine: without this, tleap failed on a real PDB (1UBQ, protein
        immediately followed on paper by crystallographic waters with no
        separator) with "FATAL: Atom .R<GLY 76>.A<OXT 8> does not have a
        type." -- the last protein residue was silently treated as a
        middle-of-chain residue.

        Only handles the two boundary kinds actually observed/common:
        a chain ID change, and the transition into a solvent residue
        (HOH/WAT/...). This is a pragmatic default for the common
        "protein (+ crystallographic waters)" case, not a general PDB
        sanitizer -- more exotic inputs (multiple ligands packed into
        one chain with no separator, etc.) may still need the user to
        pre-clean the structure (e.g. with AmberTools' own pdb4amber)
        before using this window.
    """
    with open(pdb_path, "r") as pdb_file:
        lines = pdb_file.readlines()

    out_lines = []
    previous_chain = None
    previous_resname = None
    previous_atom_line = None

    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            chain = line[21:22]
            resname = line[17:20].strip()
            is_solvent = resname in _SOLVENT_RESIDUE_NAMES
            was_solvent = previous_resname in _SOLVENT_RESIDUE_NAMES

            boundary = (
                previous_chain is not None
                and (chain != previous_chain or (is_solvent and not was_solvent))
            )
            if boundary and previous_atom_line is not None:
                out_lines.append("TER\n")

            previous_chain = chain
            previous_resname = resname
            previous_atom_line = line

        out_lines.append(line)

    if previous_atom_line is not None:
        out_lines.append("TER\n")

    with open(pdb_path, "w") as pdb_file:
        pdb_file.writelines(out_lines)


_ERRORS_RE = re.compile(r"Errors\s*=\s*(\d+)")


def run_tleap(script_text, workdir, tleap_command, output_basename="prepared_system"):
    """ Writes `script_text` to "<workdir>/tleap.in" and runs
        `tleap_command -f tleap.in` there.

        Returns a dict:
            {'ok': bool, 'stdout': str, 'stderr': str, 'returncode': int,
             'prmtop': path or None, 'inpcrd': path or None}
        ('prmtop'/'inpcrd' name what these files ARE -- an AMBER
        topology and coordinate set -- the actual files tleap writes
        are named "<output_basename>.top"/".crd", see
        build_tleap_script()'s docstring for why.)

        'ok' requires BOTH a zero exit code AND "Errors = 0" in tleap's
        own summary line -- tleap frequently exits 0 even when a
        parametrization step failed (e.g. a residue tleap doesn't
        recognize), so the exit code alone isn't trustworthy; the
        "Exiting LEaP: Errors = N; Warnings = ...; Notes = ..." line it
        always prints is the real signal.
    """
    os.makedirs(workdir, exist_ok=True)
    script_path = os.path.join(workdir, "tleap.in")
    with open(script_path, "w") as script_file:
        script_file.write(script_text)

    result = subprocess.run(
        [tleap_command, "-f", "tleap.in"],
        cwd=workdir,
        capture_output=True,
        text=True,
    )

    errors_match = _ERRORS_RE.search(result.stdout)
    n_errors = int(errors_match.group(1)) if errors_match else None
    ok = (result.returncode == 0) and (n_errors == 0)

    prmtop_path = os.path.join(workdir, output_basename + ".top")
    inpcrd_path = os.path.join(workdir, output_basename + ".crd")
    if not (os.path.isfile(prmtop_path) and os.path.isfile(inpcrd_path)):
        ok = False
        prmtop_path = None
        inpcrd_path = None

    return {
        "ok": ok,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "returncode": result.returncode,
        "prmtop": prmtop_path,
        "inpcrd": inpcrd_path,
    }
