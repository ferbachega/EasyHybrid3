#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: namd_runner
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
"""
namd_runner
===========

Backend for preparing and running NAMD molecular dynamics jobs from
AMBER prmtop/inpcrd input. Pure Python and GTK-free (like
tleap_runner.py), except for one pDynamo3 import in
compute_cell_from_coordinates()/build_restraint_pdb() -- both need
real atom/residue names and coordinates, and pDynamo3's own
ImportSystem/ImportCoordinates3 (the exact pair session.py's
load_a_new_pDynamo_system_from_dict() already uses for system_type 0)
is the simplest way to get them straight from a prmtop/inpcrd pair,
with no dependency on the pair having ever been loaded into
EasyHybrid/vismol at all.

Cross-checked (keyword-by-keyword and, for CHARMM, physically re-run)
against real, working NAMD input: the AMBER+NAMD equilibration
protocol found on this machine at
~/programs/NAMD_protocol/AMBER_protocol/ (01.namd.tcl .. 09.namd.tcl,
README.TXT, run.sh, NaMD_Object3.py, pdb_to_pdb_fixed.py -- all read in
full), and NAMD's own bundled apoa1_gpu CHARMM example. The 9-step
AMBER protocol itself is no longer built by a dedicated function here
-- prepare_namd_run.py's "Custom Pipeline" treeview supersedes it (each
row calls build_namd_config() directly, chained the same way), with
that same 9-step sequence offered there as a loadable preset instead
(its one step this module can't reproduce exactly -- the real step 2's
single continuous 100-increment Tcl heating loop -- is approximated in
that preset as several separate NVT pipeline rows at increasing target
temperatures, chained the normal way).

Used by gui/windows/setup/windows_and_dialogs/system_windows/
prepare_namd_run.py (the "Run NAMD" window).
"""

import math
import os
import re
import shutil
import subprocess


# . Applied to every axis length before turning it into a PMEGridSize --
#   plain padding, not a physical assumption. build_namd_config() emits
#   ONE PME grid per pipeline run, reused unchanged by every step (the
#   .xsc from each step still carries that step's own real, evolved
#   cellBasisVector*, only the GRID stays fixed -- see its own docstring
#   for why). Under NPT the box can grow well past its starting size
#   before the run ends; a real case on this machine grew from ~44 Å to
#   ~136 Å across just 7 chained steps (a runaway barostat response to
#   an unrelated bad `cutoff`, since fixed -- but the grid should not
#   depend on every other parameter being perfect to stay sufficient)
#   and hit NAMD's own "PMEGridSizeX N is too small for cell length"
#   FATAL ERROR because the grid was sized to the ORIGINAL ~44 Å only.
#   25% headroom is a plain safety margin, not a guarantee -- it does
#   not replace fixing whatever is actually driving the box that far in
#   the first place, and a large enough blow-up will still exceed it.
_PME_GRID_MARGIN = 1.25


def _pme_grid_size(length):
    """ PMEGridSize for one axis: `length` (Å) padded by
        _PME_GRID_MARGIN, then rounded up -- see that constant's own
        comment for why the padding exists. PMEGridSpacing is always
        1.0 Å/point in every config this module builds, so grid points
        and padded Å are numerically the same thing here.
    """
    return int(math.ceil(length * _PME_GRID_MARGIN)) + 1


def find_namd_executable():
    """ Locates the `namd3` (preferred -- confirmed installed and
        working on this machine, CUDA-enabled build) or `namd2`
        executable via PATH (shutil.which). Unlike tleap, NAMD has no
        equivalent of $AMBERHOME to fall back on -- PATH, or a path the
        user types in by hand in the window (persisted the same way as
        vm_config.gl_parameters['tleap_command'], see
        prepare_amber_system.py), is the only automatic source.
        Returns the path, or None if not found.
    """
    for name in ("namd3", "namd2"):
        command = shutil.which(name)
        if command:
            return command
    return None


def write_amber_crd(positions, output_path, title="EasyHybrid"):
    """ Writes a plain-text AMBER coordinate file (the format
        `ambercoor` reads) from a plain list/array of (x, y, z)
        tuples -- NOT via pDynamo3's own AmberCrdFileWriter, which was
        confirmed (against a real prmtop/inpcrd pair on this machine)
        to be broken in this pDynamo3 install: it calls a `self.write`
        method TextFileWriter doesn't have (every other pBabel writer
        correctly uses `self.file.write`), so it raises
        AttributeError before writing anything.

        Format verified byte-for-byte against a real, working file
        (~/programs/NAMD_protocol/AMBER_protocol/sys.crd): a title
        line, then the atom count right-justified in 6 columns, then
        the coordinates 6 values per line (2 atoms/line, last line
        short if natom is odd), each formatted "%12.7f" with no
        separator between values.

        Used for the "coordinates from a loaded vismol object" input
        option (see prepare_namd_run.py's on_amber_coords_source_changed()
        and _validate_input_files()) -- writes out whatever the
        SELECTED object's CURRENT frame coordinates are, so NAMD can
        be pointed at it exactly like any other `ambercoor` file (the
        atom count/order still has to match the accompanying prmtop,
        same requirement as a restart file already has).
    """
    positions = list(positions)
    with open(output_path, "w") as f:
        f.write("{}\n".format(title))
        f.write("{:6d}\n".format(len(positions)))
        values = [component for xyz in positions for component in xyz]
        for i in range(0, len(values), 6):
            f.write("".join("{:12.7f}".format(v) for v in values[i:i + 6]))
            f.write("\n")


def ensure_amber_top_crd(prmtop_path, inpcrd_path, workdir):
    """ pDynamo3's own system_type==0 import path
        (session.py's load_a_new_pDynamo_system_from_dict(), used by
        prepare_namd_run.py to import the finished NAMD run's original
        topology back into EasyHybrid) calls ImportSystem()/
        ImportCoordinates3() WITHOUT an explicit `format=` kwarg --
        unlike _load_amber_system() above, it relies purely on the
        file's own extension, which only recognises ".top"/".crd" (see
        tleap_runner.py's build_tleap_script() docstring for the same
        quirk -- worked around THERE by always writing ".top"/".crd"
        outputs in the first place). This window, unlike the tleap
        one, must also accept prmtop/inpcrd files it did NOT generate
        itself (the whole point of good AMBER support -- see this
        module's own docstring), which by AMBER/NAMD convention are
        normally named ".prmtop"/".inpcrd" (confirmed with this
        machine's own real ~/programs/NAMD_protocol/AMBER_protocol/
        sys.prmtop). Rather than touching the shared session.py (used
        by several other callers) to add format= there too, this
        creates ".top"/".crd"-named SYMLINKS pointing at the real
        files inside `workdir`, and returns those paths for the caller
        to pass to load_a_new_pDynamo_system_from_dict() instead of
        the originals. Falls back to a plain copy if symlinking isn't
        possible (e.g. some Windows configurations without the
        privilege).

        Returns (top_path, crd_path).
    """
    os.makedirs(workdir, exist_ok=True)
    top_path = os.path.join(workdir, "_namd_import.top")
    crd_path = os.path.join(workdir, "_namd_import.crd")
    for target, link_path in ((prmtop_path, top_path), (inpcrd_path, crd_path)):
        if os.path.lexists(link_path):
            os.remove(link_path)
        try:
            os.symlink(os.path.abspath(target), link_path)
        except OSError:
            shutil.copyfile(target, link_path)
    return top_path, crd_path


# Same set as tleap_runner.py's _SOLVENT_RESIDUE_NAMES (kept as an
# independent copy rather than a cross-module import of a
# module-private name -- both modules treat this list as an internal
# implementation detail, not a shared public constant).
_SOLVENT_RESIDUE_NAMES = {"HOH", "WAT", "TIP3", "TIP", "T3P", "T4P", "SPC"}


def _iter_pdb_atoms(pdb_path):
    """ Yields (atom_name, res_name, res_seq, x, y, z) for every
        ATOM/HETATM record in a plain-text PDB, reading the fixed PDB
        columns directly (atom name 13-16, residue name 18-20, residue
        sequence number 23-26, x/y/z 31-38/39-46/47-54 -- 1-indexed
        per the PDB spec). Used by the CHARMM
        input helpers below (compute_cell_from_charmm(),
        build_restraint_pdb_from_charmm()) instead of routing through
        pDynamo3's own CHARMM PSF/parameter reader -- see this
        module's own docstring for why: pDynamo3's CHARMMParameterFileReader
        was confirmed (against real files from this machine's own
        ~/programs/NAMD_3.0.1_.../lib/namdcph/toppar/) to silently
        parse ".str" CHARMM "stream" parameter files (the modern,
        common way water/ions parameters are distributed) as EMPTY --
        0 bonds/angles despite non-zero section counts -- which then
        crashes ImportSystem() with a KeyError the moment the PSF
        references an atom-type pair (e.g. the TIP3 water "HT-HT"
        angle) that ended up missing. NAMD ITSELF parses these files
        fine (it has its own, separate, unaffected reader) -- only
        pDynamo3's Python-side reader has this gap. Since neither of
        these two functions actually needs resolved FORCE FIELD
        parameters at all (just atom names/residue names/coordinates),
        reading the coordinates PDB directly sidesteps the whole
        problem for THEM. The one place this module still can't avoid
        pDynamo3's CHARMM reader is importing the final result back
        into EasyHybrid afterwards (see prepare_namd_run.py's
        _import_result_back()) -- that one genuinely needs a real,
        fully-parameterised pDynamo System, so it can still hit this
        same pDynamo3 limitation if the user's parameter set includes
        a ".str" file; that failure is caught and reported clearly
        there rather than crashing, since NAMD's own run already
        succeeded by that point regardless.
    """
    with open(pdb_path, "r") as pdb_file:
        for line in pdb_file:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            res_seq = line[22:26].strip()
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            yield atom_name, res_name, res_seq, x, y, z


def compute_cell_from_charmm(coordinates_path):
    """ CHARMM-input equivalent of compute_cell_from_coordinates() --
        same bounding-box-based cell (see that function's docstring for
        the full reasoning, including why cellOrigin is deliberately
        max/2 per axis, not (min+max)/2). Reads the coordinates PDB
        directly (see _iter_pdb_atoms()) rather than going through
        pDynamo3's CHARMM PSF/parameter reader.
    """
    xs, ys, zs = [], [], []
    for _atom_name, _res_name, _res_seq, x, y, z in _iter_pdb_atoms(coordinates_path):
        xs.append(x)
        ys.append(y)
        zs.append(z)
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    z_min, z_max = min(zs), max(zs)

    return {
        "cell_basis_vector1": (x_max - x_min, 0.0, 0.0),
        "cell_basis_vector2": (0.0, y_max - y_min, 0.0),
        "cell_basis_vector3": (0.0, 0.0, z_max - z_min),
        "cell_origin": (x_max / 2.0, y_max / 2.0, z_max / 2.0),
        "pme_grid_size": (
            _pme_grid_size(x_max - x_min),
            _pme_grid_size(y_max - y_min),
            _pme_grid_size(z_max - z_min),
        ),
    }


def build_restraint_pdb_from_charmm(coordinates_path, output_path):
    """ CHARMM-input equivalent of build_restraint_pdb() -- same
        occupancy/B-factor convention (occupancy: 1.00 solute / 0.00
        water; B-factor: 1.00 backbone C/CA/N of a non-water residue /
        0.00 otherwise -- see that function's docstring for the full
        `conskcol` O-vs-B reasoning, which applies identically here).
        Reads atom/residue names straight off the coordinates PDB (see
        _iter_pdb_atoms()) instead of pDynamo3's CHARMM reader.
    """
    lines = []
    for index, (atom_name, res_name, res_seq, x, y, z) in enumerate(_iter_pdb_atoms(coordinates_path)):
        is_water = res_name in _SOLVENT_RESIDUE_NAMES
        occupancy = 0.00 if is_water else 1.00
        b_factor = 1.00 if (not is_water and atom_name in ("C", "CA", "N")) else 0.00

        line = "{:<6s}{:5d} {:<4s} {:<3s} {:1s}{:>4s}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format(
            "ATOM  ", index + 1, atom_name[:4], res_name[:3], "A", res_seq[:4],
            x, y, z, occupancy, b_factor,
        )
        lines.append(line)
    lines.append("END\n")

    with open(output_path, "w") as f:
        f.writelines(lines)


def _load_amber_system(prmtop_path, inpcrd_path):
    """ ImportSystem/ImportCoordinates3 on a prmtop/inpcrd pair --
        exactly the pair session.py's load_a_new_pDynamo_system_from_dict()
        already uses for system_type 0. Returns (system, sequence),
        sequence possibly None (mirrors session.py's own
        _get_sequence_from_pdynamo_system() fallback chain, simplified
        since a prmtop-derived system always carries its own .sequence).

        Passes format="top"/"crd" EXPLICITLY (pBabel's ImportSystem/
        ImportCoordinates3 both accept a `format` kwarg that bypasses
        their normal by-extension dispatch entirely -- confirmed by
        reading ExportImport.py's GetHandler()) rather than relying on
        the file's own extension. Needed for good AMBER support: this
        window's whole point is to accept prmtop/inpcrd as they
        actually come from AmberTools/tleap/the wild (".prmtop"/
        ".inpcrd", confirmed with a real file from this machine's own
        ~/programs/NAMD_protocol/AMBER_protocol/sys.prmtop) -- pDynamo3's
        pure by-extension dispatch only recognizes ".top"/".crd" (see
        tleap_runner.py's build_tleap_script() docstring for the same
        quirk, worked around there by always writing ".top"/".crd"
        instead -- not an option here, since this function must also
        accept prmtop/inpcrd files it did NOT generate itself).
    """
    from pBabel import ImportSystem, ImportCoordinates3
    system = ImportSystem(prmtop_path, format="top")
    system.coordinates3 = ImportCoordinates3(inpcrd_path, format="crd")
    sequence = getattr(system, "sequence", None)
    return system, sequence


def read_prmtop_box_dimensions(prmtop_path):
    """ Reads the periodic box tleap itself recorded in the prmtop's own
        "%FLAG BOX_DIMENSIONS" section (present whenever the system was
        built with `setBox` -- absent for a gas-phase/non-periodic
        prmtop, in which case this returns None) -- a plain text scan,
        no pDynamo3/ImportSystem needed, so it works from the prmtop
        ALONE, before any inpcrd has even been picked.

        Format is Fortran "%FORMAT(5E16.8)": beta angle, then box
        lengths a/b/c, e.g. "9.00000000E+01  6.55170960E+01
        6.25461450E+01  5.89104020E+01" (confirmed against a real file,
        ~/programs/NAMD_protocol/AMBER_protocol/sys.prmtop).

        Returns (beta, (a, b, c)), or None if the prmtop has no
        BOX_DIMENSIONS section at all.
    """
    with open(prmtop_path) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith('%FLAG BOX_DIMENSIONS'):
            values = lines[i + 2].split()
            beta, a, b, c = (float(v) for v in values[:4])
            return beta, (a, b, c)
    return None


def cell_from_prmtop_box_dimensions(prmtop_path):
    """ Builds a cell dict (same shape as compute_cell_from_coordinates()'s
        own return value) straight from read_prmtop_box_dimensions(),
        for prefilling the "Manually set periodic cell" fields the
        moment a prmtop is picked -- before an inpcrd has necessarily
        been chosen, and without pDynamo3's heavier ImportSystem.

        Returns None when the prmtop has no BOX_DIMENSIONS section, OR
        when beta isn't ~90 degrees (a non-orthogonal/truncated-
        octahedron box -- this whole window's cell model, like
        compute_cell_from_coordinates()'s, is axis-aligned rectangular
        only, so a triclinic box can't be represented this way and is
        deliberately left for the user to set by hand instead).

        cellOrigin here is box-length/2 per axis, NOT derived from any
        actual coordinate (none read here) -- an approximation good
        enough for a freshly-solvated, centered system to look right at
        a glance; compute_cell_from_coordinates() (run automatically at
        Run time whenever "Manually set periodic cell" is left
        unchecked) still does the coordinate-accurate version from the
        real inpcrd once one is available.
    """
    result = read_prmtop_box_dimensions(prmtop_path)
    if result is None:
        return None
    beta, (a, b, c) = result
    if abs(beta - 90.0) > 0.01:
        return None

    return {
        "cell_basis_vector1": (a, 0.0, 0.0),
        "cell_basis_vector2": (0.0, b, 0.0),
        "cell_basis_vector3": (0.0, 0.0, c),
        "cell_origin": (a / 2.0, b / 2.0, c / 2.0),
        "pme_grid_size": (_pme_grid_size(a), _pme_grid_size(b), _pme_grid_size(c)),
    }


def compute_cell_from_coordinates(prmtop_path, inpcrd_path):
    """ Ports get_cell() from the user's own pdb_to_pdb_fixed.py:
        a rectangular periodic cell (cellBasisVector1/2/3 + cellOrigin)
        sized to the bounding box of the system's actual coordinates,
        plus a matching PME grid size.

        Reads coordinates directly from the prmtop/inpcrd pair via
        pDynamo3 (see _load_amber_system()) rather than parsing a PDB's
        fixed-width coordinate columns like the original script did --
        more accurate (no column truncation), and avoids requiring a
        PDB to exist at all.

        NOTE: cellOrigin is deliberately max/2 per axis (NOT
        (min+max)/2, which is what a symmetric box midpoint would
        normally use) -- this replicates the reference script's own
        get_cell() EXACTLY, since that is what actually produced the
        real, working templates on this machine. Not "fixed" here on
        purpose: this function's job is to reproduce proven-working
        behaviour, not to second-guess it.

        Only meaningful for a chain's FIRST step (or a standalone
        single run with no restart files) -- later steps get their
        cell from the previous step's own extendedSystem (.xsc) file
        instead (see prepare_namd_run.py's on_button_run_pipeline_clicked()),
        exactly like the reference protocol.
    """
    system, _ = _load_amber_system(prmtop_path, inpcrd_path)
    n = len(system.atoms.items)
    xs = [system.coordinates3[i][0] for i in range(n)]
    ys = [system.coordinates3[i][1] for i in range(n)]
    zs = [system.coordinates3[i][2] for i in range(n)]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    z_min, z_max = min(zs), max(zs)

    return {
        "cell_basis_vector1": (x_max - x_min, 0.0, 0.0),
        "cell_basis_vector2": (0.0, y_max - y_min, 0.0),
        "cell_basis_vector3": (0.0, 0.0, z_max - z_min),
        "cell_origin": (x_max / 2.0, y_max / 2.0, z_max / 2.0),
        "pme_grid_size": (
            _pme_grid_size(x_max - x_min),
            _pme_grid_size(y_max - y_min),
            _pme_grid_size(z_max - z_min),
        ),
    }


def build_restraint_pdb(prmtop_path, inpcrd_path, output_path):
    """ Ports export_PDB_const() from the user's own
        pdb_to_pdb_fixed.py: writes the companion restraint PDB
        (NAMD's `consref`/`conskfile`) that the whole 9-step
        equilibration protocol shares -- ONE file, with BOTH columns
        the reference protocol's steps switch between via `conskcol`:

          - occupancy (cols 55-60): 1.00 for every solute atom, 0.00
            for water -- read when a step sets `conskcol O` (steps
            1-4: "restrain the whole solute while water/ions settle").
          - B-factor  (cols 61-66): 1.00 ONLY for backbone atoms (atom
            name in {C, CA, N}) of non-water residues, 0.00 otherwise
            -- read when a step sets `conskcol B` (steps 5-8:
            "restrain just the backbone").

        build_namd_config()'s `restraint` dict picks which column a
        given step actually reads (see its `column` key); `constraints
        on`/`off` toggles whether this file is consulted at all (off
        for step 9). Cross-checked against all 9 real templates before
        writing this.

        Built directly from the prmtop's own atom/residue names via
        pDynamo3 (see _load_amber_system()), NOT from a vismol object
        -- so this works for ANY prmtop/inpcrd pair, including one
        never loaded into EasyHybrid at all.

        x/y/z columns hold the real atom coordinates (from inpcrd) --
        not functionally required by NAMD's constraint code here
        (which only reads the occupancy/B-factor columns off this
        file), but kept real/consistent with the actual structure so
        the file is sane to inspect by hand.
    """
    system, sequence = _load_amber_system(prmtop_path, inpcrd_path)

    lines = []
    for atom in system.atoms.items:
        xyz = system.coordinates3[atom.index]
        if sequence is not None:
            res_name, res_seq, _i_code = sequence.ParseLabel(atom.parent.label, fields=3)
        else:
            res_name, res_seq = "UNK", "1"
        res_name = res_name.strip()
        atom_name = atom.label.strip()

        is_water = res_name in _SOLVENT_RESIDUE_NAMES
        occupancy = 0.00 if is_water else 1.00
        b_factor = 1.00 if (not is_water and atom_name in ("C", "CA", "N")) else 0.00

        line = "{:<6s}{:5d} {:<4s} {:<3s} {:1s}{:>4s}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format(
            "ATOM  ", atom.index + 1, atom_name[:4], res_name[:3], "A", str(res_seq)[:4],
            xyz[0], xyz[1], xyz[2], occupancy, b_factor,
        )
        lines.append(line)
    lines.append("END\n")

    with open(output_path, "w") as f:
        f.writelines(lines)


def flag_pdb_occupancy(pdb_path, indexes, output_path):
    """ Copies `pdb_path` to `output_path`, rewriting the OCCUPANCY
        column (PDB spec columns 55-60, 1-indexed -- the SAME position
        build_restraint_pdb() and export_special_PDB() both already
        write to) of every ATOM/HETATM record: 1.00 if that record's
        0-based sequential position among ATOM/HETATM records is in
        `indexes`, 0.00 otherwise. This is the generic "flag a group"
        step NAMD's own SMDFile and fixedAtomsFile both expect
        (nonzero occupancy = "this atom is in the group").

        `indexes` are plain 0-based atom indexes, in the SAME order the
        source PDB's own records were written -- e.g. straight out of
        export_special_PDB(), whose record order is `vobject.atoms.
        keys()`, the same indexing EasyHybrid's own vismol selections
        and `p_session.psystem[id].e_selections` already use. So the
        typical caller is: export_special_PDB(vobject, -1, tmp_pdb),
        then flag_pdb_occupancy(tmp_pdb, e_selections['some name'],
        final_pdb) -- no pDynamo reload needed, unlike
        build_restraint_pdb() (which has no live vobject to draw a
        selection from in the first place, only a prmtop/inpcrd pair).
    """
    with open(pdb_path) as f:
        lines = f.readlines()

    flagged = set(indexes)
    out_lines = []
    position = 0
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            occupancy = 1.00 if position in flagged else 0.00
            line = line[:54] + "{:6.2f}".format(occupancy) + line[60:]
            position += 1
        out_lines.append(line)

    with open(output_path, "w") as f:
        f.writelines(out_lines)


def build_namd_config(
    output_basename,
    parmfile=None, ambercoor=None,
    charmm_psf=None, charmm_parameters=None, charmm_xplor=False, coordinates=None,
    bin_coordinates=None, bin_velocities=None, extended_system=None, cell=None,
    temperature=298.0, timestep=1.0, firsttimestep=0,
    cutoff=8.0, pairlistdist=12.0, watermodel="tip3",
    minimize_steps=0, run_steps=0, heating=None,
    ensemble="NPT",
    restraints=None,
    restartfreq=1000, dcdfreq=1000, outputenergies=1000, outputpressure=1000,
    langevin=True, langevin_damping=5.0, langevin_temp=None, langevin_hydrogen=False,
    flexible_cell=False,
    exclude="scaled1-4",
    switching=False, switchdist=None,
    one_four_scaling=0.833333333,
    readexclusions=True, scnb=2.0,
    rigidbonds="all", rigidtolerance=1e-08,
    vdw_force_switching=False,
    pbc=True, wrap_all=True, implicit_solvent=True,
    solvent_dielectric=78.5, ion_concentration=0.3, alpha_cutoff=15.0,
    sasa=False, surface_tension=0.005,
    fixed_atoms_file=None,
    smd=False, smd_file=None, smd_k=None, smd_k2=0.0,
    smd_vel=None, smd_dir=None, smd_output_freq=1,
):
    """ Builds the text of ONE NAMD config file (a .namd/.namd.tcl,
        NAMD's own Tcl-scripting input format). Every AMBER-mode
        keyword/value below was cross-checked against the real, working
        templates at ~/programs/NAMD_protocol/AMBER_protocol/*.namd.tcl
        -- this is a faithful port of a single step of that protocol,
        generalised into parameters instead of hardcoded per-file
        values. The CHARMM-mode keywords were cross-checked the same
        way against NAMD's own bundled real example at
        ~/programs/NAMD_3.0.1_.../apoa1_gpu/apoa1_gpuoff_{nve,npt}.namd.

        Exactly ONE of the two input modes must be used:

        AMBER mode -- parmfile, ambercoor: the AMBER prmtop/inpcrd pair
            (`amber yes` mode -- NAMD reads topology/starting
            coordinates from these two, same as the reference
            protocol's every step). AMBER-specific keywords this mode
            alone emits: `amber yes`, `readexclusions`, `scnb`
            (meaningless outside `amber yes` mode -- CHARMM mode omits
            both entirely), `oneFourScaling` (AMBER's own keyword name
            for the 1-4 electrostatic scaling factor -- see
            one_four_scaling below; CHARMM mode emits the same VALUE
            under ITS OWN keyword, `1-4scaling`, instead).

        CHARMM mode -- charmm_psf, charmm_parameters (a list of one or
            more parameter file paths -- CHARMM force fields routinely
            split parameters across several files, e.g. protein +
            water/ions), coordinates (the starting PDB). charmm_xplor
            selects which of the two parameter-file dialects NAMD
            should expect: False (default) -> `paraTypeCharmm on`
            (native CHARMM-format .prm/.par, the modern/common case);
            True -> `paraTypeCharmm off` (older X-PLOR-format .xplor
            files, matching the real apoa1_gpu example this was
            verified against). CHARMM-specific keywords this mode
            alone emits: `structure`/`parameters`/`paraTypeCharmm`
            instead of `amber`/`parmfile`.

        This function no longer hardcodes what NAMD's own AMBER-support
        docs (ug/node13.html) call the values that "characterize"
        AMBER-vs-CHARMM usage (exclude/switching/switchdist/1-4
        scaling/rigidBonds/rigidTolerance/readexclusions/scnb) -- every
        one of them is a real argument below with a default matching
        this function's PREVIOUS hardcoded AMBER behaviour, so a caller
        that passes none of them gets identical output to before; the
        "Run NAMD" window's Setup tab now exposes them as real widgets
        (a "Force Field" section) instead.
        bin_coordinates, bin_velocities, extended_system -- NAMD's own
            binary restart files (.coor/.vel/.xsc) from a PREVIOUS
            step, if this step continues one (steps 2-9 of the
            reference protocol always set all three; step 1 sets
            none). The topology/starting-structure keywords above stay
            set even when these are given -- NAMD still needs them --
            exactly like every real template (AMBER or CHARMM) does.
        bin_coordinates, bin_velocities, extended_system -- NAMD's own
            binary restart files (.coor/.vel/.xsc) from a PREVIOUS
            step, if this step continues one (steps 2-9 of the
            reference protocol always set all three; step 1 sets
            none). `ambercoor`/`parmfile` stay set even when these are
            given -- NAMD still needs them for topology -- exactly
            like every one of the real templates does.
        cell -- a dict from compute_cell_from_coordinates(), emitted as
            explicit cellBasisVector1/2/3 + cellOrigin lines. Ignored
            (and should be omitted) when extended_system is given --
            the .xsc already carries the (possibly NPT-fluctuated) cell
            from the previous step, exactly like steps 2-9 of the
            reference protocol (their cellBasisVector* lines are
            present in the template only as comments, explicitly
            noting the cell "já está escrito no xsc gerado pelo
            primeiro passo"). Also ignored entirely when `pbc` is False
            (implicit solvent has no periodic cell at all).
        pbc -- True (default) runs with a periodic cell + PME
            electrostatics, same as this function has always done.
            False drops the periodic cell entirely -- no PME, and NPT
            is meaningless (a barostat has no volume to control without
            a cell), so the langevinPiston block is skipped regardless
            of `ensemble` in that case -- and then `implicit_solvent`
            decides what replaces PME: NAMD's Generalized Born Implicit
            Solvent model (`GBIS`), or nothing at all (plain vacuum).
            Meaningless (and ignored) when `pbc` is True.
        wrap_all -- only meaningful when `pbc` is True. True (default)
            emits `wrapAll on`: every atom (not just water, unlike
            NAMD's separate `wrapWater`) gets wrapped back into the
            primary periodic image in the OUTPUT trajectory. Purely
            cosmetic (for visualization/analysis) -- it does not change
            the dynamics themselves, NAMD applies it only when writing
            coordinates out. Emitted on EVERY step that has a periodic
            cell, not just the first one that defines cellBasisVector*
            explicitly (a restart step reading its cell from a previous
            step's own .xsc still wants its own trajectory wrapped).
        implicit_solvent -- only meaningful when `pbc` is False. True
            (default) emits `GBIS on` plus the four parameters below.
            False emits neither GBIS nor PME nor a periodic cell -- a
            bare vacuum run with a plain nonbonded cutoff and nothing
            else, for when the user wants no solvent treatment at all
            rather than an implicit one.
        solvent_dielectric, ion_concentration, alpha_cutoff -- GBIS-only
            (meaningless when `pbc` is True): `solventDielectric`
            (default 78.5, water), `ionConcentration` (default 0.3,
            mol/L), `alphaCutoff` (default 15.0 Å -- GBIS's own cutoff
            for computing Born radii, distinct from and normally
            smaller than the main nonbonded `cutoff`).
        sasa, surface_tension -- GBIS-only: `SASA` (default off -- the
            optional nonpolar surface-area solvation term) and, only
            emitted when `sasa` is True, `surfaceTension` (default
            0.005 kcal/mol/Å², NAMD's own default).
        minimize_steps, run_steps -- `minimize $n` / `run $n`. Either
            may be 0 to skip (a step can be minimize-only, run-only, or
            both -- matches the reference: step 1 is minimize-only,
            steps 3/4/6/7/8/9 are run-only, step 5 is minimize-only,
            step 2 is minimize(2000, fixed) + the heating loop below).
        heating -- None, or a dict {'initial': float, 'final': float,
            'n_increments': int, 'steps_per_increment': int} -- emits
            the SAME Tcl `for` ramp loop as the reference protocol's
            step 2 (see 02.namd.tcl): steps langevinTemp from `initial`
            to `final` in `n_increments` equal increments, running
            `steps_per_increment` steps of dynamics at each rung.
        ensemble -- "NVT" (no langevinPiston -- constant volume, used
            by the reference protocol's step 2 ONLY, for the heating
            ramp) or "NPT" (langevinPiston on -- every other step, 1
            and 3 through 9, per the real templates -- verified line by
            line: only 02.namd.tcl has this block commented out).
        restraints -- None (no `constraints`/consref at all -- step 9
            of the reference protocol), or a dict {'pdb': path to a
            build_restraint_pdb() output, 'column': 'O' or 'B',
            'scaling': float} -- 'column' picks `conskcol` (steps 1-4
            use 'O': restrain the whole solute; steps 5-8 use 'B':
            backbone only -- see build_restraint_pdb()'s docstring for
            why one file serves both), 'scaling' is
            `constraintScaling` (the reference protocol's own schedule:
            100.0 for steps 1-3, 10.0 for steps 4-6, 1.0 for step 7,
            0.1 for step 8).
        langevin -- whether to run Langevin dynamics at all (the
            reference protocol always does -- every real template has
            "langevin yes"). False emits "langevin no" and omits
            langevinDamping/langevinTemp/langevinHydrogen entirely
            (meaningless without Langevin on).
        langevin_damping -- `langevinDamping` (the reference protocol's
            own value: 5, i.e. 5/ps).
        langevin_temp -- None (default) emits `langevinTemp
            $temperature` -- tied to the `temperature` argument above,
            exactly like every real template. A float instead emits a
            literal `langevinTemp <value>` -- lets the thermostat
            target a DIFFERENT temperature than the one used for the
            initial Maxwell-Boltzmann velocity assignment (an advanced,
            not-in-the-reference-protocol option, exposed for the
            window's "Advanced" tab).
        langevin_hydrogen -- `langevinHydrogen` (the reference
            protocol's own value: off -- doesn't couple the Langevin
            bath to hydrogens).
        flexible_cell -- `useFlexibleCell` (only emitted when
            ensemble == "NPT" -- meaningless otherwise). False
            (default, matches every real template) keeps the box's
            isotropic fluctuation: a/b/c scale together, ratio fixed.
            True lets each axis fluctuate independently -- needed for
            non-cubic/membrane-like systems, and worth trying when a
            run dies with NAMD's own "Periodic cell has become too
            small for original patch grid" FATAL ERROR (one of that
            message's own suggested fixes).
        exclude -- non-bonded exclusion policy: "none", "1-2", "1-3",
            "1-4" or "scaled1-4" (default -- AMBER's own recommended
            value, per NAMD's ug/node13.html; also correct for CHARMM).
        switching, switchdist -- `switching`. False (default -- AMBER's
            own convention, matches every real AMBER template on this
            machine) uses a hard cutoff with no switching function at
            all, and `switchdist` is not emitted (meaningless without
            switching on). True emits `switching on` and `switchdist`
            (CHARMM's own long-standing convention when it's the one
            enabled -- confirmed against the real apoa1 example: cutoff
            12./switchdist 10., i.e. switchdist = cutoff - 2); a
            `switchdist` of None then defaults to `cutoff - 2.0`.
        one_four_scaling -- the 1-4 electrostatic scaling factor,
            default 0.833333333 (AMBER's own fixed constant). Emitted
            under AMBER's own keyword `oneFourScaling` in AMBER mode,
            or CHARMM's own keyword `1-4scaling` in CHARMM mode (CHARMM
            force fields conventionally use 1.0 instead -- the window's
            Setup tab leaves the actual number up to the user either
            way, this is just which of the two on-disk keywords the
            SAME value ends up under).
        readexclusions, scnb -- AMBER-mode-only (silently omitted in
            CHARMM mode, where neither concept applies). Defaults yes/
            2.0 -- both match every real AMBER template on this
            machine.
        rigidbonds, rigidtolerance -- `rigidBonds`/`rigidTolerance`.
            Defaults "all"/1e-08, both match every real template
            (AMBER or CHARMM) on this machine.
        vdw_force_switching -- `vdwForceSwitching`. False (default --
            NAMD's own default, and every real template on this
            machine) uses the classic energy-switching VDW scheme;
            True switches on force-based switching instead (some newer
            force fields recommend it -- not something this codebase's
            own reference templates use, hence off by default here).
        fixed_atoms_file -- None (default, `fixedAtoms off`) or a path
            to a PDB with nonzero occupancy marking the atoms to hold
            completely fixed (`fixedAtoms on`, `fixedAtomsCol O`) --
            see flag_pdb_occupancy() for building one from an arbitrary
            set of atom indexes (e.g. a vismol/EasyHybrid selection).
            Independent of `restraints` above (harmonic, not fixed) and
            of SMD below (a fixed atom can still be an SMD anchor, or
            the two can be entirely separate atoms).
        smd, smd_file, smd_k, smd_k2, smd_vel, smd_dir, smd_output_freq
            -- Steered Molecular Dynamics (see NAMD's own ug/node48.html
            -- this window's "Run SMD" tool exists specifically to
            build these). `smd=False` (default) emits `SMD off` and
            ignores every other smd_* argument. `smd=True` requires
            `smd_file` (a PDB with nonzero occupancy marking the PULLED
            group -- same flag_pdb_occupancy() convention as
            fixed_atoms_file, built from a DIFFERENT selection), `smd_k`
            (`SMDk`, kcal/mol/Å², the pulling spring constant) and
            `smd_dir` (`SMDDir`, a (x, y, z) tuple -- NAMD normalizes it
            itself, no need to pre-normalize). `smd_vel` is `SMDVel` in
            NAMD's own native unit, Å/TIMESTEP (not Å/ps!) -- deliberately
            NOT converted here, so the caller (the "Run SMD" window)
            does that conversion itself from whatever more intuitive
            unit its own UI shows, using ITS OWN timestep value, rather
            than this function silently assuming one. `smd_k2` (`SMDk2`,
            default 0 -- NAMD's own default, no transverse restraint)
            and `smd_output_freq` (`SMDOutputFreq`, default 1 -- NAMD's
            own default) are optional.

        Returns the config text (str). Does not touch disk -- the
        caller writes it out and passes the path to run_namd().
    """
    is_charmm = charmm_psf is not None
    lines = []

    if is_charmm:
        lines.append('structure                       "{}"'.format(charmm_psf))
        lines.append('coordinates                     "{}"'.format(coordinates))
        lines.append("")
        lines.append("paraTypeCharmm                  {}".format("off" if charmm_xplor else "on"))
        for parameter_path in charmm_parameters or []:
            lines.append('parameters                      "{}"'.format(parameter_path))
    else:
        lines.append("amber                           yes")
        lines.append('parmfile                        "{}"'.format(parmfile))
        lines.append('ambercoor                        "{}"'.format(ambercoor))
    lines.append("")

    if bin_coordinates:
        lines.append('binCoordinates                  "{}"'.format(bin_coordinates))
    if bin_velocities:
        lines.append('binVelocities                   "{}"'.format(bin_velocities))
    if extended_system:
        lines.append('extendedSystem                  "{}"'.format(extended_system))
    if bin_coordinates or bin_velocities or extended_system:
        lines.append("")

    lines.append("set temperature                 {}".format(temperature))
    if not bin_velocities:
        # . Only meaningful when NOT restarting from a previous step's
        #   velocities -- matches every template: `temperature` (start
        #   from a Maxwell-Boltzmann distribution) is only ever paired
        #   with the FIRST step of a chain, later steps read binVelocities
        #   instead and never re-set `temperature`.
        lines.append("temperature                     $temperature")
    lines.append("firsttimestep                   {}".format(firsttimestep))
    if not is_charmm:
        # . `readexclusions`/`scnb` only mean anything in `amber yes`
        #   mode -- absent from the real CHARMM apoa1 example.
        lines.append("readexclusions                  {}".format("yes" if readexclusions else "no"))
        lines.append("scnb                            {}".format(scnb))
    lines.append("")

    lines.append("exclude                         {}".format(exclude))
    lines.append("switching                       {}".format("on" if switching else "off"))
    if switching:
        lines.append("switchdist                      {}".format(
            switchdist if switchdist is not None else cutoff - 2.0))
    if is_charmm:
        lines.append("1-4scaling                      {}".format(one_four_scaling))
    else:
        lines.append("oneFourScaling                  {}".format(one_four_scaling))
    lines.append("cutoff                          {}".format(cutoff))
    lines.append("watermodel                      {}".format(watermodel))
    lines.append("pairListDist                    {}".format(pairlistdist))
    if not is_charmm:
        # . Absent from the real CHARMM apoa1 example (defaults to off
        #   there) -- an AMBER-force-field-tail-correction convention,
        #   kept exactly as before for AMBER mode.
        lines.append("LJcorrection                    on")
    lines.append("rigidBonds                      {}".format(rigidbonds))
    lines.append("rigidTolerance                  {}".format(rigidtolerance))
    lines.append("vdwForceSwitching               {}".format("on" if vdw_force_switching else "off"))
    lines.append("rigidIterations                 100")
    lines.append("useSettle                       on")
    lines.append("fullElectFrequency              1")
    lines.append("nonBondedFreq                   1")
    lines.append("stepspercycle                   10")
    lines.append("timeStep                        {}".format(timestep))
    lines.append("")

    lines.append("langevin                        {}".format("yes" if langevin else "no"))
    if langevin:
        lines.append("langevinDamping                 {}".format(langevin_damping))
        lines.append("langevinTemp                    {}".format(
            "$temperature" if langevin_temp is None else langevin_temp))
        lines.append("langevinHydrogen                {}".format("on" if langevin_hydrogen else "off"))
    lines.append("")

    if pbc and ensemble == "NPT":
        # . A barostat controls PRESSURE by adjusting a periodic cell's
        #   own volume -- meaningless with no cell to adjust, so this
        #   whole block (and NPT as a concept) only applies when `pbc`
        #   is on, regardless of what `ensemble` was asked for.
        lines.append("useGroupPressure                yes")
        lines.append("useFlexibleCell                 {}".format("yes" if flexible_cell else "no"))
        lines.append("useConstantArea                 no")
        lines.append("langevinPiston                  yes")
        lines.append("langevinPistonTarget            1.01325")
        lines.append("langevinPistonPeriod            100.0")
        lines.append("langevinPistonDecay             50.0")
        lines.append("langevinPistonTemp              $temperature")
        lines.append("")

    if pbc:
        if cell and not extended_system:
            cbv1 = cell["cell_basis_vector1"]
            cbv2 = cell["cell_basis_vector2"]
            cbv3 = cell["cell_basis_vector3"]
            origin = cell["cell_origin"]
            lines.append("cellBasisVector1    {:.5f}     {:.5f}     {:.5f}".format(*cbv1))
            lines.append("cellBasisVector2    {:.5f}     {:.5f}     {:.5f}".format(*cbv2))
            lines.append("cellBasisVector3    {:.5f}     {:.5f}     {:.5f}".format(*cbv3))
            lines.append("cellOrigin          {:.5f}     {:.5f}     {:.5f}".format(*origin))
            lines.append("")

        if wrap_all:
            lines.append("wrapAll                         on")
            lines.append("")

        lines.append("PME                             yes")
        lines.append("PMEGridSpacing                  1.0")
        if cell:
            px, py, pz = cell["pme_grid_size"]
            lines.append("PMEGridSizeX                    {}".format(px))
            lines.append("PMEGridSizeY                    {}".format(py))
            lines.append("PMEGridSizeZ                    {}".format(pz))
        lines.append("PMETolerance                    1e-06")
        lines.append("PMEInterpOrder                  4")
        lines.append("")
    elif implicit_solvent:
        # . NAMD's Generalized Born Implicit Solvent model -- the
        #   window's own alternative to PME+explicit water whenever
        #   there's no periodic cell to do electrostatics against
        #   (vacuum/implicit-solvent runs). `cell` is ignored entirely
        #   in this branch, even if one was passed in.
        lines.append("GBIS                            on")
        lines.append("solventDielectric               {}".format(solvent_dielectric))
        lines.append("ionConcentration                {}".format(ion_concentration))
        lines.append("alphaCutoff                     {}".format(alpha_cutoff))
        lines.append("SASA                            {}".format("on" if sasa else "off"))
        if sasa:
            lines.append("surfaceTension                  {}".format(surface_tension))
        lines.append("")
    # . else: `pbc` is off AND `implicit_solvent` is off -- a bare
    #   vacuum run, nothing to emit here at all (no PME, no cell, no
    #   GBIS -- just the plain nonbonded cutoff already set above).

    if restraints:
        lines.append("constraints on")
        lines.append('consref   "{}"'.format(restraints["pdb"]))
        lines.append('conskfile "{}"'.format(restraints["pdb"]))
        lines.append("conskcol  {}".format(restraints["column"]))
        lines.append("constraintScaling {}".format(restraints.get("scaling", 1.0)))
        lines.append("")
    else:
        lines.append("constraints off")
        lines.append("")

    if fixed_atoms_file:
        lines.append("fixedAtoms                      on")
        lines.append('fixedAtomsFile                  "{}"'.format(fixed_atoms_file))
        lines.append("fixedAtomsCol                   O")
        lines.append("")

    lines.append("SMD                             {}".format("on" if smd else "off"))
    if smd:
        dx, dy, dz = smd_dir
        lines.append('SMDFile                         "{}"'.format(smd_file))
        lines.append("SMDk                            {}".format(smd_k))
        if smd_k2:
            lines.append("SMDk2                           {}".format(smd_k2))
        lines.append("SMDVel                          {}".format(smd_vel))
        lines.append("SMDDir                          {} {} {}".format(dx, dy, dz))
        lines.append("SMDOutputFreq                   {}".format(smd_output_freq))
    lines.append("")

    lines.append("outputName                      {}".format(output_basename))
    lines.append("restartfreq                     {}".format(restartfreq))
    lines.append("dcdfreq                         {}".format(dcdfreq))
    lines.append("outputEnergies                  {}".format(outputenergies))
    lines.append("outputPressure                  {}".format(outputpressure))
    lines.append("")

    if minimize_steps:
        lines.append("minimize                        {}".format(minimize_steps))

    if heating:
        lines.append("set temperature_initial   {}".format(heating["initial"]))
        lines.append("set temperature_final     {}".format(heating["final"]))
        lines.append("set heating_steps         {}".format(heating["n_increments"]))
        lines.append("set increment [expr {($temperature_final - $temperature_initial) / $heating_steps}]")
        lines.append("for {set i 0} {$i < $heating_steps} {incr i} {")
        lines.append("    set current_temp [expr {$temperature_initial + $i * $increment}]")
        lines.append("    langevinTemp $current_temp")
        lines.append("    run {}".format(heating["steps_per_increment"]))
        lines.append("}")
    elif run_steps:
        lines.append("run                             {}".format(run_steps))

    return "\n".join(lines) + "\n"


_ERROR_RE = re.compile(r"^ERROR:", re.MULTILINE)
_WRITING_COORDINATES_RE = re.compile(r"WRITING COORDINATES")


def run_namd(config_path, workdir, namd_command, log_path, n_procs=None):
    """ Runs `namd_command [+p n_procs] config_path`, redirecting
        stdout+stderr to `log_path` (matches the reference protocol's
        own run.sh, which does the same with plain shell redirection:
        `namd3 +p16 01.namd.tcl > log_01.txt`). Non-blocking --
        returns the subprocess.Popen object immediately; the caller
        (prepare_namd_run.py) polls process.poll() and tails log_path
        via GLib.timeout_add, same pattern already used for job
        tracking in process_manager_window.py.

        config_path is passed as a path RELATIVE to workdir when
        possible (matching run.sh's convention and keeping the log
        readable), falling back to the absolute path otherwise.
    """
    os.makedirs(workdir, exist_ok=True)

    command = [namd_command]
    if n_procs:
        command += ["+p", str(n_procs)]

    try:
        config_arg = os.path.relpath(config_path, workdir)
    except ValueError:
        config_arg = config_path
    command.append(config_arg)

    log_file = open(log_path, "w")
    process = subprocess.Popen(
        command, cwd=workdir, stdout=log_file, stderr=subprocess.STDOUT,
    )
    return process


def check_namd_log_success(log_path):
    """ Best-effort check of a finished run's log: True if NAMD's own
        "WRITING COORDINATES" (the very last thing it prints when
        writing the final -- not just periodic restart -- output
        files) is present AND no line starts with "ERROR:". A run can
        legitimately still be producing output (mid-simulation) when
        this returns False -- callers should only call this AFTER
        process.poll() confirms the process has actually exited.
    """
    if not os.path.isfile(log_path):
        return False
    with open(log_path, "r", errors="replace") as f:
        text = f.read()
    return bool(_WRITING_COORDINATES_RE.search(text)) and not _ERROR_RE.search(text)
