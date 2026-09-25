#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: AutoDock Vina docking (receptor/ligand PDBQT prep + cross-docking runner)
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#  Maintainer:
#      Fernando Bachega <ferbachega@gmail.com> or <easyhybrid3@gmail.com>
#
#  Description:
#      GTK-free (see md_analysis.py/namd_runner.py for the same
#      convention) port of ~/programs/adtVinaFR -- a small standalone
#      project of the same author's that cross-docks every receptor
#      PDBQT in a folder against every ligand PDBQT in another folder
#      via AutoDock Vina, having converted both from plain PDB with
#      OpenBabel first. "Cross-docking a folder against a folder"
#      supports RECEPTOR ENSEMBLES (e.g. several snapshots of a
#      flexible receptor -- treated as independent RIGID conformations,
#      not as AutoDock's own per-residue flexible-sidechain docking) and
#      LIGAND ENSEMBLES (e.g. several conformers of the same molecule)
#      in one run, which is adtVinaFR's whole reason to exist over
#      calling `vina` directly once.
#
#      Differences from the original adtVinaFR scripts (behaviour
#      preserved, code adapted to this project's conventions): rewritten
#      with spaces (the originals used tabs), results returned as a list
#      of dicts instead of positional lists, PDBQT conversion and the
#      docking loop both support a `cancel_event` (threading.Event) so
#      the GUI's Stop button can abort a long-running batch between (or
#      during, for the currently-running vina job) individual jobs
#      instead of only after the whole matrix finishes, and a
#      `progress_callback` for streaming status into a GUI log view --
#      neither concept existed in the original, which was a
#      run-to-completion command-line batch script.
#
import glob
import os
import re
import shlex
import shutil
import subprocess
import time

import numpy as np


# ---------------------------------------------------------------------------
#  Executable discovery
# ---------------------------------------------------------------------------

def find_vina_executable():
    """ Locates the `vina` executable via PATH (shutil.which). Returns
        the path, or None if not found.
    """
    return shutil.which("vina")


def find_obabel_executable():
    """ Locates the `obabel` (OpenBabel) executable via PATH. Returns
        the path, or None if not found.
    """
    return shutil.which("obabel")


# ---------------------------------------------------------------------------
#  PDB -> PDBQT preparation (OpenBabel)
# ---------------------------------------------------------------------------

def list_files_with_extension(folder, extension):
    """ Sorted list of full paths to every file in `folder` whose
        extension matches (case-insensitive, with or without the
        leading dot). Sorted so ensemble ordering (frame_0000,
        frame_0001, ...) is reproducible across OSes/filesystems that
        don't guarantee directory-listing order.
    """
    return list_files_with_extensions(folder, (extension,))


# . Ligand input formats accepted by prepare_ligand_ensemble() -- every
#   one of these is natively read by OpenBabel (obabel -L formats),
#   so no new parsing code is needed here, just widening the file
#   search and letting obabel infer the format from the extension
#   itself, exactly as it already does for .pdb.
LIGAND_INPUT_EXTENSIONS = ('pdb', 'mol', 'mol2', 'sdf', 'xyz')


def list_files_with_extensions(folder, extensions):
    """ Sorted list of full paths to every file in `folder` whose
        extension matches (case-insensitive, with or without the
        leading dot) ANY of `extensions`. Sorted so ensemble ordering
        is reproducible across OSes/filesystems that don't guarantee
        directory-listing order.
    """
    wanted = {e.lower().lstrip(".") for e in extensions}
    matches = []
    for name in os.listdir(folder):
        base, ext = os.path.splitext(name)
        if ext.lower().lstrip(".") in wanted:
            matches.append(os.path.join(folder, name))
    return sorted(matches)


def convert_receptor_pdb_to_pdbqt(pdb_path, pdbqt_path, obabel_bin="obabel"):
    """ Converts a receptor PDB into a RIGID PDBQT (no torsion tree --
        `-xr`) via OpenBabel, with Gasteiger partial charges. Raises
        subprocess.CalledProcessError on failure (caller decides how to
        report it -- see prepare_receptor_ensemble()'s per-file
        try/except for the batch case).

        --AddPolarH: a real receptor PDB (an X-ray structure especially)
        routinely has NO explicit hydrogens at all -- without this flag,
        OpenBabel does not invent any, so the resulting PDBQT ends up
        with ZERO "HD" (hydrogen-bond DONOR) atoms. Confirmed on a real
        protein: every backbone amide nitrogen -- a classic H-bond
        donor via its own N-H -- got typed "NA" (acceptor) instead of
        "N"+"HD" (donor) simply because OpenBabel had no hydrogen to
        tell it the amide's H was there; adding this flag correctly
        reclassified the large majority of those nitrogens from NA to
        N+HD. Without it, Vina's scoring function cannot see the
        receptor as able to DONATE any hydrogen bond at all, silently
        degrading every pose that would normally rely on one -- this
        was reported as "OpenBabel isn't generating the PDBQTs
        correctly", and this missing flag was the reason. Only POLAR
        hydrogens are added (not full protonation like the ligand's own
        `-p` below) -- matching AutoDockTools' own prepare_receptor4.py
        convention of a united-atom receptor (nonpolar C-H stays
        implicit, merged into its carbon by the PDBQT writer itself
        regardless of whether it was ever made explicit).
    """
    cmd = [obabel_bin, pdb_path, "-O", pdbqt_path, "-xr", "--AddPolarH", "--partialcharge", "gasteiger"]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def convert_ligand_file_to_pdbqt_ensemble(input_path, pdbqt_folder, ph=7.4, obabel_bin="obabel", out_name=None):
    """ Converts ONE ligand input file -- any OpenBabel-readable format
        (.pdb/.mol/.mol2/.sdf/.xyz, see LIGAND_INPUT_EXTENSIONS) -- into
        one OR MORE flexible PDBQTs (OpenBabel builds its own torsion
        tree) in `pdbqt_folder`, protonated at `ph` with Gasteiger
        partial charges.

        Uses obabel's own `-m` ("multiple output files") flag
        UNCONDITIONALLY, regardless of format: a compound-library file
        (.sdf/.mol2 with several molecules in one file) is split into
        one PDBQT PER MOLECULE this way -- confirmed this is exactly
        as safe for an ordinary single-molecule file (.pdb, or a
        single-entry .mol/.xyz/.sdf): `obabel single.pdb -O out.pdbqt
        -m` still produces exactly one file, just numbered
        "out1.pdbqt" rather than "out.pdbqt" (never left unsuffixed),
        so this one code path handles both cases uniformly with no
        format-specific branching needed.

        `out_name`, if given, overrides the output prefix (default: the
        input file's own basename) -- prepare_ligand_ensemble() passes
        one in when two input files of DIFFERENT formats share the same
        basename (e.g. "aspirin.pdb" and "aspirin.sdf" both present in
        the same folder), which would otherwise collide on the same
        output prefix and silently overwrite each other's PDBQT.

        Returns the sorted list of PDBQT paths obabel actually
        produced (more than one for a multi-molecule input) -- each
        becomes its own separate "ligand" in the Results table later,
        named "<out_name><N>" (e.g. a "compound_library.sdf" with 5
        molecules yields compound_library1..5), same as any other
        independently-prepared ligand.
    """
    name = out_name or os.path.splitext(os.path.basename(input_path))[0]
    out_prefix = os.path.join(pdbqt_folder, name)
    cmd = [obabel_bin, input_path, "-O", out_prefix + ".pdbqt", "-m",
           "-p", str(ph), "--partialcharge", "gasteiger"]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    produced = glob.glob(glob.escape(out_prefix) + "[0-9]*.pdbqt")
    produced.sort(key=lambda p: int(re.search(r'(\d+)\.pdbqt$', p).group(1)))
    return produced


def prepare_receptor_ensemble(pdb_folder, pdbqt_folder, obabel_bin="obabel"):
    """ Converts every .pdb in `pdb_folder` into a rigid .pdbqt in
        `pdbqt_folder` (created if missing). Returns
        (pdbqt_paths, failures) -- failures is a list of (pdb_path,
        error_message) for files OpenBabel could not convert; the batch
        continues past a single bad file rather than aborting (same
        "collect every failure, keep going" choice RMSF/RDF-adjacent
        batch operations in this project already make).
    """
    os.makedirs(pdbqt_folder, exist_ok=True)
    pdbqt_paths = []
    failures = []
    for pdb_path in list_files_with_extension(pdb_folder, "pdb"):
        name = os.path.splitext(os.path.basename(pdb_path))[0]
        out_path = os.path.join(pdbqt_folder, name + ".pdbqt")
        try:
            convert_receptor_pdb_to_pdbqt(pdb_path, out_path, obabel_bin=obabel_bin)
            pdbqt_paths.append(out_path)
        except subprocess.CalledProcessError as error:
            failures.append((pdb_path, error.stderr or str(error)))
    return pdbqt_paths, failures


def prepare_ligand_ensemble(ligand_folder, pdbqt_folder, ph=7.4, obabel_bin="obabel"):
    """ Converts every ligand input file in `ligand_folder` -- any mix
        of .pdb/.mol/.mol2/.sdf/.xyz, see LIGAND_INPUT_EXTENSIONS -- into
        PDBQT(s) in `pdbqt_folder` (created if missing), flexible and
        protonated at `ph`. A multi-molecule .sdf/.mol2 (e.g. a small
        compound library exported as one file) is split into ONE PDBQT
        PER MOLECULE (see convert_ligand_file_to_pdbqt_ensemble()) --
        each ends up as its own separate "ligand" in the Results table.

        Returns (pdbqt_paths, failures) -- failures is a list of
        (input_path, error_message) for files OpenBabel could not
        convert AT ALL; the batch continues past a single bad file
        rather than aborting (same "collect every failure, keep going"
        choice RMSF/RDF-adjacent batch operations in this project
        already make).
    """
    os.makedirs(pdbqt_folder, exist_ok=True)
    input_paths = list_files_with_extensions(ligand_folder, LIGAND_INPUT_EXTENSIONS)

    # Two input files with the SAME basename but different formats
    # (e.g. "aspirin.pdb" and "aspirin.sdf" both present) would
    # otherwise both convert to an "aspirin1.pdbqt" output prefix and
    # silently overwrite each other -- fold the source extension into
    # the output name, but ONLY for basenames that actually collide, so
    # the common (non-colliding) case keeps clean "aspirin1.pdbqt"-style
    # names in the Results table.
    basenames = {}
    for path in input_paths:
        base = os.path.splitext(os.path.basename(path))[0]
        basenames.setdefault(base, []).append(path)
    colliding_basenames = {base for base, paths in basenames.items() if len(paths) > 1}

    pdbqt_paths = []
    failures = []
    for input_path in input_paths:
        base = os.path.splitext(os.path.basename(input_path))[0]
        if base in colliding_basenames:
            ext = os.path.splitext(input_path)[1].lstrip('.')
            out_name = '{}_{}'.format(base, ext)
        else:
            out_name = base
        try:
            pdbqt_paths.extend(convert_ligand_file_to_pdbqt_ensemble(
                input_path, pdbqt_folder, ph=ph, obabel_bin=obabel_bin, out_name=out_name))
        except subprocess.CalledProcessError as error:
            failures.append((input_path, error.stderr or str(error)))
    return pdbqt_paths, failures


# ---------------------------------------------------------------------------
#  Docking box
# ---------------------------------------------------------------------------

def compute_centroid(coords):
    """ The plain geometric centroid (arithmetic mean position) of
        `coords` (array-like, shape (N, 3)) -- deliberately different
        from compute_box_from_coordinates()'s centre (the midpoint of
        the bounding box), which the two coincide on only for a
        symmetric arrangement of atoms. For a lopsided selection (e.g.
        a pocket with residues bulging more to one side) the centroid
        sits where the atoms' MASS actually concentrates, while the
        bounding-box midpoint sits wherever the two most extreme atoms
        happen to be, which can be pulled well away from the rest of
        the selection by a single outlying atom.

        Returns a dict with only 'center_x'/'center_y'/'center_z' --
        unlike compute_box_from_coordinates(), a centroid has no
        associated size, so the caller is expected to keep whatever Size
        is already set.

        Raises ValueError if `coords` is empty.
    """
    coords = np.asarray(coords, dtype=np.float64)
    if coords.shape[0] == 0:
        raise ValueError("compute_centroid() needs at least one atom.")
    center = coords.mean(axis=0)
    return {"center_x": float(center[0]), "center_y": float(center[1]), "center_z": float(center[2])}


def compute_box_from_coordinates(coords, padding=5.0):
    """ A docking search box tightly enclosing `coords` (array-like,
        shape (N, 3)), expanded by `padding` Angstrom on EVERY side (so
        each box dimension grows by 2*padding total) -- the same
        "bounding box of a selection, plus a margin for the ligand to
        move around in" convention AutoDockTools' own grid-box-from-
        selection helper uses. Returns a box dict with the same 6 keys
        run_vina_docking()'s `box` argument expects.

        Raises ValueError if `coords` is empty (no box can be computed
        from zero atoms -- this is the caller's job to catch, e.g. an
        empty 3D-view selection).
    """
    coords = np.asarray(coords, dtype=np.float64)
    if coords.shape[0] == 0:
        raise ValueError("compute_box_from_coordinates() needs at least one atom.")
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)
    center = (mins + maxs) / 2.0
    size = (maxs - mins) + 2.0 * padding
    return {
        "center_x": float(center[0]), "center_y": float(center[1]), "center_z": float(center[2]),
        "size_x": float(size[0]), "size_y": float(size[1]), "size_z": float(size[2]),
    }


# ---------------------------------------------------------------------------
#  Running Vina
# ---------------------------------------------------------------------------

def read_pdb_coordinates(pdb_path):
    """ Reads the x/y/z of every ATOM/HETATM record of a plain PDB
        (fixed columns 31-38/39-46/47-54, the standard PDB convention --
        the same one extract_pose_as_pdb()'s own output follows), in
        file order. Returns an (N, 3) float64 array.

        Used to grow an ALREADY-LOADED vobject's own .frames with one
        more pose without re-importing a whole new pDynamo System for
        it -- valid ONLY when every pose really does share the same
        atom count/ORDER as frame 0 (see read_pdb_elements() below,
        used to check that BEFORE trusting this function's raw
        coordinates -- two poses of the same ligand are not guaranteed
        to have matching atom order just because they came from Vina:
        see that function's own docstring for the confirmed case where
        they don't).
    """
    coords = []
    with open(pdb_path, "r") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return np.asarray(coords, dtype=np.float64)


def read_pdb_elements(pdb_path):
    """ Reads the element symbol of every ATOM/HETATM record of a plain
        PDB, in file order (same row-for-row correspondence as
        read_pdb_coordinates()'s own output) -- from columns 77-78 when
        present, else derived from the atom name. Returns a list of
        1-2 character strings, e.g. ['C', 'C', 'O', 'C', 'Cl'].

        Used to VERIFY, before trusting read_pdb_coordinates() to grow
        an existing vobject's .frames, that a later pose's atoms are
        really in the same order as frame 0's -- confirmed NOT to be
        guaranteed in general: every pose of one Vina job (one
        receptor/ligand PDBQT pair, every MODE in its output) does
        share the same order (Vina never reorders atoms itself), but
        DIFFERENT ligand PDBQT files (e.g. several conformers of one
        ligand exported from a loaded Object's multiple frames, each
        converted by OpenBabel SEPARATELY) are only guaranteed to share
        the same connectivity, not the same OpenBabel-assigned internal
        atom order -- OpenBabel's own torsion-tree root/branch selection
        can depend on which atom happens to anchor the largest rigid
        fragment, which is not always conformation-independent for a
        ligand with more than one rotatable bond. Silently trusting
        matching atom COUNT alone (as an earlier version of this tool
        did) let exactly this slip through: a correct first pose, and
        every later one built from a differently-ordered file quietly
        assigned to the wrong atoms -- "chemically correct first frame,
        completely distorted geometry every frame after it".
    """
    elements = []
    with open(pdb_path, "r") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            element = line[76:78].strip() or "".join(c for c in line[12:16] if c.isalpha())
            elements.append(element.capitalize() if len(element) > 1 else element.upper())
    return elements


def extract_pose_as_pdb(obabel_bin, docked_pdbqt, mode, out_pdb):
    """ Extracts ONE binding mode (Vina numbers them from 1) out of a
        multi-MODEL docked PDBQT into its own plain PDB, via OpenBabel's
        -f/-l (first/last model) flags, then rewrites its atom names to
        be unique (see _make_pdb_atom_names_unique() below) before
        returning.

        Why the rewrite is needed: OpenBabel's PDB writer, converting a
        ligand that had no real per-atom names to begin with (typical
        for a docking ligand built from a bare .xyz/.pdbqt with no
        naming convention), names every atom after its own element ALONE
        ("C", "C", "C", "O", "CL", ...) with no disambiguating number.
        That is invalid PDB (atom names must be unique within a residue)
        and was confirmed, against a real docked pose, to make
        pDynamo3's own strict PDB reader (used by EasyHybrid's
        vm_session.load()) log "Duplicate ATOM/HETATM record" warnings
        and SILENTLY DROP every atom past the first one sharing a name
        -- a 5-heavy-atom ligand (3 carbons all plainly named "C")
        reloaded into EasyHybrid with only 3 atoms, no error shown to
        the user. Renaming every atom to element+running-index (C1, C2,
        C3, O1, CL1, ...) before handing the file to vm_session.load()
        fixes this at the source instead of leaving it as a silent trap
        for whoever next needed a docked pose reloaded into the viewer.
    """
    cmd = [obabel_bin, docked_pdbqt, "-f", str(mode), "-l", str(mode), "-O", out_pdb]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    _make_pdb_atom_names_unique(out_pdb)
    return out_pdb


def _make_pdb_atom_names_unique(pdb_path):
    """ Rewrites every ATOM/HETATM record's atom-name field (PDB
        columns 13-16) in place so no two atoms of the same residue
        share a name -- see extract_pose_as_pdb()'s docstring for why
        this matters. Renames strictly by first-appearance order within
        each (chain, resSeq, resName) residue, as element symbol +
        1-based running count for that element in that residue (e.g.
        the 2nd carbon of residue UNL 1 becomes "C2"). The element is
        read from columns 77-78 when present (what OpenBabel writes),
        falling back to the existing name stripped of digits/spaces.
    """
    with open(pdb_path, "r") as handle:
        lines = handle.readlines()

    counters = {}  # (chain, resseq, resname) -> {element: count}
    new_lines = []
    for line in lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            new_lines.append(line)
            continue

        chain = line[21]
        resseq = line[22:26]
        resname = line[17:20]
        element = line[76:78].strip() or "".join(c for c in line[12:16] if c.isalpha())
        element = element.capitalize() if len(element) > 1 else element.upper()

        key = (chain, resseq, resname)
        element_counts = counters.setdefault(key, {})
        element_counts[element] = element_counts.get(element, 0) + 1
        new_name = "{}{}".format(element, element_counts[element])

        # . PDB's own right/left-justification convention: a 1-character
        #   element's name field starts at column 14 (index 13); a
        #   2-character element's starts at column 13 (index 12).
        if len(element) == 1:
            name_field = (" " + new_name).ljust(4)
        else:
            name_field = new_name.ljust(4)

        new_lines.append(line[:12] + name_field + line[16:])

    with open(pdb_path, "w") as handle:
        handle.writelines(new_lines)


_BOX_TO_VINA_FLAG = {
    "center_x": "--center_x", "center_y": "--center_y", "center_z": "--center_z",
    "size_x": "--size_x", "size_y": "--size_y", "size_z": "--size_z",
}
_PARAM_TO_VINA_FLAG = {
    "cpu": "--cpu", "seed": "--seed", "exhaustiveness": "--exhaustiveness",
    "num_modes": "--num_modes", "energy_range": "--energy_range",
}


def build_vina_command(vina_bin, receptor_pdbqt, ligand_pdbqt, out_pdbqt, box, params=None):
    """ Builds the full `vina` command line for one receptor/ligand
        pair. `box` must have all 6 center_*/size_* keys (see
        compute_box_from_coordinates()); `params` (optional) may have
        any of cpu/seed/exhaustiveness/num_modes/energy_range -- a
        missing or None-valued key is simply omitted, letting Vina fall
        back to its own default (same convention the original
        adtVinaFR's `parameters` dict used).
    """
    cmd = [vina_bin]
    for key, flag in _BOX_TO_VINA_FLAG.items():
        cmd += [flag, str(box[key])]
    for key, flag in _PARAM_TO_VINA_FLAG.items():
        value = (params or {}).get(key)
        if value is not None:
            cmd += [flag, str(value)]
    cmd += ["--receptor", receptor_pdbqt, "--ligand", ligand_pdbqt, "--out", out_pdbqt]
    return cmd


def write_vina_config(conf_path, receptor_pdbqt, ligand_pdbqt, out_pdbqt, box, params=None):
    """ Writes a native AutoDock Vina config file (plain "key = value"
        lines, one per parameter) that reproduces one docking job
        exactly -- confirmed directly that `vina --config <this file>`
        runs identically to the equivalent CLI flags (same key names as
        _BOX_TO_VINA_FLAG/_PARAM_TO_VINA_FLAG, just without the leading
        "--" and with "=" instead of a space). Lets a user inspect,
        tweak, or rerun a single job entirely outside EasyHybrid, with
        every parameter sitting in one plain-text file instead of a
        long command line.
    """
    lines = ['receptor = {}'.format(receptor_pdbqt), 'ligand = {}'.format(ligand_pdbqt)]
    for key in ('center_x', 'center_y', 'center_z', 'size_x', 'size_y', 'size_z'):
        lines.append('{} = {}'.format(key, box[key]))
    for key in ('exhaustiveness', 'num_modes', 'energy_range', 'cpu', 'seed'):
        value = (params or {}).get(key)
        if value is not None:
            lines.append('{} = {}'.format(key, value))
    lines.append('out = {}'.format(out_pdbqt))
    with open(conf_path, 'w') as handle:
        handle.write('\n'.join(lines) + '\n')


def write_job_shell_script(sh_path, vina_bin, conf_path):
    """ Writes a standalone, executable shell script that reproduces
        ONE docking job outside EasyHybrid by running
        `<vina_bin> --config <conf_path>` -- the .conf file (see
        write_vina_config()) already has every parameter -- receptor,
        ligand, box, search settings, output path -- baked in, so this
        script needs no arguments and can be copied/moved/rerun on its
        own (as long as the .conf's own absolute paths stay valid).
    """
    with open(sh_path, 'w') as handle:
        handle.write('#!/bin/sh\n')
        handle.write('# Reproduces this exact AutoDock Vina docking job outside EasyHybrid.\n')
        handle.write('{} --config {}\n'.format(shlex.quote(vina_bin), shlex.quote(conf_path)))
    os.chmod(sh_path, 0o755)


def write_batch_shell_script(sh_path, job_sh_paths):
    """ Writes a standalone, executable shell script that reproduces
        an ENTIRE docking batch (every receptor x ligand pair actually
        attempted) by running each job's own script (see
        write_job_shell_script()) in the same order the batch itself
        ran them. Does NOT stop at the first failed job (no `set -e`)
        -- matches this tool's own "collect every failure, keep going"
        batch philosophy (see run_vina_docking()'s per-job try/except).
    """
    with open(sh_path, 'w') as handle:
        handle.write('#!/bin/sh\n')
        handle.write('# Reproduces this entire AutoDock Vina docking batch outside EasyHybrid --\n')
        handle.write('# runs every job below in the same order the original batch did.\n')
        for job_sh_path in job_sh_paths:
            handle.write('sh {}\n'.format(shlex.quote(job_sh_path)))
    os.chmod(sh_path, 0o755)


def parse_vina_output_log(log_path):
    """ Parses Vina's own results table out of its captured stdout (see
        run_single_docking_job() -- Vina >= 1.2 dropped the old --log
        flag and prints its results table to stdout instead, so the
        caller redirects stdout into this file itself). Returns a list
        of dicts, one per binding mode: {'mode': int, 'affinity':
        float, 'rmsd_lb': float, 'rmsd_ub': float}. Lines that are not
        exactly "<int> <float> <float> <float>" (headers, banners,
        warnings, ...) are silently skipped, matching the original
        adtVinaFR parser's own tolerant behaviour.
    """
    results = []
    with open(log_path, "r", errors="replace") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                mode = int(parts[0])
                affinity = float(parts[1])
                rmsd_lb = float(parts[2])
                rmsd_ub = float(parts[3])
            except ValueError:
                continue
            results.append({"mode": mode, "affinity": affinity, "rmsd_lb": rmsd_lb, "rmsd_ub": rmsd_ub})
    return results


def count_pdbqt_models(pdbqt_path):
    """ Counts "MODEL" records in a PDBQT file. A file with no MODEL
        lines at all (a single-pose PDBQT) counts as 1 -- same
        convention pdbqt_parser.parse_pdbqt_models() uses for a
        headerless single-pose file.
    """
    count = 0
    with open(pdbqt_path, "r", errors="replace") as handle:
        for line in handle:
            if line.startswith("MODEL"):
                count += 1
    return count if count > 0 else 1


def run_single_docking_job(vina_bin, receptor_pdbqt, ligand_pdbqt, out_pdbqt, log_path,
                            box, params=None, cancel_event=None, poll_interval=0.2):
    """ Runs ONE Vina job as a polled subprocess (rather than a blocking
        subprocess.run()) so `cancel_event` (a threading.Event) can
        terminate it mid-run, not just between jobs -- a single
        exhaustive Vina search can take minutes, and a GUI's Stop
        button should not have to wait for the current job to finish on
        its own.

        Returns (returncode, results) -- results is
        parse_vina_output_log(log_path)'s own return, FILTERED down to
        only the modes Vina actually wrote to `out_pdbqt` (see below),
        or [] if the process was cancelled or failed before producing a
        table. returncode is None if cancelled before the process
        exited (the caller should not treat that as success OR as a
        normal Vina failure).

        Vina's own printed results table can list MORE modes than it
        actually writes to --out: confirmed directly (AutoDock Vina
        1.2.3, real receptor/ligand pair, default parameters) -- a
        severely clashing pose with an absurd +142 kcal/mol affinity
        was printed as "mode 9" in the table, but the output PDBQT only
        ever contained MODEL 1 through 8 (a contiguous prefix, the
        worst mode simply missing) -- Vina evidently drops that one
        pose during its own write-out step for reasons of its own. This
        is what a user reported as "the results table shows 10
        different affinities but the PDBQT only has one mode/pose" --
        a MUCH more severe case of the exact same mismatch (most of the
        table's rows had no corresponding MODEL at all). Without this
        filter, the Results table shows rows nothing can actually be
        imported from, and importing them raised UNCAUGHT StopIteration
        crashes (confirmed) in the docking windows' own pose-import
        code, which assumed every logged mode has a real MODEL behind
        it. Filtering here, once, at the source, means every OTHER
        consumer of these results (the Results table, "Import all
        poses", DockingResult metadata) never has to know this Vina
        quirk exists.
    """
    cmd = build_vina_command(vina_bin, receptor_pdbqt, ligand_pdbqt, out_pdbqt, box, params)
    with open(log_path, "w") as logfile:
        process = subprocess.Popen(cmd, stdout=logfile, stderr=subprocess.STDOUT)
        try:
            while True:
                returncode = process.poll()
                if returncode is not None:
                    break
                if cancel_event is not None and cancel_event.is_set():
                    process.terminate()
                    try:
                        process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        process.kill()
                    return None, []
                time.sleep(poll_interval)
        except Exception:
            process.kill()
            raise

    if returncode != 0:
        return returncode, []
    results = parse_vina_output_log(log_path)
    if not os.path.isfile(out_pdbqt):
        return returncode, []
    n_models = count_pdbqt_models(out_pdbqt)
    results = [r for r in results if r["mode"] <= n_models]
    return returncode, results


def run_vina_docking(vina_bin, receptor_pdbqt_list, ligand_pdbqt_list, output_folder,
                      box, params=None, progress_callback=None, cancel_event=None):
    """ Cross-docks every receptor in `receptor_pdbqt_list` against
        every ligand in `ligand_pdbqt_list` (the full N x M matrix,
        same as the original adtVinaFR) into `output_folder`
        (<ligand>_<receptor>_docked.pdbqt + _vina.log per pair, plus one
        aggregate affinity_logs.txt tab-separated table, both matching
        the original's own file-naming/columns for anyone used to
        reading them).

        progress_callback(event, info), if given, is called with:
            ('job_start',    {'receptor':.., 'ligand':.., 'index':.., 'total':..})
            ('job_done',     {'receptor':.., 'ligand':.., 'rows': [...]})
            ('job_failed',   {'receptor':.., 'ligand':.., 'error': str})
            ('cancelled',    {})
        -- meant to be used from a background thread, e.g.
        `lambda event, info: GLib.idle_add(handler, event, info)`; this
        function itself does no GUI work and has no GTK/vismol import.

        Returns the full list of result rows across every pair, each a
        dict: {'receptor':.., 'ligand':.., 'mode':.., 'affinity':..,
        'rmsd_lb':.., 'rmsd_ub':.., 'docked_pdbqt':.., 'log':..}.
    """
    os.makedirs(output_folder, exist_ok=True)

    affinity_log_path = os.path.join(output_folder, "affinity_logs.txt")
    with open(affinity_log_path, "w") as fulllog:
        fulllog.write("ligand\treceptor\tmode\taffinity\trmsd_l.b.\trmsd_u.b.\n")

    # . job_scripts accumulates every per-job reproducer .sh actually
    #   written (regardless of whether that job succeeded, failed, or
    #   the batch was later cancelled) -- write_batch_shell_script() is
    #   called at every exit point below so "rerun this batch outside
    #   EasyHybrid" always reflects exactly what's really on disk, even
    #   for a cancelled/partial run.
    job_scripts = []
    batch_sh_path = os.path.join(output_folder, "run_all.sh")

    all_rows = []
    total = len(receptor_pdbqt_list) * len(ligand_pdbqt_list)
    index = 0
    for receptor_pdbqt in receptor_pdbqt_list:
        receptor_name = os.path.splitext(os.path.basename(receptor_pdbqt))[0]
        for ligand_pdbqt in ligand_pdbqt_list:
            if cancel_event is not None and cancel_event.is_set():
                if progress_callback:
                    progress_callback("cancelled", {})
                write_batch_shell_script(batch_sh_path, job_scripts)
                return all_rows

            ligand_name = os.path.splitext(os.path.basename(ligand_pdbqt))[0]
            index += 1
            if progress_callback:
                progress_callback("job_start", {
                    "receptor": receptor_name, "ligand": ligand_name,
                    "index": index, "total": total,
                })

            out_pdbqt = os.path.join(output_folder, "{}_{}_docked.pdbqt".format(ligand_name, receptor_name))
            log_path = os.path.join(output_folder, "{}_{}_vina.log".format(ligand_name, receptor_name))
            conf_path = os.path.join(output_folder, "{}_{}.conf".format(ligand_name, receptor_name))
            job_sh_path = os.path.join(output_folder, "{}_{}_run.sh".format(ligand_name, receptor_name))

            # . Written BEFORE running the job (not just on success) so
            #   a failed or cancelled job still leaves behind a way to
            #   manually retry it outside EasyHybrid.
            write_vina_config(conf_path, receptor_pdbqt, ligand_pdbqt, out_pdbqt, box, params=params)
            write_job_shell_script(job_sh_path, vina_bin, conf_path)
            job_scripts.append(job_sh_path)

            try:
                returncode, results = run_single_docking_job(
                    vina_bin, receptor_pdbqt, ligand_pdbqt, out_pdbqt, log_path,
                    box, params=params, cancel_event=cancel_event)
            except OSError as error:
                if progress_callback:
                    progress_callback("job_failed", {
                        "receptor": receptor_name, "ligand": ligand_name, "error": str(error)})
                continue

            if returncode is None:
                if progress_callback:
                    progress_callback("cancelled", {})
                write_batch_shell_script(batch_sh_path, job_scripts)
                return all_rows

            if returncode != 0:
                if progress_callback:
                    progress_callback("job_failed", {
                        "receptor": receptor_name, "ligand": ligand_name,
                        "error": "vina exited with code {} -- see {}".format(returncode, log_path)})
                continue

            rows = []
            with open(affinity_log_path, "a") as fulllog:
                for result in results:
                    row = {
                        "receptor": receptor_name, "ligand": ligand_name,
                        "mode": result["mode"], "affinity": result["affinity"],
                        "rmsd_lb": result["rmsd_lb"], "rmsd_ub": result["rmsd_ub"],
                        "docked_pdbqt": out_pdbqt, "log": log_path,
                    }
                    rows.append(row)
                    all_rows.append(row)
                    fulllog.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        ligand_name, receptor_name, result["mode"],
                        result["affinity"], result["rmsd_lb"], result["rmsd_ub"]))

            if progress_callback:
                progress_callback("job_done", {"receptor": receptor_name, "ligand": ligand_name, "rows": rows})

    write_batch_shell_script(batch_sh_path, job_scripts)
    return all_rows


def parse_affinity_log(output_folder):
    """ Reads a PREVIOUS run's own affinity_logs.txt (written by
        run_vina_docking() above) and reconstructs full result rows --
        the exact same shape run_vina_docking() itself returns --
        pointing at the `<ligand>_<receptor>_docked.pdbqt`/`_vina.log`
        files that should already exist alongside it (same naming
        convention run_vina_docking() itself uses). Lets "Import
        previous results" add an earlier run's results to the Results
        table without re-running anything.

        Returns (rows, missing_pairs) -- `rows` only includes a row
        when its own `docked_pdbqt` file still actually exists at that
        path (e.g. not moved/deleted since the run); `missing_pairs` is
        a sorted list of (ligand, receptor) pairs for which it did NOT
        (so the caller can report exactly what could not be
        recovered, instead of silently dropping rows).
    """
    log_path = os.path.join(output_folder, "affinity_logs.txt")
    rows = []
    missing_pairs = set()
    with open(log_path, "r", errors="replace") as handle:
        handle.readline()  # header
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 6:
                continue
            ligand, receptor, mode, affinity, rmsd_lb, rmsd_ub = parts
            docked_pdbqt = os.path.join(output_folder, "{}_{}_docked.pdbqt".format(ligand, receptor))
            log_file = os.path.join(output_folder, "{}_{}_vina.log".format(ligand, receptor))
            if not os.path.isfile(docked_pdbqt):
                missing_pairs.add((ligand, receptor))
                continue
            rows.append({
                "receptor": receptor, "ligand": ligand, "mode": int(mode),
                "affinity": float(affinity), "rmsd_lb": float(rmsd_lb), "rmsd_ub": float(rmsd_ub),
                "docked_pdbqt": docked_pdbqt, "log": log_file,
            })
    return rows, sorted(missing_pairs)
