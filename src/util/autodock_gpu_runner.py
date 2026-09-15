#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: AutoDock-GPU docking (AutoGrid4 grid-map generation + AutoDock-GPU runner)
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
#      GTK-free (see vina_runner.py for the same convention) wrapping of
#      AutoDock-GPU (~/programs/adtGPU/adgpu-v1.6_linux_x64_cuda11_128wi
#      -- a CUDA build, needs an NVIDIA GPU) and AutoGrid4. Unlike Vina,
#      AutoDock-GPU does NOT compute its own search grid per job -- it
#      reads PRECOMPUTED affinity grid maps (one per atom type, plus
#      electrostatics/desolvation) that AutoGrid4 has to generate first
#      from a GPF (Grid Parameter File) text file, once per receptor
#      (every ligand docked against that receptor reuses the same
#      maps). Receptor/ligand PDBQT preparation and the docking box
#      itself are NOT duplicated here -- both are engine-agnostic and
#      already live in vina_runner.py (prepare_receptor_ensemble(),
#      prepare_ligand_ensemble(), compute_box_from_coordinates(),
#      compute_centroid()); this module only adds what is genuinely
#      specific to AutoDock-GPU: writing a GPF, running AutoGrid4,
#      running AutoDock-GPU itself, and parsing its .dlg text output.
#
#      AutoDock-GPU's own .dlg output wraps each docked pose's PDBQT
#      lines with a "DOCKED: " prefix (MODEL/ATOM/ENDMDL, one MODEL per
#      GA run) -- dlg_to_docked_pdbqt() strips that prefix into a plain
#      multi-MODEL PDBQT, deliberately so the result is byte-for-byte
#      the same shape as one of Vina's own *_docked.pdbqt files: every
#      pose-extraction/import helper in vina_runner.py
#      (extract_pose_as_pdb(), read_pdb_coordinates(),
#      read_pdb_elements()) then works on an AutoDock-GPU result
#      completely unmodified.
#
import math
import os
import re
import shlex
import shutil
import subprocess


# ---------------------------------------------------------------------------
#  Executable discovery
# ---------------------------------------------------------------------------

def find_autogrid_executable():
    """ Locates the `autogrid4` executable: PATH first (shutil.which),
        then the user-local pip/conda-style install location this was
        confirmed installed at on this machine (~/.local/bin/autogrid4)
        -- autogrid4 has no equivalent of $AMBERHOME to fall back on.
        Returns the path, or None if not found.
    """
    found = shutil.which("autogrid4")
    if found:
        return found
    local = os.path.expanduser("~/.local/bin/autogrid4")
    if os.path.isfile(local) and os.access(local, os.X_OK):
        return local
    return None


def find_autodock_gpu_executable():
    """ AutoDock-GPU ships as a single standalone binary with a
        version/build-specific name (e.g.
        "adgpu-v1.6_linux_x64_cuda11_128wi") -- there is no standard
        PATH name to guess the way `vina`/`obabel` have. Only a plain
        PATH search is attempted; the window falls back to asking the
        user to browse for it (same as namd_runner.py's own
        find_namd_executable() has to for a similarly non-standardised
        binary name, minus the couple of conventional names NAMD does
        have).
    """
    for name in ("autodock_gpu", "autodock-gpu", "adgpu"):
        found = shutil.which(name)
        if found:
            return found
    return None


# ---------------------------------------------------------------------------
#  Grid Parameter File (GPF) + AutoGrid4
# ---------------------------------------------------------------------------

def parse_pdbqt_atom_types(pdbqt_path):
    """ The set of AutoDock atom types actually present in a PDBQT (its
        ATOM/HETATM records' own last whitespace-separated column),
        sorted for a deterministic GPF. Needed for a GPF's own
        `receptor_types`/`ligand_types` lines -- AutoGrid4 only
        computes a map for a type it is explicitly told about.
    """
    types = set()
    with open(pdbqt_path, "r") as handle:
        for line in handle:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                types.add(line.split()[-1])
    return sorted(types)


def compute_npts(size, spacing):
    """ AutoGrid4 requires an EVEN number of grid points per axis --
        the smallest even integer whose (npts * spacing) covers at
        least `size` Angstrom.
    """
    n = int(math.ceil(size / spacing))
    if n % 2 != 0:
        n += 1
    return n


def write_gpf(gpf_path, receptor_pdbqt_name, receptor_types, ligand_types,
              center, size, spacing=0.375, smooth=0.5, dielectric=-0.1465):
    """ Writes a Grid Parameter File for AutoGrid4. `receptor_pdbqt_name`
        is used BARE (no directory) -- AutoGrid4 resolves every filename
        in a GPF relative to its OWN current working directory, so the
        caller is expected to run it with cwd set to the folder holding
        both the GPF and the receptor PDBQT (see run_autogrid() below).

        center, size: (x, y, z) tuples in Angstrom -- the SAME box
        convention vina_runner.compute_box_from_coordinates()'s dict
        uses (just unpacked into a plain tuple by the caller).

        Returns the "gridfld" filename this GPF declares (e.g.
        "receptor.maps.fld") -- what run_autodock_gpu()'s own `--ffile`
        needs afterward.
    """
    base = os.path.splitext(receptor_pdbqt_name)[0]
    fld_name = base + ".maps.fld"
    npts = tuple(compute_npts(s, spacing) for s in size)

    lines = [
        "npts {} {} {}".format(*npts),
        "gridfld {}".format(fld_name),
        "spacing {}".format(spacing),
        "receptor_types {}".format(" ".join(receptor_types)),
        "ligand_types {}".format(" ".join(ligand_types)),
        "receptor {}".format(receptor_pdbqt_name),
        "gridcenter {} {} {}".format(*center),
        "smooth {}".format(smooth),
    ]
    for ligand_type in ligand_types:
        lines.append("map {}.{}.map".format(base, ligand_type))
    lines.append("elecmap {}.e.map".format(base))
    lines.append("dsolvmap {}.d.map".format(base))
    lines.append("dielectric {}".format(dielectric))

    with open(gpf_path, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    return fld_name


def run_autogrid(autogrid_bin, gpf_name, log_name, cwd):
    """ Runs AutoGrid4 with cwd set to the folder containing the GPF
        (and the receptor PDBQT it references) -- both `gpf_name` and
        `log_name` are used BARE, matching write_gpf()'s own
        filename-resolution convention.

        BARE names + a short `cwd` are not just a style choice: AutoGrid4
        (confirmed directly, this version) crashes with "*** buffer
        overflow detected ***" -- a stack-smashing abort from an
        internal fixed-size filename buffer -- when a GPF's own
        `receptor`/`gridfld`/`map` entries are long absolute paths
        (~140 characters was enough to trigger it). Always run it with
        cwd set to a dedicated, short-named folder and reference every
        file in the GPF by its bare filename, never a full path.

        Returns (returncode, stdout, stderr); check_autogrid_log_success()
        on the resulting log file is the real success test (a returncode
        of 0 is necessary but AutoGrid4 has been seen to still return 0
        after warnings worth surfacing).
    """
    cmd = build_autogrid_command(autogrid_bin, gpf_name, log_name)
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    return result.returncode, result.stdout, result.stderr


def build_autogrid_command(autogrid_bin, gpf_name, log_name):
    """ The AutoGrid4 command line alone (no `cwd` -- that's a separate
        subprocess.run()/shell-`cd` concern, not part of the argv
        itself). Factored out of run_autogrid() so both it and
        write_autogrid_shell_script() build the exact same command.
    """
    return [autogrid_bin, "-p", gpf_name, "-l", log_name]


def check_autogrid_log_success(log_path):
    """ True iff AutoGrid4's own log ends with its "Successful
        Completion" banner -- the same "trust the tool's own final
        banner, not just returncode" convention this project's
        namd_runner.check_namd_log_success() already uses.
    """
    if not os.path.isfile(log_path):
        return False
    with open(log_path, "r", errors="replace") as handle:
        return "Successful Completion" in handle.read()


# ---------------------------------------------------------------------------
#  Running AutoDock-GPU
# ---------------------------------------------------------------------------

def run_autodock_gpu(adgpu_bin, ligand_pdbqt_path, fld_path,
                      nrun=20, seed=None, devnum=None, resnam=None,
                      autostop=None, heurmax=None, nev=None, ngen=None,
                      lsit=None, psize=None):
    """ Runs AutoDock-GPU. Unlike run_autogrid()/AutoGrid4 (which
        crashes -- a fixed-size internal buffer, "buffer overflow
        detected" -- on a long/absolute path, confirmed directly),
        AutoDock-GPU itself has no such issue: `ligand_pdbqt_path`,
        `fld_path`, and `resnam` may all be full absolute paths,
        confirmed working directly, so no particular `cwd` is needed
        here at all.

        `resnam`, if given, is used as the FULL output path prefix for
        the resulting .dlg/.xml (default: the ligand's own basename, in
        the current directory) -- always pass one when docking more
        than one ligand so their results don't overwrite each other.

        `autostop`/`heurmax`/`nev`/`ngen`/`lsit`/`psize` map directly to
        AutoDock-GPU's own --autostop/--heurmax/--nev/--ngen/--lsit/
        --psize search-parameter flags (see `adgpu_bin --help`); each is
        left off the command line (letting AutoDock-GPU use its own
        built-in default) when passed as None.

        Returns (returncode, stdout, stderr). AutoDock-GPU prints "All
        jobs ran without errors." on success; combined with returncode
        == 0 and the .dlg file actually existing, that is the success
        test the caller should use (no separate log-file convention
        here the way NAMD/AutoGrid4 have -- stdout IS the log).
    """
    cmd = build_autodock_gpu_command(
        adgpu_bin, ligand_pdbqt_path, fld_path, nrun=nrun, seed=seed, devnum=devnum, resnam=resnam,
        autostop=autostop, heurmax=heurmax, nev=nev, ngen=ngen, lsit=lsit, psize=psize)
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode, result.stdout, result.stderr


def build_autodock_gpu_command(adgpu_bin, ligand_pdbqt_path, fld_path,
                                nrun=20, seed=None, devnum=None, resnam=None,
                                autostop=None, heurmax=None, nev=None, ngen=None,
                                lsit=None, psize=None):
    """ The AutoDock-GPU command line alone. Factored out of
        run_autodock_gpu() so both it and write_autodock_gpu_shell_script()
        build the exact same command -- see run_autodock_gpu()'s own
        docstring for what each parameter means.
    """
    cmd = [adgpu_bin, "--lfile", ligand_pdbqt_path, "--ffile", fld_path, "--nrun", str(nrun)]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    if devnum is not None:
        cmd += ["--devnum", str(devnum)]
    if autostop is not None:
        cmd += ["--autostop", "1" if autostop else "0"]
    if heurmax is not None:
        cmd += ["--heurmax", str(heurmax)]
    if nev is not None:
        cmd += ["--nev", str(nev)]
    if ngen is not None:
        cmd += ["--ngen", str(ngen)]
    if lsit is not None:
        cmd += ["--lsit", str(lsit)]
    if psize is not None:
        cmd += ["--psize", str(psize)]
    if resnam:
        cmd += ["--resnam", resnam]
    return cmd


def write_autogrid_shell_script(sh_path, autogrid_bin, gpf_name, log_name, cwd):
    """ Writes a standalone, executable shell script that reproduces
        ONE receptor's AutoGrid4 grid-map generation outside EasyHybrid
        -- `cd`s into `cwd` first, same requirement run_autogrid()
        itself has (AutoGrid4 crashes on long/absolute paths in its GPF
        -- see run_autogrid()'s own docstring), so this script only
        works correctly if the whole `cwd` folder (with its GPF and
        receptor PDBQT) is still there, or is moved along with it.
    """
    cmd = build_autogrid_command(autogrid_bin, gpf_name, log_name)
    with open(sh_path, 'w') as handle:
        handle.write('#!/bin/sh\n')
        handle.write('# Reproduces this AutoGrid4 grid-map generation outside EasyHybrid.\n')
        handle.write('cd {}\n'.format(shlex.quote(cwd)))
        handle.write(' '.join(shlex.quote(part) for part in cmd) + '\n')
    os.chmod(sh_path, 0o755)


def write_autodock_gpu_shell_script(sh_path, adgpu_bin, ligand_pdbqt_path, fld_path,
                                     nrun=20, seed=None, devnum=None, resnam=None,
                                     autostop=None, heurmax=None, nev=None, ngen=None,
                                     lsit=None, psize=None):
    """ Writes a standalone, executable shell script that reproduces
        ONE AutoDock-GPU docking job outside EasyHybrid -- no `cwd`
        needed (see run_autodock_gpu()'s own docstring: unlike
        AutoGrid4, AutoDock-GPU has no issue with absolute paths), so
        this script works from anywhere as long as the referenced
        ligand PDBQT and .maps.fld (and its .map files alongside it)
        still exist at their original paths.
    """
    cmd = build_autodock_gpu_command(
        adgpu_bin, ligand_pdbqt_path, fld_path, nrun=nrun, seed=seed, devnum=devnum, resnam=resnam,
        autostop=autostop, heurmax=heurmax, nev=nev, ngen=ngen, lsit=lsit, psize=psize)
    with open(sh_path, 'w') as handle:
        handle.write('#!/bin/sh\n')
        handle.write('# Reproduces this exact AutoDock-GPU docking job outside EasyHybrid.\n')
        handle.write(' '.join(shlex.quote(part) for part in cmd) + '\n')
    os.chmod(sh_path, 0o755)


def write_batch_shell_script(sh_path, job_sh_paths):
    """ Writes a standalone, executable shell script that reproduces an
        ENTIRE AutoDock-GPU docking batch (every grid-generation and
        docking job actually attempted) by running each job's own
        script (see write_autogrid_shell_script()/
        write_autodock_gpu_shell_script()) in the same order the batch
        itself ran them -- a receptor's grid script always precedes its
        own ligand jobs, matching run_autodock_gpu_docking()'s own
        per-receptor-then-per-ligand loop order. Does NOT stop at the
        first failed job (no `set -e`) -- matches this tool's own
        "collect every failure, keep going" batch philosophy.
    """
    with open(sh_path, 'w') as handle:
        handle.write('#!/bin/sh\n')
        handle.write('# Reproduces this entire AutoDock-GPU docking batch outside EasyHybrid --\n')
        handle.write('# runs every grid-generation and docking job below in the same order the\n')
        handle.write('# original batch did.\n')
        for job_sh_path in job_sh_paths:
            handle.write('sh {}\n'.format(shlex.quote(job_sh_path)))
    os.chmod(sh_path, 0o755)


# ---------------------------------------------------------------------------
#  Parsing .dlg results
# ---------------------------------------------------------------------------

# . The final "RMSD TABLE" section's own per-pose summary line, e.g.:
#   "   1      1      7       -1.80      0.00     37.29           RANKING"
#   columns: overall Rank, Sub-Rank (position within its cluster), Run
#   number, Binding Energy (kcal/mol), Cluster RMSD (to that cluster's
#   own best pose), Reference RMSD (to --xraylfile, meaningless/large
#   here since none was given) -- every such line ends with the literal
#   word "RANKING", which is what this regex anchors on.
_RANKING_LINE_RE = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+RANKING\s*$")


def parse_dlg_ranking(dlg_path):
    """ Parses AutoDock-GPU's own final clustering/ranking table out of
        a .dlg file. Returns a list of dicts: {'rank':.., 'subrank':..,
        'run':.., 'energy':.., 'cluster_rmsd':.., 'ref_rmsd':..}. Note
        this table is grouped BY CLUSTER first (each cluster's members
        together, strongest cluster first) and only sorted by energy
        WITHIN a cluster -- it is NOT a single global sort by energy
        (a rank-2 cluster's best pose can score better than a rank-1
        cluster's weakest member) -- sort the returned list by 'energy'
        yourself if a strict best-first order is wanted regardless of
        clustering (e.g. before combining poses into one trajectory,
        matching vina_runner-based tools' own convention).
    """
    rows = []
    with open(dlg_path, "r", errors="replace") as handle:
        for line in handle:
            match = _RANKING_LINE_RE.match(line)
            if not match:
                continue
            rank, subrank, run, energy, cluster_rmsd, ref_rmsd = match.groups()
            rows.append({
                "rank": int(rank), "subrank": int(subrank), "run": int(run),
                "energy": float(energy), "cluster_rmsd": float(cluster_rmsd),
                "ref_rmsd": float(ref_rmsd),
            })
    return rows


def dlg_to_docked_pdbqt(dlg_path, out_pdbqt):
    """ Extracts every "DOCKED: "-prefixed line from a .dlg file,
        strips that prefix, and writes the result to `out_pdbqt` -- a
        plain multi-MODEL PDBQT (MODEL number == GA run number) with
        the exact same shape as one of Vina's own *_docked.pdbqt
        files. Lets every pose-extraction helper in vina_runner.py
        (extract_pose_as_pdb() etc.) work on an AutoDock-GPU result
        completely unmodified -- pass a ranking row's own 'run' value
        as the `mode` argument.
    """
    with open(dlg_path, "r", errors="replace") as handle:
        lines = [line[len("DOCKED: "):] for line in handle if line.startswith("DOCKED: ")]
    with open(out_pdbqt, "w") as handle:
        handle.writelines(lines)
    return out_pdbqt


# ---------------------------------------------------------------------------
#  Orchestration: grid maps (once per receptor) + docking (per ligand)
# ---------------------------------------------------------------------------

def run_autodock_gpu_docking(adgpu_bin, autogrid_bin, receptor_pdbqt_list, ligand_pdbqt_list,
                              grids_folder, output_folder, box, spacing=0.375,
                              nrun=20, seed=None, devnum=None,
                              autostop=None, heurmax=None, nev=None, ngen=None,
                              lsit=None, psize=None,
                              progress_callback=None, cancel_event=None):
    """ Cross-docks every receptor in `receptor_pdbqt_list` against
        every ligand in `ligand_pdbqt_list`, same N x M matrix
        vina_runner.run_vina_docking() computes -- but unlike Vina,
        AutoDock-GPU needs a grid map PER RECEPTOR generated up front
        (see write_gpf()/run_autogrid()), reused for every ligand
        docked against that same receptor rather than recomputed per
        job. Grid generation happens once per receptor, inside
        `grids_folder` (each receptor gets its own short-named
        subfolder there -- required by run_autogrid()'s own
        long-path-crashes-AutoGrid4 constraint, see its docstring); the
        `ligand_types` given to AutoGrid4 for a receptor is the UNION
        of every ligand's own atom types, computed once, so the SAME
        maps can score every ligand regardless of which one happens to
        run first.

        `box`: a dict with center_x/y/z, size_x/y/z (the same shape
        vina_runner.compute_box_from_coordinates() returns).

        progress_callback(event, info), if given, is called with:
            ('grid_start',  {'receptor':.., 'index':.., 'total':..})
            ('grid_done',   {'receptor':..})
            ('grid_failed', {'receptor':.., 'error': str})
            ('job_start',   {'receptor':.., 'ligand':.., 'index':.., 'total':..})
            ('job_done',    {'receptor':.., 'ligand':.., 'rows': [...]})
            ('job_failed',  {'receptor':.., 'ligand':.., 'error': str})
            ('cancelled',   {})
        -- meant to be used from a background thread, e.g.
        `lambda event, info: GLib.idle_add(handler, event, info)`; this
        function itself does no GUI work and has no GTK/vismol import.
        A receptor whose grid generation fails has every one of its
        ligand jobs skipped (there is nothing to dock against).

        Returns the full list of result rows across every pair, each a
        dict: {'receptor':.., 'ligand':.., 'run':.., 'energy':..,
        'cluster_rank':.., 'cluster_rmsd':.., 'docked_pdbqt':..,
        'dlg':..}.
    """
    os.makedirs(grids_folder, exist_ok=True)
    os.makedirs(output_folder, exist_ok=True)

    # . One aggregate, tab-separated log across the whole batch -- same
    #   convention vina_runner.run_vina_docking()'s own affinity_logs.txt
    #   uses (there, ligand/receptor/mode/affinity/rmsd_lb/rmsd_ub; here,
    #   this engine's own row shape instead) -- lets a later session
    #   reopen these results via parse_affinity_log() without re-running
    #   anything, exactly like the Vina side.
    affinity_log_path = os.path.join(output_folder, "affinity_logs.txt")
    with open(affinity_log_path, "w") as fulllog:
        fulllog.write("ligand\treceptor\trun\tenergy\tcluster_rank\tcluster_rmsd\n")

    center = (box["center_x"], box["center_y"], box["center_z"])
    size = (box["size_x"], box["size_y"], box["size_z"])

    ligand_types = sorted(set().union(*(set(parse_pdbqt_atom_types(l)) for l in ligand_pdbqt_list)))

    # . job_scripts accumulates every per-job reproducer .sh actually
    #   written (grid-generation AND docking, regardless of whether
    #   that step succeeded) -- write_batch_shell_script() is called at
    #   every exit point below so "rerun this batch outside EasyHybrid"
    #   always reflects exactly what's really on disk, even for a
    #   cancelled/partial run (same convention as vina_runner.py's own
    #   run_vina_docking()).
    job_scripts = []
    batch_sh_path = os.path.join(output_folder, "run_all.sh")

    all_rows = []
    total = len(receptor_pdbqt_list) * len(ligand_pdbqt_list)
    index = 0

    for receptor_index, receptor_pdbqt in enumerate(receptor_pdbqt_list):
        if cancel_event is not None and cancel_event.is_set():
            if progress_callback:
                progress_callback("cancelled", {})
            write_batch_shell_script(batch_sh_path, job_scripts)
            return all_rows

        receptor_name = os.path.splitext(os.path.basename(receptor_pdbqt))[0]
        # . A dedicated, SHORT subfolder per receptor -- not named after
        #   the receptor's own (possibly long) filename, precisely to
        #   avoid run_autogrid()'s long-path crash regardless of what a
        #   real receptor happens to be called.
        receptor_grid_dir = os.path.join(grids_folder, "receptor_{:03d}".format(receptor_index))
        os.makedirs(receptor_grid_dir, exist_ok=True)
        local_receptor = os.path.join(receptor_grid_dir, "receptor.pdbqt")
        shutil.copy(receptor_pdbqt, local_receptor)

        if progress_callback:
            progress_callback("grid_start", {
                "receptor": receptor_name, "index": receptor_index + 1, "total": len(receptor_pdbqt_list)})

        receptor_types = parse_pdbqt_atom_types(local_receptor)
        gpf_path = os.path.join(receptor_grid_dir, "receptor.gpf")
        fld_name = write_gpf(gpf_path, "receptor.pdbqt", receptor_types, ligand_types,
                              center=center, size=size, spacing=spacing)
        log_path = os.path.join(receptor_grid_dir, "autogrid.log")

        # . Written BEFORE running AutoGrid4 (not just on success) so a
        #   failed grid generation still leaves behind a way to
        #   manually retry it outside EasyHybrid. Lives in
        #   output_folder (alongside the docking jobs' own scripts),
        #   not receptor_grid_dir, so every reproducer script for one
        #   run ends up in one place.
        grid_sh_path = os.path.join(output_folder, "{}_grid_run.sh".format(receptor_name))
        write_autogrid_shell_script(grid_sh_path, autogrid_bin, "receptor.gpf", "autogrid.log", receptor_grid_dir)
        job_scripts.append(grid_sh_path)

        _rc, _out, err = run_autogrid(autogrid_bin, "receptor.gpf", "autogrid.log", cwd=receptor_grid_dir)

        if not check_autogrid_log_success(log_path):
            if progress_callback:
                progress_callback("grid_failed", {
                    "receptor": receptor_name,
                    "error": err or "see {}".format(log_path)})
            index += len(ligand_pdbqt_list)  # every ligand for this receptor is being skipped
            continue

        if progress_callback:
            progress_callback("grid_done", {"receptor": receptor_name})

        fld_path = os.path.join(receptor_grid_dir, fld_name)

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
                    "index": index, "total": total})

            resnam = os.path.join(output_folder, "{}_{}".format(ligand_name, receptor_name))

            # . Written BEFORE running AutoDock-GPU (not just on
            #   success) so a failed job still leaves behind a way to
            #   manually retry it outside EasyHybrid.
            job_sh_path = os.path.join(output_folder, "{}_{}_run.sh".format(ligand_name, receptor_name))
            write_autodock_gpu_shell_script(
                job_sh_path, adgpu_bin, ligand_pdbqt, fld_path, nrun=nrun, seed=seed, devnum=devnum,
                resnam=resnam, autostop=autostop, heurmax=heurmax, nev=nev, ngen=ngen, lsit=lsit, psize=psize)
            job_scripts.append(job_sh_path)

            rc, out, err = run_autodock_gpu(
                adgpu_bin, ligand_pdbqt, fld_path, nrun=nrun, seed=seed, devnum=devnum, resnam=resnam,
                autostop=autostop, heurmax=heurmax, nev=nev, ngen=ngen, lsit=lsit, psize=psize)
            dlg_path = resnam + ".dlg"

            if rc != 0 or "All jobs ran without errors" not in out or not os.path.isfile(dlg_path):
                if progress_callback:
                    progress_callback("job_failed", {
                        "receptor": receptor_name, "ligand": ligand_name,
                        "error": err or out[-800:] or "AutoDock-GPU exited with code {}".format(rc)})
                continue

            ranking = parse_dlg_ranking(dlg_path)
            docked_pdbqt = resnam + "_docked.pdbqt"
            dlg_to_docked_pdbqt(dlg_path, docked_pdbqt)

            rows = []
            with open(affinity_log_path, "a") as fulllog:
                for ranked in ranking:
                    row = {
                        "receptor": receptor_name, "ligand": ligand_name,
                        "run": ranked["run"], "energy": ranked["energy"],
                        "cluster_rank": ranked["rank"], "cluster_rmsd": ranked["cluster_rmsd"],
                        "docked_pdbqt": docked_pdbqt, "dlg": dlg_path,
                    }
                    rows.append(row)
                    all_rows.append(row)
                    fulllog.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        ligand_name, receptor_name, ranked["run"],
                        ranked["energy"], ranked["rank"], ranked["cluster_rmsd"]))

            if progress_callback:
                progress_callback("job_done", {"receptor": receptor_name, "ligand": ligand_name, "rows": rows})

    write_batch_shell_script(batch_sh_path, job_scripts)
    return all_rows


def parse_affinity_log(output_folder):
    """ Reads a PREVIOUS run's own affinity_logs.txt (written by
        run_autodock_gpu_docking() above) and reconstructs full result
        rows -- the exact same shape run_autodock_gpu_docking() itself
        returns -- pointing at the `<ligand>_<receptor>_docked.pdbqt`/
        `<ligand>_<receptor>.dlg` files that should already exist
        alongside it (same naming convention run_autodock_gpu_docking()
        itself uses). Lets "Import previous results" add an earlier
        run's results to the Results table without re-running anything.

        Returns (rows, missing_pairs) -- same contract as
        vina_runner.parse_affinity_log(): `rows` only includes a row
        when its own `docked_pdbqt` file still actually exists;
        `missing_pairs` is a sorted list of (ligand, receptor) pairs
        for which it did not.
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
            ligand, receptor, run, energy, cluster_rank, cluster_rmsd = parts
            docked_pdbqt = os.path.join(output_folder, "{}_{}_docked.pdbqt".format(ligand, receptor))
            dlg_path = os.path.join(output_folder, "{}_{}.dlg".format(ligand, receptor))
            if not os.path.isfile(docked_pdbqt):
                missing_pairs.add((ligand, receptor))
                continue
            rows.append({
                "receptor": receptor, "ligand": ligand, "run": int(run),
                "energy": float(energy), "cluster_rank": int(cluster_rank), "cluster_rmsd": float(cluster_rmsd),
                "docked_pdbqt": docked_pdbqt, "dlg": dlg_path,
            })
    return rows, sorted(missing_pairs)
