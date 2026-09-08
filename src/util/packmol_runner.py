#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Packmol input generation and execution
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
#      Pure-Python (no GTK) helpers for prepare_packmol_run.py: builds a
#      Packmol .inp file from a list of components and runs the `packmol`
#      executable as a subprocess, mirroring namd_runner.py's/
#      tleap_runner.py's own conventions (non-blocking Popen, GLib polling
#      done by the caller, log-based success check).
#
#      Constraint keyword syntax (box/sphere/fixed argument order and
#      count) confirmed against the real Packmol Fortran source shipped
#      with this machine's AmberTools install (getinp.f90's restriction
#      parser, ~/programs/ambertools25/lib/python3.12/site-packages/
#      packmol_memgen/lib/packmol/getinp.f90) rather than assumed from
#      memory of the docs alone.
#
import os
import shutil
import subprocess


# . Each entry: (constraint_type, packmol keyword line template, number of
#   float params it takes, whether "number" is meaningful for this type).
#   "fixed" places exactly one copy at a given position/orientation --
#   Packmol itself errors out if `number` > 1 is given for a fixed
#   structure (see packmol.f90's own check), so the caller always emits
#   "number 1" for it rather than omitting the line.
CONSTRAINT_TYPES = {
    "inside_box":     {"label": "Inside box",    "n_params": 6},
    "outside_box":    {"label": "Outside box",   "n_params": 6},
    "inside_sphere":  {"label": "Inside sphere", "n_params": 4},
    "outside_sphere": {"label": "Outside sphere","n_params": 4},
    "fixed":          {"label": "Fixed position/orientation", "n_params": 6},
}


def find_packmol_executable():
    """ Locates the `packmol` executable via PATH (shutil.which) -- same
        auto-detect-then-let-the-user-override pattern as
        namd_runner.find_namd_executable()/tleap_runner.find_tleap_executable().
        Returns the path, or None if not found.
    """
    return shutil.which("packmol")


def _constraint_line(constraint_type, params):
    """ Builds the single Packmol restriction line (everything after the
        `structure`/`number` lines, before `end structure`) for one
        component, given its constraint_type (a CONSTRAINT_TYPES key) and
        its raw float params in the exact order Packmol expects:

          inside_box / outside_box  -> x1 y1 z1 x2 y2 z2  (two opposite
                                        corners of an axis-aligned box)
          inside_sphere/outside_sphere -> cx cy cz r
          fixed                     -> x y z alpha beta gamma (degrees)
    """
    if constraint_type == "inside_box":
        return "inside box {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}".format(*params)
    if constraint_type == "outside_box":
        return "outside box {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}".format(*params)
    if constraint_type == "inside_sphere":
        return "inside sphere {:.4f} {:.4f} {:.4f} {:.4f}".format(*params)
    if constraint_type == "outside_sphere":
        return "outside sphere {:.4f} {:.4f} {:.4f} {:.4f}".format(*params)
    if constraint_type == "fixed":
        return "fixed {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}".format(*params)
    raise ValueError("Unknown Packmol constraint type: {!r}".format(constraint_type))


def build_packmol_input(input_path, output_path, components, tolerance=2.0,
                         seed=None, filetype="pdb"):
    """ Writes a Packmol .inp file at `input_path`.

        components: list of dicts, each shaped like
            {"structure_path": str, "number": int,
             "constraint_type": one of CONSTRAINT_TYPES, "params": tuple}
        one `structure ... end structure` block per entry, in list order
        (Packmol packs earlier structures first, which matters when a
        "fixed" solute needs to already occupy space before solvent is
        packed around it -- so the caller/UI should let the user order
        this list, fixed components first).

        tolerance is Packmol's own `tolerance` directive (Angstrom,
        default 2.0 -- the minimum distance accepted between atoms of
        different molecules).

        seed=None omits the `seed` directive entirely (Packmol then uses
        its own default/random seed); pass an int to make a packing run
        reproducible.
    """
    lines = []
    lines.append("tolerance {:.4f}".format(tolerance))
    lines.append("filetype {}".format(filetype))
    lines.append('output "{}"'.format(output_path))
    if seed is not None:
        lines.append("seed {}".format(int(seed)))
    lines.append("")

    for component in components:
        structure_path = component["structure_path"]
        constraint_type = component["constraint_type"]
        params = component["params"]

        lines.append('structure "{}"'.format(structure_path))
        if constraint_type == "fixed":
            lines.append("  number 1")
        else:
            lines.append("  number {}".format(int(component["number"])))
        lines.append("  " + _constraint_line(constraint_type, params))
        lines.append("end structure")
        lines.append("")

    with open(input_path, "w") as f:
        f.write("\n".join(lines) + "\n")


def run_packmol(input_path, workdir, packmol_command, log_path):
    """ Runs `packmol_command < input_path`, redirecting stdout+stderr to
        `log_path` -- Packmol reads its entire configuration from STDIN
        (unlike NAMD/tleap, it takes no config-file command-line argument
        at all), so the input file is wired up as the subprocess's stdin
        rather than appended to `command`.

        Non-blocking -- returns the subprocess.Popen object immediately;
        the caller polls process.poll() and tails log_path via
        GLib.timeout_add, same pattern as namd_runner.run_namd().
    """
    os.makedirs(workdir, exist_ok=True)

    stdin_file = open(input_path, "r")
    log_file = open(log_path, "w")
    process = subprocess.Popen(
        [packmol_command], cwd=workdir, stdin=stdin_file,
        stdout=log_file, stderr=subprocess.STDOUT,
    )
    return process


def check_packmol_log_success(log_path):
    """ Best-effort check of a finished run's log: True if Packmol's own
        "Solution written to file" line (its definitive "packing done and
        a structure was actually written" message -- confirmed against
        packmol.f90's own write statements) is present AND no line starts
        with "ERROR". A run can legitimately still be mid-optimization
        (writing intermediate "Current point written to file" snapshots)
        when this returns False -- callers should only call this AFTER
        process.poll() confirms the process has actually exited.
    """
    if not log_path or not os.path.isfile(log_path):
        return False
    with open(log_path, "r", errors="replace") as f:
        text = f.read()
    if any(line.strip().startswith("ERROR") for line in text.splitlines()):
        return False
    return "Solution written to file" in text
