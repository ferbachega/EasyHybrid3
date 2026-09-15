#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Generic plain-PDB file fixups (GTK/pDynamo-free)
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
#      Small, dependency-free PDB text fixups shared across import paths
#      -- kept separate from vina_runner.py because this one is a
#      general PDB-import concern, not specific to docking.
#
import shutil


def dedupe_pdb_atom_names(pdb_path, output_path=None):
    """ Rewrites ATOM/HETATM atom names (PDB columns 13-16) that
        DUPLICATE another atom's name within the same residue (chain +
        resSeq + resName), so every atom in a residue ends up with a
        unique name.

        Why this matters: pDynamo3's own PDB reader treats two atoms
        sharing a name within one residue as a "Duplicate ATOM/HETATM
        record" and SILENTLY DROPS every atom past the first one --
        confirmed against a real ligand PDB (an OpenBabel-generated
        small molecule with atoms plainly named "C", "C", "C", "O",
        "Cl" -- no OpenBabel PDB output numbers same-element atoms by
        default when the input had no real per-atom names to begin
        with, e.g. a molecule built from a bare .xyz): importing it into
        EasyHybrid via ImportSystem() silently kept only 4 of the
        molecule's real atom count, with no error show to the user --
        this is what actually produced "a PDBQT with many errors" when
        later preparing that already-atom-deficient loaded Object as a
        ligand for AutoDock Vina.

        ONLY actual duplicates are touched: a residue with no repeated
        atom name is left COMPLETELY untouched (byte-identical), which
        is what makes this safe to run unconditionally on every PDB
        import in the app, not just docking-related ones -- a normal
        protein residue's atom names (N, CA, C, O, CB, ...) never
        collide with each other, so no protein import is ever affected;
        this only ever fires for the genuinely ambiguous case.

        Every atom of a residue that has ANY name collision gets
        renamed to element + a running per-residue count (e.g. all
        three carbons of a colliding "UNL" residue become C1/C2/C3,
        even the one that happened to be unique already) -- simpler and
        just as safe as trying to preserve only the non-colliding
        names, since this whole branch only ever runs for a residue
        that already had a duplicate in it.

        output_path: written there if given; otherwise pdb_path itself
        is rewritten in place. Returns the path actually written to (or
        pdb_path unchanged, if there was nothing to fix and output_path
        was not given).
    """
    with open(pdb_path, "r") as handle:
        lines = handle.readlines()

    seen_per_residue = {}
    residues_with_duplicates = set()
    for line in lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        key = (line[21], line[22:26], line[17:20])
        name = line[12:16]
        seen = seen_per_residue.setdefault(key, set())
        if name in seen:
            residues_with_duplicates.add(key)
        seen.add(name)

    if not residues_with_duplicates:
        if output_path and output_path != pdb_path:
            shutil.copy(pdb_path, output_path)
            return output_path
        return pdb_path

    per_residue_element_counts = {}
    new_lines = []
    for line in lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            new_lines.append(line)
            continue
        key = (line[21], line[22:26], line[17:20])
        if key not in residues_with_duplicates:
            new_lines.append(line)
            continue

        element = line[76:78].strip() or "".join(c for c in line[12:16] if c.isalpha())
        element = element.capitalize() if len(element) > 1 else element.upper()

        counts = per_residue_element_counts.setdefault(key, {})
        counts[element] = counts.get(element, 0) + 1
        new_name = "{}{}".format(element, counts[element])

        # . PDB's own right/left-justification convention: a 1-character
        #   element's name field starts at column 14 (index 13); a
        #   2-character element's starts at column 13 (index 12).
        if len(element) == 1:
            name_field = (" " + new_name).ljust(4)
        else:
            name_field = new_name.ljust(4)

        new_lines.append(line[:12] + name_field + line[16:])

    out = output_path or pdb_path
    with open(out, "w") as handle:
        handle.writelines(new_lines)
    return out
