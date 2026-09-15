#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Ramachandran (backbone phi/psi) analysis
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
#      Pure NumPy (no SciPy, matching this project's own convention --
#      see md_analysis.py) computation of the protein backbone dihedral
#      angles phi/psi used by Ramachandran plots. GTK-free, so it is
#      independently testable and reusable outside
#      ramachandran_analysis_window.py.
#
import numpy as np

# A real peptide bond (C(i-1)-N(i)) is ~1.33 A. Anything much longer than
# that means the two residues are only adjacent by RESIDUE NUMBERING, not
# by an actual covalent bond -- a chain break (missing/unresolved
# residues), two separate chains sharing one chain letter, or a
# non-sequential numbering scheme. Without this check, phi/psi computed
# across such a gap would be a numerically valid but chemically
# meaningless angle (the same class of problem as the Reaction Path
# tool's own "no convergence check" caveat -- see
# docs/normal_modes_symmetry_tutorial.md's sibling caveats elsewhere in
# this project for the same philosophy: skip and say so, don't silently
# report nonsense).
_PEPTIDE_BOND_MAX_DISTANCE = 2.0

_GLYCINE_NAMES = {"GLY"}
_PROLINE_NAMES = {"PRO"}


def _signed_dihedral(r1, r2, r3, r4):
    """ Signed dihedral angle (degrees, -180 to 180) for four points, in
        the same sign convention as util/geometric_analysis.get_dihedral
        (which matches pDynamo3's own convention).

        Returns None instead of raising on degenerate geometry (collinear
        or coincident atoms) -- this is called in bulk over a whole
        structure/trajectory, and one degenerate case (e.g. two atoms
        that truly overlap in a broken structure) should not abort the
        rest of the analysis.
    """
    r1 = np.asarray(r1, dtype=np.float64)
    r2 = np.asarray(r2, dtype=np.float64)
    r3 = np.asarray(r3, dtype=np.float64)
    r4 = np.asarray(r4, dtype=np.float64)

    b1 = r2 - r1
    b2 = r3 - r2
    b3 = r4 - r3

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    b2_norm = np.linalg.norm(b2)

    if b2_norm == 0.0 or np.linalg.norm(n1) == 0.0 or np.linalg.norm(n2) == 0.0:
        return None

    b2_unit = b2 / b2_norm
    m1 = np.cross(n1, b2_unit)

    x = np.dot(n1, n2)
    y = np.dot(m1, n2)

    return -float(np.degrees(np.arctan2(y, x)))


def _find_backbone_atom(residue, name):
    for atom in residue.atoms.values():
        if atom.name.strip() == name:
            return atom
    return None


def compute_phi_psi(vismol_object, frame_indices=None):
    """ Computes the backbone phi/psi dihedral angles for every protein
        residue of `vismol_object`, chain by chain.

        frame_indices: list of frame indices (int) to compute over -- one
        point per (residue, frame) pair is returned, i.e. for a
        trajectory this is a SCATTER over conformations, not an average.
        Phi/psi are circular quantities; arithmetically averaging them
        across frames would not be physically meaningful (e.g. averaging
        +179 and -179 should give ~180, not ~0). If None (default), a
        single point per residue is computed from each atom's own
        "current frame" (Atom.coords()'s default -- see
        graphics_engine/.../model/atom.py).

        Residues are identified as protein via `residue.is_protein`
        (vismol's own classifier) AND by requiring the three backbone
        atoms N, CA, C to actually be present by name -- a residue
        missing any of them (e.g. an unresolved backbone atom in an
        X-ray structure) is silently skipped rather than crashing.

        phi needs the preceding residue's C atom; psi needs the
        following residue's N atom. Before using either, the distance to
        that neighbouring atom is checked against
        _PEPTIDE_BOND_MAX_DISTANCE -- residues that are only "adjacent"
        by residue numbering (a chain break, a gap of unresolved
        residues, ...) do not get a phi/psi computed across that gap.
        The first and last residue of each chain therefore always lack a
        phi or a psi respectively, as is standard practice for
        Ramachandran plots (e.g. PROCHECK excludes both termini).

        Returns a list of dicts, each with:
            'chain'   : chain name (str, as read from the file)
            'resname' : residue name (str, e.g. "ALA")
            'resnum'  : residue number (as read from the file)
            'frame'   : the frame index this point was computed from
                        (None if frame_indices was left as None)
            'phi'     : float, degrees
            'psi'     : float, degrees
            'category': one of 'glycine', 'proline', 'pre-proline',
                        'general' -- the four groups conventionally shown
                        with separate colours/regions on a Ramachandran
                        plot, since glycine (no side chain) and proline
                        (its ring locks phi) have markedly different
                        allowed regions than a generic residue, and the
                        residue immediately before a proline ("pre-Pro")
                        is also measurably more restricted.
        Only residues where BOTH phi and psi could be computed are
        included -- a Ramachandran point needs both.
    """
    if frame_indices is None:
        frame_indices = [None]

    results = []

    for chain_name, chain in vismol_object.chains.items():
        ordered_residues = [r for _, r in sorted(chain.residues.items(), key=lambda kv: kv[0])]

        backbone = []
        for residue in ordered_residues:
            if not residue.is_protein:
                backbone.append(None)
                continue
            n_atom  = _find_backbone_atom(residue, "N")
            ca_atom = _find_backbone_atom(residue, "CA")
            c_atom  = _find_backbone_atom(residue, "C")
            if n_atom is None or ca_atom is None or c_atom is None:
                backbone.append(None)
            else:
                backbone.append((n_atom, ca_atom, c_atom))

        n_residues = len(ordered_residues)
        for i, residue in enumerate(ordered_residues):
            this = backbone[i]
            if this is None:
                continue
            n_i, ca_i, c_i = this
            prev_bb = backbone[i - 1] if i > 0 else None
            next_bb = backbone[i + 1] if i < n_residues - 1 else None

            resname = residue.name.strip().upper()
            if resname in _GLYCINE_NAMES:
                category = "glycine"
            elif resname in _PROLINE_NAMES:
                category = "proline"
            else:
                next_residue = ordered_residues[i + 1] if i < n_residues - 1 else None
                if next_residue is not None and next_residue.name.strip().upper() in _PROLINE_NAMES:
                    category = "pre-proline"
                else:
                    category = "general"

            for frame in frame_indices:
                phi = None
                psi = None

                if prev_bb is not None:
                    r_c_prev = np.asarray(prev_bb[2].coords(frame=frame), dtype=np.float64)
                    r_n      = np.asarray(n_i.coords(frame=frame), dtype=np.float64)
                    if np.linalg.norm(r_n - r_c_prev) <= _PEPTIDE_BOND_MAX_DISTANCE:
                        phi = _signed_dihedral(r_c_prev, r_n,
                                                ca_i.coords(frame=frame),
                                                c_i.coords(frame=frame))

                if next_bb is not None:
                    r_c      = np.asarray(c_i.coords(frame=frame), dtype=np.float64)
                    r_n_next = np.asarray(next_bb[0].coords(frame=frame), dtype=np.float64)
                    if np.linalg.norm(r_n_next - r_c) <= _PEPTIDE_BOND_MAX_DISTANCE:
                        psi = _signed_dihedral(n_i.coords(frame=frame),
                                                ca_i.coords(frame=frame),
                                                r_c, r_n_next)

                if phi is None or psi is None:
                    continue

                results.append({
                    "chain": chain_name,
                    "resname": residue.name,
                    "resnum": residue.index,
                    "frame": frame,
                    "phi": phi,
                    "psi": psi,
                    "category": category,
                })

    return results
