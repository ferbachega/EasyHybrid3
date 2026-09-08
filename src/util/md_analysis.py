#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Trajectory analysis (RDF, RMSF)
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
#      Pure NumPy (no SciPy, no MDAnalysis/MDTraj -- this project deliberately
#      avoids SciPy, see the "no scipy required" commit) trajectory-analysis
#      routines: the pair radial distribution function g(r) and the
#      per-atom root-mean-square fluctuation (RMSF). GTK-free, so these are
#      independently testable and reusable outside the GUI windows that
#      call them (rdf_analysis_window.py, rmsf_analysis_window.py).
#
import numpy as np


def compute_rmsf(coords_per_frame):
    """ Per-atom RMSF (root-mean-square fluctuation around the mean
        position over the given frames), in the same length unit as the
        input coordinates.

        coords_per_frame: array-like, shape (n_frames, n_atoms, 3) -- the
        ALREADY-SELECTED subset of atoms/frames to analyse (frame range,
        step size, and atom selection are all the caller's responsibility;
        this function has no notion of "the whole trajectory").

        No superposition/alignment is performed here -- RMSF is only
        physically meaningful once rigid-body translation/rotation has
        already been removed from the trajectory (e.g. via EasyHybrid's
        own "Align Trajectory" tool). Computing it on an unaligned
        trajectory silently mixes in whole-molecule motion; this function
        does not attempt to detect that case.

        Returns a 1D array, shape (n_atoms,).
    """
    coords = np.asarray(coords_per_frame, dtype=np.float64)
    mean = coords.mean(axis=0)
    displacement = coords - mean[np.newaxis, :, :]
    mean_square = np.sum(displacement * displacement, axis=-1).mean(axis=0)
    return np.sqrt(mean_square)


def _minimum_image(delta, box):
    """ Orthorhombic minimum-image convention: box is (a, b, c). Matches
        this project's own box model elsewhere (e.g. namd_runner.py,
        vismol_object's cell_coordinates) -- rectangular cells only, no
        triclinic support.
    """
    return delta - box * np.round(delta / box)


def compute_rdf(coords_a_per_frame, coords_b_per_frame, box_per_frame,
                 same_selection, r_max=15.0, bin_width=0.1):
    """ The pair radial distribution function g(r) between group A and
        group B (which may be the same group -- see `same_selection`),
        averaged over a set of frames, plus the running coordination
        number N(r) (the average number of B atoms found within r of an
        A atom).

        coords_a_per_frame, coords_b_per_frame: array-like, each shape
            (n_frames, n_X_atoms, 3) -- the already-selected coordinates
            for every frame to include (same frame count/order for both).
        box_per_frame: array-like, shape (n_frames, 3), the (a, b, c)
            orthorhombic box for each of those same frames. A periodic
            box is required to normalise g(r) by a bulk density; there is
            no meaningful non-periodic g(r) here.
        same_selection: True when A and B are literally the same group of
            atoms (self-RDF, e.g. O-O) -- excludes an atom from being its
            own neighbour. Do not set this for two selections that merely
            happen to overlap partially; only an EXACT same-group RDF is
            handled specially.
        r_max, bin_width: histogram range/resolution, same length unit as
            the input coordinates (Angstrom, throughout this codebase).

        Returns (r_centres, g_r, coordination_number), each a 1D array of
        length ceil(r_max / bin_width).

        Density convention: the bulk number density of B in a frame is
        taken as N_B / V (not (N_B-1) / V for the self-RDF case) -- the
        standard large-N approximation used by common MD analysis tools;
        the resulting O(1/N) bias is negligible for any system large
        enough for a radial distribution function to be a meaningful
        thing to compute in the first place.
    """
    n_frames = len(coords_a_per_frame)
    if n_frames == 0:
        raise ValueError("compute_rdf() needs at least one frame.")

    n_bins = int(np.ceil(r_max / bin_width))
    edges = np.arange(n_bins + 1, dtype=np.float64) * bin_width
    histogram_total = np.zeros(n_bins, dtype=np.float64)

    sum_n_a = 0
    sum_density_b = 0.0

    for frame_index in range(n_frames):
        coords_a = np.asarray(coords_a_per_frame[frame_index], dtype=np.float64)
        coords_b = np.asarray(coords_b_per_frame[frame_index], dtype=np.float64)
        box = np.asarray(box_per_frame[frame_index], dtype=np.float64)

        delta = coords_a[:, np.newaxis, :] - coords_b[np.newaxis, :, :]
        delta = _minimum_image(delta, box)
        distances = np.sqrt(np.sum(delta * delta, axis=-1))

        if same_selection:
            np.fill_diagonal(distances, np.inf)

        in_range = distances[distances <= r_max]
        histogram, _ = np.histogram(in_range, bins=edges)
        histogram_total += histogram

        volume = box[0] * box[1] * box[2]
        sum_n_a += coords_a.shape[0]
        sum_density_b += coords_b.shape[0] / volume

    mean_n_a = sum_n_a / n_frames
    mean_density_b = sum_density_b / n_frames

    r_centres = 0.5 * (edges[:-1] + edges[1:])
    shell_volumes = (4.0 / 3.0) * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)

    ideal_counts = mean_n_a * mean_density_b * shell_volumes * n_frames
    with np.errstate(divide="ignore", invalid="ignore"):
        g_r = np.where(ideal_counts > 0.0, histogram_total / ideal_counts, 0.0)

    coordination_number = np.cumsum(histogram_total) / (mean_n_a * n_frames)

    return r_centres, g_r, coordination_number
