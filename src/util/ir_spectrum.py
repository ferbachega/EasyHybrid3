#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Infrared spectrum construction from computed normal modes
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
#      Turns a list of (frequency, IR intensity) stick data -- from
#      pSimulation.NormalModes.NormalModes_InfraredIntensities(), saved
#      per-mode by p_methods/normal_modes.py -- into a continuous,
#      broadened spectrum for plotting, the same way a real IR
#      instrument would report one. Pure NumPy, no external dependency,
#      matching every other analysis tool added this session (RDF, RMSF).
#
import numpy as np

# [EN] BUG FIX (found testing a pure-MM/DYFF water system with real
# charges assigned): the near-zero-frequency "modes" left over from an
# imperfectly projected translation/rotation basis are NOT physically
# meaningless the way a simple `frequency <= 0` filter assumes -- their
# finite-difference dipole DERIVATIVE can come out numerically huge
# (verified: a leftover rotational mode at +4.4 cm^-1 reported an
# "intensity" of ~380, larger than every real vibrational band in the
# same molecule) because rotating a dipole changes its Cartesian
# components a lot even though no real absorption is happening. A
# frequency-sign filter alone lets these through as spurious giant
# peaks near zero. Filtering by |frequency| below a threshold instead
# discards both the negative (imaginary-mode) and small-positive (RT
# residue) cases uniformly.
_MINIMUM_REAL_FREQUENCY = 50.0  # cm^-1


def broaden_spectrum(frequencies, intensities, x_min=0.0, x_max=4000.0,
                      num_points=2000, fwhm=20.0, lineshape='lorentzian',
                      minimum_frequency=_MINIMUM_REAL_FREQUENCY):
    """ Convolves stick (frequency, intensity) data into a continuous curve.

        frequencies/intensities: equal-length sequences (cm^-1 / km mol^-1).
        Non-real vibrations are skipped: imaginary modes (frequency < 0)
        and translation/rotation residues (|frequency| below
        minimum_frequency) -- see _MINIMUM_REAL_FREQUENCY's docstring
        above for why the residues need a magnitude threshold, not just
        a sign check.

        fwhm: full width at half maximum of each individual band, in
        cm^-1. Not a computed quantity -- a display choice, same as
        every real IR-plotting program lets you pick (instrumental/
        homogeneous broadening is not something a single-point harmonic
        calculation predicts).

        Returns (x, y): x the frequency grid (cm^-1, num_points values
        from x_min to x_max), y the summed intensity at each grid point.
    """
    frequencies, intensities = stick_data(frequencies, intensities, minimum_frequency=minimum_frequency)

    x = np.linspace(x_min, x_max, num_points)
    y = np.zeros_like(x)

    if frequencies.size == 0:
        return x, y

    halfWidth = fwhm / 2.0
    for freq, inten in zip(frequencies, intensities):
        if lineshape == 'gaussian':
            sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
            y += inten * np.exp(-0.5 * ((x - freq) / sigma) ** 2)
        else:
            y += inten * (halfWidth ** 2) / ((x - freq) ** 2 + halfWidth ** 2)

    return x, y


def stick_data(frequencies, intensities, minimum_frequency=_MINIMUM_REAL_FREQUENCY):
    """ Returns (frequencies, intensities) filtered to only the real
        vibrations -- frequency >= minimum_frequency -- for plotting as
        a stick spectrum (e.g. one vertical line per mode) alongside the
        broadened curve. Excludes imaginary modes AND small-magnitude
        translation/rotation residues (see module docstring): both can
        report finite-difference "intensities" with no physical meaning.
    """
    frequencies = np.asarray(frequencies, dtype=float)
    intensities = np.asarray(intensities, dtype=float)
    keep = frequencies >= minimum_frequency
    return frequencies[keep], intensities[keep]
