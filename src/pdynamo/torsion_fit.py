#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: torsion fitting of a dihedral parameter (DYFF, OPLS) against an xTB scan
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      1. Scan: one dihedral i-j-k-l is driven over 360 degrees with xtb,
#         relaxed ($constrain + $scan, one xtb run, geometries in
#         xtbscan.log) or rigid (the k side of the xtb-optimized structure
#         rotated about j-k). Every scan
#         geometry is then re-evaluated by an unconstrained xtb single point,
#         so that the reference energies carry no restraint bias.
#      2. MM energies at the same geometries with the fitted parameter set to
#         zero (every other term, charges and Lennard-Jones unchanged), given
#         by the force field's side (gui/windows/builder/dyff_parameters.py,
#         pdynamo/opls/torsion_fit.py).
#      3. Fit: the series E = sum_n c_n cos ( n phi ) (DYFF's own form; OPLS
#         k ( 1 + cos ( n phi - delta ) ) with delta 0/180 is the same up to a
#         constant: k = |c_n|, delta = 0 if c_n > 0 else 180), linear in the c_n. With all the terms that share the parameter key
#         ( d = 1..D, around the scanned bond and elsewhere ),
#
#             E_QM ( s ) - E_MM,0 ( s ) = a + sum_n c_n sum_d cos ( n phi_d ( s ) )
#
#         is solved by weighted linear least squares ( a = free offset ).
#         Weights damp the high-energy points (clashes): w = 1 up to
#         WEIGHT_ONSET above the QM minimum, exp ( -( dE - onset ) / WEIGHT_KT )
#         above it. Periods whose column (almost) vanishes ( e.g. n = 1, 2
#         for the three equivalent terms of a methyl group ) are left out.
#      4. Check: the MM energies are recomputed by pDynamo with the fitted
#         parameter (not by the formula) and compared with the QM curve
#         (RMSE after the best offset).
#
#      Units: kJ/mol and degrees throughout.
#

import math
import os
import shutil
import subprocess
import tempfile
from collections import defaultdict

import numpy as np

HARTREE_TO_KJMOL = 2625.499639
WEIGHT_ONSET     = 25.0     # kJ/mol above the QM minimum with full weight
WEIGHT_KT        = 10.0     # kJ/mol, decay of the weight above the onset
ZERO_TERM        = [ ( 0.0, 0 ) ]
RELATIVE_SPREAD  = 0.10     # periods whose column varies less than this x the largest are left out


class ScanCancelled ( Exception ):
    pass


#=============================================================================
#  Geometry
#=============================================================================
def dihedral_angle ( coordinates, quad ):
    """ Dihedral i-j-k-l in degrees, ( -180, 180 ]. """
    a, b, c, d = ( np.asarray ( coordinates[i], dtype = float ) for i in quad )
    b0, b1, b2 = a - b, c - b, d - c
    b1n = b1 / np.linalg.norm ( b1 )
    v = b0 - np.dot ( b0, b1n ) * b1n
    w = b2 - np.dot ( b2, b1n ) * b1n
    return math.degrees ( math.atan2 ( float ( np.dot ( np.cross ( b1n, v ), w ) ), float ( np.dot ( v, w ) ) ) )

def wrap180 ( x ):
    return ( x + 180.0 ) % 360.0 - 180.0

def adjacency_from_bonds ( n_atoms, bonds ):
    adjacency = defaultdict ( set )
    for ( a, b ) in bonds:
        adjacency[a].add ( b ); adjacency[b].add ( a )
    return adjacency

def moving_side ( adjacency, j, k ):
    """ Atoms on the k side of the j-k bond. ValueError when j-k is in a ring. """
    seen, stack = { k }, [ k ]
    while stack:
        u = stack.pop ( )
        for v in adjacency[u]:
            if u == k and v == j: continue
            if v == j:
                raise ValueError ( "The central bond is in a ring: it cannot be rotated over 360 degrees." )
            if v not in seen:
                seen.add ( v ); stack.append ( v )
    return seen

def rotate_about_bond ( coordinates, j, k, side, angle ):
    """ Copy of the coordinates with the atoms of side rotated by angle
        (degrees) about the j -> k axis. """
    xyz = np.array ( coordinates, dtype = float )
    origin = xyz[j]
    axis = xyz[k] - origin
    axis /= np.linalg.norm ( axis )
    t = math.radians ( angle )
    c, s = math.cos ( t ), math.sin ( t )
    K = np.array ( [ [ 0.0, -axis[2], axis[1] ], [ axis[2], 0.0, -axis[0] ], [ -axis[1], axis[0], 0.0 ] ] )
    R = np.eye ( 3 ) + s * K + ( 1.0 - c ) * ( K @ K )
    for a in side:
        xyz[a] = origin + R @ ( xyz[a] - origin )
    return xyz


#=============================================================================
#  xTB
#=============================================================================
def _write_xyz ( path, symbols, coordinates, title = "scan" ):
    with open ( path, "w" ) as handle:
        handle.write ( "{}\n{}\n".format ( len ( symbols ), title ) )
        for s, r in zip ( symbols, coordinates ):
            handle.write ( "{} {:.8f} {:.8f} {:.8f}\n".format ( s, r[0], r[1], r[2] ) )

def read_xyz_frames ( path ):
    """ [ ( comment, N x 3 array ) ] of a multi-frame xyz file. """
    frames = [ ]
    with open ( path ) as handle:
        lines = handle.read ( ).splitlines ( )
    p = 0
    while p < len ( lines ):
        if not lines[p].strip ( ): p += 1; continue
        n = int ( lines[p].split ( )[0] )
        comment = lines[p + 1]
        xyz = np.array ( [ [ float ( x ) for x in lines[p + 2 + a].split ( )[1:4] ] for a in range ( n ) ] )
        frames.append ( ( comment, xyz ) )
        p += n + 2
    return frames

def _energy_from_comment ( comment ):
    parts = comment.split ( )
    return float ( parts[parts.index ( "energy:" ) + 1] )

class XTBRunner:
    """ Runs xtb in a scratch folder; cancel() kills the running process
        (callable from another thread). """

    def __init__ ( self, xtb = None, gfn = 2, charge = 0, multiplicity = 1, solvent = None, threads = 4 ):
        from pdynamo.opls.ligand import find_xtb
        self.xtb = xtb or find_xtb ( )
        if self.xtb is None:
            raise RuntimeError ( "xtb was not found (set XTBHOME or put xtb in the PATH)." )
        self.gfn, self.charge, self.multiplicity, self.solvent = int ( gfn ), int ( charge ), int ( multiplicity ), solvent
        self.threads   = int ( threads )
        self.process   = None
        self.cancelled = False

    def cancel ( self ):
        self.cancelled = True
        if self.process is not None and self.process.poll ( ) is None:
            self.process.kill ( )

    def options ( self ):
        options = [ "--gfn", str ( self.gfn ), "--chrg", str ( self.charge ) ]
        if self.multiplicity > 1: options += [ "--uhf", str ( self.multiplicity - 1 ) ]
        if self.solvent: options += [ "--alpb", self.solvent ]
        return options

    def run ( self, arguments, cwd ):
        if self.cancelled: raise ScanCancelled ( )
        env = dict ( os.environ, OMP_NUM_THREADS = str ( self.threads ), OMP_STACKSIZE = os.environ.get ( "OMP_STACKSIZE", "1G" ) )
        self.process = subprocess.Popen ( [ self.xtb ] + arguments + self.options ( ), cwd = cwd, env = env,
                                          stdout = subprocess.PIPE, stderr = subprocess.STDOUT, text = True )
        output, _ = self.process.communicate ( )
        code = self.process.returncode
        self.process = None
        if self.cancelled: raise ScanCancelled ( )
        if code != 0 or "normal termination of xtb" not in output:
            raise RuntimeError ( "xtb failed:\n" + output[-1500:] )
        return output

    def single_point ( self, symbols, coordinates, scratch ):
        """ Energy in kJ/mol. """
        _write_xyz ( os.path.join ( scratch, "sp.xyz" ), symbols, coordinates )
        output = self.run ( [ "sp.xyz", "--sp" ], scratch )
        for line in output.splitlines ( ):
            if "TOTAL ENERGY" in line:
                return float ( line.split ( )[3] ) * HARTREE_TO_KJMOL
        raise RuntimeError ( "No TOTAL ENERGY in the xtb output." )

    def optimize ( self, symbols, coordinates, scratch ):
        """ Unconstrained xtb optimization; N x 3 array. """
        _write_xyz ( os.path.join ( scratch, "opt.xyz" ), symbols, coordinates )
        self.run ( [ "opt.xyz", "--opt" ], scratch )
        return read_xyz_frames ( os.path.join ( scratch, "xtbopt.xyz" ) )[-1][1]

    def relaxed_scan ( self, symbols, coordinates, quad, start, stop, n_points, scratch ):
        """ [ N x 3 arrays ] of the constrained optimizations. """
        _write_xyz ( os.path.join ( scratch, "start.xyz" ), symbols, coordinates )
        with open ( os.path.join ( scratch, "scan.inp" ), "w" ) as handle:
            handle.write ( "$constrain\n   force constant=1.0\n   dihedral: {},{},{},{},auto\n"
                           "$scan\n   1: {:.4f},{:.4f},{:d}\n$end\n".format ( *[ i + 1 for i in quad ], start, stop, n_points ) )
        for name in ( "xtbscan.log", "xtbrestart" ):
            if os.path.exists ( os.path.join ( scratch, name ) ): os.remove ( os.path.join ( scratch, name ) )
        self.run ( [ "start.xyz", "--opt", "--input", "scan.inp" ], scratch )
        frames = read_xyz_frames ( os.path.join ( scratch, "xtbscan.log" ) )
        if len ( frames ) != n_points:
            raise RuntimeError ( "xtb scan returned {} of {} points.".format ( len ( frames ), n_points ) )
        return [ xyz for ( comment, xyz ) in frames ]


def run_scan ( symbols, coordinates, bonds, quad, step = 15.0, relaxed = True, both_directions = False,
               runner = None, progress = None ):
    """ xTB torsion scan of quad = ( i, j, k, l ) over 360 degrees.
        Returns { "geometries": [ N x 3 ], "phi": [ deg ], "e_qm": [ kJ/mol, absolute ],
                  "settings": { ... } }. progress ( fraction, text ) is optional. """
    i, j, k, l = quad
    adjacency = adjacency_from_bonds ( len ( symbols ), bonds )
    side = moving_side ( adjacency, j, k )
    n = max ( 4, int ( round ( 360.0 / float ( step ) ) ) )
    step = 360.0 / n
    phi0 = dihedral_angle ( coordinates, quad )
    say = progress or ( lambda f, t: None )
    scratch = tempfile.mkdtemp ( prefix = "dyff_scan_" )
    try:
        geometries = [ ]
        if relaxed:
            say ( 0.0, "xtb relaxed scan ({} points)...".format ( n ) )
            forward = runner.relaxed_scan ( symbols, coordinates, quad, phi0, phi0 + 360.0 - step, n, scratch )
            geometries = list ( forward )
            if both_directions:
                say ( 0.25, "xtb relaxed scan, backward..." )
                backward = runner.relaxed_scan ( symbols, coordinates, quad, phi0, phi0 - 360.0 + step, n, scratch )
                # . backward[m] is at phi0 - m step = forward[( n - m ) % n]
                geometries = geometries + [ backward[( n - m ) % n] for m in range ( n ) ]
        else:
            # . rigid: rotation of the xTB-optimized structure (a non-optimized
            #   start would put its own strain into every point)
            say ( 0.0, "xtb optimization of the starting structure..." )
            start = runner.optimize ( symbols, coordinates, scratch )
            geometries = [ rotate_about_bond ( start, j, k, side, m * step ) for m in range ( n ) ]
        energies = [ ]
        for m, xyz in enumerate ( geometries ):
            say ( 0.5 + 0.5 * m / len ( geometries ), "xtb single points {}/{}...".format ( m + 1, len ( geometries ) ) )
            energies.append ( runner.single_point ( symbols, xyz, scratch ) )
        if relaxed and both_directions:
            # . the lower of the two runs at every angle (hysteresis)
            keep = [ m if energies[m] <= energies[m + n] else m + n for m in range ( n ) ]
            geometries = [ geometries[m] for m in keep ]
            energies = [ energies[m] for m in keep ]
    finally:
        shutil.rmtree ( scratch, ignore_errors = True )
    phi = [ dihedral_angle ( xyz, quad ) for xyz in geometries ]
    return { "geometries": geometries, "phi": phi, "e_qm": energies,
             "settings": { "quad": tuple ( quad ), "step": step, "relaxed": bool ( relaxed ),
                           "both_directions": bool ( both_directions ), "gfn": runner.gfn,
                           "charge": runner.charge, "multiplicity": runner.multiplicity, "solvent": runner.solvent } }


#=============================================================================
#  MM side
#=============================================================================
def mm_energies ( system, geometries ):
    """ pDynamo energies (kJ/mol) of system at every geometry. """
    saved = np.array ( [ [ system.coordinates3[a, c] for c in range ( 3 ) ] for a in range ( len ( system.atoms ) ) ] )
    energies = [ ]
    try:
        for xyz in geometries:
            for a in range ( len ( system.atoms ) ):
                for c in range ( 3 ): system.coordinates3[a, c] = float ( xyz[a][c] )
            energies.append ( float ( system.Energy ( log = None ) ) )
    finally:
        for a in range ( len ( system.atoms ) ):
            for c in range ( 3 ): system.coordinates3[a, c] = float ( saved[a, c] )
    return energies


#=============================================================================
#  Fit
#=============================================================================
def weights_of ( e_qm_relative ):
    e = np.asarray ( e_qm_relative, dtype = float )
    return np.where ( e <= WEIGHT_ONSET, 1.0, np.exp ( -( e - WEIGHT_ONSET ) / WEIGHT_KT ) )

def best_offset_rmse ( reference, model, weights = None ):
    """ RMSE between two curves after the (weighted) best constant offset;
        returns ( rmse, weighted rmse, offset ). """
    r, m = np.asarray ( reference, dtype = float ), np.asarray ( model, dtype = float )
    w = np.ones_like ( r ) if weights is None else np.asarray ( weights, dtype = float )
    offset = float ( np.sum ( w * ( r - m ) ) / np.sum ( w ) )
    d = r - m - offset
    return float ( np.sqrt ( np.mean ( d * d ) ) ), float ( np.sqrt ( np.sum ( w * d * d ) / np.sum ( w ) ) ), offset

def fit_periods ( geometries, e_qm, e_mm_zero, instances, max_period = 3, min_coefficient = 1.0e-3 ):
    """ Weighted linear least squares of the cosine series. Returns
        ( pairs [ ( c_n, n ) ], dropped periods, weights ). """
    e_qm = np.asarray ( e_qm, dtype = float ); e_mm_zero = np.asarray ( e_mm_zero, dtype = float )
    target = ( e_qm - e_qm.min ( ) ) - ( e_mm_zero - e_mm_zero.min ( ) )
    weights = weights_of ( e_qm - e_qm.min ( ) )
    angles = np.radians ( [ [ dihedral_angle ( xyz, quad ) for quad in instances ] for xyz in geometries ] )
    # . a period whose column barely changes along the scan (equivalent terms
    #   that cancel, e.g. n = 1, 3 for the two ortho carbons of a phenyl) is
    #   left out: fitting it would only fit the noise of the relaxed geometries
    spread, columns = { }, { }
    for n in range ( 1, int ( max_period ) + 1 ):
        column = np.cos ( n * angles ).sum ( axis = 1 )
        centred = column - np.average ( column, weights = weights )
        columns[n], spread[n] = column, float ( np.sqrt ( np.sum ( weights * centred * centred ) / np.sum ( weights ) ) )
    largest = max ( spread.values ( ) )
    periods = [ n for n in columns if largest > 1.0e-6 and spread[n] >= RELATIVE_SPREAD * largest ]
    dropped = [ n for n in columns if n not in periods ]
    columns = [ columns[n] for n in periods ]
    n_points = len ( geometries )
    if not periods:
        raise ValueError ( "None of the periods 1-{} changes along this scan (symmetric rotor).".format ( max_period ) )
    A = np.column_stack ( columns + [ np.ones ( n_points ) ] )
    sw = np.sqrt ( weights )
    solution, *_ = np.linalg.lstsq ( A * sw[:, None], target * sw, rcond = None )
    pairs = [ ( float ( c ), n ) for c, n in zip ( solution[:-1], periods ) if abs ( c ) >= min_coefficient ]
    if not pairs: pairs = list ( ZERO_TERM )
    return pairs, dropped, weights


def format_cosine_pairs ( pairs ):
    """ "c n; c n" (kJ/mol, E = sum c cos ( n phi )) -- the DYFF form. """
    return "; ".join ( "{:.3f} {:d}".format ( c, int ( n ) ) for c, n in pairs )


def fit_dihedral_key ( scan, key, instances, energies_for, current_text = "", format_value = format_cosine_pairs,
                       max_period = 3, force_field = "DYFF" ):
    """ Fit of the parameter key (all its terms in the scanned molecule:
        instances, indices of the scan geometries) to the scan.
        energies_for ( pairs or None ) -> MM energies (kJ/mol) at the scan
        geometries with every other parameter as now and this key replaced by
        pairs [ ( c_n, n ) ] (kJ/mol, E = sum c_n cos ( n phi ); None = as now).
        format_value ( pairs ) -> the fitted parameter in the force field's
        own notation. Returns the curves (relative to the QM minimum) and
        the statistics. """
    geometries = scan["geometries"]
    e_qm = np.asarray ( scan["e_qm"], dtype = float )
    e_before = np.asarray ( energies_for ( None ) )
    e_zero   = np.asarray ( energies_for ( ZERO_TERM ) )
    pairs, dropped, weights = fit_periods ( geometries, e_qm, e_zero, instances, max_period = max_period )
    e_after  = np.asarray ( energies_for ( pairs ) )
    qm = e_qm - e_qm.min ( )
    def aligned ( curve ):
        rmse, wrmse, offset = best_offset_rmse ( qm, curve, weights )
        return ( np.asarray ( curve ) + offset ), rmse, wrmse
    before, rmse_before, wrmse_before = aligned ( e_before )
    after,  rmse_after,  wrmse_after  = aligned ( e_after )
    j, k = scan["settings"]["quad"][1:3]
    around = [ q for q in instances if { q[1], q[2] } == { j, k } ]
    order = np.argsort ( scan["phi"] )
    sort = lambda v: [ float ( v[m] ) for m in order ]
    warnings = [ ]
    if len ( around ) < len ( instances ):
        warnings.append ( "{} of the {} terms with this key are on other bonds: the new parameter changes them too.".format (
                          len ( instances ) - len ( around ), len ( instances ) ) )
    if dropped:
        warnings.append ( "Period(s) {} left out: no change along this scan (equivalent terms cancel).".format (
                          ", ".join ( str ( n ) for n in dropped ) ) )
    if float ( np.min ( weights ) ) < 0.5:
        warnings.append ( "Points more than {:.0f} kJ/mol above the minimum were down-weighted (clashes).".format ( WEIGHT_ONSET ) )
    return { "key": key, "force_field": force_field, "pairs": pairs,
             "old_text": current_text, "new_text": format_value ( pairs ),
             "phi": sort ( scan["phi"] ), "qm": sort ( qm ), "before": sort ( before ), "after": sort ( after ),
             "target": sort ( qm - ( e_zero - e_zero.min ( ) ) ), "weights": sort ( weights ),
             "rmse_before": rmse_before, "rmse_after": rmse_after,
             "wrmse_before": wrmse_before, "wrmse_after": wrmse_after,
             "barrier_qm": float ( qm.max ( ) ), "n_terms": len ( instances ), "n_around": len ( around ),
             "dropped": dropped, "warnings": warnings, "settings": dict ( scan["settings"] ) }


def report_text ( result ):
    """ Plain-text summary (also saved with the plot data). """
    s = result["settings"]
    ff = result.get ( "force_field", "DYFF" )
    lines = [ "{} torsion fit: {}".format ( ff, result["key"] ),
              "Scan: dihedral {} ({}), step {:.1f} deg, GFN{}-xTB, charge {:+d}, multiplicity {}{}".format (
                  "-".join ( str ( i ) for i in s["quad"] ), "relaxed" + ( ", both directions" if s["both_directions"] else "" )
                  if s["relaxed"] else "rigid", s["step"], s["gfn"], s["charge"], s["multiplicity"],
                  ", ALPB " + s["solvent"] if s["solvent"] else "" ),
              "Terms with this key: {} ({} around the scanned bond)".format ( result["n_terms"], result["n_around"] ),
              "Old: " + result["old_text"],
              "New: " + result["new_text"],
              "Fitted cosine series (kJ/mol, E = sum c cos(n phi)): " + format_cosine_pairs ( result["pairs"] ),
              "RMSE vs xTB: {:.2f} -> {:.2f} kJ/mol (weighted {:.2f} -> {:.2f}); xTB barrier {:.2f} kJ/mol".format (
                  result["rmse_before"], result["rmse_after"], result["wrmse_before"], result["wrmse_after"], result["barrier_qm"] ) ]
    lines += [ "Note: " + w for w in result["warnings"] ]
    lines.append ( "" )
    lines.append ( "phi(deg)  E_xTB  E_{0}_before  E_{0}_after  weight   (kJ/mol)".format ( ff ) )
    for row in zip ( result["phi"], result["qm"], result["before"], result["after"], result["weights"] ):
        lines.append ( "{:8.2f} {:8.3f} {:8.3f} {:8.3f} {:6.3f}".format ( *row ) )
    return "\n".join ( lines )
