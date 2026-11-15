#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: bond and angle force constants from a QM Hessian (modified Seminario)
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      Modified Seminario method (J. M. Seminario, Int. J. Quantum Chem. 1996;
#      A. E. A. Allen, M. C. Payne, D. J. Cole, J. Chem. Theory Comput. 14, 274,
#      2018): the 3x3 interatomic block -H_AB of the Cartesian Hessian at a
#      minimum, projected on the A-B direction, gives the bond force constant;
#      projected on the in-plane perpendiculars of an angle A-B-C (with the
#      scaling of Allen et al. for the other angles that share each bond),
#      the angle force constant. Only the Hessian, the coordinates and the
#      list of bonds are needed: no atom typing, no RDKit -- metal centres
#      are treated like any other atom (the method behind AMBER's MCPB.py).
#
#      Conventions: E = 1/2 k ( r - r0 )^2 and E = 1/2 k ( theta - theta0 )^2;
#      k in kJ/mol/A^2 and kJ/mol/rad^2; r0 in A, theta0 in radians (the QM
#      geometry). Hessians in kJ/mol/A^2 (read_xtb_hessian converts xtb's
#      hartree/bohr^2).
#
#      DYFF forms (pDynamo, docs/dyff_to_namd_methodology.md):
#        bond  E = k ( b - b0 )^2                -> k_DYFF = k / 2
#        angle E = sum_n c_n cos ( n theta ), the UFF general form
#              K [ C0 + C1 cos theta + C2 cos 2 theta ], C2 = 1 / ( 4 sin^2 theta0 ),
#              C1 = -4 C2 cos theta0, C0 = C2 ( 2 cos^2 theta0 + 1 ): its
#              curvature at theta0 is exactly K -> K = k; near-linear angles
#              use K ( 1 + cos theta ) (curvature K at 180 degrees).
#

import math
import os
import shutil
import subprocess
import tempfile
from collections import defaultdict

import numpy as np

HARTREE_PER_BOHR2_TO_KJ_PER_A2 = 2625.499639 / 0.529177210903 ** 2
LINEAR_ANGLE = math.radians ( 175.0 )
_MASS = { "H": 1.008, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Si": 28.085, "P": 30.974,
          "S": 32.06, "Cl": 35.45, "Br": 79.904, "I": 126.90, "Fe": 55.845, "Cu": 63.546, "Zn": 65.38, "Co": 58.933,
          "Ni": 58.693, "Mn": 54.938, "Mg": 24.305, "Ca": 40.078, "Na": 22.990, "K": 39.098, "Ru": 101.07,
          "Pd": 106.42, "Pt": 195.08, "Rh": 102.91, "Ir": 192.22, "Cr": 51.996, "V": 50.942, "Ti": 47.867, "Mo": 95.95 }


#=============================================================================
#  xtb
#=============================================================================
def read_xtb_hessian ( path, n_atoms ):
    """ Turbomole-format `hessian` file (hartree/bohr^2) -> ( 3N, 3N ) in kJ/mol/A^2. """
    values = [ ]
    with open ( path ) as handle:
        for line in handle:
            if line.startswith ( "$" ): continue
            values += [ float ( x ) for x in line.split ( ) ]
    n = 3 * n_atoms
    if len ( values ) < n * n:
        raise ValueError ( "hessian file has {} values, {} expected".format ( len ( values ), n * n ) )
    hessian = np.array ( values[:n * n] ).reshape ( n, n )
    return 0.5 * ( hessian + hessian.T ) * HARTREE_PER_BOHR2_TO_KJ_PER_A2


def xtb_hessian ( symbols, coordinates, charge = 0, multiplicity = 1, gfn = 2, xtb = None, threads = 2, solvent = None, log_dir = None ):
    """ xtb --ohess: ( optimized coordinates A, Hessian kJ/mol/A^2, frequencies cm-1 from xtb ).
        gfn = 0/1/2 (GFN-xTB) or "gxtb" (g-xTB through tblite, xtb >= 6.7.1; gas phase
        only: tblite has no implicit solvent, so `solvent` is ignored with g-xTB). """
    if xtb is None:
        from glob import glob
        candidates = [ shutil.which ( "xtb" ) ] + sorted ( glob ( os.path.expanduser ( "~/programs/xtb-*/bin/xtb" ) ), reverse = True )
        xtb = next ( ( c for c in candidates if c and os.access ( c, os.X_OK ) ), None )
    if xtb is None:
        raise RuntimeError ( "xtb was not found" )
    tmp = tempfile.mkdtemp ( prefix = "seminario_" )
    try:
        with open ( os.path.join ( tmp, "m.xyz" ), "w" ) as handle:
            handle.write ( "{}\n\n".format ( len ( symbols ) ) )
            for s, p in zip ( symbols, coordinates ):
                handle.write ( "{} {:.8f} {:.8f} {:.8f}\n".format ( s, *p ) )
        method = [ "--gxtb" ] if str ( gfn ).lower ( ) in ( "gxtb", "g-xtb" ) else [ "--gfn", str ( gfn ) ]
        command = [ xtb, "m.xyz", "--ohess" ] + method + [ "--chrg", str ( charge ) ]
        if multiplicity > 1: command += [ "--uhf", str ( multiplicity - 1 ) ]
        if solvent and method[0] != "--gxtb": command += [ "--alpb", solvent ]
        result = subprocess.run ( command, cwd = tmp, capture_output = True, text = True, timeout = 7200,
                                  env = dict ( os.environ, OMP_NUM_THREADS = str ( threads ), OMP_STACKSIZE = "1G" ) )
        from pdynamo import param_logs
        param_logs.keep ( tmp, log_dir, [ "m.xyz", "xtbopt.xyz", "vibspectrum" ],
                          { "xtb_ohess.out": "$ " + " ".join ( command[1:] ) + "\n" + result.stdout + result.stderr } )
        if "normal termination of xtb" not in result.stdout + result.stderr:
            raise RuntimeError ( "xtb --ohess failed: " + ( result.stdout + result.stderr )[-800:] )
        lines = open ( os.path.join ( tmp, "xtbopt.xyz" ) ).read ( ).splitlines ( )
        xyz = np.array ( [ [ float ( x ) for x in l.split ( )[1:4] ] for l in lines[2:2 + len ( symbols )] ] )
        hessian = read_xtb_hessian ( os.path.join ( tmp, "hessian" ), len ( symbols ) )
        freqs = [ ]
        text = result.stdout
        start = text.find ( "projected vibrational frequencies" )
        for line in text[start:].splitlines ( )[1:]:
            if not line.startswith ( "eigval" ): break
            freqs += [ float ( x ) for x in line.split ( ":" )[1].split ( ) ]
        return xyz, hessian, np.array ( freqs )
    finally:
        shutil.rmtree ( tmp, ignore_errors = True )


#=============================================================================
#  Frequencies (validation)
#=============================================================================
def frequencies ( hessian, symbols, masses = None, n_remove = 6 ):
    """ Harmonic frequencies (cm-1) of a Hessian in kJ/mol/A^2, the n_remove
        smallest |values| (translations/rotations) left out; imaginary as negative. """
    m = np.array ( masses if masses is not None else [ _MASS.get ( s, 12.0 ) for s in symbols ] )
    w = 1.0 / np.sqrt ( np.repeat ( m, 3 ) )
    eig = np.linalg.eigvalsh ( hessian * w[:, None] * w[None, :] )
    factor = 1.0e26                                               # kJ mol-1 A-2 amu-1 -> s-2
    nu = np.sign ( eig ) * np.sqrt ( np.abs ( eig ) * factor ) / ( 2.0 * math.pi * 2.99792458e10 )
    order = np.argsort ( np.abs ( nu ) )
    return np.sort ( nu[order[n_remove:]] )


#=============================================================================
#  Modified Seminario
#=============================================================================
def _block ( hessian, a, b ):
    return -hessian[3 * a:3 * a + 3, 3 * b:3 * b + 3]


def _eig ( block ):
    values, vectors = np.linalg.eig ( block )
    return np.real ( values ), np.real ( vectors )


def bond_constants ( hessian, coordinates, bonds ):
    """ { ( a, b ): ( r0 A, k kJ/mol/A^2 ) } -- average of the A-B and B-A projections. """
    x = np.asarray ( coordinates, dtype = float )
    out = { }
    for a, b in bonds:
        a, b = int ( a ), int ( b )
        r = x[b] - x[a]
        r0 = float ( np.linalg.norm ( r ) )
        u = r / r0
        ks = [ ]
        for p, q, v in ( ( a, b, u ), ( b, a, -u ) ):
            lam, vec = _eig ( _block ( hessian, p, q ) )
            ks.append ( sum ( lam[i] * abs ( np.dot ( v, vec[:, i] ) ) for i in range ( 3 ) ) )
        out[( min ( a, b ), max ( a, b ) )] = ( r0, float ( 0.5 * ( ks[0] + ks[1] ) ) )
    return out


def _angles ( bonds ):
    neighbours = defaultdict ( set )
    for a, b in bonds:
        neighbours[int ( a )].add ( int ( b ) ); neighbours[int ( b )].add ( int ( a ) )
    angles = [ ]
    for centre, ns in neighbours.items ( ):
        ns = sorted ( ns )
        for i in range ( len ( ns ) ):
            for j in range ( i + 1, len ( ns ) ):
                angles.append ( ( ns[i], centre, ns[j] ) )
    return angles, neighbours


def angle_constants ( hessian, coordinates, bonds ):
    """ { ( a, b, c ): ( theta0 rad, k kJ/mol/rad^2 ) }, b central, a < c (modified Seminario). """
    x = np.asarray ( coordinates, dtype = float )
    angles, neighbours = _angles ( bonds )

    def unit ( v ): return v / np.linalg.norm ( v )

    def perpendicular ( a, b, c ):
        """ in-plane unit vector at a, perpendicular to a-b, pointing into the angle a-b-c. """
        u_ab, u_cb = unit ( x[a] - x[b] ), unit ( x[c] - x[b] )
        normal = np.cross ( u_cb, u_ab )
        if np.linalg.norm ( normal ) < 1.0e-6:                         # linear: any perpendicular
            trial = np.array ( [ 1.0, 0.0, 0.0 ] ) if abs ( u_ab[0] ) < 0.9 else np.array ( [ 0.0, 1.0, 0.0 ] )
            normal = np.cross ( u_ab, trial )
        normal = unit ( normal )
        return unit ( np.cross ( normal, u_ab ) )

    # . in-plane perpendiculars of every ( terminal, centre, other ) for the scaling factors
    perps = { }
    for a, b, c in angles:
        perps[( a, b, c )] = perpendicular ( a, b, c )
        perps[( c, b, a )] = perpendicular ( c, b, a )
    out = { }
    for a, b, c in angles:
        theta0 = math.acos ( max ( -1.0, min ( 1.0, float ( np.dot ( unit ( x[a] - x[b] ), unit ( x[c] - x[b] ) ) ) ) ) )
        terms = [ ]
        for p, q in ( ( a, c ), ( c, a ) ):
            u_p = perps[( p, b, q )]
            lam, vec = _eig ( _block ( hessian, p, b ) )
            proj = sum ( lam[i] * abs ( np.dot ( u_p, vec[:, i] ) ) for i in range ( 3 ) )
            r2 = float ( np.dot ( x[p] - x[b], x[p] - x[b] ) )
            others = [ perps[( p, b, o )] for o in neighbours[b] if o not in ( p, q ) ]
            scaling = 1.0 + ( sum ( abs ( np.dot ( u_p, w ) ) ** 2 for w in others ) / len ( others ) if others else 0.0 )
            terms.append ( scaling / ( r2 * proj ) if proj > 0.0 else float ( "inf" ) )
        k = 1.0 / sum ( terms ) if all ( math.isfinite ( t ) for t in terms ) and sum ( terms ) > 0 else 0.0
        out[( min ( a, c ), b, max ( a, c ) )] = ( theta0, float ( k ) )
    return out


#=============================================================================
#  DYFF parameters
#=============================================================================
def dyff_bond ( r0, k ):
    """ ( b0, k ) for DYFF's E = k ( b - b0 )^2. """
    return float ( r0 ), float ( k ) / 2.0


def dyff_angle_pairs ( theta0, k ):
    """ [ ( c_n, n ) ] for DYFF's E = sum c_n cos ( n theta ) with minimum at theta0
        and curvature k there. """
    if theta0 >= LINEAR_ANGLE:
        return [ ( float ( k ), 0 ), ( float ( k ), 1 ) ]
    s2 = math.sin ( theta0 ) ** 2
    c2 = 1.0 / ( 4.0 * s2 )
    c1 = -4.0 * c2 * math.cos ( theta0 )
    c0 = c2 * ( 2.0 * math.cos ( theta0 ) ** 2 + 1.0 )
    return [ ( float ( k * c0 ), 0 ), ( float ( k * c1 ), 1 ), ( float ( k * c2 ), 2 ) ]


# . a Seminario constant below these is not physical for a bonded term (e.g. the
#   C-COO- bond of a gas-phase carboxylate in GFN2-xTB comes out NEGATIVE: the
#   isolated anion is close to decarboxylation); such terms keep the DYFF value
MIN_BOND_K  = 100.0     # kJ/mol/A^2 (1/2 k convention)
MIN_ANGLE_K = 10.0      # kJ/mol/rad^2


def dyff_edits_from_hessian ( hessian, coordinates, bonds, skipped = None, calibration = None ):
    """ DYFF per-term edits (pdynamo/dyff_edits.py: "bond_terms", "angle_terms")
        from a Hessian at a minimum. Terms with a non-physical constant (below
        MIN_BOND_K / MIN_ANGLE_K) are left to DYFF and listed in `skipped`
        (a list, if given) as ( kind, atoms, k ). calibration ( kind, atoms ) ->
        ( shift, ratio, level ) or None (pdynamo/bond_library.HessianCalibration):
        r0/theta0 + shift, k * ratio, applied before the guard. """
    edits = { "bond_terms": { }, "angle_terms": { } }
    def calibrate ( kind, atoms, v0, k ):
        c = calibration ( kind, atoms ) if calibration is not None else None
        return ( v0 + c[0], k * c[1] ) if c else ( v0, k )
    for ( a, b ), ( r0, k ) in bond_constants ( hessian, coordinates, bonds ).items ( ):
        r0, k = calibrate ( "bond", ( a, b ), r0, k )
        if k < MIN_BOND_K:
            if skipped is not None: skipped.append ( ( "bond", ( a, b ), k ) )
            continue
        edits["bond_terms"]["{}-{}".format ( a, b )] = list ( dyff_bond ( r0, k ) )
    for ( a, b, c ), ( theta0, k ) in angle_constants ( hessian, coordinates, bonds ).items ( ):
        theta0, k = calibrate ( "angle", ( a, b, c ), theta0, k )
        if k < MIN_ANGLE_K:
            if skipped is not None: skipped.append ( ( "angle", ( a, b, c ), k ) )
            continue
        edits["angle_terms"]["{}-{}-{}".format ( a, b, c )] = [ list ( p ) for p in dyff_angle_pairs ( theta0, k ) ]
    return edits


#=============================================================================
#  Out-of-plane stiffness from the Hessian (2026-10-07)
#=============================================================================
#  Seminario has no out-of-plane analogue, and the pyramidalization stiffness
#  of a centre is shared by DYFF's out-of-plane terms AND the dihedrals about
#  the centre's bonds. The probe is the local curvature along d_c, the central
#  atom j moved along the normal of its three neighbours' plane:
#
#      k_c = d_c^T H d_c          (kJ/mol/A^2)
#
#  For a planar centre this displacement changes no bond length or angle to
#  first order, so it measures exactly the out-of-plane + torsional stiffness
#  (a pyramidal group-15 centre also feels its angles, already Seminario-fitted).
#  The MM curvature is linear in per-centre scales of DYFF's series,
#
#      k_c ( s ) = r_c + sum_m s_m A_cm ,  A_cm = d_c^T H_m d_c
#
#  (r_c: every other term, H_m: centre m's default out-of-plane terms; centres
#  sharing atoms couple), solved for k_c ( s ) = k_c ( QM ) with the scales
#  bounded. The functional form and equilibrium of DYFF's series are kept, only
#  their strength changes. Result: per-centre edits "outofplane_centres"
#  (dyff_edits.py). A relaxed-compliance target (b^T H+ b) was tried first and
#  dropped: the MM Hessian at the QM minimum is indefinite, its pseudoinverse
#  gave negative compliances.
OOP_SCALE_RANGE = ( 0.05, 25.0 )
OOP_STEP        = 0.01      # A, finite-difference step of the curvature

def _oop_cos ( p, i, j, k, l ):
    n = np.cross ( p[k] - p[j], p[l] - p[j] )
    u = p[i] - p[j]
    return float ( np.dot ( n, u ) / ( np.linalg.norm ( n ) * np.linalg.norm ( u ) ) )

def _oop_centre_energy ( p, j, neighbours, series ):
    """ The three DYFF terms of centre j (pDynamo CosineOutOfPlane convention). """
    a, b, c = neighbours
    e = 0.0
    for i, k, l in ( ( a, b, c ), ( b, a, c ), ( c, a, b ) ):
        phi = math.acos ( max ( -1.0, min ( 1.0, _oop_cos ( p, i, j, k, l ) ) ) )
        e  += sum ( fc * math.cos ( n * phi ) for fc, n in series )
    return e

def _plane_normal_displacement ( coordinates, j, neighbours ):
    x = np.asarray ( coordinates, dtype = float )
    a, b, c = neighbours
    n = np.cross ( x[b] - x[a], x[c] - x[a] )
    d = np.zeros_like ( x )
    d[j] = n / np.linalg.norm ( n )
    return d

def _curvature ( energy, coordinates, d, step = OOP_STEP ):
    x = np.asarray ( coordinates, dtype = float )
    return ( energy ( x + step * d ) - 2.0 * energy ( x ) + energy ( x - step * d ) ) / step ** 2

def oop_centres ( system ):
    """ { j: ( ( a, b, c ), series ) } of the DYFF out-of-plane terms of a system
        (the three terms of a centre share the series). """
    centres = { }
    for container in system.mmState.mmTerms:
        if type ( container ).__name__ != "CosineOutOfPlaneContainer": continue
        state = container.__getstate__ ( )
        for ( i, j, k, l ), parameter, active in state["terms"]:
            if not active: continue
            series = [ ( float ( fc ), int ( n ) ) for fc, n in state["parameters"][parameter] ]
            centres.setdefault ( int ( j ), ( tuple ( sorted ( ( int ( i ), int ( k ), int ( l ) ) ) ), series ) )
    return centres

def _system_energy ( system ):
    import contextlib, io
    from pScientific.Geometry3 import Coordinates3
    def energy ( x ):
        c3 = Coordinates3.WithExtent ( len ( x ) )
        for a in range ( len ( x ) ):
            for k in range ( 3 ): c3[a, k] = float ( x[a, k] )
        system.coordinates3 = c3
        with contextlib.redirect_stdout ( io.StringIO ( ) ):
            return float ( system.Energy ( log = None ) )
    return energy

def outofplane_edits ( hessian, coordinates, make_system, base_edits = None, report = None ):
    """ Per-centre out-of-plane edits ( { "outofplane_centres": { "j": series } } )
        matching the out-of-plane curvature of a QM Hessian (kJ/mol/A^2, at the
        QM minimum `coordinates`). make_system ( edits ) returns the DYFF pDynamo
        system; base_edits (e.g. the Seminario bonds and angles) are kept.
        `report` (a dict, if given) receives per centre { "k_qm", "k_dyff",
        "k_fit", "scale" } (kJ/mol/A^2). """
    base = dict ( base_edits or { } )
    centres = oop_centres ( make_system ( base ) )
    if not centres: return { "outofplane_centres": { } }
    x = np.asarray ( coordinates, dtype = float )
    h = np.asarray ( hessian, dtype = float )
    labels = sorted ( centres )
    rest = _system_energy ( make_system ( dict ( base, outofplane_centres = { str ( j ): [ ] for j in labels } ) ) )
    full = _system_energy ( make_system ( base ) )
    d = [ _plane_normal_displacement ( x, j, centres[j][0] ) for j in labels ]
    k_qm = np.array ( [ float ( v.ravel ( ) @ h @ v.ravel ( ) ) for v in d ] )
    r = np.array ( [ _curvature ( rest, x, v ) for v in d ] )
    k_dyff = np.array ( [ _curvature ( full, x, v ) for v in d ] )
    a = np.zeros ( ( len ( labels ), len ( labels ) ) )
    for m, j in enumerate ( labels ):
        neighbours, series = centres[j]
        atoms = { j, *neighbours }
        energy = lambda p, j = j, n = neighbours, s = series: _oop_centre_energy ( p, j, n, s )
        for c, v in enumerate ( d ):
            if labels[c] in atoms: a[c, m] = _curvature ( energy, x, v )
    lo, hi = OOP_SCALE_RANGE
    scales = np.ones ( len ( labels ) )
    for _ in range ( 200 ):                                     # projected Gauss-Seidel
        previous = scales.copy ( )
        for c in range ( len ( labels ) ):
            if a[c, c] <= 0.0: continue
            other = r[c] + a[c] @ scales - a[c, c] * scales[c]
            scales[c] = min ( hi, max ( lo, ( k_qm[c] - other ) / a[c, c] ) )
        if np.max ( np.abs ( scales - previous ) ) < 1.0e-6: break
    k_fit = r + a @ scales
    edits = { }
    for c, j in enumerate ( labels ):
        ok = k_qm[c] > 0.0 and a[c, c] > 0.0
        if report is not None:
            report[j] = { "k_qm": float ( k_qm[c] ), "k_dyff": float ( k_dyff[c] ), "k_fit": float ( k_fit[c] ),
                          "scale": float ( scales[c] ) if ok else None }
        if ok: edits[str ( j )] = [ [ float ( scales[c] ) * fc, n ] for fc, n in centres[j][1] ]
    return { "outofplane_centres": edits }


def mm_hessian ( system, coordinates, step = 0.005 ):
    """ Cartesian Hessian (kJ/mol/A^2) of a pDynamo system by central differences
        of the analytic gradient. """
    import contextlib, io
    from pScientific.Geometry3 import Coordinates3
    x = np.array ( coordinates, dtype = float )
    n = len ( x )
    def gradient ( y ):
        c3 = Coordinates3.WithExtent ( n )
        for a in range ( n ):
            for k in range ( 3 ): c3[a, k] = float ( y[a, k] )
        system.coordinates3 = c3
        with contextlib.redirect_stdout ( io.StringIO ( ) ):
            system.Energy ( log = None, doGradients = True )
        g = system.scratch.gradients3
        return np.array ( [ g[a, k] for a in range ( n ) for k in range ( 3 ) ] )
    h = np.zeros ( ( 3 * n, 3 * n ) )
    for d in range ( 3 * n ):
        a, k = divmod ( d, 3 )
        x[a, k] += step; gp = gradient ( x )
        x[a, k] -= 2 * step; gm = gradient ( x )
        x[a, k] += step
        h[d] = ( gp - gm ) / ( 2 * step )
    gradient ( x )
    return 0.5 * ( h + h.T )


#  DFT statistics instead of a per-molecule fit: GFN2-xTB is 10-15 % softer than
#  B3LYP-D3BJ/DZVP for pyramidalization (median k_xtb / k_dft 0.85-0.92 over 300
#  QCArchive Hessians), and since the dihedrals carry most of the stiffness that
#  small bias becomes a large error in the fitted scale (C: 0.16 from xTB, 0.69
#  from DFT), which made the low frequencies worse. The median DFT scales per
#  centre kind (fitted on the 300 smallest molecules of the set, after the
#  Seminario bonds/angles) are used instead. N (amide/aniline) and sp3 group 15
#  hit the lower bound: their excess stiffness is in the dihedrals (torsion fit).
OOP_DFT_SCALES = { "C": 0.69, "C=O": 0.59, "N": 0.05, "N=O": 0.59, "Nsp3": 0.05 }

def oop_kind ( symbols, j, series ):
    """ "C", "N", "C=O"/"N=O" (DYFF's sp2 centre next to an sp2 O) or "<element>sp3". """
    kind = symbols[j]
    if len ( series ) == 3: return kind + "sp3"
    if series and abs ( series[0][0] ) > 30.0: return kind + "=O"
    return kind

def outofplane_statistical_edits ( system, symbols, scales = None ):
    """ Per-centre edits scaling DYFF's out-of-plane series by the DFT median
        of the centre kind (OOP_DFT_SCALES). """
    scales = OOP_DFT_SCALES if scales is None else scales
    edits = { }
    for j, ( _, series ) in oop_centres ( system ).items ( ):
        s = scales.get ( oop_kind ( symbols, j, series ) )
        if s is not None: edits[str ( j )] = [ [ s * fc, n ] for fc, n in series ]
    return { "outofplane_centres": edits }
