#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: per-bond torsion correction of DYFF fitted to a relaxed QM scan
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      The correction of a rotatable bond j-k is ONE cosine series that
#      replaces every DYFF dihedral term around j-k (dyff_edits
#      "dihedral_bonds"). Fit (shared by the DYFF Parameters window and
#      fraglib):
#
#        1. MM side: DYFF (with whatever other edits the molecule has, e.g.
#           Hessian bonds/angles) with the dihedrals around j-k removed;
#           restrained relaxed scan starting from every QM geometry (harmonic
#           restraint on the driven dihedral, conjugate gradient), energies
#           without the restraint -> E0 ( phi ), relaxed geometries;
#        2. target E_QM - E0; design sum_d cos ( n phi_d ) over every
#           quadruple i-j-k-l at the relaxed geometries; weighted linear least
#           squares (torsion_fit.fit_periods: periods 1..max_period, cancelling
#           periods dropped, high-energy points down-weighted);
#        3. self-consistency (ffTK / ForceBalance practice): the correction
#           changes how the rest of the molecule relaxes, so the scan is relaxed
#           again WITH the current correction, the correction's own energy is
#           subtracted at the new geometries (-> E0 at consistent geometries)
#           and the series is refitted; FIT_ITERATIONS rounds. Without this the
#           fitted residual looked good but the re-relaxed profile was often
#           worse than plain DYFF (tested on QCArchive pairs).
#
#      make_system ( edits ) -> pDynamo system with the DYFF model (those
#      edits) and an NB model; xyz in angstrom; energies in kJ/mol.
#

import contextlib
import io
import math
import threading

import numpy as np

HARTREE_TO_KJ = 2625.499639
RESTRAINT_FORCE = 1000.0          # kJ/mol/rad^2
MINIMIZATION_STEPS = 2000
FIT_ITERATIONS = 3

# . redirect_stdout swaps the global sys.stdout: fits running in several threads
#   (residue_library.scan_flexible_torsions, one fragment per bond) take this lock
STDOUT_LOCK = threading.RLock ( )


def bond_label ( j, k ):
    return "{}-{}".format ( min ( int ( j ), int ( k ) ), max ( int ( j ), int ( k ) ) )


def quads_about ( bonds, j, k ):
    """ Every quadruple i-j-k-l of the bond j-k (bonds: pairs of atom indices). """
    neighbours = { }
    for a, b in bonds:
        neighbours.setdefault ( int ( a ), set ( ) ).add ( int ( b ) ); neighbours.setdefault ( int ( b ), set ( ) ).add ( int ( a ) )
    return [ ( i, j, k, l ) for i in neighbours.get ( j, ( ) ) if i != k for l in neighbours.get ( k, ( ) ) if l != j and l != i ]


def dihedral ( x, q ):
    b0, b1, b2 = x[q[0]] - x[q[1]], x[q[2]] - x[q[1]], x[q[3]] - x[q[2]]
    b1 = b1 / np.linalg.norm ( b1 )
    v = b0 - np.dot ( b0, b1 ) * b1; w = b2 - np.dot ( b2, b1 ) * b1
    return math.atan2 ( float ( np.dot ( np.cross ( b1, v ), w ) ), float ( np.dot ( v, w ) ) )


def relaxed_scan ( system, geometries, quad, targets ):
    """ ( energies kJ/mol, relaxed geometries ): restrained minimization at every
        point, starting from the given geometry; energy without the restraint. """
    from pMolecule import RestraintModel, RestraintDihedral, RestraintEnergyModel
    from pScientific.Geometry3 import Coordinates3
    from pSimulation import LBFGSMinimize_SystemGeometry
    n = len ( system.atoms )
    energies, relaxed = [ ], [ ]
    with STDOUT_LOCK, contextlib.redirect_stdout ( io.StringIO ( ) ):
        for xyz, target in zip ( geometries, targets ):
            c3 = Coordinates3.WithExtent ( n )
            for a in range ( n ):
                for c in range ( 3 ): c3[a, c] = float ( xyz[a][c] )
            system.coordinates3 = c3
            restraints = RestraintModel ( )
            system.DefineRestraintModel ( restraints )
            restraints["torsion"] = RestraintDihedral.WithOptions (
                energyModel = RestraintEnergyModel.Harmonic ( float ( target ), RESTRAINT_FORCE, period = 360.0 ),
                point1 = int ( quad[0] ), point2 = int ( quad[1] ), point3 = int ( quad[2] ), point4 = int ( quad[3] ) )
            LBFGSMinimize_SystemGeometry ( system, log = None, maximumIterations = MINIMIZATION_STEPS,
                                           rmsGradientTolerance = 0.01 )
            system.DefineRestraintModel ( None )
            energies.append ( float ( system.Energy ( log = None ) ) )
            relaxed.append ( np.array ( [ [ system.coordinates3[a, c] for c in range ( 3 ) ] for a in range ( n ) ] ) )
    return np.array ( energies ), relaxed


def fit_bond ( make_system, base_edits, bonds, quad, geometries, targets, e_qm, max_period = 4, ridge = 0.0 ):
    """ Correction of the bond quad[1]-quad[2] fitted to a relaxed QM scan
        ( geometries A, targets degrees, e_qm kJ/mol ). Returns ( pairs [ ( c, n ) ],
        in-sample RMSE kJ/mol (best offset, weighted points), dropped periods ). """
    from pdynamo.torsion_fit import fit_periods, best_offset_rmse
    j, k = int ( quad[1] ), int ( quad[2] )
    quads = quads_about ( bonds, j, k )
    e_qm = np.asarray ( e_qm, dtype = float )
    correction = lambda x, pairs: sum ( c * math.cos ( n * dihedral ( x, q ) ) for q in quads for c, n in pairs )
    pairs, start = [ ], geometries
    for iteration in range ( FIT_ITERATIONS ):
        system = make_system ( _with_bond ( base_edits, j, k, pairs ) )
        e_total, relaxed = relaxed_scan ( system, start, quad, targets )
        e0 = e_total - np.array ( [ correction ( x, pairs ) for x in relaxed ] )     # MM without the bond's torsion
        pairs, dropped, weights = fit_periods ( relaxed, e_qm, e0, quads, max_period = max_period, ridge = ridge )
        start = relaxed                                                              # continue from the relaxed geometries
    predicted = e0 + np.array ( [ correction ( x, pairs ) for x in relaxed ] )
    rmse, wrmse, _ = best_offset_rmse ( e_qm - e_qm.min ( ), predicted - predicted.min ( ), weights )
    return pairs, wrmse, dropped


def _with_bond ( base_edits, j, k, pairs ):
    import copy
    edits = copy.deepcopy ( base_edits or { } )
    edits.setdefault ( "dihedral_bonds", { } )[bond_label ( j, k )] = [ list ( p ) for p in pairs ]
    return edits


def profile_error ( make_system, edits, quad, geometries, targets, e_ref ):
    """ RMSE (kJ/mol, all points, relative to each profile's minimum) and low-energy
        RMSE (reference < 20.9 kJ/mol = 5 kcal/mol) of a relaxed DYFF scan with these
        edits against a reference profile; plus min_ok (the MM minimum lies within
        2.09 kJ/mol = 0.5 kcal/mol of the reference minimum, on the reference). """
    system = make_system ( edits )
    e_mm, _ = relaxed_scan ( system, geometries, quad, targets )
    ref = np.asarray ( e_ref, dtype = float ); ref = ref - ref.min ( )
    mm = e_mm - e_mm.min ( )
    low = ref < 20.92
    return { "rmse": float ( np.sqrt ( np.mean ( ( mm - ref ) ** 2 ) ) ),
             "rmse_low": float ( np.sqrt ( np.mean ( ( mm[low] - ref[low] ) ** 2 ) ) ) if low.any ( ) else float ( "nan" ),
             "min_ok": bool ( ref[int ( np.argmin ( mm ) )] <= 2.092 ), "mm": mm.tolist ( ) }
