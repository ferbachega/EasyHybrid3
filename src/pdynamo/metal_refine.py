#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  metal_refine.py
#
#  Fast, direct per-system refinement of the metal site of a DYFF model (2026-10-13; DYFF_validation,
#  tmQM study, "CM5 refinement" step 4b).
#
#  Goal: the DYFF (+ charges, normally CM5) optimised geometry of a complex reproduces a reference
#  geometry of the metal site (GFN2-xTB optimum, DFT or X-ray), with a few DYFF optimisations and no
#  Hessian. Every metal-centred bond M-L and angle L-M-L becomes a per-instance term
#  (dyff_edits "bond_terms" / "angle_terms"):
#      start   r0 = d_ref, theta0 = theta_ref
#      repeat  optimise; r0 += d_ref - d_opt; theta0 += theta_ref - theta_opt
#  The force constants keep the UFF formulas of the assigned types (bond: E2A zM zL / r0^3; angle:
#  general form with the new theta0, linear form above 175 deg). What it corrects: trans influence,
#  Jahn-Teller and chelate distortions, multiple-bond lengths -- what one radius per metal type cannot.
#  What it does not do: it does not touch the ligands' own terms (DYFF + bond/torsion libraries) nor
#  the charges.
#
#  refine ( build_system, reference, metals = None, iterations = 3 ) -> ( edits, report )
#    build_system ( edits ) -> a pDynamo System with the DYFF model (EH set) built with these edits,
#      its charges set and its coordinates at the reference (the caller owns charges/NB model);
#    reference: Nx3 coordinates (Angstrom).

import math
import numpy as np

STEPS   = 2000
TOL_KJ  = 0.04184                         # . 0.01 kcal/mol/A rms gradient
LINEAR  = 175.0


def _angle ( x, i, j, k ):
    u = x[i] - x[j]; v = x[k] - x[j]
    c = float ( np.dot ( u, v ) / ( np.linalg.norm ( u ) * np.linalg.norm ( v ) ) )
    return math.degrees ( math.acos ( max ( -1.0, min ( 1.0, c ) ) ) )


def metal_terms ( system, metals = None ):
    """ ( bonds [(m, l)], angles [(l1, m, l2)] ) of the metal centres of a system. """
    from pdynamo.metal_geometry import TRANSITION_METALS
    nb = { }
    for i, j in system.connectivity.bondIndices:
        nb.setdefault ( i, set ( ) ).add ( j ); nb.setdefault ( j, set ( ) ).add ( i )
    if metals is None: metals = [ i for i, a in enumerate ( system.atoms ) if a.atomicNumber in TRANSITION_METALS ]
    bonds, angles = [ ], [ ]
    for m in metals:
        ligands = sorted ( nb.get ( m, ( ) ) )
        bonds  += [ ( m, l ) for l in ligands ]
        angles += [ ( ligands[a], m, ligands[b] ) for a in range ( len ( ligands ) ) for b in range ( a + 1, len ( ligands ) ) ]
    return bonds, angles


def _model_types ( system ):
    st = system.mmState
    labels = [ st.atomTypes[i] for i in st.atomTypeIndices ]
    model = system.mmModel
    atom_types = getattr ( model, "_eh_atom_types", None )
    return labels, atom_types


def _bond_k ( atom_types, ti, tj, r0 ):
    from pScientific import Units
    zi = atom_types.GetItem ( ti ).GetProperty ( "Effective Charge" )
    zj = atom_types.GetItem ( tj ).GetProperty ( "Effective Charge" )
    return Units.Energy_E2Angstroms_To_Kilojoules_Per_Mole * zi * zj / r0 ** 3


def _angle_terms ( atom_types, ti, tj, tk, r_ij, r_jk, theta_deg ):
    """ UFF angle with natural angle theta (general form; linear form above LINEAR). """
    from pScientific import Units
    theta = math.radians ( min ( theta_deg, 179.9 ) )
    zi = atom_types.GetItem ( ti ).GetProperty ( "Effective Charge" )
    zk = atom_types.GetItem ( tk ).GetProperty ( "Effective Charge" )
    sinT, cosT = math.sin ( theta ), math.cos ( theta )
    r_ik = math.sqrt ( r_ij ** 2 + r_jk ** 2 - 2.0 * r_ij * r_jk * cosT )
    fC0 = 2.0 * Units.Energy_E2Angstroms_To_Kilojoules_Per_Mole * zi * zk * ( 3.0 * r_ij * r_jk * ( 1.0 - cosT ** 2 ) - r_ik ** 2 * cosT ) / r_ik ** 5
    if theta_deg >= LINEAR:
        return [ [ fC0, 0 ], [ fC0, 1 ] ]
    fC0 /= ( 4.0 * sinT ** 2 )
    return [ [ fC0 * ( 1.0 + 2.0 * cosT ** 2 ), 0 ], [ - fC0 * 4.0 * cosT, 1 ], [ fC0, 2 ] ]


def _optimise ( system, steps = STEPS ):
    import contextlib, io
    from pSimulation import LBFGSMinimize_SystemGeometry
    with contextlib.redirect_stdout ( io.StringIO ( ) ):
        LBFGSMinimize_SystemGeometry ( system, log = None, maximumIterations = steps, rmsGradientTolerance = TOL_KJ )
    return np.array ( [ [ system.coordinates3[i, k] for k in range ( 3 ) ] for i in range ( len ( system.atoms ) ) ] )


def refine ( build_system, reference, metals = None, iterations = 3, base_edits = None ):
    """ Per-instance metal terms reproducing the reference metal site. Returns ( edits, report ). """
    ref = np.asarray ( reference, dtype = float )
    edits = { k: dict ( v ) for k, v in ( base_edits or { } ).items ( ) }
    edits.setdefault ( "bond_terms", { } ); edits.setdefault ( "angle_terms", { } )
    s = build_system ( edits )
    bonds, angles = metal_terms ( s, metals )
    labels, atom_types = _model_types ( s )
    if atom_types is None: raise ValueError ( "metal_refine needs the EasyHybrid DYFF model (dyff_edits.MMModelDYFFEdited)" )
    r0 = { b: float ( np.linalg.norm ( ref[b[0]] - ref[b[1]] ) ) for b in bonds }
    t0 = { a: _angle ( ref, *a ) for a in angles }
    report = { "iterations": [ ] }
    for it in range ( iterations + 1 ):
        for ( m, l ), r in r0.items ( ):
            key = "{}-{}".format ( min ( m, l ), max ( m, l ) )
            edits["bond_terms"][key] = [ r, _bond_k ( atom_types, labels[m], labels[l], r ) ]
        for ( i, m, k ), t in t0.items ( ):
            a, b = ( i, k ) if i < k else ( k, i )
            edits["angle_terms"]["{}-{}-{}".format ( a, m, b )] = _angle_terms ( atom_types, labels[i], labels[m], labels[k],
                                                                                  r0[( m, i )], r0[( m, k )], t )
        s = build_system ( edits )
        x = _optimise ( s )
        dr = { b: float ( np.linalg.norm ( ref[b[0]] - ref[b[1]] ) - np.linalg.norm ( x[b[0]] - x[b[1]] ) ) for b in bonds }
        dt = { a: _angle ( ref, *a ) - _angle ( x, *a ) for a in angles }
        report["iterations"].append ( { "bond_MAE": float ( np.mean ( np.abs ( list ( dr.values ( ) ) ) ) ) if dr else 0.0,
                                        "angle_MAE": float ( np.mean ( np.abs ( list ( dt.values ( ) ) ) ) ) if dt else 0.0 } )
        if it == iterations: break
        for b in bonds:  r0[b] += dr[b]
        for a in angles: t0[a] = min ( 179.9, max ( 1.0, t0[a] + dt[a] ) )
    report["final_coordinates"] = x
    return edits, report
