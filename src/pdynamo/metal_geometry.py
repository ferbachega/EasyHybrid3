#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  metal_geometry.py
#
#  Coordination geometry of transition-metal atoms from the COORDINATES, for DYFF typing
#  (2026-10-13; DYFF_validation, tmQM study).
#
#  pDynamo's DetermineAtomGeometry gives a transition metal its geometry from the number of
#  neighbours only (4 -> Tet, 5 -> TBp, 6 -> Oct), so a square-planar Pt(II)/Pd(II)/Ni(II) is
#  typed as tetrahedral and the dyff-1.0 SPl types never match; it keeps a geometry that is
#  already set on the atom. MMModelDYFFEdited.BuildModel (dyff_edits) calls perceive() with the
#  system's coordinates and presets the metal geometries before typing.
#
#  Rules (validated on 108 541 tmQM complexes, DYFF_validation/scripts/tmqm_geometry.py):
#    * coordination SITES: the metal's bonded neighbours; neighbours bonded to each other (an
#      eta-n ligand: Cp, arene, alkene, allyl) form ONE site at their centroid;
#    * 2 sites: Lin (angle >= 160 deg) / bent; 3 sites: Tri (planar) / T-shaped / pyramidal;
#    * 4 sites: tau4 = (360 - a - b)/141 (Yang et al. 2007): SPl <= 0.35, Tet >= 0.55 (cut-offs
#      from the valley of the tmQM tau4 histogram), see-saw in between (typed Tet);
#    * 5 sites: tau5 = (b - a)/60 (Addison et al. 1984): TBp >= 0.5, else SPy;
#    * 6 sites: Oct (3 trans angles >= 150 deg), trigonal prism (none) and distorted Oct -> Oct;
#    * 7+ sites: PBp (7), Cub (8), TTP (9), else undefined.
#  ideal_angles(): for geometries the UFF angle forms cannot describe with one natural angle
#  (TBp: 90/120/180; PBp: 72/90/144/180), every L-M-L angle gets the ideal angle nearest to the
#  observed one, used by the per-instance angle factory of dyff_edits.

import math
import numpy as np

TRANSITION_METALS = { 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
                      57, 72, 73, 74, 75, 76, 77, 78, 79, 80 }

TAU4_SPL, TAU4_TET = 0.35, 0.55
IDEAL = { "TBp": ( 90.0, 120.0, 180.0 ), "PBp": ( 72.0, 90.0, 144.0, 180.0 ) }


def _angle ( u, v ):
    c = float ( np.dot ( u, v ) / ( np.linalg.norm ( u ) * np.linalg.norm ( v ) ) )
    return math.degrees ( math.acos ( max ( -1.0, min ( 1.0, c ) ) ) )


def classify ( vectors ):
    """ ( geometry tag, details ) from the metal -> site vectors. """
    n = len ( vectors )
    A = sorted ( ( _angle ( vectors[i], vectors[j] ) for i in range ( n ) for j in range ( i + 1, n ) ), reverse = True )
    if n == 0: return "Iso", { }
    if n == 1: return "Ter", { }
    if n == 2: return ( "Lin" if A[0] >= 160.0 else "Ben" ), { }
    if n == 3:
        if A[0] >= 160.0: return "TSh", { }
        return ( "Tri" if sum ( A ) >= 350.0 else "Pyr" ), { }
    if n == 4:
        t4 = ( 360.0 - A[0] - A[1] ) / 141.0
        return ( "SPl" if t4 <= TAU4_SPL else ( "Tet" if t4 >= TAU4_TET else "SSw" ) ), { "tau4": t4 }
    if n == 5:
        t5 = ( A[0] - A[1] ) / 60.0
        return ( "TBp" if t5 >= 0.5 else "SPy" ), { "tau5": t5 }
    if n == 6:
        trans = sum ( 1 for a in A if a >= 150.0 )
        return ( "Oct" if trans >= 3 else ( "TPr" if trans == 0 else "Oct-d" ) ), { "trans": trans }
    return { 7: "PBp", 8: "Cub", 9: "TTP" }.get ( n, "C{}".format ( n ) ), { }


# . coordinate class -> DYFF / pDynamo geometry tag of the metal type
DYFF_TAG = { "Lin": "Lin", "Ben": "Lin", "Tri": "Tri", "TSh": "Tri", "Pyr": "Tri", "Tet": "Tet", "SSw": "Tet", "SPl": "SPl",
             "TBp": "TBp", "SPy": "SPy", "Oct": "Oct", "Oct-d": "Oct", "TPr": "Oct", "PBp": "PBp", "Cub": "Cub", "TTP": "TTP" }


def sites ( metal, neighbours, coordinates, bonded ):
    """ Coordination sites of a metal: lists of neighbour indices (one per site; eta groups merged).
        bonded ( i, j ) -> True when two neighbours are bonded to each other. """
    groups, seen = [ ], set ( )
    for i in neighbours:
        if i in seen: continue
        stack, group = [ i ], [ ]
        seen.add ( i )
        while stack:
            k = stack.pop ( ); group.append ( k )
            for j in neighbours:
                if j not in seen and bonded ( k, j ):
                    seen.add ( j ); stack.append ( j )
        groups.append ( sorted ( group ) )
    return groups


# . eta-n ligands (2026-10-13): coordination sites occupied by a hapto ligand (eta2 1, eta3/eta4 2, eta5-eta8 3)
#   -> "effective" coordination number of the metal, used for the angles between its classic ligands
#   (a piano-stool CpM(CO)3 is pseudo-octahedral: OC-M-CO ~ 90 deg).
HAPTO_SITES = { 2: 1, 3: 2, 4: 2, 5: 3, 6: 3, 7: 3, 8: 3 }
IDEAL_EFFECTIVE = { 6: ( 90.0, 180.0 ), 5: ( 90.0, 120.0, 180.0 ), 4: ( 109.47, ), 7: ( 72.0, 90.0, 144.0, 180.0 ) }


def coordinated_waters ( atomic_numbers, coordinates, bonds, metal, radius = None ):
    """ Water oxygens (O bonded to exactly two H and nothing else) in contact with the metal but NOT bonded to
        it: EasyHybrid's protein preparation (dyff_prep) keeps coordinated waters without a metal-O bond, bound by
        electrostatics only. Contact: d <= r_M + 0.66 + 0.5 A (dyff_prep.sites_from_system, same rule). """
    xyz = np.asarray ( coordinates, dtype = float )
    nb = { }
    for i, j in bonds:
        nb.setdefault ( i, set ( ) ).add ( j ); nb.setdefault ( j, set ( ) ).add ( i )
    if radius is None:
        try:
            from pdynamo.dyff_prep import METALS
            from pScientific import PeriodicTable
            radius = METALS.get ( PeriodicTable.Symbol ( atomic_numbers[metal] ), ( 0, 0, 1.4 ) )[2]
        except Exception:
            radius = 1.4
    out = [ ]
    for o, z in enumerate ( atomic_numbers ):
        if z != 8 or o in nb.get ( metal, ( ) ): continue
        partners = nb.get ( o, set ( ) )
        if len ( partners ) != 2 or any ( atomic_numbers[k] != 1 for k in partners ): continue
        if np.linalg.norm ( xyz[o] - xyz[metal] ) <= radius + 0.66 + 0.5: out.append ( o )
    return out


def perceive ( atomic_numbers, coordinates, bonds, include_waters = True ):
    """ { metal index: { "class", "tag", "sites", "hapto", "details", "angles" } } for every transition
        metal with at least one bond. angles: { ( i, metal, k ): ideal angle (deg) } for the metal-centred
        angles whose ideal value must be set per instance (TBp, PBp; None elsewhere). """
    xyz = np.asarray ( coordinates, dtype = float )
    nb = { }
    for i, j in bonds:
        nb.setdefault ( i, set ( ) ).add ( j ); nb.setdefault ( j, set ( ) ).add ( i )
    out = { }
    for m, z in enumerate ( atomic_numbers ):
        if z not in TRANSITION_METALS or not nb.get ( m ): continue
        neighbours = sorted ( nb[m] )
        groups = sites ( m, neighbours, xyz, lambda a, b: b in nb.get ( a, ( ) ) and a != m and b != m )
        # . coordinated waters kept without a bond (protein preparation): sites for the CLASSIFICATION only
        #   (no bond / angle term involves them); 2026-10-13, full-protein test (a Co with 3 waters looked linear)
        waters = coordinated_waters ( atomic_numbers, xyz, bonds, m ) if include_waters else [ ]
        groups = groups + [ [ o ] for o in waters ]
        vectors = [ xyz[g].mean ( axis = 0 ) - xyz[m] for g in groups ]
        cls, details = classify ( vectors )
        info = { "class": cls, "tag": DYFF_TAG.get ( cls ), "sites": groups, "hapto": [ len ( g ) for g in groups if len ( g ) > 1 ],
                 "details": details, "angles": None, "unbonded_waters": waters }
        if info["hapto"]:
            # . eta complex: hapto atoms (their bonds and angles get special treatment in dyff_edits) and the
            #   ideal angles between CLASSIC ligands from the effective coordination number
            hapto_atoms = sorted ( i for g in groups if len ( g ) > 1 for i in g )
            classic = [ g[0] for g in groups if len ( g ) == 1 ]
            eff = sum ( HAPTO_SITES.get ( len ( g ), 3 ) if len ( g ) > 1 else 1 for g in groups )
            info["hapto_atoms"] = hapto_atoms; info["effective_cn"] = eff
            ideal = IDEAL_EFFECTIVE.get ( eff )
            angles = { }
            for a in range ( len ( classic ) ):
                for b in range ( a + 1, len ( classic ) ):
                    i, k = classic[a], classic[b]
                    if ideal is None: continue
                    observed = _angle ( xyz[i] - xyz[m], xyz[k] - xyz[m] )
                    angles[( i, m, k )] = min ( ideal, key = lambda t: abs ( t - observed ) )
            for i in hapto_atoms:
                for k in neighbours:
                    if k != i: angles[( min ( i, k ), m, max ( i, k ) )] = None          # . None: no angle term (see dyff_edits)
            info["angles"] = angles
        elif cls in IDEAL:
            ideal = IDEAL[cls]; angles = { }
            for a in range ( len ( neighbours ) ):
                for b in range ( a + 1, len ( neighbours ) ):
                    i, k = neighbours[a], neighbours[b]
                    observed = _angle ( xyz[i] - xyz[m], xyz[k] - xyz[m] )
                    angles[( i, m, k )] = min ( ideal, key = lambda t: abs ( t - observed ) )
            info["angles"] = angles
        out[m] = info
    return out


def perceive_system ( system ):
    """ perceive() on a pDynamo System (bonds from its connectivity, coordinates3). """
    atomic_numbers = [ a.atomicNumber for a in system.atoms ]
    coordinates = [ [ system.coordinates3[i, 0], system.coordinates3[i, 1], system.coordinates3[i, 2] ] for i in range ( len ( system.atoms ) ) ]
    return perceive ( atomic_numbers, coordinates, list ( system.connectivity.bondIndices ) )
