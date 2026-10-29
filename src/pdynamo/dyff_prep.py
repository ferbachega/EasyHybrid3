#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: preparation of a system for DYFF + residue library (no OPLS typing)
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description (2026-10-08):
#      PDB (or any loaded system) -> a system with hydrogens, bond orders, residue
#      names, protonation variants and metal-ligand bonds, ready for "DYFF + residue
#      library". Reuses the analysis and the PDB-component build of the OPLS
#      preparation (pdynamo/opls/prep.py) but types nothing with OPLS, so it accepts
#      what OPLS refuses:
#        metals    single-atom metal residues (Zn, Cu, Fe, Ni, Co, Mn, Cd, Hg, Mg, Ca, ...)
#                  are kept. The donor atoms are found by distance (sum of covalent
#                  radii + 0.5 A). A site is either BONDED (metal-donor bonds in the
#                  connectivity; default for transition metals) or IONIC (no bonds;
#                  default for Mg and Ca). The protonation of the donor residues follows
#                  the site: a bound Cys is a thiolate, a His bound through NE2 is
#                  delta-protonated (and vice versa), Asp/Glu are deprotonated. The
#                  metal and its donor residues are left to the "DYFF + residue
#                  library" step as one group outside the library: CM5 charges and a
#                  Hessian (GFN2-xTB or g-xTB) on a capped subsystem around the site;
#                  these residues are "pending" until refined.
#        truncated residues (missing heavy atoms) can be completed: each missing
#                  atom is placed from its bonded neighbours with a standard bond length
#                  and angle, extended (anti) conformation -- a geometric guess the
#                  user is warned about (minimize / refine afterwards).
#        ligands   organic residues without a PDB component get hydrogens and bond
#                  orders from Open Babel at the given pH (opls/ligand.ligand_chemistry,
#                  chemistry only).
#      Metal-containing cofactors (heme, Fe-S clusters, ...) are reported as not yet
#      supported.
#

import math
import os
import tempfile
from collections import OrderedDict

import numpy as np

# . metals: element -> ( default formal charge, default multiplicity, covalent radius A, default site model )
METALS = OrderedDict ( [
    ( "Zn", ( 2, 1, 1.22, "bonded" ) ), ( "Cu", ( 2, 2, 1.32, "bonded" ) ), ( "Fe", ( 3, 6, 1.32, "bonded" ) ),
    ( "Ni", ( 2, 3, 1.24, "bonded" ) ), ( "Co", ( 2, 4, 1.26, "bonded" ) ), ( "Mn", ( 2, 6, 1.39, "bonded" ) ),
    ( "Cd", ( 2, 1, 1.44, "bonded" ) ), ( "Hg", ( 2, 1, 1.32, "bonded" ) ), ( "Mo", ( 4, 1, 1.54, "bonded" ) ),
    ( "W",  ( 4, 1, 1.62, "bonded" ) ), ( "V",  ( 3, 3, 1.53, "bonded" ) ), ( "Cr", ( 3, 4, 1.39, "bonded" ) ),
    ( "Pt", ( 2, 1, 1.36, "bonded" ) ), ( "Ag", ( 1, 1, 1.45, "bonded" ) ), ( "Au", ( 1, 1, 1.36, "bonded" ) ),
    ( "Mg", ( 2, 1, 1.41, "ionic" ) ),  ( "Ca", ( 2, 1, 1.76, "ionic" ) ) ] )
DONORS = { "N": 0.71, "O": 0.66, "S": 1.05, "Se": 1.20, "Cl": 1.02, "Br": 1.20, "C": 0.76 }
COORDINATION_SLACK = 0.5       # A over the sum of covalent radii
COVALENT_CONTACT = 2.0         # A: a non-library residue this close to another residue is covalently bound
_BOND_RADII = { 1: 0.31, 6: 0.76, 7: 0.71, 8: 0.66, 15: 1.07, 16: 1.05, 34: 1.20 }


def _element ( z ):
    from pScientific import PeriodicTable
    return PeriodicTable.Symbol ( z )


#=============================================================================
#  Analysis
#=============================================================================
def analyze ( system ):
    """ { "analysis" (opls.prep.PrepAnalysis), "metals" [ site ], "ligands" [ residue
        indices ], "metal_cofactors" [ residue indices ], "truncated" [ ( index, label,
        missing ) ] }. site = { "residue", "label", "element", "charge", "multiplicity",
        "model", "donors": [ ( residue index, atom name, distance ) ] }. """
    from pdynamo.opls import prep
    analysis = prep.analyze_system ( system )
    residues = analysis.residues
    metals, ligands, metal_cofactors = [ ], [ ], [ ]
    for i, r in enumerate ( residues ):
        elements = [ _element ( z ) for ( name, z, xyz ) in r.atoms ]
        heavy = [ a for a in r.atoms if a[1] != 1 ]
        if len ( heavy ) == 1 and elements and _element ( heavy[0][1] ) in METALS:
            element = _element ( heavy[0][1] )
            charge, mult, radius, model = METALS[element]
            metals.append ( { "residue": i, "label": r.label, "element": element, "charge": charge, "multiplicity": mult,
                              "model": model, "radius": radius, "xyz": np.asarray ( heavy[0][2], float ), "donors": [ ] } )
        elif r.kind == "unsupported":
            ( metal_cofactors if any ( e in METALS for e in elements ) else ligands ).append ( i )
    # . donors: N, O, S, ... of other residues within the sum of covalent radii + slack
    for site in metals:
        candidates = { }
        for j, r in enumerate ( residues ):
            if j == site["residue"]: continue
            for ( name, z, xyz ) in r.atoms:
                element = _element ( z )
                if element not in DONORS or xyz is None: continue
                d = float ( np.linalg.norm ( np.asarray ( xyz ) - site["xyz"] ) )
                if d <= site["radius"] + DONORS[element] + COORDINATION_SLACK:
                    site["donors"].append ( ( j, name, d ) )
                    candidates.setdefault ( ( j, name ), ( element, np.asarray ( xyz, float ), d ) )
        # . 2026-10-13: a carbon is not a donor when it is bonded (< 1.6 A) to a closer O/N/S donor of the same
        #   residue -- the carboxylate C of a bidentate Asp/Glu falls within the contact distance and became a
        #   metal-C bond (full-protein test, 2Y3B Co-Glu63:CD). Carbon donors of CO / CN ligands are kept.
        drop = set ( )
        for ( j, name ), ( element, xyz, d ) in candidates.items ( ):
            if element != "C": continue
            for ( j2, name2 ), ( element2, xyz2, d2 ) in candidates.items ( ):
                if j2 == j and element2 in ( "O", "N", "S" ) and d2 < d and np.linalg.norm ( xyz - xyz2 ) < 1.6:
                    drop.add ( ( j, name ) ); break
        site["donors"] = [ t for t in site["donors"] if ( t[0], t[1] ) not in drop ]
        site["donors"].sort ( key = lambda t: t[2] )
    truncated = [ ( i, r.label, list ( r.missing_heavy ) ) for i, r in enumerate ( residues ) if r.missing_heavy ]
    # . residues without a library component bonded to another residue (terminal caps ACE/NME/NH2 of
    #   peptides, covalent ligands): the preparation cannot link them -- they must be removed
    covalent = [ ]
    for i in ligands + metal_cofactors:
        mine = [ np.asarray ( xyz ) for ( n, z, xyz ) in residues[i].atoms if z != 1 and xyz is not None ]
        if not mine: continue
        for j, r in enumerate ( residues ):
            if j == i or r.kind == "water": continue
            if any ( z != 1 and xyz is not None and min ( np.linalg.norm ( np.asarray ( xyz ) - m ) for m in mine ) < COVALENT_CONTACT
                     for ( n, z, xyz ) in r.atoms ):
                covalent.append ( i ); break
    return { "analysis": analysis, "metals": metals, "ligands": ligands, "metal_cofactors": metal_cofactors, "truncated": truncated,
             "covalent": covalent, "removed": [ ] }


def remove_residues ( info, indices ):
    """ Drops residues from the build (e.g. the covalently bound ones): they are not
        written to the structure that is built. """
    indices = set ( indices )
    info["removed"] = sorted ( set ( info.get ( "removed", [ ] ) ) | indices )
    info["ligands"] = [ i for i in info["ligands"] if i not in indices ]
    info["metal_cofactors"] = [ i for i in info["metal_cofactors"] if i not in indices ]
    info["covalent"] = [ i for i in info["covalent"] if i not in indices ]
    info["truncated"] = [ t for t in info["truncated"] if t[0] not in indices ]
    for site in info["metals"]:
        site["donors"] = [ d for d in site["donors"] if d[0] not in indices ]


def apply_site_protonation ( info ):
    """ Protonation of the donor residues of BONDED sites (see the module docstring).
        Returns [ ( residue label, variant ) ] changed. """
    residues = info["analysis"].residues
    changed = [ ]
    for site in info["metals"]:
        if site["model"] != "bonded": continue
        for ( j, name, d ) in site["donors"]:
            r = residues[j]
            variant = None
            if   r.library_name == "CYS" and name == "SG":  variant = "Deprotonated"
            elif r.library_name == "HIS" and name == "NE2": variant = "Delta Protonated"
            elif r.library_name == "HIS" and name == "ND1": variant = "Epsilon Protonated"
            elif r.library_name in ( "ASP", "GLU" ):         variant = "Deprotonated"
            if variant is not None and r.variant != variant:
                r.variant, r.variant_source = variant, "bound to {}".format ( site["label"] )
                changed.append ( ( r.label, variant ) )
    return changed


#=============================================================================
#  Completion of missing heavy atoms (geometric guess)
#=============================================================================
def complete_heavy_atoms ( system ):
    """ Coordinates for undefined heavy atoms from their defined bonded neighbours:
        bond length = sum of covalent radii, angle 111 degrees (sp3) or 120 (atoms with a
        multiple bond), dihedral 180 (extended). Returns the indices completed. """
    from pScientific import PeriodicTable
    c3 = system.coordinates3
    undefined = set ( c3.undefined or [ ] )
    heavy = [ a.index for a in system.atoms if a.atomicNumber != 1 and a.index in undefined ]
    if not heavy: return [ ]
    nb = { a.index: [ ] for a in system.atoms }
    multiple = set ( )
    for b in system.connectivity.bonds:
        i, j = b.node1.index, b.node2.index
        nb[i].append ( j ); nb[j].append ( i )
        if getattr ( b.type, "name", "Single" ) != "Single" or getattr ( b, "isAromatic", False ): multiple.update ( ( i, j ) )
    x = lambda i: np.array ( [ c3[i, 0], c3[i, 1], c3[i, 2] ] )
    def radius ( i ): return _BOND_RADII.get ( system.atoms[i].atomicNumber, 0.77 )
    done = [ ]
    pending = list ( heavy )
    while pending:
        progress = False
        for a in list ( pending ):
            defined = [ n for n in nb[a] if n not in undefined ]
            if not defined: continue
            j = defined[0]
            others = [ n for n in nb[j] if n != a and n not in undefined ]
            length = radius ( a ) + radius ( j )
            if len ( others ) == 2 and j not in multiple:
                # . sp3 centre with two bonds: tetrahedral position, out of their plane
                u1, u2 = [ ( x ( o ) - x ( j ) ) / np.linalg.norm ( x ( o ) - x ( j ) ) for o in others ]
                m, n = u1 + u2, np.cross ( u1, u2 )
                m = m / np.linalg.norm ( m ) if np.linalg.norm ( m ) > 1e-6 else np.array ( [ 1.0, 0.0, 0.0 ] )
                n = n / np.linalg.norm ( n ) if np.linalg.norm ( n ) > 1e-6 else np.cross ( m, [ 0.0, 0.0, 1.0 ] )
                alpha = math.radians ( 54.74 )
                p = x ( j ) + length * ( -math.cos ( alpha ) * m + math.sin ( alpha ) * n )
            elif len ( others ) >= 2:
                u = sum ( ( x ( o ) - x ( j ) ) / np.linalg.norm ( x ( o ) - x ( j ) ) for o in others )
                direction = -u / np.linalg.norm ( u ) if np.linalg.norm ( u ) > 1e-6 else np.array ( [ 1.0, 0.0, 0.0 ] )
                p = x ( j ) + length * direction
            elif len ( others ) == 1:
                k = others[0]
                angle = math.radians ( 120.0 if j in multiple else 111.0 )
                ls = [ n for n in nb[k] if n != j and n not in undefined ]
                kj = ( x ( j ) - x ( k ) ) / np.linalg.norm ( x ( j ) - x ( k ) )
                ref = ( x ( ls[0] ) - x ( k ) ) if ls else np.array ( [ 0.0, 0.0, 1.0 ] )
                perp = ref - np.dot ( ref, kj ) * kj
                if np.linalg.norm ( perp ) < 1e-6: perp = np.cross ( kj, [ 1.0, 0.0, 0.0 ] )
                perp /= np.linalg.norm ( perp )
                # . angle k-j-a = angle, dihedral l-k-j-a = 180 (a leaves opposite to l)
                p = x ( j ) + length * ( -math.cos ( angle ) * kj - math.sin ( angle ) * perp )
            else:
                p = x ( j ) + length * np.array ( [ 1.0, 0.0, 0.0 ] )
            for k3 in range ( 3 ): c3[a, k3] = float ( p[k3] )
            c3.FlagCoordinateAsDefined ( a )
            undefined.discard ( a ); pending.remove ( a ); done.append ( a ); progress = True
        if not progress: break
    return done


def place_completed_side_chains ( system, completed, clash = 3.0 ):
    """ 2026-10-15: rotamer choice for the heavy atoms placed by complete_heavy_atoms
        (extended conformation, blind to the surroundings: in 4oas the completed Gln18
        side chain landed 1.6 A from the Arg97 guanidinium and the minimization then
        distorted the Arg). For every group of completed atoms, the dihedrals about the
        bonds that move them (chi1, chi2, ...) are scanned -- 60/180/300 deg for sp3
        bonds, every 30 deg for a bond to a trigonal atom -- and the combination whose
        smallest distance from the moved atoms to every heavy atom further than three
        bonds away is the largest is kept. Returns [ ( residue path, closest contact
        before, after ) ]. """
    import itertools
    completed = set ( completed )
    if not completed: return [ ]
    c3 = system.coordinates3
    n_atoms = len ( system.atoms )
    x = np.array ( [ [ c3[i, 0], c3[i, 1], c3[i, 2] ] for i in range ( n_atoms ) ] )
    z = np.array ( [ a.atomicNumber for a in system.atoms ] )
    nb = { i: set ( ) for i in range ( n_atoms ) }
    for i, j in system.connectivity.bondIndices: nb[i].add ( j ); nb[j].add ( i )
    heavy = np.flatnonzero ( z != 1 )
    report = [ ]
    # . groups: connected completed atoms
    seen = set ( )
    for start in sorted ( completed ):
        if start in seen: continue
        group, stack = set ( ), [ start ]
        while stack:
            a = stack.pop ( )
            if a in group: continue
            group.add ( a ); stack.extend ( n for n in nb[a] if n in completed and n not in group )
        seen |= group
        # . root: the completed atom bonded to a defined one
        roots = [ ( a, d ) for a in group for d in nb[a] if d not in completed and z[d] != 1 ]
        if len ( roots ) != 1: continue                       # . a ring closure or several anchors: left as built
        root, anchor = roots[0]
        before_anchor = [ d for d in nb[anchor] if d not in completed and d != root and z[d] != 1 ]
        if not before_anchor: continue
        # . tree of the group from the root; pivots ( u, v ): rotation about u-v moves the subtree of v
        parent, order, stack = { root: anchor }, [ ], [ root ]
        while stack:
            a = stack.pop ( ); order.append ( a )
            for n in nb[a]:
                if n in group and n not in parent: parent[n] = a; stack.append ( n )
        subtree = { a: { a } for a in order }
        for a in reversed ( order ):
            if parent[a] in subtree and parent[a] != anchor: subtree[parent[a]] |= subtree[a]
        pivots = [ ( before_anchor[0], anchor, root, subtree[root] ) ]
        for a in order:
            kids = [ k for k in order if parent.get ( k ) == a ]
            if kids and len ( subtree[a] ) > 1:
                moved = set ( ).union ( *( subtree[k] for k in kids ) )
                pivots.append ( ( parent[a], a, kids[0], moved ) )
        def steps ( v ):
            # . hydrogens are not built yet: three heavy neighbours = a trigonal atom (amide
            #   C, Arg CZ) or a branch point -- a finer grid; otherwise the sp3 rotamers
            return [ math.radians ( d ) for d in ( range ( 0, 360, 30 ) if len ( nb[v] ) == 3 else ( 60, 180, 300 ) ) ]
        # . atoms to test against: heavy atoms more than three bonds from the group
        near = set ( group )
        for _ in range ( 3 ):
            near |= set ( ).union ( *( nb[a] for a in near ) )
        others = np.array ( [ h for h in heavy if h not in near ] )
        moving = np.array ( sorted ( a for a in group if z[a] != 1 ) )
        def score ( coords ):
            if len ( others ) == 0: return 99.0
            d = np.linalg.norm ( coords[moving][:, None, :] - coords[others][None, :, :], axis = 2 )
            return float ( d.min ( ) )
        def dihedral ( p0, p1, p2, p3 ):
            b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
            b1 = b1 / np.linalg.norm ( b1 )
            v = b0 - np.dot ( b0, b1 ) * b1; w = b2 - np.dot ( b2, b1 ) * b1
            return math.atan2 ( np.dot ( np.cross ( b1, v ), w ), np.dot ( v, w ) )
        def rotate ( coords, u, v, atoms, angle ):
            axis = coords[v] - coords[u]; axis /= np.linalg.norm ( axis )
            K = np.array ( [ [ 0, -axis[2], axis[1] ], [ axis[2], 0, -axis[0] ], [ -axis[1], axis[0], 0 ] ] )
            R = np.eye ( 3 ) + math.sin ( angle ) * K + ( 1.0 - math.cos ( angle ) ) * K @ K
            idx = sorted ( atoms )
            coords[idx] = ( coords[idx] - coords[v] ) @ R.T + coords[v]
        before = score ( x )
        best, best_score = x[sorted ( group )].copy ( ), before
        grids = [ steps ( v ) for ( u, v, w, moved ) in pivots ]
        for combo in itertools.product ( *grids ):
            y = x.copy ( )
            for ( u, v, w, moved ), target in zip ( pivots, combo ):
                ref = [ n for n in nb[u] if n != v ]
                if not ref: continue
                now = dihedral ( y[ref[0]], y[u], y[v], y[w] )
                rotate ( y, u, v, moved, target - now )
            sc = score ( y )
            if sc > best_score + 1.0e-6:
                best, best_score = y[sorted ( group )].copy ( ), sc
        for a, p in zip ( sorted ( group ), best ):
            for k in range ( 3 ): c3[a, k] = float ( p[k] )
            x[a] = p
        report.append ( ( system.atoms[root].parent.path if system.atoms[root].parent is not None else str ( root ),
                          round ( before, 2 ), round ( best_score, 2 ) ) )
    return report


_XH = { 6: 1.09, 7: 1.01, 8: 0.96, 16: 1.34 }

def repair_hydrogens ( system, close_parent = 0.9, close_any = 1.0 ):
    """ Rebuilds the hydrogens of every heavy atom that has one too close to it (< close_parent)
        or to any other atom (< close_any) -- pDynamo's builder puts a hydrogen almost on top of
        its parent when the parent's heavy neighbours are nearly coplanar (distorted crystal
        geometry, completed residues). Ideal directions from the heavy neighbours; a nearly
        planar centre uses the plane normal. Returns the number of hydrogens moved. """
    c3 = system.coordinates3
    n_atoms = len ( system.atoms )
    x = np.array ( [ [ c3[i, 0], c3[i, 1], c3[i, 2] ] for i in range ( n_atoms ) ] )
    z = [ a.atomicNumber for a in system.atoms ]
    nb = { i: [ ] for i in range ( n_atoms ) }
    for i, j in system.connectivity.bondIndices: nb[i].append ( j ); nb[j].append ( i )
    moved = 0
    for p in range ( n_atoms ):
        if z[p] == 1: continue
        hs = [ h for h in nb[p] if z[h] == 1 ]
        if not hs: continue
        def bad ( h ):
            if np.linalg.norm ( x[h] - x[p] ) < close_parent: return True
            d = np.linalg.norm ( x - x[h], axis = 1 ); d[h] = 9.0; d[p] = 9.0
            return d.min ( ) < close_any
        if not any ( bad ( h ) for h in hs ): continue
        heavy = [ ( x[q] - x[p] ) / np.linalg.norm ( x[q] - x[p] ) for q in nb[p] if z[q] != 1 ]
        length = _XH.get ( z[p], 1.0 )
        dirs = [ ]
        if len ( heavy ) >= 3 and len ( hs ) == 1:
            u = -sum ( heavy )
            if np.linalg.norm ( u ) < 0.3:
                n = np.cross ( heavy[1] - heavy[0], heavy[2] - heavy[0] ); u = n if np.dot ( n, u ) >= 0 else -n
            dirs = [ u ]
        elif len ( heavy ) == 2:
            m = heavy[0] + heavy[1]; n = np.cross ( heavy[0], heavy[1] )
            m = m / np.linalg.norm ( m ) if np.linalg.norm ( m ) > 1e-6 else np.cross ( heavy[0], [ 0.0, 0.0, 1.0 ] )
            n = n / np.linalg.norm ( n ) if np.linalg.norm ( n ) > 1e-6 else np.cross ( m, [ 0.0, 0.0, 1.0 ] )
            if len ( hs ) == 1: dirs = [ -m ]
            else:
                a = math.radians ( 54.74 ); dirs = [ -math.cos ( a ) * m + math.sin ( a ) * n, -math.cos ( a ) * m - math.sin ( a ) * n ]
        elif len ( heavy ) == 1:
            u = heavy[0]; ref = np.cross ( u, [ 1.0, 0.0, 0.0 ] )
            if np.linalg.norm ( ref ) < 1e-3: ref = np.cross ( u, [ 0.0, 1.0, 0.0 ] )
            ref /= np.linalg.norm ( ref ); other = np.cross ( u, ref )
            t = math.radians ( 109.47 if len ( hs ) != 2 else 120.0 )
            for k in range ( len ( hs ) ):
                phi = 2.0 * math.pi * k / max ( len ( hs ), 1 )
                dirs.append ( math.cos ( t ) * u + math.sin ( t ) * ( math.cos ( phi ) * ref + math.sin ( phi ) * other ) )
        elif not heavy and len ( hs ) == 2:
            # . a water (no heavy neighbour): the ideal orientation farthest from the other atoms
            rng = np.random.default_rng ( p )
            others = np.delete ( np.arange ( n_atoms ), [ p ] + hs )
            near = others[ np.linalg.norm ( x[others] - x[p], axis = 1 ) < 4.0 ]
            best, best_score = None, -1.0
            for _ in range ( 48 ):
                q, _r = np.linalg.qr ( rng.normal ( size = ( 3, 3 ) ) )
                half = math.radians ( 104.5 / 2.0 )
                cand = [ q @ np.array ( [ math.sin ( half ), math.cos ( half ), 0.0 ] ), q @ np.array ( [ -math.sin ( half ), math.cos ( half ), 0.0 ] ) ]
                pts = [ x[p] + length * c for c in cand ]
                score = min ( np.linalg.norm ( x[near] - pt, axis = 1 ).min ( ) if len ( near ) else 9.0 for pt in pts )
                if score > best_score: best, best_score = cand, score
            dirs = best
        for h, d in zip ( hs, dirs ):
            d = d / np.linalg.norm ( d )
            pos = x[p] + length * d
            for k in range ( 3 ): c3[h, k] = float ( pos[k] )
            x[h] = pos; moved += 1
    return moved


def planarize_conjugated_nh2 ( system, length = 1.01, built = None ):
    """ 2026-10-15: places the two hydrogens of every NH2 bonded to a trigonal carbon
        (amides of Asn/Gln, Arg NH1/NH2, amidines, anilines) in the plane of that
        carbon's substituents, at 120 deg. pDynamo's builder treats such N as a
        generic centre and turns the NH2 about the C-N bond arbitrarily (57 deg out of
        the amide plane on average in 4oas), which the minimization then only partly
        undoes. Each hydrogen goes to the planar position nearest to it (names kept).
        Only hydrogens in `built` (indices placed by the builder) are moved when it
        is given -- hydrogens read from the input file are kept as they are.
        Returns the number of NH2 groups placed. """
    c3 = system.coordinates3
    n_atoms = len ( system.atoms )
    x = np.array ( [ [ c3[i, 0], c3[i, 1], c3[i, 2] ] for i in range ( n_atoms ) ] )
    z = [ a.atomicNumber for a in system.atoms ]
    nb = { i: [ ] for i in range ( n_atoms ) }
    for i, j in system.connectivity.bondIndices: nb[i].append ( j ); nb[j].append ( i )
    placed = 0
    for n in range ( n_atoms ):
        if z[n] != 7 or len ( nb[n] ) != 3: continue
        hs = [ h for h in nb[n] if z[h] == 1 ]
        heavy = [ q for q in nb[n] if z[q] != 1 ]
        if len ( hs ) != 2 or len ( heavy ) != 1: continue
        if built is not None and not all ( h in built for h in hs ): continue
        c = heavy[0]
        others = [ q for q in nb[c] if q != n ]
        if z[c] != 6 or len ( nb[c] ) != 3 or len ( others ) != 2: continue
        normal = np.cross ( x[others[0]] - x[c], x[others[1]] - x[c] )
        if np.linalg.norm ( normal ) < 1.0e-6: continue
        normal /= np.linalg.norm ( normal )
        u = x[n] - x[c]; u -= np.dot ( u, normal ) * normal; u /= np.linalg.norm ( u )
        w = np.cross ( normal, u )
        sites = [ x[n] + length * ( 0.5 * u + math.sqrt ( 3.0 ) / 2.0 * w ), x[n] + length * ( 0.5 * u - math.sqrt ( 3.0 ) / 2.0 * w ) ]
        d = lambda h, p: float ( np.linalg.norm ( x[h] - p ) )
        if d ( hs[0], sites[0] ) + d ( hs[1], sites[1] ) > d ( hs[0], sites[1] ) + d ( hs[1], sites[0] ):
            sites = sites[::-1]
        for h, pos in zip ( hs, sites ):
            for k in range ( 3 ): c3[h, k] = float ( pos[k] )
            x[h] = pos
        placed += 1
    return placed


#=============================================================================
#  Build
#=============================================================================
def build ( info, label = None, complete_truncated = True, ph = 7.4, random_seed = 1 ):
    """ The prepared system and a report. info: from analyze(), with the sites'
        "model" / "charge" / "multiplicity" possibly edited (and apply_site_protonation()
        already called). Raises ValueError for metal-containing cofactors and, when
        complete_truncated is False, for truncated residues. """
    from pBabel import PDBFileReader, PDBModelLink, PDBModelVariant
    from pMolecule import BondType, System
    from pScientific.RandomNumbers import RandomNumberGenerator
    from pSimulation import BuildHydrogenCoordinates3FromConnectivity
    from pdynamo.opls import prep
    from pdynamo.opls import ligand as L
    analysis = info["analysis"]
    residues = analysis.residues
    if info["covalent"]:
        raise ValueError ( "Covalently bound residues without a library component (terminal caps, covalent ligands): {}. "
                           "Remove them (remove_residues).".format ( ", ".join ( residues[i].label for i in info["covalent"] ) ) )
    if info["metal_cofactors"]:
        names = sorted ( { residues[i].name for i in info["metal_cofactors"] } )
        raise ValueError ( "Metal-containing cofactors are not supported yet: {}.".format ( ", ".join ( names ) ) )
    if info["truncated"] and not complete_truncated:
        raise ValueError ( "Residues with missing heavy atoms: {}.".format ( ", ".join ( l for _, l, _ in info["truncated"] ) ) )
    work = tempfile.mkdtemp ( prefix = "dyff_prep_" )
    library_root = os.path.join ( work, "pdb_components" )
    report = { "completed": [ ], "ligands": [ ], "metal_sites": [ ], "warnings": [ ] }
    # . metals: one-atom PDB components with the site's formal charge
    for site in info["metals"]:
        r = residues[site["residue"]]
        name = r.name.upper ( )[:3]
        lig = L.LigandParameters ( name )
        lig.atoms = [ ( a[0], a[1], np.asarray ( a[2], float ) ) for a in r.atoms if a[1] != 1 ]
        lig.formal_charges = [ int ( site["charge"] ) ]
        lig.write_pdb_component ( library_root )
        r.library_name, r.kind, r.missing_heavy = name, "metal", [ ]
    # . organic ligands: hydrogens and bond orders (chemistry only)
    done_names = { }
    for i in info["ligands"]:
        r = residues[i]
        if r.name not in done_names:
            lig = L.ligand_chemistry ( r.name, r.atoms, ph = ph, bonds = r.bonds, formal_charges = r.formal_charges )
            lig.write_pdb_component ( library_root )
            done_names[r.name] = lig
            report["ligands"].append ( { "name": r.name, "atoms": len ( lig.atoms ), "charge": lig.total_charge } )
        else:
            lig = L.ligand_chemistry ( r.name, r.atoms, ph = ph, bonds = r.bonds, formal_charges = r.formal_charges )
        r.atoms = [ ( a[0], a[1], np.asarray ( a[2], float ) ) for a in lig.atoms ]
        r.library_name, r.kind, r.missing_heavy = r.name.upper ( )[:3], "ligand", [ ]
    analysis.unsupported = [ ]
    for i, label_, missing in info["truncated"]:
        residues[i].missing_heavy = [ ]
    if library_root not in analysis.extra_library_paths:
        analysis.extra_library_paths.insert ( 0, library_root )
    # . build from the PDB components (as prep.build_prepared_system), plus the metal bonds
    kept = [ i for i in range ( len ( residues ) ) if i not in set ( info.get ( "removed", [ ] ) ) ]
    path = prep.write_clean_pdb ( analysis, os.path.join ( work, "prepared.pdb" ), residue_indices = kept )
    if info.get ( "removed" ):
        report["warnings"].append ( "Removed: {}.".format ( ", ".join ( residues[i].label for i in info["removed"] ) ) )
    raw   = PDBFileReader.PathToPDBModel ( path, log = None, useSegmentEntityLabels = True )
    model = PDBFileReader.PathToPDBModel ( path, log = None, useSegmentEntityLabels = True )
    def component_of ( i ):
        segment, component_label = analysis.pdb_labels[i]
        entity = model.childIndex.get ( segment )
        return entity, ( entity.childIndex.get ( component_label ) if entity is not None else None )
    variants = { }
    for i, r in enumerate ( residues ):
        if r.variant is None or i not in analysis.pdb_labels: continue
        entity, component = component_of ( i )
        if component is None: continue
        entity.AddVariant ( PDBModelVariant.WithOptions ( component = component, label = r.variant ) )
        variants[r.label] = r.variant
    for ( i, j, d ) in analysis.disulfides:
        if i not in analysis.pdb_labels or j not in analysis.pdb_labels: continue
        _, left = component_of ( i ); _, right = component_of ( j )
        if left is not None and right is not None:
            model.AddLink ( PDBModelLink.WithOptions ( label = "Disulfide_Bridge", leftComponent = left, rightComponent = right ) )
    model.MakeAtomicModelFromComponentLibrary ( libraryPaths = list ( analysis.extra_library_paths ) + prep.library_paths ( ), log = None )
    model.ExtractAtomData ( raw, log = None )
    ( sequence, atom_mapping ) = model.MakeSequence ( )
    bonds = model.MakeBonds ( atom_mapping ) or [ ]
    atoms = sequence.GatherAtoms ( )
    position = { id ( a ): k for k, a in enumerate ( atoms ) }
    by_key = { ( a.parent.parent.label, a.parent.label, a.label ): k for k, a in enumerate ( atoms ) }
    tuples = [ ( position[id ( b.node1 )], position[id ( b.node2 )], b.type ) for b in bonds ]
    for site in info["metals"]:
        n_bondable = sum ( 1 for ( j, n, d ) in site["donors"] if residues[j].kind != "water" )
        if site["model"] == "bonded" and n_bondable < 2:
            site["model"] = "ionic"
            report["warnings"].append ( "{}: only {} protein donor(s) in the file (a site completed by crystal symmetry?): treated as "
                                        "ionic; build the biological assembly for a bonded site.".format ( site["label"], n_bondable ) )
        seg, comp = analysis.pdb_labels[site["residue"]]
        metal_name = residues[site["residue"]].atoms[0][0] if residues[site["residue"]].atoms else None
        m = by_key.get ( ( seg, comp, metal_name ) )
        entry = { "metal": site["label"], "element": site["element"], "model": site["model"], "charge": site["charge"],
                  "multiplicity": site["multiplicity"], "donors": [ ], "metal_index": m }
        for ( j, name, d ) in site["donors"]:
            dseg, dcomp = analysis.pdb_labels[j]
            k = by_key.get ( ( dseg, dcomp, name ) )
            water = residues[j].kind == "water"
            entry["donors"].append ( { "residue": residues[j].label, "atom": name, "distance": round ( d, 3 ), "index": k,
                                       "bonded": site["model"] == "bonded" and not water } )
            # . a water bound to a metal is not bonded to it: DYFF has no type for a three-coordinate
            #   water oxygen (a hydroxide would be: the user can change the water to OH- in the Builder)
            if water and site["model"] == "bonded":
                report["warnings"].append ( "{}: the coordinated water {} is kept without a metal-oxygen bond (DYFF has no "
                                            "three-coordinate water type); its interaction is electrostatic.".format ( site["label"], residues[j].label ) )
            if site["model"] == "bonded" and not water and m is not None and k is not None:
                tuples.append ( ( m, k, BondType.Single ) )
        report["metal_sites"].append ( entry )
    system = System.FromSequence ( sequence, bonds = tuples or None )
    system.coordinates3 = model.MakeCoordinates3 ( None )
    system.label = label or analysis.label or "system"
    if info["truncated"]:
        completed = complete_heavy_atoms ( system )
        report["completed"] = [ system.atoms[i].path for i in completed ]
        report["completed_rotamers"] = place_completed_side_chains ( system, completed )
        if completed:
            report["warnings"].append ( "{} heavy atoms of {} truncated residue(s) were placed geometrically (side-chain "
                                        "rotamer chosen away from the other atoms): minimize and check them.".format ( len ( completed ), len ( info["truncated"] ) ) )
    rng = RandomNumberGenerator.WithSeed ( random_seed ) if random_seed is not None else RandomNumberGenerator.WithRandomSeed ( )
    built = set ( system.coordinates3.undefined or [ ] )            # . hydrogens the builder will place
    BuildHydrogenCoordinates3FromConnectivity ( system, log = None, randomNumberGenerator = rng )
    report["repaired_hydrogens"] = repair_hydrogens ( system )
    report["planar_nh2"] = planarize_conjugated_nh2 ( system, built = built )
    undefined = sorted ( system.coordinates3.undefined or [ ] )
    report["type_overrides"] = type_overrides ( system, report )
    report.update ( atoms = len ( system.atoms ), variants = variants, disulfides = len ( analysis.disulfides ),
                    chain_breaks = len ( analysis.chain_breaks ), undefined = [ system.atoms[i].path for i in undefined ],
                    formal_charge = int ( sum ( getattr ( a, "formalCharge", 0 ) or 0 for a in system.atoms ) ) )
    if undefined:
        report["warnings"].append ( "{} atoms still without coordinates.".format ( len ( undefined ) ) )
    return system, report


_DYFF_LABELS = None

def _dyff_metal_labels ( ):
    """ { element: [ DYFF atom type labels ] } from the dyff-1.0 parameter set. """
    global _DYFF_LABELS
    if _DYFF_LABELS is None:
        import re
        from pdynamo.dyff_edits import dyff_parameter_dir   # EasyHybrid extended set (dyff-1.0-eh)
        path = os.path.join ( dyff_parameter_dir ( "dyff-1.0" ), "atomTypes.yaml" )
        labels = { }
        try:                                                   # rows: - [ "Cu:SPl:+2" , 29 , ... ] (a multi-document YAML)
            for line in open ( path ):
                m = re.match ( r'\s*-\s*\[\s*"([A-Za-z]+:[^"]+)"', line )
                if m: labels.setdefault ( m.group ( 1 ).split ( ":" )[0], [ ] ).append ( m.group ( 1 ) )
        except OSError:
            pass
        _DYFF_LABELS = labels
    return _DYFF_LABELS


_GEOMETRY = { 0: [ "Oct", "Tet", "SPl" ], 1: [ "Oct", "Tet" ], 2: [ "Lin", "Tet" ], 3: [ "Tri", "Tet" ], 4: [ "Tet", "SPl" ],
              5: [ "Oct", "Tet" ], 6: [ "Oct" ] }

def metal_type_label ( element, charge, n_bonds ):
    """ The DYFF type of a metal with this charge and number of bonds (geometry from the
        count; a free ion as an octahedral hydrated ion), or None. """
    labels = _dyff_metal_labels ( ).get ( element, [ ] )
    for geometry in _GEOMETRY.get ( n_bonds, [ "Oct" ] ):
        for want_charge in ( charge, ):
            for label in labels:
                parts = label.split ( ":" )
                if len ( parts ) >= 3 and parts[1] == geometry and parts[2] == "{:+d}".format ( int ( want_charge ) ): return label
    for label in labels:                                   # any type of that charge
        parts = label.split ( ":" )
        if len ( parts ) >= 3 and parts[2] == "{:+d}".format ( int ( charge ) ): return label
    return labels[0] if labels else None


def type_overrides ( system, report ):
    """ { metal atom index: DYFF label } for every metal site (the DYFF pattern typer has no
        free-ion type for most metals and none for a metal with a single bond): see
        session._apply_manual_atom_type_overrides. """
    nb = { }
    for i, j in system.connectivity.bondIndices: nb.setdefault ( i, [ ] ).append ( j ); nb.setdefault ( j, [ ] ).append ( i )
    out = { }
    for site in report.get ( "metal_sites", [ ] ):
        m = site.get ( "metal_index" )
        if m is None: continue
        # . geometry from the coordination number (donors in the file, a coordinated water
        #   included even though it is not bonded), not from the bonds of the model
        coordination = len ( site.get ( "donors", [ ] ) ) if site.get ( "model" ) == "bonded" else 0
        label = metal_type_label ( site["element"], site["charge"], max ( coordination, len ( nb.get ( m, [ ] ) ) ) )
        if label: out[m] = label
        # . donors left untyped by the extra metal bond (a carboxylate O- of Asp/Glu bonded to
        #   the metal, a hydroxide): fallback types, used only for untyped atoms (2026-10-08:
        #   every bonded site with an Asp/Glu donor failed DYFF typing)
        for k in nb.get ( m, [ ] ):
            fallback = _DONOR_FALLBACK.get ( _element ( system.atoms[k].atomicNumber ) )
            if fallback and k not in out:
                out[k] = fallback ( len ( nb.get ( k, [ ] ) ) ) if callable ( fallback ) else fallback
    return out


def site_groups ( system, report ):
    """ { metal atom index: { "charge", "multiplicity", "model" } } of the bonded sites
        -- for residue_library.rest_parameters (charge/spin of the site's group). """
    out = { }
    for s in report.get ( "metal_sites", [ ] ):
        if s.get ( "metal_index" ) is not None:
            out[s["metal_index"]] = { "charge": s["charge"], "multiplicity": s["multiplicity"], "model": s["model"] }
    return out


# . DYFF type of a metal donor left untyped, by element (bent/pyramidal O and S, N by its
#   number of bonds including the metal one)
_DONOR_FALLBACK = { "O": "O:Tet", "S": "S:Tet:+2", "Se": "Se:Tet:+2",
                    "N": lambda n: "N:Tet" if n >= 4 else "N:Res" }


def sites_from_system ( system, multiplicities = None ):
    """ ( type_overrides, site_states ) of the metal atoms of an already built system,
        read from the structure itself rather than from a preparation report -- so a
        system edited afterwards (Builder: a coordinated water made a bonded hydroxide,
        a donor added or removed) can get DYFF + residue library again.
          charge       : the atom's formal charge (the METALS default when it has none);
          multiplicity : multiplicities[atom path] (e.g. kept from the preparation), else
                         the METALS default;
          model        : "bonded" when the metal has bonds, else "ionic";
          DYFF type    : from the charge and the coordination number = bonded neighbours
                         + unbonded water oxygens in contact (as in type_overrides). """
    multiplicities = multiplicities or { }
    nb = { }
    for i, j in system.connectivity.bondIndices: nb.setdefault ( i, [ ] ).append ( j ); nb.setdefault ( j, [ ] ).append ( i )
    xyz = None
    overrides, states = { }, { }
    for atom in system.atoms:
        element = _element ( atom.atomicNumber )
        if element not in METALS: continue
        charge, mult, radius, _ = METALS[element]
        charge = int ( getattr ( atom, "formalCharge", 0 ) or 0 ) or charge
        mult   = int ( multiplicities.get ( atom.path, mult ) )
        bonded = nb.get ( atom.index, [ ] )
        if xyz is None:
            xyz = np.array ( [ [ system.coordinates3[k, c] for c in range ( 3 ) ] for k in range ( len ( system.atoms ) ) ] )
        waters = 0
        for other in system.atoms:
            if other.atomicNumber != 8 or other.index in bonded: continue
            heavy = [ k for k in nb.get ( other.index, [ ] ) if system.atoms[k].atomicNumber != 1 ]
            if heavy or len ( nb.get ( other.index, [ ] ) ) != 2: continue                  # only a water oxygen
            if np.linalg.norm ( xyz[other.index] - xyz[atom.index] ) <= radius + DONORS["O"] + COORDINATION_SLACK: waters += 1
        label = metal_type_label ( element, charge, len ( bonded ) + waters )
        if label: overrides[atom.index] = label
        # . donors the pattern typer cannot type with the extra metal bond (e.g. a
        #   hydroxide O- bonded to the metal and to its H): overrides are used only for
        #   atoms left untyped, so the typer's own type wins whenever it has one
        for k in bonded:
            donor = system.atoms[k]
            n_bonds = len ( nb.get ( k, [ ] ) )
            fallback = _DONOR_FALLBACK.get ( _element ( donor.atomicNumber ) )
            if fallback and k not in overrides:
                overrides[k] = fallback ( n_bonds ) if callable ( fallback ) else fallback
        states[atom.index] = { "charge": charge, "multiplicity": mult, "model": "bonded" if bonded else "ionic" }
    return overrides, states
