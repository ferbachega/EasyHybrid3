#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ligand_perception.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  [EN] 2026-10-14 user request: replace Open Babel ("obabel -p pH") in the
#  preparation of ligands that come without hydrogens (PDB HETATM residues),
#  so that system preparation no longer depends on an external program whose
#  version we do not control. No GUI and no RDKit: plain Python + numpy, so it
#  runs in the pDynamo worker threads and in scripts.
#
#  Pipeline (perceive_ligand):
#    1. connectivity from covalent radii (unless bonds are given);
#    2. hydrogen count of every heavy atom in its NEUTRAL form, from the local
#       geometry (bond angles, flat aromatic-length rings, terminal bond
#       lengths) -- the rules of the Builder's residue_chemistry.py
#       (2026-10-02), moved here unchanged except for the data structures;
#    3. Kekule bond orders by the project's Wang & Case perceiver (vismol),
#       the implicit hydrogens given to it as placeholder atoms;
#    4. aromatic bonds (Hueckel 4n+2, fused systems included);
#    5. protonation at the given pH by functional-group rules with typical
#       pKa values (PKA below). Acids lose a proton when pH > pKa, bases gain
#       one when pH < pKa. Every decision is reported;
#    6. formal charges from the final valences;
#    7. hydrogens placed by the Builder's closed-form VSEPR construction
#       (atom_ops._complete_vsepr_directions), planar on aromatic and
#       conjugated centres.
#
#  A molecule that already has hydrogens keeps them (bonds and orders are
#  perceived, no protonation change), like Open Babel without "-p".
# ============================================================================

import math
from collections import defaultdict

import numpy as np

# ---------------------------------------------------------------------------
#   element data
# ---------------------------------------------------------------------------
# Cordero et al. 2008 covalent radii (A), for connectivity only.
COVALENT_RADII = { "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "SI": 1.11, "P": 1.07, "S": 1.05,
                   "CL": 1.02, "AS": 1.19, "SE": 1.20, "BR": 1.20, "TE": 1.38, "I": 1.39 }
# X-H bond lengths (A) for the placed hydrogens.
XH_LENGTH = { "C": 1.09, "N": 1.01, "O": 0.96, "S": 1.34, "P": 1.42, "B": 1.19, "SI": 1.48, "SE": 1.46 }
DEFAULT_XH_LENGTH = 1.05

# Typical pKa values (aqueous) behind the protonation rules. A rule applies
# when pH > pKa (acid loses H) or pH < pKa (base gains H).
PKA = {
    "carboxylic acid":          4.2,
    "sulfonic/sulfuric acid":  -1.0,
    "phosphate, first OH":      1.5,
    "phosphate ester, second OH": 6.3,     # monoester / inorganic: dianion at pH 7.4
    "phosphonate, second OH":   7.6,       # stays a monoanion at pH 7.4
    "tetrazole NH":             4.9,
    "acylsulfonamide NH":       4.5,
    "sulfonylurea NH":          5.2,
    "amidine":                 12.0,
    "guanidine":               13.5,
}
# Aliphatic amines: base pKa by substitution, lowered by substituents on the alpha
# and beta carbons (Perrin-type increments; checks: 2-methoxyethylamine 9.4,
# 2,2,2-trifluoroethylamine 5.5, benzylamine 9.3, N-methylmorpholine 7.4,
# N-benzylmorpholine 6.4, glycinamide-like alpha-amino carbonyls ~8).
AMINE_BASE_PKA = { 2: 10.6, 1: 11.0, 0: 9.8 }          # by number of H on N: primary, secondary, tertiary
AMINE_PKA_SHIFT = {
    "aryl or vinyl on alpha carbon":   -1.3,
    "carbonyl on alpha carbon":        -2.5,
    "amide carbonyl on alpha carbon":  -2.0,     # lidocaine-type R2N-CH2-C(=O)NR 7.9
    "nitrile on alpha carbon":         -5.0,
    "O or N(sp2) on beta carbon":      -1.2,
    "F on beta carbon":                -1.7,
}
# A second basic amine 2-3 bonds from an already protonated one (ethylenediamine,
# piperazine) has its pKa lowered by about this much.
DIAMINE_PKA_SHIFT = 4.5

# ---------------------------------------------------------------------------
#   geometry helpers
# ---------------------------------------------------------------------------

def _angle ( a, b, c ):
    v1 = np.asarray ( a ) - np.asarray ( b ); v2 = np.asarray ( c ) - np.asarray ( b )
    n1 = np.linalg.norm ( v1 ); n2 = np.linalg.norm ( v2 )
    if n1 < 1e-6 or n2 < 1e-6:
        return 109.5
    return math.degrees ( math.acos ( max ( -1.0, min ( 1.0, float ( np.dot ( v1, v2 ) / ( n1 * n2 ) ) ) ) ) )


def _torsion ( a, b, c, d ):
    a, b, c, d = ( np.asarray ( x, float ) for x in ( a, b, c, d ) )
    b0, b1, b2 = b - a, c - b, d - c
    n1, n2 = np.cross ( b0, b1 ), np.cross ( b1, b2 )
    nb1 = np.linalg.norm ( b1 ) or 1.0
    m1 = np.cross ( n1, b1 / nb1 )
    return math.degrees ( math.atan2 ( float ( np.dot ( m1, n2 ) ), float ( np.dot ( n1, n2 ) ) ) )


def _dist ( p, q ):
    return float ( np.linalg.norm ( np.asarray ( p ) - np.asarray ( q ) ) )


def _pair ( a, b ):
    return tuple ( sorted ( ( a.upper ( ), b.upper ( ) ) ) )


def _key ( i, j ):
    return ( i, j ) if i < j else ( j, i )


def _sym ( s ):
    return s.strip ( ).upper ( )


def connectivity ( symbols, coords, factor = 1.20 ):
    """ Bonds (i, j) between atoms closer than factor x (sum of covalent radii);
        no H-H bonds, and every hydrogen bonded to its nearest partner only. """
    syms = [ _sym ( s ) for s in symbols ]
    xyz = np.asarray ( coords, float )
    n = len ( syms )
    bonds = set ( )
    for i in range ( n ):
        ri = COVALENT_RADII.get ( syms[i], 0.77 )
        for j in range ( i + 1, n ):
            if syms[i] == "H" and syms[j] == "H":
                continue
            d = float ( np.linalg.norm ( xyz[i] - xyz[j] ) )
            if d < 0.4:
                continue
            if d <= factor * ( ri + COVALENT_RADII.get ( syms[j], 0.77 ) ):
                bonds.add ( ( i, j ) )
    # . a hydrogen keeps only its nearest partner
    out = set ( )
    best = { }
    for ( i, j ) in bonds:
        for h, other in ( ( i, j ), ( j, i ) ):
            if syms[h] == "H":
                d = float ( np.linalg.norm ( xyz[h] - xyz[other] ) )
                if h not in best or d < best[h][0]:
                    best[h] = ( d, other )
    for ( i, j ) in bonds:
        if syms[i] == "H" and best[i][1] != j: continue
        if syms[j] == "H" and best[j][1] != i: continue
        out.add ( ( i, j ) )
    return sorted ( out )


# ---------------------------------------------------------------------------
#   step 2: neutral hydrogen counts from the geometry
#   (rules of gui/windows/builder/residue_chemistry.py, 2026-10-02)
# ---------------------------------------------------------------------------
_BASE_VALENCE = { "C": 4, "N": 3, "O": 2, "S": 2, "P": 3, "B": 3, "SI": 4, "SE": 2 }
_TERMINAL_LENGTH_RULES = {
    ( "O", "C" ): [ ( 1.30, 0 ) ],                     # C=O (carbonyl/carboxylate) else C-OH
    ( "O", "N" ): [ ( 1.30, 0 ) ],                     # nitro / N-oxide
    ( "O", "S" ): [ ( 1.60, 0 ) ],                     # sulfonyl/sulfoxide S=O
    ( "O", "P" ): [ ( 1.60, 0 ) ],                     # phosphoryl P=O
    ( "N", "C" ): [ ( 1.20, 0 ), ( 1.33, 1 ) ],        # nitrile, imine; else amine (2 H)
    ( "C", "C" ): [ ( 1.25, 1 ), ( 1.40, 2 ) ],        # alkyne, alkene; else methyl (3 H)
    ( "C", "N" ): [ ( 1.20, 1 ), ( 1.33, 2 ) ],
    ( "C", "O" ): [ ( 1.30, 2 ) ],                     # formyl-like C=O; else methoxy C (3 H)
    ( "S", "C" ): [ ( 1.72, 0 ) ],                     # thione C=S; else thiol
}
_MULTIPLE_BOND_MAX = { ( "C", "C" ): 1.43, ( "C", "N" ): 1.31, ( "C", "O" ): 1.30, ( "C", "S" ): 1.72,
                       ( "N", "N" ): 1.30, ( "N", "O" ): 1.30, ( "O", "S" ): 1.60, ( "O", "P" ): 1.60,
                       ( "N", "S" ): 1.60 }
# [EN] 2026-10-14: C-N, C-S, N-N, N-O and C-O limits widened (the Builder's values
# rejected a nitro-pyrazole C-N of 1.416 A and thiadiazole C-S of 1.82-1.85 A in
# PL-REX ligands). Saturated 5/6-rings of the same lengths are excluded first by
# the flatness test (puckered rings twist well beyond 15 degrees).
# Second widening (same day, PL-REX validation): fused aromatic C-C up to 1.444 A,
# ring C-N 1.438 A and oxadiazole C-O 1.434 A occur in real model geometries; the
# limits stay below the saturated single bonds (C-C 1.53, C-N 1.47, flat rings only).
# bond at most this long -> unmistakably a multiple bond, whatever the angles (a
# polyene C with a 112-degree crystal angle but a 1.36 A C=C, retinol in 1AQB)
_CLEAR_DOUBLE_MAX = { ( "C", "C" ): 1.38, ( "C", "N" ): 1.30, ( "C", "O" ): 1.26, ( "N", "N" ): 1.28, ( "C", "S" ): 1.66 }
_AROMATIC_BOND_MAX = { ( "C", "C" ): 1.47, ( "C", "N" ): 1.47, ( "N", "N" ): 1.44, ( "C", "O" ): 1.46,
                       ( "C", "S" ): 1.87, ( "N", "S" ): 1.72, ( "N", "O" ): 1.45, ( "C", "SE" ): 1.95 }


def small_rings ( adjacency, max_size = 7 ):
    """ Every simple cycle of up to max_size atoms (ordered atom lists), once each. """
    rings = { }
    for start in adjacency:
        stack = [ ( start, [ start ] ) ]
        while stack:
            node, path = stack.pop ( )
            if len ( path ) > max_size:
                continue
            for nb in adjacency[node]:
                if nb == start and len ( path ) >= 3:
                    key = frozenset ( path )
                    if key not in rings and len ( key ) == len ( path ):
                        rings[key] = list ( path )
                elif nb not in path and nb > start:
                    stack.append ( ( nb, path + [ nb ] ) )
    return list ( rings.values ( ) )


def _planar_rings ( syms, neighbors, coords, max_twist = 30.0 ):
    """ Flat 5- and 6-membered rings whose bonds all have aromatic lengths.
        max_twist 30 degrees (was 15 in the Builder): docked poses distort aromatic
        rings by 20-23 degrees; saturated rings are excluded by their C-C/C-N
        bonds (>= 1.50 A) whatever their twist. """
    out = [ ]
    adjacency = { a: set ( nbs ) for a, nbs in neighbors.items ( ) }
    for ring in small_rings ( adjacency, max_size = 6 ):
        n = len ( ring )
        if n not in ( 5, 6 ):
            continue
        twists = [ abs ( _torsion ( *[ coords[ring[( k + d ) % n]] for d in range ( 4 ) ] ) ) for k in range ( n ) ]
        if max ( twists ) > max_twist:
            continue
        angles = [ _angle ( coords[ring[k - 1]], coords[ring[k]], coords[ring[( k + 1 ) % n]] ) for k in range ( n ) ]
        # [EN] 2026-10-14: a regular hexagon (every interior angle 114-126) is
        # accepted with bonds 0.03 A longer (crystal ligands: uracil and pyrimidinone
        # rings with 1.45-1.50 A bonds in PLINDER 6WRE, 4X4L)
        slack = 0.03 if n == 6 and all ( 114.0 <= x <= 126.0 for x in angles ) else 0.0
        # . a flat hexagon (twist <= 10, every angle 116-124) cannot hold an sp3
        #   atom (~110 degrees) whatever its bond lengths (PLINDER 6ICQ: a xanthine
        #   refined with 1.46-1.52 A ring bonds)
        if n == 6 and max ( twists ) <= 10.0 and all ( 116.0 <= x <= 124.0 for x in angles ):
            slack = 0.10
        ok = True
        for k in range ( n ):
            a, b = ring[k], ring[( k + 1 ) % n]
            limit = _AROMATIC_BOND_MAX.get ( _pair ( syms[a], syms[b] ) )
            if limit is None or _dist ( coords[a], coords[b] ) > limit + slack:
                ok = False
                break
        if ok and n == 6:
            # . aromatic 6-rings have interior angles near 120, sp3 ring atoms near 110
            ok = sum ( angles ) / n >= 116.0 and min ( angles ) >= 108.0
        if ok:
            out.append ( ring )
    return out


def geometric_implicit_hydrogens ( syms, heavy_neighbors, coords ):
    """ {atom: hydrogen count} of every heavy atom in its neutral form, from the
        geometry. heavy_neighbors: {atom: [heavy neighbour atoms]}. Atoms of flat
        aromatic-length 5/6-rings are sp2 whatever their angle; a flat 5-ring with
        no lone-pair donor gets an N-H on the 2-connected N with the widest angle. """
    out = { }
    ring_sp2 = set ( )
    pyrrole_nh = set ( )
    for ring in _planar_rings ( syms, heavy_neighbors, coords ):
        ring_sp2.update ( ring )
        if len ( ring ) == 5:
            donors = [ a for a in ring if syms[a] in ( "O", "S", "SE" ) or ( syms[a] == "N" and len ( heavy_neighbors[a] ) == 3 ) ]
            if not donors:
                n2 = [ a for a in ring if syms[a] == "N" and len ( heavy_neighbors[a] ) == 2 ]
                if n2:
                    widest = max ( n2, key = lambda a: _angle ( coords[heavy_neighbors[a][0]], coords[a], coords[heavy_neighbors[a][1]] ) )
                    pyrrole_nh.add ( widest )
    for aid, nbs in heavy_neighbors.items ( ):
        sym = syms[aid]
        base = _BASE_VALENCE.get ( sym )
        if base is None:
            continue
        n = len ( nbs )
        if n == 0:
            out[aid] = base
            continue
        if n == 1:
            nb = nbs[0]
            d = _dist ( coords[aid], coords[nb] )
            rules = _TERMINAL_LENGTH_RULES.get ( ( sym, syms[nb] ) )
            h = None
            if rules:
                for max_len, nh in rules:
                    if d <= max_len:
                        h = nh
                        break
            # [EN] 2026-10-14: an exocyclic N on a flat aromatic ring is an amino group
            # (adenine, cytosine, aminopyridine) down to C-N 1.30 A; the crystal
            # adenine of PLINDER 6BQ4 has C6-N6 1.33 A, read as =NH by the 1.33 rule
            if sym == "N" and syms[nb] == "C" and nb in ring_sp2 and h == 1 and d > 1.30:
                h = 2
            out[aid] = h if h is not None else base - 1
            continue
        angles = [ _angle ( coords[nbs[i]], coords[aid], coords[nbs[j]] ) for i in range ( n ) for j in range ( i + 1, n ) ]
        mean = sum ( angles ) / len ( angles )
        planar = ( sum ( angles ) > 350.0 ) if n == 3 else ( mean > 114.0 )
        short_bond = any ( _dist ( coords[aid], coords[nb] ) <= _MULTIPLE_BOND_MAX.get ( _pair ( sym, syms[nb] ), 0.0 ) for nb in nbs )
        if mean > 155.0 and n == 2:
            pi = 2                     # sp (linear)
        elif aid in ring_sp2:
            pi = 1                     # atom of a flat aromatic-length ring
        elif planar and short_bond:
            pi = 1                     # sp2: planar AND carrying a short (multiple) bond
        elif sym == "C" and ( n == 2 or sum ( angles ) > 340.0 ) and \
                any ( _dist ( coords[aid], coords[nb] ) <= _CLEAR_DOUBLE_MAX.get ( _pair ( sym, syms[nb] ), 0.0 ) for nb in nbs ):
            pi = 1                     # an unmistakable C=X bond, unless the carbon is pyramidal
        elif sym == "C" and n == 3 and planar:
            pi = 1                     # [EN] 2026-10-14: a carbon with three heavy neighbours in a plane is
                                       # sp2 whatever its bond lengths (amidinium C-N 1.32-1.34 A is not
                                       # "short"; an sp3 CH with three heavy neighbours is pyramidal)
        else:
            pi = 0
        if sym == "N" and pi == 1 and n == 3:
            pi = 0                     # planar N with 3 neighbours: amide/aniline/pyrrole lone pair
        if sym in ( "O", "S", "SE" ) and n == 2:
            pi = 0
        if aid in pyrrole_nh:
            pi = 0                     # N-H of a pyrrole/imidazole/indole ring
        out[aid] = max ( 0, base - n - pi )
    return out


# ---------------------------------------------------------------------------
#   step 3: Kekule bond orders (Wang & Case perceiver of vismol)
# ---------------------------------------------------------------------------

def _wang_case ( elements, pairs ):
    try:
        from vismol.core.bond_order_perception import perceive_bond_orders
    except ImportError:
        import os, sys
        here = os.path.dirname ( os.path.abspath ( __file__ ) )
        sys.path.insert ( 0, os.path.join ( here, "..", "graphics_engine", "src" ) )
        from vismol.core.bond_order_perception import perceive_bond_orders
    order_map, _ = perceive_bond_orders ( list ( elements ), list ( pairs ) )
    return order_map


def _saturated ( a, syms, heavy_nb, h_counts ):
    """ True for atoms whose bonds are certainly all single: C with four
        substituents, O/S with two, halogens, and N with three substituents none
        of which is a terminal O/S (nitro and N-oxide N stay free). """
    s = syms[a]
    total = len ( heavy_nb[a] ) + h_counts.get ( a, 0 )
    if s in ( "F", "CL", "BR", "I" ):
        return True
    if s == "C":
        return total >= 4
    if s == "O":
        return total >= 2
    if s == "S":
        # . [EN] 2026-10-14: any 2-connected S (thiol, thioether, thiophene) has
        #   single bonds only; was restricted to S-H, which let the perceiver make
        #   a thioxanthone S+ (PLINDER 6UUX)
        return total == 2
    if s == "N":
        if total < 3 or any ( syms[b] in ( "O", "S" ) and len ( heavy_nb[b] ) == 1 for b in heavy_nb[a] ):
            return False
        # . N with three heavy neighbours, two of them sp2 carbons: may be a
        #   pyridinium-type N+ (thiamine, 7-methylguanine); left to the perceiver,
        #   which keeps it neutral whenever it can
        if h_counts.get ( a, 0 ) == 0 and sum ( 1 for c in heavy_nb[a] if syms[c] == "C" and len ( heavy_nb[c] ) + h_counts.get ( c, 0 ) == 3 ) >= 2:
            return False
        # . N-H of an amidine/guanidine carbon (a C with 3 substituents, 2+ of them
        #   N): both C-N look single in amidinium geometries, so the hydrogen count
        #   may put the C=N on either N; leave them to the perceiver
        if h_counts.get ( a, 0 ) > 0 and any ( syms[c] == "C" and len ( heavy_nb[c] ) + h_counts.get ( c, 0 ) == 3 and
                                              sum ( 1 for d in heavy_nb[c] if syms[d] == "N" ) >= 2 for c in heavy_nb[a] ):
            return False
        return True
    return False


def kekule_orders ( syms, heavy_bonds, h_counts ):
    """ {(i, j): order} for the heavy-heavy bonds, the hydrogens given as
        placeholder atoms (h_counts: {atom: n}).

        [EN] 2026-10-14: bonds touching a certainly saturated atom are fixed as
        single, and the Wang & Case perceiver runs on each remaining conjugated
        system separately, every other neighbour given as a placeholder hydrogen.
        On a whole 40-atom inhibitor in one search it returned C=N+/O- amide
        zwitterions and ring radicals (PL-REX HIV-PR models 16e/16f/22d). """
    heavy_nb = defaultdict ( list )
    for ( i, j ) in heavy_bonds:
        heavy_nb[i].append ( j ); heavy_nb[j].append ( i )
    atoms = { a for b in heavy_bonds for a in b }
    free = { a for a in atoms if not _saturated ( a, syms, heavy_nb, h_counts ) }
    orders = { _key ( i, j ): 1 for ( i, j ) in heavy_bonds }
    # . connected components of the free atoms (through bonds between free atoms)
    seen = set ( )
    for start in sorted ( free ):
        if start in seen:
            continue
        comp, stack = [ ], [ start ]
        seen.add ( start )
        while stack:
            u = stack.pop ( ); comp.append ( u )
            for v in heavy_nb[u]:
                if v in free and v not in seen:
                    seen.add ( v ); stack.append ( v )
        comp_bonds = [ _key ( i, j ) for ( i, j ) in heavy_bonds if i in seen and j in seen and i in comp and j in comp ]
        if not comp_bonds:
            continue
        local = { a: k for k, a in enumerate ( sorted ( comp ) ) }
        elements = [ ( syms[a].capitalize ( ) if len ( syms[a] ) > 1 else syms[a] ) for a in sorted ( comp ) ]
        pairs = [ ( local[i], local[j] ) for ( i, j ) in comp_bonds ]
        for a in comp:
            n_placeholders = h_counts.get ( a, 0 ) + sum ( 1 for b in heavy_nb[a] if b not in local )
            for _ in range ( n_placeholders ):
                elements.append ( "H" )
                pairs.append ( ( local[a], len ( elements ) - 1 ) )
        try:
            order_map = _wang_case ( elements, pairs )
        except Exception:
            order_map = { }
        for ( i, j ) in comp_bonds:
            orders[( i, j )] = int ( order_map.get ( _key ( local[i], local[j] ), 1 ) )
    return orders


# ---------------------------------------------------------------------------
#   step 4: aromaticity
# ---------------------------------------------------------------------------

def _pi_electrons ( sym, ring_orders, exo_double, n_bonds, formal_charge ):
    if 2 in ring_orders:
        return 1
    if sym == "C":
        if formal_charge == -1: return 2
        if formal_charge == 1:  return 0
        if exo_double:          return 0
        return None
    if sym in ( "N", "P" ):
        if formal_charge == 1 and n_bonds == 4:
            return None
        return 2
    if sym in ( "O", "S", "SE" ):
        return 2
    return None


def aromatic_bonds ( syms, orders, h_counts, formal_charges = None ):
    """ Set of heavy-atom bonds (i, j) in aromatic rings (Hueckel 4n+2), fused
        systems included. orders: {(i, j): order} of the heavy-atom bonds. """
    formal_charges = formal_charges or { }
    adjacency = defaultdict ( set )
    for ( i, j ) in orders:
        adjacency[i].add ( j ); adjacency[j].add ( i )
    rings = [ r for r in small_rings ( adjacency ) if 5 <= len ( r ) <= 7 ]

    def order ( a, b ): return orders.get ( _key ( a, b ), 1 )
    def bonds_of ( ring ):
        return { _key ( ring[k], ring[( k + 1 ) % len ( ring )] ) for k in range ( len ( ring ) ) }
    def ring_pi ( ring_atoms, ring_bonds ):
        total = 0
        for a in ring_atoms:
            in_ring = [ order ( a, b ) for b in adjacency[a] if _key ( a, b ) in ring_bonds ]
            exo = any ( order ( a, b ) == 2 for b in adjacency[a] if _key ( a, b ) not in ring_bonds )
            pi = _pi_electrons ( syms[a], in_ring, exo, len ( adjacency[a] ) + h_counts.get ( a, 0 ), formal_charges.get ( a, 0 ) )
            if pi is None:
                return None
            total += pi
        return total

    out = set ( )
    failed = [ ]
    for ring in rings:
        rb = bonds_of ( ring )
        pi = ring_pi ( set ( ring ), rb )
        if pi is not None and pi % 4 == 2:
            out |= rb
        elif pi is not None:
            failed.append ( ( set ( ring ), rb ) )
    for ring_atoms, rb in failed:
        for other in rings:
            ob = bonds_of ( other )
            if rb & ob and set ( other ) != ring_atoms:
                pi = ring_pi ( ring_atoms | set ( other ), rb | ob )
                if pi is not None and pi % 4 == 2:
                    out |= rb | ob
    return out


# ---------------------------------------------------------------------------
#   steps 5-6: protonation and formal charges
# ---------------------------------------------------------------------------

def formal_charges_from_valence ( syms, orders, h_counts ):
    """ {atom: q} from the final valence (bond orders + hydrogens). """
    valence = defaultdict ( int )
    for ( i, j ), o in orders.items ( ):
        valence[i] += o; valence[j] += o
    out = { }
    for a, s in enumerate ( syms ):
        if s == "H":
            continue
        v = valence[a] + h_counts.get ( a, 0 )
        if s == "N":            q = v - 3 if v in ( 2, 4 ) else 0
        elif s == "O":          q = v - 2 if v in ( 1, 3 ) else 0
        elif s in ( "S", "SE" ): q = -1 if v == 1 else ( 1 if v == 3 else 0 )
        elif s == "P":          q = 1 if v == 4 else 0
        elif s == "B":          q = -1 if v == 4 else 0
        else:                   q = 0
        out[a] = q
    return out


def neutralize ( syms, orders, h_counts, coords = None ):
    """ The geometric hydrogen count describes the NEUTRAL molecule; when the
        Kekule assignment still leaves charges, the count was wrong at those atoms
        (an amidine NH2 read as =NH, a carboxylic OH read as =O). Changes h_counts
        in place and returns how many atoms were changed:
          - N-, O-, S- not bonded to a positive atom (nitro and N-oxides keep
            their charge separation): one hydrogen added;
          - imidic acid tautomer O(H)-C=N (acyclic or ring) -> amide: the
            hydrogen moves from O to N. """
    charges = formal_charges_from_valence ( syms, orders, h_counts )
    nb = defaultdict ( list )
    for ( i, j ) in orders:
        nb[i].append ( j ); nb[j].append ( i )
    changed = 0
    for a, q in charges.items ( ):
        if q < 0 and syms[a] in ( "N", "O", "S" ) and not any ( charges.get ( b, 0 ) > 0 for b in nb[a] ):
            h_counts[a] = h_counts.get ( a, 0 ) + 1
            changed += 1
        # . an N+ carrying hydrogens (amidinium read with both C-N as single-bond
        #   lengths): the neutral form has one hydrogen less; the protonation rules
        #   put it back where the pH says so
        elif q > 0 and syms[a] == "N" and h_counts.get ( a, 0 ) > 0 and not any ( charges.get ( b, 0 ) < 0 for b in nb[a] ):
            # [EN] 2026-10-18: the hydrogen comes off the N of the amidine with the
            #   SHORTEST C-N bond (the C=N of the neutral form), not necessarily the
            #   one the Kekule pass made N+ (user: 2-amino-3,4-dihydroquinazoline of
            #   XX01D2VNS2BA-a1291 drawn as ring N-H / exocyclic C=NH; ring C-N 1.31,
            #   exocyclic 1.33 A -> ring C=N, NH2)
            target = a
            centre = [ c for c in nb[a] if syms[c] == "C" and orders.get ( _key ( a, c ), 1 ) == 2 ]
            if coords is not None and centre:
                c = centre[0]
                donors = [ n for n in nb[c] if syms[n] == "N" and h_counts.get ( n, 0 ) > 0 and charges.get ( n, 0 ) >= 0
                           and not any ( charges.get ( b, 0 ) < 0 for b in nb[n] ) ]
                if donors:
                    target = min ( donors, key = lambda n: _dist ( coords[c], coords[n] ) )
            h_counts[target] -= 1
            changed += 1
    for o, s in enumerate ( syms ):
        if s != "O" or h_counts.get ( o, 0 ) != 1 or len ( nb[o] ) != 1:
            continue
        c = nb[o][0]
        if syms[c] != "C":
            continue
        imine = [ n for n in nb[c] if syms[n] == "N" and orders.get ( _key ( c, n ), 1 ) == 2 ]
        # . lactim -> lactam: OH on a ring carbon next to a pyridine-type N (no H,
        #   two heavy neighbours): 2-pyridone, hypoxanthine, guanine, uracil
        lactim = [ n for n in nb[c] if syms[n] == "N" and h_counts.get ( n, 0 ) == 0 and len ( nb[n] ) == 2 ] if len ( nb[c] ) == 3 else [ ]
        target = imine or lactim
        if target:
            h_counts[o] = 0
            h_counts[target[0]] = h_counts.get ( target[0], 0 ) + 1
            changed += 1
    # . radicals: a carbon left with valence 3 was read as sp2 by a flattened
    #   crystal geometry but found no partner for a double bond
    valence = defaultdict ( int )
    for ( i, j ), o in orders.items ( ):
        valence[i] += o; valence[j] += o
    # [EN] 2026-10-14: before adding a hydrogen, (1) a terminal OH on that carbon
    # is a C=O read with a long C-O bond (PLINDER 6UUX thioxanthone); (2) a
    # neighbouring CH/CH2 whose bond angle is open (>= 114 degrees, sp3 ~110) is
    # the partner of a C=C read with a long bond (steroid enone C4=C5, PLINDER 6ITQ)
    for a, s in enumerate ( syms ):
        if s != "C" or not nb[a] or valence[a] + h_counts.get ( a, 0 ) != 3:
            continue
        # . only for a flat (sp2-looking) radical carbon; a pyramidal one (sugar C1
        #   read with a short ring C-O, PLINDER BGC) just gets its hydrogen
        flat = False
        if coords is not None and len ( nb[a] ) >= 2:
            angs = [ _angle ( coords[nb[a][i]], coords[a], coords[nb[a][j]] ) for i in range ( len ( nb[a] ) ) for j in range ( i + 1, len ( nb[a] ) ) ]
            flat = sum ( angs ) > 340.0 if len ( nb[a] ) == 3 else angs[0] >= 114.0
        # . ketone/aldehyde only: a carbon with another O/N (anomeric C1 of a
        #   sugar, PLINDER 3VIJ) is not turned into a lactone/amide
        oh = [ o for o in nb[a] if syms[o] == "O" and len ( nb[o] ) == 1 and h_counts.get ( o, 0 ) == 1 ] if flat else [ ]
        if sum ( 1 for x in nb[a] if syms[x] in ( "O", "N" ) ) > 1:
            oh = [ ]
        partners = [ ]
        if flat:
            for b in nb[a]:
                if syms[b] != "C" or h_counts.get ( b, 0 ) == 0 or len ( nb[b] ) != 2 or \
                   valence[b] + h_counts.get ( b, 0 ) != 4 or any ( orders.get ( _key ( b, x ), 1 ) > 1 for x in nb[b] ):
                    continue
                angle = _angle ( coords[nb[b][0]], coords[b], coords[nb[b][1]] )
                if angle >= 114.0 and _dist ( coords[a], coords[b] ) <= 1.52:   # . not a 1.58 A sp3 bond (PLINDER 2ZSH)
                    partners.append ( ( angle, b ) )
        if oh:
            h_counts[oh[0]] = 0
        elif partners:
            b = max ( partners )[1]
            h_counts[b] -= 1
            valence[b] -= 1                   # . claimed: not offered to another radical
        else:
            h_counts[a] = h_counts.get ( a, 0 ) + 1
        changed += 1
    return changed


class _Graph:
    """ Read-only view of the heavy-atom graph for the protonation rules. """
    def __init__ ( self, syms, orders, h_counts, aromatic ):
        self.syms = syms
        self.orders = orders
        self.h = h_counts
        self.aromatic = aromatic
        self.nb = defaultdict ( list )
        for ( i, j ) in orders:
            self.nb[i].append ( j ); self.nb[j].append ( i )
        self.arom_atoms = { a for b in aromatic for a in b }

    def order ( self, a, b ): return self.orders.get ( _key ( a, b ), 0 )
    def double_to ( self, a, element ): return [ b for b in self.nb[a] if self.syms[b] == element and self.order ( a, b ) == 2 ]
    def single_to ( self, a, element ): return [ b for b in self.nb[a] if self.syms[b] == element and self.order ( a, b ) == 1 ]
    def unsaturated ( self, a ):
        return a in self.arom_atoms or any ( self.order ( a, b ) >= 2 for b in self.nb[a] )
    def is_carbonyl_c ( self, a ):
        return self.syms[a] == "C" and bool ( self.double_to ( a, "O" ) or self.double_to ( a, "S" ) )
    def is_sulfonyl_s ( self, a ):
        return self.syms[a] == "S" and len ( self.double_to ( a, "O" ) ) >= 2
    def terminal_oh ( self, a ):
        """ O atoms singly bonded to a, with a hydrogen and no other heavy neighbour. """
        return [ o for o in self.single_to ( a, "O" ) if self.h.get ( o, 0 ) >= 1 and len ( self.nb[o] ) == 1 ]
    def path_within ( self, a, b, max_bonds ):
        frontier, seen = { a }, { a }
        for _ in range ( max_bonds ):
            frontier = { y for x in frontier for y in self.nb[x] } - seen
            if b in frontier:
                return True
            seen |= frontier
        return False


def protonate ( syms, orders, h_counts, aromatic, ph = 7.4 ):
    """ Changes h_counts in place by functional-group rules at the given pH.
        Returns the list of decisions [ (group, atom, pKa, action) ], action in
        { "deprotonated", "protonated", "kept neutral" }. """
    g = _Graph ( syms, orders, h_counts, aromatic )
    decisions = [ ]
    def acid ( group, atom, pka ):
        if ph > pka and h_counts.get ( atom, 0 ) > 0:
            h_counts[atom] -= 1
            decisions.append ( ( group, atom, pka, "deprotonated" ) )
            return True
        decisions.append ( ( group, atom, pka, "kept neutral" ) )
        return False
    def base ( group, atom, pka ):
        # . a nitrogen that already has four bonds (an N+ left by the perceiver) is
        #   never protonated again
        if sum ( g.order ( atom, b ) for b in g.nb[atom] ) + h_counts.get ( atom, 0 ) >= 4:
            decisions.append ( ( group, atom, pka, "already charged" ) )
            return False
        if ph < pka:
            h_counts[atom] = h_counts.get ( atom, 0 ) + 1
            decisions.append ( ( group, atom, pka, "protonated" ) )
            return True
        decisions.append ( ( group, atom, pka, "kept neutral" ) )
        return False

    for a, s in enumerate ( syms ):
        # . carboxylic acids (and thio analogues C(=S)OH)
        if g.is_carbonyl_c ( a ):
            ohs = g.terminal_oh ( a )
            if ohs and not g.single_to ( a, "N" ):
                acid ( "carboxylic acid", ohs[0], PKA["carboxylic acid"] )
        # . sulfonic acids, sulfuric mono-esters
        elif g.is_sulfonyl_s ( a ):
            for o in g.terminal_oh ( a ):
                acid ( "sulfonic/sulfuric acid", o, PKA["sulfonic/sulfuric acid"] )
        # . phosphoric / phosphonic acids
        elif s == "P" and g.double_to ( a, "O" ):
            ohs = g.terminal_oh ( a )
            if ohs:
                acid ( "phosphate, first OH", ohs[0], PKA["phosphate, first OH"] )
                if len ( ohs ) >= 2:
                    carbon_bonded = any ( syms[b] == "C" for b in g.nb[a] )
                    group = "phosphonate, second OH" if carbon_bonded else "phosphate ester, second OH"
                    acid ( group, ohs[1], PKA[group] )

    # . N-H acids: tetrazoles, acylsulfonamides, sulfonylureas
    adjacency = { a: set ( g.nb[a] ) for a in range ( len ( syms ) ) if syms[a] != "H" }
    for ring in small_rings ( adjacency, max_size = 5 ):
        if len ( ring ) == 5 and sum ( 1 for a in ring if syms[a] == "N" ) == 4 and all ( a in g.arom_atoms for a in ring ):
            nh = [ a for a in ring if syms[a] == "N" and h_counts.get ( a, 0 ) >= 1 ]
            if nh:
                acid ( "tetrazole NH", nh[0], PKA["tetrazole NH"] )
    for a, s in enumerate ( syms ):
        if s != "N" or h_counts.get ( a, 0 ) < 1 or a in g.arom_atoms:
            continue
        sulfonyl = [ b for b in g.nb[a] if g.is_sulfonyl_s ( b ) ]
        carbonyl = [ b for b in g.nb[a] if g.is_carbonyl_c ( b ) ]
        if sulfonyl and carbonyl:
            urea = any ( g.single_to ( c, "N" ) and len ( g.single_to ( c, "N" ) ) >= 2 for c in carbonyl )
            group = "sulfonylurea NH" if urea else "acylsulfonamide NH"
            acid ( group, a, PKA[group] )
        elif len ( sulfonyl ) >= 2:
            acid ( "acylsulfonamide NH", a, PKA["acylsulfonamide NH"] )

    # . amidines and guanidines: protonate the imine N
    for a, s in enumerate ( syms ):
        if s != "C" or a in g.arom_atoms:
            continue
        imine = [ n for n in g.double_to ( a, "N" ) if n not in g.arom_atoms ]
        amino = [ n for n in g.single_to ( a, "N" ) if n not in g.arom_atoms ]
        if len ( imine ) != 1 or not amino:
            continue
        n_im = imine[0]
        # . acylated / sulfonylated / N-oxygenated or cyano amidines are far less basic
        deactivated = any ( g.is_carbonyl_c ( x ) or g.is_sulfonyl_s ( x ) or syms[x] in ( "O", "S" ) or
                            any ( g.order ( x, y ) == 3 for y in g.nb[x] ) or ( syms[x] == "N" and x not in amino )
                            for n in imine + amino for x in g.nb[n] if x != a )
        if deactivated:
            decisions.append ( ( "amidine/guanidine (acyl- or hetero-substituted)", n_im, None, "kept neutral" ) )
            continue
        group = "guanidine" if len ( amino ) >= 2 else "amidine"
        base ( group, n_im, PKA[group] )

    # . aliphatic amines: every neighbour an sp3 carbon (or H); pKa from the base
    #   value and the alpha/beta substituents; one per diamine
    candidates = [ ]
    for a, s in enumerate ( syms ):
        if s != "N" or a in g.arom_atoms or any ( g.order ( a, b ) != 1 for b in g.nb[a] ):
            continue
        if len ( g.nb[a] ) + h_counts.get ( a, 0 ) != 3:
            continue
        if any ( syms[b] != "C" or g.unsaturated ( b ) for b in g.nb[a] ):
            continue
        candidates.append ( ( amine_pka ( g, a ), a ) )
    # . most basic first, so that a diamine protonates its more basic nitrogen
    protonated_amines = [ ]
    for pka, a in sorted ( candidates, key = lambda t: ( -t[0], t[1] ) ):
        if any ( g.path_within ( a, b, 3 ) for b in protonated_amines ):
            base ( "aliphatic amine, second of a diamine", a, pka - DIAMINE_PKA_SHIFT )
            continue
        if base ( "aliphatic amine", a, pka ):
            protonated_amines.append ( a )
    return decisions


def amine_pka ( g, n ):
    """ Estimated pKa of the aliphatic amine N n (Perrin-type increments). """
    syms = g.syms
    pka = AMINE_BASE_PKA.get ( min ( 2, g.h.get ( n, 0 ) ), 9.8 )
    def withdrawing_n ( x ):
        return syms[x] == "N" and ( x in g.arom_atoms or g.unsaturated ( x ) or any ( g.is_carbonyl_c ( y ) or g.is_sulfonyl_s ( y ) for y in g.nb[x] ) )
    for alpha in g.nb[n]:
        for x in g.nb[alpha]:                    # substituents on the alpha carbon
            if x == n:
                continue
            if syms[x] == "C" and ( x in g.arom_atoms or any ( g.order ( x, y ) == 2 and syms[y] == "C" for y in g.nb[x] ) ):
                pka += AMINE_PKA_SHIFT["aryl or vinyl on alpha carbon"]
            elif g.is_carbonyl_c ( x ) and g.single_to ( x, "N" ):
                pka += AMINE_PKA_SHIFT["amide carbonyl on alpha carbon"]
            elif g.is_carbonyl_c ( x ) or g.is_sulfonyl_s ( x ):
                pka += AMINE_PKA_SHIFT["carbonyl on alpha carbon"]
            elif syms[x] == "C" and any ( g.order ( x, y ) == 3 for y in g.nb[x] ):
                pka += AMINE_PKA_SHIFT["nitrile on alpha carbon"]
            elif syms[x] == "F":
                pka += 2 * AMINE_PKA_SHIFT["F on beta carbon"]
            elif syms[x] in ( "O", "S" ) or withdrawing_n ( x ):
                pka += 2 * AMINE_PKA_SHIFT["O or N(sp2) on beta carbon"]
            elif syms[x] == "C":                 # x is a saturated beta carbon: its substituents
                for y in g.nb[x]:
                    if y == alpha:
                        continue
                    if syms[y] == "F":
                        pka += AMINE_PKA_SHIFT["F on beta carbon"]
                    elif syms[y] in ( "O", "S" ) or withdrawing_n ( y ):
                        pka += AMINE_PKA_SHIFT["O or N(sp2) on beta carbon"]
    return pka


# ---------------------------------------------------------------------------
#   step 7: hydrogen placement (Builder atom_ops VSEPR construction, GUI-free)
# ---------------------------------------------------------------------------
_TETRAHEDRAL = [ v / np.linalg.norm ( v ) for v in (
    np.array ( [ -0.785298,  0.243518, -0.653254 ] ), np.array ( [ 0.322015, -0.981331, -0.189814 ] ),
    np.array ( [ -0.334691,  0.073016,  0.992665 ] ), np.array ( [ 0.798227,  0.665009, -0.149645 ] ) ) ]


def _perpendicular ( v ):
    trial = np.array ( [ 1.0, 0.0, 0.0 ] ) if abs ( v[0] ) < 0.9 else np.array ( [ 0.0, 1.0, 0.0 ] )
    p = np.cross ( v, trial )
    return p / np.linalg.norm ( p )


def _fan ( existing, n_new, angle ):
    if existing:
        avg = -np.sum ( existing, axis = 0 )
        norm = float ( np.linalg.norm ( avg ) )
        avg = avg / norm if norm > 1e-6 else np.array ( [ 0.0, 0.0, 1.0 ] )
    else:
        avg = np.array ( [ 0.0, 0.0, 1.0 ] )
    if n_new == 1:
        return [ avg ]
    perp = _perpendicular ( avg )
    spread = math.radians ( angle or 109.5 )
    out = [ ]
    for k in range ( n_new ):
        off = ( k - ( n_new - 1 ) / 2.0 ) * ( spread / max ( 1, n_new - 1 ) )
        d = avg * math.cos ( off ) + perp * math.sin ( off )
        out.append ( d / np.linalg.norm ( d ) )
    return out


def vsepr_directions ( existing, total, angle ):
    """ New unit vectors completing `total` substituents around an atom whose
        `existing` bond directions are fixed, at the given ideal angle (degrees);
        sp2 completions stay in the plane of the existing bonds. """
    existing = [ np.asarray ( d, float ) / np.linalg.norm ( d ) for d in existing ]
    n_ex, n_new = len ( existing ), total - len ( existing )
    if n_new <= 0:
        return [ ]
    theta = math.radians ( angle ) if angle else None
    try:
        if n_ex == 0:
            if total == 1:
                return [ np.array ( [ 0.0, 0.0, 1.0 ] ) ]
            if total == 2 and theta is not None:
                h = theta / 2.0
                return [ np.array ( [ math.sin ( h ), 0.0, math.cos ( h ) ] ), np.array ( [ -math.sin ( h ), 0.0, math.cos ( h ) ] ) ]
            if total == 3 and angle and abs ( angle - 120.0 ) < 1e-6:
                return [ np.array ( [ math.cos ( math.radians ( 120 * k ) ), math.sin ( math.radians ( 120 * k ) ), 0.0 ] ) for k in range ( 3 ) ]
            if total in ( 3, 4 ):
                return list ( _TETRAHEDRAL[:total] )
        elif n_ex == 1 and theta is not None:
            A = existing[0]
            e1 = _perpendicular ( A ); e2 = np.cross ( A, e1 ); e2 /= np.linalg.norm ( e2 )
            if n_new == 1:
                return [ math.cos ( theta ) * A + math.sin ( theta ) * e1 ]
            if n_new == 2:
                phi = math.acos ( max ( -1.0, min ( 1.0, math.cos ( theta ) / ( 1.0 + math.cos ( theta ) ) ) ) )
                return [ math.cos ( theta ) * A + math.sin ( theta ) * ( math.cos ( az ) * e1 + math.sin ( az ) * e2 ) for az in ( 0.0, phi ) ]
            if n_new == 3:
                return [ math.cos ( theta ) * A + math.sin ( theta ) * ( math.cos ( az ) * e1 + math.sin ( az ) * e2 )
                         for az in ( 0.0, 2 * math.pi / 3, 4 * math.pi / 3 ) ]
        elif n_ex == 2 and theta is not None:
            A, B = existing
            u = A + B; un = float ( np.linalg.norm ( u ) )
            c = np.cross ( A, B ); cn = float ( np.linalg.norm ( c ) )
            if un < 1e-3 or cn < 1e-3:
                raise ValueError
            u /= un; nrm = c / cn
            half = math.cos ( math.acos ( max ( -1.0, min ( 1.0, float ( np.dot ( A, B ) ) ) ) ) / 2.0 )
            if abs ( half ) < 1e-3:
                raise ValueError
            if n_new == 1 and abs ( angle - 120.0 ) < 1e-6:
                return [ -u ]                      # in the plane, bisecting the reflex angle
            if n_new == 1:
                a = math.cos ( theta ) / half
                if abs ( a ) > 1.0: raise ValueError
                return [ a * u + math.sqrt ( max ( 0.0, 1 - a * a ) ) * nrm ]
            if n_new == 2:
                p = math.cos ( theta ) / half
                if abs ( p ) > 1.0: raise ValueError
                r = math.sqrt ( max ( 0.0, 1 - p * p ) )
                return [ p * u + r * nrm, p * u - r * nrm ]
        elif n_ex == 3 and n_new == 1:
            s = -np.sum ( existing, axis = 0 ); sn = float ( np.linalg.norm ( s ) )
            if sn < 1e-3: raise ValueError
            return [ s / sn ]
    except ValueError:
        pass
    return _fan ( existing, n_new, angle )


def _ideal_angle ( total, max_order, double_count, planar ):
    if total < 2 or total > 4:
        return None
    if total == 2:
        if max_order == 3 or double_count == 2: return 180.0
        return 120.0 if ( max_order == 2 or planar ) else 109.5
    if total == 3:
        return 120.0 if ( max_order >= 2 or planar ) else 109.5
    return 109.5


def place_hydrogens ( syms, coords, orders, h_counts, aromatic ):
    """ [ (heavy atom, xyz) ] for every hydrogen in h_counts. Aromatic atoms and
        N/O conjugated with an unsaturated neighbour (amide, aniline, enol O) are
        completed in the plane. """
    g = _Graph ( syms, orders, h_counts, aromatic )
    out = [ ]
    for a in sorted ( h_counts ):
        nh = h_counts[a]
        if nh <= 0:
            continue
        pos = np.asarray ( coords[a], float )
        dirs = [ np.asarray ( coords[b], float ) - pos for b in g.nb[a] ]
        heavy_orders = [ g.order ( a, b ) for b in g.nb[a] ]
        total = len ( dirs ) + nh
        max_order = max ( heavy_orders ) if heavy_orders else 1
        doubles = sum ( 1 for o in heavy_orders if o == 2 )
        conjugated = syms[a] == "N" and total == 3 and any ( g.unsaturated ( b ) for b in g.nb[a] )
        planar = a in g.arom_atoms or conjugated
        angle = _ideal_angle ( total, max_order, doubles, planar )
        length = XH_LENGTH.get ( syms[a], DEFAULT_XH_LENGTH )
        for d in vsepr_directions ( dirs, total, angle ):
            out.append ( ( a, pos + length * np.asarray ( d ) ) )
    return out


# ---------------------------------------------------------------------------
#   driver
# ---------------------------------------------------------------------------

def perceive_ligand ( symbols, coords, bonds = None, ph = 7.4, protonate_groups = True, h_counts = None ):
    """ Chemistry of one residue (ligand) from its coordinates.

        symbols: element symbols; coords: (N, 3) A; bonds: optional [(i, j)] or
        [(i, j, order)] -- orders are ignored and re-perceived, the pairs are
        used as the connectivity.

        Heavy atoms only (no hydrogen in the input): neutral hydrogen counts from
        the geometry, Kekule orders, aromaticity, protonation at pH, hydrogens
        placed. With hydrogens in the input: they are kept, orders perceived, no
        protonation change.

        Returns { "symbols", "coords" (input atoms first, new H appended),
                  "bonds" [(i, j, order)], "aromatic" {(i, j)}, "formal_charges",
                  "total_charge", "added_hydrogens" [(heavy atom, index)],
                  "decisions" [(group, atom, pKa, action)], "notes" }.

        h_counts: {heavy atom index: hydrogen count} chosen by the user for an
        H-less input (2026-10-18, protonation review): used instead of the
        geometric count / neutral repair / pH rules; orders and formal charges
        follow from it. """
    syms = [ _sym ( s ) for s in symbols ]
    xyz = [ np.asarray ( c, float ) for c in coords ]
    n = len ( syms )
    if bonds is None:
        pairs = connectivity ( syms, xyz )
    else:
        pairs = sorted ( { _key ( int ( b[0] ), int ( b[1] ) ) for b in bonds } )
    notes = [ ]
    has_h = any ( s == "H" for s in syms )
    heavy_pairs = [ p for p in pairs if syms[p[0]] != "H" and syms[p[1]] != "H" ]
    heavy_nb = { a: [ ] for a in range ( n ) if syms[a] != "H" }
    for ( i, j ) in heavy_pairs:
        heavy_nb[i].append ( j ); heavy_nb[j].append ( i )
    override = None if ( has_h or h_counts is None ) else { int ( a ): int ( c ) for a, c in h_counts.items ( ) }
    if override is not None:
        h_counts = { a: max ( 0, override.get ( a, 0 ) ) for a in heavy_nb }
        notes.append ( "hydrogen counts chosen by the user; bond orders perceived (Wang & Case)" )
    elif has_h:
        h_counts = defaultdict ( int )
        h_of = defaultdict ( list )
        for ( i, j ) in pairs:
            for h, x in ( ( i, j ), ( j, i ) ):
                if syms[h] == "H" and syms[x] != "H":
                    h_counts[x] += 1; h_of[x].append ( h )
        h_counts = dict ( h_counts )
        notes.append ( "hydrogens from the input; bond orders perceived (Wang & Case)" )
    else:
        h_counts = geometric_implicit_hydrogens ( syms, heavy_nb, xyz )
        notes.append ( "neutral hydrogen counts from the geometry; bond orders perceived (Wang & Case)" )
    orders = kekule_orders ( syms, heavy_pairs, h_counts )
    if not has_h and override is None:
        changed = neutralize ( syms, orders, h_counts, xyz )
        if changed:
            orders = kekule_orders ( syms, heavy_pairs, h_counts )
            notes.append ( "neutral form repaired on {} atom(s) after the geometric hydrogen count".format ( changed ) )
    charges0 = formal_charges_from_valence ( syms, orders, h_counts )
    arom = aromatic_bonds ( syms, orders, h_counts, charges0 )
    decisions = [ ]
    if not has_h and protonate_groups and override is None:
        decisions = protonate ( syms, orders, h_counts, arom, ph = ph )
        notes.append ( "protonation at pH {:.1f} by group rules ({} decision(s))".format ( ph, len ( decisions ) ) )
    charges = formal_charges_from_valence ( syms, orders, h_counts )
    out_syms = list ( syms )
    out_xyz = list ( xyz )
    out_bonds = [ ( i, j, orders[( i, j )] ) for ( i, j ) in heavy_pairs ]
    added = [ ]
    if has_h:
        out_bonds += [ ( i, j, 1 ) for ( i, j ) in pairs if syms[i] == "H" or syms[j] == "H" ]
    else:
        for heavy, pos in place_hydrogens ( syms, xyz, orders, h_counts, arom ):
            out_syms.append ( "H" ); out_xyz.append ( pos )
            k = len ( out_syms ) - 1
            out_bonds.append ( ( heavy, k, 1 ) )
            added.append ( ( heavy, k ) )
    q = [ charges.get ( a, 0 ) for a in range ( len ( out_syms ) ) ]
    return { "symbols": out_syms, "coords": np.array ( out_xyz ), "bonds": out_bonds, "aromatic": arom,
             "formal_charges": q, "total_charge": int ( sum ( q ) ), "added_hydrogens": added,
             "decisions": decisions, "notes": notes }


def format_decisions ( result, names = None ):
    """ One line per protonation decision, for reports. """
    lines = [ ]
    for group, atom, pka, action in result.get ( "decisions", [ ] ):
        label = names[atom] if names else "{}{}".format ( result["symbols"][atom], atom + 1 )
        lines.append ( "{}: {} ({}{})".format ( label, action, group, "" if pka is None else ", pKa ~{:.1f}".format ( pka ) ) )
    return lines
