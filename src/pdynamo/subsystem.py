#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: capped subsystem around a selection (xTB refinement of part of a system)
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      Graph-only (no RDKit, no atom typing: metal centres are fine) version of
#      the fraglib fragmenter, for refining a SELECTION of a large system with
#      xTB (DYFF Parameters > Refine with xTB, "selected atoms only"):
#
#        region = selection + every atom within `shell` bonds of it, grown until
#                 stable so that
#                   - a bond that is not single (double, triple, aromatic) is
#                     never cut;
#                   - a ring bond is never cut (ring systems enter whole);
#                   - a metal inside the region keeps all its ligands, and a
#                     metal bonded to the region is pulled in (a donor atom
#                     is never capped with H in place of its metal);
#                   - the hydrogens of every kept heavy atom are kept;
#        caps   = every cut bond A(in)-B(out) becomes A-H along A->B, at a
#                 standard X-H length;
#        charge = sum of the formal charges of the region atoms.
#
#      Only terms whose atoms are ALL in the selection are refined (the shell
#      gives them their chemical environment; terms touching the shell or the
#      caps would carry the cut artifacts).
#

import numpy as np

CAP_LENGTH = { 6: 1.09, 7: 1.01, 8: 0.96, 16: 1.34, 15: 1.42, 5: 1.19, 14: 1.48 }
_METALS = set ( range ( 21, 31 ) ) | set ( range ( 39, 49 ) ) | set ( range ( 57, 81 ) ) | { 3, 4, 11, 12, 13, 19, 20, 37, 38, 55, 56, 81, 82, 83 }


def _neighbours ( n, bonds ):
    nb = [ set ( ) for _ in range ( n ) ]
    for b in bonds:
        a, c = int ( b[0] ), int ( b[1] )
        nb[a].add ( c ); nb[c].add ( a )
    return nb


def ring_bonds ( n, bonds, max_size = None ):
    """ Set of ( a, b ) (a < b) bonds that belong to a cycle (bridges excluded);
        with max_size, only the bonds of a cycle of at most max_size atoms. """
    nb = _neighbours ( n, bonds )
    index, low, counter, bridges = [ -1 ] * n, [ 0 ] * n, [ 0 ], set ( )
    for root in range ( n ):
        if index[root] != -1: continue
        stack = [ ( root, -1, iter ( sorted ( nb[root] ) ) ) ]
        index[root] = low[root] = counter[0]; counter[0] += 1
        while stack:
            u, parent, it = stack[-1]
            advanced = False
            for v in it:
                if v == parent: continue
                if index[v] == -1:
                    index[v] = low[v] = counter[0]; counter[0] += 1
                    stack.append ( ( v, u, iter ( sorted ( nb[v] ) ) ) ); advanced = True; break
                low[u] = min ( low[u], index[v] )
            if not advanced:
                stack.pop ( )
                if stack:
                    p = stack[-1][0]
                    low[p] = min ( low[p], low[u] )
                    if low[u] > index[p]: bridges.add ( ( min ( p, u ), max ( p, u ) ) )
    cyclic = { ( min ( int ( a ), int ( b ) ), max ( int ( a ), int ( b ) ) ) for a, b, *_ in bonds } - bridges
    if max_size is None: return cyclic
    small = set ( )
    for a, b in cyclic:
        # . shortest path a -> b without the bond itself, at most max_size - 1 bonds
        depth, frontier, seen = 0, { a }, { a }
        while frontier and depth < max_size - 1 and b not in seen:
            frontier = { v for u in frontier for v in nb[u] if not ( u == a and v == b ) and v not in seen }
            seen |= frontier; depth += 1
        if b in seen: small.add ( ( a, b ) )
    return small


#  rings that must not be cut: aromatic and ordinary rings up to this size. Larger
#  cycles are cut and capped like chains -- in a protein the disulfide bridges close
#  macrocycles through the backbone (crambin: a 2-shell region grew to 331 atoms)
SMALL_RING = 8


def build_subsystem ( atomic_numbers, coordinates, bonds, selection, shell = 2, formal_charges = None ):
    """ bonds: [ ( a, b, order, aromatic ) ] of the whole system (order 1/2/3).
        Returns a dict: atoms (parent indices, region order), parent_index (per
        subsystem atom, -1 = cap), atomic_numbers, coordinates, bonds [ ( a, b,
        order, aromatic ) ] in subsystem indices, selection_local, charge,
        n_caps, whole (True when the region is the whole system). """
    z = [ int ( x ) for x in atomic_numbers ]
    n = len ( z )
    xyz = np.asarray ( coordinates, dtype = float )
    nb = _neighbours ( n, bonds )
    info = { ( min ( int ( b[0] ), int ( b[1] ) ), max ( int ( b[0] ), int ( b[1] ) ) ):
             ( int ( b[2] ) if len ( b ) > 2 else 1, bool ( b[3] ) if len ( b ) > 3 else False ) for b in bonds }
    in_ring = ring_bonds ( n, bonds, max_size = SMALL_RING )
    selection = { int ( a ) for a in selection }
    region = set ( selection )
    frontier = set ( selection )
    for _ in range ( max ( 0, int ( shell ) ) ):
        frontier = { v for u in frontier for v in nb[u] } - region
        region |= frontier
    changed = True
    while changed:
        changed = False
        for a in list ( region ):
            for b in nb[a]:
                if b in region: continue
                key = ( min ( a, b ), max ( a, b ) )
                order, aromatic = info.get ( key, ( 1, False ) )
                # . never cut: multiple/aromatic or ring bonds, the bonds of a metal of the
                #   region, a bond to a metal (a donor is not capped with H), an X-H bond
                if order != 1 or aromatic or key in in_ring or z[a] in _METALS or z[b] in _METALS or z[b] == 1:
                    region.add ( b ); changed = True
        for a in list ( region ):
            if z[a] != 1:
                for b in nb[a]:
                    if z[b] == 1 and b not in region: region.add ( b ); changed = True
    atoms = sorted ( region )
    local = { a: k for k, a in enumerate ( atoms ) }
    sub_z = [ z[a] for a in atoms ]
    sub_xyz = [ xyz[a] for a in atoms ]
    parent = list ( atoms )
    sub_bonds = [ ]
    for ( a, b ), ( order, aromatic ) in info.items ( ):
        if a in local and b in local:
            sub_bonds.append ( ( local[a], local[b], order, aromatic ) )
    n_caps = 0
    for a in atoms:
        for b in nb[a]:
            if b in local: continue
            h = len ( sub_z )
            sub_z.append ( 1 ); parent.append ( -1 )
            d = xyz[b] - xyz[a]; d /= np.linalg.norm ( d )
            sub_xyz.append ( xyz[a] + d * CAP_LENGTH.get ( z[a], 1.0 ) )
            sub_bonds.append ( ( local[a], h, 1, False ) )
            n_caps += 1
    charge = sum ( int ( formal_charges[a] ) for a in atoms ) if formal_charges is not None else 0
    return { "atoms": atoms, "parent_index": parent, "atomic_numbers": sub_z, "coordinates": np.array ( sub_xyz ),
             "bonds": sub_bonds, "selection_local": sorted ( local[a] for a in selection ), "charge": charge,
             "n_caps": n_caps, "whole": len ( atoms ) == n }


def map_term_edits_to_parent ( edits, parent_index, selection ):
    """ Per-term edits computed on the subsystem -> keys in parent indices, keeping
        only the terms whose atoms are all in the selection (parent indices). """
    selection = set ( selection )
    out = { "bond_terms": { }, "angle_terms": { }, "dihedral_bonds": { } }
    def remap ( label ):
        idx = [ parent_index[int ( x )] for x in label.split ( "-" ) ]
        return idx if all ( i >= 0 and i in selection for i in idx ) else None
    for key, value in ( edits.get ( "bond_terms" ) or { } ).items ( ):
        idx = remap ( key )
        if idx: out["bond_terms"]["{}-{}".format ( min ( idx ), max ( idx ) )] = value
    for key, value in ( edits.get ( "angle_terms" ) or { } ).items ( ):
        idx = remap ( key )
        if idx:
            a, b, c = idx
            out["angle_terms"]["{}-{}-{}".format ( min ( a, c ), b, max ( a, c ) )] = value
    for key, value in ( edits.get ( "dihedral_bonds" ) or { } ).items ( ):
        idx = remap ( key )
        if idx: out["dihedral_bonds"]["{}-{}".format ( min ( idx ), max ( idx ) )] = value
    return out


def map_term_edits_to_subsystem ( edits, parent_index ):
    """ Existing per-term edits of the parent (parent indices) -> subsystem keys
        (terms whose atoms are all real atoms of the subsystem). """
    local = { p: k for k, p in enumerate ( parent_index ) if p >= 0 }
    out = { }
    for kind in ( "bond_terms", "angle_terms", "dihedral_bonds" ):
        table = { }
        for key, value in ( ( edits or { } ).get ( kind ) or { } ).items ( ):
            idx = [ int ( x ) for x in key.split ( "-" ) ]
            if all ( i in local for i in idx ):
                l = [ local[i] for i in idx ]
                if kind == "angle_terms":
                    a, b, c = l; key2 = "{}-{}-{}".format ( min ( a, c ), b, max ( a, c ) )
                else:
                    key2 = "{}-{}".format ( min ( l ), max ( l ) )
                table[key2] = value
        if table: out[kind] = table
    return out


def pdynamo_system ( atomic_numbers, bonds, coordinates, formal_charges = None, label = "subsystem" ):
    """ pDynamo System of a subsystem (same construction as the Builder's
        _build_pdynamo_system_from_vismol_object: element, formal charge, bond
        order and aromaticity -- what DYFF types from). """
    from pMolecule import Atom, Bond, BondType, Connectivity, ConvertInputConnectivity, System
    from pScientific.Geometry3 import Coordinates3
    connectivity = Connectivity ( )
    for k, z in enumerate ( atomic_numbers ):
        q = int ( formal_charges[k] ) if formal_charges is not None else 0
        connectivity.AddNode ( Atom.WithOptions ( atomicNumber = int ( z ), label = "A{}".format ( k + 1 ), formalCharge = q ) )
    types = { 1: BondType.Single, 2: BondType.Double, 3: BondType.Triple }
    for a, b, order, aromatic in bonds:
        connectivity.AddEdge ( Bond.WithNodes ( connectivity.nodes[int ( a )], connectivity.nodes[int ( b )],
                                                isAromatic = bool ( aromatic ), type = types.get ( int ( order ), BondType.Single ) ) )
    ConvertInputConnectivity ( connectivity, { } )
    system = System.FromConnectivity ( connectivity = connectivity )
    system.label = label
    c3 = Coordinates3.WithExtent ( len ( atomic_numbers ) )
    for a, p in enumerate ( coordinates ):
        for k in range ( 3 ): c3[a, k] = float ( p[k] )
    system.coordinates3 = c3
    return system

