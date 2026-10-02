#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  residue_chemistry.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  [EN] 2026-10-02 user request: a PDB structure imported by pDynamo and then
#  edited in the Builder must get its atoms' chemical nature right (bond
#  orders, aromatic rings, hybridization) and its hydrogens added correctly,
#  each hydrogen INSIDE its own residue. User's choices: residue templates
#  from pDynamo for standard residues; perception only for the rest.
#
#  - assign_bond_orders(): standard residues take their bond orders straight
#    from pDynamo's PDB component library (the same library util/
#    hydrogen_builder.py uses -- including its CYM/LYN additions), with the
#    protonation variant that matches the hydrogens the residue already has
#    (or the library default). Peptide C-N and disulfide S-S links are single.
#    Everything else (ligands, cofactors, Builder-made molecules) goes to the
#    Wang & Case perceiver one FRAGMENT at a time, each bond leaving the
#    fragment replaced by a placeholder hydrogen so the fragment's valences
#    stay right. (Running the perceiver on a whole protein at once used to
#    hang "Edit in Builder".)
#  - perceive_aromatic_bonds(): Hueckel (4n+2) rings, fused systems included,
#    from the Kekule bond orders above.
#  - atom_hybridization(): sp / sp2 / sp3 / aromatic per atom, with
#    conjugated N (amide, aniline, guanidine...) reported as sp2.
# ============================================================================
from collections import defaultdict


# residue names that are linked to the next/previous residue by a peptide
# bond (C of residue i, N of residue i+1) -- anything with a template
PEPTIDE_ATOMS = ( "C", "N" )


def _library ( ):
    from util.hydrogen_builder import _get_library
    return _get_library ( )


def _component_variants ( library, component_label ):
    """ Every variant the library defines specifically for this component
        (e.g. HIS: Delta/Epsilon/Doubly Protonated; ASP: Deprotonated). """
    out = [ ]
    for variant in library.variants.values ( ):
        if getattr ( variant, "componentLabel", None ) == component_label:
            out.append ( variant )
    return out


def template_bond_info ( component_label, variant_label = "__default__" ):
    """ ( {frozenset((name1, name2)): order}, {atom_name: formal_charge},
          {atom_name: atomic_number}, applied_variant_label )
        for a standard residue, or None when the library has no such
        component. variant_label: "__default__" = the component's own default
        variant (HIS -> Delta Protonated, ASP/GLU -> Deprotonated, others ->
        plain template); None = the plain template; anything else = that
        variant. Applies the variant exactly like pBabel's own
        ApplyLibraryVariant: bondsToDelete, atomsToDelete, atomsToAdd,
        bondsToAdd, then bondTypes and formalCharges overrides. """
    library = _library ( )
    component = library.GetComponent ( component_label )
    if component is None:
        return None
    # single-atom library components (ions: CL, NA, ...) have bonds = None
    atoms   = { a.label: a.atomicNumber for a in ( component.atoms or [ ] ) }
    charges = { a.label: int ( getattr ( a, "formalCharge", 0 ) or 0 ) for a in ( component.atoms or [ ] ) }
    orders  = { }
    for b in ( component.bonds or [ ] ):
        orders[frozenset ( ( b.atomLabel1, b.atomLabel2 ) )] = _bond_type_to_order ( b.bondType )

    if variant_label == "__default__":
        defaults = getattr ( component, "variants", None )
        variant_label = defaults[0] if defaults else None
    applied = None
    if variant_label is not None:
        variant = library.GetVariant ( variant_label, component_label )
        if variant is not None:
            applied = variant_label
            for pair in ( variant.bondsToDelete or [ ] ):
                orders.pop ( frozenset ( pair ), None )
            for label in ( variant.atomsToDelete or [ ] ):
                atoms.pop ( label, None ); charges.pop ( label, None )
                for key in [ k for k in orders if label in k ]:
                    orders.pop ( key )
            for atom in ( variant.atomsToAdd or [ ] ):
                atoms[atom.label] = atom.atomicNumber
                charges[atom.label] = int ( getattr ( atom, "formalCharge", 0 ) or 0 )
            for b in ( variant.bondsToAdd or [ ] ):
                orders[frozenset ( ( b.atomLabel1, b.atomLabel2 ) )] = _bond_type_to_order ( b.bondType )
            for b in ( getattr ( variant, "bondTypes", None ) or [ ] ):
                key = frozenset ( ( b.atomLabel1, b.atomLabel2 ) )
                if key in orders:
                    orders[key] = _bond_type_to_order ( b.bondType )
            for label, q in ( variant.formalCharges or { } ).items ( ):
                charges[label] = int ( q )
    return orders, charges, atoms, applied


def _bond_type_to_order ( bond_type ):
    name = str ( bond_type ).split ( "." )[-1].lower ( )
    return { "single": 1, "double": 2, "triple": 3 }.get ( name, 1 )


def best_variant_for_residue ( component_label, atom_names ):
    """ The protonation variant whose hydrogen set best matches the H names
        the residue already has (None / plain template included); with no
        hydrogens at all, the library default ("__default__"). """
    h_present = { n for n in atom_names if n.startswith ( "H" ) }
    if not h_present:
        return "__default__"
    library = _library ( )
    candidates = [ None ] + [ v.label for v in _component_variants ( library, component_label ) ]
    best, best_score = "__default__", None
    for label in candidates:
        info = template_bond_info ( component_label, label )
        if info is None:
            return "__default__"
        template_h = { n for n in info[2] if n.startswith ( "H" ) }
        score = len ( template_h ^ h_present )     # symmetric difference: lower is better
        if best_score is None or score < best_score:
            best, best_score = label, score
    return best


# ---------------------------------------------------------------------------
#   bond orders
# ---------------------------------------------------------------------------

def assign_bond_orders ( atoms, pairs, coords = None ):
    """ atoms: {atom_id: vismol Atom}; pairs: iterable of (i, j) atom_id
        pairs (i < j); coords (optional) {atom_id: (x, y, z)} -- used for
        fragments that have no hydrogens at all (geometry-derived implicit
        H, see geometric_implicit_hydrogens()). Returns ( orders {pair: 1|2|3}, formal_charges
        {atom_id: q}, report dict ) covering EVERY pair.

        report: {"template_residues": n, "perceived_fragments": n,
                 "perceived_atoms": n, "variants": {(chain, resi, resn): variant}} """
    pairs = [ ( min ( i, j ), max ( i, j ) ) for ( i, j ) in pairs ]
    residue_of = { }
    for aid, atom in atoms.items ( ):
        residue = getattr ( atom, "residue", None )
        residue_of[aid] = residue

    # templates per residue object
    template_of_residue = { }
    variants = { }
    for residue in { r for r in residue_of.values ( ) if r is not None }:
        names = [ a.name for a in residue.atoms.values ( ) ]
        variant = best_variant_for_residue ( residue.name, names )
        info = template_bond_info ( residue.name, variant )
        if info is None:
            continue
        template_of_residue[residue] = info
        chain_name = getattr ( getattr ( residue, "chain", None ), "name", "" )
        variants[( chain_name, residue.index, residue.name )] = info[3]

    orders = { }
    charges = { }
    unresolved = [ ]
    for ( i, j ) in pairs:
        ai, aj = atoms[i], atoms[j]
        ri, rj = residue_of.get ( i ), residue_of.get ( j )
        if ai.symbol == "H" or aj.symbol == "H":
            orders[( i, j )] = 1
            continue
        if ri is not None and ri is rj and ri in template_of_residue:
            order = template_of_residue[ri][0].get ( frozenset ( ( ai.name, aj.name ) ) )
            if order is not None:
                orders[( i, j )] = order
                continue
        if ( ri in template_of_residue and rj in template_of_residue and ri is not rj ):
            # peptide link, disulfide, or any other bond between two
            # template residues: a plain single bond
            orders[( i, j )] = 1
            continue
        unresolved.append ( ( i, j ) )

    for residue, info in template_of_residue.items ( ):
        names_to_charge = info[1]
        for atom in residue.atoms.values ( ):
            q = names_to_charge.get ( atom.name, 0 )
            if q and _residue_has_hydrogens ( residue ):
                # only trust a template charge when the residue's own H set
                # was used to pick the variant (no H -> protonation unknown)
                charges[atom.atom_id] = q

    perceived_fragments, perceived_atoms = _perceive_fragments ( atoms, pairs, unresolved, orders, coords )
    report = { "template_residues": len ( template_of_residue ),
               "perceived_fragments": perceived_fragments,
               "perceived_atoms": perceived_atoms,
               "variants": variants }
    return orders, charges, report


def _residue_has_hydrogens ( residue ):
    return any ( a.symbol == "H" for a in residue.atoms.values ( ) )


# [EN] 2026-10-02: geometry rules for fragments that come WITHOUT hydrogens
# (a heavy-atom-only PDB ligand): the perceiver needs every atom's real
# valence, and with no H a methyl C (1 bond) looks like it "needs" a triple
# bond. Hybridization is read from the geometry instead -- bond angles for
# atoms with 2+ heavy neighbours, bond length for terminal atoms -- and turned
# into implicit hydrogen counts, given to the perceiver as placeholder H.
_BASE_VALENCE = { "C": 4, "N": 3, "O": 2, "S": 2, "P": 3, "B": 3, "SI": 4, "SE": 2 }
# terminal atom: (symbol, neighbour symbol) -> [(max length, implicit H), ...]
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


def _angle ( a, b, c ):
    import math
    v1 = [ a[k] - b[k] for k in range ( 3 ) ]
    v2 = [ c[k] - b[k] for k in range ( 3 ) ]
    n1 = math.sqrt ( sum ( x * x for x in v1 ) ); n2 = math.sqrt ( sum ( x * x for x in v2 ) )
    if n1 < 1e-6 or n2 < 1e-6:
        return 109.5
    cosang = max ( -1.0, min ( 1.0, sum ( v1[k] * v2[k] for k in range ( 3 ) ) / ( n1 * n2 ) ) )
    return math.degrees ( math.acos ( cosang ) )


def _torsion ( a, b, c, d ):
    import math
    def sub ( u, v ): return [ u[k] - v[k] for k in range ( 3 ) ]
    def cross ( u, v ): return [ u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0] ]
    def dot ( u, v ): return sum ( u[k] * v[k] for k in range ( 3 ) )
    b0, b1, b2 = sub ( b, a ), sub ( c, b ), sub ( d, c )     # a->b, b->c, c->d
    n1, n2 = cross ( b0, b1 ), cross ( b1, b2 )
    nb1 = math.sqrt ( dot ( b1, b1 ) ) or 1.0
    m1 = cross ( n1, [ x / nb1 for x in b1 ] )
    return math.degrees ( math.atan2 ( dot ( m1, n2 ), dot ( n1, n2 ) ) )


# bond shorter than this -> evidence of a multiple bond on that atom
_MULTIPLE_BOND_MAX = { ( "C", "C" ): 1.43, ( "C", "N" ): 1.31, ( "C", "O" ): 1.30, ( "C", "S" ): 1.72,
                       ( "N", "N" ): 1.30, ( "N", "O" ): 1.30, ( "O", "S" ): 1.60, ( "O", "P" ): 1.60,
                       ( "N", "S" ): 1.60 }
# every bond of an aromatic ring is at most this long (saturated CH2-CH2,
# C-N single bonds of lactams/ureas are longer)
_AROMATIC_BOND_MAX = { ( "C", "C" ): 1.45, ( "C", "N" ): 1.40, ( "N", "N" ): 1.40, ( "C", "O" ): 1.40,
                       ( "C", "S" ): 1.78, ( "N", "S" ): 1.70, ( "N", "O" ): 1.40 }


def _pair ( a, b ):
    return tuple ( sorted ( ( a.upper ( ), b.upper ( ) ) ) )


def _dist ( p, q ):
    import math
    return math.sqrt ( sum ( ( p[k] - q[k] ) ** 2 for k in range ( 3 ) ) )


def _planar_rings ( neighbors, coords, atoms = None, max_twist = 15.0 ):
    """ 5- and 6-membered rings whose internal torsions are all within
        max_twist degrees of zero (flat) and -- when `atoms` is given --
        whose bonds all have aromatic lengths (see _AROMATIC_BOND_MAX). """
    out = [ ]
    adjacency = { a: set ( nbs ) for a, nbs in neighbors.items ( ) }
    for ring in _small_rings ( adjacency, max_size = 6 ):
        if len ( ring ) not in ( 5, 6 ):
            continue
        n = len ( ring )
        twists = [ abs ( _torsion ( *[ coords[ring[( k + d ) % n]] for d in range ( 4 ) ] ) ) for k in range ( n ) ]
        if max ( twists ) > max_twist:
            continue
        if atoms is not None:
            aromatic_lengths = True
            for k in range ( n ):
                a, b = ring[k], ring[( k + 1 ) % n]
                limit = _AROMATIC_BOND_MAX.get ( _pair ( atoms[a].symbol, atoms[b].symbol ) )
                if limit is None or _dist ( coords[a], coords[b] ) > limit:
                    aromatic_lengths = False
                    break
            if not aromatic_lengths:
                continue
        out.append ( ring )
    return out


def geometric_implicit_hydrogens ( atoms, neighbors, coords ):
    """ {atom_id: implicit H count} for heavy atoms with NO explicit
        hydrogen, from local geometry (see the rules above). coords:
        {atom_id: (x, y, z)}; neighbors: {atom_id: [heavy neighbour ids]}.
        Atoms of flat 5/6-membered rings are sp2 whatever their angle (a
        pentagon's ~108 deg would otherwise read as sp3); in a flat
        5-membered ring Hueckel needs exactly one lone-pair donor (O, S,
        3-connected N or an N-H): when the ring has none, the 2-connected N
        with the widest C-N-C angle is taken as the N-H. """
    import math
    out = { }
    ring_sp2 = set ( )
    pyrrole_nh = set ( )
    for ring in _planar_rings ( neighbors, coords, atoms ):
        ring_sp2.update ( ring )
        if len ( ring ) == 5:
            syms = { a: atoms[a].symbol.upper ( ) for a in ring }
            donors = [ a for a in ring if syms[a] in ( "O", "S", "SE" ) or ( syms[a] == "N" and len ( neighbors[a] ) == 3 ) ]
            if not donors:
                n2 = [ a for a in ring if syms[a] == "N" and len ( neighbors[a] ) == 2 ]
                if n2:
                    widest = max ( n2, key = lambda a: _angle ( coords[neighbors[a][0]], coords[a], coords[neighbors[a][1]] ) )
                    pyrrole_nh.add ( widest )
    for aid, nbs in neighbors.items ( ):
        sym = atoms[aid].symbol.upper ( )
        base = _BASE_VALENCE.get ( sym )
        if base is None:
            continue
        n = len ( nbs )
        if n == 0:
            out[aid] = base
            continue
        if n == 1:
            nb = nbs[0]
            d = math.sqrt ( sum ( ( coords[aid][k] - coords[nb][k] ) ** 2 for k in range ( 3 ) ) )
            rules = _TERMINAL_LENGTH_RULES.get ( ( sym, atoms[nb].symbol.upper ( ) ) )
            h = None
            if rules:
                for max_len, nh in rules:
                    if d <= max_len:
                        h = nh
                        break
            out[aid] = h if h is not None else base - 1
            continue
        angles = [ _angle ( coords[nbs[i]], coords[aid], coords[nbs[j]] )
                   for i in range ( n ) for j in range ( i + 1, n ) ]
        mean = sum ( angles ) / len ( angles )
        planar = ( sum ( angles ) > 350.0 ) if n == 3 else ( mean > 114.0 )
        short_bond = any ( _dist ( coords[aid], coords[nb] ) <= _MULTIPLE_BOND_MAX.get ( _pair ( sym, atoms[nb].symbol ), 0.0 )
                           for nb in nbs )
        if mean > 155.0 and n == 2:
            pi = 2                     # sp (linear)
        elif aid in ring_sp2:
            pi = 1                     # atom of an aromatic ring
        elif planar and short_bond:
            pi = 1                     # sp2: planar AND carrying a short (multiple) bond
        else:
            pi = 0                     # sp3 (or a planar atom whose bonds are all single)
        if sym == "N" and pi == 1 and n == 3:
            pi = 0                     # planar N with 3 neighbours: amide/aniline/pyrrole lone pair
        if sym in ( "O", "S", "SE" ) and n == 2:
            pi = 0
        if aid in pyrrole_nh:
            pi = 0                     # N-H of a pyrrole/imidazole/indole ring
        out[aid] = max ( 0, base - n - pi )
    return out


def _perceive_fragments ( atoms, all_pairs, unresolved, orders, coords = None ):
    """ Runs the Wang & Case perceiver on each connected fragment of the
        unresolved bonds. A bond from a fragment atom to an atom OUTSIDE the
        fragment (already resolved, e.g. a covalent ligand-protein link) is
        represented by a placeholder hydrogen so that atom's valence stays
        right. Fills `orders` in place; returns (n_fragments, n_atoms). """
    if not unresolved:
        return 0, 0
    from vismol.core.bond_order_perception import perceive_bond_orders

    parent = { }
    def find ( x ):
        while parent.setdefault ( x, x ) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for ( i, j ) in unresolved:
        ri, rj = find ( i ), find ( j )
        if ri != rj:
            parent[ri] = rj
    fragments = defaultdict ( list )
    for pair in unresolved:
        fragments[find ( pair[0] )].append ( pair )

    neighbors = defaultdict ( list )
    for ( i, j ) in all_pairs:
        neighbors[i].append ( j ); neighbors[j].append ( i )

    n_atoms_total = 0
    for frag_pairs in fragments.values ( ):
        frag_atoms = sorted ( { a for p in frag_pairs for a in p } )
        local = { aid: k for k, aid in enumerate ( frag_atoms ) }
        elements = [ atoms[aid].symbol for aid in frag_atoms ]
        local_pairs = [ ( local[i], local[j] ) for ( i, j ) in frag_pairs ]
        frag_pair_set = set ( frag_pairs )
        # heavy-atom-only fragment: implicit hydrogens from the geometry
        if coords is not None and not any ( atoms[nb].symbol == "H" for aid in frag_atoms for nb in neighbors[aid] ):
            heavy_nbs = { aid: [ nb for nb in neighbors[aid] if atoms[nb].symbol != "H" ] for aid in frag_atoms }
            for aid, n_h in geometric_implicit_hydrogens ( atoms, heavy_nbs, coords ).items ( ):
                for _ in range ( n_h ):
                    elements.append ( "H" )
                    local_pairs.append ( ( local[aid], len ( elements ) - 1 ) )
        for aid in frag_atoms:
            for nb in neighbors[aid]:
                key = ( min ( aid, nb ), max ( aid, nb ) )
                if key in frag_pair_set:
                    continue
                # bond leaving the fragment (or to a hydrogen): keep it as a
                # single bond to a placeholder H in the local problem
                elements.append ( "H" )
                local_pairs.append ( ( local[aid], len ( elements ) - 1 ) )
        try:
            order_map, _tps = perceive_bond_orders ( elements, local_pairs )
        except Exception:
            order_map = { }
        for ( i, j ) in frag_pairs:
            li, lj = local[i], local[j]
            orders[( i, j )] = int ( order_map.get ( ( min ( li, lj ), max ( li, lj ) ), 1 ) )
        n_atoms_total += len ( frag_atoms )
    return len ( fragments ), n_atoms_total


# ---------------------------------------------------------------------------
#   aromaticity / hybridization
# ---------------------------------------------------------------------------

def _small_rings ( adjacency, max_size = 7 ):
    """ Every simple cycle of up to max_size atoms (as frozensets of atom
        ids, plus the ordered atom list), found by a bounded DFS from each
        atom; rings are reported once. Enough for real chemistry (aromatic
        rings are 5-7 membered). """
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


def _pi_electrons ( atom_symbol, ring_bond_orders_of_atom, exo_double, n_bonds, formal_charge ):
    """ Hueckel pi-electron contribution of one ring atom, or None if it
        cannot be part of an aromatic ring (sp3 carbon, ...). """
    if 2 in ring_bond_orders_of_atom:
        return 1
    if atom_symbol == "C":
        if formal_charge == -1:
            return 2
        if formal_charge == 1:
            return 0
        if exo_double:          # e.g. C=O in a pyridone ring
            return 0
        return None
    if atom_symbol in ( "N", "P" ):
        if formal_charge == 1 and n_bonds == 4:
            return None
        return 2                # pyrrole-type N (lone pair)
    if atom_symbol in ( "O", "S", "Se" ):
        return 2
    return None


def perceive_aromatic_bonds ( atoms, orders, formal_charges = None ):
    """ Set of (i, j) bonds that belong to an aromatic ring (Hueckel 4n+2),
        judged on the Kekule orders in `orders` ({pair: order}). Single
        rings first; then fused systems (union of rings sharing a bond)
        whose individual rings failed only because of how the Kekule
        double bonds were distributed (e.g. indole, purine). """
    formal_charges = formal_charges or { }
    adjacency = defaultdict ( set )
    h_count = defaultdict ( int )
    for ( i, j ) in orders:
        if atoms[i].symbol == "H" or atoms[j].symbol == "H":
            if atoms[i].symbol != "H": h_count[i] += 1
            if atoms[j].symbol != "H": h_count[j] += 1
            continue
        adjacency[i].add ( j ); adjacency[j].add ( i )
    rings = [ r for r in _small_rings ( adjacency ) if 5 <= len ( r ) <= 7 ]

    def order ( a, b ):
        return orders.get ( ( min ( a, b ), max ( a, b ) ), 1 )

    def ring_pi_count ( ring_atoms, ring_bonds ):
        total = 0
        for a in ring_atoms:
            in_ring = [ order ( a, b ) for b in adjacency[a] if ( min ( a, b ), max ( a, b ) ) in ring_bonds ]
            exo = any ( order ( a, b ) == 2 for b in adjacency[a] if ( min ( a, b ), max ( a, b ) ) not in ring_bonds )
            n_bonds = len ( adjacency[a] ) + h_count[a]
            pi = _pi_electrons ( atoms[a].symbol, in_ring, exo, n_bonds, formal_charges.get ( a, getattr ( atoms[a], "formal_charge", 0 ) ) )
            if pi is None:
                return None
            total += pi
        return total

    def bonds_of ( ring ):
        return { ( min ( ring[k], ring[( k + 1 ) % len ( ring )] ), max ( ring[k], ring[( k + 1 ) % len ( ring )] ) ) for k in range ( len ( ring ) ) }

    aromatic = set ( )
    failed = [ ]
    for ring in rings:
        rb = bonds_of ( ring )
        pi = ring_pi_count ( set ( ring ), rb )
        if pi is not None and pi % 4 == 2:
            aromatic |= rb
        elif pi is not None:
            failed.append ( ( set ( ring ), rb ) )
    # fused systems: merge failed rings with any ring they share a bond with
    for ring_atoms, rb in failed:
        for other in rings:
            ob = bonds_of ( other )
            if rb & ob and set ( other ) != ring_atoms:
                union_atoms = ring_atoms | set ( other )
                perimeter = ( rb | ob ) - ( rb & ob )
                pi = ring_pi_count ( union_atoms, rb | ob )
                if pi is not None and pi % 4 == 2:
                    aromatic |= rb | ob
    return aromatic


def atom_hybridizations ( atoms, orders, aromatic_bonds ):
    """ {atom_id: "aromatic" | "sp" | "sp2" | "sp3" | ""} for every atom
        ("" for H / isolated atoms). Conjugated N and O with only single
        bonds but attached to an sp2/sp or aromatic centre (amide, aniline,
        guanidine, carboxylate O-, enol/phenol O) are reported as sp2, like
        most force fields do. """
    bonds_of = defaultdict ( list )
    for p, o in orders.items ( ):
        bonds_of[p[0]].append ( ( p, o ) ); bonds_of[p[1]].append ( ( p, o ) )

    def unsaturated ( a ):
        return any ( q in aromatic_bonds or o2 >= 2 for q, o2 in bonds_of[a] )

    out = { }
    for aid, atom in atoms.items ( ):
        my_bonds = bonds_of.get ( aid, [ ] )
        if atom.symbol == "H" or not my_bonds:
            out[aid] = ""
            continue
        if any ( p in aromatic_bonds for p, _o in my_bonds ):
            out[aid] = "aromatic"; continue
        max_order = max ( o for _p, o in my_bonds )
        n_double = sum ( 1 for _p, o in my_bonds if o == 2 )
        if max_order == 3 or n_double >= 2:
            out[aid] = "sp"; continue
        if max_order == 2:
            out[aid] = "sp2"; continue
        if atom.symbol in ( "N", "O" ) and any (
                unsaturated ( p[0] if p[1] == aid else p[1] ) for p, _o in my_bonds ):
            out[aid] = "sp2"; continue
        out[aid] = "sp3"
    return out


def valence_formal_charges ( atoms, orders, atom_ids ):
    """ [EN] 2026-10-02: {atom_id: formal charge} for the N and O atoms in
        `atom_ids`, read from their final valence (sum of bond orders, H
        included): O with one single bond and no H -> -1, O with 3 bonds ->
        +1, N with 4 bonds -> +1, N with 2 single bonds -> -1, else 0.
        Used for non-template residues after protonation: OpenBabel decides
        HOW MANY protons a ligand has, but which of two equivalent oxygens
        carries the charge must agree with the bond orders actually chosen
        (a carboxylate O- is the singly bonded one, never the C=O). """
    valence = defaultdict ( int )
    for ( i, j ), o in orders.items ( ):
        valence[i] += o; valence[j] += o
    out = { }
    for aid in atom_ids:
        sym = atoms[aid].symbol.upper ( )
        if sym == "O":
            out[aid] = valence[aid] - 2
        elif sym == "N":
            out[aid] = valence[aid] - 3
    return out
