#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: chemical environment codes without RDKit
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description (2026-10-14):
#      The r1 / r2 environment codes of pdynamo/bond_library.py (the keys of the
#      DFT torsion library and of the Hessian calibration, made by fraglib with
#      RDKit) computed in pure Python. Every atom invariant RDKit contributes
#      (aromaticity, ring membership, hybridization, and the formal charges and
#      bond orders after its "cleanUp" of nitro / P / halogen oxo groups) is a
#      port of the RDKit 2026.03.6 code that the library was made with:
#
#        Code/GraphMol/MolOps.cpp        sanitizeMol, cleanUp (nitrogensCleanup,
#                                        phosphorusCleanup, halogenCleanup)
#        Code/GraphMol/Aromaticity.cpp   default (RDKit) aromaticity model
#        Code/GraphMol/ConjugHybrid.cpp  setConjugation, setHybridization
#        Code/GraphMol/atomic_data.cpp   outer-shell electrons and valence lists
#
#      RDKit: Copyright (C) 2001-2026 Greg Landrum and other RDKit contributors,
#      BSD license (https://github.com/rdkit/rdkit/blob/master/license.txt).
#
#      Rings: RDKit uses its symmetrized SSSR; here the relevant cycles (cycles
#      that are not a sum of shorter ones), which are the same set for the
#      molecules of the libraries (checked against RDKit, DYFF_validation/
#      scripts/env_codes_check.py).
#
#      Input: atomic numbers, Kekule bonds ( i, j, order ) and formal charges,
#      ALL hydrogens explicit (as the library's molecules and bond_library.
#      rdkit_molecule). No implicit hydrogens, no radicals, no dative bonds.
#

import hashlib
from collections import defaultdict, deque
from itertools import combinations

ENV_CODE_VERSION = "1"            # . r1 / r2 strings identical to fraglib's version 1
ENV_CODE_VERSION_OWN = "2"        # . libraries keyed with this module: exact = exact_code ( ) (no stereo)

# . RDKit atomic_data.cpp (Release_2026_03_6): outer-shell electrons and
#   allowed valences by atomic number 0-118 (-1 = no limit known)
_N_OUTER = ( 0, 1, 2, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 )
_VALENCES = ( (-1,), (1,), (0,), (1, -1), (2,), (3,), (4,), (3,), (2,), (1,), (0,), (1, -1), (2, -1), (3,), (4,), (3, 5), (2, 4, 6), (1,), (0,), (1, -1), (2, -1), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (3,), (4,), (3, 5), (2, 4, 6), (1,), (0,), (1, -1), (2, -1), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (3,), (2, 4), (3, 5), (2, 4, 6), (1, 3, 5), (0, 2, 4, 6), (1,), (2, -1), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (2, 4), (3, 5), (2, 4, 6), (1, 3, 5), (0,), (1,), (2, -1), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,), (-1,) )

_MAX_FUSED_AROMATIC_RING_SIZE = 24          # . Aromaticity.cpp maxFusedAromaticRingSize
_MAX_FUSED_RINGS              = 6           # . aromaticityHelper: applyHuckelToFused ( ..., 6, ... )

# . electron donor types (Aromaticity.cpp ElectronDonorType)
_VACANT, _ONE, _TWO, _ONE_OR_TWO, _ANY, _NONE = range ( 6 )


def _default_valence ( z ):
    return _VALENCES[z][0] if 0 <= z < len ( _VALENCES ) else -1

def _n_outer ( z ):
    return _N_OUTER[z] if 0 <= z < len ( _N_OUTER ) else 0

def _more_electronegative ( z1, z2 ):
    """ PeriodicTable::moreElectroNegative: more outer electrons, then the lighter. """
    n1, n2 = _n_outer ( z1 ), _n_outer ( z2 )
    return n1 > n2 or ( n1 == n2 and z1 < z2 )


#=============================================================================
#  Molecule
#=============================================================================
class Molecule:
    """ Atoms ( atomic numbers, formal charges ) and bonds ( i, j, order ), the
        orders Kekule (1, 2, 3); after perceive ( ): .aromatic_atoms, bond types
        ( 1, 2, 3 or "a" ), .ring_atoms, .hybridization, .conjugated. """

    def __init__ ( self, atomic_numbers, bonds, formal_charges = None ):
        self.z      = [ int ( z ) for z in atomic_numbers ]
        n           = len ( self.z )
        self.charge = [ int ( q ) for q in formal_charges ] if formal_charges is not None else [ 0 ] * n
        self.bonds  = [ ]                         # . [ i, j, type ]: type 1, 2, 3 or "a"
        self.nb     = [ [ ] for _ in range ( n ) ]   # . [ ( neighbour, bond index ) ]
        self.bond_of = { }
        for i, j, order in bonds:
            i, j = int ( i ), int ( j )
            b = len ( self.bonds )
            self.bonds.append ( [ i, j, int ( order ) ] )
            self.nb[i].append ( ( j, b ) ); self.nb[j].append ( ( i, b ) )
            self.bond_of[( i, j )] = self.bond_of[( j, i )] = b
        self.perceive ( )

    # . valence helpers ------------------------------------------------------
    def degree ( self, a ):
        return len ( self.nb[a] )

    def _contrib ( self, b ):
        t = self.bonds[b][2]
        return 1.5 if t == "a" else float ( t )

    def explicit_valence ( self, a ):
        """ Atom::calcExplicitValence on the Kekule molecule (no implicit H). """
        return int ( round ( sum ( self._contrib ( b ) for _, b in self.nb[a] ) ) )

    # . perception ----------------------------------------------------------
    def perceive ( self ):
        self._clean_up ( )
        # . valences cached as RDKit's property cache: computed on the Kekule form,
        #   before aromaticity, and used by conjugation and hybridization
        self.valence = [ self.explicit_valence ( a ) for a in range ( len ( self.z ) ) ]
        self.radicals = [ self._radicals ( a ) for a in range ( len ( self.z ) ) ]
        self.rings = relevant_cycles ( len ( self.z ), [ ( i, j ) for i, j, _ in self.bonds ] )
        self.ring_bonds = set ( )
        for ring in self.rings:
            self.ring_bonds.update ( self._ring_bonds ( ring ) )
        self.ring_atoms = { a for ring in self.rings for a in ring }
        self.aromatic_atoms = set ( )
        self._aromaticity ( )
        self._conjugation ( )
        self._hybridization ( )

    def _ring_bonds ( self, ring ):
        return [ self.bond_of[( ring[k], ring[( k + 1 ) % len ( ring )] )] for k in range ( len ( ring ) ) ]

    #=========================================================================
    #  MolOps.cpp cleanUp
    #=========================================================================
    def _clean_up ( self ):
        z, q = self.z, self.charge
        # . nitrogensCleanup: neutral 5-valent N=O -> [N+]-[O-]; then N#N -> [N+]=[N-]
        five = [ a for a in range ( len ( z ) ) if z[a] == 7 and q[a] == 0 and self.explicit_valence ( a ) == 5 ]
        for a in five:
            for o, b in self.nb[a]:
                if z[o] == 8 and q[o] == 0 and self.bonds[b][2] == 2:
                    self.bonds[b][2] = 1; q[a] = 1; q[o] = -1
                    break
        for a in five:
            for o, b in self.nb[a]:
                if z[o] == 7 and q[o] == 0 and self.bonds[b][2] == 3:
                    self.bonds[b][2] = 2; q[a] = 1; q[o] = -1
                    break
        for a in range ( len ( z ) ):
            if z[a] == 15:
                # . phosphorusCleanup: C=P(=O)X -> C=[P+]([O-])X
                if q[a] == 0 and self.explicit_valence ( a ) == 5 and self.degree ( a ) == 3:
                    double_o = None; double_cn = False
                    for o, b in self.nb[a]:
                        if z[o] == 8 and q[o] == 0 and self.bonds[b][2] == 2:
                            double_o = ( o, b )
                        elif z[o] in ( 6, 7 ) and self.degree ( o ) >= 2 and self.bonds[b][2] == 2:
                            double_cn = True
                    if double_cn and double_o is not None:
                        q[double_o[0]] = -1; self.bonds[double_o[1]][2] = 1; q[a] = 1
            elif z[a] in ( 17, 35, 53 ):
                # . halogenCleanup: X(=O)n(O) -> [X+n]([O-])n(O)
                ev = self.explicit_valence ( a )
                if q[a] == 0 and ev in ( 7, 5, 3 ) and all ( z[o] == 8 for o, _ in self.nb[a] ):
                    charge = 0
                    for o, b in self.nb[a]:
                        if self.bonds[b][2] == 2:
                            self.bonds[b][2] = 1; charge += 1; q[o] = -1
                    q[a] = charge

    #=========================================================================
    #  MolOps.cpp assignRadicals (every atom has no implicit H)
    #=========================================================================
    def _radicals ( self, a ):
        z = self.z[a]
        if z == 0:
            return 0
        valens = _VALENCES[z] if z < len ( _VALENCES ) else ( -1, )
        chg = self.charge[a]; nouter = _n_outer ( z )
        if len ( valens ) != 1 or valens[0] != -1:
            total = self.valence[a]
            base = 2 if z in ( 1, 2 ) else 8
            n = base - nouter - total + chg
            if n < 0:
                n = 0
                if len ( valens ) > 1:
                    for v in valens:
                        if v - total + chg >= 0:
                            n = v - total + chg
                            break
            n2 = nouter - total - chg
            if n2 >= 0:
                n = min ( n, n2 )
            return n
        if self.degree ( a ) > 0:
            return 0
        return max ( nouter - chg, 0 ) % 2

    #=========================================================================
    #  Aromaticity.cpp (AROMATICITY_DEFAULT = aromaticityHelper ( 0, 0, fused ))
    #=========================================================================
    def _count_atom_elec ( self, a ):
        dv = _default_valence ( self.z[a] )
        if dv <= 1:
            return -1
        degree = self.degree ( a )
        if degree > 3:
            return -1
        nlp = max ( _n_outer ( self.z[a] ) - dv - self.charge[a], 0 )
        res = ( dv - degree ) + nlp - self.radicals[a]
        if res > 1 and self.valence[a] - degree > 1:
            res = 1
        return res

    def _incident_noncyclic_multiple ( self, a ):
        for o, b in self.nb[a]:
            if b not in self.ring_bonds and self._contrib ( b ) >= 2.0:
                return o
        return None

    def _incident_cyclic_multiple ( self, a ):
        return any ( b in self.ring_bonds and self._contrib ( b ) >= 2.0 for _, b in self.nb[a] )

    def _incident_multiple ( self, a ):
        return self.valence[a] != self.degree ( a )

    def _donor_type ( self, a ):
        z = self.z[a]
        if z == 0:
            return _ONE if self._incident_cyclic_multiple ( a ) else _ANY
        nelec = self._count_atom_elec ( a )
        if nelec < 0:
            return _NONE
        if nelec == 0:
            if self._incident_noncyclic_multiple ( a ) is not None:
                return _VACANT
            if self._incident_cyclic_multiple ( a ):
                return _ONE
            return _NONE
        if nelec == 1:
            who = self._incident_noncyclic_multiple ( a )
            if who is not None:
                return _VACANT if _more_electronegative ( self.z[who], z ) else _ONE
            if self._incident_multiple ( a ):
                return _ONE
            if self.charge[a] == 1:
                return _VACANT
            return _NONE
        who = self._incident_noncyclic_multiple ( a )
        if who is not None and _more_electronegative ( self.z[who], z ):
            nelec -= 1
        return _ONE if nelec % 2 == 1 else _TWO

    def _candidate ( self, a, edon ):
        z = self.z[a]
        if z > 18 and z not in ( 34, 52 ):
            return False
        if edon not in ( _VACANT, _ONE, _TWO, _ONE_OR_TWO, _ANY ):
            return False
        dv = _default_valence ( z )
        if dv > 0 and self.valence[a] > _default_valence ( z - self.charge[a] ):
            return False
        if self.radicals[a] and ( z != 6 or self.charge[a] ):
            return False
        if self.valence[a] - self.degree ( a ) > 1:
            n_mult = sum ( 1 for _, b in self.nb[a] if self.bonds[b][2] in ( 2, 3 ) )
            if n_mult > 1:
                return False
        return True

    def _aromaticity ( self ):
        edon = { }; cand = { }
        candidate_rings = [ ]
        for ring in self.rings:
            size = len ( ring )
            all_aromatic = True
            for a in ring:
                if a in edon:
                    if not cand[a]: all_aromatic = False
                    continue
                donor = self._donor_type ( a )
                # . macrocycle fix (2026.03): an ether O / S of a ring of 9+ atoms
                if donor == _TWO and size >= 9 and self.z[a] in ( 8, 16 ) and self.degree ( a ) == 2 and self.charge[a] == 0 and \
                   not any ( self.bonds[b][2] in ( 2, 3 ) for _, b in self.nb[a] ):
                    donor = _NONE
                edon[a] = donor
                cand[a] = self._candidate ( a, donor )
                if not cand[a]: all_aromatic = False
            if all_aromatic and not all ( self.z[a] == 0 for a in ring ):
                candidate_rings.append ( ring )
        bond_rings = [ self._ring_bonds ( r ) for r in candidate_rings ]
        # . neighbour map: rings of <= 24 atoms sharing exactly one bond
        neigh = { i: [ ] for i in range ( len ( candidate_rings ) ) }
        for i, j in combinations ( range ( len ( candidate_rings ) ), 2 ):
            if len ( bond_rings[i] ) > _MAX_FUSED_AROMATIC_RING_SIZE or len ( bond_rings[j] ) > _MAX_FUSED_AROMATIC_RING_SIZE:
                continue
            if len ( set ( bond_rings[i] ) & set ( bond_rings[j] ) ) == 1:
                neigh[i].append ( j ); neigh[j].append ( i )
        done = set ( )
        for start in range ( len ( candidate_rings ) ):
            if start in done: continue
            system, stack = [ ], [ start ]
            while stack:                                   # . pickFusedRings
                r = stack.pop ( )
                if r in done: continue
                done.add ( r ); system.append ( r )
                stack.extend ( x for x in neigh[r] if x not in done )
            self._huckel_fused ( candidate_rings, bond_rings, system, edon, neigh )

    def _huckel_fused ( self, rings, bond_rings, fused, edon, neigh ):
        n_ring_bonds = len ( { b for r in fused for b in bond_rings[r] } )
        done_bonds = set ( )
        for size in range ( 1, min ( len ( fused ), _MAX_FUSED_RINGS ) + 1 ):
            if len ( done_bonds ) >= n_ring_bonds:
                break
            if size == 3 and len ( fused ) > 300:
                break
            for comb in combinations ( fused, size ):
                if size > 1 and not _connected ( comb, neigh ):
                    continue
                count = defaultdict ( int )
                for r in comb:
                    for a in rings[r]:
                        count[a] += 1
                atoms = [ a for a, c in count.items ( ) if c in ( 1, 2 ) ]
                if not _huckel ( [ edon[a] for a in atoms ] ):
                    continue
                bond_count = defaultdict ( int )
                for r in comb:
                    for b in bond_rings[r]:
                        bond_count[b] += 1
                for b, c in bond_count.items ( ):
                    if c != 1: continue
                    bond = self.bonds[b]
                    if bond[2] in ( 1, 2, "a" ):
                        if bond[2] != "a":
                            bond[2] = "a"
                        self.aromatic_atoms.update ( bond[:2] )
                    done_bonds.add ( b )

    #=========================================================================
    #  ConjugHybrid.cpp
    #=========================================================================
    def _conjugation_candidate ( self, a ):
        z = self.z[a]
        vals = _VALENCES[z] if z < len ( _VALENCES ) else ( -1, )
        if not self.charge[a] and vals[0] >= 0 and self.valence[a] > vals[0]:
            return False
        nouter = _n_outer ( z )
        return ( z <= 10 or nouter not in ( 5, 6 ) or ( nouter == 6 and self.degree ( a ) < 2 ) ) and self._count_atom_elec ( a ) > 0

    def _conjugation ( self ):
        self.conjugated = { b for b, bond in enumerate ( self.bonds ) if bond[2] == "a" }
        cand = [ self._conjugation_candidate ( a ) for a in range ( len ( self.z ) ) ]
        subst = [ self.degree ( a ) if cand[a] else 0 for a in range ( len ( self.z ) ) ]
        for a in range ( len ( self.z ) ):
            if not cand[a] or subst[a] < 2 or subst[a] > 3:
                continue
            for o1, b1 in self.nb[a]:
                if self._contrib ( b1 ) < 1.5 or not cand[o1]:
                    continue
                for o2, b2 in self.nb[a]:
                    if b1 == b2 or subst[o2] > 3:
                        continue
                    if cand[o2]:
                        self.conjugated.add ( b1 ); self.conjugated.add ( b2 )

    def _bonds_plus_lone_pairs ( self, a ):
        deg = self.degree ( a )
        if self.z[a] <= 1:
            return deg
        nouter = _n_outer ( self.z[a] ); tv = self.valence[a]; chg = self.charge[a]
        free = nouter - ( tv + chg )
        if tv + nouter - chg < 8:
            rad = self.radicals[a]
            return deg + int ( ( free - rad ) / 2 ) + rad    # . C++ integer division
        return deg + int ( free / 2 )

    def _hybridization ( self ):
        self.hybridization = [ ]
        stereo = None
        for a, z in enumerate ( self.z ):
            if z == 0:
                self.hybridization.append ( "UNSPECIFIED" ); continue
            # . a tetrahedral stereo centre with 4 neighbours is SP3 (setHybridization
            #   reads the chiral tag first). RDKit's tags come from the SMILES (library)
            #   or the 3D structure (lookup: bond_library.rdkit_molecule); here a
            #   4-connected atom whose neighbours all differ is taken as tagged. It
            #   only matters where the rule below gives something else: S / Se / Te
            #   with 4 neighbours and a lone pair.
            if self.degree ( a ) == 4:
                if stereo is None:
                    stereo = _stereo_candidates ( self )
                if a in stereo:
                    self.hybridization.append ( "SP3" ); continue
            norbs = self._bonds_plus_lone_pairs ( a ) if z < 89 else self.degree ( a )
            if norbs in ( 0, 1 ):
                h = "S"
            elif norbs == 2:
                h = "SP"
            elif norbs == 3:
                h = "SP2"
            elif norbs == 4:
                h = "SP3" if ( self.degree ( a ) > 3 or not any ( b in self.conjugated for _, b in self.nb[a] ) ) else "SP2"
            elif norbs == 5:
                h = "SP3D"
            elif norbs == 6:
                h = "SP3D2"
            else:
                h = "UNSPECIFIED"
            self.hybridization.append ( h )

    #=========================================================================
    #  Invariants and codes (same strings as bond_library / fraglib)
    #=========================================================================
    def bond_symbol ( self, b ):
        t = self.bonds[b][2]
        return { 1: "1", 2: "2", 3: "3", "a": "a" }.get ( t, "?" )

    def atom_invariant ( self, a ):
        # . GetTotalNumHs ( ) is 0: every hydrogen is an explicit atom
        return "{}|{}|{}|{}|{}|{}|{}".format ( self.z[a], int ( a in self.aromatic_atoms ), self.degree ( a ), 0,
                                               self.charge[a], int ( a in self.ring_atoms ), self.hybridization[a] )


def _symmetry_classes ( mol ):
    """ Morgan-style refinement of ( Z, charge, degree, bond types ) to stability. """
    cls = [ hash ( ( z, q, mol.degree ( a ) ) ) for a, ( z, q ) in enumerate ( zip ( mol.z, mol.charge ) ) ]
    n_classes = len ( set ( cls ) )
    while True:
        new = [ hash ( ( cls[a], tuple ( sorted ( ( mol.bond_symbol ( b ), cls[o] ) for o, b in mol.nb[a] ) ) ) ) for a in range ( len ( cls ) ) ]
        n_new = len ( set ( new ) )
        cls = new
        if n_new == n_classes:
            return cls
        n_classes = n_new


def _stereo_candidates ( mol ):
    """ 4-connected atoms whose four neighbours are in four symmetry classes. """
    cls = _symmetry_classes ( mol )
    return { a for a in range ( len ( mol.z ) ) if mol.degree ( a ) == 4 and len ( { cls[o] for o, _ in mol.nb[a] } ) == 4 }


def _connected ( comb, neigh ):
    """ RingUtils::checkFused: the rings of comb form one connected system. """
    members = set ( comb ); seen = { comb[0] }; stack = [ comb[0] ]
    while stack:
        r = stack.pop ( )
        for x in neigh[r]:
            if x in members and x not in seen:
                seen.add ( x ); stack.append ( x )
    return len ( seen ) == len ( members )


def _huckel ( donors ):
    """ applyHuckel: 4n+2 within the electron range of the ring atoms. """
    low = up = 0; n_any = 0
    for d in donors:
        if d == _ANY:
            n_any += 1
            if n_any > 1: return False
        lo, hi = { _ANY: ( 1, 2 ), _ONE_OR_TWO: ( 1, 2 ), _ONE: ( 1, 1 ), _TWO: ( 2, 2 ) }.get ( d, ( 0, 0 ) )
        low += lo; up += hi
    if up >= 6:
        return any ( ( e - 2 ) % 4 == 0 for e in range ( low, up + 1 ) )
    return up == 2


#=============================================================================
#  Rings: relevant cycles
#=============================================================================
def relevant_cycles ( n_atoms, edges ):
    """ Cycles (atom lists in ring order) that are not a GF(2) sum of shorter
        cycles, from the candidates "shortest cycles through each ring bond"
        (all shortest paths). Sorted by size, then atoms. """
    adj = defaultdict ( set )
    for i, j in edges:
        adj[i].add ( j ); adj[j].add ( i )
    # . ring bonds: remove degree-1 atoms repeatedly, then bridges
    ring_edges = _cyclic_edges ( adj )
    if not ring_edges:
        return [ ]
    radj = defaultdict ( set )
    for i, j in ring_edges:
        radj[i].add ( j ); radj[j].add ( i )
    edge_index = { }
    for k, ( i, j ) in enumerate ( sorted ( ring_edges ) ):
        edge_index[( i, j )] = edge_index[( j, i )] = k
    candidates = { }
    for x, y in sorted ( ring_edges ):
        for path in _all_shortest_paths ( radj, y, x, ( x, y ) ):
            cyc = tuple ( path )
            key = frozenset ( edge_index[( cyc[k], cyc[( k + 1 ) % len ( cyc )] )] for k in range ( len ( cyc ) ) )
            if len ( key ) == len ( cyc ) and key not in candidates:
                candidates[key] = cyc
    ordered = sorted ( candidates.items ( ), key = lambda kv: ( len ( kv[1] ), sorted ( kv[1] ) ) )
    out = [ ]
    basis_shorter = { }                 # . pivot -> vector, from cycles strictly shorter
    basis_same    = { }
    current = None
    for key, cyc in ordered:
        if len ( cyc ) != current:
            basis_shorter.update ( basis_same ); basis_same = { }
            current = len ( cyc )
        vec = 0
        for e in key: vec |= 1 << e
        # . independent of strictly shorter cycles -> relevant
        if _reduce ( vec, basis_shorter ):
            out.append ( list ( cyc ) )
            red = _reduce ( vec, { **basis_shorter, **basis_same } )
            if red:
                basis_same[red.bit_length ( ) - 1] = red
    return out


def _reduce ( vec, basis ):
    while vec:
        top = vec.bit_length ( ) - 1
        if top not in basis:
            return vec
        vec ^= basis[top]
    return 0


def _cyclic_edges ( adj ):
    """ Edges on at least one cycle (non-bridges), Tarjan. """
    index = { }; low = { }; counter = [ 0 ]; bridges = set ( )
    for root in list ( adj ):
        if root in index: continue
        stack = [ ( root, None, iter ( sorted ( adj[root] ) ) ) ]
        index[root] = low[root] = counter[0]; counter[0] += 1
        while stack:
            v, parent, it = stack[-1]
            advanced = False
            for w in it:
                if w == parent: continue
                if w not in index:
                    index[w] = low[w] = counter[0]; counter[0] += 1
                    stack.append ( ( w, v, iter ( sorted ( adj[w] ) ) ) )
                    advanced = True
                    break
                low[v] = min ( low[v], index[w] )
            if not advanced:
                stack.pop ( )
                if parent is not None:
                    low[parent] = min ( low[parent], low[v] )
                    if low[v] > index[parent]:
                        bridges.add ( frozenset ( ( parent, v ) ) )
    return { ( min ( i, j ), max ( i, j ) ) for i in adj for j in adj[i] if frozenset ( ( i, j ) ) not in bridges }


def _all_shortest_paths ( adj, source, target, skip ):
    """ Every shortest path source -> target not using the edge skip. """
    dist = { source: 0 }; parents = defaultdict ( list ); queue = deque ( [ source ] )
    blocked = { skip, skip[::-1] }
    while queue:
        v = queue.popleft ( )
        if v == target: continue
        for w in adj[v]:
            if ( v, w ) in blocked: continue
            if w not in dist:
                dist[w] = dist[v] + 1; parents[w].append ( v ); queue.append ( w )
            elif dist[w] == dist[v] + 1:
                parents[w].append ( v )
    if target not in dist:
        return [ ]
    paths = [ ]
    def walk ( v, acc ):
        if len ( paths ) > 64: return
        if v == source:
            paths.append ( [ source ] + acc[::-1] ); return
        for p in parents[v]:
            walk ( p, acc + [ v ] )
    walk ( target, [ ] )
    return paths


#=============================================================================
#  Codes
#=============================================================================
def _sha ( text, n = 16 ):
    return hashlib.sha1 ( text.encode ( ) ).hexdigest ( )[:n]


def atom_identifiers ( mol, radius ):
    ids = [ _sha ( mol.atom_invariant ( a ), 12 ) for a in range ( len ( mol.z ) ) ]
    for _ in range ( radius ):
        new = [ ]
        for a in range ( len ( mol.z ) ):
            neighbours = sorted ( "{}{}".format ( mol.bond_symbol ( b ), ids[o] ) for o, b in mol.nb[a] )
            new.append ( _sha ( ids[a] + "(" + ",".join ( neighbours ) + ")", 12 ) )
        ids = new
    return ids


def term_codes_r ( mol, atoms, kind = "bond", identifiers = None ):
    """ { "r1", "r2" } of a chain term (bond_library.term_codes without "exact"). """
    atoms = [ int ( a ) for a in atoms ]
    orders = [ ]
    for a, b in zip ( atoms, atoms[1:] ):
        k = mol.bond_of.get ( ( a, b ) )
        orders.append ( mol.bond_symbol ( k ) if k is not None else "-" )
    codes = { }
    for radius in ( 1, 2 ):
        ids = identifiers[radius] if identifiers else atom_identifiers ( mol, radius )
        forward = [ ids[a] for a in atoms ]
        key = min ( forward, forward[::-1] )
        orders_k = orders[::-1] if key != forward else orders
        codes["r{}".format ( radius )] = "{}|r{}|{}".format ( kind, radius, _sha ( "{}|{}|{}".format ( kind, ",".join ( key ), "".join ( orders_k ) ) ) )
    return codes


#=============================================================================
#  "exact" level without RDKit (2026-10-14)
#=============================================================================
#  RDKit's exact code is sha1( canonical isomeric SMILES | kind | symmetry classes
#  of the term atoms | charge | multiplicity ). Here the canonical SMILES is replaced
#  by a canonical graph identifier (iterated refinement of the atom invariants and
#  bond symbols to a stable partition, then the sorted atom labels and labelled
#  edges), and the symmetry classes by the refined atom labels. Same equivalence
#  of terms, different strings: the library keys must be recomputed with it.

def _refined_labels ( mol ):
    labels = [ _sha ( mol.atom_invariant ( a ), 12 ) for a in range ( len ( mol.z ) ) ]
    n = len ( set ( labels ) )
    while True:
        new = [ _sha ( labels[a] + "(" + ",".join ( sorted ( mol.bond_symbol ( b ) + labels[o] for o, b in mol.nb[a] ) ) + ")", 12 )
                for a in range ( len ( labels ) ) ]
        m = len ( set ( new ) )
        labels = new
        if m == n:
            return labels
        n = m


def graph_identifier ( mol, labels = None ):
    labels = labels or _refined_labels ( mol )
    edges = sorted ( "{}{}{}".format ( *sorted ( ( labels[i], labels[j] ) ), mol.bond_symbol ( b ) )
                     for b, ( i, j, _ ) in enumerate ( mol.bonds ) )
    return _sha ( "|".join ( sorted ( labels ) ) + "#" + "|".join ( edges ), 24 )


def exact_code ( mol, atoms, kind = "bond", charge = 0, multiplicity = 1, labels = None, graph_id = None ):
    labels = labels or _refined_labels ( mol )
    graph_id = graph_id or graph_identifier ( mol, labels )
    classes = [ labels[int ( a )] for a in atoms ]
    classes = min ( classes, classes[::-1] )
    return _sha ( "{}|{}|{}|{}|{}".format ( graph_id, kind, "-".join ( classes ), charge, multiplicity ) )
