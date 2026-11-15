#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: force-field-agnostic ligand chemistry perception
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      Split out of pdynamo/opls/ligand.py (2026-11): that file mixes this
#      generic part (hydrogens at a given pH, bond orders, Kekulization,
#      formal charges, topological equivalence classes -- none of it tied to
#      any force field) with the OPLS-AA-specific atom-class/LJ/charge
#      pipeline built on top of it (assign_classes, class_lennard_jones,
#      cm5_charges, finalize_types, parametrize_residue, ... -- still in
#      pdynamo/opls/ligand.py, which now imports LigandParameters and
#      ligand_chemistry() from here). DYFF's own preparation
#      (dyff_prep.py, residue_library.py, torsion_fit.py) only ever needed
#      this generic part, never the OPLS-specific one.
#
#      Output: a LigandParameters object with atoms, bonds, aromatic bonds
#      and formal charges -- everything a force-field-specific layer (OPLS
#      classes, or DYFF's own typing) needs to take it from there. Monoatomic
#      metal ions get their ionic formal charge directly (no QC).
#
#      find_xtb/find_obabel: locating the xtb / Open Babel executables is
#      also force-field-agnostic (DYFF's own esp_charges.py and the OPLS
#      layer both need an xtb binary).
#

import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict, OrderedDict

import numpy as np

ELEMENT_Z = { "H": 1, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Na": 11, "Mg": 12, "Si": 14, "P": 15, "S": 16,
              "Cl": 17, "K": 19, "Ca": 20, "Zn": 30, "Br": 35, "I": 53, "Li": 3, "Rb": 37, "Sr": 38, "Cs": 55, "Ba": 56 }
Z_ELEMENT = { z: s for s, z in ELEMENT_Z.items ( ) }

# . Monoatomic ions with a known formal charge (values match the OPLS-AA 2008
#   ion parameters, the only place this project has them tabulated): element
#   -> ( formal charge, OPLS-AA Tinker type -- the latter is unused here,
#   kept only because callers in pdynamo/opls/ligand.py still want it).
OPLS_IONS = { "Li": ( 1, 348 ), "Na": ( 1, 349 ), "K": ( 1, 350 ), "Rb": ( 1, 351 ), "Cs": ( 1, 352 ),
              "Mg": ( 2, 353 ), "Ca": ( 2, 354 ), "Sr": ( 2, 355 ), "Ba": ( 2, 356 ),
              "F": ( -1, 343 ), "Cl": ( -1, 344 ), "Br": ( -1, 345 ), "I": ( -1, 346 ) }


def find_xtb ( ):
    """ xtb executable: $XTBHOME/bin/xtb, PATH, or ~/programs/xtb-*/bin/xtb. """
    candidates = [ ]
    if os.getenv ( "XTBHOME" ): candidates.append ( os.path.join ( os.getenv ( "XTBHOME" ), "bin", "xtb" ) )
    found = shutil.which ( "xtb" )
    if found: candidates.append ( found )
    import glob
    candidates += sorted ( glob.glob ( os.path.expanduser ( "~/programs/xtb-*/bin/xtb" ) ), reverse = True )
    for path in candidates:
        if os.path.isfile ( path ) and os.access ( path, os.X_OK ): return path
    return None

def find_obabel ( ):
    return shutil.which ( "obabel" )


#=============================================================================
#  Data
#=============================================================================
class LigandParameters:
    """ Everything a force-field-specific layer (pdynamo/opls/ligand.py,
        parameters.py, ...) needs for one residue name. """
    def __init__ ( self, residue_name ):
        self.residue_name   = residue_name.upper ( )[:3]
        self.atoms          = [ ]      # [ ( name, Z, xyz ) ]
        self.bonds          = [ ]      # [ ( i, j, order ) ]
        self.aromatic_bonds = set ( )  # { ( i, j ) }
        self.formal_charges = [ ]      # per atom
        self.classes        = [ ]      # force-field-specific class per atom (e.g. OPLS)
        self.labels         = [ ]      # pDynamo type label per atom
        self.charges        = [ ]      # per atom
        self.lj             = { }      # label -> ( epsilon, sigma )
        self.equivalence    = [ ]      # equivalence class per atom
        self.charge_method  = None
        self.charge_override = None    # total charge set by the user (None = sum of the formal charges)
        self.multiplicity   = 1        # spin multiplicity of the QC charge calculation
        self.notes          = [ ]
        self.warnings       = [ ]

    @property
    def formal_charge_sum ( self ):
        return int ( sum ( self.formal_charges ) )

    @property
    def total_charge ( self ):
        """ The molecule's total charge: the user's value when given,
            otherwise the sum of the perceived formal charges. """
        return int ( self.charge_override ) if self.charge_override is not None else self.formal_charge_sum

    @property
    def n_electrons ( self ):
        return int ( sum ( z for ( _, z, _ ) in self.atoms ) ) - self.total_charge

    # . Files for pDynamo -----------------------------------------------------
    def pdb_component_mapping ( self ):
        bond_type = { 1: "Single", 2: "Double", 3: "Triple" }
        return { "Atom Fields": [ "Atomic Number", "Formal Charge", "Label" ],
                 "Atoms": [ [ z, int ( q ), name ] for ( name, z, _ ), q in zip ( self.atoms, self.formal_charges ) ],
                 "Bond Fields": [ "Atom 1", "Atom 2", "Type" ],
                 "Bonds": [ [ self.atoms[i][0], self.atoms[j][0], bond_type.get ( o, "Single" ) ] for ( i, j, o ) in self.bonds ] or None,
                 "Component Class": "Non-Polymer", "Formal Charge": self.total_charge, "Is Heteroatom": True,
                 "Label": self.residue_name, "Name": "EasyHybrid OPLS ligand", "PDB Class": "Hetas" }

    def write_pdb_component ( self, library_root ):
        import yaml
        folder = os.path.join ( library_root, "components" )
        os.makedirs ( folder, exist_ok = True )
        mapping = self.pdb_component_mapping ( )
        if mapping["Bonds"] is None:
            for key in ( "Bond Fields", "Bonds" ): mapping.pop ( key )
        path = os.path.join ( folder, self.residue_name.lower ( ) + ".yaml" )
        with open ( path, "w" ) as handle:
            handle.write ( "# . !PDBComponent\n" )
            yaml.safe_dump ( mapping, handle, default_flow_style = None, sort_keys = False, width = 200 )
        return path

    def sequence_component_mapping ( self ):
        return { "Atom Fields": [ "Label", "Type", "Charge" ],
                 "Atoms": [ [ name, label, round ( float ( q ), 6 ) ] for ( name, _, _ ), label, q in zip ( self.atoms, self.labels, self.charges ) ],
                 "Label": self.residue_name }

    def atom_type_rows ( self ):
        rows, seen = [ ], set ( )
        for ( name, z, _ ), label, q in zip ( self.atoms, self.labels, self.charges ):
            if label in seen: continue
            seen.add ( label )
            rows.append ( [ label, z, round ( float ( q ), 6 ), None, "{} {} (EasyHybrid ligand)".format ( self.residue_name, name ) ] )
        return rows

    def summary ( self ):
        return { "residue": self.residue_name, "atoms": len ( self.atoms ), "total_charge": self.total_charge,
                 "formal_charge_sum": self.formal_charge_sum, "multiplicity": self.multiplicity,
                 "charge_method": self.charge_method, "classes": dict ( Counter ( self.classes ) ),
                 "charge_sum": round ( float ( sum ( self.charges ) ), 6 ), "notes": list ( self.notes ),
                 "warnings": list ( self.warnings ) }


#=============================================================================
#  Graph helpers
#=============================================================================
def _neighbors ( n, bonds ):
    adj = defaultdict ( list )
    for ( i, j, o ) in bonds:
        adj[i].append ( ( j, o ) ); adj[j].append ( ( i, o ) )
    return adj

def _rings ( n, adj, max_size = 7 ):
    """ Smallest rings through each atom (simple BFS), as sets of atom indices. """
    rings = set ( )
    for start in range ( n ):
        for ( first, _ ) in adj[start]:
            # . shortest path from first back to start not using the start-first bond
            parent = { first: None }
            queue  = [ first ]
            found  = None
            while queue and found is None:
                nxt = [ ]
                for u in queue:
                    for ( v, _ ) in adj[u]:
                        if u == first and v == start: continue
                        if v == start:
                            found = u; break
                        if v not in parent:
                            parent[v] = u; nxt.append ( v )
                    if found is not None: break
                queue = nxt
                if len ( parent ) > 200: break
            if found is not None:
                path = [ start ]
                u = found
                while u is not None:
                    path.append ( u ); u = parent[u]
                if len ( path ) <= max_size: rings.add ( frozenset ( path ) )
    return [ set ( r ) for r in rings ]

def equivalence_classes ( elements, adj, n_rounds = 8 ):
    """ Topologically equivalent atoms (iterative Morgan-like refinement on
        element + neighbour labels). """
    labels = [ elements[i] for i in range ( len ( elements ) ) ]
    for _ in range ( n_rounds ):
        new = [ ]
        for i in range ( len ( elements ) ):
            new.append ( labels[i] + "(" + ",".join ( sorted ( labels[j] + str ( o ) for ( j, o ) in adj[i] ) ) + ")" )
        mapping = { }
        compact = [ mapping.setdefault ( lab, len ( mapping ) ) for lab in new ]
        if len ( set ( compact ) ) == len ( set ( labels ) ):
            labels = [ str ( c ) for c in compact ]
            break
        labels = [ str ( c ) for c in compact ]
    mapping = { }
    return [ mapping.setdefault ( lab, len ( mapping ) ) for lab in labels ]


#=============================================================================
#  Chemistry from Open Babel
#=============================================================================
def _write_pdb ( atoms, path, residue_name ):
    with open ( path, "w" ) as handle:
        for k, ( name, z, xyz ) in enumerate ( atoms ):
            symbol = Z_ELEMENT.get ( z, "X" )
            pdb_name = name if ( len ( name ) >= 4 or len ( symbol ) == 2 ) else " " + name
            handle.write ( "HETATM{:5d} {:<4s} {:>3s} L   1    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00          {:>2s}\n".format (
                    k + 1, pdb_name[:4], residue_name[:3], xyz[0], xyz[1], xyz[2], symbol.upper ( ) ) )
        handle.write ( "END\n" )

def _parse_mol2 ( path ):
    atoms, bonds, section = [ ], [ ], None
    with open ( path ) as handle:
        for line in handle:
            if line.startswith ( "@<TRIPOS>" ):
                section = line.strip ( )[9:]; continue
            tokens = line.split ( )
            if not tokens: continue
            if section == "ATOM" and len ( tokens ) >= 6:
                element = tokens[5].split ( "." )[0]
                element = element[0].upper ( ) + element[1:].lower ( )
                atoms.append ( ( tokens[1], element, np.array ( [ float ( tokens[2] ), float ( tokens[3] ), float ( tokens[4] ) ] ), tokens[5] ) )
            elif section == "BOND" and len ( tokens ) >= 4:
                bonds.append ( ( int ( tokens[1] ) - 1, int ( tokens[2] ) - 1, tokens[3] ) )
    return atoms, bonds

def _kekulize ( elements, bonds ):
    """ ( i, j, order ) with "ar"/"am" resolved: amide = 1, aromatic by the
        project's own Wang & Case perceiver (vismol, optional -- falls back
        to alternating guesses if it is not importable), keeping the
        explicit orders. """
    explicit = { }
    aromatic = set ( )
    for ( i, j, code ) in bonds:
        key = ( min ( i, j ), max ( i, j ) )
        if code in ( "1", "2", "3" ): explicit[key] = int ( code )
        elif code == "am":            explicit[key] = 1
        else:                         aromatic.add ( key ); explicit[key] = 1
    if aromatic:
        try:
            from vismol.core.bond_order_perception import perceive_bond_orders
            order_map, _ = perceive_bond_orders ( list ( elements ), list ( explicit.keys ( ) ) )
            for key in aromatic:
                o = order_map.get ( key, order_map.get ( ( key[1], key[0] ) ) )
                if o: explicit[key] = int ( o )
        except Exception:
            pass
    return [ ( i, j, explicit[( i, j )] ) for ( i, j ) in explicit ], aromatic

def _repair_valence_deficits ( elements, bonds, xyz ):
    """ A carbon left with valence 3 after the bond-order perception (Open
        Babel does not flag cationic rings such as thiazolium as aromatic,
        and the Kekule perceiver assumes neutral ring nitrogens) gets a
        double bond to the closest singly-bonded neighbour that also lacks
        valence -- a carbon (C=C), an oxygen (C=O), a nitrogen with valence 2
        (neutral C=N) or, last, a nitrogen with valence 3, which then becomes
        N+ (iminium/pyridinium/thiazolium). """
    bonds = [ list ( b ) for b in bonds ]
    def valences ( ):
        v = defaultdict ( int )
        for ( i, j, o ) in bonds: v[i] += o; v[j] += o
        return v
    for _ in range ( len ( elements ) ):
        v = valences ( )
        best = None
        for b in bonds:
            i, j, o = b
            if o != 1: continue
            for c, other in ( ( i, j ), ( j, i ) ):
                if elements[c] != "C" or v[c] != 3: continue
                e = elements[other]
                if   e == "C" and v[other] == 3: rank = 0          # C=C
                elif e == "O" and v[other] == 1: rank = 1          # C=O
                elif e == "N" and v[other] == 2: rank = 2          # C=N, neutral imine N
                elif e == "N" and v[other] == 3: rank = 3          # C=N+, iminium
                else: continue
                d = float ( np.linalg.norm ( np.asarray ( xyz[c] ) - np.asarray ( xyz[other] ) ) )
                if best is None or ( rank, d ) < best[0]: best = ( ( rank, d ), b )
        if best is None: break
        best[1][2] = 2
    return [ tuple ( b ) for b in bonds ]

def _aromatic_bonds ( elements, bonds ):
    """ Ring bonds of planar conjugated 5/6-rings from a Kekule structure:
        every ring atom carries a double bond (in or out of the ring -- a
        Kekule structure of a fused system can leave a ring with a single
        internal double bond), except, in a 5-ring, one N/O/S with a lone
        pair (pyrrole, furan, thiophene, purine N9). Enough to tell aromatic
        from aliphatic atoms for the class rules. """
    n = len ( elements )
    adj = _neighbors ( n, bonds )
    aromatic = set ( )
    for ring in _rings ( n, adj, max_size = 6 ):
        if len ( ring ) not in ( 5, 6 ): continue
        without = [ a for a in ring if not any ( o == 2 for ( j, o ) in adj[a] ) ]
        if len ( ring ) == 6 and without: continue
        if len ( ring ) == 5 and ( len ( without ) > 1 or ( without and elements[without[0]] not in ( "N", "O", "S" ) ) ): continue
        # . exocyclic C=O on every sp2 carbon is not aromaticity (e.g. quinones)
        exo_carbonyl = sum ( 1 for a in ring for ( j, o ) in adj[a] if o == 2 and j not in ring and elements[j] in ( "O", "S" ) )
        if exo_carbonyl >= 2: continue
        for a in ring:
            for ( j, o ) in adj[a]:
                if j in ring: aromatic.add ( ( min ( a, j ), max ( a, j ) ) )
    return aromatic

def _formal_charges ( elements, bonds ):
    valence = defaultdict ( int )
    for ( i, j, o ) in bonds:
        valence[i] += o; valence[j] += o
    charges = [ ]
    for i, e in enumerate ( elements ):
        v = valence[i]
        if   e == "N": charges.append ( v - 3 if v in ( 2, 4 ) else 0 )
        elif e == "O": charges.append ( v - 2 if v in ( 1, 3 ) else 0 )
        elif e == "S": charges.append ( -1 if v == 1 else ( 1 if v == 3 else 0 ) )
        elif e == "P": charges.append ( 1 if v == 4 else 0 )        # P+ - O- representation of P=O
        elif e == "C": charges.append ( 0 )
        else:          charges.append ( 0 )
    return charges

_ALLOWED_VALENCES = { "N": ( 3, 5 ), "O": ( 2, ), "S": ( 2, 4, 6 ), "Se": ( 2, 4, 6 ), "P": ( 3, 5 ), "F": ( 1, ),
                      "Cl": ( 1, 3, 5, 7 ), "Br": ( 1, 3, 5 ), "I": ( 1, 3, 5, 7 ) }

def _valences_complete ( elements, bonds, formal_charges = None ):
    """ True when every atom's bond-order sum matches an allowed valence for its
        formal charge (given ones, else none) -- i.e. no hydrogen is missing.
        Unknown elements -> False (the caller keeps the old behaviour). """
    total = [ 0 ] * len ( elements )
    for ( i, j, order ) in bonds:
        total[i] += int ( order ); total[j] += int ( order )
    charges = list ( formal_charges ) if formal_charges is not None else [ 0 ] * len ( elements )
    for element, v, q in zip ( elements, total, charges ):
        q = int ( q or 0 )
        if element in ( "C", "Si" ):
            if v != 4 - abs ( q ): return False
        elif element == "B":
            if v != 3 - q: return False
        elif element in _ALLOWED_VALENCES:
            if not any ( v == a + q for a in _ALLOWED_VALENCES[element] ): return False
        else:
            return False
    return True


def ligand_chemistry ( residue_name, atoms, ph = 7.4, add_hydrogens = None, obabel = None, bonds = None,
                       formal_charges = None, h_counts = None ):
    """ LigandParameters with atoms (names kept for the input atoms, new
        hydrogens named H1, H2, ...), bonds, formal charges. atoms: [ ( name,
        Z, xyz ) ] -- heavy atoms with or without hydrogens. add_hydrogens:
        None = only when the input has none. bonds/formal_charges: the
        input's own chemistry ( i, j, order ) -- used as is (no Open Babel)
        when given and the molecule already has hydrogens.
        h_counts: { heavy atom name: hydrogen count } chosen by the user
        (protonation review, 2026-10-18): the input hydrogens are dropped and
        rebuilt with these counts; bond orders and formal charges follow. """
    lig = LigandParameters ( residue_name )
    if h_counts:
        return _ligand_chemistry_native ( lig, atoms, ph, True, bonds, h_counts = h_counts )
    obabel = obabel or find_obabel ( )
    heavy = [ a for a in atoms if a[1] != 1 ]
    if len ( heavy ) == 1 and len ( atoms ) == 1:
        element = Z_ELEMENT.get ( heavy[0][1] )
        if element in OPLS_IONS:
            lig.atoms = [ ( heavy[0][0], heavy[0][1], np.asarray ( heavy[0][2], float ) ) ]
            lig.formal_charges = [ OPLS_IONS[element][0] ]
            lig.notes.append ( "monoatomic ion {}{:+d}".format ( element, OPLS_IONS[element][0] ) )
            return lig
        raise ValueError ( "Monoatomic residue {} ({}) has no known ion parameters.".format ( residue_name, element ) )
    has_h = any ( a[1] == 1 for a in atoms )
    # . [EN] 2026-10-17 BUG FIX: a molecule that really has no hydrogens (SO2(N3)2
    #   from the Builder, CCl4, CO2, SF6) was taken as "heavy atoms only" and got
    #   hydrogens rebuilt from the geometry (the azide nitrogens became amines).
    #   With the input's own bond orders, an H-free molecule whose valences are
    #   already complete is used as it is; deficits still mean missing hydrogens.
    if bonds and not has_h and add_hydrogens is not True:
        complete = _valences_complete ( [ Z_ELEMENT.get ( z, "X" ) for ( _, z, _ ) in atoms ], bonds, formal_charges )
    else:
        complete = False
    if bonds and ( has_h or complete ) and add_hydrogens is not True:
        counts = Counter ( a[0] for a in atoms )
        names, used = [ ], set ( )
        for ( name, z, xyz ) in atoms:
            new = name[:4] if counts[name] == 1 and name[:4] not in used else None
            if new is None:
                element, n = Z_ELEMENT.get ( z, "X" ), 1
                while "{}{}".format ( element, n )[:4] in used or "{}{}".format ( element, n ) in counts: n += 1
                new = "{}{}".format ( element, n )[:4]
            used.add ( new ); names.append ( new )
        lig.atoms = [ ( names[k], z, np.asarray ( xyz, float ) ) for k, ( _, z, xyz ) in enumerate ( atoms ) ]
        lig.bonds = [ ( i, j, int ( o ) ) for ( i, j, o ) in bonds ]
        elements = [ Z_ELEMENT.get ( z, "X" ) for ( _, z, _ ) in atoms ]
        lig.aromatic_bonds = _aromatic_bonds ( elements, lig.bonds )
        # . input charges where set (Builder), valence-derived otherwise (the pDynamo
        #   mol2 reader leaves formalCharge = 0 even on a quaternary N+)
        derived = _formal_charges ( elements, lig.bonds )
        given   = list ( formal_charges ) if formal_charges is not None else [ 0 ] * len ( derived )
        lig.formal_charges = [ g if g else d for g, d in zip ( given, derived ) ]
        lig.notes.append ( "bond orders and formal charges from the input" )
        return lig
    if add_hydrogens is None: add_hydrogens = not has_h
    # . [EN] 2026-10-14: native chemistry (pdynamo/ligand_perception.py) by default;
    #   Open Babel only when asked for (EASYHYBRID_LIGAND_CHEMISTRY=obabel) and found.
    #   Validation on 523 prepared ligands (DYFF_validation/scripts/ligprep_validate.py):
    #   native 507 identical to the reference, Open Babel 353.
    if os.environ.get ( "EASYHYBRID_LIGAND_CHEMISTRY", "native" ).lower ( ) != "obabel" or obabel is None:
        return _ligand_chemistry_native ( lig, atoms, ph, add_hydrogens, bonds if not has_h else None )
    if add_hydrogens is None: add_hydrogens = not has_h
    with tempfile.TemporaryDirectory ( ) as scratch:
        pdb  = os.path.join ( scratch, "lig.pdb" )
        mol2 = os.path.join ( scratch, "lig.mol2" )
        _write_pdb ( atoms if not add_hydrogens else heavy, pdb, residue_name )
        command = [ obabel, pdb, "-O", mol2 ]
        if add_hydrogens: command += [ "-p", "{:.1f}".format ( ph ) ]
        subprocess.run ( command, capture_output = True, text = True, timeout = 300 )
        if not os.path.exists ( mol2 ):
            raise RuntimeError ( "Open Babel failed for residue {}.".format ( residue_name ) )
        mol2_atoms, mol2_bonds = _parse_mol2 ( mol2 )
    source = atoms if not add_hydrogens else heavy
    # . Names: the input names where they are unique, otherwise element +
    #   counter (Open Babel names every hydrogen "H"), at most 4 characters.
    source_names = Counter ( a[0] for a in source )
    names, used = [ ], set ( )
    for k, ( name, element, xyz, tripos ) in enumerate ( mol2_atoms ):
        new = None
        if k < len ( source ) and source[k][1] == ELEMENT_Z.get ( element ) and source_names[source[k][0]] == 1:
            new = source[k][0][:4]
        if new is None or new in used:
            n = 1
            while "{}{}".format ( element, n )[:4] in used or ( "{}{}".format ( element, n ) in source_names ): n += 1
            new = "{}{}".format ( element, n )[:4]
        used.add ( new ); names.append ( new )
    elements = [ a[1] for a in mol2_atoms ]
    bonds, aromatic = _kekulize ( elements, mol2_bonds )
    bonds = _repair_valence_deficits ( elements, bonds, [ a[2] for a in mol2_atoms ] )
    lig.atoms = [ ( names[k], ELEMENT_Z.get ( elements[k], 0 ), mol2_atoms[k][2] ) for k in range ( len ( mol2_atoms ) ) ]
    lig.bonds = bonds
    lig.aromatic_bonds = aromatic
    lig.formal_charges = _formal_charges ( elements, bonds )
    if add_hydrogens:
        lig.notes.append ( "hydrogens added by Open Babel at pH {:.1f}".format ( ph ) )
    if any ( z == 0 for ( _, z, _ ) in lig.atoms ):
        raise ValueError ( "Residue {} has elements this perception does not handle.".format ( residue_name ) )
    return lig


_NATIVE_Z_SYMBOL = { 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 14: "Si", 15: "P", 16: "S", 17: "Cl", 33: "As", 34: "Se",
                     35: "Br", 52: "Te", 53: "I", 1: "H" }

def _ligand_chemistry_native ( lig, atoms, ph, add_hydrogens, bonds, h_counts = None ):
    """ ligand_chemistry() without Open Babel: pdynamo/ligand_perception.py
        (geometric hydrogen counts, Wang & Case orders, protonation rules at pH,
        VSEPR hydrogens). add_hydrogens: strip the input hydrogens and rebuild
        them; otherwise the input hydrogens are kept. Names: input names where
        unique, new hydrogens element + counter, at most 4 characters. """
    from pdynamo import ligand_perception as LP
    source = [ a for a in atoms if a[1] != 1 ] if add_hydrogens else list ( atoms )
    symbols = [ _NATIVE_Z_SYMBOL.get ( int ( z ), Z_ELEMENT.get ( int ( z ), "X" ) ) for ( _, z, _ ) in source ]
    if any ( s == "X" for s in symbols ):
        raise ValueError ( "Residue {} has elements the ligand perception does not handle.".format ( lig.residue_name ) )
    pairs = [ ( int ( b[0] ), int ( b[1] ) ) for b in bonds ] if ( bonds and not add_hydrogens ) or ( bonds and len ( source ) == len ( atoms ) ) else None
    override = None
    if h_counts:
        # . the user's hydrogen counts (by atom name); the input's heavy-atom bonds
        #   re-indexed onto the heavy atoms
        override = { k: int ( h_counts.get ( a[0], 0 ) ) for k, a in enumerate ( source ) }
        if bonds and pairs is None:
            heavy_index = { i: k for k, i in enumerate ( i for i, a in enumerate ( atoms ) if a[1] != 1 ) }
            pairs = [ ( heavy_index[int ( b[0] )], heavy_index[int ( b[1] )] ) for b in bonds
                      if int ( b[0] ) in heavy_index and int ( b[1] ) in heavy_index ] or None
    result = LP.perceive_ligand ( symbols, [ np.asarray ( a[2], float ) for a in source ], bonds = pairs, ph = ph,
                                  h_counts = override )
    counts = Counter ( a[0] for a in source )
    names, used = [ ], set ( )
    for k, element in enumerate ( result["symbols"] ):
        new = None
        if k < len ( source ) and counts[source[k][0]] == 1 and source[k][0][:4] not in used:
            new = source[k][0][:4]
        if new is None:
            element = element.capitalize ( )
            n = 1
            while "{}{}".format ( element, n )[:4] in used or "{}{}".format ( element, n ) in counts: n += 1
            new = "{}{}".format ( element, n )[:4]
        used.add ( new ); names.append ( new )
    lig.atoms = [ ( names[k], ELEMENT_Z.get ( result["symbols"][k].capitalize ( ), 0 ) or ( source[k][1] if k < len ( source ) else 1 ),
                    np.asarray ( result["coords"][k], float ) ) for k in range ( len ( result["symbols"] ) ) ]
    lig.bonds = [ ( int ( i ), int ( j ), int ( o ) ) for ( i, j, o ) in result["bonds"] ]
    elements = [ s.capitalize ( ) for s in result["symbols"] ]
    lig.aromatic_bonds = _aromatic_bonds ( elements, lig.bonds )
    lig.formal_charges = [ int ( q ) for q in result["formal_charges"] ]
    lig.notes.extend ( result["notes"] )
    for line in LP.format_decisions ( result, names ):
        lig.notes.append ( "protonation: " + line )
    return lig
