#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: simplified OPLS parametrization of small molecules (phase 4)
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      Parametrizes a residue that has no PDB component (ligand, cofactor,
#      metal ion) so that prep.py can build and type it:
#
#        1. chemistry: hydrogens at a given pH and bond orders from Open
#           Babel ("obabel -p pH"), aromatic bonds Kekulized, formal
#           charges from the valences;
#        2. OPLS-AA atom classes by rules (element, hybridization,
#           aromaticity, ring size, neighbours) -- the class names of the
#           OPLS-AA files (CT, CA, C, CM, CZ, CW, CS, CR, NA, NB, NC, N,
#           NT, N2, N3, NO, NZ, O, O2, OH, OS, ON, OY, S, SH, SZ, SY, P,
#           F, Cl, Br, I, HC, HA, H, H3, HO, HS);
#        3. charges: CM5 from GFN1-xTB (the "Mulliken/CM5 charges" block),
#           unscaled by default (CM5_SCALE_NEUTRAL; the Jorgensen group's
#           1.20/1.27 factors were derived for DFT CM5 and over-polarize
#           xTB CM5, see CM5_SCALE_NEUTRAL), averaged over
#           topologically equivalent atoms and corrected to the exact
#           total charge;
#        4. Lennard-Jones per class: the most common (sigma, epsilon) of
#           that class among the OPLS-AA 2008 atom types;
#        5. bonded terms: from the class database (parameters.py) when
#           available, otherwise ESTIMATED -- equilibrium values from the
#           geometry, force constants = median of the database terms with
#           the same elements, missing torsions = 0 -- and listed in the
#           report.
#
#      Output: a LigandParameters object with what prep.py needs: a PDB
#      component file (atoms, bonds, formal charges) for building, and an
#      MM sequence component (atom name -> type, charge) + LJ + type rows
#      for the parameter set. Monoatomic metal ions get their OPLS-AA ion
#      parameters directly (no QC).
#
#      This is a SIMPLIFIED parametrization: types by rules, no torsion
#      fitting, xTB-level CM5 (CM5 was defined on DFT densities). Check the
#      report before production simulations.
#

import os
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict, OrderedDict

import numpy as np

# . OPLS-AA X-CA-CA-X ring torsion, V2 = 7.25 kcal/mol -> pDynamo k = V2 / 2.
AROMATIC_RING_TORSION_K = 3.625

# . Scale of the GFN1-xTB CM5 charges. 1.20 (the Jorgensen group's universal
#   factor) was derived for CM5 from DFT densities; xTB CM5 is already more
#   polarized: against the OPLS-AA 2008 reference charges of 14 small
#   molecules (tests/benchmarks/opls_ligand_charges.py) the MAE is 0.054 e
#   unscaled, 0.082 e with 1.20, and the least-squares optimum is 0.835.
CM5_SCALE_NEUTRAL = 1.00
CM5_SCALE_CHARGED = 1.00

ELEMENT_Z = { "H": 1, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Na": 11, "Mg": 12, "Si": 14, "P": 15, "S": 16,
              "Cl": 17, "K": 19, "Ca": 20, "Zn": 30, "Br": 35, "I": 53, "Li": 3, "Rb": 37, "Sr": 38, "Cs": 55, "Ba": 56 }
Z_ELEMENT = { z: s for s, z in ELEMENT_Z.items ( ) }

# . Monoatomic ions with OPLS-AA 2008 parameters: element -> ( formal charge, Tinker type ).
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
    """ Everything prep.py/parameters.py need for one residue name. """
    def __init__ ( self, residue_name ):
        self.residue_name   = residue_name.upper ( )[:3]
        self.atoms          = [ ]      # [ ( name, Z, xyz ) ]
        self.bonds          = [ ]      # [ ( i, j, order ) ]
        self.aromatic_bonds = set ( )  # { ( i, j ) }
        self.formal_charges = [ ]      # per atom
        self.classes        = [ ]      # OPLS class per atom
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
        project's own Wang & Case perceiver (vismol), keeping the explicit
        orders. Falls back to alternating guesses only if it fails. """
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

def ligand_chemistry ( residue_name, atoms, ph = 7.4, add_hydrogens = None, obabel = None, bonds = None,
                       formal_charges = None ):
    """ LigandParameters with atoms (names kept for the input atoms, new
        hydrogens named H1, H2, ...), bonds, formal charges. atoms: [ ( name,
        Z, xyz ) ] -- heavy atoms with or without hydrogens. add_hydrogens:
        None = only when the input has none. bonds/formal_charges: the
        input's own chemistry ( i, j, order ) -- used as is (no Open Babel)
        when given and the molecule already has hydrogens. """
    lig = LigandParameters ( residue_name )
    obabel = obabel or find_obabel ( )
    heavy = [ a for a in atoms if a[1] != 1 ]
    if len ( heavy ) == 1 and len ( atoms ) == 1:
        element = Z_ELEMENT.get ( heavy[0][1] )
        if element in OPLS_IONS:
            lig.atoms = [ ( heavy[0][0], heavy[0][1], np.asarray ( heavy[0][2], float ) ) ]
            lig.formal_charges = [ OPLS_IONS[element][0] ]
            lig.notes.append ( "monoatomic ion {}{:+d}".format ( element, OPLS_IONS[element][0] ) )
            return lig
        raise ValueError ( "Monoatomic residue {} ({}) has no OPLS ion parameters.".format ( residue_name, element ) )
    has_h = any ( a[1] == 1 for a in atoms )
    if bonds and has_h and add_hydrogens is not True:
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
    if obabel is None:
        raise RuntimeError ( "Open Babel (obabel) was not found; it is needed for ligand hydrogens and bond orders." )
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
        raise ValueError ( "Residue {} has elements without OPLS support.".format ( residue_name ) )
    return lig


#=============================================================================
#  OPLS classes
#=============================================================================
def assign_classes ( lig ):
    n = len ( lig.atoms )
    el = [ Z_ELEMENT.get ( z, "X" ) for ( _, z, _ ) in lig.atoms ]
    adj = _neighbors ( n, lig.bonds )
    rings = _rings ( n, adj )
    arom_atoms = set ( )
    for ( i, j ) in lig.aromatic_bonds: arom_atoms.update ( ( i, j ) )
    def ring_size ( i ):
        sizes = [ len ( r ) for r in rings if i in r and r <= arom_atoms ]
        return min ( sizes ) if sizes else None
    def heavy_nbrs ( i ): return [ j for ( j, o ) in adj[i] if el[j] != "H" ]
    def n_h ( i ): return sum ( 1 for ( j, o ) in adj[i] if el[j] == "H" )
    def has_double_to ( i, element ): return any ( o == 2 and el[j] == element for ( j, o ) in adj[i] )
    def max_order ( i ): return max ( [ o for ( _, o ) in adj[i] ] or [ 0 ] )
    classes = [ None ] * n
    q = lig.formal_charges
    # . Heavy atoms.
    for i in range ( n ):
        e = el[i]
        if e == "C":
            if max_order ( i ) == 3 or sum ( 1 for ( _, o ) in adj[i] if o == 2 ) == 2: classes[i] = "CZ"
            elif i in arom_atoms:
                size = ring_size ( i )
                hetero = [ j for j in heavy_nbrs ( i ) if el[j] in ( "N", "O", "S" ) and j in arom_atoms ]
                n_aromatic_rings = sum ( 1 for r in rings if i in r and r <= arom_atoms )
                if n_aromatic_rings >= 2 and any ( len ( r ) == 5 for r in rings if i in r and r <= arom_atoms ):
                    classes[i] = "CB"                       # 5/6 fusion carbon (purine C4/C5, indole C3a/C7a)
                elif size == 5:
                    ring5 = [ r for r in rings if i in r and len ( r ) == 5 and r <= arom_atoms ]
                    fused = any ( len ( r & other ) >= 2 for r in ring5 for other in rings
                                  if other is not r and other <= arom_atoms and other != r )
                    if len ( hetero ) >= 2: classes[i] = "CK" if fused else "CR"   # purine C8 / imidazole C2
                    elif hetero:
                        # . next to a pyridine-like N (2 connections, no H): CV (imidazole/oxazole/
                        #   thiazole C4); next to N-H, O or S: CW (imidazole C5, pyrrole/furan C2)
                        j = hetero[0]
                        classes[i] = "CV" if el[j] == "N" and len ( adj[j] ) == 2 else "CW"
                    else:                   classes[i] = "CS"
                elif len ( [ j for j in hetero if el[j] == "N" ] ) >= 2:
                    classes[i] = "CQ"                       # C between two ring N (pyrimidine/purine C2)
                else:
                    classes[i] = "CA"
            elif has_double_to ( i, "O" ) or has_double_to ( i, "S" ): classes[i] = "C"
            elif any ( o == 2 for ( _, o ) in adj[i] ):
                # . C=N of guanidinium/amidinium: OPLS uses CA
                if has_double_to ( i, "N" ) and sum ( 1 for j in heavy_nbrs ( i ) if el[j] == "N" ) >= 2: classes[i] = "CA"
                else: classes[i] = "CM"
            else: classes[i] = "CT"
        elif e == "N":
            hn = heavy_nbrs ( i )
            if max_order ( i ) == 3: classes[i] = "NZ"
            elif sum ( 1 for j in hn if el[j] == "O" ) >= 2: classes[i] = "NO"
            elif i in arom_atoms:
                ring_partners = { j for r in rings if i in r for j in r }
                if len ( adj[i] ) == 3 and any ( el[j] == "C" and j not in ring_partners and j not in arom_atoms for j in hn ):
                    classes[i] = "N*"                       # aromatic N bonded to an sp3 C (nucleoside N9/N1)
                elif len ( adj[i] ) == 3: classes[i] = "NA"
                else: classes[i] = "NB" if ring_size ( i ) == 5 else "NC"
            elif q[i] == 1 and len ( adj[i] ) == 4: classes[i] = "N3"
            elif any ( el[j] == "C" and any ( o == 2 and el[k] in ( "O", "S" ) for ( k, o ) in adj[j] ) for j in hn ):
                classes[i] = "N"
            elif any ( el[j] == "C" and sum ( 1 for ( k, o ) in adj[j] if el[k] == "N" ) >= 2 and any ( o2 == 2 for ( _, o2 ) in adj[j] ) for j in hn ):
                classes[i] = "N2"
            elif max_order ( i ) == 2: classes[i] = "NB"
            else: classes[i] = "NT"
        elif e == "O":
            hn = heavy_nbrs ( i )
            partner = hn[0] if hn else None
            if n_h ( i ) == 2: classes[i] = "OW"
            elif n_h ( i ) == 1: classes[i] = "OH"
            elif partner is not None and el[partner] == "N" and classes[partner] in ( None, "NO" ) and \
                 sum ( 1 for j in heavy_nbrs ( partner ) if el[j] == "O" ) >= 2: classes[i] = "ON"
            elif partner is not None and el[partner] == "S" and len ( hn ) == 1: classes[i] = "OY"
            elif partner is not None and el[partner] == "P" and len ( hn ) == 1: classes[i] = "O2"
            elif len ( hn ) == 1 and partner is not None and el[partner] == "C":
                # . carboxylate: C bonded to two terminal O, one of them charged
                terminal_o = [ j for ( j, o ) in adj[partner] if el[j] == "O" and len ( adj[j] ) == 1 ]
                if len ( terminal_o ) == 2 and any ( q[j] == -1 for j in terminal_o ): classes[i] = "O2"
                else: classes[i] = "O"
            else: classes[i] = "OS"
        elif e == "S":
            n_o = sum ( 1 for j in heavy_nbrs ( i ) if el[j] == "O" and len ( adj[j] ) == 1 )
            if   n_h ( i ) >= 1: classes[i] = "SH"
            elif n_o >= 2:       classes[i] = "SY"
            elif n_o == 1:       classes[i] = "SZ"
            else:                classes[i] = "S"
        elif e == "P":  classes[i] = "P"
        elif e in ( "F", "Cl", "Br", "I", "Si", "B" ): classes[i] = e
        elif e in OPLS_IONS: classes[i] = e
    # . Hydrogens.
    for i in range ( n ):
        if el[i] != "H": continue
        parent = adj[i][0][0] if adj[i] else None
        pe, pc = ( el[parent], classes[parent] ) if parent is not None else ( None, None )
        if   pe == "C": classes[i] = "HA" if pc in ( "CA", "CW", "CS", "CR", "CV", "CQ", "CB", "CK" ) else "HC"
        elif pe == "N": classes[i] = "H3" if pc == "N3" else "H"
        elif pe == "O": classes[i] = "HW" if pc == "OW" else "HO"
        elif pe == "S": classes[i] = "HS"
        else:           classes[i] = "HC"
    if any ( c is None for c in classes ):
        raise ValueError ( "No OPLS class for atoms {} of {}.".format (
                [ lig.atoms[i][0] for i, c in enumerate ( classes ) if c is None ], lig.residue_name ) )
    lig.classes = classes
    lig.equivalence = equivalence_classes ( el, adj )
    return classes


#=============================================================================
#  Lennard-Jones per class
#=============================================================================
_CLASS_LJ = None

def class_lennard_jones ( ):
    """ { class: ( epsilon, sigma ) } -- most common pair among the
        OPLS-AA 2008 types of the class (united-atom and GBSA types out). """
    global _CLASS_LJ
    if _CLASS_LJ is None:
        from pdynamo.opls.parameters import TINKER_FILES
        from pdynamo.opls.tinker_prm import TinkerPRM
        path = dict ( TINKER_FILES )["oplsaa08"]
        prm = TinkerPRM ( path )
        counts = defaultdict ( Counter )
        for t, data in prm.atoms.items ( ):
            if "(UA)" in data["description"] or "GBSA" in data["description"]: continue
            lj = prm.type_lj ( t )
            if lj is None: continue
            counts[data["symbol"]][( round ( lj[1], 5 ), round ( lj[0], 5 ) )] += 1
        _CLASS_LJ = { s: c.most_common ( 1 )[0][0] for s, c in counts.items ( ) }
        for element, ( _, atom_type ) in OPLS_IONS.items ( ):
            lj = prm.type_lj ( atom_type )
            if lj is not None: _CLASS_LJ["ION_" + element] = ( lj[1], lj[0] )
    return _CLASS_LJ


#=============================================================================
#  Charges
#=============================================================================
def cm5_charges ( lig, xtb = None, scale = None, keep_scratch = None ):
    """ GFN1-xTB CM5 charges (scaled, symmetrized, exact total). """
    xtb = xtb or find_xtb ( )
    if xtb is None:
        raise RuntimeError ( "xtb was not found (set XTBHOME or put xtb in the PATH); it is needed for CM5 charges." )
    total = lig.total_charge
    check_multiplicity ( lig )
    if scale is None: scale = CM5_SCALE_NEUTRAL if total == 0 else CM5_SCALE_CHARGED
    scratch = tempfile.mkdtemp ( prefix = "opls_cm5_" )
    try:
        xyz = os.path.join ( scratch, "lig.xyz" )
        with open ( xyz, "w" ) as handle:
            handle.write ( "{}\n{}\n".format ( len ( lig.atoms ), lig.residue_name ) )
            for ( name, z, r ) in lig.atoms:
                handle.write ( "{} {:.6f} {:.6f} {:.6f}\n".format ( Z_ELEMENT[z], r[0], r[1], r[2] ) )
        command = [ xtb, xyz, "--gfn", "1", "--chrg", str ( total ) ]
        if lig.multiplicity > 1: command += [ "--uhf", str ( lig.multiplicity - 1 ) ]
        result = subprocess.run ( command, cwd = scratch,
                                  capture_output = True, text = True, timeout = 1800,
                                  env = dict ( os.environ, OMP_NUM_THREADS = os.environ.get ( "OMP_NUM_THREADS", "2" ) ) )
        text = result.stdout
        start = text.find ( "Mulliken/CM5 charges" )
        if start < 0:
            raise RuntimeError ( "xtb gave no CM5 charges for {}:\n{}".format ( lig.residue_name, ( result.stdout + result.stderr )[-1500:] ) )
        charges = [ ]
        for line in text[start:].splitlines ( )[1:len ( lig.atoms ) + 1]:
            charges.append ( float ( line.split ( )[2] ) )
        if len ( charges ) != len ( lig.atoms ):
            raise RuntimeError ( "Could not read the CM5 charges of {}.".format ( lig.residue_name ) )
        if keep_scratch:
            shutil.copy ( os.path.join ( scratch, "lig.xyz" ), keep_scratch + ".xyz" )
            with open ( keep_scratch + "_xtb.log", "w" ) as handle: handle.write ( text )
    finally:
        shutil.rmtree ( scratch, ignore_errors = True )
    charges = np.array ( charges ) * scale
    # . Symmetrize over topologically equivalent atoms.
    groups = defaultdict ( list )
    for i, g in enumerate ( lig.equivalence ): groups[g].append ( i )
    for members in groups.values ( ):
        charges[members] = charges[members].mean ( )
    # . Exact total: spread the residual over all atoms.
    charges += ( total - charges.sum ( ) ) / len ( charges )
    # . 6 decimals (as written to the parameter files) with an exact total: the
    #   rounding residual (< 1e-5) goes to the atom with the largest |q|
    charges = np.round ( charges, 6 )
    charges[int ( np.argmax ( np.abs ( charges ) ) )] += round ( total - float ( charges.sum ( ) ), 6 )
    lig.charges = [ float ( c ) for c in charges ]
    lig.charge_method = "CM5 (GFN1-xTB) x {:.2f}".format ( scale )
    if lig.multiplicity > 1:
        lig.charge_method += ", multiplicity {}".format ( lig.multiplicity )
    return lig.charges


#=============================================================================
#  Labels, LJ
#=============================================================================
def finalize_types ( lig ):
    """ Type labels CLASS.RESn (one per equivalence class) and LJ. """
    lj = class_lennard_jones ( )
    if len ( lig.atoms ) == 1 and lig.notes and lig.notes[0].startswith ( "monoatomic" ):
        element = Z_ELEMENT[lig.atoms[0][1]]
        label = "{}.{}".format ( element, lig.residue_name )
        lig.classes, lig.labels = [ element ], [ label ]
        lig.charges = [ float ( lig.formal_charges[0] ) ]
        lig.lj = { label: lj["ION_" + element] }
        lig.charge_method = "formal charge (OPLS-AA ion)"
        lig.equivalence = [ 0 ]
        return
    labels = [ ]
    for i, ( c, g ) in enumerate ( zip ( lig.classes, lig.equivalence ) ):
        labels.append ( "{}.{}{}".format ( c, lig.residue_name, g + 1 ) )
    lig.labels = labels
    for label, c in zip ( labels, lig.classes ):
        if c not in lj:
            lig.warnings.append ( "no LJ for class {} -- using CT values".format ( c ) )
        lig.lj[label] = lj.get ( c, lj["CT"] )

def check_multiplicity ( lig ):
    """ Raises ValueError when the number of electrons (sum of Z minus the
        total charge) and the multiplicity are incompatible (an even number
        of electrons needs an odd multiplicity and vice versa). """
    m = int ( lig.multiplicity )
    if m < 1:
        raise ValueError ( "{}: the multiplicity must be 1 or larger.".format ( lig.residue_name ) )
    if ( lig.n_electrons + m - 1 ) % 2 != 0:
        raise ValueError ( "{}: charge {:+d} gives {} electrons, incompatible with multiplicity {} "
                           "({} electrons need an {} multiplicity).".format (
                           lig.residue_name, lig.total_charge, lig.n_electrons, m, "even" if lig.n_electrons % 2 == 0 else "odd",
                           "odd" if lig.n_electrons % 2 == 0 else "even" ) )

def apply_user_charge ( lig, total_charge = None, multiplicity = None ):
    """ Total charge / multiplicity chosen by the user. The formal charges
        of the atoms (used for building and for the atom classes) are not
        changed: a different total only changes the number of electrons of
        the QC charge calculation and the normalization of the charges. To
        change protons, change the pH or the hydrogens instead. """
    if total_charge is not None and int ( total_charge ) != lig.formal_charge_sum:
        lig.charge_override = int ( total_charge )
        lig.warnings.append ( "total charge {:+d} set by the user (perceived formal charges sum to {:+d})".format (
                lig.charge_override, lig.formal_charge_sum ) )
    if multiplicity is not None and int ( multiplicity ) != 1:
        lig.multiplicity = int ( multiplicity )
        lig.notes.append ( "spin multiplicity {} set by the user".format ( lig.multiplicity ) )
    check_multiplicity ( lig )
    return lig

def parametrize_residue ( residue_name, atoms, ph = 7.4, xtb = None, cm5_scale = None, add_hydrogens = None,
                          keep_scratch = None, bonds = None, formal_charges = None, total_charge = None,
                          multiplicity = None ):
    """ Full simplified parametrization of one residue. total_charge /
        multiplicity: the user's values (see apply_user_charge()). """
    lig = ligand_chemistry ( residue_name, atoms, ph = ph, add_hydrogens = add_hydrogens, bonds = bonds,
                             formal_charges = formal_charges )
    if len ( lig.atoms ) == 1 and ( total_charge is not None or multiplicity not in ( None, 1 ) ):
        lig.warnings.append ( "charge/multiplicity of a monoatomic ion come from the OPLS-AA ion parameters (user values ignored)" )
    elif len ( lig.atoms ) > 1:
        apply_user_charge ( lig, total_charge, multiplicity )
    if len ( lig.atoms ) > 1:
        assign_classes ( lig )
        cm5_charges ( lig, xtb = xtb, scale = cm5_scale, keep_scratch = keep_scratch )
    finalize_types ( lig )
    return lig


#=============================================================================
#  Missing bonded terms
#=============================================================================
def estimate_missing_terms ( system, types, database, base_set = None ):
    """ { term: { key: value } } for the bonds/angles/dihedrals of `system`
        that neither the base set nor the class database has: equilibrium
        values from the current geometry, force constants = median over the
        database terms with the same element tuple (fallback 300 / 50),
        torsions = no term (empty list). Also returns notes. """
    from pdynamo.opls import parameters as P
    used = P.system_terms ( system, types )
    crd  = system.coordinates3
    z    = [ a.atomicNumber for a in system.atoms ]
    element_of_class = { }
    for label, z_ in zip ( types, z ): element_of_class[P.atom_class ( label )] = z_
    def elements ( key ): return tuple ( sorted ( element_of_class.get ( P.atom_class ( l ), 0 ) for l in key ) )
    medians = { "bond": defaultdict ( list ), "angle": defaultdict ( list ) }
    for term in ( "bond", "angle" ):
        for key, entry in database.entries[term].items ( ):
            medians[term][elements ( key )].append ( entry["value"][1] )
    # . geometry per key
    geometry = { "bond": defaultdict ( list ), "angle": defaultdict ( list ) }
    c = system.connectivity
    def xyz ( i ): return np.array ( [ crd[i, 0], crd[i, 1], crd[i, 2] ] )
    for ( i, j ) in c.bondIndices:
        geometry["bond"][P.key_bond ( types[i], types[j] )].append ( float ( np.linalg.norm ( xyz ( i ) - xyz ( j ) ) ) )
    for ( i, j, k ) in c.angleIndices:
        u, v = xyz ( i ) - xyz ( j ), xyz ( k ) - xyz ( j )
        cosine = np.dot ( u, v ) / ( np.linalg.norm ( u ) * np.linalg.norm ( v ) )
        geometry["angle"][P.key_angle ( types[i], types[j], types[k] )].append ( float ( np.degrees ( np.arccos ( np.clip ( cosine, -1, 1 ) ) ) ) )
    # . dihedrals whose central bond is aromatic (both central atoms aromatic, in one ring)
    aromatic_central = set ( )
    for ( i, j, k, l ) in c.dihedralIndices:
        if getattr ( system.atoms[j], "isAromatic", False ) and getattr ( system.atoms[k], "isAromatic", False ):
            aromatic_central.add ( P.key_dihedral ( types[i], types[j], types[k], types[l] ) )
    out = { "bond": { }, "angle": { }, "dihedral": { }, "outofplane": { } }
    notes = [ ]
    # . out-of-plane terms: OPLS-AA impropers for planar centres (Tinker V / 2):
    #   carbonyl C 10.5, aromatic C/N 1.1, amide-like N 1.0; others stay 0
    central_aromatic = { }
    central_carbonyl = { }
    for ( i, j, k, l ) in c._IndicesNeighbor3 ( ):
        key = P.key_outofplane ( types[i], types[j], types[k], types[l] )
        central_aromatic[key] = bool ( getattr ( system.atoms[i], "isAromatic", False ) )
        central_carbonyl[key] = system.atoms[i].atomicNumber == 6 and any (
                system.atoms[n].atomicNumber == 8 and len ( c.adjacentNodes[system.atoms[n]] ) == 1 for n in ( j, k, l ) )
    for key in sorted ( used["outofplane"] ):
        if base_set is not None and key in base_set.terms["outofplane"]: continue
        value, _, _ = database.lookup ( "outofplane", key )
        if value is not None: continue
        central = P.atom_class ( key[0] )
        if central_carbonyl.get ( key ):      K = 10.5
        elif central_aromatic.get ( key ):    K = 1.1
        elif central in ( "N", "N2", "NA", "N*" ): K = 1.0
        else: continue
        out["outofplane"][key] = K
        notes.append ( "outofplane {} : generic OPLS improper (K {:.1f})".format ( "-".join ( key ), K ) )
    for term in ( "bond", "angle", "dihedral" ):
        for key in sorted ( used[term] ):
            if base_set is not None and key in base_set.terms[term]: continue
            value, _, _ = database.lookup ( term, key )
            if value is not None: continue
            if term == "dihedral":
                if key in aromatic_central:
                    out[term][key] = [ ( AROMATIC_RING_TORSION_K, 2, 180.0 ) ]
                    notes.append ( "dihedral {} : generic aromatic ring torsion (k {:.3f}, n 2)".format ( "-".join ( key ), AROMATIC_RING_TORSION_K ) )
                else:
                    out[term][key] = [ ]
                    notes.append ( "dihedral {} : no parameters, set to zero".format ( "-".join ( key ) ) )
                continue
            ks = medians[term].get ( elements ( key ) )
            k_value = float ( np.median ( ks ) ) if ks else ( 300.0 if term == "bond" else 50.0 )
            eq = float ( np.mean ( geometry[term][key] ) )
            out[term][key] = ( round ( eq, 4 ), round ( k_value, 2 ) )
            notes.append ( "{} {} : estimated ({:.3f}, k {:.1f})".format ( term, "-".join ( key ), eq, k_value ) )
    return out, notes
