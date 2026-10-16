#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: OPLS system preparation -- phase 0 (no GUI)
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      Turns a system loaded in EasyHybrid (typically from a PDB file,
#      with or without hydrogens) into a pDynamo system that the OPLS
#      atom typer can type:
#
#        analyze_system ( system )            -> PrepAnalysis
#        build_prepared_system ( analysis )   -> ( new System, report )
#        typing_report ( system )             -> dict
#
#      The chemistry (bond orders, formal charges, termini, protonation
#      states, proline/disulfide links) comes from pDynamo's own PDB
#      component library through its PDBModel route -- the same route the
#      pDynamo PDB examples use -- fed with a clean PDB written here:
#
#        - one SEGMENT per polymer piece: a chain is split wherever the
#          peptide C(i)-N(i+1) distance shows a break, and every piece
#          ends with a TER record (that is what makes pDynamo create the
#          linear polymer: peptide links, N/C termini, proline links);
#        - waters, ions and other non-polymer components in segments of
#          their own;
#        - disulfide bridges found by geometry (SG-SG) added as model links;
#        - protonation variants given per residue (defaults: pH 7 library
#          states, or the state implied by the hydrogens already present).
#
#      Residue/atom name aliases from AMBER/CHARMM/GROMACS files (HID, HSE,
#      CYX, ASH, WAT, TIP3, OT1/OT2, ILE CD, ...) are mapped to the PDB
#      component library names. Residues without a library component
#      (ligands, cofactors, nucleic acids) cannot go through this route:
#      they are reported as "unsupported" (to be parametrized separately).
#

import os
import tempfile
from collections import Counter, OrderedDict

import numpy as np


# . Extra PDB component library entries shipped with EasyHybrid (CYM, LYN).
EXTRA_LIBRARY_PATH = os.path.join ( os.path.dirname ( os.path.abspath ( __file__ ) ), "pdb_components_extra" )

# . Default OPLS parameter set (pDynamo3's own).
DEFAULT_PARAMETER_SET = "protein"

PEPTIDE_BOND_MAX  = 2.0      # A, C(i)-N(i+1) above this = chain break
DISULFIDE_SG_MAX  = 2.5      # A, SG-SG below this = disulfide bridge

# . Residue name aliases -> ( PDB library name, protonation variant or None ).
RESIDUE_ALIASES = {
    "HID" : ( "HIS", "Delta Protonated"   ), "HSD" : ( "HIS", "Delta Protonated"   ),
    "HIE" : ( "HIS", "Epsilon Protonated" ), "HSE" : ( "HIS", "Epsilon Protonated" ),
    "HIP" : ( "HIS", "Doubly Protonated"  ), "HSP" : ( "HIS", "Doubly Protonated"  ),
    "HSH" : ( "HIS", "Doubly Protonated"  ),
    "ASH" : ( "ASP", "Protonated"   ), "ASPP": ( "ASP", "Protonated"   ),
    "GLH" : ( "GLU", "Protonated"   ), "GLUP": ( "GLU", "Protonated"   ),
    "LYN" : ( "LYS", "Neutral"      ), "LSN" : ( "LYS", "Neutral"      ),
    "CYM" : ( "CYS", "Deprotonated" ), "CYX" : ( "CYS", None ), "CYS2": ( "CYS", None ),
    "WAT" : ( "HOH", None ), "TIP3": ( "HOH", None ), "TIP" : ( "HOH", None ),
    "SOL" : ( "HOH", None ), "H2O" : ( "HOH", None ), "SPC" : ( "HOH", None ),
    "SOD" : ( "NA", None ), "NA+" : ( "NA", None ), "Na+" : ( "NA", None ),
    "POT" : ( "K",  None ), "K+"  : ( "K",  None ),
    "CLA" : ( "CL", None ), "CL-" : ( "CL", None ), "Cl-" : ( "CL", None ),
}

# . Atom name aliases: ( residue or "*", name ) -> library name.
ATOM_ALIASES = {
    ( "*",   "OT1" ) : "O",   ( "*",   "OT2" ) : "OXT", ( "*", "OC1" ) : "O", ( "*", "OC2" ) : "OXT",
    ( "*",   "HN"  ) : "H",   ( "ILE", "CD"  ) : "CD1",
    ( "HOH", "OW"  ) : "O",   ( "HOH", "OH2" ) : "O",
    ( "HOH", "HW1" ) : "H1",  ( "HOH", "HW2" ) : "H2",
    ( "NA",  "SOD" ) : "NA",  ( "K",   "POT" ) : "K",   ( "CL", "CLA" ) : "CL",
}

# . Residues whose protonation state is a choice.
PROTONATION_CHOICES = {
    "HIS" : [ "Delta Protonated", "Epsilon Protonated", "Doubly Protonated" ],
    "ASP" : [ "Deprotonated", "Protonated" ],
    "GLU" : [ "Deprotonated", "Protonated" ],
    "LYS" : [ None, "Neutral" ],                 # None = library form (NH3+)
    "CYS" : [ None, "Deprotonated" ],            # None = library form (SH)
}


#=============================================================================
#  Data classes
#=============================================================================
class PrepResidue:
    """ One residue of the input system. """
    def __init__ ( self ):
        self.entity        = None      # original entity (chain) label
        self.name          = None      # original residue name
        self.number        = None      # residue number (int)
        self.icode         = ""
        self.library_name  = None      # PDB component library name, None = unsupported
        self.kind          = None      # "amino acid" | "water" | "ion" | "other" | "unsupported"
        self.atoms         = [ ]       # [ ( library atom name, atomic number, xyz ) ]
        self.n_hydrogens   = 0
        self.variant       = None      # protonation variant to apply (None = library default)
        self.variant_source = "library default"
        self.segment       = None      # segment label in the prepared system
        self.disulfide_partner = None  # index of the partner PrepResidue
        self.missing_heavy = [ ]       # library heavy atoms absent from the file
        self.added_atoms   = [ ]       # atoms placed by this module (e.g. OXT)
        self.bonds         = None      # [ ( i, j, order ) ] inside the residue when the input has bond orders
        self.source_atoms  = None      # the input system's atoms of this residue
        self.notes         = ""
        self.formal_charges = None     # per atom, from the input system

    @property
    def label ( self ):
        return "{}:{}.{}{}".format ( self.entity, self.name, self.number, self.icode )

    def atom_xyz ( self, name ):
        for ( atom_name, _z, xyz ) in self.atoms:
            if atom_name == name: return xyz
        return None


class PrepAnalysis:
    """ Result of analyze_system(). """
    def __init__ ( self ):
        self.label       = None
        self.residues    = [ ]                      # [PrepResidue]
        self.segments    = OrderedDict ( )          # segment -> [residue indices]
        self.polymers    = [ ]                      # segments that are linear polymers
        self.chain_breaks = [ ]                     # [ ( residue i, residue i+1, distance ) ]
        self.disulfides  = [ ]                      # [ ( residue i, residue j, distance ) ]
        self.unsupported = [ ]                      # residue indices
        self.ligands     = OrderedDict ( )          # residue name -> ligand.LigandParameters
        self.ligand_options = { }                   # residue name -> { "charge": int or None, "multiplicity": int }
        self.extra_library_paths = [ ]              # per-system PDB component folders (ligands)
        self.ligand_warnings = [ ]

    def summary ( self ):
        kinds = Counter ( r.kind for r in self.residues )
        return { "residues"     : len ( self.residues ),
                 "kinds"        : dict ( kinds ),
                 "segments"     : { s: len ( v ) for s, v in self.segments.items ( ) },
                 "polymers"     : list ( self.polymers ),
                 "chain_breaks" : [ ( self.residues[i].label, self.residues[j].label, round ( d, 2 ) ) for i, j, d in self.chain_breaks ],
                 "disulfides"   : [ ( self.residues[i].label, self.residues[j].label, round ( d, 2 ) ) for i, j, d in self.disulfides ],
                 "unsupported"  : [ self.residues[i].label for i in self.unsupported ],
                 "ligands"      : { name: lig.summary ( ) for name, lig in self.ligands.items ( ) },
                 "missing_heavy": { r.label: list ( r.missing_heavy ) for r in self.residues if r.missing_heavy },
                 "added_atoms"  : { r.label: list ( r.added_atoms ) for r in self.residues if r.added_atoms },
                 "variants"     : { self.residues[i].label: ( r.variant, r.variant_source )
                                    for i, r in enumerate ( self.residues ) if r.library_name in PROTONATION_CHOICES } }


#=============================================================================
#  Library helpers
#=============================================================================
_LIBRARY = None

def library ( ):
    """ PDB component library: EasyHybrid extras first, then pDynamo's own. """
    global _LIBRARY
    if _LIBRARY is None:
        from pBabel import PDBComponentLibrary
        standard = PDBComponentLibrary.WithOptions ( ).MakeStandardLibraryPath ( )
        _LIBRARY = PDBComponentLibrary.WithOptions ( paths = [ EXTRA_LIBRARY_PATH, standard ] )
    return _LIBRARY

def library_paths ( ):
    return list ( library ( ).libraryPaths )

def _component ( name ):
    try:
        return library ( ).GetComponent ( name )
    except Exception:
        return None


#=============================================================================
#  Analysis
#=============================================================================
def _residue_name_and_number ( component ):
    """ ( name, number, icode ) from a pDynamo sequence component label
        such as "ALA.12" or "ALA.12A". """
    label = component.label
    tokens = label.split ( "." )
    name = tokens[0]
    number, icode = 0, ""
    if len ( tokens ) > 1:
        text = tokens[1]
        digits = text.lstrip ( "-" )
        k = len ( text ) - len ( digits )
        j = 0
        while j < len ( digits ) and digits[j].isdigit ( ): j += 1
        try:
            number = int ( text[:k + j] )
        except ValueError:
            number = 0
        icode = digits[j:]
    return name, number, icode

def _residue_chemistry ( system, atoms ):
    """ ( bonds [ ( i, j, order ) ], formal charges ) of a list of atoms
        when the input connectivity has DEFINED bond types for all their
        internal bonds (mol2, Builder), else ( None, None ). """
    connectivity = getattr ( system, "connectivity", None )
    if connectivity is None or not connectivity.bonds: return None, None
    atoms = list ( atoms )
    local = { id ( a ): k for k, a in enumerate ( atoms ) }
    orders = { "Single": 1, "Double": 2, "Triple": 3 }
    bonds = [ ]
    for bond in connectivity.bonds:
        i, j = local.get ( id ( bond.node1 ) ), local.get ( id ( bond.node2 ) )
        if i is None or j is None: continue
        order = orders.get ( bond.type.name if hasattr ( bond.type, "name" ) else str ( bond.type ).split ( "." )[-1] )
        if order is None: return None, None
        bonds.append ( ( i, j, order ) )
    if not bonds: return None, None
    return bonds, [ int ( getattr ( a, "formalCharge", 0 ) or 0 ) for a in atoms ]

def _infer_variant ( residue ):
    """ Protonation variant implied by the hydrogens already present, or
        None when the residue has no hydrogens (library default). """
    names = { n for ( n, z, _ ) in residue.atoms if z == 1 }
    if not names: return None, "library default"
    name = residue.library_name
    if name == "HIS":
        d, e = "HD1" in names, "HE2" in names
        if d and e: return "Doubly Protonated", "hydrogens in file"
        if e:       return "Epsilon Protonated", "hydrogens in file"
        return "Delta Protonated", "hydrogens in file"
    if name == "ASP": return ( "Protonated" if "HD2" in names else "Deprotonated" ), "hydrogens in file"
    if name == "GLU": return ( "Protonated" if "HE2" in names else "Deprotonated" ), "hydrogens in file"
    if name == "LYS":
        n_hz = sum ( 1 for n in names if n.startswith ( "HZ" ) )
        return ( "Neutral" if n_hz == 2 else None ), "hydrogens in file"
    if name == "CYS":
        if residue.disulfide_partner is not None: return None, "disulfide"
        return ( None if "HG" in names else "Deprotonated" ), "hydrogens in file"
    return None, "library default"

_COVALENT_RADII = { 1: 0.31, 5: 0.84, 6: 0.76, 7: 0.71, 8: 0.66, 9: 0.57, 14: 1.11, 15: 1.07, 16: 1.05,
                    17: 1.02, 35: 1.20, 53: 1.39 }

def _molecules ( system, crd ):
    """ Connected molecules (lists of atoms) of a system without residue
        information: from its bonds, or by covalent radii (d < r1 + r2 +
        0.45 A) when it has none (e.g. an .xyz file). """
    atoms = list ( system.atoms )
    n = len ( atoms )
    parent = list ( range ( n ) )
    def find ( i ):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    connectivity = getattr ( system, "connectivity", None )
    if connectivity is not None and connectivity.bonds:
        index = { id ( a ): k for k, a in enumerate ( atoms ) }
        pairs = [ ( index[id ( b.node1 )], index[id ( b.node2 )] ) for b in connectivity.bonds ]
    else:
        if n > 5000:
            raise ValueError ( "A system without residues or bonds is limited to 5000 atoms." )
        xyz = np.array ( [ [ crd[i, 0], crd[i, 1], crd[i, 2] ] for i in range ( n ) ] )
        radii = np.array ( [ _COVALENT_RADII.get ( a.atomicNumber, 1.5 ) for a in atoms ] )
        d = np.linalg.norm ( xyz[:, None, :] - xyz[None, :, :], axis = 2 )
        limit = radii[:, None] + radii[None, :] + 0.45
        pairs = [ ( int ( i ), int ( j ) ) for i, j in zip ( *np.nonzero ( np.triu ( d < limit, 1 ) ) ) ]
    for ( i, j ) in pairs:
        a, b = find ( i ), find ( j )
        if a != b: parent[a] = b
    groups = OrderedDict ( )
    for k in range ( n ):
        groups.setdefault ( find ( k ), [ ] ).append ( atoms[k] )
    return list ( groups.values ( ) )

def _atom_groups ( system, crd ):
    """ [ ( entity label, residue name, number, insertion code, [ atoms ] ) ]:
        the sequence components, or -- for a system without residues (a
        molecule drawn in the Builder, an .xyz file) -- one group per
        connected molecule: an isolated water becomes HOH, an isolated
        Na/K/Cl atom an ion, anything else a ligand named LIG (or L01,
        L02, ... when there are several). """
    sequence = getattr ( system, "sequence", None )
    if sequence is not None and len ( sequence.children ) > 0:
        groups = [ ]
        for entity in sequence.children:
            for component in entity.children:
                name, number, icode = _residue_name_and_number ( component )
                groups.append ( ( entity.label, name, number, icode, list ( component.children ) ) )
        return groups
    molecules = _molecules ( system, crd )
    ligands = [ m for m in molecules if not _small_solvent ( m ) ]
    groups, k = [ ], 0
    for number, molecule in enumerate ( molecules, start = 1 ):
        solvent = _small_solvent ( molecule )
        if solvent:
            groups.append ( ( "A", solvent, number, "", molecule ) )
        else:
            k += 1
            groups.append ( ( "A", "LIG" if len ( ligands ) == 1 else "L{:02d}".format ( k ), number, "", molecule ) )
    return groups

def _small_solvent ( atoms ):
    """ "HOH" / "NA" / "K" / "CL" for an isolated water or monoatomic ion. """
    z = sorted ( a.atomicNumber for a in atoms )
    if z == [ 1, 1, 8 ]: return "HOH"
    if len ( z ) == 1 and z[0] in ( 11, 19, 17 ): return { 11: "NA", 19: "K", 17: "CL" }[z[0]]
    return None

_SOLVENT_ATOM_NAMES = { "HOH": [ "O", "H1", "H2" ], "NA": [ "NA" ], "K": [ "K" ], "CL": [ "CL" ] }

def _template_mismatch ( residue, component ):
    """ True when the residue has heavy atoms whose names the library
        component does not know (its variants only change hydrogens). """
    known = { a.label for a in component.atoms } | { "OXT" }
    return any ( z != 1 and name not in known for ( name, z, _ ) in residue.atoms )

def _merge_bonded_ligands ( system, analysis ):
    """ Fuses residues covalently bonded to a ligand (e.g. the "GLU" and
        "UNK" substructures of folate in a mol2 file) into ONE ligand
        residue. Amino acids peptide-bonded to other amino acids are never
        absorbed (a ligand bonded to a protein stays a covalent ligand,
        rejected later). Needs the input bonds (mol2, Builder). """
    connectivity = getattr ( system, "connectivity", None )
    residues = analysis.residues
    if connectivity is None or not connectivity.bonds or not any ( r.kind == "unsupported" for r in residues ):
        return
    owner = { }
    for k, r in enumerate ( residues ):
        for a in ( r.source_atoms or [ ] ): owner[id ( a )] = k
    def is_protein ( k ):
        """ a real amino acid: complete backbone (fragments of a mol2 "GLU"
            substructure are not) """
        r = residues[k]
        names = { n for ( n, z, _ ) in r.atoms }
        return r.kind == "amino acid" and { "N", "CA", "C", "O" } <= names
    links = [ ]
    for bond in connectivity.bonds:
        i, j = owner.get ( id ( bond.node1 ) ), owner.get ( id ( bond.node2 ) )
        if i is None or j is None or i == j: continue
        links.append ( ( i, j ) )
    parent = list ( range ( len ( residues ) ) )
    def find ( k ):
        while parent[k] != k:
            parent[k] = parent[parent[k]]; k = parent[k]
        return k
    for ( i, j ) in links:
        if residues[i].kind in ( "water", "ion" ) or residues[j].kind in ( "water", "ion" ): continue
        if is_protein ( i ) or is_protein ( j ): continue          # protein-ligand bond: covalent ligand
        a, b = find ( i ), find ( j )
        if a != b: parent[a] = b
    groups = OrderedDict ( )
    for k in range ( len ( residues ) ): groups.setdefault ( find ( k ), [ ] ).append ( k )
    merged, remove = { }, set ( )
    for members in groups.values ( ):
        if len ( members ) < 2 or not any ( residues[k].kind == "unsupported" for k in members ): continue
        members.sort ( )
        main = max ( ( k for k in members if residues[k].kind == "unsupported" ),
                     key = lambda k: sum ( 1 for a in residues[k].atoms if a[1] != 1 ) )
        r = residues[main]
        atoms = [ a for k in members for a in residues[k].source_atoms ]
        r.atoms = [ a for k in members for a in residues[k].atoms ]
        r.source_atoms = atoms
        r.kind, r.library_name, r.variant, r.variant_source = "unsupported", None, None, "library default"
        r.n_hydrogens = sum ( 1 for a in r.atoms if a[1] == 1 )
        r.bonds, r.formal_charges = _residue_chemistry ( system, atoms )
        r.notes = "merged covalently bonded fragments: " + ", ".join ( residues[k].label for k in members )
        remove.update ( k for k in members if k != main )
    if remove:
        analysis.residues = [ r for k, r in enumerate ( residues ) if k not in remove ]


def analyze_system ( system, coordinates3 = None, known_ligands = None, extra_library_paths = None ):
    """ Reads residues, chain breaks and disulfides of a loaded pDynamo
        system (with residues -- PDB, mol2 -- or without -- Builder, xyz:
        see _atom_groups()). The protonation variants can be edited in the
        returned analysis (residue.variant) before build_prepared_system(). """
    from pScientific import PeriodicTable
    analysis = PrepAnalysis ( )
    analysis.label = getattr ( system, "label", None )
    crd = system.coordinates3 if coordinates3 is None else coordinates3
    if len ( system.atoms ) == 0:
        raise ValueError ( "The system has no atoms." )

    for ( entity_label, name, number, icode, group_atoms ) in _atom_groups ( system, crd ):
            no_sequence = getattr ( system, "sequence", None ) is None or len ( system.sequence.children ) == 0
            if no_sequence and name in _SOLVENT_ATOM_NAMES:
                # . an isolated water/ion drawn without residue: library atom names
                ordered = sorted ( group_atoms, key = lambda a: -a.atomicNumber )
                renamed = dict ( zip ( [ id ( a ) for a in ordered ], _SOLVENT_ATOM_NAMES[name] ) )
            else:
                renamed = { }
            residue = PrepResidue ( )
            residue.entity, residue.name, residue.number, residue.icode = entity_label, name, number, icode
            residue.source_atoms = group_atoms
            library_name, alias_variant = RESIDUE_ALIASES.get ( name, ( name, None ) )
            library_component = _component ( library_name )
            residue.library_name = library_name if library_component is not None else None
            for atom in group_atoms:
                label = renamed.get ( id ( atom ), atom.label )
                atom_name = ATOM_ALIASES.get ( ( library_name, label ), ATOM_ALIASES.get ( ( "*", label ), label ) )
                xyz = None
                try:
                    xyz = np.array ( [ crd[atom.index, 0], crd[atom.index, 1], crd[atom.index, 2] ], dtype = float )
                except Exception:
                    pass
                residue.atoms.append ( ( atom_name, atom.atomicNumber, xyz ) )
            residue.n_hydrogens = sum ( 1 for a in residue.atoms if a[1] == 1 )
            if library_component is not None and _template_mismatch ( residue, library_component ):
                # . a residue name of the library with atoms the component does not
                #   have (mol2 substructures such as "GLU" in folate, Builder "UNK")
                residue.notes = "atoms do not match the {} template: treated as a ligand".format ( library_name )
                library_component, residue.library_name = None, None
            if library_component is None:
                residue.bonds, residue.formal_charges = _residue_chemistry ( system, group_atoms )
            if library_component is None and known_ligands and name in known_ligands:
                residue.library_name, residue.kind = name, "ligand"
            elif library_component is None:
                residue.kind = "unsupported"
            else:
                labels = { a.label for a in library_component.atoms }
                if { "N", "CA", "C" } <= labels:  residue.kind = "amino acid"
                elif library_name == "HOH":       residue.kind = "water"
                elif len ( labels ) == 1:         residue.kind = "ion"
                else:                             residue.kind = "other"
            if alias_variant is not None:
                residue.variant, residue.variant_source = alias_variant, "residue name {}".format ( name )
            analysis.residues.append ( residue )

    _merge_bonded_ligands ( system, analysis )
    residues = analysis.residues
    # . Disulfides (before variant inference: a bridged CYS has no HG).
    cys = [ i for i, r in enumerate ( residues ) if r.library_name == "CYS" and r.atom_xyz ( "SG" ) is not None ]
    for a in range ( len ( cys ) ):
        for b in range ( a + 1, len ( cys ) ):
            i, j = cys[a], cys[b]
            d = float ( np.linalg.norm ( residues[i].atom_xyz ( "SG" ) - residues[j].atom_xyz ( "SG" ) ) )
            if d <= DISULFIDE_SG_MAX and residues[i].disulfide_partner is None and residues[j].disulfide_partner is None:
                residues[i].disulfide_partner, residues[j].disulfide_partner = j, i
                analysis.disulfides.append ( ( i, j, d ) )
    for r in residues:
        if r.library_name in PROTONATION_CHOICES and r.variant_source == "library default":
            r.variant, r.variant_source = _infer_variant ( r )
        if r.disulfide_partner is not None:
            r.variant, r.variant_source = None, "disulfide"

    # . Segments: polymer pieces split at chain breaks; non-polymers apart.
    used = set ( )
    def new_segment ( base ):
        base = "".join ( c for c in base if c.isalnum ( ) )[:3] or "S"
        for k in range ( 0, 1000 ):
            label = base if k == 0 else "{}{}".format ( base, k )
            label = label[-4:] if len ( label ) > 4 else label
            if label not in used:
                used.add ( label )
                return label
        raise ValueError ( "Too many segments." )

    current_entity, current_segment, previous = None, None, None
    others = { }
    for i, r in enumerate ( residues ):
        if r.kind == "unsupported":
            analysis.unsupported.append ( i )
        if r.kind != "amino acid":
            key = ( r.entity, r.kind )
            if key not in others:
                tag = { "water": "W", "ion": "I" }.get ( r.kind, "X" )
                others[key] = new_segment ( "{}{}".format ( r.entity, tag ) )
            r.segment = others[key]
            current_entity = None                     # a non-polymer ends the polymer
            continue
        linked = False
        if current_entity == r.entity and previous is not None:
            c_prev, n_this = residues[previous].atom_xyz ( "C" ), r.atom_xyz ( "N" )
            if c_prev is not None and n_this is not None:
                d = float ( np.linalg.norm ( c_prev - n_this ) )
                if d <= PEPTIDE_BOND_MAX: linked = True
                else:                     analysis.chain_breaks.append ( ( previous, i, d ) )
        if not linked:
            current_segment = new_segment ( r.entity )
            analysis.polymers.append ( current_segment )
        r.segment = current_segment
        current_entity, previous = r.entity, i

    for i, r in enumerate ( residues ):
        analysis.segments.setdefault ( r.segment, [ ] ).append ( i )

    # . C-terminal OXT (often absent from PDB files): placed in the CA-C-O
    #   plane, 1.25 A from C, at 120 degrees from both CA and O.
    for segment in analysis.polymers:
        last = residues[analysis.segments[segment][-1]]
        if last.atom_xyz ( "OXT" ) is None:
            ca, c, o = last.atom_xyz ( "CA" ), last.atom_xyz ( "C" ), last.atom_xyz ( "O" )
            if ca is not None and c is not None and o is not None:
                u, v = ca - c, o - c
                u, v = u / np.linalg.norm ( u ), v / np.linalg.norm ( v )
                direction = - ( u + v )
                direction /= np.linalg.norm ( direction )
                last.atoms.append ( ( "OXT", 8, c + 1.25 * direction ) )
                last.added_atoms.append ( "OXT" )
    # . Heavy atoms the library expects but the file lacks (truncated side
    #   chains of crystal structures): reported -- see truncate_residue().
    for r in residues:
        r.missing_heavy = _missing_heavy_atoms ( r )
    if known_ligands:
        analysis.ligands.update ( known_ligands )
        analysis.extra_library_paths = list ( extra_library_paths or [ ] )
    return analysis


_TERMINAL_ONLY = { "OXT" }

def _missing_heavy_atoms ( residue ):
    if residue.library_name is None or residue.kind == "ligand": return [ ]
    component = _component ( residue.library_name )
    if component is None: return [ ]
    present = { name for ( name, z, xyz ) in residue.atoms if xyz is not None }
    expected = [ a.label for a in component.atoms if a.atomicNumber != 1 and a.label not in _TERMINAL_ONLY ]
    return [ name for name in expected if name not in present ]

def truncate_residue ( analysis, index ):
    """ Turns an amino acid with an incomplete side chain into ALA (CB
        present) or GLY, dropping the side-chain atoms beyond CB. A
        structural compromise the user must accept explicitly. """
    r = analysis.residues[index]
    if r.kind != "amino acid": raise ValueError ( "Only amino acids can be truncated." )
    names = { n for ( n, z, xyz ) in r.atoms if xyz is not None }
    target = "ALA" if "CB" in names else "GLY"
    keep = { "N", "CA", "C", "O", "OXT", "H", "H1", "H2", "H3", "HA", "HA2", "HA3" } | ( { "CB" } if target == "ALA" else set ( ) )
    r.atoms = [ a for a in r.atoms if a[0] in keep ]
    r.name, r.library_name, r.variant, r.variant_source = target, target, None, "truncated"
    r.missing_heavy = _missing_heavy_atoms ( r )
    return target


#=============================================================================
#  Ligands (phase 4)
#=============================================================================
COVALENT_CONTACT = 2.0      # A, ligand heavy atom closer than this to another residue = covalent

def parametrize_ligands ( analysis, work_dir, ph = 7.4, xtb = None, cm5_scale = None, progress = None ):
    """ Parametrizes every unsupported residue of the analysis
        (pdynamo/opls/ligand.py: Open Babel hydrogens at pH, OPLS classes by
        rules, GFN1-xTB CM5 charges) and makes it buildable: its PDB
        component goes to <work_dir>/pdb_components, the residue becomes a
        "ligand". One parametrization per residue NAME; every copy gets its
        own hydrogens. Raises ValueError for covalently bound residues. """
    from pdynamo.opls import ligand as L
    say = progress or ( lambda text: None )
    if not analysis.unsupported: return { }
    residues = analysis.residues
    heavy_others = [ ( i, xyz ) for i, r in enumerate ( residues ) for ( n, z, xyz ) in r.atoms
                     if z != 1 and xyz is not None and r.kind != "water" ]
    others_xyz = np.array ( [ x for ( _, x ) in heavy_others ] ) if heavy_others else np.zeros ( ( 0, 3 ) )
    others_id  = np.array ( [ i for ( i, _ ) in heavy_others ] ) if heavy_others else np.zeros ( 0, int )
    library_root = os.path.join ( work_dir, "pdb_components" )
    by_name = OrderedDict ( )
    for i in analysis.unsupported:
        by_name.setdefault ( residues[i].name, [ ] ).append ( i )
    for name, indices in by_name.items ( ):
        for i in indices:
            r = residues[i]
            mine = np.array ( [ xyz for ( n, z, xyz ) in r.atoms if z != 1 and xyz is not None ] )
            if len ( mine ) and len ( others_xyz ):
                mask = others_id != i
                d = np.linalg.norm ( others_xyz[mask][:, None, :] - mine[None, :, :], axis = 2 )
                if d.size and d.min ( ) < COVALENT_CONTACT and len ( r.atoms ) > 1:
                    raise ValueError ( "Residue {} is closer than {:.1f} A to another residue (covalent link?); "
                                       "covalent ligands are not supported.".format ( r.label, COVALENT_CONTACT ) )
        first = residues[indices[0]]
        say ( "Parametrizing {} (CM5 charges)...".format ( name ) )
        options = analysis.ligand_options.get ( name, { } )
        lig = L.parametrize_residue ( name, first.atoms, ph = ph, xtb = xtb, cm5_scale = cm5_scale,
                                      keep_scratch = os.path.join ( work_dir, "ligand_" + name ),
                                      bonds = first.bonds, formal_charges = first.formal_charges,
                                      total_charge = options.get ( "charge" ), multiplicity = options.get ( "multiplicity" ) )
        lig.write_pdb_component ( library_root )
        names = [ a[0] for a in lig.atoms ]
        for i in indices:
            r = residues[i]
            atoms = lig.atoms if i == indices[0] else L.ligand_chemistry ( name, r.atoms, ph = ph, bonds = r.bonds,
                                                                             formal_charges = r.formal_charges ).atoms
            if [ a[0] for a in atoms ] != names:
                raise ValueError ( "The copies of residue {} differ in atoms/protonation.".format ( name ) )
            r.atoms = [ ( a[0], a[1], np.asarray ( a[2], float ) ) for a in atoms ]
            r.library_name, r.kind, r.missing_heavy = name, "ligand", [ ]
            r.n_hydrogens = sum ( 1 for a in r.atoms if a[1] == 1 )
        analysis.ligands[name] = lig
    analysis.unsupported = [ ]
    if library_root not in analysis.extra_library_paths:
        analysis.extra_library_paths.insert ( 0, library_root )
    return dict ( analysis.ligands )


def ligand_preview ( analysis, ph = 7.4 ):
    """ { residue name: { "copies", "heavy", "hydrogens", "charge", "electrons",
          "multiplicity", "source" } } for the unsupported residues, from the
        chemistry step only (hydrogens/bond orders/formal charges, no QC) --
        what the GUI shows before the user edits charge and multiplicity.
        "charge" is None when the chemistry could not be perceived. """
    from pdynamo.opls import ligand as L
    preview = OrderedDict ( )
    for i in analysis.unsupported:
        r = analysis.residues[i]
        if r.name in preview:
            preview[r.name]["copies"] += 1
            continue
        entry = { "copies": 1, "heavy": sum ( 1 for a in r.atoms if a[1] != 1 ), "hydrogens": None,
                  "charge": None, "electrons": None, "multiplicity": 1, "source": None, "error": None }
        try:
            lig = L.ligand_chemistry ( r.name, r.atoms, ph = ph, bonds = r.bonds, formal_charges = r.formal_charges )
            entry["hydrogens"] = sum ( 1 for a in lig.atoms if a[1] == 1 )
            entry["charge"]    = lig.total_charge
            entry["electrons"] = lig.n_electrons
            entry["multiplicity"] = 1 if lig.n_electrons % 2 == 0 else 2
            entry["source"]    = lig.notes[0] if lig.notes else ""
        except Exception as error:
            entry["error"] = str ( error )
        preview[r.name] = entry
    return preview


#=============================================================================
#  Build
#=============================================================================
def _pdb_atom_name ( name, symbol ):
    """ PDB columns 13-16: one-letter elements start in column 14. """
    if len ( name ) >= 4 or len ( symbol ) == 2: return "{:<4s}".format ( name[:4] )
    return " {:<3s}".format ( name )

def write_clean_pdb ( analysis, path, residue_indices = None ):
    """ Writes the residues (all, or residue_indices) with one segment ID
        per segment and TER after every polymer segment. Residues are
        renumbered per segment when a number does not fit (> 9999). """
    from pScientific import PeriodicTable
    wanted = set ( range ( len ( analysis.residues ) ) ) if residue_indices is None else set ( residue_indices )
    lines, serial = [ ], 1
    analysis.pdb_labels = { }               # residue index -> ( segment, "NAME.number" )
    for segment, indices in analysis.segments.items ( ):
        indices = [ i for i in indices if i in wanted ]
        if not indices: continue
        renumber = any ( abs ( analysis.residues[i].number ) > 9999 for i in indices )
        for k, i in enumerate ( indices ):
            r = analysis.residues[i]
            number = ( k + 1 ) if renumber else r.number
            icode  = "" if renumber else r.icode
            analysis.pdb_labels[i] = ( segment, "{}.{}{}".format ( r.library_name, number, icode ) )
            record = "ATOM  " if r.kind == "amino acid" else "HETATM"
            for ( name, z, xyz ) in r.atoms:
                if xyz is None: continue
                symbol = PeriodicTable.Symbol ( z )
                lines.append ( "{:6s}{:5d} {:4s} {:>3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:<4s}{:>2s}\n".format (
                        record, serial % 100000, _pdb_atom_name ( name, symbol ), r.library_name[:3], " ",
                        number, icode[:1], xyz[0], xyz[1], xyz[2], 1.0, 0.0, segment, symbol.upper ( ) ) )
                serial += 1
        if segment in analysis.polymers:
            lines.append ( "TER\n" )
    lines.append ( "END\n" )
    with open ( path, "w" ) as handle:
        handle.writelines ( lines )
    return path

def _make_system ( model, label = None, altLoc = None ):
    """ PDBModel.MakeSystem() with the bonds handed to pDynamo as
        ( i, j, type ) integer tuples: Bond.FromIterable() tests
        `node in nodes` on a plain LIST for Bond objects, O(atoms) per
        bond -- 121 of 134 s for a 57 000-atom solvated protein. """
    from pMolecule import System
    ( sequence, atom_mapping ) = model.MakeSequence ( )
    bonds = model.MakeBonds ( atom_mapping ) or [ ]
    position = { id ( atom ): k for k, atom in enumerate ( sequence.GatherAtoms ( ) ) }
    tuples = [ ( position[id ( b.node1 )], position[id ( b.node2 )], b.type ) for b in bonds ]
    system = System.FromSequence ( sequence, bonds = tuples or None )
    system.coordinates3 = model.MakeCoordinates3 ( altLoc )
    if   label is not None:      system.label = label
    elif model.label is not None: system.label = model.label
    return system

def build_prepared_system ( analysis, label = None, random_seed = 1, keep_pdb = None ):
    """ Builds the prepared system. Returns ( system, report ):
          report = { "atoms", "added_hydrogens", "undefined_heavy" [paths],
                     "formal_charge", "variants" {label: variant},
                     "disulfides", "chain_breaks", "pdb" (path if kept) }
        Raises ValueError for unsupported residues (see analysis.unsupported). """
    from pBabel import PDBFileReader, PDBModelLink, PDBModelVariant
    from pSimulation import BuildHydrogenCoordinates3FromConnectivity
    from pScientific.RandomNumbers import RandomNumberGenerator

    if analysis.unsupported:
        names = sorted ( { analysis.residues[i].name for i in analysis.unsupported } )
        raise ValueError ( "Residues without a PDB component (ligands/cofactors): {}. "
                           "They must be parametrized separately.".format ( ", ".join ( names ) ) )
    incomplete = [ r.label for r in analysis.residues if r.missing_heavy ]
    if incomplete:
        raise ValueError ( "Residues with missing heavy atoms: {}. Complete them or truncate them "
                           "(truncate_residue()).".format ( ", ".join ( incomplete ) ) )
    with tempfile.TemporaryDirectory ( ) as scratch:
        path = write_clean_pdb ( analysis, os.path.join ( scratch, "prepared.pdb" ) )
        if keep_pdb:
            import shutil
            shutil.copy ( path, keep_pdb )
        raw   = PDBFileReader.PathToPDBModel ( path, log = None, useSegmentEntityLabels = True )
        model = PDBFileReader.PathToPDBModel ( path, log = None, useSegmentEntityLabels = True )

    # . Protonation variants and disulfide links on the model.
    def component_of ( i ):
        segment, component_label = analysis.pdb_labels[i]
        entity = model.childIndex.get ( segment )
        return entity, ( entity.childIndex.get ( component_label ) if entity is not None else None )
    variants = { }
    for i, r in enumerate ( analysis.residues ):
        if r.variant is None: continue
        entity, component = component_of ( i )
        if component is None: continue
        entity.AddVariant ( PDBModelVariant.WithOptions ( component = component, label = r.variant ) )
        variants[r.label] = r.variant
    for ( i, j, d ) in analysis.disulfides:
        _, left  = component_of ( i )
        _, right = component_of ( j )
        if left is not None and right is not None:
            model.AddLink ( PDBModelLink.WithOptions ( label = "Disulfide_Bridge", leftComponent = left, rightComponent = right ) )

    model.MakeAtomicModelFromComponentLibrary ( libraryPaths = list ( analysis.extra_library_paths ) + library_paths ( ), log = None )
    model.ExtractAtomData ( raw, log = None )
    system = _make_system ( model, label = label or analysis.label )
    undefined_before = system.coordinates3.numberUndefined
    rng = RandomNumberGenerator.WithSeed ( random_seed ) if random_seed is not None else RandomNumberGenerator.WithRandomSeed ( )
    BuildHydrogenCoordinates3FromConnectivity ( system, log = None, randomNumberGenerator = rng )

    undefined_heavy = [ ]
    if system.coordinates3.numberUndefined > 0:
        rows = set ( getattr ( system.coordinates3, "undefined", None ) or [ ] )
        undefined_heavy = [ a.path for a in system.atoms if a.index in rows ]
    report = { "atoms"           : len ( system.atoms ),
               "input_atoms"     : sum ( len ( r.atoms ) for r in analysis.residues ),
               "added_hydrogens" : undefined_before - len ( undefined_heavy ),
               "undefined_atoms" : undefined_heavy,
               "formal_charge"   : int ( sum ( getattr ( a, "formalCharge", 0 ) or 0 for a in system.atoms ) ),
               "variants"        : variants,
               "disulfides"      : len ( analysis.disulfides ),
               "chain_breaks"    : len ( analysis.chain_breaks ),
               "segments"        : list ( analysis.segments.keys ( ) ) }
    return system, report


#=============================================================================
#  Typing
#=============================================================================
def parameter_set_path ( parameter_set = DEFAULT_PARAMETER_SET ):
    """ A parameter set name (pDynamo's opls/<name>) or a folder path. """
    if os.path.isdir ( parameter_set ):
        return os.path.abspath ( parameter_set )
    return os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "forceFields", "opls", parameter_set )

def type_atoms ( system, parameter_set = DEFAULT_PARAMETER_SET ):
    """ ( types [label or None per atom], charges, untyped [Atom] ) without
        building an MM model. """
    from pMolecule.MMModel.MMAtomTyper import MMAtomTyper
    from pMolecule.MMModel.MMModelError import MMModelError
    typer = MMAtomTyper.FromPath ( parameter_set_path ( parameter_set ) )
    types   = [ None ] * len ( system.atoms )
    charges = [ 0.0  ] * len ( system.atoms )
    untyped = set ( system.connectivity.atoms )
    typer.TypeBySequence ( system.sequence, types, charges, untyped )
    typer.TypeByPattern  ( system.connectivity, types, charges, untyped )
    return types, charges, sorted ( untyped, key = lambda a: a.index )

def typing_report ( system, parameter_set = DEFAULT_PARAMETER_SET ):
    """ { "ok", "n_atoms", "untyped" [paths], "untyped_by_residue" {residue: [atom names]},
          "types" Counter, "charge" (sum of typed charges), "formal_charge" } """
    types, charges, untyped = type_atoms ( system, parameter_set )
    by_residue = OrderedDict ( )
    for atom in untyped:
        by_residue.setdefault ( atom.parent.path, [ ] ).append ( atom.label )
    return { "ok"                 : len ( untyped ) == 0,
             "n_atoms"            : len ( types ),
             "untyped"            : [ a.path for a in untyped ],
             "untyped_by_residue" : by_residue,
             "types"              : Counter ( t for t in types if t is not None ),
             "charge"             : float ( sum ( charges ) ),
             "formal_charge"      : int ( sum ( getattr ( a, "formalCharge", 0 ) or 0 for a in system.atoms ) ) }


#=============================================================================
#  End to end
#=============================================================================
def prepare_opls_system ( system, parameter_set = DEFAULT_PARAMETER_SET, work_dir = None, label = None,
                          analysis = None, define_nb_model = True, solvation = None, progress = None,
                          parametrize_unsupported = True, ligand_ph = 7.4, cm5_scale = None, xtb = None ):
    """ Loaded system -> prepared, typed system with an OPLS MM model.

        analysis: an edited analyze_system() result (protonation choices);
        solvation: None, or the options of solvate.solvate_coordinates()
        (shape, padding, box_size, neutralize, concentration, cation, anion, ...);
        progress: optional callable ( text ) for status messages;
        work_dir: where the per-system parameter set is written (default:
        a new temporary folder). Returns ( system, report ) where report =
        build report + "typing" + "parameters" (materialize report). The
        MM model is NOT defined when atoms stay untyped or parameters are
        missing; report["ok"] tells. """
    from pdynamo.opls import parameters
    from pMolecule.MMModel import MMModelOPLS
    from pMolecule.NBModel import NBModelCutOff
    say = progress or ( lambda text: None )
    if work_dir is None:
        work_dir = tempfile.mkdtemp ( prefix = "opls_set_" )
    os.makedirs ( work_dir, exist_ok = True )
    if analysis is None:
        say ( "Analyzing residues..." )
        analysis = analyze_system ( system )
    if analysis.unsupported and parametrize_unsupported:
        parametrize_ligands ( analysis, work_dir, ph = ligand_ph, cm5_scale = cm5_scale, progress = say, xtb = xtb )
    say ( "Building hydrogens, termini and links..." )
    new_system, report = build_prepared_system ( analysis, label = label )
    # . ligand total charges set by the user (not carried by the atoms' formal charges)
    correction = 0
    for r in analysis.residues:
        lig = analysis.ligands.get ( r.library_name ) if r.kind == "ligand" else None
        if lig is not None and len ( lig.atoms ) > 1:
            correction += lig.total_charge - lig.formal_charge_sum
    report["total_charge"] = report["formal_charge"] + correction
    if solvation is not None and not report["undefined_atoms"]:
        from pdynamo.opls import solvate
        say ( "Solvating..." )
        new_system, solvation_report = solvate.solvate_prepared_system ( new_system, label = label, ligands = dict ( analysis.ligands ),
                                                                         extra_library_paths = analysis.extra_library_paths,
                                                                         charge_correction = correction, **solvation )
        report["solvation"] = solvation_report
        report["atoms"] = solvation_report["atoms"]
        report["formal_charge"] = solvation_report["formal_charge"]
        report["total_charge"]  = solvation_report["formal_charge"] + correction
    if analysis.ligands:
        from pdynamo.opls import parameters as _parameters
        parameter_set = _parameters.extend_parameter_set ( os.path.join ( work_dir, "base_with_ligands" ),
                                                           base = parameter_set, ligands = analysis.ligands.values ( ) )
        report["ligands"] = { name: lig.summary ( ) for name, lig in analysis.ligands.items ( ) }
    say ( "Typing atoms..." )
    typing = typing_report ( new_system, parameter_set )
    report["typing"] = typing
    report["ok"] = typing["ok"] and not report["undefined_atoms"]
    if not report["ok"]:
        return new_system, report
    types, charges, _ = type_atoms ( new_system, parameter_set )
    say ( "Writing the parameter set..." )
    database = parameters.build_class_database ( parameter_set )
    fallback = None
    if analysis.ligands:
        from pdynamo.opls import ligand as L
        fallback, notes = L.estimate_missing_terms ( new_system, types, database,
                                                     parameters.PDynamoSet ( parameters.pdynamo_set_path ( parameter_set ) ) )
        report["estimated_terms"] = notes
    set_dir = os.path.join ( work_dir, "parameters" )
    materialized = parameters.materialize_parameter_set ( new_system, types, set_dir, base = parameter_set,
                                                          database = database, fallback_terms = fallback )
    report["parameters"] = materialized
    if any ( materialized["missing"].values ( ) ):
        report["ok"] = False
        return new_system, report
    if analysis.ligands:
        from pdynamo.opls import quality
        try:
            report["quality"] = quality.assess ( new_system, types, charges, materialized, set ( analysis.ligands ) )
            import json
            with open ( os.path.join ( work_dir, "quality.json" ), "w" ) as handle:
                json.dump ( report["quality"], handle, indent = 1, default = str )
        except Exception as error:
            report["quality_error"] = str ( error )
    say ( "Building the MM model..." )
    new_system.DefineMMModel ( MMModelOPLS.WithParameterSet ( os.path.abspath ( set_dir ) ), log = None )
    if define_nb_model:
        new_system.DefineNBModel ( NBModelCutOff.WithDefaults ( ) )
    report["parameter_folder"] = os.path.abspath ( set_dir )
    return new_system, report


#=============================================================================
#  Report
#=============================================================================
def report_markdown ( report, analysis = None, title = "OPLS system preparation" ):
    """ Human-readable report (Markdown) of prepare_opls_system(). """
    import datetime
    lines = [ "# {}".format ( title ), "", "Generated by EasyHybrid (pdynamo/opls), {}.".format (
              datetime.datetime.now ( ).strftime ( "%Y-%m-%d %H:%M" ) ), "" ]
    lines += [ "## Result", "",
               "- Status: **{}**".format ( "ready (OPLS MM model defined)" if report.get ( "ok" ) else "FAILED" ),
               "- Atoms: {} (input {}, hydrogens built {})".format ( report.get ( "atoms" ), report.get ( "input_atoms" ), report.get ( "added_hydrogens" ) ),
               "- Total charge: {:+d}{}".format ( report.get ( "total_charge", report.get ( "formal_charge", 0 ) ),
                    "" if report.get ( "total_charge", report.get ( "formal_charge" ) ) == report.get ( "formal_charge" )
                    else " (formal charges of the atoms: {:+d}; ligand charges set by the user)".format ( report.get ( "formal_charge", 0 ) ) ),
               "- Segments: {}".format ( ", ".join ( report.get ( "segments", [ ] ) ) ),
               "- Disulfide bridges: {}; chain breaks: {}".format ( report.get ( "disulfides", 0 ), report.get ( "chain_breaks", 0 ) ) ]
    if report.get ( "parameter_folder" ):
        lines.append ( "- Parameter set: `{}`".format ( report["parameter_folder"] ) )
    if analysis is not None:
        summary = analysis.summary ( )
        changed = { k: v for k, v in summary["variants"].items ( ) if v[1] != "library default" or v[0] }
        lines += [ "", "## Protonation states", "", "Library defaults (pH 7: ASP/GLU deprotonated, LYS/ARG protonated, HIS delta-protonated) except:", "" ]
        lines += [ "- {}: {} ({})".format ( k, v[0] or "library form", v[1] ) for k, v in changed.items ( ) ] or [ "- none" ]
        if summary["added_atoms"]:
            lines += [ "", "Atoms placed by geometry: " + ", ".join ( "{} {}".format ( k, "/".join ( v ) ) for k, v in summary["added_atoms"].items ( ) ) ]
        truncated = [ r.label for r in analysis.residues if r.variant_source == "truncated" ]
        if truncated:
            lines += [ "", "Truncated residues (incomplete side chains): " + ", ".join ( truncated ) ]
    solvation = report.get ( "solvation" )
    if solvation:
        lines += [ "", "## Solvation", "",
                   "- Shape: {}".format ( solvation["shape"] ),
                   "- Cell: {}".format ( "{:.2f} x {:.2f} x {:.2f} A".format ( *solvation["cell"][:3] ) if solvation.get ( "cell" ) else "none (droplet, radius {:.2f} A)".format ( solvation["sphere"][1] ) ),
                   "- Waters (TIP3P): {}; removed by the solute: {}".format ( solvation["waters"], solvation["removed_by_solute"] ),
                   "- Ions: {} (solute charge {:+d})".format ( ", ".join ( "{} {}".format ( n, k ) for k, n in solvation["ions"].items ( ) ) or "none", solvation["solute_charge"] ) ]
    ligands = report.get ( "ligands" )
    if ligands:
        lines += [ "", "## Ligands (simplified parametrization)", "" ]
        for name, data in ligands.items ( ):
            lines.append ( "- **{}**: {} atoms, charge {:+d}, multiplicity {}, charges: {}, classes: {}".format (
                    name, data["atoms"], data["total_charge"], data.get ( "multiplicity", 1 ), data["charge_method"],
                    ", ".join ( "{} {}".format ( c, n ) for c, n in data["classes"].items ( ) ) ) )
            for note in data["notes"] + data["warnings"]:
                lines.append ( "  - {}".format ( note ) )
        estimated = report.get ( "estimated_terms", [ ] )
        lines += [ "", "### Estimated bonded terms ({})".format ( len ( estimated ) ), "" ]
        lines += [ "- " + e for e in estimated ] or [ "None: every term came from the OPLS database." ]
    parameters = report.get ( "parameters" )
    if parameters:
        def friendly ( source ):
            if source == "base":      return "pDynamo set, exact labels"
            if source.startswith ( "pdynamo:" ): return "pDynamo {} set, by class".format ( source.split ( ":", 1 )[1] )
            return { "oplsaal": "OPLS-AA/L (oplsaal.prm)", "oplsaa08": "OPLS-AA 2008 (oplsaa08.prm)",
                     "estimated": "estimated (see above)", "zero": "zero (no OPLS improper)" }.get ( source, source )
        lines += [ "", "## Parameter sources", "", "| Term | Sources |", "|---|---|" ]
        for term, counts in parameters["sources"].items ( ):
            if counts: lines.append ( "| {} | {} |".format ( term, ", ".join ( "{} {}".format ( friendly ( k ), v ) for k, v in counts.items ( ) ) ) )
        if parameters.get ( "zero_outofplane" ):
            lines.append ( "" )
            lines.append ( "Out-of-plane terms without parameters (K = 0): {}".format ( len ( parameters["zero_outofplane"] ) ) )
        missing = { t: v for t, v in parameters["missing"].items ( ) if v }
        if missing:
            lines += [ "", "**Missing parameters:**", "" ] + [ "- {}: {}".format ( t, ", ".join ( "-".join ( k ) for k in v ) ) for t, v in missing.items ( ) ]
    if report.get ( "quality" ):
        from pdynamo.opls import quality
        lines += quality.markdown ( report["quality"] )
    typing = report.get ( "typing" )
    if typing and not typing["ok"]:
        lines += [ "", "## Untyped atoms", "" ] + [ "- {}: {}".format ( r, ", ".join ( a ) ) for r, a in typing["untyped_by_residue"].items ( ) ]
    lines += [ "", "## Before production", "",
               "Minimize, then equilibrate (NVT, then NPT for periodic boxes). Ligand parameters are a simplified",
               "OPLS-AA/CM5 model: check the estimated terms above." ]
    return "\n".join ( lines ) + "\n"
