#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: DYFF torsion corrections looked up by chemical environment
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description (2026-10-07):
#      A per-bond DYFF torsion correction ( dyff_edits "dihedral_bonds" ) taken
#      from a library of B3LYP-D3BJ/DZVP torsion scans instead of a GFN2-xTB
#      scan of the user's molecule. The library (data/dyff_bond_library.json,
#      made by fraglib/fraglib/export_bond_library.py) holds, per chemical
#      environment code of the bond, the cosine series fitted on top of DYFF +
#      modified Seminario bonds/angles (from the GFN2-xTB Hessian) -- the same
#      base the DYFF Parameters window's refinement builds, so a library
#      series is used together with the molecule's own Hessian bonds/angles.
#
#      Lookup: exact (same fragment) -> r2 -> r1 environment of the bond's two
#      atoms (ECFP-like codes, fraglib/fraglib/environment.py; the functions
#      below are a copy and MUST give identical codes, tests/test_bond_library.py).
#      Tested against DFT on 669 held-out scans (fraglib/reports/fraglib_status.md
#      section 9): median RMSE 1.29 kcal/mol where an r2 code matched (own xTB
#      scan 1.27), 1.43 at r1 (own xTB scan 1.59).
#
#      Two tables (2026-10-08): data/dyff_bond_library.json, fitted WITHOUT partial
#      charges, and data/dyff_bond_library_cm5.json, fitted with GFN1-xTB CM5 charges
#      and the polar-H LJ type (dyff_edits.soften_polar_hydrogens) -- the model a
#      molecule has once it gets CM5 or residue-library charges. A correction must
#      be used with the model it was fitted on: on the charged model the uncharged
#      table gives 1.76 kcal/mol (median vs DFT, 588 held-out scans), the CM5 one
#      1.62 (fraglib/reports/fraglib_status.md section 10). load_library ( charged ).
#
#      Codes without RDKit (2026-10-14): the codes are computed by pdynamo/env_codes.py,
#      a pure-Python port of the RDKit 2026.03.6 perception the libraries were made
#      with (r1 / r2 identical to RDKit on 18 413 molecules, DYFF_validation/scripts/
#      env_codes_check.py). Its exact level has no stereo, so the library files were
#      re-keyed (fraglib/fraglib/rekey_env_codes.py, env_code_version "2"); the files
#      keyed with RDKit are kept as *_envcode1.json. EASYHYBRID_ENV_CODES=rdkit restores
#      the RDKit codes and those files (the results produced before 2026-10-14).
#

import hashlib
import json
import os

RDKIT_CODES = os.environ.get ( "EASYHYBRID_ENV_CODES", "" ).lower ( ) == "rdkit"
_SUFFIX = "_envcode1" if RDKIT_CODES else ""
_DATA = os.path.join ( os.path.dirname ( os.path.abspath ( __file__ ) ), "data" )
HESSIAN_PATH = os.path.join ( _DATA, "dyff_hessian_library{}.json".format ( _SUFFIX ) )
DEFAULT_PATH = os.path.join ( _DATA, "dyff_bond_library{}.json".format ( _SUFFIX ) )
CHARGED_PATH = os.path.join ( _DATA, "dyff_bond_library_cm5{}.json".format ( _SUFFIX ) )
ENV_CODE_VERSION = "1" if RDKIT_CODES else "2"
LEVELS = ( "exact", "r2", "r1" )

_cache = { }

def load_library ( path = None, charged = False ):
    """ The library dict (the CM5 table when charged), or None when the file is
        missing / of another code version. """
    path = os.path.abspath ( path or ( CHARGED_PATH if charged else DEFAULT_PATH ) )
    if path not in _cache:
        library = None
        if os.path.exists ( path ):
            with open ( path ) as handle: library = json.load ( handle )
            if str ( library.get ( "env_code_version" ) ) != ENV_CODE_VERSION: library = None
        _cache[path] = library
    return _cache[path]


#=============================================================================
#  Environment codes (copy of fraglib/fraglib/environment.py, version "1")
#=============================================================================
def _sha ( text, n = 16 ):
    return hashlib.sha1 ( text.encode ( ) ).hexdigest ( )[:n]

def _bond_order_symbol ( bond ):
    from rdkit import Chem
    return { Chem.BondType.SINGLE: "1", Chem.BondType.DOUBLE: "2", Chem.BondType.TRIPLE: "3",
             Chem.BondType.AROMATIC: "a" }.get ( bond.GetBondType ( ), "?" )

def _atom_invariant ( atom ):
    return "{}|{}|{}|{}|{}|{}|{}".format ( atom.GetAtomicNum ( ), int ( atom.GetIsAromatic ( ) ), atom.GetDegree ( ),
                                           atom.GetTotalNumHs ( ), atom.GetFormalCharge ( ), int ( atom.IsInRing ( ) ),
                                           str ( atom.GetHybridization ( ) ) )

def atom_identifiers ( mol, radius ):
    ids = [ _sha ( _atom_invariant ( a ), 12 ) for a in mol.GetAtoms ( ) ]
    for _ in range ( radius ):
        new = [ ]
        for atom in mol.GetAtoms ( ):
            neighbours = sorted ( "{}{}".format ( _bond_order_symbol ( b ), ids[b.GetOtherAtomIdx ( atom.GetIdx ( ) )] )
                                  for b in atom.GetBonds ( ) )
            new.append ( _sha ( ids[atom.GetIdx ( )] + "(" + ",".join ( neighbours ) + ")", 12 ) )
        ids = new
    return ids

def term_codes ( mol, atoms, kind = "bond", charge = 0, multiplicity = 1, identifiers = None, smiles = None, ranks = None ):
    """ { "exact", "r1", "r2" } of a chain term (kind "bond": 2 atoms, "angle": 3)
        -- fraglib term_codes. identifiers / smiles / ranks: precomputed for the
        molecule (speed). """
    from rdkit import Chem
    atoms = [ int ( a ) for a in atoms ]
    if smiles is None: smiles = Chem.MolToSmiles ( Chem.RemoveHs ( mol ), isomericSmiles = True )
    if ranks is None:  ranks = list ( Chem.CanonicalRankAtoms ( mol, breakTies = False ) )
    classes = [ ranks[a] for a in atoms ]
    classes = min ( classes, classes[::-1] )
    codes = { "exact": _sha ( "{}|{}|{}|{}|{}".format ( smiles, kind, "-".join ( map ( str, classes ) ), charge, multiplicity ) ) }
    orders = [ ]
    for a, b in zip ( atoms, atoms[1:] ):
        bond = mol.GetBondBetweenAtoms ( a, b )
        orders.append ( _bond_order_symbol ( bond ) if bond is not None else "-" )
    for radius in ( 1, 2 ):
        ids = identifiers[radius] if identifiers else atom_identifiers ( mol, radius )
        forward = [ ids[a] for a in atoms ]
        key = min ( forward, forward[::-1] )
        orders_k = orders[::-1] if key != forward else orders
        codes["r{}".format ( radius )] = "{}|r{}|{}".format ( kind, radius, _sha ( "{}|{}|{}".format ( kind, ",".join ( key ), "".join ( orders_k ) ) ) )
    return codes

def bond_codes ( mol, j, k, charge = 0, multiplicity = 1, identifiers = None, smiles = None, ranks = None ):
    """ { "exact", "r1", "r2" } of the bond j-k (fraglib term_codes, kind "bond"). """
    return term_codes ( mol, ( j, k ), "bond", charge, multiplicity, identifiers, smiles, ranks )


#=============================================================================
#  Molecule -> RDKit
#=============================================================================
def rdkit_molecule ( atomic_numbers, bonds, formal_charges = None, coordinates = None ):
    """ Sanitized RDKit molecule (all hydrogens explicit) from atomic numbers,
        bonds ( a, b, order, aromatic ) (Kekule orders) and formal charges, in
        the given atom order. Aromaticity, hybridization and rings are then
        perceived by RDKit, as for the library's molecules; with coordinates
        (A), stereo centres and double bonds are assigned from 3D (the library's
        isomeric SMILES carry them: the "exact" code needs them). """
    from rdkit import Chem
    types = { 1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE }
    formal = list ( formal_charges ) if formal_charges is not None else [ 0 ] * len ( atomic_numbers )
    def build ( kekule ):
        mol = Chem.RWMol ( )
        for z, q in zip ( atomic_numbers, formal ):
            atom = Chem.Atom ( int ( z ) )
            atom.SetFormalCharge ( int ( q ) ); atom.SetNoImplicit ( True )
            mol.AddAtom ( atom )
        for a, b, order, aromatic in bonds:
            if kekule or not aromatic:
                mol.AddBond ( int ( a ), int ( b ), types.get ( int ( order ), Chem.BondType.SINGLE ) )
            else:
                mol.AddBond ( int ( a ), int ( b ), Chem.BondType.AROMATIC )
                mol.GetAtomWithIdx ( int ( a ) ).SetIsAromatic ( True ); mol.GetAtomWithIdx ( int ( b ) ).SetIsAromatic ( True )
        mol = mol.GetMol ( )
        Chem.SanitizeMol ( mol )
        if coordinates is not None:
            from rdkit.Geometry import Point3D
            conformer = Chem.Conformer ( mol.GetNumAtoms ( ) )
            for i, x in enumerate ( coordinates ):
                conformer.SetAtomPosition ( i, Point3D ( float ( x[0] ), float ( x[1] ), float ( x[2] ) ) )
            mol.AddConformer ( conformer, assignId = True )
            Chem.AssignStereochemistryFrom3D ( mol )
        return mol
    try:
        return build ( True )
    except Exception:
        return build ( False )                 # aromatic flags without a consistent Kekule form


#=============================================================================
#  Lookup
#=============================================================================
def lookup ( library, codes, levels = LEVELS ):
    """ ( level, pairs, n scans ) of the first level whose code is in the
        library, or ( None, None, 0 ). A code may be a list of candidates. """
    for level in levels:
        candidates = codes[level] if isinstance ( codes[level], ( list, tuple ) ) else [ codes[level] ]
        for code in candidates:
            entry = library["levels"].get ( level, { } ).get ( code )
            if entry and entry["pairs"]:
                return level, [ ( float ( c ), int ( n ) ) for c, n in entry["pairs"] ], int ( entry["n"] )
    return None, None, 0

def code_molecule ( atomic_numbers, bonds, formal_charges = None ):
    """ pdynamo.env_codes.Molecule from atomic numbers, bonds ( a, b, order, aromatic )
        (Kekule orders, all hydrogens explicit) and formal charges. """
    from pdynamo import env_codes as E
    orders = { 1: 1, 2: 2, 3: 3 }
    return E.Molecule ( atomic_numbers, [ ( int ( b[0] ), int ( b[1] ), orders.get ( int ( round ( float ( b[2] ) ) ), 1 ) ) for b in bonds ],
                        formal_charges )


def acyclic_single ( mol, a, b ):
    """ True for a non-aromatic single bond outside rings (mol: code_molecule, or an
        RDKit molecule with EASYHYBRID_ENV_CODES=rdkit). """
    if RDKIT_CODES:
        bond = mol.GetBondBetweenAtoms ( int ( a ), int ( b ) )
        return bond is not None and str ( bond.GetBondType ( ) ) == "SINGLE" and not bond.IsInRing ( )
    k = mol.bond_of.get ( ( int ( a ), int ( b ) ) )
    return k is not None and mol.bond_symbol ( k ) == "1" and k not in mol.ring_bonds


def library_torsions ( atomic_numbers, bonds, formal_charges, targets, coordinates = None, library = None, levels = LEVELS ):
    """ { ( j, k ): ( level, pairs, n ) } for the target bonds ( j, k ) found in
        the library (bonds ( a, b, Kekule order, aromatic ), formal charges; the
        coordinates are used only by the RDKit codes, for stereo). The exact codes
        of the library were made with charge 0 and singlet (fraglib fit_dyff.
        bond_codes): same here, the fragment already carries the formal charges. """
    library = library if library is not None else load_library ( )
    if library is None or not targets: return { }
    if not RDKIT_CODES:
        from pdynamo import env_codes as E
        mol = code_molecule ( atomic_numbers, bonds, formal_charges )
        identifiers = { r: E.atom_identifiers ( mol, r ) for r in ( 1, 2 ) }
        labels = E._refined_labels ( mol ); graph_id = E.graph_identifier ( mol, labels )
        found = { }
        for j, k in targets:
            codes = E.term_codes_r ( mol, ( j, k ), "bond", identifiers )
            codes["exact"] = E.exact_code ( mol, ( j, k ), "bond", 0, 1, labels, graph_id )
            level, pairs, n = lookup ( library, codes, levels )
            if pairs is not None: found[( j, k )] = ( level, pairs, n )
        return found
    from rdkit import Chem
    mol = rdkit_molecule ( atomic_numbers, bonds, formal_charges, coordinates )
    identifiers = { r: atom_identifiers ( mol, r ) for r in ( 1, 2 ) }
    ranks = list ( Chem.CanonicalRankAtoms ( mol, breakTies = False ) )
    # . exact: the library's SMILES carry the stereo its source specified, or
    #   none: try the 3D-assigned stereo, then none (r1/r2 ignore stereo)
    flat = Chem.Mol ( mol ); Chem.RemoveStereochemistry ( flat )
    variants = list ( dict.fromkeys ( Chem.MolToSmiles ( Chem.RemoveHs ( m ), isomericSmiles = True ) for m in ( mol, flat ) ) )
    found = { }
    for j, k in targets:
        codes = bond_codes ( mol, j, k, identifiers = identifiers, smiles = variants[0], ranks = ranks )
        codes["exact"] = [ codes["exact"] ] + [ bond_codes ( mol, j, k, identifiers = identifiers, smiles = v, ranks = ranks )["exact"]
                                                for v in variants[1:] ]
        level, pairs, n = lookup ( library, codes, levels )
        if pairs is not None: found[( j, k )] = ( level, pairs, n )
    return found


#=============================================================================
#  Calibration of the xTB Hessian bonds/angles by environment (2026-10-08)
#=============================================================================
#  data/dyff_hessian_library.json (fraglib/fraglib/run_hessian_library.py --export):
#  per environment code of a bond / angle, the mean DFT - xTB shift of the
#  equilibrium value and the mean k_DFT / k_xTB ratio of the modified-Seminario
#  constants (B3LYP-D3BJ/DZVP vs GFN2-xTB Hessians of 928 QCArchive molecules).
#  Applied to the xTB constants of the user's own molecule (seminario.
#  dyff_edits_from_hessian calibration=...): on 309 held-out molecules the
#  frequency RMSE vs DFT went 116.5 -> 53.9 cm-1 (DFT Hessian itself: 51.3),
#  bond MAE 0.0086 -> 0.0079 A, angle MAE 1.12 -> 0.86 degrees (fraglib/reports/
#  fraglib_status.md section 11). Replacing by the DFT mean instead was as good
#  for frequencies but worse for angles: the calibration keeps the molecule's
#  own geometry and corrects only the method's bias.

def load_hessian_library ( path = None ):
    path = os.path.abspath ( path or HESSIAN_PATH )
    if path not in _cache:
        table = None
        if os.path.exists ( path ):
            with open ( path ) as handle: table = json.load ( handle )
            if str ( table.get ( "env_code_version" ) ) != ENV_CODE_VERSION: table = None
        _cache[path] = table
    return _cache[path]


class HessianCalibration:
    """ calibration ( kind, atoms ) -> ( shift, ratio, level ) or None, for a molecule
        given as in library_torsions; kind "bond" ( a, b ) or "angle" ( a, b, c ). """

    def __init__ ( self, atomic_numbers, bonds, formal_charges = None, coordinates = None, table = None ):
        self.table = table if table is not None else load_hessian_library ( )
        self.used = { }
        if not RDKIT_CODES:
            from pdynamo import env_codes as E
            self.mol = code_molecule ( atomic_numbers, bonds, formal_charges )
            self.identifiers = { r: E.atom_identifiers ( self.mol, r ) for r in ( 1, 2 ) }
            self.labels = E._refined_labels ( self.mol ); self.graph_id = E.graph_identifier ( self.mol, self.labels )
            return
        from rdkit import Chem
        self.mol = rdkit_molecule ( atomic_numbers, bonds, formal_charges, coordinates )
        self.identifiers = { r: atom_identifiers ( self.mol, r ) for r in ( 1, 2 ) }
        self.ranks = list ( Chem.CanonicalRankAtoms ( self.mol, breakTies = False ) )
        flat = Chem.Mol ( self.mol ); Chem.RemoveStereochemistry ( flat )
        self.smiles = list ( dict.fromkeys ( Chem.MolToSmiles ( Chem.RemoveHs ( m ), isomericSmiles = True ) for m in ( self.mol, flat ) ) )

    def _codes ( self, kind, atoms ):
        if not RDKIT_CODES:
            from pdynamo import env_codes as E
            c = E.term_codes_r ( self.mol, atoms, kind, self.identifiers )
            c["exact"] = E.exact_code ( self.mol, atoms, kind, 0, 1, self.labels, self.graph_id )
            return [ c ]
        return [ term_codes ( self.mol, atoms, kind, identifiers = self.identifiers, smiles = s, ranks = self.ranks ) for s in self.smiles ]

    def __call__ ( self, kind, atoms ):
        if self.table is None: return None
        levels = self.table["table"].get ( kind, { } )
        codes = self._codes ( kind, atoms )
        for level in LEVELS:
            for c in ( codes if level == "exact" else codes[:1] ):
                entry = levels.get ( level, { } ).get ( c[level] )
                if entry:
                    self.used[level] = self.used.get ( level, 0 ) + 1
                    return float ( entry[0] ), float ( entry[1] ), level
        return None

    @classmethod
    def mapped ( cls, calibration, parent_index ):
        """ A calibration of the parent molecule used with subsystem indices (caps: none). """
        def f ( kind, atoms ):
            idx = [ parent_index[int ( a )] for a in atoms ]
            return None if any ( i < 0 for i in idx ) else calibration ( kind, idx )
        f.used = calibration.used
        return f
