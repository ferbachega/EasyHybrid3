#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: DYFF with user-edited parameters
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      DYFF generates every bonded parameter with "factories" (pMolecule/
#      MMModel/DYFFUtilities.py) that return ( value, key ), the key being
#      the atom types plus the bond tags (e.g. "C:Res:C:Res:C:Res:C:Res:S").
#      Its parameter files are empty, and pDynamo's reader of explicit
#      cosine rows has a bug, so edits cannot go through files. Instead,
#      MMModelDYFFEdited wraps the factories: when the generated key has a
#      user edit, the edited value is returned instead of the rule-based one.
#
#      edits = { "bond": { key: ( b0, k ) }, "angle" | "dihedral" | "outofplane":
#                { key: [ ( force constant, period ), ... ] } }   (kJ/mol units)
#
#      Per-bond dihedral edits (2026-10-06, fragment library refinement):
#
#        edits["dihedral_bonds"] = { "j-k": [ ( c, n ), ... ] }   (j < k, atom indices)
#
#      replace EVERY dihedral term whose central atoms are j and k (whatever
#      its types) by the given series; an empty list removes them. They take
#      precedence over key edits. Each such term returns its own key
#      ( key + ( "bond j-k", ) ): pDynamo shares parameters between terms
#      with the same returned key, so a per-bond value must not reuse the
#      type key of other terms.
#
#      Per-term bond and angle edits (2026-10-07, Hessian / modified Seminario,
#      pdynamo/seminario.py):
#
#        edits["bond_terms"]  = { "a-b": ( b0, k ) }               (a < b)
#        edits["angle_terms"] = { "a-b-c": [ ( c, n ), ... ] }     (b central, a < c)
#
#      replace the DYFF value of that one term (same unique-key rule).
#
#      Per-centre out-of-plane edits (2026-10-07, Hessian compliance fit,
#      pdynamo/seminario.py outofplane_edits):
#
#        edits["outofplane_centres"] = { "j": [ ( c, n ), ... ] }  (j central atom)
#
#      DYFF builds THREE out-of-plane terms per trigonal / group-15 centre (each
#      neighbour in turn as the out-of-plane atom, pDynamo's CosineOutOfPlane
#      container stores [ i, j, k, l ] with j central), all from ONE factory call
#      whose indices start with the central atom: the edit replaces the series
#      shared by the three.
#
#      The edits live on the Builder object (vismol_object.dyff_parameter_
#      edits, kept by Undo) and are used by the DYFF Parameters window,
#      define_MMModel('DYFF') and the Builder's Optimize.
#

import os
import threading

from pMolecule.MMModel import MMModelDYFF

_TYPER_LOCK = threading.RLock ( )

_FACTORY_TERM = { "harmonicBondParameters"     : "bond",
                  "cosineAngleParameters"      : "angle",
                  "cosineDihedralParameters"   : "dihedral",
                  "cosineOutOfPlaneParameters" : "outofplane" }
KEY_SEPARATOR = ":"


#=============================================================================
#  Angles at the bent resonant types of dyff-1.0-eh (2026-10-08)
#=============================================================================
#  DYFF (as UFF) gives every angle centred on a "Res"/"Tri" atom the trigonal three-fold
#  cosine form, whose minimum is 120 deg whatever the type's natural angle. Right for C, N,
#  B in six-rings; wrong for a two-connected heteroatom whose natural angle is ~87-101 deg
#  (P in phosphinines, Se in selenophenes): the validation against X-ray structures of the
#  dyff-1.0-eh types (DYFF_validation/report/new_types_report.html) found C-P-C +7.8 deg and
#  C-Se-C +3.1 deg too wide. For these types only, the angle uses the UFF GENERAL form, which
#  has its minimum at the type's natural angle; the torsions keep the resonant (Res) terms.
#  dyff-1.0 types (S:Res included) are not touched.
BENT_RESONANT_TYPES = ( "P:Res", "Se:Res" )

def _bent_resonant_angles ( atomTypes, base ):
    import math
    from pMolecule.MMModel.DYFFUtilities import DYFFBondOrder, DYFFNaturalBondLength
    from pScientific import Units
    def f ( types, indices, connectivity, key ):
        ( iT, jT, kT ) = types
        if jT not in BENT_RESONANT_TYPES:
            return base ( types, indices, connectivity, key )
        ( i, j, k ) = indices
        ( ijOrder, ijTag ) = DYFFBondOrder ( i, j, connectivity )
        ( jkOrder, jkTag ) = DYFFBondOrder ( j, k, connectivity )
        theta = math.radians ( atomTypes.GetItem ( jT ).GetProperty ( "Bond Angle" ) )
        zI = atomTypes.GetItem ( iT ).GetProperty ( "Effective Charge" )
        zK = atomTypes.GetItem ( kT ).GetProperty ( "Effective Charge" )
        rIJ = DYFFNaturalBondLength ( atomTypes, iT, jT, ijOrder )
        rJK = DYFFNaturalBondLength ( atomTypes, jT, kT, jkOrder )
        sinT, cosT = math.sin ( theta ), math.cos ( theta )
        rIK = math.sqrt ( rIJ ** 2 + rJK ** 2 - 2.0 * rIJ * rJK * cosT )
        fC0 = 2.0 * Units.Energy_E2Angstroms_To_Kilojoules_Per_Mole * zI * zK * ( 3.0 * rIJ * rJK * ( 1.0 - cosT ** 2 ) - rIK ** 2 * cosT ) / rIK ** 5
        results = [ ]
        if fC0 != 0.0:                                 # . UFF general case (natural angle theta)
            fC0 /= ( 4.0 * sinT ** 2 )
            results = [ ( fC0 * ( 1.0 + 2.0 * cosT ** 2 ), 0 ), ( - fC0 * 4.0 * cosT, 1 ), ( fC0, 2 ) ]
        return ( results, tuple ( list ( key ) + [ ijTag, jkTag, "bent" ] ) )
    return f


#  Metal-centred angles (2026-10-13, tmQM study). UFF has angle forms for Lin, Tri/Res (cos 3),
#  SPl/Oct (cos 4: minima at 90 and 180) and the general form (one natural angle). A square
#  pyramid (90/180) takes the Oct form; a trigonal bipyramid (90/120/180) or a pentagonal bipyramid
#  has no single natural angle, so each L-M-L angle of those centres gets its own ideal angle
#  (metal_geometry.perceive, set in BuildModel), in the general form -- or the linear form at 180.
OCTAHEDRAL_FORM_GEOMETRIES = ( "SPy", )

def _uff_angle_terms ( atomTypes, types, indices, connectivity, theta_deg ):
    import math
    from pMolecule.MMModel.DYFFUtilities import DYFFBondOrder, DYFFNaturalBondLength
    from pScientific import Units
    ( iT, jT, kT ) = types
    ( i, j, k ) = indices
    ( ijOrder, ijTag ) = DYFFBondOrder ( i, j, connectivity )
    ( jkOrder, jkTag ) = DYFFBondOrder ( j, k, connectivity )
    theta = math.radians ( theta_deg )
    zI = atomTypes.GetItem ( iT ).GetProperty ( "Effective Charge" )
    zK = atomTypes.GetItem ( kT ).GetProperty ( "Effective Charge" )
    rIJ = DYFFNaturalBondLength ( atomTypes, iT, jT, ijOrder )
    rJK = DYFFNaturalBondLength ( atomTypes, jT, kT, jkOrder )
    sinT, cosT = math.sin ( theta ), math.cos ( theta )
    rIK = math.sqrt ( rIJ ** 2 + rJK ** 2 - 2.0 * rIJ * rJK * cosT )
    fC0 = 2.0 * Units.Energy_E2Angstroms_To_Kilojoules_Per_Mole * zI * zK * ( 3.0 * rIJ * rJK * ( 1.0 - cosT ** 2 ) - rIK ** 2 * cosT ) / rIK ** 5
    if fC0 == 0.0: return [ ], ( ijTag, jkTag )
    if abs ( theta - math.pi ) < 1.0e-3:
        return [ ( fC0, 0 ), ( fC0, 1 ) ], ( ijTag, jkTag )
    fC0 /= ( 4.0 * sinT ** 2 )
    return [ ( fC0 * ( 1.0 + 2.0 * cosT ** 2 ), 0 ), ( - fC0 * 4.0 * cosT, 1 ), ( fC0, 2 ) ], ( ijTag, jkTag )

# . M-C(eta) bond offsets per metal (Angstrom, atomic number -> offset added to the natural length): the
#   radii are calibrated on classic (sigma) complexes; eta bonds are longer than a sigma M-C of the same
#   type. Values from DYFF_validation/scripts/tmqm_calibrate_hapto.py (tmQM TRAIN eta complexes).
HAPTO_BOND_OFFSETS = { 21: +0.142, 22: +0.281, 23: +0.251, 24: +0.392, 25: +0.359, 26: +0.247, 27: +0.132, 28: +0.053, 29: -0.064, 30: +0.076, 39: +0.030, 40: +0.268, 41: +0.310, 42: +0.233, 44: +0.302, 45: +0.192, 46: +0.082, 47: +0.027, 72: +0.303, 73: +0.294, 74: +0.174, 75: +0.321, 76: +0.313, 77: +0.205, 78: +0.084, 79: +0.173 }

# . High-spin metal ions (2026-10-13): the metal radii are calibrated on tmQM, almost entirely closed-shell
#   (low-spin) complexes; high-spin Mn(II)/Fe(II)/Co(II) bind with longer M-L bonds. Bonds of the metal atoms
#   flagged high-spin (MMModelDYFFEdited.highSpinAtoms / dyff_model ( high_spin_atoms = ... ), e.g. from the
#   site multiplicity of a protein preparation) get this offset (A, atomic number -> offset). The spin state
#   cannot be read from the structure, so nothing is applied without the flag. Values from
#   DYFF_validation/scripts/metalpdb_calibrate_hs.py (MetalPDB calibration set, 1.5-2.0 A, not the validation set).
HIGH_SPIN_BOND_OFFSETS = { 25: +0.384, 26: +0.326, 27: +0.303 }   # . Mn, Fe, Co

def _metal_bonds ( base, hapto_bonds, high_spin_bonds = None ):
    high_spin_bonds = high_spin_bonds or { }
    def f ( types, indices, connectivity, key ):
        value, new_key = base ( types, indices, connectivity, key )
        pair = ( int ( indices[0] ), int ( indices[1] ) )
        tags, offset = [ ], 0.0
        z = hapto_bonds.get ( pair )
        if z is not None:
            offset += HAPTO_BOND_OFFSETS.get ( z, 0.0 ); tags.append ( "hapto" )
        z = high_spin_bonds.get ( pair )
        if z is not None:
            offset += HIGH_SPIN_BOND_OFFSETS.get ( z, 0.0 ); tags.append ( "high-spin" )
        if not tags: return value, new_key
        if offset == 0.0: return value, tuple ( list ( new_key ) + tags )
        ( eq, fC ) = value
        return ( eq + offset, fC ), tuple ( list ( new_key ) + tags + [ "{:+.3f}".format ( offset ) ] )
    return f

def _metal_angles ( atomTypes, base, per_instance ):
    from pMolecule.MMModel.DYFFUtilities import DYFFHybridization
    def f ( types, indices, connectivity, key ):
        idx = tuple ( int ( x ) for x in indices )
        theta = per_instance.get ( idx )
        if idx in per_instance and theta is None:                  # . angle with an eta-bonded atom: no term
            return ( [ ], tuple ( list ( key ) + [ "hapto" ] ) )
        if theta is not None:
            terms, tags = _uff_angle_terms ( atomTypes, types, indices, connectivity, theta )
            return ( terms, tuple ( list ( key ) + list ( tags ) + [ "ideal {:.0f}".format ( theta ) ] ) )
        if DYFFHybridization ( types[1] ) in OCTAHEDRAL_FORM_GEOMETRIES:
            # . cos 4 form of SPl/Oct (pDynamo DYFFCosineAngleParameters): fC0/16 * ( 1 - cos 4 theta )
            import math
            from pMolecule.MMModel.DYFFUtilities import DYFFBondOrder, DYFFNaturalBondLength
            from pScientific import Units
            ( iT, jT, kT ) = types
            ( i, j, k ) = indices
            ( ijOrder, ijTag ) = DYFFBondOrder ( i, j, connectivity )
            ( jkOrder, jkTag ) = DYFFBondOrder ( j, k, connectivity )
            zI = atomTypes.GetItem ( iT ).GetProperty ( "Effective Charge" )
            zK = atomTypes.GetItem ( kT ).GetProperty ( "Effective Charge" )
            rIJ = DYFFNaturalBondLength ( atomTypes, iT, jT, ijOrder )
            rJK = DYFFNaturalBondLength ( atomTypes, jT, kT, jkOrder )
            rIK = math.sqrt ( rIJ ** 2 + rJK ** 2 )                  # . natural angle 90 deg: cos = 0
            fC0 = 2.0 * Units.Energy_E2Angstroms_To_Kilojoules_Per_Mole * zI * zK * 3.0 * rIJ * rJK / rIK ** 5
            return ( [ ( fC0 / 16.0, 0 ), ( - fC0 / 16.0, 4 ) ], tuple ( list ( key ) + [ ijTag, jkTag ] ) )
        return base ( types, indices, connectivity, key )
    return f


class MMModelDYFFEdited ( MMModelDYFF ):
    """ DYFF with parameters overridden by key (see the module docstring).

        typeOverrides (2026-10-08): { atom index: DYFF type label } for atoms the DYFF
        pattern typer leaves UNTYPED -- e.g. a free Ca2+ or Cu2+ (DYFF has no free-ion
        pattern for most metals) or a metal with an unusual coordination
        (pdynamo/dyff_prep.type_overrides). pDynamo checks for untyped atoms inside
        MMAtomTyper.TypeAtoms, before any later substitution, so the typer is wrapped
        for the duration of BuildModel. """

    _attributable = dict ( MMModelDYFF._attributable )
    _attributable.update ( { "parameterEdits" : None, "typeOverrides" : None, "highSpinAtoms" : None } )

    def BuildModel ( self, target, *args, **kwargs ):
        self._preset_metal_geometry ( target )
        overrides = self.typeOverrides or { }
        if not overrides:
            return super ( MMModelDYFFEdited, self ).BuildModel ( target, *args, **kwargs )
        import threading
        from pMolecule.MMModel.MMAtomTyper import MMAtomTyper
        with _TYPER_LOCK:
            original = MMAtomTyper.TypeAtoms
            def type_atoms ( typer, connectivity, sequence, log ):
                atom_charges = [ 0.0  for i in range ( len ( connectivity.atoms ) ) ]
                atom_types   = [ None for i in range ( len ( connectivity.atoms ) ) ]
                untyped = set ( connectivity.atoms )
                typer.TypeBySequence ( sequence, atom_types, atom_charges, untyped )
                typer.TypeByPattern  ( connectivity, atom_types, atom_charges, untyped )
                for atom in list ( untyped ):
                    label = overrides.get ( atom.index )
                    if label is None: continue
                    try:
                        known = typer.mmAtomTypes.GetItem ( label ) is not None
                    except Exception:
                        known = False
                    if known:
                        atom_types[atom.index] = label
                        untyped.discard ( atom )
                typer.CheckUntypedAtoms ( connectivity, untyped, log )
                return ( atom_types, atom_charges )
            MMAtomTyper.TypeAtoms = type_atoms
            try:
                return super ( MMModelDYFFEdited, self ).BuildModel ( target, *args, **kwargs )
            finally:
                MMAtomTyper.TypeAtoms = original

    def _preset_metal_geometry ( self, target ):
        """ Transition metals (2026-10-13): geometry from the coordinates (pdynamo.metal_geometry), preset on
            the connectivity atoms before pDynamo's DetermineAtomGeometry (which only counts neighbours and
            keeps a preset value), plus the per-instance ideal angles of TBp / PBp centres. EH set only;
            EASYHYBRID_DYFF_NO_METAL_GEOMETRY=1 restores pDynamo's neighbour-count geometry. """
        self._metal_angles = { }; self._hapto_bonds = { }; self._high_spin_bonds = { }
        if getattr ( self, "parameterSet", None ) == EH_DYFF_SET and getattr ( self, "highSpinAtoms", None ):
            try:
                for i, j in target.connectivity.bondIndices:
                    for m, l in ( ( i, j ), ( j, i ) ):
                        if m in self.highSpinAtoms:
                            self._high_spin_bonds[( m, l )] = self._high_spin_bonds[( l, m )] = target.atoms[m].atomicNumber
            except Exception:
                self._high_spin_bonds = { }
        if getattr ( self, "parameterSet", None ) != EH_DYFF_SET or os.environ.get ( "EASYHYBRID_DYFF_NO_METAL_GEOMETRY" ): return
        try:
            connectivity = target.connectivity; coordinates = target.coordinates3
            if connectivity is None or coordinates is None: return
            from pdynamo.metal_geometry import TRANSITION_METALS, perceive_system
            if not any ( a.atomicNumber in TRANSITION_METALS for a in target.atoms ): return
            undefined = set ( getattr ( coordinates, "undefined", None ) or [ ] )
            self._hapto_bonds = { }
            mode = os.environ.get ( "EASYHYBRID_DYFF_HAPTO_ANGLES", "keep" )   # . keep: better shell RMSD on tmQM (LAB_NOTEBOOK)
            for m, info in perceive_system ( target ).items ( ):
                if m in undefined or info["tag"] is None: continue
                connectivity.atoms[m].geometry = info["tag"]
                for ( i, j, k ), theta in ( info["angles"] or { } ).items ( ):
                    if theta is None and mode != "drop": continue        # . "keep": the type's own angle form
                    self._metal_angles[( i, j, k )] = self._metal_angles[( k, j, i )] = theta
                symbol = target.atoms[m].symbol if hasattr ( target.atoms[m], "symbol" ) else None
                for c in info.get ( "hapto_atoms", ( ) ):
                    self._hapto_bonds[( m, c )] = self._hapto_bonds[( c, m )] = target.atoms[m].atomicNumber
        except Exception:
            self._metal_angles = { }

    def ParameterFactories ( self, atomTypes ):
        self._eh_atom_types = atomTypes                      # . kept for pdynamo.metal_refine (UFF force constants)
        factories = super ( MMModelDYFFEdited, self ).ParameterFactories ( atomTypes )
        if getattr ( self, "parameterSet", None ) == EH_DYFF_SET and "cosineAngleParameters" in factories:
            factories["cosineAngleParameters"] = _bent_resonant_angles ( atomTypes, factories["cosineAngleParameters"] )
            factories["cosineAngleParameters"] = _metal_angles ( atomTypes, factories["cosineAngleParameters"], getattr ( self, "_metal_angles", { } ) )
        if getattr ( self, "parameterSet", None ) == EH_DYFF_SET and "harmonicBondParameters" in factories and \
           ( getattr ( self, "_hapto_bonds", None ) or getattr ( self, "_high_spin_bonds", None ) ):
            factories["harmonicBondParameters"] = _metal_bonds ( factories["harmonicBondParameters"], getattr ( self, "_hapto_bonds", { } ) or { },
                                                                 getattr ( self, "_high_spin_bonds", { } ) )
        edits = self.parameterEdits or { }
        bonds = edits.get ( "dihedral_bonds" ) or { }
        per_term = { "bond": edits.get ( "bond_terms" ) or { }, "angle": edits.get ( "angle_terms" ) or { },
                     "outofplane": edits.get ( "outofplane_centres" ) or { } }
        wrapped = { }
        for name, factory in factories.items ( ):
            term  = _FACTORY_TERM.get ( name )
            table = edits.get ( term ) or { }
            if term in per_term and per_term[term]:
                def per_atoms ( types, indices, connectivity, key, factory = factory, table = table, term = term, terms = per_term[term] ):
                    idx = [ int ( i ) for i in indices ]
                    if term == "angle" and idx[0] > idx[2]: idx = idx[::-1]
                    if term == "bond" and idx[0] > idx[1]: idx = idx[::-1]
                    label = str ( idx[0] ) if term == "outofplane" else "-".join ( map ( str, idx ) )
                    override = terms.get ( label )
                    value, new_key = factory ( types, indices, connectivity, key )
                    if override is not None:
                        if term == "bond":
                            return ( float ( override[0] ), float ( override[1] ) ), tuple ( new_key ) + ( "term " + label, )
                        return [ ( float ( fc ), int ( n ) ) for ( fc, n ) in override ], tuple ( new_key ) + ( "term " + label, )
                    edit = table.get ( KEY_SEPARATOR.join ( new_key ) )
                    if edit is None: return value, new_key
                    if term == "bond": return ( float ( edit[0] ), float ( edit[1] ) ), new_key
                    return [ ( float ( fc ), int ( n ) ) for ( fc, n ) in edit ], new_key
                wrapped[name] = per_atoms
                continue
            if term == "dihedral" and bonds:
                def per_bond ( types, indices, connectivity, key, factory = factory, table = table, bonds = bonds ):
                    j, k = int ( indices[1] ), int ( indices[2] )
                    label = "{}-{}".format ( min ( j, k ), max ( j, k ) )
                    override = bonds.get ( label )
                    if override is not None:
                        _, new_key = factory ( types, indices, connectivity, key )
                        return [ ( float ( fc ), int ( n ) ) for ( fc, n ) in override ], tuple ( new_key ) + ( "bond " + label, )
                    value, new_key = factory ( types, indices, connectivity, key )
                    edit = table.get ( KEY_SEPARATOR.join ( new_key ) )
                    if edit is None: return value, new_key
                    return [ ( float ( fc ), int ( n ) ) for ( fc, n ) in edit ], new_key
                wrapped[name] = per_bond
                continue
            if not table:
                wrapped[name] = factory
                continue
            def edited ( types, indices, connectivity, key, factory = factory, table = table, term = term ):
                value, new_key = factory ( types, indices, connectivity, key )
                override = table.get ( KEY_SEPARATOR.join ( new_key ) )
                if override is None:
                    return value, new_key
                if term == "bond":
                    return ( float ( override[0] ), float ( override[1] ) ), new_key
                return [ ( float ( fc ), int ( n ) ) for ( fc, n ) in override ], new_key
            wrapped[name] = edited
        return wrapped


#=============================================================================
#  Parameter set: EasyHybrid's extension of dyff-1.0 (2026-10-08)
#=============================================================================
#  src/pdynamo/data/forcefields/dyff/dyff-1.0-eh = pDynamo3 dyff-1.0 + ADDED types/patterns
#  (aromatic P and Se, aromatic B, borane clusters; see the headers of its atomTypes.yaml and
#  patterns.yaml). Every atom typed by dyff-1.0 keeps its type. "dyff-1.0" requests are served
#  by the extension; EASYHYBRID_DYFF_ORIGINAL=1 restores pDynamo's own set.
EH_DYFF_ROOT = os.path.join ( os.path.dirname ( os.path.abspath ( __file__ ) ), "data", "forcefields", "dyff" )
EH_DYFF_SET  = "dyff-1.0-eh"
# . 2026-10-21: explicit request for pDynamo's own dyff-1.0 (Prepare DYFF System's version choice)
PDYNAMO_DYFF_SET = "dyff-1.0-pdynamo"

def dyff_parameter_set ( parameter_set = "dyff-1.0" ):
    """ ( name, directory containing it ) of the DYFF parameter set to use. """
    if parameter_set == PDYNAMO_DYFF_SET:
        return "dyff-1.0", os.path.join ( os.environ.get ( "PDYNAMO3_PARAMETERS", "" ), "forceFields", "dyff" )
    if parameter_set in ( "dyff-1.0", EH_DYFF_SET ) and not os.environ.get ( "EASYHYBRID_DYFF_ORIGINAL" ) \
       and os.path.isdir ( os.path.join ( EH_DYFF_ROOT, EH_DYFF_SET ) ):
        return EH_DYFF_SET, EH_DYFF_ROOT
    return parameter_set, os.path.join ( os.environ.get ( "PDYNAMO3_PARAMETERS", "" ), "forceFields", "dyff" )

def dyff_parameter_dir ( parameter_set = "dyff-1.0" ):
    """ Full path of the parameter set directory (atomTypes.yaml, patterns.yaml, ...). """
    name, root = dyff_parameter_set ( parameter_set )
    return os.path.join ( root, name )


def dyff_model ( parameter_set = "dyff-1.0", edits = None, type_overrides = None, high_spin_atoms = None ):
    """ MMModelDYFF, or MMModelDYFFEdited when there are edits or type overrides (call
        applicable_edits() first when the edits may hold per-term values). """
    edits = { term: dict ( table ) for term, table in ( edits or { } ).items ( )
              if table and isinstance ( table, dict ) }
    name, root = dyff_parameter_set ( parameter_set )
    if not edits and not type_overrides and name != EH_DYFF_SET:
        return MMModelDYFF.WithParameterSet ( name, path = root )
    model = MMModelDYFFEdited.WithParameterSet ( name, path = root )
    model.parameterEdits = edits
    model.typeOverrides = { int ( k ): v for k, v in ( type_overrides or { } ).items ( ) }
    model.highSpinAtoms = set ( int ( i ) for i in ( high_spin_atoms or ( ) ) ) or None
    return model


#=============================================================================
#  Per-term edits are tied to atom indices: signature of the molecule
#=============================================================================
PER_TERM_KEYS = ( "bond_terms", "angle_terms", "dihedral_bonds", "outofplane_centres" )

def term_signature ( n_atoms, bonds ):
    """ Atom count + bond list: per-term edits are valid only for this topology. """
    import hashlib
    pairs = sorted ( ( min ( int ( a ), int ( b ) ), max ( int ( a ), int ( b ) ) ) for a, b in bonds )
    return "{}:{}".format ( n_atoms, hashlib.sha1 ( repr ( pairs ).encode ( ) ).hexdigest ( )[:16] )

def applicable_edits ( edits, n_atoms = None, bonds = None ):
    """ The edits that apply to a molecule with n_atoms atoms and these bonds:
        per-term edits (indexed by atoms) are dropped when the molecule changed
        since they were made (edits["term_signature"] differs). Returns
        ( edits, dropped ). """
    if not edits: return edits, False
    edits = dict ( edits )
    signature = edits.pop ( "term_signature", None )
    if signature is None or n_atoms is None or bonds is None:
        return edits, False
    if signature == term_signature ( n_atoms, bonds ):
        return edits, False
    dropped = any ( edits.get ( k ) for k in PER_TERM_KEYS )
    for k in PER_TERM_KEYS: edits.pop ( k, None )
    return edits, dropped


#=============================================================================
#  Polar hydrogens with partial charges (2026-10-07)
#=============================================================================
#  DYFF gives EVERY hydrogen the same Lennard-Jones type (sigma 2.571 A,
#  epsilon 0.184 kJ/mol). With partial charges (CM5 / residue library) a
#  hydrogen bond puts a polar H ~1.9 A from its acceptor, deep inside that
#  wall: on crambin the DYFF LJ forces were 187 kJ/mol/A (heavy-atom RMS)
#  against 37 for OPLS, which gives polar H no LJ at all (the charges keep
#  them off the acceptor). A small polar-H type (as CHARMM) fixes it: force
#  RMSE vs GFN2-xTB 306 -> 247 kJ/mol/A, the same as swapping in OPLS's whole
#  nonbonded model (246); sigma anywhere in 0-1.6 A gives the same result.
#  NOT applied without charges: nothing would then keep a polar H off an
#  acceptor.
#  2026-10-08 (DYFF_validation Phase IV, NCI Atlas): the type is also given to
#  H bonded to F/Cl/Br/I (HF, HCl ... as H-bond donors): with the UFF H their
#  H-bond curves were off by 25-380 kcal/mol; NCI validation equilibrium MAE
#  3.10 -> 2.38 kcal/mol, no new parameter. The fit of the polar-H sigma /
#  epsilon and of a global LJ scale did not transfer to protein-ligand
#  complexes and were not adopted (report: DYFF_validation/report/phase4_report.html).
POLAR_H_SIGMA   = 1.0       # A
POLAR_H_EPSILON = 0.192     # kJ/mol (CHARMM polar H, 0.046 kcal/mol)
POLAR_H_LABEL   = "H.polar"
POLAR_H_PARTNERS = ( 7, 8, 16, 9, 17, 35, 53 )   # N, O, S, F, Cl, Br, I

def polar_hydrogens ( system ):
    """ Indices of the hydrogens bonded to N, O, S or a halogen (F, Cl, Br, I). """
    numbers = [ a.atomicNumber for a in system.atoms ]
    polar = set ( )
    for i, j in system.connectivity.bondIndices:
        for h, x in ( ( i, j ), ( j, i ) ):
            if numbers[h] == 1 and numbers[x] in POLAR_H_PARTNERS: polar.add ( h )
    return sorted ( polar )

def soften_polar_hydrogens ( system, sigma = POLAR_H_SIGMA, epsilon = POLAR_H_EPSILON ):
    """ Gives the polar hydrogens of a DYFF system their own small LJ type (see
        above), in the full and the 1-4 tables (same 1-4 scale as the model).
        Call after DefineMMModel, only when the system has partial charges.
        Returns the number of hydrogens changed. """
    from pMolecule.MMModel import LJParameterContainer
    state = system.mmState
    polar = polar_hydrogens ( system )
    if not polar or state.ljParameters is None: return 0
    lj = state.ljParameters.__getstate__ ( )
    keys = list ( lj["parameterKeys"] )
    if POLAR_H_LABEL in keys: return 0
    scale14 = float ( getattr ( system.mmModel, "lennardJonesScale14", 0.5 ) )
    keys.append ( POLAR_H_LABEL )
    epsilons = list ( lj["epsilons"] ) + [ float ( epsilon ) ]
    sigmas   = list ( lj["sigmas"] )   + [ float ( sigma ) ]
    state.ljParameters   = LJParameterContainer.FromEpsilonSigma ( len ( keys ), epsilons, sigmas,
                                                                    analyticForm = lj["analyticForm"], parameterKeys = keys )
    state.ljParameters14 = LJParameterContainer.FromEpsilonSigma ( len ( keys ), [ scale14 * e for e in epsilons ], sigmas,
                                                                    analyticForm = lj["analyticForm"], parameterKeys = keys )
    for h in polar: state.ljTypeIndices[h] = len ( keys ) - 1
    return len ( polar )


#  TIP3P water oxygen LJ (2026-10-08, residue_library.solvent_and_ions): water with TIP3P
#  charges needs the TIP3P oxygen LJ (sigma 3.15061 A, epsilon 0.6364 kJ/mol) -- DYFF's own
#  O (3.118, 0.251) would collapse the liquid. Hydrogens: the polar-H type.
WATER_O_SIGMA, WATER_O_EPSILON, WATER_O_LABEL = 3.15061, 0.6364, "OW.tip3p"

def water_oxygens ( system ):
    numbers = [ a.atomicNumber for a in system.atoms ]
    neighbours = [ [ ] for _ in numbers ]
    for i, j in system.connectivity.bondIndices:
        neighbours[i].append ( j ); neighbours[j].append ( i )
    return [ a for a, z in enumerate ( numbers ) if z == 8 and len ( neighbours[a] ) == 2 and
             all ( numbers[h] == 1 and len ( neighbours[h] ) == 1 for h in neighbours[a] ) ]

def water_lj ( system ):
    """ TIP3P oxygen LJ type for the water oxygens (full and 1-4 tables). Call after
        DefineMMModel, with partial charges. Returns the number of oxygens changed. """
    from pMolecule.MMModel import LJParameterContainer
    state = system.mmState
    oxygens = water_oxygens ( system )
    if not oxygens or state.ljParameters is None: return 0
    lj = state.ljParameters.__getstate__ ( )
    keys = list ( lj["parameterKeys"] )
    if WATER_O_LABEL in keys: return 0
    scale14 = float ( getattr ( system.mmModel, "lennardJonesScale14", 0.5 ) )
    keys.append ( WATER_O_LABEL )
    epsilons = list ( lj["epsilons"] ) + [ WATER_O_EPSILON ]
    sigmas   = list ( lj["sigmas"] )   + [ WATER_O_SIGMA ]
    state.ljParameters   = LJParameterContainer.FromEpsilonSigma ( len ( keys ), epsilons, sigmas, analyticForm = lj["analyticForm"], parameterKeys = keys )
    state.ljParameters14 = LJParameterContainer.FromEpsilonSigma ( len ( keys ), [ scale14 * e for e in epsilons ], sigmas,
                                                                    analyticForm = lj["analyticForm"], parameterKeys = keys )
    for o in oxygens: state.ljTypeIndices[o] = len ( keys ) - 1
    return len ( oxygens )
