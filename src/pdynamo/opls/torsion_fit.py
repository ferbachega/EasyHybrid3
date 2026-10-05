#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: OPLS side of the xTB torsion fit (pdynamo/torsion_fit.py)
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      An OPLS system is usually a ligand inside a protein and/or solvent,
#      so the scan is done on the MOLECULE that holds the dihedral (its
#      connected component), which must be small enough for xTB. That
#      molecule is cut out of the system with System.Prune: types, charges
#      (CM5 of the ligand, overrides) and every MM term stay exactly those
#      of the system; the NB model is replaced by all pairs, no cutoff (gas
#      phase, like the xTB reference) and restraints are dropped.
#
#      MM energy with the parameter of the key replaced:
#
#          E ( s ) = E_fragment ( s ) - E_key,now ( s ) + sum_d sum_n c_n cos ( n phi_d ( s ) )
#
#      where E_key,now is the pDynamo Fourier form of the current parameter,
#      sum k ( 1 + cos ( n phi - delta ) ) (kcal/mol, converted to kJ/mol)
#      over the terms of the key in the molecule (test_opls_torsion_fit
#      checks it against a real rebuild of the MM model from edited files).
#      The fitted series is written back in the OPLS form, k = |c_n| (kcal/mol),
#      delta = 0 if c_n > 0 else 180 (the constant k is irrelevant), with
#      editor.apply_edits() -- the same path as a manual edit.
#

import math
from collections import defaultdict

import numpy as np

from pdynamo.torsion_fit import dihedral_angle, mm_energies

KCAL_TO_KJ    = 4.184
MAX_MOLECULE  = 250          # atoms; above this an xTB scan of the whole molecule is not offered


def connected_component ( system, start ):
    adjacency = defaultdict ( set )
    for ( i, j ) in system.connectivity.bondIndices:
        adjacency[i].add ( j ); adjacency[j].add ( i )
    seen, stack = { start }, [ start ]
    while stack:
        u = stack.pop ( )
        for v in adjacency[u]:
            if v not in seen:
                seen.add ( v ); stack.append ( v )
    return sorted ( seen )


def fourier_energy ( value, phi_degrees ):
    """ pDynamo OPLS dihedral energy (kJ/mol) of value [ ( k kcal/mol, n, phase deg ) ]. """
    phi = math.radians ( phi_degrees )
    return KCAL_TO_KJ * sum ( float ( k ) * ( 1.0 + math.cos ( int ( n ) * phi - math.radians ( float ( p ) ) ) )
                              for ( k, n, p ) in ( value or [ ] ) )


def pairs_to_fourier ( pairs ):
    """ [ ( c kJ/mol, n ) ] of sum c cos ( n phi ) -> OPLS [ ( k kcal/mol, n, phase ) ]. """
    terms = [ ( round ( abs ( c ) / KCAL_TO_KJ, 6 ), int ( n ), 0.0 if c > 0.0 else 180.0 ) for ( c, n ) in pairs if int ( n ) > 0 and c != 0.0 ]
    return terms or [ ( 0.0, 1, 0.0 ) ]


def format_fourier ( value ):
    return "; ".join ( "{:.4f} {:d} {:.1f}".format ( k, int ( n ), p ) for ( k, n, p ) in value if float ( k ) != 0.0 ) or "none"


class OPLSTorsionTarget:
    """ What the torsion fit window needs from an OPLS system: the molecule
        holding the dihedral, the parameter of a Dihedrals row of
        editor.build_tables() and how to apply a fit. coordinates: N x 3 of
        the whole system (default: system.coordinates3). """

    force_field = "OPLS"
    units_note  = "k n phase (kcal/mol, deg)"

    def __init__ ( self, system, row, quad, coordinates = None ):
        from pCore import Selection
        from pMolecule.NBModel import NBModelFull
        from pScientific import PeriodicTable
        if getattr ( system, "qcModel", None ) is not None:
            raise ValueError ( "This system has a QC model: pDynamo removes it when the MM model is rebuilt with the fitted "
                               "parameter. Remove the QC region first (define it again afterwards)." )
        self.system, self.row = system, row
        self.working_folder = getattr ( system, "e_working_folder", None )
        self.key = tuple ( row["key"] )
        self.key_text = " - ".join ( self.key )
        quad = tuple ( int ( i ) for i in quad )
        molecule = connected_component ( system, quad[0] )
        if len ( molecule ) > MAX_MOLECULE:
            raise ValueError ( "This dihedral belongs to a molecule of {} atoms (a polymer?): the xTB scan is offered for "
                               "molecules up to {} atoms (ligands).".format ( len ( molecule ), MAX_MOLECULE ) )
        local = { a: k for k, a in enumerate ( molecule ) }
        self.molecule = molecule
        self.n_system_atoms = len ( system.atoms )
        if coordinates is None:
            coordinates = np.array ( [ [ system.coordinates3[a, c] for c in range ( 3 ) ] for a in range ( len ( system.atoms ) ) ] )
        coordinates = np.asarray ( coordinates, dtype = float )
        atoms = list ( system.atoms )
        self.symbols = [ PeriodicTable.Symbol ( atoms[a].atomicNumber ) for a in molecule ]
        self.coordinates = np.array ( [ coordinates[a] for a in molecule ] )
        self.bonds = [ ( local[i], local[j] ) for ( i, j ) in system.connectivity.bondIndices if i in local and j in local ]
        self.quad = tuple ( local[a] for a in quad )
        self.atom_names = [ atoms[a].label for a in quad ]
        all_instances = [ tuple ( q ) for q in row["instances"] ]
        self.instances = [ tuple ( local[a] for a in q ) for q in all_instances if all ( a in local for a in q ) ]
        self.n_outside = len ( all_instances ) - len ( self.instances )
        self.current = [ tuple ( t ) for t in ( row["value"] or [ ] ) ]
        self.current_text = format_fourier ( self.current )
        self.default_charge = int ( round ( sum ( float ( system.mmState.charges[a] ) for a in molecule ) ) )
        # . the molecule alone: same MM terms and charges, gas phase
        fragment = system.Prune ( Selection.FromIterable ( molecule ) )
        if getattr ( fragment, "restraintModel", None ) is not None:
            fragment.DefineRestraintModel ( None )
        fragment.DefineNBModel ( NBModelFull.WithDefaults ( ) )
        self.fragment = fragment

    def notes ( self ):
        notes = [ ]
        if self.n_outside:
            notes.append ( "{} term(s) with these types are in other molecules (other copies of the ligand, protein): "
                           "the fitted parameter changes them too.".format ( self.n_outside ) )
        notes.append ( "Scanned molecule: {} atoms, cut out of the system with its OPLS types and charges; "
                       "MM side in vacuum (all pairs).".format ( len ( self.molecule ) ) )
        return notes

    def _key_energy ( self, value_or_pairs, geometries, pairs ):
        energies = [ ]
        for xyz in geometries:
            phis = [ dihedral_angle ( xyz, q ) for q in self.instances ]
            if pairs:
                energies.append ( sum ( c * math.cos ( n * math.radians ( p ) ) for p in phis for ( c, n ) in value_or_pairs ) )
            else:
                energies.append ( sum ( fourier_energy ( value_or_pairs, p ) for p in phis ) )
        return np.array ( energies )

    def energies_for ( self, pairs, geometries ):
        base = np.asarray ( mm_energies ( self.fragment, geometries ) )
        if pairs is None: return base
        return base - self._key_energy ( self.current, geometries, False ) + self._key_energy ( pairs, geometries, True )

    def format_value ( self, pairs ):
        return format_fourier ( pairs_to_fourier ( pairs ) )

    def changed ( self ):
        if len ( self.system.atoms ) != self.n_system_atoms:
            return "The system changed during the scan: run it again."
        return None

    def apply ( self, pairs ):
        from pdynamo.opls import editor
        summary = editor.apply_edits ( self.system, { "dihedral": { self.key: pairs_to_fourier ( pairs ) } } )
        return "Fitted parameter written to the OPLS parameter set ({} term(s) in the system); MM model rebuilt.".format (
               len ( self.instances ) + self.n_outside ) if summary.get ( "rebuilt" ) else "Applied."
