#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: DYFF parameters of a Builder object (backend of the DYFF Parameters window)
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      DYFF builds every bonded parameter from RULES (UFF-like: atom type,
#      hybridization, bond order) -- its parameter files are empty. So there
#      is no "database vs estimated" here; what tells a good assignment from
#      a bad one is whether the force field AGREES WITH THE STRUCTURE. For
#      every term the strain energy at the current geometry is computed
#      (energy above the minimum of the term's own analytic form):
#
#          green   < 2 kJ/mol        yellow   < 8 kJ/mol        red   >= 8 kJ/mol
#
#      A large strain on a freshly built/optimized molecule usually means a
#      wrong atom type (hybridization, aromaticity) on one of its atoms.
#      Exception: an angle whose three atoms lie in one 3-, 4- or 5-membered
#      ring is constrained by the ring closure (DYFF's natural angles are
#      109.5/120 degrees): it is judged by its deviation from the polygon
#      angle (60/90/108 degrees, tolerance 10) instead of by its energy; so
#      are the exocyclic angles at a planar ring atom, ( 360 - polygon ) / 2.
#      Atoms: untyped -> red, manual type -> blue, every DYFF charge being
#      zero -> yellow ("no partial charges": no electrostatics).
#
#      The parameters are read from the MM model pDynamo actually builds
#      (mmState.mmTerms: parameter keys, values and the atoms of every term),
#      on a scratch pDynamo system made from the Builder object with the same
#      manual type overrides and guanidinium correction as
#      session.define_MMModel().
#
#      Partial charges: DYFF has none (all zero). compute_cm5_charges() puts
#      GFN1-xTB CM5 charges (pdynamo/opls/ligand.py) on atom.mm_charge; the
#      Builder keeps them through Undo and define_MMModel('DYFF') applies them.
#

import math
from collections import defaultdict, OrderedDict

import numpy as np

STRAIN_GREEN  = 2.0       # kJ/mol
RING_ANGLE_TOLERANCE = 10.0   # degrees, angle inside a 3/4/5-membered ring vs the polygon angle
STRAIN_YELLOW = 8.0       # kJ/mol
LEVEL_COLORS = { "green": "#4caf50", "yellow": "#f2c94c", "red": "#e5534b", "blue": "#4a90d9", "gray": "#9e9e9e" }
TERM_OF_CONTAINER = { "HarmonicBondContainer": "bond", "CosineAngleContainer": "angle",
                      "CosineDihedralContainer": "dihedral", "CosineOutOfPlaneContainer": "outofplane" }


#=============================================================================
#  Scratch system with the DYFF model
#=============================================================================
_AS_OBJECT = object ( )

def soften_polar_hydrogens ( system ):
    """ pdynamo.dyff_edits.soften_polar_hydrogens (imported late: pDynamo order). """
    from pdynamo.dyff_edits import soften_polar_hydrogens as soften
    return soften ( system )

def build_dyff_system ( vismol_object, parameter_set = "dyff-1.0", edits = _AS_OBJECT ):
    """ ( system or None, error text ). Same typing as define_MMModel('DYFF');
        the parameter edits are the object's unless given. """
    from gui.windows.builder.empty_object import _build_pdynamo_system_from_vismol_object
    from pdynamo.pDynamo2EasyHybrid.session import ( _apply_dyff_guanidinium_correction,
                                                     _apply_manual_atom_type_overrides )
    from pdynamo.dyff_edits import dyff_model
    system = _build_pdynamo_system_from_vismol_object ( vismol_object, label = vismol_object.name )
    overrides = getattr ( vismol_object, "manual_atom_type_overrides", None )
    added = [ ]
    try:
        _apply_dyff_guanidinium_correction ( system )
        added = _apply_manual_atom_type_overrides ( system, overrides )
        if edits is _AS_OBJECT: edits = getattr ( vismol_object, "dyff_parameter_edits", None )
        from pdynamo.dyff_edits import applicable_edits
        edits, _ = applicable_edits ( edits, len ( vismol_object.atoms ), getattr ( vismol_object, "manual_bonds", None ) or ( ) )
        system.DefineMMModel ( dyff_model ( parameter_set, edits ), log = None )
    except Exception as error:
        return system, str ( error ).split ( "\n" )[0]
    finally:
        if added:
            from pBabel import MOL2AtomTypesToDYFF
            for key in added: MOL2AtomTypesToDYFF.TriposMOL2AtomTypes.pop ( key, None )
    return system, None


#=============================================================================
#  Strain of one term at a geometry
#=============================================================================
def _angle ( a, b, c ):
    u, v = a - b, c - b
    return math.acos ( max ( -1.0, min ( 1.0, float ( np.dot ( u, v ) / ( np.linalg.norm ( u ) * np.linalg.norm ( v ) ) ) ) ) )

def _dihedral ( a, b, c, d ):
    b0, b1, b2 = a - b, c - b, d - c
    b1n = b1 / np.linalg.norm ( b1 )
    v = b0 - np.dot ( b0, b1n ) * b1n
    w = b2 - np.dot ( b2, b1n ) * b1n
    return math.atan2 ( float ( np.dot ( np.cross ( b1n, v ), w ) ), float ( np.dot ( v, w ) ) )

def _oop_gamma ( i, j, k, l ):
    """ Angle between the j-i bond and the normal of the j-k-l plane
        (pDynamo/DYFF out-of-plane convention, j central). """
    normal = np.cross ( k - j, l - j )
    u = i - j
    c = float ( np.dot ( u, normal ) / ( np.linalg.norm ( u ) * np.linalg.norm ( normal ) ) )
    return math.acos ( max ( -1.0, min ( 1.0, c ) ) )

# . angles: search the minimum in [60, 180] degrees -- a trigonal centre's
#   c0 ( 1 - cos 3 theta ) also has a (meaningless) minimum at 0 degrees
_GRID = np.linspace ( math.pi / 3.0, math.pi, 481 )
_GRID_FULL = np.linspace ( -math.pi, math.pi, 1441 )

def _cosine_energy ( pairs, x ):
    return sum ( fc * np.cos ( n * x ) for ( fc, n ) in pairs )

def term_strain ( term, parameters, points ):
    """ ( internal coordinate, strain in kJ/mol ) of one term. """
    p = [ np.asarray ( x, dtype = float ) for x in points ]
    if term == "bond":
        b0, k = parameters
        b = float ( np.linalg.norm ( p[0] - p[1] ) )
        return b, k * ( b - b0 ) ** 2
    if term == "angle":
        x = _angle ( p[0], p[1], p[2] )
        energy = _cosine_energy ( parameters, x ) - float ( np.min ( _cosine_energy ( parameters, _GRID ) ) )
        return math.degrees ( x ), float ( energy )
    if term == "dihedral":
        x = _dihedral ( *p )
        energy = _cosine_energy ( parameters, x ) - float ( np.min ( _cosine_energy ( parameters, _GRID_FULL ) ) )
        return math.degrees ( x ), float ( energy )
    if term == "outofplane":
        x = _oop_gamma ( *p )
        energy = _cosine_energy ( parameters, x ) - float ( np.min ( _cosine_energy ( parameters, _GRID ) ) )
        return math.degrees ( x ), float ( energy )
    return None, 0.0

def strain_level ( energy ):
    if energy < STRAIN_GREEN:  return "green"
    if energy < STRAIN_YELLOW: return "yellow"
    return "red"

def equilibrium_text ( term, parameters ):
    """ Readable parameters: b0/k for bonds, minimum of the cosine form for
        angles, the Fourier pairs otherwise. """
    if term == "bond":
        return "{:.4f}".format ( parameters[0] ), "{:.1f}".format ( parameters[1] )
    if term == "angle":
        minimum = _GRID[int ( np.argmin ( _cosine_energy ( parameters, _GRID ) ) )]
        return "{:.1f}".format ( math.degrees ( minimum ) ), "; ".join ( "{:.2f} {:d}".format ( fc, int ( n ) ) for fc, n in parameters )
    return "; ".join ( "{:.3f} {:d}".format ( fc, int ( n ) ) for fc, n in parameters ), ""


#=============================================================================
#  Tables
#=============================================================================
def _hybridizations ( vismol_object ):
    try:
        from gui.windows.builder.residue_chemistry import perceive_aromatic_bonds, atom_hybridizations
        obj = vismol_object
        orders = { }
        manual_orders = getattr ( obj, "manual_bond_orders", None ) or { }
        for pair in ( getattr ( obj, "manual_bonds", None ) or set ( ) ):
            if pair[0] in obj.atoms and pair[1] in obj.atoms:
                orders[pair] = int ( manual_orders.get ( pair, 1 ) )
        charges = { aid: int ( getattr ( a, "formal_charge", 0 ) ) for aid, a in obj.atoms.items ( ) }
        aromatic = perceive_aromatic_bonds ( obj.atoms, orders, charges )
        return atom_hybridizations ( obj.atoms, orders, aromatic )
    except Exception:
        return { }

def _small_rings ( adjacency, max_size = 5 ):
    """ Rings of 3-5 atoms (as sets of atom ids), by shortest cycles. """
    rings = set ( )
    for start in adjacency:
        for first in adjacency[start]:
            parent, queue, found = { first: None }, [ first ], None
            while queue and found is None:
                nxt = [ ]
                for u in queue:
                    for v in adjacency[u]:
                        if u == first and v == start: continue
                        if v == start: found = u; break
                        if v not in parent: parent[v] = u; nxt.append ( v )
                    if found is not None: break
                queue = nxt
                if len ( parent ) > 64: break
            if found is not None:
                path, u = [ start ], found
                while u is not None: path.append ( u ); u = parent[u]
                if 3 <= len ( path ) <= max_size: rings.add ( frozenset ( path ) )
    return [ set ( r ) for r in rings ]

def _key_text ( key ):
    """ "C:Tet:C:Res:S" -> "C:Tet - C:Res (S)": DYFF keys are type labels
        joined with ':' plus bond tags; shown as the atom types of the term. """
    return key

def build_tables ( vismol_object ):
    """ { "error", "atoms": [ rows ], "bond"|"angle"|"dihedral"|"outofplane": [ rows ],
          "lj": [ rows ], "total_formal", "total_partial", "has_charges" }. """
    from gui.windows.builder.atom_types import compute_atom_types
    tables = { "error": None }
    system, error = build_dyff_system ( vismol_object )
    tables["error"] = error
    types_by_id = compute_atom_types ( vismol_object )
    overrides = getattr ( vismol_object, "manual_atom_type_overrides", None ) or { }
    hybrid = _hybridizations ( vismol_object )
    xyz = { aid: np.asarray ( vismol_object.frames[0, aid], dtype = float ) for aid in vismol_object.atoms }
    has_charges = all ( getattr ( a, "mm_charge", None ) is not None for a in vismol_object.atoms.values ( ) )
    # . small rings (angles constrained by ring closure)
    adjacency = defaultdict ( set )
    for ( a, b ) in ( getattr ( vismol_object, "manual_bonds", None ) or set ( ) ):
        adjacency[a].add ( b ); adjacency[b].add ( a )
    small_rings = _small_rings ( adjacency )
    def ring_angle ( inst ):
        """ ( expected angle, text ) for an angle shaped by a small ring:
            inside the ring -> polygon angle; at a planar ring atom (3
            neighbours) between a ring atom and an outside atom -> the
            exocyclic angle ( 360 - polygon ) / 2; else None. """
        i, j, k = inst
        for ring in small_rings:
            polygon = 180.0 * ( len ( ring ) - 2 ) / len ( ring )
            if { i, j, k } <= ring:
                return polygon, "{}-ring angle (polygon {:.0f} deg)".format ( len ( ring ), polygon )
        for ring in small_rings:
            polygon = 180.0 * ( len ( ring ) - 2 ) / len ( ring )
            if j in ring and len ( adjacency[j] ) == 3 and ( ( i in ring ) != ( k in ring ) ):
                expected = ( 360.0 - polygon ) / 2.0
                return expected, "exocyclic angle at a {}-ring atom ({:.0f} deg)".format ( len ( ring ), expected )
        return None
    # . terms
    atom_worst = defaultdict ( lambda: ( 0.0, None ) )
    for term in ( "bond", "angle", "dihedral", "outofplane" ):
        tables[term] = [ ]
    edits = getattr ( vismol_object, "dyff_parameter_edits", None ) or { }
    if error is None:
        for container in system.mmState.mmTerms:
            term = TERM_OF_CONTAINER.get ( type ( container ).__name__ )
            if term is None: continue
            state = container.__getstate__ ( )
            keys, parameters = state["parameterKeys"], state["parameters"]
            grouped = defaultdict ( list )
            for row in state["terms"]:
                if term == "bond": indices, p = ( row[0], row[1] ), row[2]
                else:              indices, p = tuple ( row[0] ), row[1]
                grouped[p].append ( tuple ( int ( i ) for i in indices ) )
            for p, instances in grouped.items ( ):
                values = parameters[p]
                values = tuple ( values ) if term == "bond" else [ ( float ( fc ), int ( n ) ) for fc, n in values ]
                worst, worst_value, ring_note = 0.0, None, None
                for inst in instances:
                    measured, energy = term_strain ( term, values, [ xyz[i] for i in inst ] )
                    ring = ring_angle ( inst ) if term == "angle" else None
                    if ring is not None:
                        expected, ring_note = ring
                        deviation = abs ( measured - expected )
                        energy = 0.0 if deviation <= RING_ANGLE_TOLERANCE else STRAIN_YELLOW * deviation / RING_ANGLE_TOLERANCE
                    if energy >= worst: worst, worst_value = energy, measured
                    for i in inst:
                        if energy > atom_worst[i][0]: atom_worst[i] = ( energy, term )
                v0, v1 = equilibrium_text ( term, values )
                example = "-".join ( vismol_object.atoms[i].name for i in instances[0] )
                edited = keys[p] in ( edits.get ( term ) or { } ) or "term " in keys[p] or "bond " in keys[p]
                tables[term].append ( { "key": keys[p], "types": _key_text ( keys[p] ), "value": values, "v0": v0, "v1": v1,
                                        "edited": edited,
                                        "instances": instances, "count": len ( instances ), "example": example,
                                        "strain": worst, "measured": worst_value, "level": strain_level ( worst ),
                                        "quality": ( "{} ok".format ( ring_note ) if ring_note and worst == 0.0
                                                     else "{:.2f} kJ/mol".format ( worst ) + ( "  [{}]".format ( ring_note ) if ring_note else "" ) ) } )
            for row in tables[term]:
                if row["edited"]:
                    row["level"], row["quality"] = "blue", "edited -- " + row["quality"]
            tables[term].sort ( key = lambda r: -r["strain"] )
    # . atoms
    rows = [ ]
    for aid in sorted ( vismol_object.atoms ):
        atom = vismol_object.atoms[aid]
        perceived, effective = types_by_id.get ( aid, ( None, None ) )
        q = getattr ( atom, "mm_charge", None )
        worst, worst_term = atom_worst[aid]
        if effective is None:            level, text = "red", "untyped"
        elif aid in overrides:           level, text = "blue", "manual type"
        else:
            level = strain_level ( worst )
            text = "max strain {:.1f} kJ/mol ({})".format ( worst, worst_term ) if worst_term else "ok"
            if not has_charges and level == "green": level, text = "yellow", "no partial charges"
        residue = "{}{}".format ( atom.residue.name, atom.residue.index ) if atom.residue is not None else ""
        rows.append ( { "atom_id": aid, "symbol": atom.symbol, "name": atom.name, "residue": residue,
                        "type": effective or "?", "perceived": perceived or "?", "override": aid in overrides,
                        "hybridization": hybrid.get ( aid, "" ), "formal": int ( getattr ( atom, "formal_charge", 0 ) or 0 ),
                        "charge": q, "instances": [ ( aid, ) ], "level": level, "quality": text } )
    tables["atoms"] = rows
    # . Lennard-Jones of the types present
    tables["lj"] = [ ]
    if error is None:
        state = system.mmState
        types = [ state.atomTypes[i] for i in state.atomTypeIndices ]
        lj = _dyff_lj ( )
        for label in sorted ( set ( types ) ):
            members = [ i for i, t in enumerate ( types ) if t == label ]
            e, s = lj.get ( label, ( None, None ) )
            tables["lj"].append ( { "key": label, "types": label, "v0": "" if e is None else "{:.4f}".format ( e ),
                                    "v1": "" if s is None else "{:.4f}".format ( s ), "instances": [ ( i, ) for i in members ],
                                    "count": len ( members ), "example": label, "level": "green" if e is not None else "gray",
                                    "quality": "DYFF" if e is not None else "no LJ" } )
    tables["total_formal"] = sum ( r["formal"] for r in rows )
    tables["total_partial"] = sum ( r["charge"] for r in rows if r["charge"] is not None )
    tables["has_charges"] = has_charges
    return tables

_LJ = None

def _dyff_lj ( ):
    global _LJ
    if _LJ is None:
        import os, yaml
        from pdynamo.dyff_edits import dyff_parameter_dir   # EasyHybrid extended set (dyff-1.0-eh)
        path = os.path.join ( dyff_parameter_dir ( "dyff-1.0" ), "lennardJonesParameters.yaml" )
        with open ( path ) as handle: data = yaml.safe_load ( handle )
        fields = data["Parameter Fields"]
        ie, isg = fields.index ( "Epsilon" ), fields.index ( "Sigma" )
        _LJ = { str ( r[0] ): ( float ( r[ie] ), float ( r[isg] ) ) for r in data["Parameter Values"] }
    return _LJ


#=============================================================================
#  Partial charges
#=============================================================================
def compute_cm5_charges ( vismol_object, total_charge = None, multiplicity = 1, scale = None, xtb = None ):
    """ GFN1-xTB CM5 charges for the whole Builder object (the same
        protocol as the OPLS ligands: symmetrized, exact total) stored on
        atom.mm_charge. total_charge defaults to the sum of the formal
        charges. Returns the list of charges. """
    from pdynamo.opls import ligand as L
    from gui.windows.builder.atom_ops import push_undo_snapshot
    ids = sorted ( vismol_object.atoms )
    lig = L.LigandParameters ( "MOL" )
    lig.atoms = [ ( "{}{}".format ( vismol_object.atoms[i].symbol, k + 1 ), L.ELEMENT_Z.get ( vismol_object.atoms[i].symbol, 0 ),
                    np.asarray ( vismol_object.frames[0, i], dtype = float ) ) for k, i in enumerate ( ids ) ]
    if any ( z == 0 for ( _, z, _ ) in lig.atoms ):
        raise ValueError ( "CM5 charges: elements without xTB/CM5 support in this molecule." )
    lig.formal_charges = [ int ( getattr ( vismol_object.atoms[i], "formal_charge", 0 ) or 0 ) for i in ids ]
    position = { i: k for k, i in enumerate ( ids ) }
    orders = getattr ( vismol_object, "manual_bond_orders", None ) or { }
    lig.bonds = [ ( position[a], position[b], int ( orders.get ( ( a, b ), 1 ) ) )
                  for ( a, b ) in ( getattr ( vismol_object, "manual_bonds", None ) or set ( ) ) if a in position and b in position ]
    from pdynamo.opls.ligand import equivalence_classes, _neighbors, Z_ELEMENT
    lig.equivalence = equivalence_classes ( [ Z_ELEMENT[z] for ( _, z, _ ) in lig.atoms ], _neighbors ( len ( lig.atoms ), lig.bonds ) )
    L.apply_user_charge ( lig, total_charge, multiplicity )
    charges = L.cm5_charges ( lig, xtb = xtb, scale = scale )
    push_undo_snapshot ( vismol_object )
    for i, q in zip ( ids, charges ):
        vismol_object.atoms[i].mm_charge = float ( q )
    return charges

def clear_partial_charges ( vismol_object ):
    from gui.windows.builder.atom_ops import push_undo_snapshot
    push_undo_snapshot ( vismol_object )
    for atom in vismol_object.atoms.values ( ):
        atom.mm_charge = None

def set_partial_charge ( vismol_object, atom_id, value ):
    """ One partial charge (atoms without one get 0.0 so that the set is
        complete). The total is NOT forced here; the window shows it. """
    for atom in vismol_object.atoms.values ( ):
        if getattr ( atom, "mm_charge", None ) is None: atom.mm_charge = 0.0
    vismol_object.atoms[atom_id].mm_charge = float ( value )


#=============================================================================
#  Parameter edits (pdynamo/dyff_edits.py)
#=============================================================================
def parse_pairs ( text ):
    """ "c n; c n" -> [ ( c, n ) ]. """
    pairs = [ ]
    for chunk in text.replace ( ",", " " ).split ( ";" ):
        parts = chunk.split ( )
        if not parts: continue
        if len ( parts ) != 2: raise ValueError ( 'Use "force constant period; ..." (e.g. "13.07 0; -13.07 2").' )
        pairs.append ( ( float ( parts[0] ), int ( parts[1] ) ) )
    if not pairs: raise ValueError ( "No terms given." )
    return pairs

def set_parameter_edit ( vismol_object, term, key, value ):
    """ Stores ( b0, k ) for a bond or [ ( c, n ) ] for the cosine terms;
        value None removes the edit (back to the DYFF rules). One undo step. """
    from gui.windows.builder.atom_ops import push_undo_snapshot
    push_undo_snapshot ( vismol_object )
    edits = getattr ( vismol_object, "dyff_parameter_edits", None )
    if edits is None:
        edits = vismol_object.dyff_parameter_edits = { }
    table = edits.setdefault ( term, { } )
    if value is None: table.pop ( key, None )
    else:             table[key] = [ list ( v ) for v in value ] if term != "bond" else list ( value )


#=============================================================================
#  Torsion fitting (pdynamo/torsion_fit.py, gui/windows/builder/torsion_fit_window.py)
#=============================================================================
def scan_inputs ( vismol_object ):
    """ ( symbols, coordinates N x 3, bonds [ ( a, b ) ] ) in system index order. """
    ids = sorted ( vismol_object.atoms )
    symbols = [ vismol_object.atoms[i].symbol for i in ids ]
    coordinates = np.array ( [ vismol_object.frames[0, i] for i in ids ], dtype = float )
    bonds = [ tuple ( pair ) for pair in ( getattr ( vismol_object, "manual_bonds", None ) or set ( ) ) ]
    return symbols, coordinates, bonds

def make_fit_system_factory ( vismol_object, term, key ):
    """ make_system ( value ) for DYFFTorsionTarget: the
        DYFF system of the object (its edits and partial charges, all pairs
        nonbonded interactions: gas phase, like the xTB reference) with the
        key set to value ( None = as now ). """
    import copy
    from pMolecule.NBModel import NBModelFull
    base = copy.deepcopy ( getattr ( vismol_object, "dyff_parameter_edits", None ) or { } )
    charges = [ getattr ( vismol_object.atoms[i], "mm_charge", None ) for i in sorted ( vismol_object.atoms ) ]
    def make_system ( value ):
        edits = copy.deepcopy ( base )
        if value is not None:
            edits.setdefault ( term, { } )[key] = [ list ( v ) for v in value ]
        system, error = build_dyff_system ( vismol_object, edits = edits )
        if error: raise RuntimeError ( "DYFF model not built: " + error )
        if all ( q is not None for q in charges ):
            for i, q in enumerate ( charges ): system.mmState.charges[i] = float ( q )
            soften_polar_hydrogens ( system )                     # as define_MMModel does with charges
        system.DefineNBModel ( NBModelFull.WithDefaults ( ) )
        return system
    return make_system


class DYFFTorsionTarget:
    """ What the torsion fit window needs from DYFF: the Builder molecule,
        the dihedral parameter of a Dihedrals row and how to apply a fit. """

    force_field = "DYFF"
    units_note  = "c n (kJ/mol)"

    def __init__ ( self, vismol_object, row, quad ):
        self.vismol_object, self.row = vismol_object, row
        self.key, self.key_text = row["key"], row["key"]
        self.symbols, self.coordinates, self.bonds = scan_inputs ( vismol_object )
        self.quad = tuple ( int ( i ) for i in quad )
        self.instances = [ tuple ( q ) for q in row["instances"] ]
        self.n_atoms = len ( vismol_object.atoms )
        self.current_text = row["v0"]
        self.atom_names = [ vismol_object.atoms[i].name for i in self.quad ]
        self.default_charge = sum ( int ( getattr ( a, "formal_charge", 0 ) or 0 ) for a in vismol_object.atoms.values ( ) )
        self._make_system = make_fit_system_factory ( vismol_object, "dihedral", self.key )

    def notes ( self ):
        notes = [ ]
        if not all ( getattr ( a, "mm_charge", None ) is not None for a in self.vismol_object.atoms.values ( ) ):
            notes.append ( "No partial charges: DYFF electrostatics are zero, so the fitted torsion will also absorb "
                           "the missing electrostatics. Compute CM5 charges first (Atoms tab) for a transferable parameter." )
        return notes

    def energies_for ( self, pairs, geometries ):
        from pdynamo.torsion_fit import mm_energies
        return mm_energies ( self._make_system ( pairs ), geometries )

    def format_value ( self, pairs ):
        from pdynamo.torsion_fit import format_cosine_pairs
        return format_cosine_pairs ( pairs )

    def changed ( self ):
        if len ( self.vismol_object.atoms ) != self.n_atoms:
            return "The molecule changed during the scan: run it again."
        return None

    def apply ( self, pairs ):
        set_parameter_edit ( self.vismol_object, "dihedral", self.key, pairs )
        return "Fitted parameter applied to {} term(s) (Undo reverts it).".format ( len ( self.instances ) )


#=============================================================================
#  Closing the Builder: incorporate the system with DYFF (builder_sidebar)
#=============================================================================
def builder_customizations ( vismol_object ):
    """ What the Builder session customized for DYFF: { "edits": { term: n },
        "charges": bool, "total_partial", "overrides": n, "formal": n atoms
        with a formal charge }. """
    edits = getattr ( vismol_object, "dyff_parameter_edits", None ) or { }
    atoms = list ( vismol_object.atoms.values ( ) )
    charges = bool ( atoms ) and all ( getattr ( a, "mm_charge", None ) is not None for a in atoms )
    return { "edits": { term: len ( table ) for term, table in edits.items ( ) if table },
             "charges": charges,
             "total_partial": sum ( a.mm_charge for a in atoms ) if charges else 0.0,
             "overrides": len ( getattr ( vismol_object, "manual_atom_type_overrides", None ) or { } ),
             "formal": sum ( 1 for a in atoms if int ( getattr ( a, "formal_charge", 0 ) or 0 ) ) }

def customizations_text ( summary ):
    names = { "bond": "bond", "angle": "angle", "dihedral": "dihedral", "outofplane": "improper" }
    lines = [ ]
    for term, n in summary["edits"].items ( ):
        lines.append ( "{} edited {} parameter(s)".format ( n, names.get ( term, term ) ) )
    lines.append ( "partial charges: {}".format ( "set (total {:+.3f})".format ( summary["total_partial"] )
                                                  if summary["charges"] else "none (DYFF: all zero)" ) )
    if summary["overrides"]: lines.append ( "{} manual atom type(s)".format ( summary["overrides"] ) )
    if summary["formal"]:    lines.append ( "{} atom(s) with a formal charge".format ( summary["formal"] ) )
    return lines

def carry_builder_customizations ( source, target ):
    """ Copies the DYFF customizations of a Builder clone onto the object it
        is folded back into (same atom count and order): edited parameters,
        type overrides, formal and partial charges. """
    import copy
    target.dyff_parameter_edits = copy.deepcopy ( getattr ( source, "dyff_parameter_edits", None ) or { } )
    target.manual_atom_type_overrides = dict ( getattr ( source, "manual_atom_type_overrides", None ) or { } )
    for i, atom in source.atoms.items ( ):
        if i in target.atoms:
            target.atoms[i].formal_charge = int ( getattr ( atom, "formal_charge", 0 ) or 0 )
            target.atoms[i].mm_charge = getattr ( atom, "mm_charge", None )

def incorporate_with_dyff ( main, vismol_object ):
    """ Assigns DYFF, with the Builder customizations of vismol_object, to its
        pDynamo system (session.define_MMModel picks the edits, overrides and
        partial charges up from the vobject by e_id). ( ok, message ). """
    p_session = main.p_session
    system = p_session.psystem.get ( getattr ( vismol_object, "e_id", None ) )
    if system is None:
        return False, "No pDynamo system for '{}'.".format ( vismol_object.name )
    if getattr ( system, "qcModel", None ) is not None:
        return False, "'{}' has a QC model: DYFF was not assigned.".format ( system.label )
    ok, message = p_session.define_MMModel ( force_field = "DYFF", system = system )
    if ok and getattr ( vismol_object, "dyff_parameter_edits", None ):
        n = sum ( len ( t ) for t in vismol_object.dyff_parameter_edits.values ( ) )
        if n: message += " ({} edited DYFF parameter(s).)".format ( n )
    return ok, message


#=============================================================================
#  Bonds and angles from an xTB Hessian (modified Seminario)
#=============================================================================
def molecule_arrays ( vismol_object ):
    """ ( atomic numbers, bonds [ ( a, b, order, aromatic ) ], formal charges, coordinates )
        in system index order (atoms sorted by id). """
    ids = sorted ( vismol_object.atoms )
    index = { aid: k for k, aid in enumerate ( ids ) }
    orders = getattr ( vismol_object, "manual_bond_orders", None ) or { }
    aromatic = getattr ( vismol_object, "manual_aromatic_bonds", None ) or set ( )
    bonds = [ ( index[a], index[b], int ( orders.get ( ( a, b ), orders.get ( ( b, a ), 1 ) ) ), ( a, b ) in aromatic or ( b, a ) in aromatic )
              for ( a, b ) in ( getattr ( vismol_object, "manual_bonds", None ) or set ( ) ) if a in index and b in index ]
    z = [ _atomic_number ( vismol_object.atoms[aid].symbol ) for aid in ids ]
    formal = [ int ( getattr ( vismol_object.atoms[aid], "formal_charge", 0 ) or 0 ) for aid in ids ]
    xyz = np.array ( [ vismol_object.frames[0, aid] for aid in ids ], dtype = float )
    return z, bonds, formal, xyz

def hessian_calibration ( vismol_object ):
    """ pdynamo/bond_library.HessianCalibration of the molecule (DFT/xTB correction
        of the Seminario constants by environment), or None (no table, no RDKit,
        unsanitizable molecule). """
    try:
        from pdynamo import bond_library as BL
        if BL.load_hessian_library ( ) is None: return None
        return BL.HessianCalibration ( *molecule_arrays ( vismol_object ) )
    except Exception:
        return None

def hessian_bonds_angles ( vismol_object, total_charge = None, multiplicity = 1, adopt_geometry = False, threads = 4,
                           solvent = None, progress = None ):
    """ xtb --ohess on the Builder molecule (its own atoms and bonds: no
        typing, metals included), modified Seminario on the Hessian
        (pdynamo/seminario.py) and per-term DYFF edits "bond_terms" /
        "angle_terms" (b0 and theta0 from the xTB minimum), tied to the
        molecule's topology (term_signature). One Undo step. adopt_geometry:
        the Builder coordinates become the xTB minimum.
        The xTB part can run in a worker thread (compute_hessian_edits);
        store_hessian_edits() must run in the GTK thread. """
    result = compute_hessian_edits ( vismol_object, total_charge, multiplicity, threads, solvent )
    return store_hessian_edits ( vismol_object, result, adopt_geometry )

def compute_hessian_edits ( vismol_object, total_charge = None, multiplicity = 1, threads = 4, solvent = None, calibrate = True ):
    from pdynamo import seminario as SM
    symbols, coordinates, bonds = scan_inputs ( vismol_object )
    if total_charge is None:
        total_charge = sum ( int ( getattr ( a, "formal_charge", 0 ) or 0 ) for a in vismol_object.atoms.values ( ) )
    xyz, hessian, freqs = SM.xtb_hessian ( symbols, coordinates, charge = total_charge, multiplicity = multiplicity,
                                           threads = threads, solvent = solvent )
    skipped = [ ]
    calibration = hessian_calibration ( vismol_object ) if calibrate else None
    edits = SM.dyff_edits_from_hessian ( hessian, xyz, bonds, skipped, calibration = calibration )
    n_imag = int ( sum ( 1 for f in freqs if f < -50.0 ) )
    return { "edits": edits, "xyz": xyz, "n_atoms": len ( symbols ), "bonds": bonds, "n_imaginary": n_imag,
             "frequencies": freqs, "skipped": skipped, "solvent": solvent,
             "calibrated": dict ( calibration.used ) if calibration is not None else { } }

def store_hessian_edits ( vismol_object, result, adopt_geometry = False ):
    from gui.windows.builder.atom_ops import push_undo_snapshot
    from pdynamo.dyff_edits import term_signature
    if len ( vismol_object.atoms ) != result["n_atoms"]:
        raise ValueError ( "The molecule changed while xTB was running: run it again." )
    push_undo_snapshot ( vismol_object )
    edits = getattr ( vismol_object, "dyff_parameter_edits", None )
    if edits is None:
        edits = vismol_object.dyff_parameter_edits = { }
    edits["bond_terms"] = result["edits"]["bond_terms"]
    edits["angle_terms"] = result["edits"]["angle_terms"]
    edits["term_signature"] = term_signature ( result["n_atoms"], result["bonds"] )
    if adopt_geometry:
        for i, aid in enumerate ( sorted ( vismol_object.atoms ) ):
            vismol_object.frames[0, aid] = result["xyz"][i]
        from gui.windows.builder.empty_object import sync_pdynamo_system
        sync_pdynamo_system ( vismol_object )
    return { "bonds": len ( edits["bond_terms"] ), "angles": len ( edits["angle_terms"] ),
             "n_imaginary": result["n_imaginary"], "skipped": result.get ( "skipped", [ ] ) }

def clear_hessian_edits ( vismol_object ):
    from gui.windows.builder.atom_ops import push_undo_snapshot
    edits = getattr ( vismol_object, "dyff_parameter_edits", None ) or { }
    if not ( edits.get ( "bond_terms" ) or edits.get ( "angle_terms" ) or edits.get ( "dihedral_bonds" ) ):
        return False
    push_undo_snapshot ( vismol_object )
    for key in ( "bond_terms", "angle_terms", "dihedral_bonds", "term_signature" ): edits.pop ( key, None )
    return True



#=============================================================================
#  Rotatable-bond torsions fitted to xTB scans of the molecule itself
#=============================================================================
def rotatable_bonds ( vismol_object ):
    """ [ ( ( j, k ), quad ) ] -- single, acyclic bonds between two atoms that
        both have other neighbours; quad: the heavy-atom dihedral i-j-k-l with
        the heaviest outer atoms (deterministic). Indices = system order. """
    ids = sorted ( vismol_object.atoms )
    index = { aid: k for k, aid in enumerate ( ids ) }
    orders = getattr ( vismol_object, "manual_bond_orders", None ) or { }
    aromatic = getattr ( vismol_object, "manual_aromatic_bonds", None ) or set ( )
    bonds = [ ( index[a], index[b] ) for ( a, b ) in ( getattr ( vismol_object, "manual_bonds", None ) or set ( ) )
              if a in index and b in index ]
    neighbours = defaultdict ( set )
    for a, b in bonds:
        neighbours[a].add ( b ); neighbours[b].add ( a )
    z = [ _atomic_number ( vismol_object.atoms[aid].symbol ) for aid in ids ]
    def in_ring ( a, b ):
        seen, stack = { a }, [ a ]
        while stack:
            u = stack.pop ( )
            for v in neighbours[u]:
                if ( u, v ) in ( ( a, b ), ( b, a ) ): continue
                if v == b: return True
                if v not in seen: seen.add ( v ); stack.append ( v )
        return False
    out = [ ]
    for ( a, b ) in ( getattr ( vismol_object, "manual_bonds", None ) or set ( ) ):
        if a not in index or b not in index: continue
        if int ( orders.get ( ( a, b ), orders.get ( ( b, a ), 1 ) ) ) != 1 or ( a, b ) in aromatic or ( b, a ) in aromatic: continue
        j, k = index[a], index[b]
        if len ( neighbours[j] ) < 2 or len ( neighbours[k] ) < 2 or in_ring ( j, k ): continue
        outer_j = [ n for n in neighbours[j] if n != k ]
        outer_k = [ n for n in neighbours[k] if n != j ]
        i = max ( outer_j, key = lambda n: ( z[n], -n ) )
        l = max ( outer_k, key = lambda n: ( z[n], -n ) )
        if z[i] == 1 and z[l] == 1 and all ( z[n] == 1 for n in outer_j ) and all ( z[n] == 1 for n in outer_k ):
            continue                                       # e.g. H3C-CH3: nothing but hydrogens around
        out.append ( ( ( min ( j, k ), max ( j, k ) ), ( i, j, k, l ) ) )
    return sorted ( out )


def _atomic_number ( symbol ):
    from pScientific import PeriodicTable
    return PeriodicTable.AtomicNumber ( symbol )


def system_factory ( vismol_object ):
    """ make_system ( edits ) for pdynamo/dyff_bond_fit.py: the object's DYFF system
        with these edits, its partial charges and all-pairs nonbonded (gas phase). """
    import pBabel                                    # noqa: F401 -- pDynamo import order
    from pMolecule.NBModel import NBModelFull
    charges = [ getattr ( vismol_object.atoms[i], "mm_charge", None ) for i in sorted ( vismol_object.atoms ) ]
    def make ( edits ):
        system, error = build_dyff_system ( vismol_object, edits = edits )
        if error: raise RuntimeError ( "DYFF model not built: " + error )
        if all ( q is not None for q in charges ):
            for i, q in enumerate ( charges ): system.mmState.charges[i] = float ( q )
            soften_polar_hydrogens ( system )                     # as define_MMModel does with charges
        system.DefineNBModel ( NBModelFull.WithDefaults ( ) )
        return system
    return make


def torsion_library_available ( ):
    from pdynamo import bond_library as BL
    if BL.RDKIT_CODES:                                    # . RDKit codes only when asked for (EASYHYBRID_ENV_CODES=rdkit)
        try:
            import rdkit                                  # noqa: F401
        except ImportError:
            return False
    return BL.load_library ( ) is not None


def library_torsion_results ( vismol_object, targets, say = None ):
    """ Per-bond torsions from the DFT library (pdynamo/bond_library.py) for the
        targets ( ( j, k ), quad ) (system indices) whose bond environment is in
        it. Returns ( results, remaining targets ); results like _scan_and_fit's,
        with "source" instead of the fit statistics. """
    from pdynamo import bond_library as BL
    ids = sorted ( vismol_object.atoms )
    # . the table fitted on the same model: with partial charges (CM5, residue library) or without
    charged = all ( getattr ( vismol_object.atoms[aid], "mm_charge", None ) is not None for aid in ids )
    library = BL.load_library ( charged = charged )
    if library is None and charged:
        library = BL.load_library ( )
        if library is not None and say: say ( "no torsion table for charged molecules: using the uncharged one" )
    if library is None or not targets: return [ ], list ( targets )
    index = { aid: k for k, aid in enumerate ( ids ) }
    orders = getattr ( vismol_object, "manual_bond_orders", None ) or { }
    aromatic = getattr ( vismol_object, "manual_aromatic_bonds", None ) or set ( )
    bonds = [ ( index[a], index[b], int ( orders.get ( ( a, b ), orders.get ( ( b, a ), 1 ) ) ), ( a, b ) in aromatic or ( b, a ) in aromatic )
              for ( a, b ) in ( getattr ( vismol_object, "manual_bonds", None ) or set ( ) ) if a in index and b in index ]
    z = [ _atomic_number ( vismol_object.atoms[aid].symbol ) for aid in ids ]
    formal = [ int ( getattr ( vismol_object.atoms[aid], "formal_charge", 0 ) or 0 ) for aid in ids ]
    xyz = np.array ( [ vismol_object.frames[0, aid] for aid in ids ], dtype = float )
    try:
        found = BL.library_torsions ( z, bonds, formal, [ bond for bond, quad in targets ], coordinates = xyz, library = library )
    except Exception as error:                                   # e.g. RDKit cannot sanitize the molecule
        if say: say ( "torsion library not used: {}".format ( str ( error ).split ( "\n" )[0][:120] ) )
        return [ ], list ( targets )
    results, remaining = [ ], [ ]
    for bond, quad in targets:
        if bond in found:
            level, pairs, n = found[bond]
            results.append ( { "bond": bond, "quad": quad, "pairs": pairs, "source": "library", "level": level, "n_scans": n,
                               "table": "CM5" if charged and library is BL.load_library ( charged = True ) else "uncharged" } )
        else:
            remaining.append ( ( bond, quad ) )
    if say: say ( "{} of {} rotatable bond(s) from the DFT torsion library".format ( len ( results ), len ( targets ) ) )
    return results, remaining


def refine_torsions ( vismol_object, total_charge = None, multiplicity = 1, solvent = None, step = 15.0, threads = 1,
                      workers = 4, base_edits = None, progress = None, cancel = None, use_library = True ):
    """ For every rotatable bond: the DFT torsion library when the bond's
        environment is in it (use_library, library_torsion_results), else a
        relaxed GFN2-xTB scan of the WHOLE molecule (pdynamo/torsion_fit.run_scan,
        scans run in parallel) and the self-consistent per-bond fit
        (pdynamo/dyff_bond_fit.fit_bond) on top of base_edits (e.g. the Hessian
        bonds/angles). Worker-thread safe (reads the object, writes nothing).
        Returns [ result per bond ]. """
    say = progress or ( lambda text: None )
    symbols, coordinates, bonds = scan_inputs ( vismol_object )
    ids = sorted ( vismol_object.atoms )
    index = { aid: k for k, aid in enumerate ( ids ) }
    bonds = [ ( index[a], index[b] ) for a, b in bonds ]
    if total_charge is None:
        total_charge = sum ( int ( getattr ( a, "formal_charge", 0 ) or 0 ) for a in vismol_object.atoms.values ( ) )
    targets = rotatable_bonds ( vismol_object )
    library = [ ]
    if use_library: library, targets = library_torsion_results ( vismol_object, targets, say )
    return library + _scan_and_fit ( symbols, coordinates, bonds, targets, system_factory ( vismol_object ), base_edits, total_charge,
                                     multiplicity, solvent, step, threads, workers, say, cancel )


def _scan_and_fit ( symbols, coordinates, bonds, targets, make, base_edits, charge, multiplicity, solvent, step, threads,
                    workers, say, cancel = None, report_bond = None ):
    """ Relaxed xTB scans of the targets ( ( bond, quad ) ) in parallel, then the
        self-consistent per-bond fit on the MM system given by make ( edits ).
        report_bond ( bond ) -> the bond as reported (e.g. mapped to parent indices). """
    import concurrent.futures
    from pdynamo import torsion_fit as TF
    from pdynamo import dyff_bond_fit as BF
    if not targets: return [ ]
    say ( "{} rotatable bond(s): xTB scans...".format ( len ( targets ) ) )
    def scan ( item ):
        bond, quad = item
        if cancel is not None and cancel ( ): return bond, quad, None
        runner = TF.XTBRunner ( gfn = 2, charge = charge, multiplicity = multiplicity, solvent = solvent, threads = threads )
        return bond, quad, TF.run_scan ( symbols, coordinates, bonds, quad, step = step, relaxed = True, runner = runner )
    scans = [ ]
    with concurrent.futures.ThreadPoolExecutor ( max_workers = max ( 1, workers ) ) as pool:
        for n, ( bond, quad, result ) in enumerate ( pool.map ( scan, targets ), 1 ):
            scans.append ( ( bond, quad, result ) )
            say ( "xTB scans {}/{}".format ( n, len ( targets ) ) )
    results = [ ]
    for n, ( bond, quad, scan_result ) in enumerate ( scans, 1 ):
        if scan_result is None: continue
        if cancel is not None and cancel ( ): break
        reported = report_bond ( bond ) if report_bond else bond
        say ( "fitting torsion {}/{} (bond {}-{})".format ( n, len ( scans ), *reported ) )
        e_qm = np.array ( scan_result["e_qm"], dtype = float )            # torsion_fit.run_scan: already kJ/mol
        broken = TF.broken_bonds ( symbols, bonds, scan_result["geometries"] )
        if broken:
            p, ( a, b ), d = broken[0]
            results.append ( { "bond": reported, "quad": quad, "error": "the xTB scan changed the bonding (atoms {}-{} at {:.2f} A: "
                               "proton transfer?) -- try implicit solvent (ALPB) or check charges/protonation".format ( a + 1, b + 1, d ) } )
            continue
        try:
            pairs, wrmse, dropped = BF.fit_bond ( make, base_edits, bonds, quad, scan_result["geometries"], scan_result["phi"], e_qm )
            before = BF.profile_error ( make, base_edits, quad, scan_result["geometries"], scan_result["phi"], e_qm )
            results.append ( { "bond": reported, "quad": quad, "pairs": pairs, "fit_rmse": wrmse / 4.184, "before_rmse": before["rmse"] / 4.184,
                               "barrier_xtb": float ( ( e_qm.max ( ) - e_qm.min ( ) ) / 4.184 ), "dropped": dropped } )
        except Exception as error:
            results.append ( { "bond": reported, "quad": quad, "error": str ( error ).split ( "\n" )[0][:160] } )
    return results


def store_torsion_edits ( vismol_object, results, hessian_result = None ):
    """ One Undo step: per-bond torsions (+ Hessian bonds/angles when given). """
    from gui.windows.builder.atom_ops import push_undo_snapshot
    from pdynamo.dyff_edits import term_signature
    symbols, coordinates, bonds = scan_inputs ( vismol_object )
    push_undo_snapshot ( vismol_object )
    edits = getattr ( vismol_object, "dyff_parameter_edits", None )
    if edits is None: edits = vismol_object.dyff_parameter_edits = { }
    if hessian_result is not None:
        edits["bond_terms"] = hessian_result["edits"]["bond_terms"]
        edits["angle_terms"] = hessian_result["edits"]["angle_terms"]
    table = edits.setdefault ( "dihedral_bonds", { } )
    n = 0
    for r in results:
        if "pairs" in r:
            table["{}-{}".format ( *r["bond"] )] = [ list ( p ) for p in r["pairs"] ]; n += 1
    edits["term_signature"] = term_signature ( len ( vismol_object.atoms ), bonds )     # same as store_hessian_edits
    return n


#=============================================================================
#  Refinement of a SELECTION (capped subsystem, pdynamo/subsystem.py)
#=============================================================================
SELECTION_SHELL = 3          # bonds of chemical environment around a selection (see fraglib/run_selection_test.py)


def selection_subsystem ( vismol_object, selection_ids, shell = SELECTION_SHELL ):
    """ Capped subsystem around the selected atoms (atom ids). Indices in the
        returned dict: subsystem order; "parent_index" maps to system order. """
    from pdynamo import subsystem as SS
    ids = sorted ( vismol_object.atoms )
    index = { aid: k for k, aid in enumerate ( ids ) }
    orders = getattr ( vismol_object, "manual_bond_orders", None ) or { }
    aromatic = getattr ( vismol_object, "manual_aromatic_bonds", None ) or set ( )
    bonds = [ ( index[a], index[b], int ( orders.get ( ( a, b ), orders.get ( ( b, a ), 1 ) ) ), ( a, b ) in aromatic or ( b, a ) in aromatic )
              for ( a, b ) in ( getattr ( vismol_object, "manual_bonds", None ) or set ( ) ) if a in index and b in index ]
    z = [ _atomic_number ( vismol_object.atoms[aid].symbol ) for aid in ids ]
    xyz = np.array ( [ vismol_object.frames[0, aid] for aid in ids ], dtype = float )
    formal = [ int ( getattr ( vismol_object.atoms[aid], "formal_charge", 0 ) or 0 ) for aid in ids ]
    selection = [ index[a] for a in selection_ids if a in index ]
    sub = SS.build_subsystem ( z, xyz, bonds, selection, shell = shell, formal_charges = formal )
    sub["selection_parent"] = selection
    sub["formal"] = [ formal[p] if p >= 0 else 0 for p in sub["parent_index"] ]
    return sub


def subsystem_factory ( vismol_object, sub ):
    """ make_system ( edits ) on the capped subsystem: DYFF with the object's manual
        types and partial charges (caps: 0), all-pairs nonbonded. """
    import pBabel                                    # noqa: F401 -- pDynamo import order
    from pMolecule.NBModel import NBModelFull
    from pdynamo import subsystem as SS
    from pdynamo.dyff_edits import dyff_model
    from pdynamo.pDynamo2EasyHybrid.session import _apply_dyff_guanidinium_correction, _apply_manual_atom_type_overrides
    ids = sorted ( vismol_object.atoms )
    local = { p: k for k, p in enumerate ( sub["parent_index"] ) if p >= 0 }
    overrides = { local[ids.index ( aid )]: label for aid, label in
                  ( getattr ( vismol_object, "manual_atom_type_overrides", None ) or { } ).items ( ) if aid in vismol_object.atoms
                  and ids.index ( aid ) in local }
    charges = [ getattr ( vismol_object.atoms[aid], "mm_charge", None ) for aid in ids ]
    has_charges = all ( q is not None for q in charges )
    def make ( edits ):
        import contextlib, io
        with contextlib.redirect_stdout ( io.StringIO ( ) ):
            system = SS.pdynamo_system ( sub["atomic_numbers"], sub["bonds"], sub["coordinates"], sub["formal"] )
            added = [ ]
            try:
                _apply_dyff_guanidinium_correction ( system )
                added = _apply_manual_atom_type_overrides ( system, overrides ) if overrides else [ ]
                system.DefineMMModel ( dyff_model ( "dyff-1.0", edits ), log = None )
            finally:
                if added:
                    from pBabel import MOL2AtomTypesToDYFF
                    for key in added: MOL2AtomTypesToDYFF.TriposMOL2AtomTypes.pop ( key, None )
            if has_charges:
                for k, p in enumerate ( sub["parent_index"] ):
                    system.mmState.charges[k] = float ( charges[p] ) if p >= 0 else 0.0
                soften_polar_hydrogens ( system )
            system.DefineNBModel ( NBModelFull.WithDefaults ( ) )
        return system
    return make


def refine_selection ( vismol_object, selection_ids, multiplicity = 1, solvent = None, shell = SELECTION_SHELL, do_hessian = True,
                       step = 15.0, workers = 4, progress = None, cancel = None, use_library = True ):
    """ xTB refinement of the selected atoms only, computed on a capped subsystem
        (selection + `shell` bonds): Hessian bonds/angles and the rotatable bonds
        whose two central atoms are selected; only terms with ALL atoms in the
        selection are returned, in system indices. Worker-thread safe. """
    import copy
    from pScientific import PeriodicTable
    from pdynamo import seminario as SM
    from pdynamo import subsystem as SS
    say = progress or ( lambda text: None )
    sub = selection_subsystem ( vismol_object, selection_ids, shell = shell )
    symbols = [ PeriodicTable.Symbol ( z ) for z in sub["atomic_numbers"] ]
    pairs = [ ( a, b ) for a, b, o, ar in sub["bonds"] ]
    selection = set ( sub["selection_parent"] )
    say ( "subsystem: {} atoms ({} capping H), charge {:+d}".format ( len ( symbols ), sub["n_caps"], sub["charge"] ) )
    # . existing edits of the object, carried into the subsystem (type keys apply by themselves)
    edits = getattr ( vismol_object, "dyff_parameter_edits", None ) or { }
    base = { k: copy.deepcopy ( v ) for k, v in edits.items ( ) if k in ( "bond", "angle", "dihedral", "outofplane" ) and v }
    base.update ( SS.map_term_edits_to_subsystem ( edits, sub["parent_index"] ) )
    hessian_parent, skipped = None, [ ]
    if do_hessian:
        say ( "xTB Hessian of the subsystem..." )
        xyz, hessian, freqs = SM.xtb_hessian ( symbols, sub["coordinates"], charge = sub["charge"], multiplicity = multiplicity,
                                               threads = 2, solvent = solvent )
        calibration = hessian_calibration ( vismol_object )
        if calibration is not None:
            from pdynamo.bond_library import HessianCalibration
            calibration = HessianCalibration.mapped ( calibration, sub["parent_index"] )
        sem = SM.dyff_edits_from_hessian ( hessian, xyz, pairs, skipped, calibration = calibration )
        base["bond_terms"] = dict ( base.get ( "bond_terms", { } ), **sem["bond_terms"] )
        base["angle_terms"] = dict ( base.get ( "angle_terms", { } ), **sem["angle_terms"] )
        hessian_parent = SS.map_term_edits_to_parent ( sem, sub["parent_index"], selection )
    # . rotatable bonds of the SYSTEM with both central atoms selected; the
    #   library is searched with the codes of the whole molecule (caps would
    #   change the environment near the cut)
    local = { p: k for k, p in enumerate ( sub["parent_index"] ) if p >= 0 }
    selected_bonds = [ ( bond, quad ) for bond, quad in rotatable_bonds ( vismol_object ) if bond[0] in selection and bond[1] in selection ]
    library = [ ]
    if use_library: library, selected_bonds = library_torsion_results ( vismol_object, selected_bonds, say )
    targets = [ ]
    for bond, quad in selected_bonds:
        if all ( a in local for a in quad ):
            q = tuple ( local[a] for a in quad )
            targets.append ( ( ( min ( q[1], q[2] ), max ( q[1], q[2] ) ), q ) )
    parent = sub["parent_index"]
    results = library + _scan_and_fit ( symbols, sub["coordinates"], pairs, targets, subsystem_factory ( vismol_object, sub ), base,
                                        sub["charge"], multiplicity, solvent, step, 1, workers, say, cancel,
                                        report_bond = lambda b: ( min ( parent[b[0]], parent[b[1]] ), max ( parent[b[0]], parent[b[1]] ) ) )
    return { "subsystem": { "n_atoms": len ( symbols ), "n_caps": sub["n_caps"], "charge": sub["charge"],
                            "n_selected": len ( selection ), "whole": sub["whole"] },
             "hessian": hessian_parent, "skipped": skipped, "results": results }


def store_selection_refinement ( vismol_object, refinement ):
    """ One Undo step; MERGES the selection's terms into the existing edits
        (other regions keep theirs). """
    from gui.windows.builder.atom_ops import push_undo_snapshot
    from pdynamo.dyff_edits import term_signature
    symbols, coordinates, bonds = scan_inputs ( vismol_object )
    push_undo_snapshot ( vismol_object )
    edits = getattr ( vismol_object, "dyff_parameter_edits", None )
    if edits is None: edits = vismol_object.dyff_parameter_edits = { }
    if refinement.get ( "hessian" ):
        edits.setdefault ( "bond_terms", { } ).update ( refinement["hessian"]["bond_terms"] )
        edits.setdefault ( "angle_terms", { } ).update ( refinement["hessian"]["angle_terms"] )
    n = 0
    for r in refinement["results"]:
        if "pairs" in r:
            edits.setdefault ( "dihedral_bonds", { } )["{}-{}".format ( *r["bond"] )] = [ list ( p ) for p in r["pairs"] ]; n += 1
    edits["term_signature"] = term_signature ( len ( vismol_object.atoms ), bonds )
    h = refinement.get ( "hessian" ) or { }
    return { "torsions": n, "bonds": len ( h.get ( "bond_terms", { } ) ), "angles": len ( h.get ( "angle_terms", { } ) ) }

