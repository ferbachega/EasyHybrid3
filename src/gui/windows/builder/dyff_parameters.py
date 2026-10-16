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
                edited = keys[p] in ( edits.get ( term ) or { } )
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
        path = os.path.join ( os.environ["PDYNAMO3_PARAMETERS"], "forceFields", "dyff", "dyff-1.0", "lennardJonesParameters.yaml" )
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
