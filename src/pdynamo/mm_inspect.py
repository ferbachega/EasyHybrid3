#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: read-only tables of the MM parameters of any pDynamo system
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description (2026-10-21):
#      The Force Field Parameters window inspects any MM model, whatever its
#      origin (OPLS, AMBER, CHARMM, DYFF...). The parameters are read from
#      the term containers of system.mmState (their __getstate__ gives the
#      terms and parameters of every container type in the same layout), so
#      no parameter file is needed. One row per parameter of a container,
#      with every term that uses it.
#
#      Units are pDynamo's: kJ/mol and Angstroms; angles shown in degrees.
#      pDynamo harmonic terms are E = k (x - x0)^2 (no 1/2).
#

import math
from collections import defaultdict

FORCE_FIELD_NAMES = { "MMModelOPLS": "OPLS-AA", "MMModelAMBER": "AMBER", "MMModelCHARMM": "CHARMM",
                      "MMModelDYFF": "DYFF", "MMModelDYFFEdited": "DYFF (EasyHybrid edits)" }

# . kind of geometry measured for the atoms of one term of a container
_MEASURE = { "HarmonicBondContainer": "bond", "HarmonicAngleContainer": "angle", "HarmonicImproperContainer": "dihedral",
             "FourierDihedralContainer": "dihedral", "CosineDihedralContainer": "dihedral", "CosineAngleContainer": "angle",
             "CosineOutOfPlaneContainer": "outofplane" }


def force_field_name ( system ):
    model = getattr ( system, "mmModel", None )
    if model is None: return None
    name = type ( model ).__name__
    return FORCE_FIELD_NAMES.get ( name, name.replace ( "MMModel", "" ) or name )


def _deg ( value ): return math.degrees ( float ( value ) )


def _term_atoms_and_parameter ( container_class, term ):
    """ ( atom indices, parameter index, active ) of one entry of __getstate__ ( )["terms"].
        Cosine containers store [ indices, parameter, active ]; the other ones
        [ atom1, ..., parameter, active ] (HarmonicAngle's termFields lists an atom4
        that its terms do not have, so the fields are not used). """
    if isinstance ( term[0], ( list, tuple ) ):
        return tuple ( int ( i ) for i in term[0] ), int ( term[1] ), bool ( term[2] )
    return tuple ( int ( i ) for i in term[:-2] ), int ( term[-2] ), bool ( term[-1] )


def _parameter_columns ( container_class, label ):
    """ [ ( title, function ( parameter ) -> text ) ] for one container type. """
    f = lambda v, d: ( "{:.%df}" % d ).format ( float ( v ) )
    if container_class == "HarmonicBondContainer":
        return [ ( "r0 (Å)", lambda p: f ( p[0], 4 ) ), ( "k (kJ/mol/Å²)", lambda p: f ( p[1], 2 ) ) ]
    if container_class in ( "HarmonicAngleContainer", "HarmonicImproperContainer" ):
        return [ ( "θ0 (°)", lambda p: f ( _deg ( p[0] ), 2 ) ), ( "k (kJ/mol/rad²)", lambda p: f ( p[1], 2 ) ) ]
    if container_class == "FourierDihedralContainer":
        return [ ( "V (kJ/mol)", lambda p: f ( p[0], 4 ) ), ( "n", lambda p: str ( int ( p[1] ) ) ),
                 ( "phase (°)", lambda p: f ( _deg ( p[2] ), 1 ) ) ]
    if container_class.startswith ( "Cosine" ):
        return [ ( "Σ c·cos(nθ): c (kJ/mol) n", lambda p: "; ".join ( "{:.4f} {:d}".format ( float ( c ), int ( n ) ) for c, n in p ) or "none" ) ]
    if container_class == "CMAPDihedralContainer":
        return [ ( "Grid", lambda p: "{} x {} map".format ( len ( p[0] ), len ( p[1] ) ) ) ]
    return [ ( "Parameters", lambda p: str ( p ) ) ]


def build_tables ( system ):
    """ { "force_field", "kinds": [ ( key, title ) ], "columns": { key: [ ( title, func ( row ) ) ] },
          "measure": { key: bond/angle/dihedral/outofplane/None }, key: [ rows ], "total_charge" }.
        Every row: "key" (type tuple), "value", "count", "instances", "example", "source". """
    state = getattr ( system, "mmState", None )
    if getattr ( system, "mmModel", None ) is None or state is None:
        raise ValueError ( "The system has no MM model." )
    atoms = list ( system.atoms )
    n = len ( atoms )
    if state.atomTypes is not None and state.atomTypeIndices is not None:
        types = [ str ( state.atomTypes[state.atomTypeIndices[i]] ) for i in range ( n ) ]
    else:
        types = [ a.label for a in atoms ]
    residue = [ a.parent.path if a.parent is not None else "" for a in atoms ]

    def example ( indices ):
        a0 = atoms[indices[0]]
        prefix = a0.parent.path if a0.parent is not None else ""
        return "{} {}".format ( prefix, "-".join ( atoms[i].label for i in indices ) ).strip ( )

    tables = { "force_field": force_field_name ( system ), "kinds": [ ( "charges", "Charges" ) ], "columns": { }, "measure": { } }

    # . charges
    charges = [ float ( q ) for q in state.charges ] if state.charges is not None else [ 0.0 ] * n
    tables["charges"] = [ { "index": i, "key": ( types[i], ), "value": q, "instances": [ ( i, ) ], "count": 1,
                            "residue": residue[i], "name": atoms[i].label, "type": types[i], "example": example ( ( i, ) ),
                            "source": "" } for i, q in enumerate ( charges ) ]
    tables["columns"]["charges"] = [ ( "Atom", lambda r: "{} {}".format ( r["residue"], r["name"] ).strip ( ) ),
                                     ( "Type", lambda r: r["type"] ), ( "Charge (e)", lambda r: "{:+.4f}".format ( r["value"] ) ) ]
    tables["measure"]["charges"] = None
    tables["total_charge"] = sum ( charges )

    # . bonded terms, one kind per container
    titles = defaultdict ( int )
    for k, container in enumerate ( state.mmTerms or [ ] ):
        cls = type ( container ).__name__
        data = container.__getstate__ ( )
        label = data.get ( "label" ) or cls.replace ( "Container", "" )
        titles[label] += 1
        title = label if titles[label] == 1 else "{} ({})".format ( label, titles[label] )
        keys = data.get ( "parameterKeys" ) or [ ]
        parameters = data["parameters"]
        groups = defaultdict ( list )
        for term in data["terms"]:
            indices, p, active = _term_atoms_and_parameter ( cls, term )
            groups[p].append ( indices )
        rows = [ ]
        for p, instances in groups.items ( ):
            tuples = defaultdict ( int )
            for t in instances: tuples[tuple ( types[i] for i in t )] += 1
            ordered = sorted ( tuples, key = lambda t: - tuples[t] )
            rows.append ( { "key": ordered[0], "other_keys": ordered[1:], "value": parameters[p] if p < len ( parameters ) else None,
                            "count": len ( instances ), "instances": instances, "example": example ( instances[0] ),
                            "source": "parameter {}".format ( keys[p] ) if p < len ( keys ) and keys[p] else "parameter #{}".format ( p ) } )
        rows.sort ( key = lambda r: r["key"] )
        name = "term{}".format ( k )
        tables[name] = rows
        tables["kinds"].append ( ( name, title ) )
        value_columns = _parameter_columns ( cls, label )
        def types_text ( r ):
            text = " - ".join ( r["key"] )
            return text + ( "  (+{} more)".format ( len ( r["other_keys"] ) ) if r["other_keys"] else "" )
        tables["columns"][name] = [ ( "Types", types_text ) ] + \
                                  [ ( t, ( lambda r, fn = fn: fn ( r["value"] ) if r["value"] is not None else "" ) ) for t, fn in value_columns ] + \
                                  [ ( "Terms", lambda r: str ( r["count"] ) ) ]
        tables["measure"][name] = _MEASURE.get ( cls )
        if cls == "HarmonicBondContainer" and "Urey" in label: tables["measure"][name] = "bond"

    # . Lennard-Jones, one row per LJ type
    lj = state.ljParameters
    if lj is not None and state.ljTypeIndices is not None:
        data = lj.__getstate__ ( )
        data14 = state.ljParameters14.__getstate__ ( ) if state.ljParameters14 is not None else None
        lj_keys = data.get ( "parameterKeys" ) or [ ]
        members = defaultdict ( list )
        for i in range ( n ): members[int ( state.ljTypeIndices[i] )].append ( i )
        form = str ( data.get ( "analyticForm", "" ) ).split ( "." )[-1]
        rows = [ ]
        for t, atoms_t in sorted ( members.items ( ) ):
            name = lj_keys[t] if t < len ( lj_keys ) else ", ".join ( sorted ( { types[i] for i in atoms_t } ) )
            value = ( data["epsilons"][t], data["sigmas"][t] ) if "epsilons" in data else None
            value14 = ( data14["epsilons"][t], data14["sigmas"][t] ) if data14 is not None and "epsilons" in data14 else None
            rows.append ( { "key": ( str ( name ), ), "value": value, "value14": value14, "count": len ( atoms_t ),
                            "instances": [ ( i, ) for i in atoms_t ], "example": example ( ( atoms_t[0], ) ),
                            "source": "{} form".format ( form ) if form else "" } )
        tables["lj"] = rows
        tables["kinds"].append ( ( "lj", "Lennard-Jones" ) )
        # . the second parameter is sigma for the OPLS form and Rmin (= 2 Rmin/2) for the AMBER form
        second = "σ (Å)" if form.upper ( ) == "OPLS" else "Rmin (Å)"
        columns = [ ( "LJ type", lambda r: r["key"][0] ),
                    ( "ε (kJ/mol)", lambda r: "{:.5f}".format ( r["value"][0] ) if r["value"] else "table" ),
                    ( second, lambda r: "{:.4f}".format ( r["value"][1] ) if r["value"] else "table" ) ]
        if data14 is not None and any ( r["value14"] != r["value"] for r in rows ):
            columns += [ ( "ε 1-4", lambda r: "{:.5f}".format ( r["value14"][0] ) if r["value14"] else "" ),
                         ( second.replace ( " (Å)", " 1-4" ), lambda r: "{:.4f}".format ( r["value14"][1] ) if r["value14"] else "" ) ]
        columns.append ( ( "Atoms", lambda r: str ( r["count"] ) ) )
        tables["columns"]["lj"] = columns
        tables["measure"]["lj"] = None
    return tables


def measure ( kind, points ):
    """ Bond length (Å), angle or dihedral (degrees), or the distance (Å) of the
        central (second) atom from the plane of the other three. """
    import numpy as np
    p = [ np.array ( x, dtype = float ) for x in points ]
    try:
        if kind == "bond" and len ( p ) >= 2:
            return float ( np.linalg.norm ( p[1] - p[0] ) )
        if kind == "angle" and len ( p ) >= 3:
            a, b = p[0] - p[1], p[2] - p[1]
            return math.degrees ( math.acos ( max ( -1.0, min ( 1.0, float ( a @ b / ( np.linalg.norm ( a ) * np.linalg.norm ( b ) ) ) ) ) ) )
        if kind == "dihedral" and len ( p ) >= 4:
            b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
            b1n = b1 / np.linalg.norm ( b1 )
            v = b0 - ( b0 @ b1n ) * b1n
            w = b2 - ( b2 @ b1n ) * b1n
            return math.degrees ( math.atan2 ( float ( np.cross ( b1n, v ) @ w ), float ( v @ w ) ) )
        if kind == "outofplane" and len ( p ) >= 4:
            # . pDynamo's out-of-plane terms list the central atom second
            normal = np.cross ( p[2] - p[0], p[3] - p[0] )
            return abs ( float ( ( p[1] - p[0] ) @ normal / np.linalg.norm ( normal ) ) )
    except ( ValueError, ZeroDivisionError, FloatingPointError ):
        return None
    return None
