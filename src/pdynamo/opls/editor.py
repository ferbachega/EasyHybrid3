#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: OPLS parameter tables, quality levels and manual editing
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      Backend of the "OPLS Parameters" window (no GTK here):
#
#        build_tables ( system )         -> rows of every tab (charges, bonds,
#                                           angles, dihedrals, impropers, LJ)
#                                           with a quality level per row;
#        apply_edits ( system, edits )    -> writes the edits and updates the
#                                           system.
#
#      Quality level of a parameter (key), from sources.json of the system's
#      parameter set (written by parameters.materialize_parameter_set):
#        green   from the OPLS database (pDynamo exact or by class, OPLS-AA/L,
#                OPLS-AA 2008);
#        yellow  generic: aromatic ring torsion, generic or zero improper,
#                zero torsion that is NOT a rotatable bond; for charges: out
#                of the OPLS-AA range of the class;
#        red     estimated bond/angle (geometry), or a rotatable bond whose
#                torsion has no parameters;
#        blue    edited by the user;
#        gray    unknown origin (parameter set without sources.json).
#
#      Editing:
#        - charges are changed IN PLACE in system.mmState.charges (QC, NB and
#          restraint models untouched). Ligand atoms (MM sequence components)
#          are also written to components/<res>.yaml -- every copy of the
#          ligand gets the new charge; other atoms are kept in
#          charge_overrides.json and re-applied after every rebuild. The
#          total charge must stay integral (pDynamo tolerance 1e-4);
#        - bonds, angles, dihedrals, impropers and Lennard-Jones are written to
#          the parameter files and the MM model is rebuilt. pDynamo drops the
#          NB, QC and restraint models when the MM model is redefined: NB and
#          restraints are put back, a system with a QC model is refused
#          (remove the QC region, edit, define it again).
#

import json
import math
import os
from collections import defaultdict, OrderedDict

import numpy as np
import yaml

from pdynamo.opls import parameters as P

LEVEL_COLORS = { "green": "#4caf50", "yellow": "#f2c94c", "red": "#e5534b", "blue": "#4a90d9", "gray": "#9e9e9e" }
LEVEL_TEXT   = { "green": "database", "yellow": "generic", "red": "estimated", "blue": "edited", "gray": "unknown" }
CHARGE_OVERRIDES_FILE = "charge_overrides.json"
_DATABASE_SOURCES = ( "base", "oplsaal", "oplsaa08" )


def parameter_folder ( system ):
    """ The system's own OPLS parameter folder, or None. """
    folder = ( getattr ( system, "e_input_files", None ) or { } ).get ( "opls_parameters" )
    if folder and os.path.isdir ( folder ): return folder
    model = getattr ( system, "mmModel", None )
    path = getattr ( model, "parameterSetPath", None )
    if path and os.path.isdir ( path ) and getattr ( model, "forceField", "" ).upper ( ) == "OPLS": return path
    return None

def system_types ( system ):
    state = system.mmState
    return [ state.atomTypes[i] for i in state.atomTypeIndices ]

def ligand_names ( folder ):
    components = os.path.join ( folder, "components" )
    if not os.path.isdir ( components ): return set ( )
    return { os.path.splitext ( f )[0].upper ( ) for f in os.listdir ( components ) if f.endswith ( ".yaml" ) }


def _has_torsion ( value ):
    """ A dihedral value with at least one nonzero force constant (a zero
        torsion is written as one k = 0 row). """
    return bool ( value ) and any ( float ( k ) != 0.0 for ( k, n, phase ) in value )

def _level ( term, source, value, rotatable ):
    if source is None:                         return "gray"
    if source == "edited":                     return "blue"
    if source in _DATABASE_SOURCES or source.startswith ( "pdynamo:" ): return "green"
    if source == "zero":                       return "yellow"
    if source == "estimated":
        if term in ( "bond", "angle" ):        return "red"
        if term == "dihedral":
            if _has_torsion ( value ):         return "yellow"     # generic aromatic ring torsion
            return "red" if rotatable else "yellow"
        return "yellow"                                              # generic improper
    return "gray"

def _source_text ( source ):
    if source is None: return "unknown"
    if source == "base": return "pDynamo (exact)"
    if source.startswith ( "pdynamo:" ): return "pDynamo (class)"
    return { "oplsaal": "OPLS-AA/L", "oplsaa08": "OPLS-AA 2008", "estimated": "estimated",
             "zero": "zero", "edited": "edited" }.get ( source, source )


#=============================================================================
#  Tables
#=============================================================================
def build_tables ( system ):
    """ { "folder", "ligands", "charges": [ rows ], "bond": [ rows ], "angle", "dihedral",
          "outofplane", "lj": [ rows ], "total_charge" } -- every row a dict with
        "level", "quality" (text), "ligand" (bool) and the editable values. """
    from pdynamo.opls import quality
    folder = parameter_folder ( system )
    if folder is None:
        raise ValueError ( "The system has no OPLS parameter set prepared by EasyHybrid (Extras > OPLS > Prepare OPLS System)." )
    pset    = P.PDynamoSet ( folder )
    sources = P.read_sources ( folder )
    types   = system_types ( system )
    charges = system.mmState.charges
    ligands = ligand_names ( folder )
    atoms   = list ( system.atoms )
    residue = [ a.parent.genericLabel if a.parent is not None else "" for a in atoms ]
    is_lig  = [ r.upper ( ) in ligands for r in residue ]
    c = system.connectivity
    adjacency = defaultdict ( set )
    bond_type = { }
    for b in c.bonds:
        i, j = b.node1.index, b.node2.index
        adjacency[i].add ( j ); adjacency[j].add ( i )
        bond_type[( min ( i, j ), max ( i, j ) )] = b.type.name if hasattr ( b.type, "name" ) else str ( b.type )

    instances = { t: defaultdict ( list ) for t in P.TERMS }
    for ( i, j ) in c.bondIndices:            instances["bond"][P.key_bond ( types[i], types[j] )].append ( ( i, j ) )
    for ( i, j, k ) in c.angleIndices:        instances["angle"][P.key_angle ( types[i], types[j], types[k] )].append ( ( i, j, k ) )
    for ( i, j, k, l ) in c.dihedralIndices:  instances["dihedral"][P.key_dihedral ( types[i], types[j], types[k], types[l] )].append ( ( i, j, k, l ) )
    for ( i, j, k, l ) in c._IndicesNeighbor3 ( ):
        instances["outofplane"][P.key_outofplane ( types[i], types[j], types[k], types[l] )].append ( ( i, j, k, l ) )

    def rotatable ( quad ):
        i, j, k, l = quad
        if bond_type.get ( ( min ( j, k ), max ( j, k ) ) ) != "Single": return False
        if len ( adjacency[j] ) < 2 or len ( adjacency[k] ) < 2: return False
        return not quality._in_ring ( adjacency, j, k, set ( range ( len ( atoms ) ) ) )

    def example ( tuple_ ):
        a0 = atoms[tuple_[0]]
        prefix = a0.parent.path if a0.parent is not None else ""
        return "{} {}".format ( prefix, "-".join ( atoms[a].label for a in tuple_ ) )

    tables = { "folder": folder, "ligands": sorted ( ligands ) }
    for term in ( "bond", "angle", "dihedral", "outofplane" ):
        rows = [ ]
        for key, tuples in sorted ( instances[term].items ( ) ):
            value  = pset.terms[term].get ( key )
            source = sources[term].get ( key )
            rot    = term == "dihedral" and any ( rotatable ( t ) for t in tuples )
            level  = _level ( term, source, value, rot )
            label  = LEVEL_TEXT[level]
            if term == "dihedral" and source == "estimated":
                label = "generic ring torsion" if _has_torsion ( value ) else ( "no torsion, rotatable bond" if rot else "zero torsion" )
            elif term == "outofplane" and source == "estimated": label = "generic improper"
            elif term == "outofplane" and source == "zero":      label = "zero improper"
            elif source == "estimated":                           label = "estimated (geometry)"
            rows.append ( { "key": key, "value": value, "count": len ( tuples ), "example": example ( tuples[0] ),
                            "instances": [ tuple ( int ( a ) for a in t ) for t in tuples ],
                            "source": _source_text ( source ), "level": level, "quality": label,
                            "ligand": any ( is_lig[a] for t in tuples for a in t ) } )
        tables[term] = rows
    # . Lennard-Jones: the types present.
    counts = defaultdict ( int )
    for t in types: counts[t] += 1
    lj_rows = [ ]
    for label in sorted ( counts ):
        e, s = pset.lj.get ( label, ( None, None ) )
        lig = any ( is_lig[i] for i, t in enumerate ( types ) if t == label )
        edited = sources.get ( "lj", { } ).get ( ( label, ) ) == "edited" if "lj" in sources else False
        level = "blue" if edited else ( "gray" if e is None else "green" )
        lj_rows.append ( { "key": ( label, ), "value": ( e, s ), "count": counts[label], "example": label,
                           "instances": [ ( i, ) for i, t in enumerate ( types ) if t == label ],
                           "source": "OPLS-AA 2008 class" if lig else "pDynamo", "level": level,
                           "quality": "edited" if edited else ( "class LJ" if lig else "database" ), "ligand": lig } )
    tables["lj"] = lj_rows
    # . Charges.
    ranges = quality.class_charge_ranges ( )
    overrides = _read_overrides ( folder )
    charge_rows = [ ]
    for i, a in enumerate ( atoms ):
        cls = P.atom_class ( types[i] )
        q = float ( charges[i] )
        low, high, _ = ranges.get ( cls, ( None, None, None ) )
        outlier = low is not None and ( q < low - quality.CHARGE_MARGIN or q > high + quality.CHARGE_MARGIN )
        edited = str ( i ) in overrides or ( is_lig[i] and _ligand_charge_edited ( folder, residue[i], a.label ) )
        level = "blue" if edited else ( "yellow" if outlier else "green" )
        charge_rows.append ( { "index": i, "instances": [ ( i, ) ], "residue": a.parent.path if a.parent is not None else "", "name": a.label,
                               "type": types[i], "value": q, "source": "ligand (CM5)" if is_lig[i] else "OPLS pattern",
                               "range": "{:+.3f} .. {:+.3f}".format ( low, high ) if low is not None else "",
                               "level": level, "quality": "edited" if edited else ( "outside class range" if outlier else "in range" ),
                               "ligand": is_lig[i] } )
    tables["charges"] = charge_rows
    tables["total_charge"] = float ( sum ( charges ) )
    return tables


#=============================================================================
#  Editing
#=============================================================================
def _read_overrides ( folder ):
    path = os.path.join ( folder, CHARGE_OVERRIDES_FILE )
    if not os.path.exists ( path ): return { }
    with open ( path ) as handle: return json.load ( handle )

def _ligand_charge_edited ( folder, residue_name, atom_name ):
    path = os.path.join ( folder, "edited_ligand_charges.json" )
    if not os.path.exists ( path ): return False
    with open ( path ) as handle: data = json.load ( handle )
    return atom_name in data.get ( residue_name.upper ( ), [ ] )

def _component_path ( folder, residue_name ):
    return os.path.join ( folder, "components", residue_name.lower ( ) + ".yaml" )

def check_total_charge ( charges, tolerance = 1.0e-4 ):
    total = float ( sum ( charges ) )
    if math.fabs ( total - round ( total ) ) > tolerance:
        raise ValueError ( "The total charge would be {:.4f}: it must be an integer (tolerance {:g}). "
                           "Adjust other charges first.".format ( total, tolerance ) )
    return total

def apply_edits ( system, edits, rebuild_callback = None ):
    """ edits: { "charges": { atom index: q }, "bond": { key: ( b0, k ) },
        "angle": { key: ( t0, k ) }, "dihedral": { key: [ ( k, n, phase ) ] },
        "outofplane": { key: K }, "lj": { ( label, ): ( epsilon, sigma ) } }.
        Returns a summary dict. """
    folder = parameter_folder ( system )
    if folder is None: raise ValueError ( "The system has no EasyHybrid OPLS parameter set." )
    bonded = { t: v for t, v in edits.items ( ) if t != "charges" and v }
    if bonded and getattr ( system, "qcModel", None ) is not None:
        raise ValueError ( "This system has a QC model: pDynamo removes it when the MM model is rebuilt. "
                           "Remove the QC region, edit the parameters, then define the QC region again "
                           "(charges can be edited with the QC model in place)." )
    types = system_types ( system )
    summary = { "charges": 0, "terms": 0, "rebuilt": False }
    # . Validate charges first.
    new_charges = None
    charge_edits = { int ( i ): float ( q ) for i, q in ( edits.get ( "charges" ) or { } ).items ( ) }
    ligands = ligand_names ( folder )
    atoms = list ( system.atoms )
    if charge_edits:
        new_charges = [ float ( q ) for q in system.mmState.charges ]
        # . a ligand atom edit applies to every copy of the ligand (same residue name and atom name)
        expanded = dict ( charge_edits )
        for i, q in charge_edits.items ( ):
            res = atoms[i].parent.genericLabel if atoms[i].parent is not None else ""
            if res.upper ( ) in ligands:
                for k, a in enumerate ( atoms ):
                    if a.parent is not None and a.parent.genericLabel == res and a.label == atoms[i].label: expanded[k] = q
        for i, q in expanded.items ( ): new_charges[i] = q
        check_total_charge ( new_charges )
        charge_edits = expanded
    # . Bonded / LJ: write the parameter files.
    if bonded:
        pset = P.PDynamoSet ( folder )
        sources = P.read_sources ( folder )
        files = { "bond": "harmonicBondParameters.yaml", "angle": "harmonicAngleParameters.yaml",
                  "dihedral": "fourierDihedralParameters.yaml", "outofplane": "fourierOutOfPlaneParameters.yaml",
                  "lj": "lennardJonesParameters.yaml" }
        for term, changes in bonded.items ( ):
            path = os.path.join ( folder, files[term] )
            with open ( path ) as handle: mapping = yaml.safe_load ( handle )
            rows = mapping["Parameter Values"]
            if term == "dihedral":
                keep = [ r for r in rows if P.key_dihedral ( *[ str ( x ) for x in r[:4] ] ) not in changes ]
                for key, terms in changes.items ( ):
                    for ( k, n, phase ) in ( terms or [ ( 0.0, 0, 0.0 ) ] ):
                        keep.append ( [ key[0], key[1], key[2], key[3], float ( k ), int ( n ), float ( phase ) ] )
                rows = keep
            else:
                for r in rows:
                    if term == "bond":       key = P.key_bond ( str ( r[0] ), str ( r[1] ) )
                    elif term == "angle":    key = P.key_angle ( str ( r[0] ), str ( r[1] ), str ( r[2] ) )
                    elif term == "outofplane": key = P.key_outofplane ( *[ str ( x ) for x in r[:4] ] )
                    else:                    key = ( str ( r[0] ), )
                    if key in changes:
                        new = changes[key]
                        if term == "bond":         r[2], r[3] = float ( new[0] ), float ( new[1] )
                        elif term == "angle":      r[3], r[4] = float ( new[0] ), float ( new[1] )
                        elif term == "outofplane": r[4] = float ( new )
                        else:                      r[1], r[2] = float ( new[0] ), float ( new[1] )
            mapping["Parameter Values"] = rows
            P._dump ( path, mapping )
            sources.setdefault ( term, { } )
            for key in changes: sources[term][tuple ( key )] = "edited"
            summary["terms"] += len ( changes )
        P.write_sources ( folder, sources )
    # . Persist the charges.
    if charge_edits:
        overrides = _read_overrides ( folder )
        edited_lig = { }
        path_edited = os.path.join ( folder, "edited_ligand_charges.json" )
        if os.path.exists ( path_edited ):
            with open ( path_edited ) as handle: edited_lig = json.load ( handle )
        components = defaultdict ( dict )
        for i, q in charge_edits.items ( ):
            res = atoms[i].parent.genericLabel if atoms[i].parent is not None else ""
            if res.upper ( ) in ligands:
                components[res][atoms[i].label] = q
            else:
                overrides[str ( i )] = q
        for res, by_name in components.items ( ):
            path = _component_path ( folder, res )
            with open ( path ) as handle: mapping = yaml.safe_load ( handle )
            for row in mapping["Atoms"]:
                if row[0] in by_name: row[2] = round ( by_name[row[0]], 6 )
            with open ( path, "w" ) as handle:
                handle.write ( "# . !MMSequenceComponent\n" )
                yaml.safe_dump ( mapping, handle, default_flow_style = None, sort_keys = False, width = 200 )
            edited_lig.setdefault ( res.upper ( ), [ ] )
            edited_lig[res.upper ( )] = sorted ( set ( edited_lig[res.upper ( )] ) | set ( by_name ) )
        with open ( os.path.join ( folder, CHARGE_OVERRIDES_FILE ), "w" ) as handle: json.dump ( overrides, handle, indent = 0 )
        with open ( path_edited, "w" ) as handle: json.dump ( edited_lig, handle, indent = 0 )
        summary["charges"] = len ( charge_edits )
    # . Update the system.
    if bonded:
        rebuild_mm_model ( system, folder )
        summary["rebuilt"] = True
    elif charge_edits:
        for i, q in charge_edits.items ( ): system.mmState.charges[i] = q
        system.scratch.Clear ( ) if hasattr ( system.scratch, "Clear" ) else None
    summary["total_charge"] = float ( sum ( system.mmState.charges ) )
    return summary

def rebuild_mm_model ( system, folder ):
    """ Redefines the OPLS MM model from `folder`, putting back the NB and
        restraint models and the charge overrides of non-ligand atoms. """
    from pMolecule.MMModel import MMModelOPLS
    nb_model = getattr ( system, "nbModel", None )
    restraints = getattr ( system, "restraintModel", None )
    system.DefineMMModel ( MMModelOPLS.WithParameterSet ( os.path.abspath ( folder ) ), log = None )
    overrides = _read_overrides ( folder )
    for i, q in overrides.items ( ): system.mmState.charges[int ( i )] = float ( q )
    if nb_model is not None: system.DefineNBModel ( nb_model )
    if restraints is not None: system.DefineRestraintModel ( restraints )


#=============================================================================
#  Geometry of one term (shown next to the parameter)
#=============================================================================
def measure ( term, points ):
    """ The internal coordinate of a term from its atoms' coordinates
        (points in the term's atom order): bond length (A), angle (deg),
        dihedral (deg); for an out-of-plane ( central, j, k, l ): the
        height of the central atom above the j-k-l plane (A). None for
        charges and Lennard-Jones. """
    p = [ np.asarray ( x, dtype = float ) for x in points ]
    if term == "bond":
        return float ( np.linalg.norm ( p[0] - p[1] ) )
    if term == "angle":
        u, v = p[0] - p[1], p[2] - p[1]
        return float ( np.degrees ( np.arccos ( np.clip ( np.dot ( u, v ) / ( np.linalg.norm ( u ) * np.linalg.norm ( v ) ), -1.0, 1.0 ) ) ) )
    if term == "dihedral":
        b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
        b1n = b1 / np.linalg.norm ( b1 )
        v = b0 - np.dot ( b0, b1n ) * b1n
        w = b2 - np.dot ( b2, b1n ) * b1n
        return float ( np.degrees ( np.arctan2 ( np.dot ( np.cross ( b1n, v ), w ), np.dot ( v, w ) ) ) )
    if term == "outofplane":
        normal = np.cross ( p[2] - p[1], p[3] - p[1] )
        return float ( abs ( np.dot ( p[0] - p[1], normal ) ) / np.linalg.norm ( normal ) )
    return None
