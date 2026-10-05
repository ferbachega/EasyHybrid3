#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: OPLS parameter database and per-system parameter sets
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      pDynamo's OPLS parameter sets key EVERY parameter by the full atom
#      type label ("HC.METHANOL", "N.PRIMARYAMIDE", ...) and its angle and
#      out-of-plane containers have no fallback, so a set only works for
#      the label combinations written in it. This module:
#
#        1. reads a pDynamo OPLS set (default: pDynamo's opls/protein);
#        2. builds a CLASS-level database (class = label before the dot)
#           merging, in this order of precedence:
#              pDynamo opls/protein  >  oplsaal.prm (OPLS-AA/L)  >  oplsaa08.prm (OPLS-AA 2008)
#           every entry remembers its source and any different value
#           another source gives (conflict report);
#        3. materializes a parameter set FOR ONE SYSTEM: the terms the
#           system really has (bonds, angles, dihedrals, out-of-planes of
#           its connectivity, with its atom type labels) are looked up
#           exactly in the base set first, then by class in the database,
#           and written as a pDynamo parameter folder usable with
#           MMModelOPLS.WithParameterSet ( folder ).
#
#      The Tinker files live in data/tinker/ (copied from the user's
#      reference files; see data/tinker/README.md).
#

import os
import json
import shutil
from collections import OrderedDict, Counter

import yaml

from pdynamo.opls.tinker_prm import TinkerPRM, WILD

DATA_DIR    = os.path.join ( os.path.dirname ( os.path.abspath ( __file__ ) ), "data" )
TINKER_DIR  = os.path.join ( DATA_DIR, "tinker" )
TINKER_FILES = [ ( "oplsaal",  os.path.join ( TINKER_DIR, "oplsaal.prm"  ) ),
                 ( "oplsaa08", os.path.join ( TINKER_DIR, "oplsaa08.prm" ) ) ]

TERMS = ( "bond", "angle", "dihedral", "outofplane", "ureybradley" )

_FILES = { "atomTypes"   : "atomTypes.yaml",
           "patterns"    : "patterns.yaml",
           "lj"          : "lennardJonesParameters.yaml",
           "bond"        : "harmonicBondParameters.yaml",
           "angle"       : "harmonicAngleParameters.yaml",
           "dihedral"    : "fourierDihedralParameters.yaml",
           "outofplane"  : "fourierOutOfPlaneParameters.yaml",
           "ureybradley" : "ureyBradleyParameters.yaml" }


#=============================================================================
#  Keys -- identical to pDynamo's MakeKey() of each container
#=============================================================================
def key_bond ( a, b ):
    return ( a, b ) if a >= b else ( b, a )

def key_angle ( a, b, c ):
    return ( a, b, c ) if a >= c else ( c, b, a )

def key_dihedral ( a, b, c, d ):
    if   b > c: return ( a, b, c, d )
    elif b < c: return ( d, c, b, a )
    return ( a, b, c, d ) if a > d else ( d, c, b, a )

def key_outofplane ( central, p1, p2, p3 ):
    return tuple ( [ central ] + sorted ( [ p1, p2, p3 ] ) )

KEY = { "bond": key_bond, "angle": key_angle, "dihedral": key_dihedral,
        "outofplane": key_outofplane, "ureybradley": key_angle }

def atom_class ( label ):
    """ OPLS class of a pDynamo type label: the part before the dot. """
    return label.split ( "." )[0]


#=============================================================================
#  pDynamo parameter sets
#=============================================================================
class PDynamoSet:
    """ A pDynamo OPLS parameter folder read into plain dictionaries. """
    def __init__ ( self, path ):
        self.path = path
        self.raw  = { }
        for name, filename in _FILES.items ( ):
            full = os.path.join ( path, filename )
            if os.path.exists ( full ):
                with open ( full ) as handle:
                    self.raw[name] = yaml.safe_load ( handle )
        self.atom_types = OrderedDict ( )        # label -> row (list)
        for row in self.raw["atomTypes"]["Parameter Values"]:
            self.atom_types[str ( row[0] )] = row
        self.lj = OrderedDict ( ( str ( r[0] ), ( float ( r[1] ), float ( r[2] ) ) )
                                for r in self.raw["lj"]["Parameter Values"] )     # label -> ( epsilon, sigma )
        self.terms = { t: OrderedDict ( ) for t in TERMS }
        for r in self.raw.get ( "bond", { } ).get ( "Parameter Values", [ ] ):
            self.terms["bond"][key_bond ( str ( r[0] ), str ( r[1] ) )] = ( float ( r[2] ), float ( r[3] ) )
        for r in self.raw.get ( "angle", { } ).get ( "Parameter Values", [ ] ):
            self.terms["angle"][key_angle ( str ( r[0] ), str ( r[1] ), str ( r[2] ) )] = ( float ( r[3] ), float ( r[4] ) )
        for r in self.raw.get ( "dihedral", { } ).get ( "Parameter Values", [ ] ):
            k = key_dihedral ( *[ str ( x ) for x in r[:4] ] )
            self.terms["dihedral"].setdefault ( k, [ ] ).append ( ( float ( r[4] ), int ( r[5] ), float ( r[6] ) ) )
        for r in self.raw.get ( "outofplane", { } ).get ( "Parameter Values", [ ] ):
            self.terms["outofplane"][key_outofplane ( *[ str ( x ) for x in r[:4] ] )] = float ( r[4] )
        for r in self.raw.get ( "ureybradley", { } ).get ( "Parameter Values", [ ] ):
            self.terms["ureybradley"][key_angle ( str ( r[0] ), str ( r[1] ), str ( r[2] ) )] = ( float ( r[3] ), float ( r[4] ) )

def pdynamo_set_path ( name = "protein" ):
    if os.path.isdir ( name ): return os.path.abspath ( name )
    return os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "forceFields", "opls", name )


#=============================================================================
#  Class database
#=============================================================================
def _same ( v1, v2, tolerance = 1.0e-6 ):
    if isinstance ( v1, ( list, tuple ) ) and isinstance ( v2, ( list, tuple ) ):
        if len ( v1 ) != len ( v2 ): return False
        return all ( _same ( a, b, tolerance ) for a, b in zip ( v1, v2 ) )
    try:
        return abs ( float ( v1 ) - float ( v2 ) ) <= tolerance * max ( 1.0, abs ( float ( v1 ) ) )
    except ( TypeError, ValueError ):
        return v1 == v2

def _normalize_dihedral ( terms ):
    """ Sorted ( k, n, phase ) without zero force constants. """
    return sorted ( ( round ( k, 6 ), int ( n ), round ( p % 360.0, 3 ) ) for ( k, n, p ) in terms if k != 0.0 )

class ClassDatabase:
    """ entries[term][key] = { "value", "source", "others": [ ( source, value ) ] } """
    def __init__ ( self ):
        self.entries = { t: OrderedDict ( ) for t in TERMS }
        self.lj_types = OrderedDict ( )      # ( source, type ) -> { symbol, description, z, sigma, epsilon, charge }
        self.sources = [ ]

    def add ( self, term, key, value, source ):
        entry = self.entries[term].get ( key )
        if entry is None:
            self.entries[term][key] = { "value": value, "source": source, "others": [ ] }
            return
        reference = _normalize_dihedral ( entry["value"] ) if term == "dihedral" else entry["value"]
        candidate = _normalize_dihedral ( value ) if term == "dihedral" else value
        if not _same ( reference, candidate ):
            if not any ( s == source and _same ( v, value ) for ( s, v ) in entry["others"] ):
                entry["others"].append ( ( source, value ) )

    def lookup ( self, term, labels ):
        """ ( value, source, key ) for a tuple of TYPE LABELS, by class;
            exact class key first, then Tinker wild cards. """
        classes = [ atom_class ( l ) for l in labels ]
        key = KEY[term] ( *classes )
        entry = self.entries[term].get ( key )
        if entry is None and term == "dihedral":
            key = key_dihedral ( WILD, classes[1], classes[2], WILD )
            entry = self.entries[term].get ( key )
            if entry is None:
                for a, d in ( ( classes[0], WILD ), ( WILD, classes[3] ) ):
                    key = key_dihedral ( a, classes[1], classes[2], d )
                    entry = self.entries[term].get ( key )
                    if entry is not None: break
        if entry is None and term == "outofplane":
            central, peripherals = classes[0], classes[1:]
            for n_wild in ( 1, 2 ):
                from itertools import combinations
                for wild in combinations ( range ( 3 ), n_wild ):
                    trial = [ WILD if i in wild else p for i, p in enumerate ( peripherals ) ]
                    key = key_outofplane ( central, *trial )
                    entry = self.entries[term].get ( key )
                    if entry is not None: break
                if entry is not None: break
        if entry is None: return None, None, None
        return entry["value"], entry["source"], key

    def conflicts ( self ):
        out = { }
        for term in TERMS:
            out[term] = [ ( key, e["source"], e["value"], e["others"] ) for key, e in self.entries[term].items ( ) if e["others"] ]
        return out

def build_class_database ( base = "protein", tinker_files = None, include_base = True ):
    """ Merged class-level database (see module docstring). include_base=False:
        Tinker files only (used to cross-check the conversions). """
    db = ClassDatabase ( )
    # . 1. pDynamo set: only its plain class labels (no dot) are class-level data.
    if include_base:
        pset = PDynamoSet ( pdynamo_set_path ( base ) )
        origin = os.path.join ( pset.path, "SOURCE" )
        name = open ( origin ).read ( ).strip ( ) if os.path.exists ( origin ) else os.path.basename ( pset.path )
        source = "pdynamo:{}".format ( name )
        db.sources.append ( source )
        for term in TERMS:
            for key, value in pset.terms[term].items ( ):
                if all ( "." not in label for label in key ):
                    db.add ( term, key, value, source )
    # . 2. Tinker files.
    for ( name, path ) in ( tinker_files or TINKER_FILES ):
        if not os.path.exists ( path ): continue
        prm = TinkerPRM ( path )
        db.sources.append ( name )
        for term, items in prm.converted_terms ( ).items ( ):
            for labels, value in items:
                if any ( l.startswith ( "?" ) for l in labels ): continue
                db.add ( term, KEY[term] ( *labels ), value, name )
        for atom_type, data in prm.atoms.items ( ):
            lj = prm.type_lj ( atom_type )
            if lj is None: continue
            db.lj_types[( name, atom_type )] = { "symbol": data["symbol"], "description": data["description"],
                                                 "z": data["z"], "sigma": lj[0], "epsilon": lj[1],
                                                 "charge": prm.charges.get ( atom_type ) }
    return db


#=============================================================================
#  Checks
#=============================================================================
def compare_with_pdynamo ( base = "protein", tinker_files = None ):
    """ For class-level keys present both in the pDynamo set and in each
        Tinker file: { source: { term: ( n_equal, n_different, [examples] ) } }.
        A high agreement validates the unit conversions of tinker_prm. """
    pset = PDynamoSet ( pdynamo_set_path ( base ) )
    report = { }
    for ( name, path ) in ( tinker_files or TINKER_FILES ):
        if not os.path.exists ( path ): continue
        prm = TinkerPRM ( path )
        per_term = { }
        for term, items in prm.converted_terms ( ).items ( ):
            first = OrderedDict ( )
            for labels, value in items:
                first.setdefault ( KEY[term] ( *labels ), value )
            equal, different, examples = 0, 0, [ ]
            for key, value in first.items ( ):
                if key not in pset.terms[term]: continue
                mine = pset.terms[term][key]
                a = _normalize_dihedral ( mine ) if term == "dihedral" else mine
                b = _normalize_dihedral ( value ) if term == "dihedral" else value
                if _same ( a, b, 1.0e-3 ): equal += 1
                else:
                    different += 1
                    if len ( examples ) < 6: examples.append ( ( key, mine, value ) )
            per_term[term] = ( equal, different, examples )
        report[name] = per_term
    return report


#=============================================================================
#  Per-system parameter sets
#=============================================================================
def system_terms ( system, types ):
    """ { term: set of label keys } actually present in the system. """
    c = system.connectivity
    used = { t: set ( ) for t in TERMS }
    for ( i, j ) in c.bondIndices:
        used["bond"].add ( key_bond ( types[i], types[j] ) )
    for ( i, j, k ) in c.angleIndices:
        used["angle"].add ( key_angle ( types[i], types[j], types[k] ) )
    for ( i, j, k, l ) in c.dihedralIndices:
        used["dihedral"].add ( key_dihedral ( types[i], types[j], types[k], types[l] ) )
    for ( i, j, k, l ) in c._IndicesNeighbor3 ( ):
        used["outofplane"].add ( key_outofplane ( types[i], types[j], types[k], types[l] ) )
    return used

def _dump ( path, mapping ):
    with open ( path, "w" ) as handle:
        handle.write ( "---\n" )
        yaml.safe_dump ( mapping, handle, default_flow_style = None, sort_keys = False, width = 200 )

def extend_parameter_set ( out_dir, base = "protein", ligands = ( ) ):
    """ A full copy of the base set plus, for each ligand.LigandParameters,
        its atom types, LJ and MM sequence component (components/<res>.yaml,
        read by pDynamo's TypeBySequence before the patterns). Used both to
        type a system with ligands and as the base of its materialized set. """
    pset = PDynamoSet ( pdynamo_set_path ( base ) )
    if os.path.exists ( out_dir ): shutil.rmtree ( out_dir )
    shutil.copytree ( pset.path, out_dir )
    atom_types = dict ( pset.raw["atomTypes"] )
    atom_types["Parameter Values"] = list ( pset.raw["atomTypes"]["Parameter Values"] )
    lj = dict ( pset.raw["lj"] )
    lj["Parameter Values"] = [ [ label, e, s ] for label, ( e, s ) in pset.lj.items ( ) ]
    components = os.path.join ( out_dir, "components" )
    for ligand in ligands:
        atom_types["Parameter Values"] += ligand.atom_type_rows ( )
        for label, ( e, s ) in ligand.lj.items ( ):
            lj["Parameter Values"].append ( [ label, e, s ] )
        os.makedirs ( components, exist_ok = True )
        mapping = ligand.sequence_component_mapping ( )
        with open ( os.path.join ( components, ligand.residue_name.lower ( ) + ".yaml" ), "w" ) as handle:
            handle.write ( "# . !MMSequenceComponent\n" )
            yaml.safe_dump ( mapping, handle, default_flow_style = None, sort_keys = False, width = 200 )
    _dump ( os.path.join ( out_dir, _FILES["atomTypes"] ), atom_types )
    _dump ( os.path.join ( out_dir, _FILES["lj"] ), lj )
    with open ( os.path.join ( out_dir, "SOURCE" ), "w" ) as handle:
        handle.write ( os.path.basename ( pset.path ) + "\n" )
    return out_dir

def materialize_parameter_set ( system, types, out_dir, base = "protein", database = None,
                                prefer_base = True, extra_atom_types = None, extra_patterns = None,
                                extra_lj = None, fallback_terms = None ):
    """ Writes a pDynamo OPLS parameter folder with exactly the terms of
        `system` (types = one label per atom, e.g. from prep.type_atoms()).

        Lookup per term: the base set with the full labels (prefer_base),
        then the class database. Out-of-planes absent from every source
        get K = 0 (OPLS only defines impropers for some centres; pDynamo
        nevertheless wants an entry for each one). Bonds, angles and
        dihedrals that no source has are reported as missing.

        extra_atom_types: rows to append to atomTypes.yaml ([label, Z, charge, hydrogenType, description]);
        extra_patterns:   pattern mappings to append to patterns.yaml;
        extra_lj:         { label: ( epsilon, sigma ) }.

        Returns a report { "path", "sources": {term: Counter}, "missing": {term: [keys]},
                           "zero_outofplane": [keys] }. """
    database = database or build_class_database ( base )
    pset = PDynamoSet ( pdynamo_set_path ( base ) )
    used = system_terms ( system, types )
    os.makedirs ( out_dir, exist_ok = True )

    report = { "path": out_dir, "sources": { t: Counter ( ) for t in TERMS },
               "missing": { t: [ ] for t in TERMS }, "zero_outofplane": [ ],
               "key_source": { t: { } for t in TERMS }, "key_value": { t: { } for t in TERMS } }
    values = { t: OrderedDict ( ) for t in TERMS }
    for term in ( "bond", "angle", "dihedral", "outofplane" ):
        for key in sorted ( used[term] ):
            value, source = None, None
            if prefer_base and key in pset.terms[term]:
                value, source = pset.terms[term][key], "base"
            if value is None:
                value, source, _ = database.lookup ( term, key )
            if value is None and fallback_terms and key in fallback_terms.get ( term, { } ):
                value, source = fallback_terms[term][key], "estimated"
            if value is None and term == "outofplane":
                value, source = 0.0, "zero"
                report["zero_outofplane"].append ( key )
            if value is None:
                report["missing"][term].append ( key )
                continue
            values[term][key] = value
            report["sources"][term][source] += 1
            report["key_source"][term][key] = source
            report["key_value"][term][key] = value
    # . Urey-Bradley terms (water) only from the base set.
    for ( i, j, k ) in system.connectivity.angleIndices:
        key = key_angle ( types[i], types[j], types[k] )
        if key in pset.terms["ureybradley"]:
            values["ureybradley"][key] = pset.terms["ureybradley"][key]

    # . Atom types and patterns: the base ones plus extras.
    atom_types = dict ( pset.raw["atomTypes"] )
    atom_types["Parameter Values"] = list ( pset.raw["atomTypes"]["Parameter Values"] ) + list ( extra_atom_types or [ ] )
    _dump ( os.path.join ( out_dir, _FILES["atomTypes"] ), atom_types )
    patterns = dict ( pset.raw["patterns"] )
    patterns["Patterns"] = list ( pset.raw["patterns"].get ( "Patterns", [ ] ) ) + list ( extra_patterns or [ ] )
    _dump ( os.path.join ( out_dir, _FILES["patterns"] ), patterns )
    # . Lennard-Jones: every base label (cheap) plus extras.
    lj = dict ( pset.raw["lj"] )
    rows = [ [ label, e, s ] for label, ( e, s ) in pset.lj.items ( ) ]
    for label, ( e, s ) in ( extra_lj or { } ).items ( ):
        rows.append ( [ label, e, s ] )
    lj["Parameter Values"] = rows
    _dump ( os.path.join ( out_dir, _FILES["lj"] ), lj )
    # . Bonded terms.
    def container ( name, rows ):
        mapping = dict ( pset.raw[name] )
        mapping["Parameter Values"] = rows
        _dump ( os.path.join ( out_dir, _FILES[name] ), mapping )
    container ( "bond",  [ [ k[0], k[1], v[0], v[1] ] for k, v in values["bond"].items ( ) ] )
    container ( "angle", [ [ k[0], k[1], k[2], v[0], v[1] ] for k, v in values["angle"].items ( ) ] )
    dihedral_rows = [ ]
    for k, terms in values["dihedral"].items ( ):
        for ( fc, n, phase ) in ( terms or [ ( 0.0, 0, 0.0 ) ] ):
            dihedral_rows.append ( [ k[0], k[1], k[2], k[3], fc, n, phase ] )
    container ( "dihedral", dihedral_rows )
    container ( "outofplane", [ [ k[0], k[1], k[2], k[3], v ] for k, v in values["outofplane"].items ( ) ] )
    if "ureybradley" in pset.raw:
        container ( "ureybradley", [ [ k[0], k[1], k[2], v[0], v[1] ] for k, v in values["ureybradley"].items ( ) ] )
    # . MM sequence components (ligands) of the base set.
    components = os.path.join ( pset.path, "components" )
    if os.path.isdir ( components ) and os.path.abspath ( components ) != os.path.abspath ( os.path.join ( out_dir, "components" ) ):
        shutil.copytree ( components, os.path.join ( out_dir, "components" ), dirs_exist_ok = True )
    write_sources ( out_dir, report["key_source"] )
    with open ( os.path.join ( out_dir, "provenance.json" ), "w" ) as handle:
        json.dump ( { "base": pset.path, "sources": { t: dict ( c ) for t, c in report["sources"].items ( ) },
                      "missing": { t: [ list ( k ) for k in v ] for t, v in report["missing"].items ( ) },
                      "zero_outofplane": [ list ( k ) for k in report["zero_outofplane"] ] }, handle, indent = 1 )
    return report


#=============================================================================
#  Database report
#=============================================================================
def database_report_markdown ( database = None, base = "protein" ):
    """ Markdown text: coverage per source and every conflict. """
    database = database or build_class_database ( base )
    pset = PDynamoSet ( pdynamo_set_path ( base ) )
    lines = [ "# OPLS class-level parameter database", "",
              "Generated by `pdynamo/opls/parameters.py` (`database_report_markdown()`).", "",
              "Precedence: " + "  >  ".join ( database.sources ) + ". The first source that has a key wins;",
              "different values from the other sources are listed below as conflicts (they are NOT used).", "",
              "## Coverage (class-level keys)", "",
              "| Term | pDynamo `{}` (all labels) | Merged database | from each source |".format ( os.path.basename ( pset.path ) ),
              "|---|---:|---:|---|" ]
    for term in TERMS:
        per_source = Counter ( e["source"] for e in database.entries[term].values ( ) )
        lines.append ( "| {} | {} | {} | {} |".format ( term, len ( pset.terms[term] ), len ( database.entries[term] ),
                       ", ".join ( "{} {}".format ( s, n ) for s, n in per_source.items ( ) ) ) )
    lines += [ "", "Lennard-Jones/charge data of {} Tinker atom types are kept for ligand typing.".format ( len ( database.lj_types ) ), "" ]
    conflicts = database.conflicts ( )
    lines += [ "## Conflicts", "" ]
    for term in TERMS:
        items = conflicts[term]
        lines.append ( "### {} ({} keys)".format ( term, len ( items ) ) )
        lines.append ( "" )
        if not items:
            lines += [ "None.", "" ]
            continue
        lines.append ( "| Key | Used (source) | Other values |" )
        lines.append ( "|---|---|---|" )
        for key, source, value, others in items:
            other_text = "; ".join ( "{}: {}".format ( s, v ) for s, v in others )
            lines.append ( "| {} | {} ({}) | {} |".format ( "-".join ( key ), value, source, other_text ) )
        lines.append ( "" )
    return "\n".join ( lines )


#=============================================================================
#  Per-key provenance file (read by the parameter quality window)
#=============================================================================
SOURCES_FILE = "sources.json"

def write_sources ( folder, key_source ):
    """ sources.json: { term: [ [ [ label, ... ], source ], ... ] }. """
    data = { term: [ [ list ( key ), source ] for key, source in items.items ( ) ] for term, items in key_source.items ( ) }
    with open ( os.path.join ( folder, SOURCES_FILE ), "w" ) as handle:
        json.dump ( data, handle, indent = 0 )

def read_sources ( folder ):
    """ { term: { key tuple: source } } (empty when the file is absent). """
    path = os.path.join ( folder, SOURCES_FILE )
    if not os.path.exists ( path ): return { t: { } for t in TERMS }
    with open ( path ) as handle:
        data = json.load ( handle )
    out = { t: { } for t in TERMS }
    for term, items in data.items ( ):
        out.setdefault ( term, { } )
        for key, source in items: out[term][tuple ( key )] = source
    return out
