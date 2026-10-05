#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: quality assessment of the OPLS parameters of ligands
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      For every atom of a parametrized ligand (prep.prepare_opls_system()),
#      from the per-key provenance of the materialized parameter set:
#
#        - terms: how many bonded terms contain the atom and how many of
#          them were not found in the OPLS database (estimated, generic or
#          zero out-of-plane);
#        - rotatable torsions without parameters: dihedrals left at zero
#          whose central bond is single, acyclic and between two non-
#          terminal atoms -- the candidates for a QC torsion scan;
#        - charge: outside the range of the charges of its class among the
#          OPLS-AA 2008 types (with a margin) -> outlier;
#        - level: "red"    estimated bond/angle, or central atom of a
#                          rotatable torsion without parameters;
#                 "yellow" only generic terms (aromatic ring torsion,
#                          generic or zero improper), or a charge outlier;
#                 "green"  every term from the database, charge in range.
#
#      Per ligand: level counts, torsion scan candidates, estimated
#      bonds/angles, charge outliers and plain-text recommendations.
#

from collections import defaultdict, OrderedDict

import numpy as np

CHARGE_MARGIN = 0.15        # e, tolerance outside the class charge range
LEVELS = ( "green", "yellow", "red" )

_CLASS_CHARGES = None

def class_charge_ranges ( ):
    """ { class: ( min, max, mean ) } of the OPLS-AA 2008 type charges
        (united-atom and GBSA types out). """
    global _CLASS_CHARGES
    if _CLASS_CHARGES is None:
        from pdynamo.opls.parameters import TINKER_FILES
        from pdynamo.opls.tinker_prm import TinkerPRM
        prm = TinkerPRM ( dict ( TINKER_FILES )["oplsaa08"] )
        values = defaultdict ( list )
        for t, data in prm.atoms.items ( ):
            if "(UA)" in data["description"] or "GBSA" in data["description"]: continue
            if t in prm.charges: values[data["symbol"]].append ( prm.charges[t] )
        _CLASS_CHARGES = { c: ( min ( v ), max ( v ), float ( np.mean ( v ) ) ) for c, v in values.items ( ) }
    return _CLASS_CHARGES


def ligand_atom_groups ( system, ligand_names ):
    """ { residue path: [ atom indices ] } of the components whose name is a
        parametrized ligand (monoatomic ions left out). """
    groups = OrderedDict ( )
    for entity in system.sequence.children:
        for component in entity.children:
            if component.genericLabel in ligand_names and len ( component.children ) > 1:
                groups[component.path] = [ a.index for a in component.children ]
    return groups


def _in_ring ( adjacency, j, k, allowed ):
    """ True when j and k stay connected without the j-k bond. """
    seen, stack = { j }, [ j ]
    while stack:
        u = stack.pop ( )
        for v in adjacency[u]:
            if ( u == j and v == k ) or ( u == k and v == j ) or v not in allowed: continue
            if v == k: return True
            if v not in seen:
                seen.add ( v ); stack.append ( v )
    return False


def assess ( system, types, charges, materialized, ligand_names ):
    """ Quality report (dict) of the ligands of a prepared system. """
    from pdynamo.opls import parameters as P
    connectivity = system.connectivity
    atoms = list ( system.atoms )
    adjacency = defaultdict ( set )
    bond_type = { }
    for bond in connectivity.bonds:
        i, j = bond.node1.index, bond.node2.index
        adjacency[i].add ( j ); adjacency[j].add ( i )
        bond_type[( min ( i, j ), max ( i, j ) )] = bond.type.name if hasattr ( bond.type, "name" ) else str ( bond.type )
    source = materialized["key_source"]
    value  = materialized["key_value"]
    groups = ligand_atom_groups ( system, ligand_names )
    ranges = class_charge_ranges ( )
    report = OrderedDict ( )
    for path, indices in groups.items ( ):
        members = set ( indices )
        per_atom = { i: { "terms": 0, "estimated": defaultdict ( int ) } for i in indices }
        def count ( term, key, tuple_ ):
            src = source[term].get ( key )
            flag = None
            if src == "estimated":
                flag = term
                if term == "dihedral":
                    flag = "dihedral_ring" if value[term][key] else "dihedral_zero"
            elif src == "zero":
                flag = "outofplane_zero"
            for a in tuple_:
                if a in per_atom:
                    per_atom[a]["terms"] += 1
                    if flag: per_atom[a]["estimated"][flag] += 1
        for ( i, j ) in connectivity.bondIndices:
            if i in members or j in members: count ( "bond", P.key_bond ( types[i], types[j] ), ( i, j ) )
        for ( i, j, k ) in connectivity.angleIndices:
            if members & { i, j, k }: count ( "angle", P.key_angle ( types[i], types[j], types[k] ), ( i, j, k ) )
        zero_torsions = defaultdict ( list )        # central bond -> dihedrals
        for ( i, j, k, l ) in connectivity.dihedralIndices:
            if not members & { i, j, k, l }: continue
            key = P.key_dihedral ( types[i], types[j], types[k], types[l] )
            count ( "dihedral", key, ( i, j, k, l ) )
            if source["dihedral"].get ( key ) == "estimated" and not value["dihedral"][key]:
                central = ( min ( j, k ), max ( j, k ) )
                rotatable = ( bond_type.get ( central ) == "Single" and len ( adjacency[j] ) > 1 and len ( adjacency[k] ) > 1
                              and not _in_ring ( adjacency, j, k, members ) )
                if rotatable: zero_torsions[central].append ( ( i, j, k, l ) )
        for ( i, j, k, l ) in connectivity._IndicesNeighbor3 ( ):
            if i in members: count ( "outofplane", P.key_outofplane ( types[i], types[j], types[k], types[l] ), ( i, j, k, l ) )
        torsion_atoms = { a for central in zero_torsions for a in central }
        atoms_out = OrderedDict ( )
        levels = defaultdict ( int )
        outliers = [ ]
        for i in indices:
            cls = P.atom_class ( types[i] )
            q = float ( charges[i] )
            low, high, mean = ranges.get ( cls, ( None, None, None ) )
            outlier = low is not None and ( q < low - CHARGE_MARGIN or q > high + CHARGE_MARGIN )
            est = per_atom[i]["estimated"]
            if est.get ( "bond" ) or est.get ( "angle" ) or i in torsion_atoms: level = "red"
            elif est or outlier: level = "yellow"
            else: level = "green"
            levels[level] += 1
            name = atoms[i].label
            if outlier:
                outliers.append ( "{} ({}, q {:+.3f}; OPLS-AA {} range {:+.3f} to {:+.3f})".format ( name, cls, q, cls, low, high ) )
            atoms_out[i] = { "name": name, "class": cls, "type": types[i], "charge": round ( q, 4 ),
                             "terms": per_atom[i]["terms"], "estimated": dict ( est ),
                             "charge_outlier": bool ( outlier ), "level": level }
        def label ( t ): return "-".join ( atoms[a].label for a in t )
        scans = [ { "bond": label ( central ), "dihedral": label ( dihedrals[0] ), "n_dihedrals": len ( dihedrals ) }
                  for central, dihedrals in zero_torsions.items ( ) ]
        est_terms = sorted ( { "{} {}".format ( term, "-".join ( key ) )
                               for term in ( "bond", "angle" ) for key, src in source[term].items ( )
                               if src == "estimated" and any ( P.atom_class ( l ) for l in key ) } )
        recommendations = [ ]
        if scans:
            recommendations.append ( "{} rotatable bond(s) without torsion parameters ({}): run a QC dihedral scan "
                                     "and fit V1-V3 before conformational sampling.".format (
                                     len ( scans ), ", ".join ( s["bond"] for s in scans ) ) )
        n_est = sum ( 1 for a in atoms_out.values ( ) if a["estimated"].get ( "bond" ) or a["estimated"].get ( "angle" ) )
        if n_est:
            recommendations.append ( "{} atom(s) in estimated bonds/angles (equilibrium values from the starting "
                                     "geometry): optimize the ligand at the QC level (e.g. GFN2-xTB) and compare "
                                     "the MM geometry; refit if the deviations are large.".format ( n_est ) )
        if outliers:
            recommendations.append ( "{} charge(s) outside the OPLS-AA range of their class: check the protonation "
                                     "state and the atom classes.".format ( len ( outliers ) ) )
        if any ( a["estimated"].get ( "dihedral_ring" ) for a in atoms_out.values ( ) ):
            recommendations.append ( "Ring torsions use the generic aromatic term (V2 = 7.25 kcal/mol): adequate for "
                                     "planarity, not for ring puckering." )
        if not recommendations:
            recommendations.append ( "All bonded terms come from the OPLS database and the charges are within the "
                                     "OPLS-AA ranges of their classes." )
        report[path] = { "atoms": atoms_out, "levels": { lv: levels.get ( lv, 0 ) for lv in LEVELS },
                         "torsion_scans": scans, "charge_outliers": outliers,
                         "recommendations": recommendations }
    return report


def markdown ( quality ):
    """ Markdown section for prep.report_markdown(). """
    lines = [ "", "## Parameter quality (ligands)", "",
              "Levels: **green** every term from the OPLS database and the charge within its class range; "
              "**yellow** generic terms (aromatic ring torsion, generic/zero improper) or a charge outlier; "
              "**red** estimated bonds/angles or a rotatable bond without torsion parameters. "
              "In EasyHybrid every parameter, with its level, is listed (and can be edited) in Extras > OPLS > OPLS Parameters.", "" ]
    for path, data in quality.items ( ):
        lv = data["levels"]
        total = sum ( lv.values ( ) ) or 1
        lines.append ( "### {}: {} green, {} yellow, {} red ({:.0f} % green)".format (
                path, lv["green"], lv["yellow"], lv["red"], 100.0 * lv["green"] / total ) )
        lines.append ( "" )
        red = [ "{} ({})".format ( a["name"], a["class"] ) for a in data["atoms"].values ( ) if a["level"] == "red" ]
        if red: lines.append ( "- Red atoms: " + ", ".join ( red ) )
        for scan in data["torsion_scans"]:
            lines.append ( "- QC scan candidate: bond {} (e.g. dihedral {}, {} dihedral(s) with zero parameters)".format (
                    scan["bond"], scan["dihedral"], scan["n_dihedrals"] ) )
        for outlier in data["charge_outliers"]:
            lines.append ( "- Charge outlier: " + outlier )
        lines.append ( "" )
        lines.append ( "Recommendations:" )
        lines += [ "- " + r for r in data["recommendations"] ]
        lines.append ( "" )
    return lines
