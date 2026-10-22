#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: carry atom-indexed session data over to a rebuilt system
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      [EN] 2026-10-18, user's request ("se eu precisar reeditar um sistema
#      ela sempre cria um sistema novo"): a tool that rebuilds a system
#      (Prepare OPLS System) can now replace it in place. The EasyHybrid data
#      that stores atom indices (QC region, fixed atoms, selections,
#      restraints, custom colors) is translated with an old -> new index map
#      instead of being dropped; whatever cannot be translated is reported.
#
#      New atoms (hydrogens, OXT) join a set when every atom they are bonded
#      to is in it -- the hydrogens of a QC residue stay in the QC region.
#
#      system_history() tells whether a system already has results that a
#      replacement would orphan (jobs, trajectories, other coordinate sets).
#

import copy

# . e_* attributes that hold atom indices and are NOT translated (reset)
RESET_INDEX_ATTRIBUTES = { 'e_color_segments': dict, 'e_bonds': lambda: None, 'e_sequence': lambda: None,
                           'e_charges_backup': list }


def _neighbours ( system ):
    """ { index: [ bonded indices ] } from system.connectivity. """
    nb = { }
    connectivity = getattr ( system, 'connectivity', None )
    for bond in ( getattr ( connectivity, 'bonds', None ) or [ ] ):
        i, j = bond.node1.index, bond.node2.index
        nb.setdefault ( i, [ ] ).append ( j )
        nb.setdefault ( j, [ ] ).append ( i )
    return nb


class IndexRemapper:
    """ Translates old atom indices to the new system. """

    def __init__ ( self, mapping, new_system ):
        self.mapping = dict ( mapping )                       # old index -> new index
        mapped = set ( self.mapping.values ( ) )
        nb = _neighbours ( new_system )
        # . new atoms (not the image of an old one) and the atoms they are bonded to
        self.followers = { }
        for atom in new_system.atoms:
            k = atom.index
            if k in mapped: continue
            partners = nb.get ( k, [ ] )
            if partners: self.followers[k] = partners

    def index ( self, old ):
        return self.mapping.get ( old )

    def index_set ( self, old_indices, follow = True ):
        """ ( sorted new indices, [ lost old indices ] ) """
        new, lost = set ( ), [ ]
        for i in old_indices:
            j = self.mapping.get ( i )
            if j is None: lost.append ( i )
            else:         new.add ( j )
        if follow and new:
            for k, partners in self.followers.items ( ):
                if all ( p in new for p in partners ):
                    new.add ( k )
        return sorted ( new ), lost


def remap_session_data ( old_system, new_system, mapping ):
    """ Copies every e_* attribute of old_system onto new_system, translating
        the atom-indexed ones (old_system is not changed). Returns a report:
          { "qc": ( n kept, n lost ), "fixed": ..., "selections": { name: ( kept, lost ) },
            "restraints_kept": [names], "restraints_lost": [names], "colors": ( kept, lost ),
            "qc_model": bool (old system had a QC model -- must be defined again) } """
    remap = IndexRemapper ( mapping, new_system )
    report = { 'selections': { }, 'restraints_kept': [ ], 'restraints_lost': [ ],
               'qc_model': getattr ( old_system, 'qcModel', None ) is not None }
    translated = { 'e_qc_table', 'e_qc_residue_table', 'e_fixed_table', 'e_selections',
                   'e_restraints_dict', 'e_custom_colors', 'e_id' } | set ( RESET_INDEX_ATTRIBUTES )
    for name, value in vars ( old_system ).items ( ):
        if name.startswith ( 'e_' ) and name not in translated:
            setattr ( new_system, name, value )
    for name, make in RESET_INDEX_ATTRIBUTES.items ( ):
        setattr ( new_system, name, make ( ) )
    try:
        # . same as pDynamoSession.add_new_system_to_psession()
        new_system.e_charges_backup = list ( new_system.AtomicCharges ( ) ).copy ( )
    except Exception:
        pass

    # . QC region (+ per-residue table, rebuilt from the new sequence)
    old_qc = list ( getattr ( old_system, 'e_qc_table', None ) or [ ] )
    qc, lost = remap.index_set ( old_qc )
    new_system.e_qc_table = qc
    residue_table = { }
    for k in qc:
        # . keyed by residue number, like eSession.build_index_list_from_atom_selection()
        label = getattr ( getattr ( new_system.atoms[k], 'parent', None ), 'label', '' ) or ''
        try:    key = int ( ''.join ( c for c in label.split ( '.' )[-1] if c.isdigit ( ) or c == '-' ) )
        except ValueError: key = 0
        residue_table.setdefault ( key, [ ] ).append ( k )
    new_system.e_qc_residue_table = residue_table if qc else { }
    report['qc'] = ( len ( qc ), len ( lost ) )

    fixed, lost = remap.index_set ( getattr ( old_system, 'e_fixed_table', None ) or [ ] )
    new_system.e_fixed_table = fixed
    report['fixed'] = ( len ( fixed ), len ( lost ) )

    selections = { }
    for name, indices in ( getattr ( old_system, 'e_selections', None ) or { } ).items ( ):
        kept, lost = remap.index_set ( indices )
        if kept: selections[name] = kept
        report['selections'][name] = ( len ( kept ), len ( lost ) )
    new_system.e_selections = selections

    restraints = { }
    for name, entry in ( getattr ( old_system, 'e_restraints_dict', None ) or { } ).items ( ):
        entry = copy.copy ( list ( entry ) )
        kind = entry[2] if len ( entry ) > 2 else None
        atoms = entry[3] if len ( entry ) > 3 else None
        if kind == 'position':
            kept, lost = remap.index_set ( atoms or [ ] )
            ok = bool ( kept )
            entry[3] = kept
        else:
            new_atoms = [ remap.index ( a ) for a in ( atoms or [ ] ) ]
            ok = bool ( new_atoms ) and all ( a is not None for a in new_atoms )
            entry[3] = new_atoms
        if ok:
            restraints[name] = entry
            report['restraints_kept'].append ( name )
        else:
            report['restraints_lost'].append ( name )
    new_system.e_restraints_dict = restraints

    colors, n_kept, n_lost = [ ], 0, 0
    for item in ( getattr ( old_system, 'e_custom_colors', None ) or [ ] ):
        item = dict ( item )
        kept, lost = remap.index_set ( item.get ( 'indexes', [ ] ) )
        if kept:
            item['indexes'] = kept
            colors.append ( item )
            n_kept += 1
        if lost: n_lost += 1
    new_system.e_custom_colors = colors
    report['colors'] = ( n_kept, n_lost )
    return report


def remap_summary ( report ):
    """ Lines describing what was carried over / lost (empty when nothing
        index-based existed). """
    lines = [ ]
    kept, lost = report.get ( 'qc', ( 0, 0 ) )
    if kept or lost:
        lines.append ( 'QC region: {} atoms kept{}.'.format ( kept, ', {} removed'.format ( lost ) if lost else '' ) )
    if report.get ( 'qc_model' ):
        lines.append ( 'The QC model must be defined again (the region is kept).' )
    kept, lost = report.get ( 'fixed', ( 0, 0 ) )
    if kept or lost:
        lines.append ( 'Fixed atoms: {} kept{}.'.format ( kept, ', {} removed'.format ( lost ) if lost else '' ) )
    for name, ( kept, lost ) in report.get ( 'selections', { } ).items ( ):
        if lost or not kept:
            lines.append ( 'Selection "{}": {} atoms kept, {} removed{}.'.format (
                    name, kept, lost, ' (selection dropped)' if not kept else '' ) )
    if report.get ( 'restraints_kept' ):
        lines.append ( 'Restraints kept: {}.'.format ( ', '.join ( report['restraints_kept'] ) ) )
    if report.get ( 'restraints_lost' ):
        lines.append ( 'Restraints removed (atoms no longer exist): {}.'.format ( ', '.join ( report['restraints_lost'] ) ) )
    kept, lost = report.get ( 'colors', ( 0, 0 ) )
    if lost:
        lines.append ( 'Custom colors: {} kept, {} lost atoms.'.format ( kept, lost ) )
    return lines


def report_has_losses ( report ):
    return bool ( report.get ( 'qc', ( 0, 0 ) )[1] or report.get ( 'fixed', ( 0, 0 ) )[1]
                  or report.get ( 'restraints_lost' ) or report.get ( 'qc_model' )
                  or any ( lost for ( _k, lost ) in report.get ( 'selections', { } ).values ( ) ) )


def system_history ( main, e_id ):
    """ What a replacement of psystem[e_id] would leave behind:
          { "jobs": n entries of e_job_history, "objects": [ names of the OTHER
            coordinate objects of the system (trajectories, scans, results) ],
            "free": no jobs and no other objects } """
    system = main.p_session.psystem.get ( e_id )
    jobs = len ( getattr ( system, 'e_job_history', None ) or { } ) if system is not None else 0
    objects = [ ]
    for vobject in main.vm_session.vm_objects_dic.values ( ):
        if getattr ( vobject, 'e_id', None ) != e_id or getattr ( vobject, 'is_surface', False ):
            continue
        objects.append ( vobject )
    others = [ getattr ( v, 'name', '?' ) for v in objects[1:] ]
    return { 'jobs': jobs, 'objects': others, 'n_objects': len ( objects ),
             'free': jobs == 0 and len ( objects ) <= 1 }
