#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Build EasyHybrid-registered VismolObjects directly from PDBQT
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Maintainer:
#      Fernando Bachega <ferbachega@gmail.com> or <easyhybrid3@gmail.com>
#
#  Description:
#      Loads PDBQT (AutoDock/Vina) receptors and ligand-docking poses
#      directly into EasyHybrid's VisMol view -- no temporary PDB file
#      on disk, no OpenBabel re-invocation, no pDynamo ImportSystem()
#      call reading a file.
#
#      IMPORTANT ARCHITECTURE NOTE (discovered while building this):
#      EasyHybrid3's own VismolSession subclass (src/gui/eSession.py,
#      _add_vismol_object()) requires EVERY VismolObject to carry an
#      `e_id` pointing to a REAL pDynamo System registered in
#      p_session.psystem -- unlike the vendored vismol library's own
#      generic, pDynamo-free VismolObject construction path
#      (vismol/utils/parser.py:_load_pdb_file), which is NOT reachable
#      from this app's actual runtime session. So "no temp file" here
#      means building a pDynamo System ENTIRELY IN MEMORY (proven
#      precedent: src/gui/windows/builder/empty_object.py's
#      _build_pdynamo_system_from_vismol_object(), which builds a
#      System from a VismolObject's own atoms with
#      System.FromConnectivity() -- but that path only produces a
#      dummy single-residue "UNK.1" sequence, losing chain/residue
#      identity, which the user explicitly requires be preserved here).
#      Instead this module uses pMolecule.Sequence.FromAtomPaths(),
#      which builds a REAL chain/residue/atom hierarchy from plain
#      "chain:resName.resSeq:atomName" path strings -- no file I/O,
#      confirmed working against a real 2746-atom receptor PDBQT.
#      System.FromSequence() then wraps it; the resulting System has no
#      mmModel/qcModel (confirmed this is NOT required for display --
#      System.AtomicCharges() returns None gracefully, and
#      eSession._add_vismol_object()/update_restaint_representation()/
#      _apply_custom_colors_to_vobject() only touch plain 'e_*'
#      attributes, not any pDynamo energy-model internals), so the
#      existing, already-tested EasyHybrid pipeline
#      (p_session._add_vismol_object_to_easyhybrid_session(), which
#      calls _build_vobject_from_pdynamo_system() internally) builds
#      the actual VismolObject/Chain/Residue/Atom tree AND falls back
#      to VismolObject.find_bonded_and_nonbonded_atoms() (distance/
#      covalent-radius bond perception) for bonds, since there is no
#      mmState connectivity to read bonds from -- exactly right for a
#      receptor/ligand with no force field attached.
#
import os
from dataclasses import dataclass

import numpy as np


@dataclass
class DockingResult:
    """ One imported docking pose, with its receptor pairing made
        explicit -- ligand_frame and receptor_frame are independent
        indices into two DIFFERENT VismolObjects (see VismolObject.
        pin_frame()/unpin_frame()) and must never be assumed equal.
        receptor_object/receptor_frame are None/-1 if the receptor was
        never loaded into VisMol (load_receptor_pdbqt()/
        load_receptor_ensemble_pdbqt() were never called for this
        docking run) -- receptor_id is always available regardless,
        straight from the docking results table.
    """
    ligand_object: object
    ligand_frame: int
    receptor_object: object
    receptor_frame: int
    pose_rank: int
    score: float
    receptor_id: str


def _build_pdynamo_system_from_pdbqt_topology(pdbqt_parser_module, topology, label, dedupe_policy='drop'):
    """ Builds a pDynamo System entirely in memory from a pdbqt_parser
        topology list -- no file I/O. See this module's docstring for
        why this route (Sequence.FromAtomPaths), rather than
        System.FromConnectivity() (which loses chain/residue identity)
        or ImportSystem() (which structurally requires a real
        filesystem path -- confirmed: pDynamo3's TextFile base class
        calls os.access()/open() on a path, no stream/StringIO support
        anywhere in that chain).

        `dedupe_policy`: 'drop' (default, for RECEPTORS -- see
        pdbqt_parser.drop_duplicate_atom_names()) or 'rename' (for
        LIGANDS -- see pdbqt_parser.rename_duplicate_atom_names()).
        These two real-world duplicate-atom-name cases need OPPOSITE
        handling: a receptor's collision is redundant/overlapping data
        to discard, while a ligand's collision is often several
        genuinely distinct, generically-labeled atoms that must all be
        kept -- confirmed against real fixtures for both (see either
        function's own docstring). Passing the wrong policy either
        leaves a receptor's bond perception confused by overlapping
        duplicate atoms, or silently deletes real ligand atoms.
    """
    # Imported here (not at module top) so this module itself stays
    # importable without pDynamo3 already being on sys.path/loaded --
    # by the time any caller actually calls this function, EasyHybrid
    # has already imported pMolecule/pScientific itself.
    from pMolecule import Atom, System
    from pMolecule.Sequence import Sequence
    from pScientific import PeriodicTable
    from pScientific.Geometry3 import Coordinates3

    # Work on a copy of the per-atom dicts -- both dedupe policies
    # mutate `topology` (dropping atoms and/or renaming 'name'), and
    # the caller's own topology (e.g. the reference used for a later
    # pdbqt_topology_names_match() check) should not be silently
    # changed out from under it.
    topology = [dict(a) for a in topology]
    if dedupe_policy == 'rename':
        pdbqt_parser_module.rename_duplicate_atom_names(topology)
    else:
        pdbqt_parser_module.drop_duplicate_atom_names(topology)

    atoms = [Atom.WithOptions(atomicNumber=PeriodicTable.AtomicNumber(a['element']), label=a['name'])
             for a in topology]
    atom_paths = []
    for a in topology:
        chain = a['chain'] if a['chain'].strip() else 'A'
        atom_paths.append('{}:{}.{}:{}'.format(chain, a['resname'], a['resseq'], a['name']))

    sequence = Sequence.FromAtomPaths(atom_paths, atoms=atoms)
    system = System.FromSequence(sequence)
    system.label = label

    # [EN] IMPORTANT: System.FromSequence() REASSIGNS every atom's own
    # .index attribute to match its position in the resulting
    # system.atoms container -- which is TREE order (atoms grouped by
    # residue/chain), NOT the file/atomPaths order Sequence.FromAtomPaths()
    # originally assigned. Confirmed directly: for a residue whose atoms
    # are NOT contiguous in the source file (e.g. a receptor PDBQT where
    # a residue's sidechain reappears later in the file, interleaved
    # with other residues -- a real, confirmed case in
    # adtVinaFR/7CO0/flexible_receptor_PDBQT/receptor.pdbqt), the same
    # Atom object's .index value visibly changes between right after
    # FromAtomPaths() and after FromSequence() (e.g. 5 -> 3 in a minimal
    # repro). Since EasyHybrid's own _build_vobject_from_pdynamo_system()
    # (session.py) fetches each atom's coordinates via
    # `system.coordinates3[atom.index]`, that FINAL, POST-FromSequence
    # index is the one that matters -- writing coordinates by ORIGINAL
    # topology position instead (as an earlier version of this function
    # did) put each atom's real xyz in the WRONG row whenever its final
    # index differed from its file-order position, silently swapping
    # that atom's displayed position with whichever OTHER atom ended up
    # sharing its old row -- this is what a user reported as "the
    # structure shows wrong bonds, like they're out of order" (verified
    # live: an ARG214 atom was rendered at MET77's own coordinates).
    # `atoms[i]` is the exact same object as `topology[i]` throughout
    # (never replaced/copied), so `atoms[i].index` always reflects
    # whatever the atom's CURRENT, authoritative index is.
    coords = Coordinates3.WithExtent(len(topology))
    for i, a in enumerate(topology):
        row = atoms[i].index
        coords[row, 0] = a['x']
        coords[row, 1] = a['y']
        coords[row, 2] = a['z']
    system.coordinates3 = coords
    return system


def load_pdbqt_as_vobject(p_session, pdbqt_parser_module, topology, name, tag='pdbqt', dedupe_policy='drop'):
    """ Builds a pDynamo System from `topology` (already parsed --
        see pdbqt_parser.parse_pdbqt_topology()/parse_pdbqt_models()),
        registers it (p_session.add_new_system_to_psession()), and
        builds+registers the corresponding VismolObject through
        EasyHybrid's existing, already-tested pipeline
        (p_session._add_vismol_object_to_easyhybrid_session()) -- the
        SAME pipeline every file-based system import already uses, so
        this object behaves identically in the treeview/selection/
        representation machinery to any other loaded system.

        `dedupe_policy`: passed straight through to
        _build_pdynamo_system_from_pdbqt_topology() -- 'drop' for
        receptors (the default here), 'rename' for ligand poses. Ligand
        pose callers (on_button_import_pose_clicked in the docking
        windows) MUST pass dedupe_policy='rename' explicitly.

        Returns the new VismolObject (vm_object.metadata is left as
        an empty dict -- the caller sets it, e.g. {'type': 'receptor',
        ...} or {'type': 'receptor_ensemble', 'frames': {...}}).
    """
    system = _build_pdynamo_system_from_pdbqt_topology(pdbqt_parser_module, topology, name, dedupe_policy=dedupe_policy)
    p_session.add_new_system_to_psession(system=system, name=name, tag=tag)
    # Every other caller of _add_vismol_object_to_easyhybrid_session()
    # in session.py first registers the system with the main treeview
    # (add_new_system_to_treeview() sets up main.system_treeview_iters
    # [e_id], which add_vismol_object_to_treeview() -- called inside
    # _add_vismol_object_to_easyhybrid_session() -- looks up immediately;
    # skipping this step raises KeyError there, confirmed directly).
    p_session.main.main_treeview.add_new_system_to_treeview(system)
    vm_object = p_session._add_vismol_object_to_easyhybrid_session(system, True, name=name)
    vm_object.metadata = {}
    return vm_object


def append_coord_frame(vm_object, coords):
    """ Appends one more coordinate frame to an already-built
        VismolObject, e.g. a second docking pose of the same ligand or
        a second receptor-ensemble structure.

        No topology check here -- the caller is responsible for having
        already verified topology match (pdbqt_parser.
        pdbqt_topology_names_match()) before calling this. Does NOT
        rebuild the pDynamo System or re-run bond perception -- bonds
        stay derived from frame 0 only, matching this codebase's
        existing contract for appended trajectory frames (see project
        memory: vismol dynamic_bonds bug -- do not apply any QC/
        dynamic representation to this object until ALL of its frames
        are final).
    """
    vm_object.frames = np.vstack(
        (vm_object.frames, coords[np.newaxis, :, :].astype(vm_object.frames.dtype)))


def load_receptor_pdbqt(p_session, pdbqt_parser_module, pdbqt_path, name=None):
    """ Loads a single receptor PDBQT file directly as a 1-frame
        VismolObject (frame_number 0), named after the receptor unless
        `name` is given.
    """
    topology = pdbqt_parser_module.parse_pdbqt_topology(pdbqt_path)
    name = name or os.path.splitext(os.path.basename(pdbqt_path))[0]
    vm_object = load_pdbqt_as_vobject(p_session, pdbqt_parser_module, topology, name, tag='receptor_pdbqt')
    vm_object.metadata = {"type": "receptor", "source": pdbqt_path}
    return vm_object


def load_receptor_ensemble_pdbqt(p_session, pdbqt_parser_module, pdbqt_paths, name=None):
    """ Loads an ordered list of per-frame receptor PDBQT files as ONE
        VismolObject, one frame per file, in the given order. The
        first file's topology is authoritative; every subsequent file
        must match it atom-for-atom (name/resname/chain/resseq, AFTER
        drop_duplicate_atom_names() is applied to each independently
        -- a deterministic, structure-driven drop, so two files of
        the genuinely same receptor still match after it) -- ensembles
        are, by construction, the SAME receptor (NMR/MD snapshots
        etc.), so a mismatch here is a genuine user error, not
        something to silently skip (unlike ligand-pose mismatches in
        the docking-results import path, which legitimately get
        skipped-and-reported since different docked files there can be
        genuinely different ligands/conformers).
    """
    topologies = []
    for path in pdbqt_paths:
        topo = pdbqt_parser_module.parse_pdbqt_topology(path)
        topo = [dict(a) for a in topo]
        pdbqt_parser_module.drop_duplicate_atom_names(topo)
        topologies.append(topo)

    for i, t in enumerate(topologies[1:], start=1):
        if not pdbqt_parser_module.pdbqt_topology_names_match(topologies[0], t):
            raise ValueError('Receptor ensemble atom mismatch in {}'.format(pdbqt_paths[i]))

    name = name or 'receptor_ensemble'
    vm_object = load_pdbqt_as_vobject(p_session, pdbqt_parser_module, topologies[0], name, tag='receptor_ensemble')
    frames_meta = {0: {
        'receptor_id': os.path.splitext(os.path.basename(pdbqt_paths[0]))[0],
        'source': pdbqt_paths[0]}}
    for i, (path, topo) in enumerate(zip(pdbqt_paths[1:], topologies[1:]), start=1):
        coords = np.array([[a['x'], a['y'], a['z']] for a in topo], dtype=vm_object.frames.dtype)
        append_coord_frame(vm_object, coords)
        frames_meta[i] = {
            'receptor_id': os.path.splitext(os.path.basename(path))[0],
            'source': path}

    vm_object.metadata = {
        'type': 'receptor_ensemble',
        'source': os.path.dirname(pdbqt_paths[0]),
        'frames': frames_meta,
    }
    return vm_object


def create_or_update_docking_box(p_session, vm_session, pdbqt_parser_module, vm_object, box, name,
                                  color=(0.2, 0.85, 0.3)):
    """ Creates (if `vm_object` is None) or updates (otherwise) a
        wireframe box in VisMol representing an AutoDock Vina/
        AutoDock-GPU docking search box -- reuses vismol's OWN "unit
        cell" mechanism (VismolObject.set_cell() +
        VismolSession.show_cell()), confirmed to already do exactly
        what's needed here (build 8 box-corner vertices + 12 wireframe
        edges from a/b/c/alpha/beta/gamma, GL_LINES rendering) rather
        than writing new box-drawing code.

        `box`: a dict with center_x/y/z, size_x/y/z (the same shape
        vina_runner.compute_box_from_coordinates() returns).

        Pass the SAME `vm_object` back in on a later call (instead of
        None) to move/resize the existing box in place rather than
        creating a new one every time the user changes Center/Size --
        confirmed live precedent (src/gui/windows/setup/edit_cell.py):
        calling set_cell() again on an already-cell'd object, then
        show_cell() again, is the established way to update one.

        Why this needs a real (if minimal) VismolObject rather than a
        lighter-weight "geometric object": confirmed this app's render
        loop (vismol_glcore.py's render()) only draws representations
        for objects living in vm_session.vm_objects_dic --
        vm_geometric_object_dic (picking/dash markers) is never drawn
        generically, only for a few hardcoded picking-mode
        representation names, so a box object living there would never
        actually appear on screen. And every object in vm_objects_dic
        in THIS app (not the generic vismol library) needs a real
        `e_id` backed by a pDynamo System (see this module's own
        docstring on load_pdbqt_as_vobject()) -- so the box object
        carries exactly ONE dummy marker atom (element 'X', pDynamo
        atomicNumber -1, confirmed handled gracefully -- EasyHybrid's
        own _get_atom_info_from_pdynamo_atom_obj() falls back to
        symbol "X" on any lookup failure) positioned at the box's own
        center, purely so it's a valid, registerable structure. Its
        default "lines"/"nonbonded" representations are deactivated
        right after creation so only the wireframe cell box itself is
        visible, never the marker atom.
    """
    center = (box['center_x'], box['center_y'], box['center_z'])
    size = (box['size_x'], box['size_y'], box['size_z'])

    if vm_object is None:
        topology = [{
            'serial': 1, 'name': 'X', 'resname': 'BOX', 'chain': 'X', 'resseq': 1,
            'x': center[0], 'y': center[1], 'z': center[2],
            'occupancy': 0.0, 'bfactor': 0.0, 'charge': 0.0,
            'autodock_type': '', 'element': 'X',
        }]
        vm_object = load_pdbqt_as_vobject(p_session, pdbqt_parser_module, topology, name, tag='docking_box')
        for rep_name in ('lines', 'nonbonded'):
            representation = vm_object.representations.get(rep_name)
            if representation is not None:
                representation.active = False

    # set_cell() always builds the box with one corner at the LOCAL
    # origin (0,0,0) -- there is no center parameter (confirmed
    # directly in vismol_object.py). Shifting cell_coordinates by
    # (center - size/2) afterward re-centers it at the actual docking
    # box center instead of leaving a corner pinned to the origin.
    vm_object.set_cell(size[0], size[1], size[2], 90, 90, 90, color=list(color))
    offset = np.array(
        [center[0] - size[0] / 2.0, center[1] - size[1] / 2.0, center[2] - size[2] / 2.0], dtype=np.float32)
    vm_object.cell_coordinates[0] += offset

    vm_session.show_cell(vm_object)
    return vm_object
