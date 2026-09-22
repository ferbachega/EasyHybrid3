#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Molecule Builder -- atom-level editing operations
#
#  Description:
#      Second building block of the "Builder" tool. add_atom() below adds a
#      single atom to a VismolObject, whether it currently has zero atoms
#      (a builder-created empty object, see empty_object.py) or already has
#      some (a real, file-loaded molecule -- this function does not care
#      which; it only touches the same internal structures either way).
#
#      Design notes / what this reuses from the existing codebase (all
#      read from source, not guessed):
#        - vismol_object.frames is the single source of truth for atom
#          coordinates: Atom.coords() reads self.vm_object.frames[frame,
#          self.atom_id] directly, NOT some per-atom cached position (the
#          Atom.pos attribute is described in atom.py itself as just "the
#          coordinates of the first frame", i.e. a convenience snapshot,
#          not the live value). Adding an atom therefore means growing
#          .frames along its atom axis (axis=1) for every existing frame,
#          not just setting atom.pos.
#        - Chain / Residue containers are created lazily, matching the
#          exact same pattern used by the real file loaders (e.g.
#          vismol.utils.XYZFiles.load_xyz_file): if the requested chain_id
#          / residue index don't exist yet on this object, create them
#          with sensible defaults ("A" / 1 / "UNK") before creating the
#          Atom itself.
#        - Atom.unique_id / _generate_atom_unique_color_id() are what the
#          picking/selection system uses (an RGB-encoded unique ID, see
#          _generate_atom_unique_color_id in atom.py) -- every atom loaded
#          from a file gets these assigned via vm_session.atom_id_counter,
#          so builder-created atoms do the same to stay selectable/
#          pickable like any other atom.
#        - Bond (re)detection reuses VismolObject.find_bonded_and_
#          nonbonded_atoms() UNCHANGED -- the same grid/covalent-radius
#          based auto-detection already used "at the object's genesis"
#          when loading a real file (see that method's own docstring).
#          This also populates self.non_bonded_atoms, which is what the
#          "nonbonded" representation (dots/spheres for atoms with no
#          bonds -- necessary for a single freshly-placed atom to be
#          visible at all, since "lines"/"sticks" only draw BONDS) needs.
#          self.cov_radii_array / self.electronegativity_array / .
#          index_bonds are reset to None first so that method recomputes
#          them fresh at the new atom count instead of using stale,
#          wrong-sized arrays from before this atom was added.
#
#      NOT implemented yet (left for later steps): removing atoms, moving
#      an existing atom, adding an explicit bond without relying on
#      distance-based auto-detection, undo/redo.
#
import numpy as np
from vismol.model.atom import Atom
from vismol.model.chain import Chain
from vismol.model.residue import Residue
from util.debug import dprint


def add_atom ( vismol_object, symbol, x, y, z, name = None,
               chain_id = "A", resi = 1, resn = "UNK",
               bonded_to = None, bond_order = 1, aromatic = False,
               recompute_bonds = True, update_representation = True ):
    """ Adds a single atom to vismol_object at position (x, y, z), in the
    same coordinate units/frame convention already used by the rest of
    the object (Angstrom, matching every file-loaded VismolObject).

    Works whether vismol_object currently has zero atoms (e.g. one just
    created by empty_object.create_empty_vismol_object()) or already has
    atoms -- either way the new atom is simply appended after the
    existing ones (atom_id = current atom count).

    Parameters
    ----------
    vismol_object : VismolObject
    symbol   : element symbol, e.g. "C", "N", "O", "H" -- used to look up
               default color/radii/electronegativity via the periodic
               table (see Atom._init_color / _init_vdw_rad / etc.).
    x, y, z  : position of the new atom, in Angstrom.
    name     : atom name (e.g. "CA", "OXT"...); defaults to `symbol` if
               not given, which is the common case for a freshly-placed,
               not-yet-part-of-a-residue-template atom.
    chain_id, resi, resn : which Chain/Residue this atom belongs to.
               Created on the fly, with sensible defaults, if they don't
               exist yet on this object -- same pattern the file loaders
               already use.
    bonded_to : atom_id, or list of atom_ids, this new atom should be
               EXPLICITLY bonded to (via the same manual_bonds mechanism
               add_bond() uses -- see _reapply_manual_bonds()). This is
               now the ONLY way a newly-added atom ends up bonded to
               anything: distance-based auto-detection was intentionally
               turned OFF (see below), so an atom placed near another one
               by coincidence does NOT become bonded to it just from
               proximity anymore.
    bond_order : order (1/2/3) recorded in manual_bond_orders for EVERY
               pair in `bonded_to` -- ignored if `bonded_to` is None.
               Defaults to 1 (plain single bond), matching the previous,
               unconditional behaviour before this parameter existed.
    aromatic : if True, ALSO records every `bonded_to` pair in
               manual_aromatic_bonds (see add_bond()'s own docstring for
               why aromaticity is tracked as a flag separate from
               bond_order -- mirrors pDynamo3's own Bond model, where
               "Aromatic bond types are not defined. Instead bonds
               themselves are flagged if they form part of an aromatic
               system." -- pMolecule/Bond.py's own comment).
    recompute_bonds : if True (default), rebuilds bonds/topology/molecule
               grouping from vismol_object.manual_bonds (via
               _reapply_manual_bonds()) after adding the atom. Set False
               only if you're about to add several atoms in a row and
               want to pay for this once, at the end (call
               _reapply_manual_bonds(vismol_object) yourself afterwards
               in that case).
    update_representation : if True (default), (re)creates the "sticks"
               and "nonbonded" representations so the new atom actually
               becomes visible immediately -- "sticks" (not "lines") is
               what the Builder actually draws with, specifically because
               it's bond-order-aware (single/double/triple render as
               1/2/3 parallel cylinders; "lines" only ever draws one
               plain segment per bond regardless of order). Same
               reasoning as recompute_bonds for batching -- skip and do
               it yourself once at the end if adding many atoms in a row.

    [EN] DESIGN CHANGE: distance-based bond auto-detection (vismol_
    object.find_bonded_and_nonbonded_atoms(), the same grid/covalent-
    radius heuristic used when loading a real file) used to run here on
    every call. Turned OFF entirely for the Builder (confirmed, by
    request, after live testing kept surfacing correctness issues that
    all traced back to it -- see _reapply_manual_bonds()'s own docstring
    for the two concrete bugs this caused: a manually-created bond
    silently vanishing when this ran again later for an unrelated atom,
    and, after fixing that, visibly DUPLICATED lines, because that
    native detection's own raw output isn't deduplicated by itself). The
    Builder is a deliberate, atom-by-atom sketchpad -- every bond should
    come from an explicit action (dragging, the 'b' key, or the auto-
    hydrogenation convenience, which now also bonds explicitly -- see
    click_mode.py), never from two atoms merely ending up close to each
    other. Bonds now live ENTIRELY in vismol_object.manual_bonds, rebuilt
    by _reapply_manual_bonds() (bond-object/topology/molecule-grouping
    bookkeeping only -- no more distance queries at all).

    Returns the new Atom object.
    """
    vm_session = vismol_object.vm_session

    if name is None:
        name = symbol

    # --- Chain / Residue (create lazily, same pattern as the file loaders) ---
    if chain_id not in vismol_object.chains:
        vismol_object.chains[chain_id] = Chain ( vismol_object, name = chain_id )
    chain = vismol_object.chains[chain_id]

    if resi not in chain.residues:
        chain.residues[resi] = Residue ( vismol_object, name = resn, index = resi, chain = chain )
    residue = chain.residues[resi]

    atom_id = len ( vismol_object.atoms )   # next free index, sequential

    atom = Atom (
        vismol_object = vismol_object,
        name          = name,
        index         = atom_id + 1,   # index e 1-based por convencao (ver comentario em Atom.__init__)
        residue       = residue,
        chain         = chain,
        symbol        = symbol,
        atom_id       = atom_id,
    )
    atom.pos = np.array ( [ x, y, z ], dtype = np.float32 )
    atom.unique_id = vm_session.atom_id_counter
    atom._generate_atom_unique_color_id ( )
    vm_session.atom_id_counter += 1
    # [EN] REQUIRED for picking/selection to find this atom: the picking
    # system decodes the RGB colour under the mouse click back into
    # pickedID (== atom.unique_id) and looks it up in
    # vm_session.atom_dic_id[pickedID] to get the actual Atom object (see
    # vismol_glcore.py's _pick()/_selection_box_pick()). Every atom
    # normally gets registered here by eSession._add_vismol_object() (for
    # atoms present when an object is FIRST added to the session) -- but
    # register_builder_object() (empty_object.py) intentionally skips
    # that whole method (see its own docstring for why), and an object
    # is only ever added to the session ONCE, while atoms are added one
    # at a time afterwards -- so each new atom needs this registration
    # done here individually instead. Missing this was a real, confirmed
    # bug (atoms added via add_atom() were not selectable at all).
    vm_session.atom_dic_id[atom.unique_id] = atom

    residue.atoms[atom_id]        = atom
    vismol_object.atoms[atom_id]  = atom

    # --- grows vismol_object.frames by +1 atom, for ALL frames already
    # existing (the new atom enters with the SAME position (x,y,z) in each
    # frame -- coerente pra um objeto builder, que tipicamente tem so 1
    # frame mesmo; se este objeto ja tiver uma trajetoria de verdade, isso
    # still works, except the new atom stays "still" across all
    # frames ate que algo mais sofisticado seja implementado). ---
    n_frames = vismol_object.frames.shape[0]
    new_frames = np.zeros ( (n_frames, atom_id + 1, 3), dtype = np.float32 )
    if atom_id > 0:
        new_frames[:, :atom_id, :] = vismol_object.frames
    new_frames[:, atom_id, :] = [ x, y, z ]
    vismol_object.frames = new_frames

    vismol_object.mass_center = np.mean ( vismol_object.frames[0], axis = 0 )

    # recalcula cores/vdw/tamanhos de ponto pra TODOS os atomos -- simples
    # is correct; the colors_id_start parameter is not even used inside the
    # method (checked by reading the source), so the exact value passed
    # here does not matter.
    vismol_object._generate_color_vectors ( vm_session.atom_id_counter )

    # [EN] Register any explicitly-requested bond(s) for this new atom --
    # the ONLY way it ends up bonded to anything now that distance-based
    # auto-detection is off (see this function's own docstring above).
    if bonded_to is not None:
        if not hasattr ( vismol_object, "manual_bonds" ) or vismol_object.manual_bonds is None:
            vismol_object.manual_bonds = set ( )
        if not hasattr ( vismol_object, "manual_bond_orders" ) or vismol_object.manual_bond_orders is None:
            vismol_object.manual_bond_orders = { }
        if not hasattr ( vismol_object, "manual_aromatic_bonds" ) or vismol_object.manual_aromatic_bonds is None:
            vismol_object.manual_aromatic_bonds = set ( )
        targets = bonded_to if isinstance ( bonded_to, ( list, tuple, set ) ) else [ bonded_to ]
        for other_id in targets:
            pair = ( min ( atom_id, other_id ), max ( atom_id, other_id ) )
            vismol_object.manual_bonds.add ( pair )
            vismol_object.manual_bond_orders[pair] = int ( bond_order )
            if aromatic:
                vismol_object.manual_aromatic_bonds.add ( pair )
            else:
                vismol_object.manual_aromatic_bonds.discard ( pair )

    if recompute_bonds:
        # [EN] No more find_bonded_and_nonbonded_atoms() call here (that
        # was the native, distance-based detection -- see this function's
        # docstring for why it was removed). _reapply_manual_bonds() alone
        # is enough: it rebuilds self.bonds/topology/molecule-grouping
        # purely from vismol_object.manual_bonds (which now includes the
        # bonded_to pair(s) just registered above, if any), no distance
        # queries involved at all.
        vismol_object.cov_radii_array = None
        vismol_object.electronegativity_array = None
        vismol_object.index_bonds = None
        vismol_object.bonds = None
        vismol_object.non_bonded_atoms = None
        _reapply_manual_bonds ( vismol_object )

    if update_representation:
        vismol_object.create_representation ( rep_type = "sticks" )
        vismol_object.create_representation ( rep_type = "stick_spheres" )
        _activate_new_sphere_representation ( vismol_object, "stick_spheres" )
        vismol_object.create_representation ( rep_type = "nonbonded" )
        # [EN] REQUIRED for the new atom to be selectable/pickable via
        # mouse click, in addition to the atom_dic_id registration above:
        # vm_glcore.render() only ever calls build_core_representations()
        # (which is what creates "picking_dots", the representation the
        # click-picking system actually reads pixel colours from) ONCE --
        # guarded by "if core_representations['picking_dots'] is None".
        # After that first build, adding more atoms does NOT automatically
        # grow/rebuild it (was_rep_coord_modified only tells the renderer
        # to re-upload COORDINATES for the atom count it already knows
        # about, not that the atom count itself changed). Resetting both
        # core_representations entries to None here forces
        # build_core_representations() to run again, fresh, the next time
        # render() fires -- picking up the current, correct atom count.
        vismol_object.core_representations["picking_dots"] = None
        vismol_object.core_representations["picking_text"] = None
        # [EN] Any OTHER representation this object might carry (e.g.
        # 'cartoon', inherited by an "Edit in Builder" clone of a normal
        # system) needs to either be resized or deactivated here too, or
        # its own stale self.indexes crashes the next picking pass -- see
        # _refresh_bond_dependent_representations()'s own (extended)
        # docstring for the full story.
        _refresh_bond_dependent_representations ( vismol_object )
        if getattr ( vm_session, "vm_glcore", None ) is not None:
            vm_session.vm_glcore.queue_draw ( )

    return atom


def set_atom_element ( vismol_object, atom_id, symbol, name = None,
                        recompute_bonds = True, update_representation = True ):
    """ [EN] Changes the chemical ELEMENT of an already-existing atom, IN

    PLACE -- same atom_id, same position, same bonds-list slot. Used by
    the Builder's "add" tool when a plain click lands ON an atom that
    already exists (see click_mode.handle_click_to_place_atom()):
    instead of stacking a new, essentially-overlapping atom on top of
    it, the existing atom is turned into the newly selected element.

    Unlike add_atom()/remove_atom(), this does NOT touch vismol_object.
    frames (position is unchanged) or vm_session.atom_dic_id / atom.
    unique_id (identity for picking is unchanged -- it's still "the same
    atom", just a different element now) -- only the element-derived
    per-atom attributes are recomputed, exactly the same way Atom.
    __init__ computes them the first time (self.symbol/self.name is set
    first, THEN each _init_*() is called, since every one of them reads
    self.symbol/self.name at call time -- see atom.py):
      - atom.color             (_init_color)
      - atom.vdw_rad            (_init_vdw_rad)
      - atom.cov_rad            (_init_cov_rad)   <- affects bond detection
      - atom.ball_rad           (_init_ball_rad)
      - atom.electronegativity  (_init_electronegativity)

    recompute_bonds : if True (default), rebuilds bonds/topology from
               vismol_object.manual_bonds (via _reapply_manual_bonds())
               after the change. [EN] DESIGN CHANGE: this parameter used
               to matter a lot more -- changing an atom's element used to
               trigger a fresh DISTANCE-based find_bonded_and_nonbonded_
               atoms() call, since cov_rad (which that detection uses)
               changes with the element (e.g. C -> O). Distance-based
               auto-detection is now turned OFF entirely for the Builder
               (see add_atom()'s docstring for why) -- bonds live purely
               in vismol_object.manual_bonds, which changing an element
               never touches. Kept as a parameter mostly for interface
               symmetry with add_atom()/remove_atom(); the rebuild it
               triggers now is just topology/molecule-grouping
               bookkeeping (harmless, cheap), not a real bond
               recalculation.
    update_representation : if True (default), rebuilds "sticks"/
               "nonbonded" and forces a re-build of the picking
               representations, same reasoning as add_atom()/
               remove_atom() (the new colour otherwise wouldn't show up
               until something else forced a full rebuild).

    Returns the (mutated) Atom object.
    """
    if atom_id not in vismol_object.atoms:
        raise ValueError ( "set_atom_element: atom_id {} does not exist in this object.".format ( atom_id ) )

    atom = vismol_object.atoms[atom_id]

    if name is None:
        name = symbol

    atom.symbol = symbol
    atom.name   = name

    # Recompute every element-derived attribute -- same _init_*() methods
    # Atom.__init__ itself calls, now reading the just-updated symbol/name.
    atom.color             = atom._init_color ( )
    atom.vdw_rad           = atom._init_vdw_rad ( )
    atom.cov_rad           = atom._init_cov_rad ( )
    atom.ball_rad          = atom._init_ball_rad ( )
    atom.electronegativity = atom._init_electronegativity ( )

    vm_session = vismol_object.vm_session

    # same call add_atom()/remove_atom() make -- colors_id_start isn't
    # actually used inside the method (confirmed reading the source), so
    # the exact value passed doesn't matter, only that ALL atoms' colour
    # vectors get regenerated (this atom's new colour included).
    vismol_object._generate_color_vectors ( vm_session.atom_id_counter )

    if recompute_bonds:
        _reapply_manual_bonds ( vismol_object )

    if update_representation:
        vismol_object.create_representation ( rep_type = "sticks" )
        vismol_object.create_representation ( rep_type = "stick_spheres" )
        _activate_new_sphere_representation ( vismol_object, "stick_spheres" )
        vismol_object.create_representation ( rep_type = "nonbonded" )
        vismol_object.core_representations["picking_dots"] = None
        vismol_object.core_representations["picking_text"] = None
        # [EN] Any OTHER representation this object might carry (e.g.
        # 'cartoon', inherited by an "Edit in Builder" clone of a normal
        # system) needs to either be resized or deactivated here too, or
        # its own stale self.indexes crashes the next picking pass -- see
        # _refresh_bond_dependent_representations()'s own (extended)
        # docstring for the full story.
        _refresh_bond_dependent_representations ( vismol_object )
        if getattr ( vm_session, "vm_glcore", None ) is not None:
            vm_session.vm_glcore.queue_draw ( )

    return atom


def move_atom ( vismol_object, atom_id, x, y, z, update_representation = True ):
    """ [EN] Repositions an EXISTING atom to (x, y, z) -- element, bonds,
    atom_id, everything else about it stays untouched. Used by the
    Builder's click-and-drag-to-create-a-bonded-atom interaction (see
    click_mode.py's start_bond_drag()/update_bond_drag()/
    finish_bond_drag()) to move the still-being-dragged atom on EVERY
    mouse_motion event, without paying for a full bond-detection
    recompute on every pixel of mouse movement -- that only happens
    ONCE, when the drag ends (finish_bond_drag()).

    Deliberately does NOT touch bonds/topology/cov_radii_array/etc, and
    does NOT call vismol_object.find_bonded_and_nonbonded_atoms() --
    unlike add_atom()/remove_atom(), a simple reposition can't change
    the ATOM COUNT (the thing that method's reset-then-recompute dance
    exists to handle safely), so none of that is needed here; it would
    just make dragging noticeably less smooth for no benefit, run on
    every single mouse pixel of movement.

    update_representation (default True): instead of the FULL
    create_representation() rebuild add_atom()/remove_atom() do (needed
    there because the atom COUNT changed), this uses the SAME cheap
    "just re-upload coordinates, not the whole VBO" mechanism already
    used for trajectory/MD-frame playback -- see forward_frame() /
    reverse_frame() in vismol_session.py -- by setting
    vm_glcore.updated_coords = True (consumed once per frame inside
    render(), which marks every ACTIVE representation's coordinates
    dirty for re-upload). Safe to call every motion event during a drag;
    finish_bond_drag() turns updated_coords back off once the drag ends.

    Returns the (repositioned) Atom object.
    """
    if atom_id not in vismol_object.atoms:
        raise ValueError ( "move_atom: atom_id {} does not exist in this object.".format ( atom_id ) )

    vismol_object.frames[:, atom_id, :] = [ x, y, z ]
    vismol_object.mass_center = np.mean ( vismol_object.frames[0], axis = 0 )

    if update_representation:
        vm_session = vismol_object.vm_session

        # [EN] Mark THIS object's representations dirty DIRECTLY (bug
        # fixed after the user reported the dragged atom only appearing
        # at its new position after releasing the mouse button, instead
        # of following the cursor live): the earlier version relied
        # solely on vm_glcore.updated_coords (the same flag trajectory/MD
        # playback uses -- see forward_frame()/reverse_frame() in
        # vismol_session.py), which is only turned into per-representation
        # rep.was_rep_coord_modified = True INSIDE render()'s own
        # "if self.updated_coords" block, on whatever vm_objects_dic
        # currently holds -- one extra layer of indirection between
        # "coordinates changed" and "this rep will actually re-upload
        # them" that wasn't reliably keeping up every single motion
        # event during a live drag. Setting the flag directly on this
        # object's own representations removes that indirection
        # entirely -- correctness no longer depends on render()'s global
        # per-frame sweep happening to run with the flag still set.
        for rep in vismol_object.representations.values ( ):
            if rep is not None and rep.active:
                rep.was_rep_coord_modified = True
                rep.was_sel_coord_modified = True
        picking_dots = vismol_object.core_representations.get ( "picking_dots" )
        if picking_dots is not None:
            picking_dots.was_rep_coord_modified = True

        if getattr ( vm_session, "vm_glcore", None ) is not None:
            # Kept in addition to the direct marking above (belt-and-
            # suspenders, matches the existing global convention) --
            # harmless if redundant, and covers any OTHER representation
            # that might read this flag in a way not accounted for above.
            vm_session.vm_glcore.updated_coords = True
            vm_session.vm_glcore.queue_draw ( )

    return vismol_object.atoms[atom_id]


def remove_atom ( vismol_object, atom_id ):
    """ [EN] Fourth building block of the Builder (delete-atom tool, 'd'
    key). Removes a single atom from vismol_object by its atom_id.

    Non-trivial because atom_id doubles as a dense array INDEX into
    vismol_object.frames (shape (n_frames, n_atoms, 3)) -- removing atom
    K means deleting frames[:, K, :] AND renumbering every atom with
    atom_id > K down by one, everywhere that number is used as a dict
    key (vismol_object.atoms, residue.atoms) or stored on the Atom object
    itself (atom.atom_id). Getting this renumbering wrong would silently
    corrupt frames/atoms alignment for every atom after the removed one.

    Also cleans up the two places keyed by identity/unique_id rather
    than atom_id (so they don't need renumbering, just removal of the
    deleted atom's own entry): vm_session.atom_dic_id (picking -- a
    stale entry pointing at a removed Atom would be a dangling
    reference) and every active selection's selected_atoms set (so a
    deleted atom can't stay "selected", which would be a real problem
    for add_bond() below if it silently included a since-deleted atom).

    Bonds/topology: rebuilt from vismol_object.manual_bonds (already
    remapped to the new atom numbering just above) via
    _reapply_manual_bonds() -- no distance queries involved. Distance-
    based auto-detection is turned OFF entirely for the Builder now (see
    add_atom()'s docstring for why); a manually-added bond DOES survive
    remove_atom(), as long as it doesn't involve the atom actually being
    removed (which has nothing sensible to remap to, so that specific
    pair is dropped -- see the id_map remapping above). """
    if atom_id not in vismol_object.atoms:
        raise ValueError ( "remove_atom: atom_id {} does not exist in this object.".format ( atom_id ) )

    vm_session = vismol_object.vm_session
    removed_atom = vismol_object.atoms[atom_id]

    # tira do dicionario global de picking (senao fica uma referencia
    # dangling pointing to an Atom that no longer exists in the object)
    if removed_atom.unique_id in vm_session.atom_dic_id:
        del vm_session.atom_dic_id[removed_atom.unique_id]

    # removes it from any active selection (otherwise a subsequent add_bond()
    # could mistakenly include an atom that was just deleted)
    for sel in vm_session.selections.values ( ):
        sel.selected_atoms.discard ( removed_atom )

    # tira do dict de atomos do residuo
    if removed_atom.residue is not None:
        removed_atom.residue.atoms.pop ( atom_id, None )

    del vismol_object.atoms[atom_id]

    # removes the corresponding row from ALL frames (axis=1 = atoms axis)
    vismol_object.frames = np.delete ( vismol_object.frames, atom_id, axis = 1 )

    # renumera atom_id de todo mundo com atom_id > removido, subtraindo 1 --
    # em vismol_object.atoms, em residue.atoms, e no proprio atom.atom_id
    old_atoms = vismol_object.atoms
    new_atoms = {}
    id_map = {}   # old_id -> new_id, para tambem remapear manual_bonds abaixo
    for old_id in sorted ( old_atoms.keys ( ) ):
        atom = old_atoms[old_id]
        new_id = old_id - 1 if old_id > atom_id else old_id
        id_map[old_id] = new_id
        if new_id != old_id:
            atom.atom_id = new_id
            if atom.residue is not None:
                atom.residue.atoms.pop ( old_id, None )
                atom.residue.atoms[new_id] = atom
        new_atoms[new_id] = atom
    vismol_object.atoms = new_atoms

    # [EN] REQUIRED, not optional: atom_id is a dense array index, so
    # removing one RENUMBERS every atom after it (see id_map above) --
    # manual_bonds stores atom_id PAIRS, so without remapping them the
    # SAME way, a bond like (3, 5) would silently start pointing at
    # whatever atoms happen to occupy ids 3/5 AFTER the shift, not the
    # atoms it was actually meant for. A pair involving the atom being
    # removed itself has nothing sensible to remap to, so it's dropped
    # instead (that connection genuinely no longer exists).
    if hasattr ( vismol_object, "manual_bonds" ) and vismol_object.manual_bonds:
        remapped_manual_bonds = set ( )
        for a, b in vismol_object.manual_bonds:
            if a == atom_id or b == atom_id:
                continue
            new_pair = ( min ( id_map[a], id_map[b] ), max ( id_map[a], id_map[b] ) )
            remapped_manual_bonds.add ( new_pair )
        vismol_object.manual_bonds = remapped_manual_bonds

    # [EN] Same remapping, same reasoning, for manually-set bond ORDERS
    # (click_mode.apply_selected_bond_order()) -- these are ALSO keyed by atom_id
    # pairs, so they'd silently end up describing the wrong bond after a
    # removal-triggered renumbering otherwise.
    if hasattr ( vismol_object, "manual_bond_orders" ) and vismol_object.manual_bond_orders:
        remapped_manual_bond_orders = { }
        for ( a, b ), order in vismol_object.manual_bond_orders.items ( ):
            if a == atom_id or b == atom_id:
                continue
            new_pair = ( min ( id_map[a], id_map[b] ), max ( id_map[a], id_map[b] ) )
            remapped_manual_bond_orders[new_pair] = order
        vismol_object.manual_bond_orders = remapped_manual_bond_orders

    # [EN] Same remapping, same reasoning, for the aromatic-flag set --
    # see add_bond()'s own docstring for why aromaticity is tracked
    # separately from bond_order (mirrors pDynamo3's own Bond model).
    if hasattr ( vismol_object, "manual_aromatic_bonds" ) and vismol_object.manual_aromatic_bonds:
        remapped_manual_aromatic_bonds = set ( )
        for a, b in vismol_object.manual_aromatic_bonds:
            if a == atom_id or b == atom_id:
                continue
            new_pair = ( min ( id_map[a], id_map[b] ), max ( id_map[a], id_map[b] ) )
            remapped_manual_aromatic_bonds.add ( new_pair )
        vismol_object.manual_aromatic_bonds = remapped_manual_aromatic_bonds

    vismol_object.mass_center = ( np.mean ( vismol_object.frames[0], axis = 0 )
                                  if vismol_object.frames.shape[1] > 0
                                  else np.zeros ( 3, dtype = np.float32 ) )
    vismol_object._generate_color_vectors ( vm_session.atom_id_counter )

    # bonds/topology now come entirely from vismol_object.manual_bonds
    # (ja remapeado acima) -- ver nota no docstring sobre deteccao por
    # distancia ter sido desligada.
    vismol_object.cov_radii_array = None
    vismol_object.electronegativity_array = None
    vismol_object.index_bonds = None
    vismol_object.bonds = None
    vismol_object.non_bonded_atoms = None
    if len ( vismol_object.atoms ) > 0:
        _reapply_manual_bonds ( vismol_object )
    else:
        vismol_object.index_bonds = []
        vismol_object.bonds = {}
        vismol_object.non_bonded_atoms = []

    vismol_object.create_representation ( rep_type = "sticks" )
    vismol_object.create_representation ( rep_type = "stick_spheres" )
    _activate_new_sphere_representation ( vismol_object, "stick_spheres" )
    vismol_object.create_representation ( rep_type = "nonbonded" )
    # ver nota em add_atom() -- forca build_core_representations() a
    # reconstruir do zero na proxima renderizacao, com a contagem de
    # atomos (agora menor) correta.
    vismol_object.core_representations["picking_dots"] = None
    vismol_object.core_representations["picking_text"] = None

    # [EN] Atom count just changed -- see _refresh_bond_dependent_
    # representations()'s own (extended) docstring: any OTHER
    # representation this object carries (e.g. 'cartoon', inherited by
    # an "Edit in Builder" clone of a normal system) needs to either be
    # resized or deactivated here too, or its own stale self.indexes
    # crashes the next picking pass.
    _refresh_bond_dependent_representations ( vismol_object )

    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.queue_draw ( )


def _reapply_manual_bonds ( vismol_object ):
    """ [EN] Re-merges every explicitly-added bond (this module's own
    add_bond()) recorded in vismol_object.manual_bonds back into the
    CURRENT index_bonds/bonds/topology -- i.e. restores manual bonds
    that a from-scratch find_bonded_and_nonbonded_atoms() recompute may
    have just silently dropped (that recompute is purely distance-
    based; a manual bond exists specifically because the two atoms are
    NOT within auto-detection range, so it can never reconstruct one on
    its own).

    BUG FIXED #1 (caught by live testing, not just reading): dragging
    out a SECOND bonded atom (click_mode.py's start_bond_drag(), which
    calls add_atom() to create the new atom) was silently erasing the
    bond created by a PREVIOUS drag -- add_atom()'s own find_bonded_and_
    nonbonded_atoms() call has no memory of bonds add_bond() added
    manually in an earlier, unrelated step. This was exactly the "KNOWN
    LIMITATION" add_bond() already documented in its own docstring, just
    not yet acted on anywhere. Fixed generically (not special-cased just
    for the drag feature) by calling this helper from every call site
    that does a from-scratch recompute: add_atom(), remove_atom(),
    set_atom_element(), and click_mode.finish_bond_drag() -- so ANY
    manually-added bond survives ANY of them, not just the one pair the
    bug happened to be reported for.

    BUG FIXED #2 (caught by live testing right after #1 -- dragging a
    SECOND bonded atom out was then showing duplicated/spurious lines
    instead of erasing anything): the native, purely distance-based
    find_bonded_and_nonbonded_atoms() (vismol_object.py's own
    cdist.get_atomic_bonds_from_grid()) stores each detected bond
    potentially MORE THAN ONCE in self.index_bonds -- both [i,j] and
    [j,i] for the SAME physical bond, and evidently not even always
    symmetrically (confirmed live: 3 atoms, all 3 pairs within bonding
    distance, ended up with index_bonds holding 10 flat integers -- 5
    raw (i,j) entries -- for what should be at most 3 UNIQUE bonds).
    That native call's own _bonds_from_pair_of_indexes_list() (which
    find_bonded_and_nonbonded_atoms() triggers internally, BEFORE this
    function ever runs) does NOT deduplicate: it creates a fresh Bond()
    object for EVERY raw occurrence and appends it to BOTH atoms' own
    .bonds lists every time (self.bonds, the OBJECT-level dict, happens
    to self-correct since it's keyed by a normalized (min,max) tuple and
    a later duplicate simply overwrites the same key -- but self.
    index_bonds itself, and each ATOM's .bonds list, keep every
    duplicate, and THAT's what the "lines" representation actually
    renders from). The first version of this function only rebuilt
    (deduplicating as a side effect) when a genuinely NEW manual bond
    needed adding -- so any call where nothing NEW needed adding (e.g.
    the pair was already auto-detected) skipped the rebuild entirely and
    let the native call's own un-deduplicated data pass straight
    through untouched, duplicates and all. Fixed by ALWAYS rebuilding
    from the deduplicated pair set below whenever there's at least one
    bond at all (not only when something NEW was added) -- this function
    now doubles as "restore manual bonds" AND "sanitize whatever the
    native call just produced", every single time it's called.

    Returns True if there was at least one bond ((re)built normally),
    False if the object has zero bonds at all (an isolated atom / all
    atoms unbonded) -- in that case self.bonds/self.non_bonded_atoms/
    self.index_bonds are still properly initialised to valid EMPTY
    structures ({}, all-atoms-non-bonded, and an empty array
    respectively), not left as None -- downstream representation code
    expects to be able to iterate them regardless of whether there are
    any bonds. Silently forgets (removes from vismol_object.manual_bonds)
    any recorded pair where an atom no longer exists (e.g. removed via
    remove_atom() since the bond was added) -- nothing sensible to
    restore there. """
    if not hasattr ( vismol_object, "manual_bonds" ) or vismol_object.manual_bonds is None:
        vismol_object.manual_bonds = set ( )
    manual_bonds = vismol_object.manual_bonds

    if vismol_object.index_bonds is not None and len ( vismol_object.index_bonds ) > 0:
        flat = np.asarray ( vismol_object.index_bonds ).reshape ( -1, 2 )
        existing_pairs = set ( tuple ( sorted ( p ) ) for p in flat.tolist ( ) )
    else:
        existing_pairs = set ( )

    for pair in list ( manual_bonds ):
        a, b = pair
        if a not in vismol_object.atoms or b not in vismol_object.atoms:
            manual_bonds.discard ( pair )

    merged_pairs = existing_pairs | manual_bonds   # union -- already deduplicated, since both sides are sets of normalized (min,max) tuples

    # [EN] BUG FIXED: this used to return early here (before rebuilding
    # anything) whenever there were zero bonds at all -- leaving self.
    # bonds/self.non_bonded_atoms/self.index_bonds sitting at None
    # (whatever the caller had just reset them to), instead of the valid,
    # empty-but-well-formed structures create_representation() actually
    # expects to iterate. Harmless for a lone atom added by itself, but
    # confirmed to break the VERY NEXT add_atom()/remove_atom() call on
    # the same object once distance-based auto-detection (which used to
    # paper over this by always producing a real, if empty, result) was
    # turned off. Fixed by always running the full rebuild below,
    # including the zero-bonds case -- _bonds_from_pair_of_indexes_list()
    # and friends handle an empty index_bonds array just fine (same
    # "zero bonds" state a freshly-created single atom always had
    # anyway, even back when auto-detection was still on).

    new_index_bonds = []
    for i, j in sorted ( merged_pairs ):
        new_index_bonds.append ( i )
        new_index_bonds.append ( j )

    vismol_object.bonds = None
    vismol_object.non_bonded_atoms = None
    vismol_object.index_bonds = np.array ( new_index_bonds, dtype = np.int64 )

    # [EN] Aligns with _bonds_from_pair_of_indexes_list()'s own
    # "Convention B": external_orders[k] must be the order of the k-th
    # bond that SURVIVES that method's own exclude_list filter (default
    # [['H','H']]) -- not the k-th raw pair.
    #
    # [EN] BUG FIX: this used to build external_orders with
    # manual_bond_orders.get((i,j), 1) -- i.e. DEFAULT TO A PLAIN SINGLE
    # BOND for every pair the user hadn't explicitly cycled/set an order
    # for. Combined with _bonds_from_pair_of_indexes_list() ALSO ignoring
    # external_orders entirely (separately fixed there -- see that
    # method's own updated comment), the net effect used to be "always
    # auto-perceive, manual overrides silently discarded". Fixing ONLY
    # the other method without fixing this default would have flipped
    # the bug the other way: "always plain single bonds, auto-perception
    # (aromatic/conjugated rings etc. -- see perceive_bond_order_for_pairs)
    # never runs for anything drawn in the Builder". Fixed by computing a
    # BASELINE via the same automatic perception used everywhere else
    # (perceive_bond_order_for_pairs -- exact maximum-matching heuristic,
    # handles conjugated rings correctly, see that method's own docstring
    # for the full history), then overriding ONLY the specific pairs the
    # user explicitly set in manual_bond_orders (via 'bond order=N' or
    # Ctrl+click apply_selected_bond_order()) on top of that baseline -- matching
    # apply_selected_bond_order()'s own stated intent ("persist the new order...
    # REQUIRED... without persisting it... changing a bond's order here
    # would get silently overwritten back to the default").
    manual_bond_orders = getattr ( vismol_object, "manual_bond_orders", None ) or { }
    exclude_list = [ [ 'H', 'H' ] ]

    surviving_pairs = [ ]
    for i, j in sorted ( merged_pairs ):
        symbol_i = vismol_object.atoms[i].symbol
        symbol_j = vismol_object.atoms[j].symbol
        is_excluded = any ( symbol_i in pair and symbol_j in pair for pair in exclude_list )
        if is_excluded:
            continue
        surviving_pairs.append ( ( i, j ) )

    if surviving_pairs:
        flat_pairs = [ idx for pair in surviving_pairs for idx in pair ]
        baseline_orders = vismol_object.perceive_bond_order_for_pairs ( flat_pairs ).tolist ( )
    else:
        baseline_orders = [ ]

    external_orders = [ ]
    for k, ( i, j ) in enumerate ( surviving_pairs ):
        external_orders.append ( manual_bond_orders.get ( ( i, j ), baseline_orders[k] ) )

    vismol_object._bonds_from_pair_of_indexes_list ( external_orders = external_orders )
    vismol_object._get_non_bonded_from_bonded_list ( )
    vismol_object._generate_topology_from_index_bonds ( )
    vismol_object.define_molecules ( )
    vismol_object.define_Calpha_backbone ( )
    return bool ( merged_pairs )


def add_bond ( vismol_object, atom_id_a, atom_id_b, bond_order = 1, aromatic = False ):
    """ [EN] Fifth building block of the Builder ('b' key -- add a bond
    between exactly two currently-selected atoms). Unlike the automatic,
    distance-based bond detection find_bonded_and_nonbonded_atoms() does
    (used internally by add_atom()/remove_atom() above), this creates an
    EXPLICIT bond regardless of the distance between the two atoms --
    necessary for cases automatic detection would never produce on its
    own (e.g. deliberately closing a ring across a gap the covalent-
    radius/distance heuristic wouldn't recognise as bonded).

    Implementation: the pair is remembered PERMANENTLY in
    vismol_object.manual_bonds (created lazily), regardless of whether
    it happens to already be auto-detected at the moment this is called
    -- then _reapply_manual_bonds() (see its own docstring) merges
    manual_bonds into the current index_bonds and does the same
    downstream bookkeeping find_bonded_and_nonbonded_atoms() itself does
    (_bonds_from_pair_of_indexes_list, _get_non_bonded_from_bonded_list,
    _generate_topology_from_index_bonds, define_molecules,
    define_Calpha_backbone). Recording it unconditionally (not only when
    it turns out NOT to already be auto-detected) is what makes it
    survive add_atom()/remove_atom()/set_atom_element() recomputing
    everything from scratch LATER, even if at the moment of THIS call it
    happened to already be in bonding distance (fixes the KNOWN
    LIMITATION this function used to have, and still describes below for
    context -- confirmed via live testing to have actually been hit by
    the Builder's click-and-drag-to-create-a-bonded-atom feature: a
    SECOND drag's add_atom() call was silently erasing the FIRST drag's
    bond).

    [EN] BUG FIX: bond_order used to be accepted as a parameter here but
    was NEVER ACTUALLY USED anywhere in this function's body -- calling
    add_bond(..., bond_order=2) silently produced a plain single bond
    every time, exactly the same class of dead-parameter bug already
    found and fixed once in click_mode.apply_selected_bond_order() (see that
    function's own docstring: "_bonds_from_pair_of_indexes_list()'s
    external_orders parameter already existed but its actual assignment
    ... was commented out"). Fixed the same way apply_selected_bond_order()
    already does it: persist the requested order in
    vismol_object.manual_bond_orders (keyed by the normalized (min,max)
    pair), which _reapply_manual_bonds() below feeds into
    _bonds_from_pair_of_indexes_list() as external_orders -- required,
    not optional, since bonds get rebuilt FROM SCRATCH (fresh Bond()
    objects) on every structural edit, so anything not persisted there
    is silently lost on the very next add_atom()/remove_atom()/add_bond()
    call.

    [EN] aromatic: mirrors bond_order's own persistence, but in a
    SEPARATE set (manual_aromatic_bonds) rather than folded into
    bond_order's own 1/2/3 range -- matches pDynamo3's own Bond model
    (pMolecule/Bond.py), whose BondType enum has only Single/Double/
    Triple ("Aromatic bond types are not defined. Instead bonds
    themselves are flagged if they form part of an aromatic system."). """
    if atom_id_a == atom_id_b:
        raise ValueError ( "add_bond: cannot bond an atom to itself." )
    if atom_id_a not in vismol_object.atoms or atom_id_b not in vismol_object.atoms:
        raise ValueError ( "add_bond: atom_id_a={} or atom_id_b={} does not exist in this object.".format (
                            atom_id_a, atom_id_b ) )
    if bond_order not in ( 1, 2, 3 ):
        raise ValueError ( "add_bond: bond_order must be 1 (single), 2 (double) or 3 (triple), got {!r}.".format (
                            bond_order ) )

    pair = ( min ( atom_id_a, atom_id_b ), max ( atom_id_a, atom_id_b ) )

    if not hasattr ( vismol_object, "manual_bonds" ) or vismol_object.manual_bonds is None:
        vismol_object.manual_bonds = set ( )
    vismol_object.manual_bonds.add ( pair )

    if not hasattr ( vismol_object, "manual_bond_orders" ) or vismol_object.manual_bond_orders is None:
        vismol_object.manual_bond_orders = { }
    vismol_object.manual_bond_orders[pair] = int ( bond_order )

    if not hasattr ( vismol_object, "manual_aromatic_bonds" ) or vismol_object.manual_aromatic_bonds is None:
        vismol_object.manual_aromatic_bonds = set ( )
    if aromatic:
        vismol_object.manual_aromatic_bonds.add ( pair )
    else:
        vismol_object.manual_aromatic_bonds.discard ( pair )

    changed = _reapply_manual_bonds ( vismol_object )
    if not changed:
        return False   # object with no bonds at all (neither this new one nor any other) -- nothing to do (rare/defensive case)

    _refresh_bond_dependent_representations ( vismol_object )
    vismol_object.core_representations["picking_dots"] = None
    vismol_object.core_representations["picking_text"] = None

    vm_session = vismol_object.vm_session
    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.queue_draw ( )

    return True


def remove_bond ( vismol_object, atom_id_a, atom_id_b ):
    """ [EN] Sixth building block of the Builder (right-click on a bond,
    Avogadro-style -- see the mouse_released() hook in vismol_glcore.py).
    Removes the single bond between atom_id_a and atom_id_b -- the atoms
    THEMSELVES are untouched, only the connection between them.

    Since bonds are ENTIRELY explicit now (vismol_object.manual_bonds --
    distance-based auto-detection is off, see add_atom()'s docstring),
    this is simply: forget the pair (and any custom order recorded for
    it), then rebuild. Nothing can silently bring it back later the way
    the old, now-removed, distance auto-detector could have (if the two
    atoms still happened to be close together).

    Returns True if a bond was actually removed, False if that pair
    wasn't bonded to begin with (nothing to do -- e.g. double right-click
    on the same spot after the bond is already gone). """
    pair = ( min ( atom_id_a, atom_id_b ), max ( atom_id_a, atom_id_b ) )

    manual_bonds = getattr ( vismol_object, "manual_bonds", None )
    if not manual_bonds or pair not in manual_bonds:
        return False

    manual_bonds.discard ( pair )

    manual_bond_orders = getattr ( vismol_object, "manual_bond_orders", None )
    if manual_bond_orders:
        manual_bond_orders.pop ( pair, None )
    manual_aromatic_bonds = getattr ( vismol_object, "manual_aromatic_bonds", None )
    if manual_aromatic_bonds:
        manual_aromatic_bonds.discard ( pair )

    vismol_object.cov_radii_array = None
    vismol_object.electronegativity_array = None
    vismol_object.index_bonds = None
    vismol_object.bonds = None
    vismol_object.non_bonded_atoms = None
    _reapply_manual_bonds ( vismol_object )

    _refresh_bond_dependent_representations ( vismol_object )
    vismol_object.core_representations["picking_dots"] = None
    vismol_object.core_representations["picking_text"] = None

    vm_session = vismol_object.vm_session
    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.queue_draw ( )

    return True


def _activate_new_sphere_representation ( vismol_object, rep_type ):
    """ [EN] REQUIRED right after vismol_object.create_representation()
    for any Sphere-family rep_type ("spheres"/"stick_spheres"/
    "vdw_spheres"/"picking_spheres") -- found via live testing (a real
    screenshot showing nothing at all where sphere joints should be,
    right after standardising the Builder on ball-and-stick rendering).

    Representations.SpheresRepresentation's OWN _make_gl_representation_
    vao_and_vbos() (unlike Lines/SticksRepresentation's, which upload
    REAL coordinate/colour data immediately) fills its colour/radius/
    per-instance-coordinate GPU buffers with np.zeros(...) PLACEHOLDERS
    at construction time, and only uploads the real data the first time
    draw_representation() sees was_rep_ind_modified (or was_col_modified/
    was_rep_coord_modified) set True. In the main app, VismolSession.
    show_or_hide() is what sets that flag (right after computing which
    atoms should be shown, based on each atom's own .spheres/.
    stick_spheres/... boolean) -- the Builder has no equivalent per-atom
    show/hide toggle (it always wants ALL current atoms shown), so
    without this, a freshly created Sphere-family representation here
    would silently render as zero-radius, zero-colour spheres stacked at
    the origin: technically "active" and drawn every frame, but
    completely invisible. """
    rep = vismol_object.representations.get ( rep_type )
    if rep is not None:
        rep.was_rep_ind_modified = True


def _refresh_bond_dependent_representations ( vismol_object ):
    """ [EN] Recreates every bond-dependent representation that is
    CURRENTLY PRESENT on this object (i.e. was created at some point --
    not necessarily the one the user happens to be looking at right now,
    but any that exist), preserving each one's active/inactive state.

    [EN] BUG FIX (found via live testing): set_bond_order()/unset_bond()
    below (and, before them, add_bond()/remove_bond()) only ever called
    vismol_object.create_representation() for rep_type='lines' and
    'nonbonded' -- NEVER 'sticks', even though 'sticks' is the ONLY
    representation that actually draws double/triple bonds distinctly
    (see representations.py's SticksRepresentation._get_bond_order_per_
    bond() / the geometry shader's u_bond_order_tbo). Confirmed live:
    'unbond' correctly removed the bond from vismol_object.bonds (a
    second 'unbond' on the same pair correctly reported "no bond
    existed"), but the on-screen STICKS view kept showing the old bond --
    because that representation object's OWN self.indexes (its private
    copy of the atom-pair list actually used for the GPU element buffer,
    set once at creation time by define_new_indexes_to_vbo()) was never
    refreshed; only vismol_object.index_bonds itself was updated, and
    nothing told the existing SticksRepresentation object to rebuild from
    it. 'lines' was fine only because it happened to be the one
    representation type these functions did already recreate.

    Only rep_types that are NOT None in vismol_object.representations get
    touched (i.e. only ones actually in use for this object) -- creating
    'sticks' from scratch for an object that never had it would silently
    turn sticks rendering ON for something the user never asked to see in
    that style. metal_dash refreshes as a side effect of recreating
    'lines'/'sticks' (see create_representation()'s own _ensure_metal_
    dash() calls), so it doesn't need its own entry here.

    [EN] EXTENDED (BUG FIX, found by live user testing -- a real crash,
    "IndexError: index 18 is out of bounds for axis 0 with size 18",
    clicking to place an atom on a system opened via "Edit in Builder"):
    ALSO deactivates every OTHER currently-active representation whose
    rep_type is NOT one of the 3 rebuilt above. Root cause: "Edit in
    Builder" (empty_object.begin_editing_existing_system()) clones a
    NORMALLY-loaded system (via p_session.clone_system() ->
    _add_vismol_object_to_easyhybrid_session()), which -- unlike a
    from-scratch Builder object (create_empty_vismol_object(), which
    only ever has the 3 rep_types above) -- can carry a much bigger set
    of representations active by default (e.g. 'cartoon' for a protein).
    Every structural edit that changes the ATOM COUNT (add_atom(),
    remove_atom(), attach_fragment_at_hydrogen()) resizes vismol_object.
    frames, but this function (and those callers) only ever rebuilt the
    3 Builder-relevant rep_types -- any OTHER already-active
    representation kept its OLD, now too-large-or-small self.indexes
    array (each Representation's OWN private copy, set once at creation
    time -- same staleness class already described above for 'sticks'
    pre-fix), and vm_glcore's picking pass (click_mode.
    _read_depth_and_atom_at_pixel(), used by EVERY Builder click, not
    just this one) iterates and draws every ACTIVE representation of
    EVERY object in the scene to read pixel IDs, crashing the instant it
    reached one of these stale ones. Deactivating them (rather than
    trying to correctly rebuild every possible representation type, some
    of which -- 'cartoon' -- need extra state this module has no
    business recomputing) sidesteps the crash entirely and also matches
    this Builder's own already-established preference for 'sticks' over
    'lines' as its one true bond-order-aware view (see the "Builder
    sticks representation" plan entry) -- a stale-but-inactive
    representation is picked up again, correctly, once editing ends
    (finish_editing_existing_system() rebuilds a normal object's usual
    representations from scratch via the SAME channel a plain reload
    would use, not by reactivating these disabled ones as-is). """
    # [EN] EXTENDED AGAIN (user's own explicit follow-up request --
    # "vamos desativar a representacao de linhas para os objetos que
    # estao sendo editados"): 'lines' is rebuilt (so its own self.indexes
    # never goes stale either, in case something re-activates it later)
    # but is now ALWAYS forced INACTIVE here, never preserving whatever
    # its previous active state was -- the Builder's own established
    # preference is 'sticks' (bond-order-aware) over 'lines' (a single
    # plain segment regardless of order, see the "Builder sticks
    # representation" plan entry) for ANYTHING it draws itself, but an
    # "Edit in Builder" clone of a normal system can still have 'lines'
    # active (inherited from the system it was cloned from) alongside
    # 'sticks' -- both drawing the same bonds on top of each other is
    # redundant at best and a Z-fighting risk at worst (see the
    # "aromatic double-bond position error" diagnosis elsewhere in this
    # project's history, a real Z-fighting bug between two near-
    # duplicate objects, not this exact case, but the same visual
    # symptom class). 'sticks'/'nonbonded' keep the existing preserve-
    # active-state behaviour unchanged.
    # [EN] BUG FOUND + FIXED (user's own follow-up request -- "vamos
    # padronizar a representacao de fragmentos e estruturas na forma de
    # esferas e bastoes"): 'stick_spheres' used to be MISSING from
    # rebuilt_types below -- every single caller of this function
    # (add_atom(), attach_fragment_at_hydrogen(), add_structure_at_
    # position(), ...) already does its OWN vismol_object.create_
    # representation(rep_type='stick_spheres') + _activate_new_sphere_
    # representation() call right before reaching here (see e.g. add_
    # atom()'s own tail), correctly leaving it ACTIVE -- but the second
    # loop below ("deactivate anything not in rebuilt_types", added LATER
    # for the Edit-in-Builder stale-representation crash fix, see this
    # function's own docstring) then immediately swept 'stick_spheres' up
    # as if it were a foreign, possibly-stale representation (like
    # 'cartoon') and forced it back OFF -- silently undoing the ball-and-
    # stick standardisation on EVERY Builder mutation, not just fragments/
    # structures. Confirmed live: vismol_object.representations['stick_
    # spheres'].active was False immediately after add_structure_at_
    # position() returned, despite that function's own explicit
    # _activate_new_sphere_representation() call a few lines earlier.
    # Fixed by treating 'stick_spheres' as a 4th REBUILT (not swept)
    # type -- same "preserve whatever active state it already had"
    # handling as 'sticks'/'nonbonded' -- plus re-running _activate_new_
    # sphere_representation() on the freshly-recreated copy, since ANY
    # newly-constructed SpheresRepresentation needs that same GPU-buffer-
    # fill kick (see that function's own docstring), not just the very
    # first one ever created for a given object.
    reps = getattr ( vismol_object, "representations", None ) or { }
    rebuilt_types = ( "lines", "sticks", "nonbonded", "stick_spheres" )
    for rep_type in rebuilt_types:
        rep = reps.get ( rep_type )
        if rep is None:
            continue
        was_active = False if rep_type == "lines" else getattr ( rep, "active", True )
        vismol_object.create_representation ( rep_type = rep_type )
        new_rep = vismol_object.representations.get ( rep_type )
        if new_rep is not None:
            new_rep.active = was_active
            if rep_type == "stick_spheres":
                _activate_new_sphere_representation ( vismol_object, "stick_spheres" )

    for rep_type, rep in reps.items ( ):
        if rep is None or rep_type in rebuilt_types:
            continue
        if getattr ( rep, "active", False ):
            rep.active = False


def _sync_index_bonds_and_order_list_from_bonds_dict ( vismol_object ):
    """ [EN] Small helper shared by set_bond_order()/unset_bond() below:
    regenerates vismol_object.index_bonds and vismol_object.bond_order_list
    DIRECTLY from vismol_object.bonds (the dict, iterated in its current
    order -- Python dicts preserve insertion order, so this stays stable
    across repeated calls). Both arrays end up paired 1:1, same index k in
    both. Does NOT touch self.bonds itself, self.atoms[*].bonds, manual_
    bonds, topology, non_bonded_atoms, or anything else -- purely a
    flatten of whatever self.bonds already is. """
    flat = [ ]
    orders = [ ]
    for ( i, j ), bond in vismol_object.bonds.items ( ):
        flat.append ( i )
        flat.append ( j )
        orders.append ( int ( bond.bond_order ) )
    vismol_object.index_bonds = np.array ( flat, dtype = np.int64 )
    vismol_object.bond_order_list = orders


def set_bond_order ( vismol_object, atom_id_a, atom_id_b, bond_order = 1, aromatic = False ):
    """ [EN] Adds a bond between atom_id_a/atom_id_b if none exists yet,

    or updates its order if one already does -- WITHOUT rebuilding the
    object's entire bond set from scratch, unlike add_bond()/
    _reapply_manual_bonds().

    WHY THIS EXISTS (found via live testing, not just reading): add_bond()
    /_reapply_manual_bonds() are built around the BUILDER's specific
    design, where vismol_object.manual_bonds is meant to be the ONLY
    source of truth for connectivity (see add_atom()'s own docstring:
    "distance-based auto-detection is off... Bonds now live ENTIRELY in
    vismol_object.manual_bonds"). That design is correct FOR THE BUILDER
    (an empty object grown one atom/bond at a time, where every bond really
    did come from an explicit add_bond() call) -- but breaks badly on a
    NORMALLY LOADED structure (e.g. a PDB file read via
    find_bonded_and_nonbonded_atoms()): those bonds live in self.bonds/
    self.index_bonds from distance-based auto-detection and were NEVER
    registered in self.manual_bonds. Calling add_bond() there merges
    self.index_bonds (which, AT THAT MOMENT, still has everything) with
    manual_bonds and rebuilds via _reapply_manual_bonds() -- confirmed
    live to end up DROPPING every other bond in the structure, keeping
    only the just-added one. (remove_bond() is worse: it explicitly resets
    index_bonds = None BEFORE rebuilding, so uses ONLY manual_bonds --
    guaranteed to drop everything not explicitly added via add_bond() in
    the current session.)

    This function sidesteps all of that: it mutates vismol_object.bonds
    (the actual dict everything else already reads from -- get_bond(),
    representations.py's _get_bond_order_per_bond(), etc.) directly, for
    ONLY the one pair involved, then regenerates index_bonds/
    bond_order_list FROM that dict (see
    _sync_index_bonds_and_order_list_from_bonds_dict() above) -- every
    OTHER bond, wherever it originally came from, is left completely
    untouched.

    Also keeps manual_bonds/manual_bond_orders in sync (adds this pair)
    purely so this bond SURVIVES if the Builder's own add_atom()/
    remove_atom()/add_bond() run later on this same object -- without that,
    a bond created this way would vanish the next time any of those does
    its own from-scratch, manual_bonds-only rebuild.

    Returns True if a NEW bond was created, False if an existing bond's
    order was updated (or left the same -- still returns False, since no
    NEW bond was created either way). """
    if atom_id_a == atom_id_b:
        raise ValueError ( "set_bond_order: cannot bond an atom to itself." )
    if atom_id_a not in vismol_object.atoms or atom_id_b not in vismol_object.atoms:
        raise ValueError ( "set_bond_order: atom_id_a={} or atom_id_b={} does not exist in this object.".format (
                            atom_id_a, atom_id_b ) )
    if bond_order not in ( 1, 2, 3 ):
        raise ValueError ( "set_bond_order: bond_order must be 1 (single), 2 (double) or 3 (triple), got {!r}.".format (
                            bond_order ) )

    if vismol_object.bonds is None:
        vismol_object.bonds = { }

    key = ( min ( atom_id_a, atom_id_b ), max ( atom_id_a, atom_id_b ) )
    existing = vismol_object.bonds.get ( key )

    if existing is not None:
        existing.bond_order = int ( bond_order )
        created = False
    else:
        from vismol.model.bond import Bond
        atom_a = vismol_object.atoms[atom_id_a]
        atom_b = vismol_object.atoms[atom_id_b]
        bond = Bond ( atom_i = atom_a, atom_index_i = atom_id_a,
                       atom_j = atom_b, atom_index_j = atom_id_b )
        bond.bond_order = int ( bond_order )
        vismol_object.bonds[key] = bond
        atom_a.bonds.append ( bond )
        atom_b.bonds.append ( bond )
        atom_a.nbonds = len ( atom_a.bonds )
        atom_b.nbonds = len ( atom_b.bonds )
        created = True

    _sync_index_bonds_and_order_list_from_bonds_dict ( vismol_object )

    if not hasattr ( vismol_object, "manual_bonds" ) or vismol_object.manual_bonds is None:
        vismol_object.manual_bonds = set ( )
    vismol_object.manual_bonds.add ( key )
    if not hasattr ( vismol_object, "manual_bond_orders" ) or vismol_object.manual_bond_orders is None:
        vismol_object.manual_bond_orders = { }
    vismol_object.manual_bond_orders[key] = int ( bond_order )

    if not hasattr ( vismol_object, "manual_aromatic_bonds" ) or vismol_object.manual_aromatic_bonds is None:
        vismol_object.manual_aromatic_bonds = set ( )
    if aromatic:
        vismol_object.manual_aromatic_bonds.add ( key )
    else:
        vismol_object.manual_aromatic_bonds.discard ( key )

    if created:
        # Topologia/non-bonded/moleculas so' precisam ser refeitas quando a
        # CONNECTIVITY changes (new bond) -- a simple bond-order change on an
        # already existing bond does not affect any of these (identical graph).
        vismol_object.non_bonded_atoms = None
        vismol_object._get_non_bonded_from_bonded_list ( )
        vismol_object._generate_topology_from_index_bonds ( )
        vismol_object.define_molecules ( )
        vismol_object.define_Calpha_backbone ( )

    _refresh_bond_dependent_representations ( vismol_object )
    vismol_object.core_representations["picking_dots"] = None
    vismol_object.core_representations["picking_text"] = None

    vm_session = vismol_object.vm_session
    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.queue_draw ( )

    return created


def unset_bond ( vismol_object, atom_id_a, atom_id_b ):
    """ [EN] Companion to set_bond_order() above, for REMOVAL -- same
    reasoning as that function's own docstring for why this exists
    instead of reusing remove_bond() (which is Builder/manual_bonds-only
    and, per that function's own docstring, drops every bond not
    explicitly added via add_bond() when called on a normally-loaded
    structure).

    Works on ANY bond present in vismol_object.bonds, regardless of
    whether it came from automatic distance-based detection (a loaded
    file) or was added explicitly (add_bond()/set_bond_order() above) --
    removes it directly from that dict and from both atoms' own .bonds
    lists, then regenerates index_bonds/bond_order_list from what
    remains. No other bond is touched.

    Returns True if a bond was actually removed, False if that pair
    wasn't bonded to begin with. """
    if vismol_object.bonds is None:
        return False

    key = ( min ( atom_id_a, atom_id_b ), max ( atom_id_a, atom_id_b ) )
    bond = vismol_object.bonds.get ( key )
    if bond is None:
        return False

    del vismol_object.bonds[key]

    atom_a = vismol_object.atoms.get ( atom_id_a )
    atom_b = vismol_object.atoms.get ( atom_id_b )
    if atom_a is not None and bond in atom_a.bonds:
        atom_a.bonds.remove ( bond )
        atom_a.nbonds = len ( atom_a.bonds )
    if atom_b is not None and bond in atom_b.bonds:
        atom_b.bonds.remove ( bond )
        atom_b.nbonds = len ( atom_b.bonds )

    manual_bonds = getattr ( vismol_object, "manual_bonds", None )
    if manual_bonds:
        manual_bonds.discard ( key )
    manual_bond_orders = getattr ( vismol_object, "manual_bond_orders", None )
    if manual_bond_orders:
        manual_bond_orders.pop ( key, None )
    manual_aromatic_bonds = getattr ( vismol_object, "manual_aromatic_bonds", None )
    if manual_aromatic_bonds:
        manual_aromatic_bonds.discard ( key )

    _sync_index_bonds_and_order_list_from_bonds_dict ( vismol_object )

    vismol_object.non_bonded_atoms = None
    vismol_object._get_non_bonded_from_bonded_list ( )
    vismol_object._generate_topology_from_index_bonds ( )
    vismol_object.define_molecules ( )
    vismol_object.define_Calpha_backbone ( )

    _refresh_bond_dependent_representations ( vismol_object )
    vismol_object.core_representations["picking_dots"] = None
    vismol_object.core_representations["picking_text"] = None

    vm_session = vismol_object.vm_session
    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.queue_draw ( )

    return True


# =====================================================================================
#   Dynamic Bonds (representacao POR FRAME) -- 'bond'/'unbond' com frame=...
#   ------------------------------------------------------------------------------
#   Tudo acima (set_bond_order/unset_bond) edita a TOPOLOGIA ESTATICA do
#   object (vismol_object.bonds/index_bonds) -- applies to ALL frames and
#   e' o que fica gravado se o objeto for salvo/exportado.
#
#   As funcoes abaixo, em vez disso, editam vismol_object.dynamic_bonds[f]
#   -- a lista de pares POR FRAME usada pela representacao "Dynamic Bonds"
#   (tipicamente a regiao QC de uma trajetoria QM/MM, recalculada
#   automaticamente por distancia a cada frame -- ver VismolSession.
#   define_dynamic_bonds() / VismolObject.find_bonded_and_nonbonded_atoms()).
#
#   *** AVISO IMPORTANTE, repetido no docstring de cmd_bond/cmd_unbond no
#   terminal: using frame=... here does NOT create a real chemical bond.
#   E' PURAMENTE uma edicao da REPRESENTACAO/VISUALIZACAO para aquele(s)
#   frame(s) -- does not touch anything in the pDynamo system (topology, force
#   field, bond force constants, charges, etc.). To create a
#   real bond (with force-field parameters), the system's
#   sistema pDynamo precisa ser editada por outros meios -- isso aqui e'
#   so' para ajustar o que aparece na tela enquanto se inspeciona/prepara
#   a trajectory (e.g. visually forcing a bond that the automatic
#   distance-based detection missed in a specific frame).
# =====================================================================================

def resolve_frame_arg ( vismol_object, frame ):
    """ [EN] Converts the 'frame' argument accepted by 'bond'/'unbond'
    (terminal cmd_bond/cmd_unbond, or the picking-based equivalents in
    click_mode.py) into either:
      - None, meaning "no frame given -- edit the STATIC topology
        (vismol_object.bonds) via set_bond_order()/unset_bond() above,
        NOT Dynamic Bonds". This is the default or when frame= is simply
        not part of the command -- keeps existing 'bond'/'unbond' usage
        completely unchanged.
      - a sorted list of int frame indices to apply a Dynamic Bonds edit
        to (set_dynamic_bond_order()/unset_dynamic_bond() below).

    Accepted forms for `frame` (already run through the terminal DSL's
    own type coercion -- see Command._coerce() in easyhybrid_terminal.py):
      None                            -> None (static topology path)
      True (from 'frame=true'/'yes'/'on') -> [current frame only]
      5 (int, from 'frame=5')         -> [5]
      'all' (any case, 'frame=all')   -> every frame in vismol_object.frames
      '1:5' (str, 'frame=1:5')        -> frames 1..5 INCLUSIVE (both ends
                                          kept -- matches this codebase's
                                          existing 'resi=10-20' filter
                                          convention elsewhere, also
                                          inclusive; note the SEPARATOR
                                          here is ':' not '-', to avoid
                                          colliding with negative frame
                                          indices/typos)
      '7' (plain digit string, in case DSL coercion didn't catch it)
                                       -> [7]
    """
    if frame is None:
        return None

    vm_session = getattr ( vismol_object, "vm_session", None )
    frames_arr = getattr ( vismol_object, "frames", None )
    n_frames = int ( frames_arr.shape[0] ) if frames_arr is not None else 0

    if isinstance ( frame, bool ):
        current_frame = int ( vm_session.get_frame ( ) ) if ( vm_session is not None
                              and hasattr ( vm_session, "get_frame" ) ) else 0
        return [ current_frame ]

    if isinstance ( frame, ( int, float ) ):
        return [ int ( frame ) ]

    if isinstance ( frame, str ):
        s = frame.strip ( ).lower ( )
        if s in ( "all", "*" ):
            return list ( range ( n_frames ) )
        if ":" in s:
            a_str, b_str = s.split ( ":", 1 )
            try:
                a = int ( a_str.strip ( ) )
                b = int ( b_str.strip ( ) )
            except ValueError:
                raise ValueError ( "frame: intervalo invalido {!r} (use algo como '1:5').".format ( frame ) )
            if a > b:
                a, b = b, a
            return list ( range ( a, b + 1 ) )
        try:
            return [ int ( s ) ]
        except ValueError:
            raise ValueError ( "frame: valor invalido {!r} (use um inteiro, 'A:B' ou 'all').".format ( frame ) )

    raise ValueError ( "frame: tipo invalido {!r} (use um inteiro, 'A:B' ou 'all').".format ( frame ) )


def _refresh_dynamic_bond_representations ( vismol_object, affected_frames ):
    """ [EN] Marks any 'is_dynamic' representation (the Dynamic Bonds
    sticks/lines, see representations.py's Representation.__init__
    is_dynamic flag) as needing its index buffer reloaded, so a Dynamic
    Bonds edit made by set_dynamic_bond_order()/unset_dynamic_bond()
    below becomes visible immediately on screen -- WITHOUT this, the
    change is still correctly stored in vismol_object.dynamic_bonds[f]
    (and would show up correctly if the user steps away to another frame
    and back), but wouldn't redraw right away for a frame that's already
    on screen.

    Only bothers marking anything if the CURRENTLY DISPLAYED frame is
    among affected_frames -- editing frame 500 while looking at frame 0
    doesn't need an immediate redraw; that frame's geometry gets rebuilt
    naturally, the same way normal frame-stepping already does, whenever
    the user actually navigates there (is_dynamic representations always
    read vismol_object.dynamic_bonds[f] fresh when their index buffer is
    reloaded -- see representations.py's _load_ind_vbo()). """
    vm_session = getattr ( vismol_object, "vm_session", None )
    if vm_session is None or not hasattr ( vm_session, "get_frame" ):
        return
    current_frame = int ( vm_session.get_frame ( ) )
    if current_frame not in affected_frames:
        return
    for rep in ( getattr ( vismol_object, "representations", None ) or { } ).values ( ):
        if rep is not None and getattr ( rep, "is_dynamic", False ):
            rep.was_rep_ind_modified = True
    vm_glcore = getattr ( vm_session, "vm_glcore", None )
    if vm_glcore is not None:
        vm_glcore.queue_draw ( )


def set_dynamic_bond_order ( vismol_object, atom_id_a, atom_id_b, bond_order = 1, frames = None ):
    """ [EN] Dynamic Bonds equivalent of set_bond_order() above -- adds
    the pair to vismol_object.dynamic_bonds[f] (if not already present)
    and/or forces its order via vismol_object.dynamic_manual_bond_orders
    (see get_dynamic_bond_order_for_frame()'s own updated docstring),
    for every frame index in `frames`.

    *** Purely a REPRESENTATION edit -- see this section's own banner
    comment above. Does NOT touch vismol_object.bonds, the pDynamo
    system's topology, force-field parameters, or anything used for an
    actual QM/MM calculation.

    frames: iterable of int frame indices (already resolved -- see
    resolve_frame_arg() above; None is NOT accepted here, the caller is
    responsible for routing to set_bond_order() instead when frame=None).

    Returns the number of frames where a NEW pair had to be added (as
    opposed to a pair that already existed there and just had its order
    updated/confirmed). """
    if atom_id_a == atom_id_b:
        raise ValueError ( "set_dynamic_bond_order: cannot bond an atom to itself." )
    if atom_id_a not in vismol_object.atoms or atom_id_b not in vismol_object.atoms:
        raise ValueError ( "set_dynamic_bond_order: atom_id_a={} or atom_id_b={} does not exist in this object.".format (
                            atom_id_a, atom_id_b ) )
    if bond_order not in ( 1, 2, 3 ):
        raise ValueError ( "set_dynamic_bond_order: bond_order must be 1 (single), 2 (double) or 3 (triple), got {!r}.".format (
                            bond_order ) )
    if vismol_object.dynamic_bonds is None or len ( vismol_object.dynamic_bonds ) == 0:
        raise ValueError ( "set_dynamic_bond_order: this object has no Dynamic Bonds defined yet "
                            "(define a Dynamic Bonds selection first)." )
    if not frames:
        raise ValueError ( "set_dynamic_bond_order: no frame(s) given." )

    key = ( min ( atom_id_a, atom_id_b ), max ( atom_id_a, atom_id_b ) )

    if not hasattr ( vismol_object, "dynamic_manual_bond_orders" ) or vismol_object.dynamic_manual_bond_orders is None:
        vismol_object.dynamic_manual_bond_orders = { }

    n_created = 0
    touched = [ ]
    for f in frames:
        if f < 0 or f >= len ( vismol_object.dynamic_bonds ):
            continue
        flat = list ( np.asarray ( vismol_object.dynamic_bonds[f] ).ravel ( ).tolist ( ) )
        pairs_in_frame = { ( min ( flat[k], flat[k + 1] ), max ( flat[k], flat[k + 1] ) )
                            for k in range ( 0, len ( flat ), 2 ) }
        if key not in pairs_in_frame:
            flat.append ( key[0] )
            flat.append ( key[1] )
            vismol_object.dynamic_bonds[f] = flat
            n_created += 1

        vismol_object.dynamic_manual_bond_orders.setdefault ( f, { } )[key] = int ( bond_order )

        if vismol_object.dynamic_bond_orders is not None and f < len ( vismol_object.dynamic_bond_orders ):
            vismol_object.dynamic_bond_orders[f] = None  # invalida o cache -- ver get_dynamic_bond_order_for_frame

        touched.append ( f )

    _refresh_dynamic_bond_representations ( vismol_object, touched )

    return n_created


def unset_dynamic_bond ( vismol_object, atom_id_a, atom_id_b, frames = None ):
    """ [EN] Dynamic Bonds equivalent of unset_bond() above -- removes the
    pair from vismol_object.dynamic_bonds[f] (if present) and forgets any
    forced order for it, for every frame index in `frames`.

    *** Purely a REPRESENTATION edit -- see this section's own banner
    comment above.

    frames: iterable of int frame indices (already resolved -- see
    resolve_frame_arg() above; None is NOT accepted here).

    Returns the number of frames where the pair was actually present (and
    got removed). """
    if vismol_object.dynamic_bonds is None or len ( vismol_object.dynamic_bonds ) == 0:
        return 0
    if not frames:
        return 0

    key = ( min ( atom_id_a, atom_id_b ), max ( atom_id_a, atom_id_b ) )
    dmb = getattr ( vismol_object, "dynamic_manual_bond_orders", None )

    n_removed = 0
    touched = [ ]
    for f in frames:
        if f < 0 or f >= len ( vismol_object.dynamic_bonds ):
            continue
        flat = list ( np.asarray ( vismol_object.dynamic_bonds[f] ).ravel ( ).tolist ( ) )
        new_flat = [ ]
        removed_here = False
        for k in range ( 0, len ( flat ), 2 ):
            pair = ( min ( flat[k], flat[k + 1] ), max ( flat[k], flat[k + 1] ) )
            if pair == key and not removed_here:
                removed_here = True
                continue
            new_flat.append ( flat[k] )
            new_flat.append ( flat[k + 1] )

        if removed_here:
            vismol_object.dynamic_bonds[f] = new_flat
            n_removed += 1

        if dmb and f in dmb:
            dmb[f].pop ( key, None )

        if vismol_object.dynamic_bond_orders is not None and f < len ( vismol_object.dynamic_bond_orders ):
            vismol_object.dynamic_bond_orders[f] = None  # invalida o cache

        touched.append ( f )

    _refresh_dynamic_bond_representations ( vismol_object, touched )

    return n_removed


# =====================================================================================
#   Undo
#   ------------------------------------------------------------------------------
#   Snapshot-based, not command-based: before each "logical" user action
#   (place/replace an atom, a whole click-and-drag-to-create-a-bond
#   gesture, deleting an atom or bond, cycling a bond's order...), the
#   CALLER (click_mode.py / vismol_glcore.py's right-click handler) calls
#   push_undo_snapshot() once, capturing everything needed to reconstruct
#   the object's current state. undo() just pops the last one and rebuilds
#   the object from it via add_atom() -- reusing the exact same function
#   that already correctly handles chains/residues/frames/atom_dic_id/
#   representations, rather than trying to clone live Atom/Chain/Residue
#   Python objects directly (those hold back-references to each other and
#   to vismol_object itself, which is a much easier thing to get subtly
#   wrong -- and impossible for me to fully verify without a live GL
#   context to test against).
#
#   Deliberately does NOT snapshot on every single mouse_motion event
#   during a drag (would flood the stack with one entry per pixel of
#   movement, and undo-ing would need dozens of clicks to get anywhere) --
#   only once per gesture, at the point the caller considers the action
#   "started" (e.g. start_bond_drag(), not update_bond_drag()).
# =====================================================================================

def _snapshot_builder_state ( vismol_object ):
    """ [EN] Captures everything needed to reconstruct vismol_object's
    CURRENT state via _restore_builder_state() below: every atom's
    element/name/chain/residue/position, plus the explicit bond set and
    any custom bond orders (manual_bonds/manual_bond_orders -- see
    add_atom()'s docstring for why bonds are entirely explicit now).
    Does NOT snapshot representations/VAOs/etc -- those get rebuilt fresh
    on restore anyway, same as every other atom_ops function already
    does after a mutation. """
    atoms_snapshot = [ ]
    for atom_id in sorted ( vismol_object.atoms.keys ( ) ):
        atom = vismol_object.atoms[atom_id]
        pos  = vismol_object.frames[0, atom_id]
        atoms_snapshot.append ( {
            'symbol'  : atom.symbol,
            'name'    : atom.name,
            'chain_id': atom.chain.name if atom.chain is not None else 'A',
            'resi'    : atom.residue.index if atom.residue is not None else 1,
            'resn'    : atom.residue.name if atom.residue is not None else 'UNK',
            'x'       : float ( pos[0] ), 'y': float ( pos[1] ), 'z': float ( pos[2] ),
        } )

    return {
        'atoms'                 : atoms_snapshot,
        'manual_bonds'          : set ( getattr ( vismol_object, 'manual_bonds', None ) or set ( ) ),
        'manual_bond_orders'    : dict ( getattr ( vismol_object, 'manual_bond_orders', None ) or { } ),
        'manual_aromatic_bonds' : set ( getattr ( vismol_object, 'manual_aromatic_bonds', None ) or set ( ) ),
    }


def push_undo_snapshot ( vismol_object, max_depth = 50 ):
    """ Records the object's CURRENT state onto its undo stack, BEFORE
    the caller goes on to mutate it. Call this once per logical user
    action (see the module-level note above), never per low-level
    primitive (a single "replace element" click might internally call
    set_atom_element() once -- one snapshot; a drag calls add_atom() +
    several move_atom()s + add_bond() -- still just ONE snapshot, pushed
    by start_bond_drag() before any of that starts).

    max_depth caps memory/undo-stack growth for a long editing session --
    oldest snapshots are dropped first once exceeded. """
    if not hasattr ( vismol_object, 'undo_stack' ) or vismol_object.undo_stack is None:
        vismol_object.undo_stack = [ ]

    vismol_object.undo_stack.append ( _snapshot_builder_state ( vismol_object ) )

    if len ( vismol_object.undo_stack ) > max_depth:
        vismol_object.undo_stack.pop ( 0 )


def _restore_builder_state ( vismol_object, snapshot ):
    """ Rebuilds vismol_object from scratch to match `snapshot`, reusing
    add_atom() for every atom (recompute_bonds=False,
    update_representation=False while looping, to avoid paying for a
    full rebuild after EVERY individual atom -- one single rebuild at the
    end is enough, same batching add_atom()'s own docstring already
    recommends for adding several atoms in a row). """
    # [EN] BUG FIX: every atom about to be discarded here was registered in
    # vm_session.atom_dic_id (by add_atom(), keyed by unique_id -- see
    # remove_atom()'s own docstring above for why that dict exists). This
    # rebuild replaces vismol_object.atoms wholesale and re-adds fresh atoms
    # via add_atom() below, each getting a brand-new unique_id -- but until
    # now nothing ever popped the OLD atoms' entries first, so every single
    # undo() leaked one dangling atom_dic_id entry per atom that existed
    # before the undo (unbounded growth over a long editing+undo session).
    vm_session = vismol_object.vm_session
    for old_atom in vismol_object.atoms.values ( ):
        if old_atom.unique_id in vm_session.atom_dic_id:
            del vm_session.atom_dic_id[old_atom.unique_id]

    vismol_object.atoms  = { }
    vismol_object.chains = { }
    vismol_object.frames = np.zeros ( ( 1, 0, 3 ), dtype = np.float32 )
    vismol_object.manual_bonds          = set ( )
    vismol_object.manual_bond_orders    = { }
    vismol_object.manual_aromatic_bonds = set ( )

    for atom_data in snapshot['atoms']:
        add_atom ( vismol_object,
                   symbol   = atom_data['symbol'],
                   x        = atom_data['x'], y = atom_data['y'], z = atom_data['z'],
                   name     = atom_data['name'],
                   chain_id = atom_data['chain_id'],
                   resi     = atom_data['resi'],
                   resn     = atom_data['resn'],
                   recompute_bonds       = False,
                   update_representation = False )

    vismol_object.manual_bonds          = set ( snapshot['manual_bonds'] )
    vismol_object.manual_bond_orders    = dict ( snapshot['manual_bond_orders'] )
    vismol_object.manual_aromatic_bonds = set ( snapshot.get ( 'manual_aromatic_bonds', ( ) ) or ( ) )

    vismol_object.cov_radii_array = None
    vismol_object.electronegativity_array = None
    vismol_object.index_bonds = None
    vismol_object.bonds = None
    vismol_object.non_bonded_atoms = None
    _reapply_manual_bonds ( vismol_object )

    vismol_object.create_representation ( rep_type = "sticks" )
    vismol_object.create_representation ( rep_type = "stick_spheres" )
    _activate_new_sphere_representation ( vismol_object, "stick_spheres" )
    vismol_object.create_representation ( rep_type = "nonbonded" )
    vismol_object.core_representations["picking_dots"] = None
    vismol_object.core_representations["picking_text"] = None

    # [EN] Atom count just changed -- see _refresh_bond_dependent_
    # representations()'s own (extended) docstring: any OTHER
    # representation this object carries (e.g. 'cartoon', inherited by
    # an "Edit in Builder" clone of a normal system) needs to either be
    # resized or deactivated here too, or its own stale self.indexes
    # crashes the next picking pass.
    _refresh_bond_dependent_representations ( vismol_object )

    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.queue_draw ( )

    # [EN] undo() rebuilds vismol_object completely from scratch (every
    # atom re-added via add_atom() above) -- the linked pDynamo system
    # would otherwise be left describing the PRE-undo structure. Import
    # done here (not at module level) to avoid a circular import:
    # empty_object.py's sync_pdynamo_system() itself calls back into
    # THIS module (_reapply_manual_bonds() indirectly, via add_atom()).
    from gui.windows.builder.empty_object import sync_pdynamo_system
    sync_pdynamo_system ( vismol_object )


def undo ( vismol_object ):
    """ Pops and restores the most recent snapshot (see
    push_undo_snapshot()), undoing the last logical action. Returns True
    if something was actually undone, False if the stack was already
    empty (nothing left to undo -- caller should probably disable/grey
    out the Undo button in that case, see the sidebar wiring). """
    undo_stack = getattr ( vismol_object, 'undo_stack', None )
    if not undo_stack:
        return False

    snapshot = undo_stack.pop ( )
    _restore_builder_state ( vismol_object, snapshot )
    return True


# =====================================================================================
#   Automatic hydrogen adjustment
#   ------------------------------------------------------------------------------
#   Generalises the ad-hoc auto-hydrogenation that used to live ONLY in
#   click_mode.handle_click_to_place_atom() (the `tmp` dict of C/N/O
#   direction templates, only ever applied to a brand-new, completely
#   isolated atom). adjust_hydrogens() instead works for ANY atom at ANY
#   point -- newly placed, freshly bonded to something via a drag,
#   replaced to a different element, had a bond's order cycled, or lost a
#   bond/neighbour to a deletion -- by comparing its CURRENT total bond
#   order (to non-H neighbours) against its standard valence, and adding
#   or removing H atoms to match. Called explicitly by click_mode.py /
#   vismol_glcore.py after whichever operation may have changed an atom's
#   bonding (see each call site's own comment for which atom(s) it's
#   called on and why).
# =====================================================================================

# Default valence (maximum sum of bond order) per element -- the
# poucos elementos que este primeiro conjunto de ferramentas do Builder
# realmente produz ate agora (C/N/O/H, mais alguns halogenios/S/P comuns
# so as not to leave it manual once the element selector grows). Elements
# outside this table simply do not get their hydrogens adjusted (see
# adjust_hydrogens() -- retorna sem fazer nada nesse caso, em vez de
# guessing a wrong valence).
STANDARD_VALENCE = {
    'H' : 1,
    'C' : 4,
    'N' : 3,
    'O' : 2,
    'F' : 1,
    'CL': 1,
    'BR': 1,
    'I' : 1,
    'P' : 3,
    'S' : 2,

    # [EN] EXTENDED (user's own explicit request, alongside the new full
    # periodic-table picker -- periodic_table_dialog.py -- "use os
    # parametros do UFF como referencia"): every OTHER per-element datum
    # this Builder needs (covalent/van der Waals radius, colour,
    # electronegativity -- read by Atom._init_cov_rad()/_init_vdw_rad()/
    # _init_color()/_init_electronegativity() in vismol/model/atom.py)
    # ALREADY comes from vismol's own PeriodicTable table for the WHOLE
    # periodic table, and that table's own radius/electronegativity
    # columns are themselves UFF-derived (confirmed directly against
    # pDynamo3's own DYFF force field, whose atomTypes.yaml explicitly
    # says "Atom types copied from UFF" and cites the exact same Rappe
    # et al. 1992 source) -- so bond-length rendering and later DYFF
    # assignment already work correctly for ANY element the new picker
    # offers, with NO changes needed here. This table alone (STANDARD_
    # VALENCE, used only by adjust_hydrogens()'s automatic H add/remove)
    # was the one piece still narrower than that -- extended here with
    # the remaining MAIN-GROUP elements that have a single, well-defined
    # "typical" covalent valence in ordinary closed-shell compounds (the
    # same classic 8-N-rule valence UFF itself assumes for its own
    # generic sp3/sp2/sp atom types), same reasoning as the original 10.
    #
    # Deliberately NOT extended to transition metals, lanthanides,
    # actinides, or noble gases: none of those have one single "typical"
    # covalent valence the way a p-block or alkali/alkaline-earth element
    # does (oxidation state varies by compound) -- adjust_hydrogens()
    # already handles "not in this table" gracefully (returns without
    # guessing, per this module's own docstring above), so the picker
    # still lets the user PLACE those atoms freely, they just don't get
    # automatic hydrogen bookkeeping, exactly like any other unsupported
    # element today.
    'B' : 3,
    'AL': 3,
    'SI': 4,
    'GE': 4,
    'SN': 4,
    'PB': 4,
    'AS': 3,
    'SB': 3,
    'BI': 3,
    'SE': 2,
    'TE': 2,
    'AT': 1,
    'LI': 1,
    'NA': 1,
    'K' : 1,
    'RB': 1,
    'CS': 1,
    'MG': 2,
    'CA': 2,
    'SR': 2,
    'BA': 2,
}


def _bond_order_sum_excluding_symbol ( vismol_object, atom_id, exclude_symbol ):
    """ Sum of bond_order over every bond touching atom_id, EXCLUDING
    bonds to a neighbour whose symbol == exclude_symbol (used to exclude
    hydrogens -- see adjust_hydrogens(), which needs "how much valence is
    already used up by NON-hydrogen neighbours" to know how many
    hydrogens are left to place). """
    total = 0
    for bond in vismol_object.bonds.values ( ):
        if bond.atom_index_i != atom_id and bond.atom_index_j != atom_id:
            continue
        other_id = bond.atom_index_j if bond.atom_index_i == atom_id else bond.atom_index_i
        if vismol_object.atoms[other_id].symbol == exclude_symbol:
            continue
        total += bond.bond_order
    return total


def _bonded_neighbors_with_symbol ( vismol_object, atom_id, symbol ):
    """ atom_ids of every neighbour of atom_id whose symbol == symbol
    (used to find atom_id's CURRENT hydrogens -- see adjust_hydrogens()). """
    neighbor_ids = [ ]
    for bond in vismol_object.bonds.values ( ):
        if bond.atom_index_i != atom_id and bond.atom_index_j != atom_id:
            continue
        other_id = bond.atom_index_j if bond.atom_index_i == atom_id else bond.atom_index_i
        if vismol_object.atoms[other_id].symbol == symbol:
            neighbor_ids.append ( other_id )
    return neighbor_ids


def _ideal_angle_rule ( n, max_order, double_count = 0 ):
    """ [EN] Pure "which angle" decision -- the ONE place this simple,
    deliberately-not-a-real-hybridisation-calculation VSEPR-style rule
    lives, shared by BOTH _ideal_angle_for_atom() (heavy-atom skeleton
    relaxation, reads n/max_order/double_count off vismol_object's
    CURRENT bonds) and adjust_hydrogens() (decides the geometry for
    hydrogens not yet placed -- reads n/max_order/double_count off the
    atom's HEAVY bonds plus its PROSPECTIVE hydrogen count, since
    hydrogens are always single bonds and never raise max_order/
    double_count themselves). Sharing this one function is what makes
    the heavy-atom skeleton and its hydrogens agree on the SAME
    hybridisation for the same atom, instead of two independently
    drifting copies of this rule.

    n : total number of substituents (bonded neighbour ATOMS -- each
        counted once regardless of that bond's order).
    max_order : highest bond order among this atom's HEAVY bonds (1 if
        it has none).
    double_count : how many of this atom's HEAVY bonds are order 2
        (only matters for the n==2 allene-like case below).

    Returns the ideal angle in DEGREES, or None if n is outside [2, 4]
    (no angle to speak of, or more substituents than this Builder's
    supported elements ever need). """
    if n < 2 or n > 4:
        return None
    if n == 2:
        if max_order == 3 or double_count == 2:
            return 180.0
        elif max_order == 2:
            return 120.0
        else:
            return 109.5
    elif n == 3:
        if max_order >= 2:
            return 120.0
        else:
            return 109.5
    else:
        return 109.5


# Regular-tetrahedron unit vectors (pairwise ~109.47 degrees apart) --
# the "no existing neighbours yet" starting orientation for a
# 4-substituent (sp3) centre, and (a 3-vector subset of it) for a
# 3-substituent PYRAMIDAL (non-planar) centre -- see
# _complete_vsepr_directions() below.
_TETRAHEDRAL_TEMPLATE = [ v / np.linalg.norm ( v ) for v in (
    np.array ( [-0.785298,  0.243518, -0.653254] ), np.array ( [ 0.322015, -0.981331, -0.189814] ),
    np.array ( [-0.334691,  0.073016,  0.992665] ), np.array ( [ 0.798227,  0.665009, -0.149645] ),
) ]


def _heavy_neighbor_directions_and_orders ( vismol_object, atom_id ):
    """ [EN] Unit-vector directions (from atom_id's CURRENT position) and
    bond orders to every bonded neighbour that is NOT a hydrogen -- the
    "fixed anchors" adjust_hydrogens() completes a VSEPR arrangement
    around (see _complete_vsepr_directions()). Hydrogens are deliberately
    excluded here: adjust_hydrogens() always recomputes ALL of atom_id's
    OWN hydrogens fresh (see its own docstring), so an existing
    hydrogen's (possibly stale/wrong) direction must never be treated as
    a fixed anchor for that recomputation. """
    atom_pos = vismol_object.frames[0, atom_id]
    dirs, orders = [ ], [ ]
    for bond in vismol_object.bonds.values ( ):
        if bond.atom_index_i != atom_id and bond.atom_index_j != atom_id:
            continue
        other_id = bond.atom_index_j if bond.atom_index_i == atom_id else bond.atom_index_i
        if vismol_object.atoms[other_id].symbol == 'H':
            continue
        vec  = vismol_object.frames[0, other_id] - atom_pos
        norm = float ( np.linalg.norm ( vec ) )
        if norm > 1e-6:
            dirs.append ( vec / norm )
            orders.append ( bond.bond_order )
    return dirs, orders


def _arbitrary_orthonormal_basis ( v ):
    """ Two unit vectors (e1, e2), both perpendicular to `v` and to each
    other -- used to place new VSEPR directions at a specific angle
    around a single fixed axis `v` (see _complete_vsepr_directions()'s
    n_existing==1 cases). Which exact e1/e2 pair is picked is arbitrary
    (there's no "correct" azimuth around a single bond, chemically) --
    only that they're orthonormal to `v` and to each other, which is all
    the geometry below actually needs. """
    e1 = _arbitrary_perpendicular ( v )
    e2 = np.cross ( v, e1 )
    return e1, e2 / np.linalg.norm ( e2 )


def _fallback_fan_directions ( existing_dirs, n_new, ideal_angle_deg ):
    """ [EN] The ORIGINAL, simpler heuristic this Builder used before
    _complete_vsepr_directions() below existed: places n_new directions
    fanned out, symmetrically, around the single direction AWAY from the
    (normalised sum of the) existing ones -- doesn't try to reproduce
    exact VSEPR angles, just avoids new directions landing on top of
    existing bonds. Only reached now as a degenerate-case fallback (see
    _complete_vsepr_directions()'s own docstring) -- e.g. existing
    directions numerically collinear/opposite, where the closed-form
    construction's divisions become unstable. """
    if existing_dirs:
        avg_dir = -np.sum ( existing_dirs, axis = 0 )
        norm = float ( np.linalg.norm ( avg_dir ) )
        avg_dir = ( avg_dir / norm ) if norm > 1e-6 else np.array ( [ 0.0, 0.0, 1.0 ] )
    else:
        avg_dir = np.array ( [ 0.0, 0.0, 1.0 ] )

    if n_new == 1:
        return [ avg_dir ]

    perp = np.cross ( avg_dir, np.array ( [ 0.0, 0.0, 1.0 ] ) )
    if float ( np.linalg.norm ( perp ) ) < 1e-6:
        perp = np.cross ( avg_dir, np.array ( [ 0.0, 1.0, 0.0 ] ) )
    perp = perp / np.linalg.norm ( perp )

    spread = np.deg2rad ( ideal_angle_deg if ideal_angle_deg else 109.5 )
    directions = [ ]
    for k in range ( n_new ):
        offset = ( k - ( n_new - 1 ) / 2.0 ) * ( spread / max ( 1, n_new - 1 ) )
        d = avg_dir * np.cos ( offset ) + perp * np.sin ( offset )
        directions.append ( d / np.linalg.norm ( d ) )
    return directions


def _complete_vsepr_directions ( existing_dirs, total_count, ideal_angle_deg ):
    """ [EN] Returns (total_count - len(existing_dirs)) NEW unit-vector
    directions that complete a symmetric VSEPR-style arrangement of
    `total_count` substituents around a central atom, given the
    directions ALREADY fixed (existing_dirs -- 0 to 3 of them, e.g. this
    atom's existing heavy-neighbour bonds) and the ideal pairwise angle
    for that arrangement (from _ideal_angle_rule()).

    [EN] REPLACES the old adjust_hydrogens()/_new_hydrogen_directions()
    heuristic ("place away from whatever's there, fanned out by an
    arbitrary fixed 50-degree angle, regardless of hybridisation") --
    this is closed-form geometry (not a generic numerical solver, and
    not a real force field either) for every (len(existing_dirs),
    total_count) combination this Builder's STANDARD_VALENCE table can
    actually produce (total_count 1-4, existing 0-3), each verified
    numerically to reproduce the EXACT ideal_angle_deg pairwise, both
    between a new direction and an existing one, and between two new
    directions:

      existing=0: an arbitrarily-oriented canonical template for that
        total_count/angle (nothing to align to yet).
      existing=1 (axis A): new directions D_k = cos(theta)*A +
        sin(theta)*(cos(az_k) e1 + sin(az_k) e2), theta=ideal_angle,
        (e1, e2) an arbitrary orthonormal basis perpendicular to A.
        Azimuths az_k are 0/phi for 2 new directions (phi solved from
        cos(phi) = cos(theta) / (1 + cos(theta)), so the 2 new
        directions are ALSO ideal_angle apart from EACH OTHER, not just
        from A), or 0/120/240 degrees for 3 new directions (a special
        case that's only self-consistent because ideal_angle_deg is
        always 109.5 whenever 3 new directions are needed here).
      existing=2 (A, B): decomposed in the orthonormal basis
        {u=normalize(A+B), n=normalize(cross(A,B))} -- by the mirror
        symmetry between A and B, every new direction has a ZERO
        component along the 3rd axis (perpendicular to both u and n),
        so this reduces to 2 unknowns solved directly from "new
        direction must be ideal_angle from A" (and, by that same
        symmetry, automatically ideal_angle from B too).
      existing=3: the classic "4th regular-tetrahedron vertex is minus
        the sum of the other 3, already unit length" identity (only
        ever reached when total_count==4, i.e. ideal_angle_deg==109.5).

    Falls back to _fallback_fan_directions() above -- correctness isn't
    guaranteed there, but it degrades gracefully instead of dividing by
    zero -- for any input this closed-form approach can't handle cleanly
    (existing directions numerically collinear/opposite, or an
    existing_count/total_count combination outside the cases above,
    e.g. an atom with more bonds than any element this Builder supports
    should chemically have). """
    n_existing = len ( existing_dirs )
    n_new      = total_count - n_existing
    if n_new <= 0:
        return [ ]

    existing_dirs = [ np.asarray ( d, dtype = np.float64 ) / np.linalg.norm ( np.asarray ( d, dtype = np.float64 ) )
                       for d in existing_dirs ]
    theta = np.radians ( ideal_angle_deg ) if ideal_angle_deg is not None else None

    try:
        if n_existing == 0:
            if theta is None:
                if total_count == 1:
                    return [ np.array ( [ 0.0, 0.0, 1.0 ] ) ]
            elif total_count == 2:
                half = theta / 2.0
                return [ np.array ( [  np.sin ( half ), 0.0, np.cos ( half ) ] ),
                         np.array ( [ -np.sin ( half ), 0.0, np.cos ( half ) ] ) ]
            elif total_count == 3 and abs ( ideal_angle_deg - 120.0 ) < 1e-6:
                return [ np.array ( [ np.cos ( np.radians ( 120.0 * k ) ),
                                       np.sin ( np.radians ( 120.0 * k ) ), 0.0 ] ) for k in range ( 3 ) ]
            elif total_count == 3:
                return list ( _TETRAHEDRAL_TEMPLATE[:3] )
            elif total_count == 4:
                return list ( _TETRAHEDRAL_TEMPLATE )

        elif n_existing == 1 and theta is not None:
            A = existing_dirs[0]
            e1, e2 = _arbitrary_orthonormal_basis ( A )
            if n_new == 1:
                return [ np.cos ( theta ) * A + np.sin ( theta ) * e1 ]
            elif n_new == 2:
                cos_phi = np.clip ( np.cos ( theta ) / ( 1.0 + np.cos ( theta ) ), -1.0, 1.0 )
                phi = np.arccos ( cos_phi )
                return [ np.cos ( theta ) * A + np.sin ( theta ) * ( np.cos ( az ) * e1 + np.sin ( az ) * e2 )
                         for az in ( 0.0, phi ) ]
            elif n_new == 3:
                return [ np.cos ( theta ) * A + np.sin ( theta ) * ( np.cos ( az ) * e1 + np.sin ( az ) * e2 )
                         for az in ( 0.0, 2 * np.pi / 3.0, 4 * np.pi / 3.0 ) ]

        elif n_existing == 2 and theta is not None:
            A, B = existing_dirs
            cos_ab   = float ( np.clip ( np.dot ( A, B ), -1.0, 1.0 ) )
            u_raw    = A + B
            u_norm   = float ( np.linalg.norm ( u_raw ) )
            cross_ab = np.cross ( A, B )
            n_norm   = float ( np.linalg.norm ( cross_ab ) )
            if u_norm < 1e-3 or n_norm < 1e-3:
                raise ValueError ( "existing directions too close to collinear/opposite" )
            u = u_raw / u_norm
            n = cross_ab / n_norm
            c = float ( np.cos ( np.arccos ( cos_ab ) / 2.0 ) )   # = dot(A, u) = dot(B, u)
            if abs ( c ) < 1e-3:
                raise ValueError ( "existing directions too close to a 180-degree pair" )
            if n_new == 1 and ideal_angle_deg is not None and abs ( ideal_angle_deg - 120.0 ) < 1e-6:
                # [EN] BUG FIX (user's own explicit request -- "para o
                # cleanUP, nao esta ajustando corretamente a geometria de
                # H de carbonos SP2, use a disposicao dos outros atomos
                # ligados ao carbono sp2 para colocar o H no mesmo
                # plano"): the GENERAL cone formula below (still used for
                # every OTHER n_existing==2/n_new==1 case, e.g. completing
                # a PYRAMIDAL 109.5-degree centre like a plain sp3 amine
                # nitrogen, where the 3rd substituent SHOULD pucker out of
                # the A/B plane) finds a direction that's EXACTLY
                # ideal_angle_deg from BOTH A and B -- which only happens
                # to land in-plane (b=0 below) when A and B themselves are
                # ALREADY exactly 120 degrees apart. For any REAL,
                # imperfectly-relaxed sp2 centre (A/B not exactly 120
                # degrees apart -- the normal case), that formula instead
                # returns one of the family's two OUT-OF-PLANE mirror
                # solutions, visibly wrong for something that's supposed
                # to be planar (confirmed live: clean_up_structure()'s own
                # new planarity term, added earlier this session, flattens
                # the HEAVY skeleton correctly, but this H-placement
                # formula was still pushing the hydrogen back out of that
                # same plane on every Phase-2 pass).
                #
                # Fixed by NOT trying to hit ideal_angle_deg from each of
                # A/B exactly -- instead directly using the same "new
                # direction = negative sum of the existing ones,
                # normalised" identity the n_existing==3 tetrahedral case
                # below already uses for its own 4th-vertex completion:
                # -normalize(A+B) is, by construction, a pure LINEAR
                # COMBINATION of A and B (zero component along their
                # shared plane's normal n) -- i.e. ALWAYS exactly
                # coplanar with A and B, for ANY angle between them, not
                # just the ideal one -- while still pointing directly
                # "away" from both, exactly bisecting the reflex angle
                # between them (the chemically correct trigonal
                # completion direction).
                s = -( A + B )
                s_norm = float ( np.linalg.norm ( s ) )
                if s_norm < 1e-3:
                    raise ValueError ( "existing directions sum to ~zero (A/B ~180 degrees apart)" )
                return [ s / s_norm ]
            elif n_new == 1:
                a = np.cos ( theta ) / c
                if abs ( a ) > 1.0:
                    # [EN] existing A/B are too far from ideal_angle apart
                    # from EACH OTHER for a valid single completion at
                    # EXACTLY ideal_angle to both to even exist (found via
                    # live testing: a freshly-built, not-yet-relaxed chain
                    # -- e.g. an as-placed C-C-C with a ~169-degree angle,
                    # nowhere near the 109.5 target -- makes |a| > 1 here,
                    # which would otherwise silently produce a NON-unit
                    # direction) -- fall back instead of returning garbage.
                    raise ValueError ( "no valid single completion for this A/B pair" )
                b = np.sqrt ( max ( 0.0, 1.0 - a * a ) )
                return [ a * u + b * n ]
            elif n_new == 2:
                p = np.cos ( theta ) / c
                if abs ( p ) > 1.0:
                    # [EN] same reasoning as the n_new==1 guard above --
                    # additionally, this is the case that actually
                    # surfaced the bug: |p| > 1 collapses D1/D2 to the
                    # SAME (non-unit, wrong-length) direction (r would be
                    # computed as 0 via the max(0, ...) clamp below,
                    # silently), so both new hydrogens would otherwise
                    # land on top of each other instead of falling back.
                    raise ValueError ( "no valid 2-direction completion for this A/B pair" )
                r = np.sqrt ( max ( 0.0, 1.0 - p * p ) )
                return [ p * u + r * n, p * u - r * n ]

        elif n_existing == 3 and n_new == 1:
            s = -np.sum ( existing_dirs, axis = 0 )
            norm = float ( np.linalg.norm ( s ) )
            if norm < 1e-3:
                raise ValueError ( "existing directions sum to ~zero" )
            return [ s / norm ]

    except ValueError:
        pass

    return _fallback_fan_directions ( existing_dirs, n_new, ideal_angle_deg )


def _hydrogen_bond_length ( vismol_object, heavy_atom ):
    """ Covalent-radius-sum ideal length for a (single) bond from
    heavy_atom to a hydrogen -- same model _ideal_bond_length() below
    already uses for heavy-heavy bonds, just looked up directly since
    the new H atom doesn't exist yet at the point this is needed (there's
    no Atom object yet to read .cov_rad off of). Replaces the old flat
    1.05-Angstrom-for-everything default (e.g. a C-H bond and an O-H bond
    are chemically quite different lengths -- ~1.09 vs ~0.96 Angstrom;
    this model gives ~1.11 and ~1.03 respectively, much closer than a
    single constant for both). """
    h_cov_rad = vismol_object.vm_session.periodic_table.elements_by_symbol['H'][7]
    return float ( heavy_atom.cov_rad ) + float ( h_cov_rad )


def adjust_hydrogens ( vismol_object, atom_id ):
    """ Recomputes ALL of atom_id's own hydrogens -- adding, removing,
    AND REPOSITIONING as needed -- so that (a) their COUNT matches
    atom_id's STANDARD_VALENCE (as before), and (b) their POSITIONS form
    a proper VSEPR-style arrangement together with atom_id's existing
    HEAVY-atom bonds (_heavy_neighbor_directions_and_orders() /
    _complete_vsepr_directions() / _ideal_angle_rule() above) --
    respecting the SAME sp/sp2/sp3-ish hybridisation rule the heavy-atom
    skeleton relaxation uses (_ideal_angle_for_atom()), instead of the
    old "just place away from whatever's already there, fanned out by an
    arbitrary fixed angle" heuristic.

    [EN] BEHAVIOUR CHANGE: this used to return early (doing nothing) as
    soon as the CURRENT hydrogen count already matched the target --
    correct for COUNT, but meant a hydrogen's POSITION, once placed
    (e.g. under the OLD heuristic, or before a neighbouring bond's order
    later changed this atom's hybridisation), was never revisited even
    when it no longer made geometric sense. Now ALWAYS recomputes every
    hydrogen's position from the atom's CURRENT heavy-neighbour geometry
    on every call -- called explicitly after any operation that may have
    changed atom_id's bonding OR its neighbours' positions (see the
    module-level note above for the full list of call sites), and now
    also called directly by clean_up_structure()'s own hydrogen-
    repositioning phase (see that function's docstring).

    No-op if atom_id's element isn't in STANDARD_VALENCE (unknown
    elements are left alone rather than guessed at), or if it needs (and
    already has) exactly zero hydrogens.

    Returns the number of hydrogens added (positive) or removed
    (negative), or 0 if the count didn't change (positions may still
    have been adjusted even when this returns 0). """
    if atom_id not in vismol_object.atoms:
        return 0

    atom = vismol_object.atoms[atom_id]
    target_valence = STANDARD_VALENCE.get ( atom.symbol.upper ( ) )
    if target_valence is None:
        return 0

    heavy_bond_order_sum = _bond_order_sum_excluding_symbol ( vismol_object, atom_id, exclude_symbol = 'H' )
    needed_h = max ( 0, target_valence - heavy_bond_order_sum )

    current_h_ids = sorted ( _bonded_neighbors_with_symbol ( vismol_object, atom_id, 'H' ) )
    initial_h_count = len ( current_h_ids )

    if needed_h == 0 and initial_h_count == 0:
        return 0

    # [EN] Removes any EXCESS hydrogens FIRST (highest atom_id downward:
    # remove_atom() renumbers every id ABOVE the one it removes, so
    # removing highest-first means none of the other hydrogens still
    # pending removal/reposition below are affected by an earlier
    # removal -- only atom_id itself might shift, tracked explicitly).
    if initial_h_count > needed_h:
        for h_id in sorted ( current_h_ids[needed_h:], reverse = True ):
            remove_atom ( vismol_object, h_id )
            if h_id < atom_id:
                atom_id -= 1
        current_h_ids = current_h_ids[:needed_h]
        atom = vismol_object.atoms[atom_id]

    heavy_dirs, heavy_orders = _heavy_neighbor_directions_and_orders ( vismol_object, atom_id )
    total_count  = len ( heavy_dirs ) + needed_h
    max_order    = max ( heavy_orders ) if heavy_orders else 1
    double_count = sum ( 1 for o in heavy_orders if o == 2 )
    ideal_angle  = _ideal_angle_rule ( total_count, max_order, double_count )

    fresh_dirs = _complete_vsepr_directions ( heavy_dirs, total_count, ideal_angle ) if needed_h > 0 else [ ]
    bond_length = _hydrogen_bond_length ( vismol_object, atom )
    pos = vismol_object.frames[0, atom_id]

    # [EN] Repositions every SURVIVING hydrogen to one of the freshly
    # computed directions (arbitrary correspondence -- they're identical
    # H atoms, so it doesn't matter which specific one gets which
    # direction), then adds brand-new ones for any remaining directions.
    for k, h_id in enumerate ( current_h_ids ):
        offset = fresh_dirs[k] * bond_length
        move_atom ( vismol_object, h_id,
                    float ( pos[0] + offset[0] ), float ( pos[1] + offset[1] ), float ( pos[2] + offset[2] ) )

    for offset_dir in fresh_dirs[len ( current_h_ids ):]:
        offset = offset_dir * bond_length
        add_atom ( vismol_object, symbol = "H",
                   x = float ( pos[0] + offset[0] ),
                   y = float ( pos[1] + offset[1] ),
                   z = float ( pos[2] + offset[2] ),
                   bonded_to = atom_id )

    return needed_h - initial_h_count


# =====================================================================================
#   Structure clean-up (simple geometric relaxation -- NOT a force field)
#   ------------------------------------------------------------------------------
#   Deliberately NOT UFF (or any other real force field): no energy, no
#   gradient, no van der Waals/electrostatics/torsions -- just two simple
#   geometric corrections, applied iteratively:
#     1. Bond LENGTHS pulled toward an ideal length (covalent-radius sum,
#        shortened a bit for double/triple bonds).
#     2. Bond ANGLES at each atom pulled toward an ideal angle, chosen
#        from a simple VSEPR-style rule based on how many neighbours the
#        atom has and whether any of its bonds are double/triple.
#   Each iteration computes ALL corrections first (from the CURRENT,
#   not-yet-updated positions) and only applies them at the end (a
#   Jacobi-style update) -- avoids the order-dependence/oscillation that
#   applying corrections one bond/angle at a time, immediately, would
#   cause when several of them share an atom.
# =====================================================================================

# Fator de encurtamento aplicado a soma dos raios covalentes conforme a
# bond order -- approximate values, just to give a reasonable
# geometry (does not intend to reproduce the literature precisely).
_BOND_ORDER_LENGTH_FACTOR = { 1: 1.00, 2: 0.87, 3: 0.78 }


def _ideal_bond_length ( atom_i, atom_j, bond_order ):
    """ Soma dos raios covalentes dos dois atomos, encurtada conforme a
    ordem da ligacao (dupla/tripla mais curtas que simples). """
    base   = float ( atom_i.cov_rad ) + float ( atom_j.cov_rad )
    factor = _BOND_ORDER_LENGTH_FACTOR.get ( bond_order, 1.00 )
    return base * factor


def bond_order_from_distance ( atom_i, atom_j, distance ):
    """ [EN] Classifies `distance` (the CURRENT/final distance between
    atom_i and atom_j, e.g. mid-drag or at drag-release) as whichever
    bond order (1, 2 or 3) has the CLOSEST _ideal_bond_length() to it --
    short drag snaps to triple, medium to double, long to single, with
    the exact thresholds simply being the midpoints between each pair of
    consecutive ideal lengths for THESE TWO SPECIFIC elements (not a
    single universal cutoff -- a C-C and an O-H pair have very different
    absolute bond lengths, so the classification has to be relative to
    each element pair's own covalent radii, reusing the exact same model
    _ideal_bond_length() already uses for clean_up_structure()'s bond-
    length relaxation, for consistency).

    Used by click_mode's drag-to-bond-by-distance feature (Builder
    roadmap section 4) -- see start_bond_drag()/update_bond_drag()/
    finish_bond_drag() for where this is actually wired in, and their
    own docstrings for how this interacts with the sidebar's explicit
    Single/Double/Triple/Aromatic selection (distance only decides the
    order when the sidebar is at its default, Single/non-aromatic --
    an explicit sidebar choice always wins over the drag distance). """
    candidates = [ ( order, _ideal_bond_length ( atom_i, atom_j, order ) ) for order in ( 1, 2, 3 ) ]
    return min ( candidates, key = lambda pair: abs ( pair[1] - distance ) )[0]


def _ideal_angle_for_atom ( vismol_object, atom_id ):
    """ [EN] Gathers atom_id's CURRENT neighbour count/bond orders and
    delegates the actual "which angle" decision to _ideal_angle_rule()
    above -- the SAME rule adjust_hydrogens() uses when placing this
    atom's hydrogens, so a heavy atom's own skeleton angle and the angle
    its hydrogens get placed at always agree (one shared hybridisation
    decision, not two independently-drifting copies of the same logic).

    Returns the ideal angle in DEGREES, or None if atom_id has fewer
    than 2 neighbours (no angle to speak of) or more than 4 (unusual for
    the elements this Builder currently supports -- left alone rather
    than guessed at). """
    neighbor_bonds = [ bond for bond in vismol_object.bonds.values ( )
                       if bond.atom_index_i == atom_id or bond.atom_index_j == atom_id ]
    n = len ( neighbor_bonds )
    if n < 2 or n > 4:
        return None

    max_order    = max ( bond.bond_order for bond in neighbor_bonds )
    double_count = sum ( 1 for bond in neighbor_bonds if bond.bond_order == 2 )
    return _ideal_angle_rule ( n, max_order, double_count )


def _rotate_vector_rodrigues ( v, axis, angle_rad ):
    """ Rotaciona o vetor `v` em torno de `axis` (unitario) por
    `angle_rad`, via formula de Rodrigues -- preserva o comprimento de
    `v`, so muda a direcao. """
    cos_a = np.cos ( angle_rad )
    sin_a = np.sin ( angle_rad )
    return ( v * cos_a
             + np.cross ( axis, v ) * sin_a
             + axis * np.dot ( axis, v ) * ( 1.0 - cos_a ) )


def _arbitrary_perpendicular ( v ):
    """ Um vetor unitario qualquer, perpendicular a `v` -- usado quando
    dois vizinhos de um mesmo atomo estao quase colineares (produto
    vetorial proximo de zero), caso em que o eixo de rotacao "natural"
    (perpendicular aos dois) fica mal-definido. """
    axis = np.cross ( v, np.array ( [ 0.0, 0.0, 1.0 ] ) )
    if float ( np.linalg.norm ( axis ) ) < 1e-6:
        axis = np.cross ( v, np.array ( [ 0.0, 1.0, 0.0 ] ) )
    return axis / np.linalg.norm ( axis )


def clean_up_structure ( vismol_object, atom_ids = None, n_iterations = 40, step_scale = 0.7,
                          adjust_hydrogen_count = True ):
    """ Fixes up a molecule's geometry in TWO PHASES, heavy atoms first,
    then hydrogens -- NOT a real force field like UFF (see the module-
    level note above), just simple, deliberately cheap geometric
    correction.

    [EN] TWO-PHASE REDESIGN (this used to be a single phase that moved
    every atom -- heavy AND hydrogen -- together in the same Jacobi-style
    correction loop; found, via the user's own testing, to give visibly
    bad results): the OLD single-phase design let hydrogens (whose
    "ideal" direction is comparatively arbitrary before their heavy
    neighbour's OWN geometry has settled) pull on/skew their heavy
    neighbour's bond angles every iteration, and vice-versa -- both
    categories fighting over the same shared correction budget instead
    of converging to a sensible geometry. Splitting into two ordered
    phases removes that interference entirely:

      Phase 1 relaxes ONLY the heavy-atom skeleton (bond lengths/angles
      between non-hydrogen atoms) -- hydrogens are excluded completely
      here, neither moved nor used as a bond/angle target.

      Phase 2 then calls adjust_hydrogens() (see its own docstring -- now
      hybridisation-aware: builds a proper VSEPR arrangement from each
      heavy atom's bond pattern via _complete_vsepr_directions() /
      _ideal_angle_rule(), instead of the old fixed-angle fan heuristic)
      for every heavy atom in scope, using the skeleton phase 1 just
      finalised as fixed anchors -- i.e. exactly the "first fix the
      heavy atoms, then place H according to the bond pattern/
      hybridisation" order requested -- unless adjust_hydrogen_count is
      False (see its own parameter doc below), in which case this whole
      phase is skipped entirely.

    [EN] PLANARITY TERM (user's own follow-up request -- "os carbonos
    SP2 nao ficam planos"): Phase 1's angle correction, on its own, only
    pulls PAIRS of a trigonal centre's bond vectors toward 120 degrees --
    that alone doesn't pin down the group's third rotational degree of
    freedom (all 3 individually ~120 degrees apart from each other, yet
    still twisted out of their own shared plane as a "tripod"). Added a
    THIRD term to the same Phase-1 iteration (see the inline comment at
    its own call site, in the "angle term" loop below) that fits a plane
    through every 120-degree (trigonal) heavy centre with >=3 heavy
    neighbours and pulls each neighbour's out-of-plane component back to
    zero -- genuinely flattens sp2 centres now, not just their pairwise
    angles. Centres with fewer than 3 heavy neighbours (e.g. a terminal
    >C=CH2 alkene carbon) still rely on Phase 2's VSEPR placement alone,
    unchanged, since a plane can't be meaningfully fit through only 1-2
    heavy positions (hydrogens are excluded from Phase 1 entirely).

    [EN] Phase 1's n_iterations=40 / step_scale=0.7 defaults are
    unchanged from before this redesign (verified with a standalone
    numpy simulation predating the hydrogen-exclusion change -- a
    deliberately worst-case, badly distorted 5-atom test converged to
    within +/-0.7 degrees of the ideal tetrahedral angle and correct
    bond lengths to 3 decimal places). Excluding hydrogens from the
    correction loop only REMOVES simultaneous constraints competing for
    the same atom's displacement budget, so convergence with these same
    defaults is at least as good, not worse.

    atom_ids : if None (default), relaxes the WHOLE molecule. If a
               collection of atom_ids is given instead, only those atoms
               (PLUS their directly bonded neighbours -- moving an
               atom's bonds/angles without letting its neighbours
               respond too wouldn't converge to anything sensible) are
               in scope: free to move in phase 1 (if heavy), and have
               their own hydrogens recomputed in phase 2 (if heavy).
               Every OTHER atom stays fixed and acts as an anchor for
               bonds/angles that touch the in-scope set.
    n_iterations : how many phase-1 correction passes to run. More
               iterations converge closer to the ideal geometry but take
               proportionally longer -- 40 is a reasonable default for
               small Builder molecules (see the verification note above).
    step_scale : how much of each phase-1 iteration's computed
               correction is actually applied (0-1). Less than 1.0
               (default 0.7) acts as a damping factor, trading
               convergence speed for stability (a full, undamped
               correction every iteration can overshoot and oscillate,
               especially once bond and angle corrections on the same
               atom start interacting).
    adjust_hydrogen_count : if True (default), Phase 2 runs as described
               above -- every heavy atom in scope gets its hydrogens
               fully recomputed (repositioned, AND added/removed to
               match standard valence -- see adjust_hydrogens()'s own
               docstring). If False (new -- sidebar checkbox, user's own
               request), Phase 2 is skipped ENTIRELY: existing hydrogens
               are left exactly where they were (not even repositioned
               to follow a heavy atom Phase 1 just moved), and none are
               added or removed. For a molecule whose current hydrogen
               count is already intentional/correct (e.g. just imported
               from a source the Builder shouldn't second-guess), this
               avoids Clean Up silently changing which atoms exist, not
               just where they sit.

    Rebuilds representations once at the very end (this is a one-off,
    on-demand action -- e.g. the sidebar's "Clean Up" button -- not a
    per-frame update, so there's no reason to batch/optimise beyond
    that single rebuild; phase 2's own adjust_hydrogens() calls also
    each trigger their own rebuild internally whenever they add/remove a
    hydrogen -- see that function's docstring -- this final rebuild just
    also covers phase 1's own, otherwise-unreflected position changes). """
    if atom_ids is None:
        scope_ids = set ( vismol_object.atoms.keys ( ) )
    else:
        scope_ids = set ( atom_ids )
        for aid in list ( atom_ids ):
            for bond in vismol_object.bonds.values ( ):
                if bond.atom_index_i == aid:
                    scope_ids.add ( bond.atom_index_j )
                elif bond.atom_index_j == aid:
                    scope_ids.add ( bond.atom_index_i )

    if not scope_ids:
        return

    # =========================== Phase 1: heavy-atom skeleton ===========================
    movable_ids = { aid for aid in scope_ids if vismol_object.atoms[aid].symbol != 'H' }

    positions = { aid: np.array ( vismol_object.frames[0, aid], dtype = np.float64 )
                  for aid in vismol_object.atoms.keys ( ) }

    for _iteration in range ( n_iterations ) if movable_ids else ( ):
        displacement = { aid: np.zeros ( 3, dtype = np.float64 ) for aid in movable_ids }
        contributions = { aid: 0 for aid in movable_ids }

        # --- distance term (heavy-heavy bonds only -- hydrogens excluded) ---
        for bond in vismol_object.bonds.values ( ):
            i, j = bond.atom_index_i, bond.atom_index_j
            if vismol_object.atoms[i].symbol == 'H' or vismol_object.atoms[j].symbol == 'H':
                continue
            movable_i = i in movable_ids
            movable_j = j in movable_ids
            if not movable_i and not movable_j:
                continue

            vec  = positions[j] - positions[i]
            dist = float ( np.linalg.norm ( vec ) )
            if dist < 1e-6:
                continue

            ideal = _ideal_bond_length ( vismol_object.atoms[i], vismol_object.atoms[j], bond.bond_order )
            direction = vec / dist
            delta = direction * ( ideal - dist )   # o quanto J precisaria se mover pra afastar/aproximar de I

            if movable_i and movable_j:
                displacement[i] += -0.5 * delta ; contributions[i] += 1
                displacement[j] +=  0.5 * delta ; contributions[j] += 1
            elif movable_j:
                displacement[j] += delta ; contributions[j] += 1
            else:
                displacement[i] += -delta ; contributions[i] += 1

        # --- angle term (heavy neighbours only, at heavy centres only) ---
        for atom_id in vismol_object.atoms.keys ( ):
            if vismol_object.atoms[atom_id].symbol == 'H':
                continue
            ideal_angle = _ideal_angle_for_atom ( vismol_object, atom_id )
            if ideal_angle is None:
                continue

            neighbor_ids = [ ]
            for bond in vismol_object.bonds.values ( ):
                if bond.atom_index_i == atom_id and vismol_object.atoms[bond.atom_index_j].symbol != 'H':
                    neighbor_ids.append ( bond.atom_index_j )
                elif bond.atom_index_j == atom_id and vismol_object.atoms[bond.atom_index_i].symbol != 'H':
                    neighbor_ids.append ( bond.atom_index_i )

            center_pos = positions[atom_id]

            # --- planarity term (trigonal/sp2 heavy centres with >=3
            # heavy neighbours only -- see this function's own extended
            # docstring: pairwise ANGLE correction alone, right below,
            # can converge to 3 substituents that are each individually
            # ~120 degrees apart yet still out of plane as a group (a
            # "tripod", not a flat trigonal centre) -- angles between
            # PAIRS of vectors don't pin down the third rotational
            # degree of freedom (how the whole trio twists out of its
            # own average plane). This term adds that missing
            # constraint directly: fit a plane through the centre atom
            # using all of its heavy neighbours' current directions
            # (a sum-of-pairwise-cross-products normal -- the natural
            # generalisation of the single cross-product normal a plain
            # 3-neighbour case would use, but also well-defined/robust
            # for a 4th, non-ideal heavy neighbour if one happens to be
            # present), then pulls each neighbour's OUT-OF-PLANE
            # component back toward zero. Left for Phase 2's own VSEPR-
            # based adjust_hydrogens() to handle -- unchanged -- whenever
            # fewer than 3 heavy neighbours exist (e.g. a terminal
            # >C=CH2 alkene carbon), since a plane can't be meaningfully
            # fit through only 1-2 heavy positions here (hydrogens are
            # deliberately excluded from Phase 1 entirely, see the
            # module docstring's own two-phase reasoning).
            if ideal_angle == 120.0 and len ( neighbor_ids ) >= 3:
                unit_vecs = { }
                for nid in neighbor_ids:
                    v = positions[nid] - center_pos
                    vnorm = float ( np.linalg.norm ( v ) )
                    if vnorm > 1e-6:
                        unit_vecs[nid] = v / vnorm
                if len ( unit_vecs ) >= 3:
                    plane_normal = np.zeros ( 3 )
                    ids = list ( unit_vecs.keys ( ) )
                    for idx_a in range ( len ( ids ) ):
                        for idx_b in range ( idx_a + 1, len ( ids ) ):
                            plane_normal += np.cross ( unit_vecs[ids[idx_a]], unit_vecs[ids[idx_b]] )
                    normal_len = float ( np.linalg.norm ( plane_normal ) )
                    if normal_len > 1e-6:
                        plane_normal = plane_normal / normal_len
                        for nid in ids:
                            if nid not in movable_ids:
                                continue
                            out_of_plane = float ( np.dot ( positions[nid] - center_pos, plane_normal ) )
                            if abs ( out_of_plane ) < 1e-4:
                                continue
                            displacement[nid] += -out_of_plane * plane_normal * step_scale
                            contributions[nid] += 1

            for idx_a in range ( len ( neighbor_ids ) ):
                for idx_b in range ( idx_a + 1, len ( neighbor_ids ) ):
                    j = neighbor_ids[idx_a]
                    k = neighbor_ids[idx_b]
                    if j not in movable_ids and k not in movable_ids:
                        continue

                    v1 = positions[j] - center_pos
                    v2 = positions[k] - center_pos
                    n1 = float ( np.linalg.norm ( v1 ) )
                    n2 = float ( np.linalg.norm ( v2 ) )
                    if n1 < 1e-6 or n2 < 1e-6:
                        continue

                    cos_theta = np.clip ( np.dot ( v1, v2 ) / ( n1 * n2 ), -1.0, 1.0 )
                    theta = float ( np.arccos ( cos_theta ) )
                    diff  = np.deg2rad ( ideal_angle ) - theta
                    if abs ( diff ) < 1e-4:
                        continue

                    axis = np.cross ( v1, v2 )
                    axis_norm = float ( np.linalg.norm ( axis ) )
                    if axis_norm < 1e-6:
                        axis = _arbitrary_perpendicular ( v1 )
                    else:
                        axis = axis / axis_norm

                    half = ( diff / 2.0 ) * step_scale
                    new_v1 = _rotate_vector_rodrigues ( v1, axis, -half )
                    new_v2 = _rotate_vector_rodrigues ( v2, axis,  half )

                    if j in movable_ids:
                        displacement[j] += ( center_pos + new_v1 ) - positions[j]
                        contributions[j] += 1
                    if k in movable_ids:
                        displacement[k] += ( center_pos + new_v2 ) - positions[k]
                        contributions[k] += 1

        # --- aplica a media de todas as correcoes acumuladas nesta iteracao ---
        for aid in movable_ids:
            if contributions[aid] > 0:
                positions[aid] = positions[aid] + displacement[aid] / contributions[aid]

    for aid in movable_ids:
        vismol_object.frames[0, aid] = positions[aid].astype ( np.float32 )

    # =========================== Phase 2: hydrogens ===========================
    # [EN] For every HEAVY atom in scope, recompute ALL of its hydrogens
    # fresh from the now-finalised heavy skeleton (adjust_hydrogens() --
    # see its own docstring for why it's now safe/correct to call
    # unconditionally: it always repositions, not just count-adjusts).
    # Uses the ATOM OBJECTS themselves (not cached integer ids): an
    # earlier iteration's adjust_hydrogens() call can remove excess
    # hydrogens, which renumbers every atom_id ABOVE the removed one
    # (see remove_atom()'s own docstring) -- reading .atom_id fresh off
    # each heavy atom's own Atom object, right before using it, is what
    # finish_bond_drag() already does for the exact same reason, and
    # stays correct regardless of what an earlier iteration shifted.
    if adjust_hydrogen_count:
        heavy_atoms_in_scope = [ vismol_object.atoms[aid] for aid in scope_ids
                                  if vismol_object.atoms[aid].symbol != 'H' ]
        for heavy_atom in heavy_atoms_in_scope:
            adjust_hydrogens ( vismol_object, heavy_atom.atom_id )

    vismol_object.mass_center = np.mean ( vismol_object.frames[0], axis = 0 )

    vismol_object.create_representation ( rep_type = "sticks" )
    vismol_object.create_representation ( rep_type = "stick_spheres" )
    _activate_new_sphere_representation ( vismol_object, "stick_spheres" )
    vismol_object.create_representation ( rep_type = "nonbonded" )
    vismol_object.core_representations["picking_dots"] = None
    vismol_object.core_representations["picking_text"] = None

    vm_session = vismol_object.vm_session
    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.queue_draw ( )


# =====================================================================================
#   DYFF quick geometry optimisation -- a REAL force-field minimisation (unlike
#   clean_up_structure() above, which is a cheap geometric heuristic that never
#   actually evaluates/minimises an energy), run on a disposable, throwaway
#   pDynamo System built from the CURRENT Builder state. The scratch System is
#   never registered with p_session (no e_id, no psystem entry, no treeview
#   row) -- only the optimised coordinates are copied back, then it is simply
#   left to be garbage-collected.
# =====================================================================================

def optimize_geometry_dyff ( vismol_object, maximum_iterations = 100, rms_gradient_tolerance = 0.5 ):
    """ [EN] "Optimize (DYFF)" button -- user's own request: "vamos colocar a
    possibilidade de otimizar a geometria com o DYFF, vamos criar um sistema
    provisorio no background, otimizar por 100 passos e atualizar as
    coordenadas, depois podemos descartar o sistema provisorio." -- a quick,
    REAL force-field relaxation from inside the Builder, as an alternative to
    clean_up_structure()'s purely geometric fix.

    1. Builds a scratch pDynamo System from vismol_object's CURRENT atoms/
       bonds/positions via empty_object._build_pdynamo_system_from_vismol_
       object() -- the exact same builder every real Builder-to-pDynamo sync
       already uses (bond orders/aromaticity come from manual_bonds/
       manual_bond_orders/manual_aromatic_bonds, same as always).
    2. Assigns DYFF to it via p_session.define_MMModel() -- reusing that
       method's own guanidinium auto-correction and (see below) the Atom
       Types window's manual per-atom overrides, instead of duplicating any
       of that logic here.
    3. Runs ConjugateGradientMinimize_SystemGeometry() capped at
       `maximum_iterations` steps (default 100, per the user's own request),
       silently (log=None) -- this is a quick in-Builder convenience, not a
       logged/trajectoried job through the Process Manager.
    4. Copies the resulting coordinates straight into vismol_object.
       frames[0] and discards the scratch System -- it is never added to
       p_session.psystem, never gets its own e_id, no new treeview row.

    Manual Atom-Types-window overrides: define_MMModel() only looks up
    per-atom overrides by matching `system.e_id` against vm_objects_dic
    (see its own docstring) -- the scratch System has no e_id of its own,
    so it's temporarily tagged with vismol_object's OWN e_id (if it already
    has one, from a prior Builder sync) right before calling
    define_MMModel(), purely so that lookup succeeds; never written back
    anywhere, no effect beyond this one call.

    Only positions change (bonds/atom count/representation TYPES are
    untouched), so this reuses move_atom()'s cheap "just re-upload
    coordinates" refresh (was_rep_coord_modified/was_sel_coord_modified on
    every active representation) rather than clean_up_structure()'s full
    create_representation() rebuild, which isn't needed here.

    Returns (True, message) on success, (False, message) on failure -- e.g.
    DYFF can't type some atom's environment (same MMModelError define_
    MMModel() itself already turns into a clean message rather than a raw
    traceback). The molecule is left COMPLETELY UNTOUCHED on failure: the
    scratch System is discarded before any coordinates are copied back, and
    a minimiser exception part-way through still copies back whatever
    partial progress it made rather than losing it. """
    main      = vismol_object.vm_session.main
    p_session = main.p_session

    from gui.windows.builder import empty_object
    scratch_system = empty_object._build_pdynamo_system_from_vismol_object (
        vismol_object, label = "builder_dyff_optimisation_scratch" )

    existing_e_id = getattr ( vismol_object, "e_id", None )
    if existing_e_id is not None:
        scratch_system.e_id = existing_e_id

    ok, message = p_session.define_MMModel ( force_field = 'DYFF', system = scratch_system )
    if not ok:
        return False, message

    from pSimulation import ConjugateGradientMinimize_SystemGeometry
    try:
        ConjugateGradientMinimize_SystemGeometry ( scratch_system,
                                                     maximumIterations    = maximum_iterations,
                                                     rmsGradientTolerance = rms_gradient_tolerance,
                                                     log                  = None )
    except Exception as exc:
        dprint ( "WARNING atom_ops.optimize_geometry_dyff: minimisation raised part-way through -- "
                 "keeping whatever geometry it reached so far. Error: {}".format ( exc ) )

    n_atoms = len ( vismol_object.atoms )
    for atom_id in range ( n_atoms ):
        xyz = scratch_system.coordinates3[atom_id]
        vismol_object.frames[0, atom_id] = [ float ( xyz[0] ), float ( xyz[1] ), float ( xyz[2] ) ]

    vismol_object.mass_center = np.mean ( vismol_object.frames[0], axis = 0 )

    vm_session = vismol_object.vm_session
    for rep in vismol_object.representations.values ( ):
        if rep is not None and rep.active:
            rep.was_rep_coord_modified = True
            rep.was_sel_coord_modified = True
    picking_dots = vismol_object.core_representations.get ( "picking_dots" )
    if picking_dots is not None:
        picking_dots.was_rep_coord_modified = True

    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.updated_coords = True
        vm_session.vm_glcore.queue_draw ( )

    return True, "DYFF geometry optimisation finished ({} iteration(s) max).".format ( maximum_iterations )


# =====================================================================================
#   Fragment library -- attach a small pre-built molecule (fragment_library.
#   load_fragment()'s return value) onto a hydrogen atom of the object being
#   edited, replacing that H with the fragment's own "root" (attachment
#   point) atom. See fragment_library.py's own top-of-file docstring for
#   the fragment-file convention (exactly one atom left one bond short of
#   its standard valence = the attachment point, auto-detected there, not
#   here) and click_mode.handle_click_to_attach_fragment() for the click-
#   driven UI this is wired into.
# =====================================================================================

def attach_fragment_at_hydrogen ( vismol_object, target_h_atom_id, fragment ):
    """ [EN] Removes target_h_atom_id (which MUST be a hydrogen) and
    attaches `fragment` (a dict from fragment_library.load_fragment()) in
    its place -- the fragment's own "root" atom (fragment["root_index"])
    ends up bonded, via a plain single bond, to whatever heavy atom
    target_h_atom_id was bonded to, oriented along the EXACT SAME
    direction that hydrogen already occupied (that direction was already
    placed correctly by adjust_hydrogens()'s own VSEPR logic, so it's
    reused directly here rather than recomputed).

    Geometry (see fragment_library.py's own docstring for why fragments
    are authored with completely arbitrary global orientation -- only
    the fragment's OWN internal relative angles matter):
      1. Gather the fragment root atom's existing bonded-neighbour
         directions, IN THE FRAGMENT's OWN LOCAL COORDINATES, and
         complete the missing (attachment) direction via the SAME VSEPR
         machinery adjust_hydrogens() itself uses (_ideal_angle_rule() /
         _complete_vsepr_directions()) -- this is the fragment's own
         "local outward direction" for the new bond.
      2. Build a single rigid rotation mapping that local outward
         direction onto the target's own outward direction (the removed
         hydrogen's direction from its parent), via _rotate_vector_
         rodrigues() -- same helper clean_up_structure()'s angle term
         already uses for the same kind of vector-to-vector rotation.
      3. Apply that ONE rotation to every fragment atom (relative to the
         root atom's own local position), then translate the whole
         fragment so the root atom lands at parent_position +
         outward_dir * (parent.cov_rad + root_element_cov_rad) -- the
         same cov-radius-sum single-bond-length model _ideal_bond_
         length()/_hydrogen_bond_length() already use, just generalised
         to the fragment root's OWN element instead of always H.

    Insertion: batch-inserts every fragment atom via add_atom(
    recompute_bonds=False, update_representation=False) (same batching
    convention add_atom()'s own docstring recommends for adding several
    atoms in a row), tracking a local-index -> new-atom_id map; registers
    every INTERNAL fragment bond plus the ONE new cross-bond directly into
    manual_bonds/manual_bond_orders (same direct-population approach
    _restore_builder_state() already uses for a from-scratch rebuild);
    then ONE _reapply_manual_bonds() and ONE full representation rebuild
    (matching remove_atom()'s own tail exactly, since the atom COUNT
    changes here same as there); then adjust_hydrogens() on both the
    target's parent atom and the fragment's newly-attached root atom
    (their valence usage just changed: parent lost one H but gained one
    new substituent of the SAME order, and the fragment root gained the
    one bond it was deliberately authored one short of).

    Does NOT push an undo snapshot or call sync_pdynamo_system() itself
    -- matches this file's existing convention (see add_atom()/
    remove_atom()'s own docstrings): that is the CALLER's job (see
    click_mode.handle_click_to_attach_fragment()).

    Raises ValueError if target_h_atom_id doesn't exist, isn't a
    hydrogen, or has no parent bond (should never happen for a real atom
    in a well-formed object, but checked defensively since this is
    reachable directly from a user click). """
    if target_h_atom_id not in vismol_object.atoms:
        raise ValueError ( "attach_fragment_at_hydrogen: atom_id {} does not exist.".format ( target_h_atom_id ) )
    h_atom = vismol_object.atoms[target_h_atom_id]
    if h_atom.symbol != 'H':
        raise ValueError ( "attach_fragment_at_hydrogen: atom_id {} ('{}') is not a hydrogen.".format (
                            target_h_atom_id, h_atom.symbol ) )

    parent_atom = None
    for bond in vismol_object.bonds.values ( ):
        if bond.atom_index_i == target_h_atom_id:
            parent_atom = vismol_object.atoms[bond.atom_index_j] ; break
        elif bond.atom_index_j == target_h_atom_id:
            parent_atom = vismol_object.atoms[bond.atom_index_i] ; break
    if parent_atom is None:
        raise ValueError ( "attach_fragment_at_hydrogen: atom_id {} has no parent bond.".format ( target_h_atom_id ) )

    h_pos      = np.array ( vismol_object.frames[0, target_h_atom_id], dtype = np.float64 )
    parent_pos = np.array ( vismol_object.frames[0, parent_atom.atom_id], dtype = np.float64 )
    outward_dir = h_pos - parent_pos
    outward_norm = float ( np.linalg.norm ( outward_dir ) )
    if outward_norm < 1e-6:
        raise ValueError ( "attach_fragment_at_hydrogen: degenerate (zero-length) H-parent bond." )
    outward_dir = outward_dir / outward_norm

    frag_atoms  = fragment["atoms"]
    frag_bonds  = fragment["bonds"]
    root_index  = fragment["root_index"]

    root_symbol, rx, ry, rz = frag_atoms[root_index]
    root_local_pos = np.array ( [ rx, ry, rz ], dtype = np.float64 )

    existing_dirs   = [ ]
    existing_orders = [ ]
    for ( i, j, order ) in frag_bonds:
        other = j if i == root_index else ( i if j == root_index else None )
        if other is None:
            continue
        _sym, ox, oy, oz = frag_atoms[other]
        vec = np.array ( [ ox, oy, oz ], dtype = np.float64 ) - root_local_pos
        vnorm = float ( np.linalg.norm ( vec ) )
        if vnorm > 1e-6:
            existing_dirs.append ( vec / vnorm )
            existing_orders.append ( order )

    total_count  = len ( existing_dirs ) + 1
    max_order    = max ( existing_orders ) if existing_orders else 1
    double_count = sum ( 1 for o in existing_orders if o == 2 )
    ideal_angle  = _ideal_angle_rule ( total_count, max_order, double_count )
    local_outward = _complete_vsepr_directions ( existing_dirs, total_count, ideal_angle )[0]
    local_outward = local_outward / np.linalg.norm ( local_outward )

    # [EN] Rigid rotation mapping local_outward (fragment space) onto
    # outward_dir (target space) -- standard axis/angle construction,
    # with the 2 degenerate cases (already aligned, or exactly opposite)
    # handled explicitly since cross(a, a) / cross(a, -a) are both zero
    # vectors and can't supply a rotation axis on their own.
    cos_theta = float ( np.clip ( np.dot ( local_outward, outward_dir ), -1.0, 1.0 ) )
    axis = np.cross ( local_outward, outward_dir )
    axis_norm = float ( np.linalg.norm ( axis ) )
    if axis_norm < 1e-8:
        if cos_theta > 0:
            rotation_axis, rotation_angle = np.array ( [ 0.0, 0.0, 1.0 ] ), 0.0
        else:
            rotation_axis, rotation_angle = _arbitrary_perpendicular ( local_outward ), np.pi
    else:
        rotation_axis, rotation_angle = axis / axis_norm, float ( np.arccos ( cos_theta ) )

    root_cov_rad = vismol_object.vm_session.periodic_table.elements_by_symbol[root_symbol][7]
    bond_length  = float ( parent_atom.cov_rad ) + float ( root_cov_rad )
    new_root_pos = parent_pos + outward_dir * bond_length

    # [EN] Remove the target hydrogen FIRST (captures parent_atom's own
    # .atom_id fresh afterwards -- remove_atom() renumbers every id above
    # the removed one, and mutates every SURVIVING atom's .atom_id IN
    # PLACE, so re-reading it off the live Atom object is always correct,
    # same convention finish_bond_drag()/clean_up_structure() already use).
    remove_atom ( vismol_object, target_h_atom_id )

    # [EN] Batch-insert every fragment atom -- world position = rotate
    # (this atom's position relative to the fragment's own root atom) by
    # the single rotation computed above, then translate to new_root_pos.
    # The root atom's OWN relative position is (0,0,0) by construction,
    # so it lands exactly at new_root_pos, as intended.
    local_index_to_new_id = { }
    for local_index, ( symbol, x, y, z ) in enumerate ( frag_atoms ):
        local_pos = np.array ( [ x, y, z ], dtype = np.float64 ) - root_local_pos
        rotated   = _rotate_vector_rodrigues ( local_pos, rotation_axis, rotation_angle )
        world_pos = new_root_pos + rotated
        new_atom = add_atom ( vismol_object, symbol = symbol,
                               x = float ( world_pos[0] ), y = float ( world_pos[1] ), z = float ( world_pos[2] ),
                               recompute_bonds = False, update_representation = False )
        local_index_to_new_id[local_index] = new_atom.atom_id

    if not hasattr ( vismol_object, "manual_bonds" ) or vismol_object.manual_bonds is None:
        vismol_object.manual_bonds = set ( )
    if not hasattr ( vismol_object, "manual_bond_orders" ) or vismol_object.manual_bond_orders is None:
        vismol_object.manual_bond_orders = { }

    for ( i, j, order ) in frag_bonds:
        pair = ( min ( local_index_to_new_id[i], local_index_to_new_id[j] ),
                 max ( local_index_to_new_id[i], local_index_to_new_id[j] ) )
        vismol_object.manual_bonds.add ( pair )
        vismol_object.manual_bond_orders[pair] = order

    root_new_id = local_index_to_new_id[root_index]
    cross_pair = ( min ( parent_atom.atom_id, root_new_id ), max ( parent_atom.atom_id, root_new_id ) )
    vismol_object.manual_bonds.add ( cross_pair )
    vismol_object.manual_bond_orders[cross_pair] = 1

    _reapply_manual_bonds ( vismol_object )

    vismol_object.create_representation ( rep_type = "sticks" )
    vismol_object.create_representation ( rep_type = "stick_spheres" )
    _activate_new_sphere_representation ( vismol_object, "stick_spheres" )
    vismol_object.create_representation ( rep_type = "nonbonded" )
    vismol_object.core_representations["picking_dots"] = None
    vismol_object.core_representations["picking_text"] = None

    # [EN] Atom count just changed -- see _refresh_bond_dependent_
    # representations()'s own (extended) docstring: any OTHER
    # representation this object carries (e.g. 'cartoon', inherited by
    # an "Edit in Builder" clone of a normal system) needs to either be
    # resized or deactivated here too, or its own stale self.indexes
    # crashes the next picking pass.
    _refresh_bond_dependent_representations ( vismol_object )

    adjust_hydrogens ( vismol_object, parent_atom.atom_id )
    adjust_hydrogens ( vismol_object, root_new_id )

    vm_session = vismol_object.vm_session
    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.queue_draw ( )


# =====================================================================================
#   Structure library -- user's own request: "adicionar estruturas prontas ...
#   semelhante a ferramenta de fragmentos, mas nao tem necessidade de ter um H
#   como referencia." Places a COMPLETE molecule (structure_library.
#   load_structure()'s return value) as a whole at a given position, unlike
#   attach_fragment_at_hydrogen() above -- no attachment atom, no VSEPR
#   rotation/alignment, no hydrogen-count adjustment: the structure file
#   already IS a finished molecule, atoms/bonds are inserted exactly as
#   authored, just translated to the click location.
# =====================================================================================

def add_structure_at_position ( vismol_object, structure, x, y, z ):
    """ [EN] Inserts every atom/bond of `structure` (a dict from structure_
    library.load_structure()) into vismol_object, translated so the
    structure's own CENTROID (the plain average of its atoms' file
    coordinates -- not atom #1, so a structure authored off-center still
    ends up centred under the click) lands at (x, y, z) -- the Builder's
    "Add Structure" tool's click_mode.handle_click_to_add_structure() call
    site computes (x, y, z) the exact same way handle_click_to_place_atom()
    already does for a plain new atom on empty space (world_pos_from_mouse()
    + the target object's own inverse model_mat).

    Unlike attach_fragment_at_hydrogen() (which removes a hydrogen, rotates
    the fragment onto that bond's direction via VSEPR, and re-adjusts both
    ends' hydrogens afterward): no atom is removed, no rotation is applied
    (the structure keeps its own file authored orientation -- a real,
    complete molecule's internal geometry is exactly what should be
    preserved, unlike a fragment's arbitrary local frame), and adjust_
    hydrogens() is never called (the structure's own hydrogens, exactly as
    the file declares them, ARE the finished molecule -- there is no
    "missing valence" to complete, unlike a fragment's one intentionally
    open attachment point).

    Batch-inserts every atom via add_atom(recompute_bonds=False,
    update_representation=False) (same batching convention attach_fragment_
    at_hydrogen() already uses for the same reason: one full representation
    rebuild at the end, not N of them), tracking a local-index -> new-
    atom_id map, then registers every structure bond directly into manual_
    bonds/manual_bond_orders/manual_aromatic_bonds (same direct-population
    approach attach_fragment_at_hydrogen() and _restore_builder_state() both
    already use) before ONE _reapply_manual_bonds() + ONE full
    representation rebuild.

    Does NOT push an undo snapshot or call sync_pdynamo_system() itself --
    matches this file's existing convention (see add_atom()/attach_
    fragment_at_hydrogen()'s own docstrings): that is the CALLER's job (see
    click_mode.handle_click_to_add_structure()).

    Returns the list of newly-created Atom objects, in the structure's own
    file order. Raises ValueError if `structure` has no atoms (should never
    happen for a structure that passed structure_library.load_structure()'s
    own validation, checked defensively since this is reachable directly
    from a user click). """
    struct_atoms = structure["atoms"]
    struct_bonds = structure["bonds"]
    if not struct_atoms:
        raise ValueError ( "add_structure_at_position: structure '{}' has no atoms.".format (
                structure.get ( "name", "?" ) ) )

    centroid = np.mean ( np.array ( [ ( ax, ay, az ) for ( _sym, ax, ay, az ) in struct_atoms ], dtype = np.float64 ), axis = 0 )
    target   = np.array ( [ x, y, z ], dtype = np.float64 )

    local_index_to_new_id = { }
    new_atoms = [ ]
    for local_index, ( symbol, ax, ay, az ) in enumerate ( struct_atoms ):
        local_pos = np.array ( [ ax, ay, az ], dtype = np.float64 ) - centroid
        world_pos = target + local_pos
        new_atom = add_atom ( vismol_object, symbol = symbol,
                               x = float ( world_pos[0] ), y = float ( world_pos[1] ), z = float ( world_pos[2] ),
                               recompute_bonds = False, update_representation = False )
        local_index_to_new_id[local_index] = new_atom.atom_id
        new_atoms.append ( new_atom )

    if not hasattr ( vismol_object, "manual_bonds" ) or vismol_object.manual_bonds is None:
        vismol_object.manual_bonds = set ( )
    if not hasattr ( vismol_object, "manual_bond_orders" ) or vismol_object.manual_bond_orders is None:
        vismol_object.manual_bond_orders = { }
    if not hasattr ( vismol_object, "manual_aromatic_bonds" ) or vismol_object.manual_aromatic_bonds is None:
        vismol_object.manual_aromatic_bonds = set ( )

    for ( i, j, order, is_aromatic ) in struct_bonds:
        pair = ( min ( local_index_to_new_id[i], local_index_to_new_id[j] ),
                 max ( local_index_to_new_id[i], local_index_to_new_id[j] ) )
        vismol_object.manual_bonds.add ( pair )
        vismol_object.manual_bond_orders[pair] = order
        if is_aromatic:
            vismol_object.manual_aromatic_bonds.add ( pair )

    _reapply_manual_bonds ( vismol_object )

    vismol_object.create_representation ( rep_type = "sticks" )
    vismol_object.create_representation ( rep_type = "stick_spheres" )
    _activate_new_sphere_representation ( vismol_object, "stick_spheres" )
    vismol_object.create_representation ( rep_type = "nonbonded" )
    vismol_object.core_representations["picking_dots"] = None
    vismol_object.core_representations["picking_text"] = None

    # [EN] Atom count just changed -- see _refresh_bond_dependent_
    # representations()'s own (extended) docstring: any OTHER
    # representation this object carries (e.g. 'cartoon', inherited by an
    # "Edit in Builder" clone of a normal system) needs to either be
    # resized or deactivated here too, or its own stale self.indexes
    # crashes the next picking pass.
    _refresh_bond_dependent_representations ( vismol_object )

    vm_session = vismol_object.vm_session
    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.queue_draw ( )

    return new_atoms


# ============================================================================
#  Dihedral rotation tool + "edit an existing system" bootstrap -- see the
#  plan file for the full design (`click_mode.handle_click_to_pick_dihedral_
#  atom()` for the click-picking side, `empty_object.begin_editing_existing_
#  system()` for the other caller of bootstrap_manual_bonds_from_existing()
#  below).
# ============================================================================

def bootstrap_manual_bonds_from_existing ( vismol_object, source_vobject = None ):
    """ [EN] Populates manual_bonds/manual_bond_orders/manual_aromatic_
    bonds (the Builder's own explicit bond bookkeeping -- see add_atom()'s
    docstring for why every Builder mutation is driven from these instead
    of distance-based auto-detection) from an object that was NOT created
    by the Builder -- e.g. a real file-loaded molecule being opened via
    "Edit in Builder" for the first time. Without this, the first
    structural edit's own _reapply_manual_bonds() call would see an EMPTY
    manual_bonds and silently drop every real bond the molecule already
    had.

    source_vobject (new parameter -- BUG FIX, see below): if given, bonds
    are read from ITS state instead of vismol_object's own. empty_object.
    begin_editing_existing_system() passes the ORIGINAL vobject here, not
    the freshly-cloned temp vobject being edited -- confirmed by live user
    testing that using the clone's OWN bonds was silently dropping real
    bonds. Root cause: clone_system() clones the pDynamo SYSTEM, then
    _build_vobject_from_pdynamo_system() (pDynamo2EasyHybrid/session.py)
    rebuilds a BRAND NEW VismolObject from it by reading the System's own
    mmState 'Harmonic Bond' term -- or, if that term/mmState doesn't
    exist (common for a system that was loaded but never had a force
    field set up), falling back to a purely DISTANCE-based
    find_bonded_and_nonbonded_atoms() redetection. Neither path is
    guaranteed to reproduce whatever the ORIGINAL vismol_object's own
    .bonds already correctly showed (which is the real ground truth the
    user is editing) -- the clone's bond set can silently diverge from
    it. Reading from source_vobject instead sidesteps the whole System-
    roundtrip and uses the ACTUAL displayed topology directly. Atom ids
    correspond 1:1 between a clone_system() clone and the object it was
    cloned from (same atom order, deepcopied from the same underlying
    System), so no remapping is needed.

    Prefers source_vobject.manual_bonds (if already non-empty -- i.e. a
    PREVIOUS "Edit in Builder" round already recorded explicit bond
    orders/aromaticity onto it, see empty_object.
    finish_editing_existing_system()) over its raw .bonds, so those
    earlier edits survive a second round; falls back to source_vobject.
    bonds (a dict keyed by a normalised (min, max) atom-id tuple -> Bond
    object) otherwise -- CONNECTIVITY only (WHICH atoms are bonded) is
    taken from that dict; the ORDER is always re-derived fresh via
    vismol.core.bond_order_perception.perceive_bond_orders() (see BUG
    FIX below), never read from Bond.bond_order directly.

    [EN] BUG FIX (found by live user testing -- "carbonos que sao SP2
    estao com geometria de SP3 quando associados ao DYFF"): this used to
    default EVERY bond in the fallback branch to a plain single order
    (Bond.bond_order defaults to 1 on the class, and for the DOMINANT
    real-world case -- a system loaded through pDynamo2EasyHybrid/
    session.py's _build_vobject_from_pdynamo_system(), used by BOTH
    normal file loading and .easy session loading, see io_data.py -- that
    default is genuinely ALL that's stored there; nothing upstream ever
    computes a real per-bond order for that path, confirmed by reading
    it). Feeding every carbonyl/peptide-bond/aromatic-ring bond into
    DYFF as a plain single bond doesn't just leave a FEW bonds wrong: it
    silently turns EVERY sp2/sp carbon (and N/O) in the whole molecule
    tetrahedral in DetermineAtomGeometry's eyes, confirmed live with an
    N-methylacetamide test molecule (a stand-in for a protein backbone's
    own peptide bond) -- with orders lost this way, the carbonyl C
    (should be "C:Tri"), its O, and even the amide N all came out
    "*:Tet" instead; going through perceive_bond_orders() first (see
    below) restored the carbonyl C correctly to "C:Tri" in the same test.

    Fixed by reusing vismol.core.bond_order_perception.
    perceive_bond_orders(elements, bond_pairs) -- an already-proven,
    already-shipped Wang & Case-style (Antechamber 'bondtype') valence-
    penalty search, purely from element symbols + WHICH atoms are
    bonded, no 3D geometry needed (already used elsewhere in this
    codebase's own vismol_object.py for exactly this "no reliable
    external bond orders" situation -- e.g. loading a plain XYZ) --
    instead of defaulting every order to 1. Handles rings/conjugation
    (its own docstring cites a 36-atom, 6-fused-ring coronene test) and
    common charged functional groups (carboxylate, nitro, sulfonium, ...)
    via its own penalty model, so this ALSO makes formal-charge-sensitive
    DYFF patterns (e.g. "Guanidinium Cation", see [[project_dyff_force_
    field_assignment]]) more likely to match correctly, though this
    function still only sets bond orders, never atom.formalCharge itself
    -- no aromaticity flag is set either: pDynamo3's own
    ConvertInputConnectivity()/DetermineAromaticity() (already run inside
    empty_object._build_pdynamo_system_from_vismol_object()) correctly
    re-derives ring aromaticity from a properly Kekulized (alternating
    single/double) bond pattern on its own, no explicit flag needed here.

    No-op (returns without changing anything) if vismol_object already
    HAS a non-empty manual_bonds of its own. """
    existing = getattr ( vismol_object, "manual_bonds", None )
    if existing:
        return

    vismol_object.manual_bonds          = set ( )
    vismol_object.manual_bond_orders    = { }
    vismol_object.manual_aromatic_bonds = set ( )

    reference = source_vobject if source_vobject is not None else vismol_object

    reference_manual_bonds = getattr ( reference, "manual_bonds", None )
    if reference_manual_bonds:
        vismol_object.manual_bonds          = set ( reference_manual_bonds )
        vismol_object.manual_bond_orders    = dict ( getattr ( reference, "manual_bond_orders", None ) or { } )
        vismol_object.manual_aromatic_bonds = set ( getattr ( reference, "manual_aromatic_bonds", None ) or set ( ) )
        return

    bonds = getattr ( reference, "bonds", None ) or { }
    if not bonds:
        return

    pairs = list ( bonds.keys ( ) )
    try:
        from vismol.core.bond_order_perception import perceive_bond_orders
        elements = [ None ] * ( max ( max ( i, j ) for ( i, j ) in pairs ) + 1 )
        for atom_id, atom in reference.atoms.items ( ):
            if atom_id < len ( elements ):
                elements[atom_id] = atom.symbol
        order_map, _tps = perceive_bond_orders ( elements, pairs )
    except Exception as exc:
        # [EN] Never let a perception failure (unknown element, disabled
        # module, ...) block editing entirely -- falls back to the OLD,
        # known-imperfect "everything single" behaviour rather than
        # crashing the whole "Edit in Builder" entry point.
        order_map = { }
        dprint ( "WARNING atom_ops.bootstrap_manual_bonds_from_existing: bond order "
                "perception failed ({}) -- every bond defaulted to a plain single "
                "order instead.".format ( exc ) )

    for pair in pairs:
        vismol_object.manual_bonds.add ( pair )
        vismol_object.manual_bond_orders[pair] = int ( order_map.get ( pair, 1 ) )


def compute_dihedral_rotation_subgroup ( vismol_object, atom2_id, atom3_id ):
    """ [EN] For a dihedral picked as atom1-atom2-atom3-atom4, decides
    which atoms move when the atom2-atom3 bond is rotated: every atom
    reachable from atom3 WITHOUT crossing back through atom2 (atom2's own
    side -- including atom1 -- stays fixed as the pivot). A plain BFS over
    vismol_object.bonds' adjacency (same (min,max)-keyed dict
    bootstrap_manual_bonds_from_existing() reads above).

    Returns the subgroup as a set of atom_ids (always includes atom3_id
    itself), or None if atom2_id turns out to be reachable from atom3's
    side through some OTHER path too -- i.e. the atom2-atom3 bond is part
    of a ring, so "rotate only one side" isn't a well-defined operation
    (the two sides are still connected via the rest of the ring after the
    rotation, which would tear the structure apart). Callers must check
    for None and reject with a clear message instead of rotating. """
    adjacency = { }
    for ( i, j ) in vismol_object.bonds.keys ( ):
        adjacency.setdefault ( i, set ( ) ).add ( j )
        adjacency.setdefault ( j, set ( ) ).add ( i )

    visited  = { atom3_id }
    frontier = [ atom3_id ]
    is_ring  = False

    while frontier:
        current = frontier.pop ( )
        for neighbor in adjacency.get ( current, ( ) ):
            if current == atom3_id and neighbor == atom2_id:
                continue  # the pivot bond itself -- never cross it
            if neighbor == atom2_id:
                is_ring = True
                continue
            if neighbor not in visited:
                visited.add ( neighbor )
                frontier.append ( neighbor )

    return None if is_ring else visited


def rotate_atoms_around_bond ( vismol_object, pivot_atom_id, axis_atom_id, subgroup_atom_ids,
                                angle_rad, frame = 0, update_representation = True ):
    """ [EN] Rigidly rotates every atom in subgroup_atom_ids by angle_rad
    around the axis running from pivot_atom_id to axis_atom_id (anchored
    at pivot_atom_id's own position, which does not move). Reuses
    _rotate_vector_rodrigues() (same rotation primitive
    attach_fragment_at_hydrogen()/clean_up_structure() already use) rather
    than any of the OTHER rotate_bond()/rotate_dihedral() implementations
    already in this codebase (util/geometric_analysis.py, util/topol.py,
    eSession.py) -- those operate on inconsistent conventions (global
    picking-selection state, a foreign list-based molecule representation,
    or are confirmed dead/broken code) that don't match how the Builder
    itself always addresses vismol_object.frames[frame, atom_id] directly.

    update_representation mirrors move_atom()'s own flag/behaviour exactly
    (mark every active representation's coordinates dirty for re-upload,
    not a full rebuild -- rotation never changes atom count or bonds). """
    pivot_pos = np.array ( vismol_object.frames[frame, pivot_atom_id], dtype = np.float64 )
    axis_pos  = np.array ( vismol_object.frames[frame, axis_atom_id],  dtype = np.float64 )
    axis      = axis_pos - pivot_pos
    axis_norm = float ( np.linalg.norm ( axis ) )
    if axis_norm < 1e-6:
        raise ValueError ( "rotate_atoms_around_bond: degenerate (zero-length) rotation axis." )
    axis = axis / axis_norm

    for atom_id in subgroup_atom_ids:
        local_pos = np.array ( vismol_object.frames[frame, atom_id], dtype = np.float64 ) - pivot_pos
        new_pos   = _rotate_vector_rodrigues ( local_pos, axis, angle_rad ) + pivot_pos
        vismol_object.frames[frame, atom_id] = new_pos.astype ( np.float32 )

    vismol_object.mass_center = np.mean ( vismol_object.frames[0], axis = 0 )

    if update_representation:
        vm_session = vismol_object.vm_session
        for rep in vismol_object.representations.values ( ):
            if rep is not None and rep.active:
                rep.was_rep_coord_modified = True
                rep.was_sel_coord_modified = True
        picking_dots = vismol_object.core_representations.get ( "picking_dots" )
        if picking_dots is not None:
            picking_dots.was_rep_coord_modified = True
        if getattr ( vm_session, "vm_glcore", None ) is not None:
            vm_session.vm_glcore.updated_coords = True
            vm_session.vm_glcore.queue_draw ( )


# ============================================================================
#  Transform Selection tool (user's own request: "transladar e rotacionar
#  (x, y, z) uma selecao de atomos, usando o viewing selection") -- see
#  gui/windows/builder/transform_selection_window.py for the live-slider UI
#  this feeds. Unlike rotate_atoms_around_bond() above (dihedral rotation,
#  axis defined by 2 SPECIFIC atoms, one of them the fixed pivot), this
#  operates on an ARBITRARY atom_ids set with no inherent axis/pivot of its
#  own -- translation is a plain per-axis offset, rotation is anchored at
#  the selection's own CENTROID (selection_centroid() below, recomputed
#  FRESH every call, never cached -- stays correct across interleaved
#  translate/rotate gestures) around the vismol_object's own LOCAL X/Y/Z
#  axes (same local space vismol_object.frames itself already lives in,
#  not world/view space -- consistent with every other Builder transform
#  operating directly on frames[frame, atom_id]).
# ============================================================================

def selection_centroid ( vismol_object, atom_ids, frame = 0 ):
    """ Mean CURRENT position of atom_ids at `frame` -- the rotation pivot
    transform_selection_window.py's rotate sliders use, recomputed fresh on
    every rotation step (not captured once) so it stays correct even after
    an earlier translate/rotate step already moved the group. Raises
    ValueError for an empty atom_ids (nothing to centre on). """
    positions = [ vismol_object.frames[frame, atom_id] for atom_id in atom_ids ]
    if not positions:
        raise ValueError ( "selection_centroid: atom_ids is empty." )
    return np.mean ( np.array ( positions, dtype = np.float64 ), axis = 0 )


def translate_atoms ( vismol_object, atom_ids, dx, dy, dz, frame = 0, update_representation = True ):
    """ Rigidly translates every atom in atom_ids by (dx, dy, dz) -- plain
    per-axis offset, no pivot/rotation involved. update_representation
    mirrors move_atom()/rotate_atoms_around_bond()'s own flag exactly (mark
    every active representation's coordinates dirty for re-upload, not a
    full rebuild -- translation never changes atom count or bonds). """
    delta = np.array ( [ dx, dy, dz ], dtype = np.float64 )
    for atom_id in atom_ids:
        current = np.array ( vismol_object.frames[frame, atom_id], dtype = np.float64 )
        vismol_object.frames[frame, atom_id] = ( current + delta ).astype ( np.float32 )

    vismol_object.mass_center = np.mean ( vismol_object.frames[0], axis = 0 )

    if update_representation:
        vm_session = vismol_object.vm_session
        for rep in vismol_object.representations.values ( ):
            if rep is not None and rep.active:
                rep.was_rep_coord_modified = True
                rep.was_sel_coord_modified = True
        picking_dots = vismol_object.core_representations.get ( "picking_dots" )
        if picking_dots is not None:
            picking_dots.was_rep_coord_modified = True
        if getattr ( vm_session, "vm_glcore", None ) is not None:
            vm_session.vm_glcore.updated_coords = True
            vm_session.vm_glcore.queue_draw ( )


def rotate_atoms_around_point ( vismol_object, atom_ids, pivot, axis, angle_rad,
                                 frame = 0, update_representation = True ):
    """ Rigidly rotates every atom in atom_ids by angle_rad around `axis`
    (any nonzero vector -- normalized here, so (1,0,0)/(0,1,0)/(0,0,1) for
    the object's own local X/Y/Z axes work directly), anchored at `pivot` --
    a plain (x, y, z) point, NOT necessarily another atom's position (unlike
    rotate_atoms_around_bond(), which is always anchored at a specific
    pivot ATOM). Reuses the same _rotate_vector_rodrigues() primitive every
    other Builder rotation already shares. """
    axis = np.array ( axis, dtype = np.float64 )
    axis_norm = float ( np.linalg.norm ( axis ) )
    if axis_norm < 1e-9:
        raise ValueError ( "rotate_atoms_around_point: degenerate (zero-length) rotation axis." )
    axis = axis / axis_norm
    pivot = np.array ( pivot, dtype = np.float64 )

    for atom_id in atom_ids:
        local_pos = np.array ( vismol_object.frames[frame, atom_id], dtype = np.float64 ) - pivot
        new_pos   = _rotate_vector_rodrigues ( local_pos, axis, angle_rad ) + pivot
        vismol_object.frames[frame, atom_id] = new_pos.astype ( np.float32 )

    vismol_object.mass_center = np.mean ( vismol_object.frames[0], axis = 0 )

    if update_representation:
        vm_session = vismol_object.vm_session
        for rep in vismol_object.representations.values ( ):
            if rep is not None and rep.active:
                rep.was_rep_coord_modified = True
                rep.was_sel_coord_modified = True
        picking_dots = vismol_object.core_representations.get ( "picking_dots" )
        if picking_dots is not None:
            picking_dots.was_rep_coord_modified = True
        if getattr ( vm_session, "vm_glcore", None ) is not None:
            vm_session.vm_glcore.updated_coords = True
            vm_session.vm_glcore.queue_draw ( )


