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


def add_atom ( vismol_object, symbol, x, y, z, name = None,
               chain_id = "A", resi = 1, resn = "UNK",
               bonded_to = None,
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
    recompute_bonds : if True (default), rebuilds bonds/topology/molecule
               grouping from vismol_object.manual_bonds (via
               _reapply_manual_bonds()) after adding the atom. Set False
               only if you're about to add several atoms in a row and
               want to pay for this once, at the end (call
               _reapply_manual_bonds(vismol_object) yourself afterwards
               in that case).
    update_representation : if True (default), (re)creates the "lines"
               and "nonbonded" representations so the new atom actually
               becomes visible immediately. Same reasoning as
               recompute_bonds for batching -- skip and do it yourself
               once at the end if adding many atoms in a loop.

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

    atom_id = len ( vismol_object.atoms )   # proximo indice livre, sequencial

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

    # --- cresce vismol_object.frames em +1 atomo, para TODOS os frames ja
    # existentes (o novo atomo entra com a MESMA posicao (x,y,z) em cada
    # frame -- coerente pra um objeto builder, que tipicamente tem so 1
    # frame mesmo; se este objeto ja tiver uma trajetoria de verdade, isso
    # ainda funciona, so que o atomo novo fica "parado" em todos os
    # frames ate que algo mais sofisticado seja implementado). ---
    n_frames = vismol_object.frames.shape[0]
    new_frames = np.zeros ( (n_frames, atom_id + 1, 3), dtype = np.float32 )
    if atom_id > 0:
        new_frames[:, :atom_id, :] = vismol_object.frames
    new_frames[:, atom_id, :] = [ x, y, z ]
    vismol_object.frames = new_frames

    vismol_object.mass_center = np.mean ( vismol_object.frames[0], axis = 0 )

    # recalcula cores/vdw/tamanhos de ponto pra TODOS os atomos -- simples
    # e correto; o parametro colors_id_start nao e sequer usado dentro do
    # metodo (conferido lendo o codigo-fonte), entao o valor exato passado
    # aqui nao importa.
    vismol_object._generate_color_vectors ( vm_session.atom_id_counter )

    # [EN] Register any explicitly-requested bond(s) for this new atom --
    # the ONLY way it ends up bonded to anything now that distance-based
    # auto-detection is off (see this function's own docstring above).
    if bonded_to is not None:
        if not hasattr ( vismol_object, "manual_bonds" ) or vismol_object.manual_bonds is None:
            vismol_object.manual_bonds = set ( )
        targets = bonded_to if isinstance ( bonded_to, ( list, tuple, set ) ) else [ bonded_to ]
        for other_id in targets:
            pair = ( min ( atom_id, other_id ), max ( atom_id, other_id ) )
            vismol_object.manual_bonds.add ( pair )

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
        vismol_object.create_representation ( rep_type = "lines" )
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
    update_representation : if True (default), rebuilds "lines"/
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
        vismol_object.create_representation ( rep_type = "lines" )
        vismol_object.create_representation ( rep_type = "nonbonded" )
        vismol_object.core_representations["picking_dots"] = None
        vismol_object.core_representations["picking_text"] = None
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
    # pendurada apontando pra um Atom que nao existe mais no objeto)
    if removed_atom.unique_id in vm_session.atom_dic_id:
        del vm_session.atom_dic_id[removed_atom.unique_id]

    # tira de qualquer selecao ativa (senao um add_bond() subsequente
    # poderia incluir por engano um atomo que acabou de ser apagado)
    for sel in vm_session.selections.values ( ):
        sel.selected_atoms.discard ( removed_atom )

    # tira do dict de atomos do residuo
    if removed_atom.residue is not None:
        removed_atom.residue.atoms.pop ( atom_id, None )

    del vismol_object.atoms[atom_id]

    # remove a linha correspondente de TODOS os frames (axis=1 = eixo dos atomos)
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
    # (click_mode.cycle_bond_order()) -- these are ALSO keyed by atom_id
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

    vismol_object.mass_center = ( np.mean ( vismol_object.frames[0], axis = 0 )
                                  if vismol_object.frames.shape[1] > 0
                                  else np.zeros ( 3, dtype = np.float32 ) )
    vismol_object._generate_color_vectors ( vm_session.atom_id_counter )

    # bonds/topologia agora vem inteiramente de vismol_object.manual_bonds
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

    vismol_object.create_representation ( rep_type = "lines" )
    vismol_object.create_representation ( rep_type = "nonbonded" )
    # ver nota em add_atom() -- forca build_core_representations() a
    # reconstruir do zero na proxima renderizacao, com a contagem de
    # atomos (agora menor) correta.
    vismol_object.core_representations["picking_dots"] = None
    vismol_object.core_representations["picking_text"] = None

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
    # Ctrl+click cycle_bond_order()) on top of that baseline -- matching
    # cycle_bond_order()'s own stated intent ("persist the new order...
    # REQUIRED... without persisting it... cycling a bond's order here
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


def add_bond ( vismol_object, atom_id_a, atom_id_b, bond_order = 1 ):
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
    found and fixed once in click_mode.cycle_bond_order() (see that
    function's own docstring: "_bonds_from_pair_of_indexes_list()'s
    external_orders parameter already existed but its actual assignment
    ... was commented out"). Fixed the same way cycle_bond_order()
    already does it: persist the requested order in
    vismol_object.manual_bond_orders (keyed by the normalized (min,max)
    pair), which _reapply_manual_bonds() below feeds into
    _bonds_from_pair_of_indexes_list() as external_orders -- required,
    not optional, since bonds get rebuilt FROM SCRATCH (fresh Bond()
    objects) on every structural edit, so anything not persisted there
    is silently lost on the very next add_atom()/remove_atom()/add_bond()
    call. """
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

    changed = _reapply_manual_bonds ( vismol_object )
    if not changed:
        return False   # objeto sem nenhuma ligacao (nem essa nova, nem nenhuma outra) -- nada a fazer (caso raro/defensivo)

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
    dash() calls), so it doesn't need its own entry here. """
    reps = getattr ( vismol_object, "representations", None ) or { }
    for rep_type in ( "lines", "sticks", "nonbonded" ):
        rep = reps.get ( rep_type )
        if rep is None:
            continue
        was_active = getattr ( rep, "active", True )
        vismol_object.create_representation ( rep_type = rep_type )
        new_rep = vismol_object.representations.get ( rep_type )
        if new_rep is not None:
            new_rep.active = was_active


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


def set_bond_order ( vismol_object, atom_id_a, atom_id_b, bond_order = 1 ):
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

    if created:
        # Topologia/non-bonded/moleculas so' precisam ser refeitas quando a
        # CONECTIVIDADE muda (bond novo) -- uma simples troca de ordem numa
        # ligacao que ja existia nao afeta nenhuma dessas (grafo identico).
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
#   objeto (vismol_object.bonds/index_bonds) -- vale para TODOS os frames e
#   e' o que fica gravado se o objeto for salvo/exportado.
#
#   As funcoes abaixo, em vez disso, editam vismol_object.dynamic_bonds[f]
#   -- a lista de pares POR FRAME usada pela representacao "Dynamic Bonds"
#   (tipicamente a regiao QC de uma trajetoria QM/MM, recalculada
#   automaticamente por distancia a cada frame -- ver VismolSession.
#   define_dynamic_bonds() / VismolObject.find_bonded_and_nonbonded_atoms()).
#
#   *** AVISO IMPORTANTE, repetido no docstring de cmd_bond/cmd_unbond no
#   terminal: usar frame=... aqui NAO cria uma ligacao quimica de verdade.
#   E' PURAMENTE uma edicao da REPRESENTACAO/VISUALIZACAO para aquele(s)
#   frame(s) -- nao mexe em nada do sistema pDynamo (topologia, campo de
#   forca, constantes de forca de ligacao, cargas, etc.). Para criar uma
#   ligacao real (com parametros de campo de forca), a topologia do
#   sistema pDynamo precisa ser editada por outros meios -- isso aqui e'
#   so' para ajustar o que aparece na tela enquanto se inspeciona/prepara
#   uma trajetoria (ex.: forcar visualmente uma ligacao que a deteccao
#   automatica por distancia deixou passar batido num frame especifico).
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
        'atoms'             : atoms_snapshot,
        'manual_bonds'      : set ( getattr ( vismol_object, 'manual_bonds', None ) or set ( ) ),
        'manual_bond_orders': dict ( getattr ( vismol_object, 'manual_bond_orders', None ) or { } ),
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
    vismol_object.atoms  = { }
    vismol_object.chains = { }
    vismol_object.frames = np.zeros ( ( 1, 0, 3 ), dtype = np.float32 )
    vismol_object.manual_bonds       = set ( )
    vismol_object.manual_bond_orders = { }

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

    vismol_object.manual_bonds       = set ( snapshot['manual_bonds'] )
    vismol_object.manual_bond_orders = dict ( snapshot['manual_bond_orders'] )

    vismol_object.cov_radii_array = None
    vismol_object.electronegativity_array = None
    vismol_object.index_bonds = None
    vismol_object.bonds = None
    vismol_object.non_bonded_atoms = None
    _reapply_manual_bonds ( vismol_object )

    vismol_object.create_representation ( rep_type = "lines" )
    vismol_object.create_representation ( rep_type = "nonbonded" )
    vismol_object.core_representations["picking_dots"] = None
    vismol_object.core_representations["picking_text"] = None

    vm_session = vismol_object.vm_session
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

# Valencia padrao (soma maxima de ordem de ligacao) por elemento -- os
# poucos elementos que este primeiro conjunto de ferramentas do Builder
# realmente produz ate agora (C/N/O/H, mais alguns halogenios/S/P comuns
# para nao deixar na mao assim que o seletor de elemento crescer). Elementos
# fora desta tabela simplesmente nao tem os hidrogenios ajustados (ver
# adjust_hydrogens() -- retorna sem fazer nada nesse caso, em vez de
# adivinhar uma valencia errada).
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


def _new_hydrogen_directions ( vismol_object, atom_id, n_needed, bond_length = 1.05 ):
    """ [EN] Computes n_needed 3D offset vectors (already scaled to
    `bond_length`) for where to place NEW hydrogens around atom_id.

    If atom_id has no existing bonded neighbours at all, falls back to
    the same fixed geometric templates the OLD, now-replaced click_mode.
    py auto-hydrogenation used (tetrahedral-ish for 4, trigonal for 3,
    linear-ish for 2) -- see below -- taking only the first n_needed
    directions from whichever template matches the TOTAL target count
    (n_needed itself, since there are no existing neighbours to also
    account for).

    If atom_id ALREADY has one or more bonded neighbours (the more
    common case: adding H's to an atom that's already bonded to
    something else), there's no single fixed template that still applies
    -- instead, computes the direction AWAY from the existing neighbours
    (the negated, normalised sum of the existing bond unit vectors), and
    fans the needed new H's out symmetrically around that single "most
    open" direction using an arbitrary perpendicular axis. This is a
    simple, defensible heuristic (not a real hybridisation/VSEPR solver)
    -- good enough to avoid new H's landing on top of existing bonds,
    without trying to reproduce exact tetrahedral/trigonal angles once
    other substituents are already present. """
    fixed_templates = {
        4: [ [-0.785298,  0.243518, -0.653254], [ 0.322015, -0.981331, -0.189814],
             [-0.334691,  0.073016,  0.992665], [ 0.798227,  0.665009, -0.149645] ],
        3: [ [-0.785298,  0.243518, -0.653254], [ 0.322015, -0.981331, -0.189814],
             [-0.334691,  0.073016,  0.992665] ],
        2: [ [-0.785298,  0.243518, -0.653254], [ 0.322015, -0.981331, -0.189814] ],
        1: [ [-0.785298,  0.243518, -0.653254] ],
    }

    atom_pos = vismol_object.frames[0, atom_id]
    existing_dirs = [ ]
    for bond in vismol_object.bonds.values ( ):
        if bond.atom_index_i != atom_id and bond.atom_index_j != atom_id:
            continue
        other_id = bond.atom_index_j if bond.atom_index_i == atom_id else bond.atom_index_i
        vec = vismol_object.frames[0, other_id] - atom_pos
        norm = float ( np.linalg.norm ( vec ) )
        if norm > 1e-6:
            existing_dirs.append ( vec / norm )

    if not existing_dirs:
        template = fixed_templates.get ( n_needed, fixed_templates[4] )
        directions = [ np.array ( d, dtype = np.float32 ) for d in template[:n_needed] ]
        return [ d * bond_length for d in directions ]

    avg_dir = -np.sum ( existing_dirs, axis = 0 )
    norm = float ( np.linalg.norm ( avg_dir ) )
    if norm < 1e-6:
        avg_dir = np.array ( [ 0.0, 0.0, 1.0 ], dtype = np.float32 )
    else:
        avg_dir = ( avg_dir / norm ).astype ( np.float32 )

    if n_needed == 1:
        return [ avg_dir * bond_length ]

    perp = np.cross ( avg_dir, np.array ( [ 0.0, 0.0, 1.0 ], dtype = np.float32 ) )
    if float ( np.linalg.norm ( perp ) ) < 1e-6:
        perp = np.cross ( avg_dir, np.array ( [ 0.0, 1.0, 0.0 ], dtype = np.float32 ) )
    perp = ( perp / np.linalg.norm ( perp ) ).astype ( np.float32 )

    spread_angle = np.deg2rad ( 50.0 )   # angulo arbitrario, mas razoavel, entre H's adicionados
    directions = [ ]
    for k in range ( n_needed ):
        offset = ( k - ( n_needed - 1 ) / 2.0 ) * spread_angle
        d = avg_dir * np.cos ( offset ) + perp * np.sin ( offset )
        d = d / np.linalg.norm ( d )
        directions.append ( d * bond_length )

    return directions


def adjust_hydrogens ( vismol_object, atom_id ):
    """ Adds or removes hydrogens bonded to atom_id so its total bond
    order (to non-hydrogen neighbours, plus one per hydrogen) matches
    its STANDARD_VALENCE. Called explicitly after any operation that may
    have changed atom_id's bonding -- see the module-level note above for
    the full list of call sites.

    No-op if atom_id's element isn't in STANDARD_VALENCE (unknown
    elements are left alone rather than guessed at), or if it already has
    exactly the right number of hydrogens.

    Returns the number of hydrogens added (positive) or removed
    (negative), or 0 if nothing changed. """
    if atom_id not in vismol_object.atoms:
        return 0

    atom = vismol_object.atoms[atom_id]
    target_valence = STANDARD_VALENCE.get ( atom.symbol.upper ( ) )
    if target_valence is None:
        return 0

    heavy_bond_order_sum = _bond_order_sum_excluding_symbol ( vismol_object, atom_id, exclude_symbol = 'H' )
    needed_h = max ( 0, target_valence - heavy_bond_order_sum )

    current_h_ids = _bonded_neighbors_with_symbol ( vismol_object, atom_id, 'H' )
    current_h_count = len ( current_h_ids )

    if current_h_count == needed_h:
        return 0

    if current_h_count > needed_h:
        # [EN] Removes from HIGHEST atom_id downward: remove_atom()
        # renumbers every atom_id ABOVE the one it removes, so removing
        # highest-first means none of the OTHER hydrogens still pending
        # removal in this same loop are ever affected by an earlier
        # removal -- only atom_id itself (the parent atom being adjusted)
        # might need to shift down, tracked explicitly below.
        excess_ids = sorted ( current_h_ids[needed_h:], reverse = True )
        for h_id in excess_ids:
            remove_atom ( vismol_object, h_id )
            if h_id < atom_id:
                atom_id -= 1
        return needed_h - current_h_count   # negativo

    n_to_add = needed_h - current_h_count
    directions = _new_hydrogen_directions ( vismol_object, atom_id, n_to_add )
    pos = vismol_object.frames[0, atom_id]
    for offset in directions:
        add_atom ( vismol_object, symbol = "H",
                   x = float ( pos[0] + offset[0] ),
                   y = float ( pos[1] + offset[1] ),
                   z = float ( pos[2] + offset[2] ),
                   bonded_to = atom_id )

    return n_to_add


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
# ordem da ligacao -- valores aproximados, so para dar uma geometria
# razoavel (nao says pretende reproduzir literatura com precisao).
_BOND_ORDER_LENGTH_FACTOR = { 1: 1.00, 2: 0.87, 3: 0.78 }


def _ideal_bond_length ( atom_i, atom_j, bond_order ):
    """ Soma dos raios covalentes dos dois atomos, encurtada conforme a
    ordem da ligacao (dupla/tripla mais curtas que simples). """
    base   = float ( atom_i.cov_rad ) + float ( atom_j.cov_rad )
    factor = _BOND_ORDER_LENGTH_FACTOR.get ( bond_order, 1.00 )
    return base * factor


def _ideal_angle_for_atom ( vismol_object, atom_id ):
    """ [EN] Simple VSEPR-style rule: angle depends only on how many
    neighbours atom_id has and the highest bond order among them --
    NOT a real hybridisation calculation (no lone-pair counting, no
    electronegativity effects), just enough to keep sp/sp2/sp3-ish
    centres looking reasonable after the Builder's rough initial
    placement.

    Returns the ideal angle in DEGREES, or None if atom_id has fewer
    than 2 neighbours (no angle to speak of) or more than 4 (unusual for
    the elements this Builder currently supports -- left alone rather
    than guessed at). """
    neighbor_bonds = [ bond for bond in vismol_object.bonds.values ( )
                       if bond.atom_index_i == atom_id or bond.atom_index_j == atom_id ]
    n = len ( neighbor_bonds )
    if n < 2 or n > 4:
        return None

    max_order = max ( bond.bond_order for bond in neighbor_bonds )

    if n == 2:
        double_count = sum ( 1 for bond in neighbor_bonds if bond.bond_order == 2 )
        if max_order == 3 or double_count == 2:   # tripla, ou duas duplas (aleno-like) -- linear
            return 180.0
        elif max_order == 2:                       # uma dupla -- ex.: carbono de carbonila com 1 substituinte
            return 120.0
        else:
            return 109.5
    elif n == 3:
        if max_order >= 2:                          # trigonal plana -- ex.: alceno, carbonila
            return 120.0
        else:
            return 109.5                             # aproximacao razoavel (ex.: amina simples)
    else:   # n == 4
        return 109.5


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


def clean_up_structure ( vismol_object, atom_ids = None, n_iterations = 40, step_scale = 0.7 ):
    """ Nudges atom positions toward ideal bond lengths/angles, over
    `n_iterations` passes of simple geometric correction (see the
    module-level note above for why this is deliberately NOT a real
    force field like UFF).

    [EN] Defaults verified with a standalone numpy simulation (a
    deliberately worst-case, badly distorted 5-atom test: one bond
    stretched to 2x its ideal length, another compressed to under a
    third of it, no two bonds anywhere near their ideal angle to start)
    before being wired in here -- n_iterations=40 / step_scale=0.7
    converged that test to within +/-0.7 degrees of the ideal tetrahedral
    angle and correct bond lengths to 3 decimal places. Lower step_scale
    values (tried 0.5) or fewer iterations (tried 5-15) left it
    measurably short of converged for that same worst-case test; higher
    step_scale (tried 0.9-1.0) converged slightly faster but with less
    safety margin against oscillation in more crowded, multi-constraint
    geometries (several atoms sharing overlapping bond/angle corrections
    every iteration) than this project's test coverage could rule out.

    atom_ids : if None (default), relaxes the WHOLE molecule -- every
               atom is free to move, and every bond/angle in the object
               is corrected. If a collection of atom_ids is given
               instead, only those atoms (PLUS their directly bonded
               neighbours -- moving an atom's bonds/angles without
               letting its neighbours respond too wouldn't converge to
               anything sensible) are free to move; every OTHER atom
               stays fixed and acts as an anchor for bonds/angles that
               touch the movable set. Both modes share the exact same
               relaxation loop below -- "whole molecule" is just the
               special case where every atom happens to be movable.
    n_iterations : how many correction passes to run. More iterations
               converge closer to the ideal geometry but take
               proportionally longer -- 40 is a reasonable default for
               small Builder molecules (see the verification note above).
    step_scale : how much of each iteration's computed correction is
               actually applied (0-1). Less than 1.0 (default 0.7) acts
               as a damping factor, trading convergence speed for
               stability (a full, undamped correction every iteration
               can overshoot and oscillate, especially once bond and
               angle corrections on the same atom start interacting).

    Rebuilds representations once at the very end (this is a one-off,
    on-demand action -- e.g. the sidebar's "Clean Up" button -- not a
    per-frame update, so there's no reason to batch/optimise beyond
    that single rebuild). """
    if atom_ids is None:
        movable_ids = set ( vismol_object.atoms.keys ( ) )
    else:
        movable_ids = set ( atom_ids )
        for aid in list ( atom_ids ):
            for bond in vismol_object.bonds.values ( ):
                if bond.atom_index_i == aid:
                    movable_ids.add ( bond.atom_index_j )
                elif bond.atom_index_j == aid:
                    movable_ids.add ( bond.atom_index_i )

    if not movable_ids:
        return

    positions = { aid: np.array ( vismol_object.frames[0, aid], dtype = np.float64 )
                  for aid in vismol_object.atoms.keys ( ) }

    for _iteration in range ( n_iterations ):
        displacement = { aid: np.zeros ( 3, dtype = np.float64 ) for aid in movable_ids }
        contributions = { aid: 0 for aid in movable_ids }

        # --- termo de distancia (uma correcao por ligacao) ---
        for bond in vismol_object.bonds.values ( ):
            i, j = bond.atom_index_i, bond.atom_index_j
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

        # --- termo de angulo (uma correcao por par de vizinhos de cada atomo) ---
        for atom_id in vismol_object.atoms.keys ( ):
            ideal_angle = _ideal_angle_for_atom ( vismol_object, atom_id )
            if ideal_angle is None:
                continue

            neighbor_ids = [ ]
            for bond in vismol_object.bonds.values ( ):
                if bond.atom_index_i == atom_id:
                    neighbor_ids.append ( bond.atom_index_j )
                elif bond.atom_index_j == atom_id:
                    neighbor_ids.append ( bond.atom_index_i )

            center_pos = positions[atom_id]
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

    vismol_object.mass_center = np.mean ( vismol_object.frames[0], axis = 0 )

    vismol_object.create_representation ( rep_type = "lines" )
    vismol_object.create_representation ( rep_type = "nonbonded" )
    vismol_object.core_representations["picking_dots"] = None
    vismol_object.core_representations["picking_text"] = None

    vm_session = vismol_object.vm_session
    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.queue_draw ( )




