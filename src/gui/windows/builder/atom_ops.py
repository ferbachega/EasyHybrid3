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
    recompute_bonds : if True (default), re-runs
               vismol_object.find_bonded_and_nonbonded_atoms() after
               adding the atom, so bonds/non-bonded status stay correct
               for the WHOLE object (not just the new atom). Set False
               if you are about to add several atoms in a row and only
               want to pay for bond detection once, at the end (call
               vismol_object.find_bonded_and_nonbonded_atoms() yourself
               afterwards in that case).
    update_representation : if True (default), (re)creates the "lines"
               and "nonbonded" representations so the new atom actually
               becomes visible immediately. Same reasoning as
               recompute_bonds for batching -- skip and do it yourself
               once at the end if adding many atoms in a loop.

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

    if recompute_bonds:
        # forca recalculo com o novo tamanho, em vez de reusar arrays de
        # tamanho antigo (que dariam erro de shape ou dados errados).
        # Resetar TODOS os campos abaixo antes de chamar
        # find_bonded_and_nonbonded_atoms() de novo e necessario porque
        # os metodos internos que ela chama tem asserts explicitos
        # exigindo estado ainda-nao-inicializado (confirmado rodando de
        # verdade, nao so por leitura: a 2a chamada de add_atom() dava
        # "AssertionError" em _bonds_from_pair_of_indexes_list()
        # [assert self.bonds is None] e teria dado o mesmo em
        # _get_non_bonded_from_bonded_list() [assert self.non_bonded_atoms
        # is None] logo em seguida). index_bonds tambem e zerado so pra
        # nao disparar o aviso (nao um erro) "ja existe informacao de
        # contatos" que o metodo imprime quando index_bonds != None.
        vismol_object.cov_radii_array = None
        vismol_object.electronegativity_array = None
        vismol_object.index_bonds = None
        vismol_object.bonds = None
        vismol_object.non_bonded_atoms = None
        vismol_object.find_bonded_and_nonbonded_atoms ( )

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

    Bonds/topology: recomputed from scratch via the same
    find_bonded_and_nonbonded_atoms() used by add_atom() -- correct
    because bond detection is purely distance-based; the remaining
    atoms' positions haven't changed, so re-running detection on the
    smaller, renumbered atom set gives the right answer directly. KNOWN
    LIMITATION shared with add_bond() below: a manually-added bond is
    NOT preserved across a remove_atom() call, for the same reason it
    isn't guaranteed to survive add_atom() either -- see add_bond()'s
    docstring. """
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
    for old_id in sorted ( old_atoms.keys ( ) ):
        atom = old_atoms[old_id]
        new_id = old_id - 1 if old_id > atom_id else old_id
        if new_id != old_id:
            atom.atom_id = new_id
            if atom.residue is not None:
                atom.residue.atoms.pop ( old_id, None )
                atom.residue.atoms[new_id] = atom
        new_atoms[new_id] = atom
    vismol_object.atoms = new_atoms

    vismol_object.mass_center = ( np.mean ( vismol_object.frames[0], axis = 0 )
                                  if vismol_object.frames.shape[1] > 0
                                  else np.zeros ( 3, dtype = np.float32 ) )
    vismol_object._generate_color_vectors ( vm_session.atom_id_counter )

    # recalcula ligacoes do zero (deteccao por distancia -- ver nota no
    # docstring sobre bonds manuais nao sobreviverem a isso)
    vismol_object.cov_radii_array = None
    vismol_object.electronegativity_array = None
    vismol_object.index_bonds = None
    vismol_object.bonds = None
    vismol_object.non_bonded_atoms = None
    if len ( vismol_object.atoms ) > 0:
        vismol_object.find_bonded_and_nonbonded_atoms ( )
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


def add_bond ( vismol_object, atom_id_a, atom_id_b, bond_order = 1 ):
    """ [EN] Fifth building block of the Builder ('b' key -- add a bond
    between exactly two currently-selected atoms). Unlike the automatic,
    distance-based bond detection find_bonded_and_nonbonded_atoms() does
    (used internally by add_atom()/remove_atom() above), this creates an
    EXPLICIT bond regardless of the distance between the two atoms --
    necessary for cases automatic detection would never produce on its
    own (e.g. deliberately closing a ring across a gap the covalent-
    radius/distance heuristic wouldn't recognise as bonded).

    Implementation: appends the new pair to the EXISTING index_bonds
    (preserving whatever bonds -- auto-detected or previously manually
    added -- were already there) rather than re-running distance-based
    detection from scratch, then re-runs just the downstream bookkeeping
    steps (_bonds_from_pair_of_indexes_list, _get_non_bonded_from_bonded_list,
    _generate_topology_from_index_bonds, define_molecules,
    define_Calpha_backbone) that also run inside
    find_bonded_and_nonbonded_atoms() -- same reset-before-recompute
    pattern as add_atom() (self.bonds / self.non_bonded_atoms reset to
    None first, required by the asserts inside those methods -- see
    add_atom()'s own comment for where that was discovered).

    KNOWN LIMITATION: a manually-added bond is only preserved as long as
    nothing calls find_bonded_and_nonbonded_atoms() FROM SCRATCH
    afterwards -- and add_atom()/remove_atom() both do exactly that
    (needed for THEIR correctness, since a newly added/removed atom can
    change what's genuinely in bonding distance of everything else).
    Concretely: add a manual out-of-distance bond with add_bond(), then
    add or remove a different atom elsewhere in the same molecule, and
    the manual bond will be silently dropped (overwritten by fresh
    distance-based detection). Not fixed here -- would need a persistent
    "manual bonds" list on vismol_object that gets re-merged back in
    every time add_atom()/remove_atom() recompute from scratch; flagged
    as a follow-up rather than attempted in this pass, to keep this
    step's scope contained. """
    if atom_id_a == atom_id_b:
        raise ValueError ( "add_bond: cannot bond an atom to itself." )
    if atom_id_a not in vismol_object.atoms or atom_id_b not in vismol_object.atoms:
        raise ValueError ( "add_bond: atom_id_a={} or atom_id_b={} does not exist in this object.".format (
                            atom_id_a, atom_id_b ) )

    pair = ( min ( atom_id_a, atom_id_b ), max ( atom_id_a, atom_id_b ) )

    # [EN] BUG FIXED (caught by testing, not just reading): the auto-
    # detected index_bonds already stores each bond TWICE, once per
    # direction (e.g. a bond between atoms 0 and 1 appears as both
    # [0,1] and [1,0] in the flat list) -- confirmed live: an object with
    # one auto-detected bond had index_bonds == [0,1,1,0], four entries
    # for ONE bond. Converting each raw pair to a sorted tuple collapses
    # both directions to the SAME (0,1) -- but without deduplicating,
    # that duplicate carried straight through, and the pair being added
    # here got appended on top of it, producing an index_bonds with a
    # repeated pair (confirmed: a second bond ended up as
    # [0,1,0,1,0,2] instead of the correct [0,1,0,2]). Fixed by
    # deduplicating via a set before re-flattening.
    if vismol_object.index_bonds is not None and len ( vismol_object.index_bonds ) > 0:
        flat = np.asarray ( vismol_object.index_bonds ).reshape ( -1, 2 )
        existing_pairs = set ( tuple ( sorted ( p ) ) for p in flat.tolist ( ) )
    else:
        existing_pairs = set ( )

    if pair in existing_pairs:
        return False   # ja existe esse bond -- nada a fazer

    existing_pairs.add ( pair )
    new_index_bonds = []
    for i, j in sorted ( existing_pairs ):
        new_index_bonds.append ( i )
        new_index_bonds.append ( j )

    vismol_object.bonds = None
    vismol_object.non_bonded_atoms = None
    vismol_object.index_bonds = np.array ( new_index_bonds, dtype = np.int64 )

    vismol_object._bonds_from_pair_of_indexes_list ( )
    vismol_object._get_non_bonded_from_bonded_list ( )
    vismol_object._generate_topology_from_index_bonds ( )
    vismol_object.define_molecules ( )
    vismol_object.define_Calpha_backbone ( )

    vismol_object.create_representation ( rep_type = "lines" )
    vismol_object.create_representation ( rep_type = "nonbonded" )
    vismol_object.core_representations["picking_dots"] = None
    vismol_object.core_representations["picking_text"] = None

    vm_session = vismol_object.vm_session
    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.queue_draw ( )

    return True
