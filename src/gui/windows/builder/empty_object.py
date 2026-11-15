#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: Molecule Builder -- empty object creation + pDynamo system sync
#
#  Description:
#      First building block of the "Builder" tool (draw molecules from
#      scratch). Creates a VismolObject AND, since the user explicitly
#      asked for the two to be kept in sync, a matching pDynamo System
#      linked to it via vismol_object.e_id -- REVERSING an earlier,
#      explicit decision (kept here for context, not because it's still
#      true): "the builder starts as a pure visual sketchpad; promoting a
#      drawn molecule into a real pDynamo system is a separate, later
#      step". That later step is THIS one.
#
#      sync_pdynamo_system() (see its own docstring) rebuilds the linked
#      System from scratch to match the VismolObject's CURRENT atoms/
#      bonds/positions -- called by atom_ops.adjust_hydrogens() every
#      time it runs (which itself already runs after every structural
#      edit -- place/replace an atom, a bond-drag, cycling a bond's
#      order, deleting an atom/bond -- see each of those call sites' own
#      comments in click_mode.py/vismol_glcore.py). Same "always fully
#      rebuild rather than incrementally patch" philosophy already used
#      throughout atom_ops.py for the VismolObject's own bonds/topology
#      (_reapply_manual_bonds()) -- NOT wired into move_atom() (called on
#      every single mouse-motion event during a live drag), which would
#      rebuild a whole pDynamo System dozens of times a second for no
#      benefit.
#
#      [EN] BUG FIX / DESIGN CHANGE: create_empty_vismol_object() used to
#      ALSO call sync_pdynamo_system() immediately, right at creation --
#      meaning even opening the Builder sidebar and closing it again
#      without drawing anything already permanently committed an empty
#      system + treeview row to the project (and, separately, caused
#      main_treeview.refresh() -- called after every structural edit --
#      to accumulate one stale duplicate row per edit in
#      main_window.system_liststore, since that ListStore was never
#      cleared between refreshes; both fixed together). Now the object
#      stays a pure, uncommitted, builder-only VisMol object (e_id=None)
#      until the FIRST real structural edit's own call to
#      sync_pdynamo_system() promotes it -- see
#      create_empty_vismol_object()'s own docstring, and
#      discard_builder_object_if_unused() for cleanly removing it again
#      if editing ends with no promotion ever having happened.
#
#      [EN] The underlying pDynamo API used here (Connectivity,
#      Atom.WithOptions, Bond.WithNodes, ConvertInputConnectivity,
#      System.FromConnectivity) was NOT guessed at -- it was read
#      directly out of util/extras/MOL2FileReader.py's own ToSystem()
#      method (already vendored in this repo), which builds a System the
#      exact same way when importing a real MOL2 file. The zero-atom case
#      (creating the system for a brand new, still-empty Builder object)
#      was verified against pDynamo3's own upstream source (github.com/
#      pdynamo/pDynamo3, cloned and read directly, not executed -- no
#      live pDynamo environment available here either) rather than
#      assumed: ConvertInputConnectivity has its own explicit
#      `if len(self.nodes) > 0:` guard (a no-op for zero atoms, not an
#      error), and AtomContainer._SetItemsFromIterable's handling of an
#      empty iterable is ordinary Python (empty list in, empty list out,
#      no special-casing needed). This is still the LEAST-verified part
#      of the whole Builder feature set so far -- everything else in
#      this project could at least be reasoned about against this same
#      codebase's own, already-working patterns; this one relies on
#      reading a separate, external library's source with no live
#      execution feedback at all. Test this specific piece first,
#      in isolation, before relying on the rest.
#
#      Because of the (former) separate-from-pDynamo choice, this module
#      still does NOT use eSession._add_vismol_object() for the
#      VismolObject side of things (see register_builder_object()'s own
#      docstring for why) -- only the pDynamo SYSTEM side now goes
#      through the normal channel (p_session.add_new_system_to_psession(),
#      the exact same method real, file-loaded systems use).
#
#      TREEVIEW DISPLAY: main_treeview.add_vismol_object_to_treeview()
#      IS used now (unlike an earlier version of this module/comment,
#      before the user asked "shouldn't this show up in the treeview?")
#      -- fixed to handle e_id is None gracefully instead of skipping it
#      entirely. That fix was needed for a real reason beyond just
#      "make it visible", found by tracing through
#      main_treeview.refresh() (which already iterates and calls
#      add_vismol_object_to_treeview() on EVERY vm_objects_dic entry
#      unconditionally): the very first refresh() after creating a
#      Builder object would have raised KeyError(None) there, before
#      this fix. This e_id-is-None path is now the NORMAL state for a
#      whole Builder session that hasn't made a structural edit yet (see
#      the module-level note above about create_empty_vismol_object() no
#      longer promoting to a real system immediately) -- not just a rare
#      failure fallback as it was when this comment was first written.
#
#      Verified safe against a zero-atom object (this is exactly what an
#      empty builder object is) by reading, but not executing (no live
#      pDynamo/GTK environment available), the relevant code paths:
#        - VismolObject.__init__: self.atoms = {} by default, no crash risk.
#        - VismolObject._generate_color_vectors(): explicitly returns False
#          early when len(self.atoms) == 0 (already handled upstream, not
#          something added for the builder).
#        - VismolObject.build_core_representations() /
#          create_representation(): pass indexes=self.index_bonds /
#          indexes=self.non_bonded_atoms, both None by default for an atom-
#          less object -- previously a crash (see the Representation.__init__
#          fix earlier in this project's history: "indexes pode chegar como
#          None"), now safely normalised to an empty array there.
#        - vm_glcore.render(): iterates self.vm_session.vm_objects_dic.values()
#          directly, calling build_core_representations() and looping over
#          vm_object.representations.values() (skipping None entries) -- no
#          dependency on p_session/psystem anywhere in this loop.
#
from util.debug import dprint
import numpy as np
from vismol.core.vismol_object import VismolObject


def register_builder_object ( vm_session, vismol_object, show_molecule = False ):
    """ Minimal, psystem-free registration for a builder-created object.
    Mirrors ONLY the bookkeeping parts of eSession._add_vismol_object()
    that do not require vismol_object.e_id to reference a real pDynamo
    system: unique index assignment, key6 tag, unique name -- PLUS a
    treeview row (main_treeview.add_vismol_object_to_treeview(), now
    that it handles e_id is None -- see module docstring). Deliberately
    still skips: e_id/e_tag handling and psystem lookups (see module
    docstring above for why). """
    vismol_object.key6 = vm_session.gen_random_tag_string ( length = 6 )

    vismol_object.index = vm_session.vm_object_counter
    vm_session.vm_objects_dic[vismol_object.index] = vismol_object
    vm_session.vm_object_counter += 1

    vismol_object.name = vismol_object.name.replace ( " ", "_" )
    while vismol_object.name in vm_session.vobject_names:
        vismol_object.name = "{}_X".format ( vismol_object.name )
    vm_session.vobject_names[vismol_object.name] = vismol_object

    # Flag used elsewhere (now and in future code) to recognise a
    # builder-only object that has no pDynamo system behind it yet --
    # e.g. to skip/guard any code path that would otherwise assume
    # vismol_object.e_id is a valid psystem key.
    vismol_object.is_builder_only = True
    vismol_object.e_id = None

    # [EN] 2026-10-09 BUG FIX: same model matrix as every other object (as
    # pDynamo2EasyHybrid's own object creation does). Left at identity, a
    # Builder object created after the scene had been centred/panned (e.g.
    # the DYFF window centring on a previous, since-deleted object) was drawn
    # off the pk1-pk4 markers and atom labels, which use vm_glcore.model_mat.
    vm_glcore = getattr ( vm_session, "vm_glcore", None )
    if vm_glcore is not None:
        vismol_object.set_model_matrix ( vm_glcore.model_mat )

    if show_molecule:
        vismol_object.create_representation ( rep_type = "sticks" )
        vismol_object.create_representation ( rep_type = "stick_spheres" )
        from gui.windows.builder.atom_ops import _activate_new_sphere_representation
        _activate_new_sphere_representation ( vismol_object, "stick_spheres" )
        vismol_object.create_representation ( rep_type = "nonbonded" )

    # vm_session.main is set by main_window.py during normal app startup
    # (self.vm_session.main = self) -- guarded with getattr because it is
    # NOT set in the headless test harness used to validate this module
    # (no real main window there), so tests can still exercise everything
    # above this line without needing a GTK main window at all.
    main = getattr ( vm_session, "main", None )
    if main is not None and getattr ( main, "main_treeview", None ) is not None:
        main.main_treeview.add_vismol_object_to_treeview ( vismol_object )

    return vismol_object


def create_empty_vismol_object ( vm_session, name = "new_molecule" ):
    """ Creates and registers a brand-new, completely empty (zero atoms)
    VismolObject -- the starting point of the molecule Builder. Does not
    show any representation yet (nothing to show with zero atoms); once
    atom-placement is implemented, representations can be (re)created as
    atoms are added.

    [EN] BUG FIX / DESIGN CHANGE: this used to call sync_pdynamo_system()
    immediately here, which -- since sync_pdynamo_system() treats "no
    e_id yet" as "first call, register a real pDynamo system now" --
    meant simply opening the Builder sidebar and closing it again,
    WITHOUT drawing anything, already permanently committed an empty
    system + a treeview row to the project. Not calling it here at all
    keeps the freshly-created object a pure, uncommitted, builder-only
    VisMol object (e_id=None, is_builder_only=True, both already set by
    register_builder_object() above) -- the FIRST real structural edit
    (add_atom()/remove_atom()/add_bond()/set_atom_element(), all of which
    already end in a call to sync_pdynamo_system() via adjust_hydrogens()
    or their own explicit call -- see click_mode.py's various handlers)
    is what promotes it to a real, tracked system, reusing
    sync_pdynamo_system()'s own existing first-call logic unchanged. See
    discard_builder_object_if_unused() below for the other half of this:
    cleanly removing the object again if editing ends with no promotion
    ever having happened. """
    vismol_object = VismolObject ( vismol_session = vm_session, name = name, active = True )

    # frames shape convention: (n_frames, n_atoms, 3) -- see e.g.
    # generate_new_empty_vismol_object() in pDynamo2EasyHybrid/session.py,
    # which uses np.empty([0, n_atoms, 3]) for an existing system awaiting
    # a trajectory. Here we have neither frames nor atoms yet, so
    # (1, 0, 3): one (empty) frame, zero atoms -- keeps `len(vismol_object.frames)`
    # and `.shape[0]`-based frame-count code elsewhere well-defined (1, not 0),
    # while still being unambiguously "no atoms".
    vismol_object.frames = np.zeros ( (1, 0, 3), dtype = np.float32 )

    register_builder_object ( vm_session, vismol_object, show_molecule = False )

    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.queue_draw ( )

    return vismol_object


def discard_builder_object_if_unused ( vismol_object ):
    """ Called when Builder editing mode ends (see builder_sidebar.py's
    close_window()) -- cleanly removes `vismol_object` from the session
    if it was NEVER promoted to a real pDynamo system (its own
    `is_builder_only` flag is still True, meaning sync_pdynamo_system()
    was never called with a structural edit behind it -- see
    create_empty_vismol_object()'s own docstring above). A no-op if a
    real edit DID happen (is_builder_only is False): the object stays
    exactly as today, a normal part of the project.

    Deliberately does NOT touch p_session.psystem/system_liststore/the
    treestore's system row directly -- an unpromoted builder object never
    got one (see sync_pdynamo_system(): the "else" branch that calls
    add_new_system_to_psession() only ever runs on a REAL edit). Only
    vm_objects_dic/vobject_names (registered by register_builder_object())
    need cleaning up, followed by a plain treeview refresh so its
    (vobject-only, e_id=None) row disappears along with it. """
    if vismol_object is None:
        return
    if getattr ( vismol_object, "is_builder_only", False ) is not True:
        return

    vm_session = vismol_object.vm_session
    index      = getattr ( vismol_object, "index", None )
    name       = getattr ( vismol_object, "name", None )

    if index is not None:
        vm_session.vm_objects_dic.pop ( index, None )
    if name is not None:
        vm_session.vobject_names.pop ( name, None )

    main = getattr ( vm_session, "main", None )
    if main is not None and getattr ( main, "main_treeview", None ) is not None:
        main.main_treeview.refresh ( )

    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.queue_draw ( )


def _build_pdynamo_system_from_vismol_object ( vismol_object, label = None ):
    """ [EN] Builds a brand-new pDynamo System matching vismol_object's
    CURRENT atoms/bonds/positions -- see this module's own docstring
    (top of file) for where this exact sequence of calls came from
    (util/extras/MOL2FileReader.py's ToSystem()) and how the zero-atom
    case was verified against pDynamo3's own upstream source.

    Bond orders: BondType.Single/Double/Triple, from vismol_object.
    manual_bond_orders (defaulting to Single for any pair not explicitly
    recorded there -- see atom_ops.add_bond()'s own docstring for why
    EVERY bond, regardless of order, ends up in manual_bonds).

    Aromaticity: isAromatic, from vismol_object.manual_aromatic_bonds --
    a SEPARATE flag from the numeric bond type, mirroring pDynamo3's own
    Bond model (pMolecule/Bond.py: "Aromatic bond types are not defined.
    Instead bonds themselves are flagged if they form part of an
    aromatic system.").

    Returns the new System (not yet registered with any eSession/
    treeview -- see sync_pdynamo_system(), which calls this and then
    handles registration). """
    return _build_system_from_snapshot ( _pdynamo_snapshot ( vismol_object, label ) )


#-----------------------------------------------------------------------------
#  [EN] 2026-10-17, user's report ("editing a big system is still very slow"):
#  the pDynamo system is rebuilt from scratch on every Builder edit (~1.3 s at
#  20k atoms). The build is split in two: _pdynamo_snapshot() reads the
#  vismol object (cheap, GTK thread) and _build_system_from_snapshot() makes
#  the System from plain data only -- so sync_pdynamo_system() can run it in a
#  worker thread for big objects (see _BackgroundSync below).
#-----------------------------------------------------------------------------
def _pdynamo_snapshot ( vismol_object, label = None ):
    """ Plain data for _build_system_from_snapshot(): per atom ( symbol, name,
        formal charge, ( chain name, residue name, residue index ) or None ),
        bonds ( i, j, order, aromatic ) in atom positions, coordinates, cell. """
    atom_ids = sorted ( vismol_object.atoms.keys ( ) )
    atoms = [ ]
    for atom_id in atom_ids:
        atom    = vismol_object.atoms[atom_id]
        residue = getattr ( atom, "residue", None )
        chain   = getattr ( atom, "chain", None ) or getattr ( residue, "chain", None )
        where   = None
        if residue is not None and chain is not None:
            where = ( str ( getattr ( chain, "name", "" ) or "" ), residue.name, residue.index )
        atoms.append ( ( atom.symbol, atom.name, int ( getattr ( atom, "formal_charge", 0 ) or 0 ), where ) )
    position_of = { atom_id: k for k, atom_id in enumerate ( atom_ids ) }
    orders   = getattr ( vismol_object, "manual_bond_orders", None ) or { }
    aromatic = getattr ( vismol_object, "manual_aromatic_bonds", None ) or set ( )
    bonds = [ ( position_of[i], position_of[j], orders.get ( ( i, j ), 1 ), ( i, j ) in aromatic )
              for ( i, j ) in ( getattr ( vismol_object, "manual_bonds", None ) or set ( ) )
              if i in position_of and j in position_of ]
    coordinates = np.array ( vismol_object.frames[0][atom_ids], dtype = np.float64 ) if atom_ids else np.zeros ( ( 0, 3 ) )
    return { "label": label if label else vismol_object.name, "atoms": atoms, "bonds": bonds,
             "coordinates": coordinates, "cell": getattr ( vismol_object, "builder_cell", None ) }


def _build_system_from_snapshot ( snapshot ):
    """ The System of a snapshot: with a sequence when possible
        (_build_sequenced_from_snapshot()), otherwise from the connectivity
        alone; the Builder cell as its symmetry. Plain data only: safe in a
        worker thread. """
    from pMolecule             import Atom, Bond, BondType, Connectivity, ConvertInputConnectivity, System
    from pScientific            import PeriodicTable
    from pScientific.Geometry3  import Coordinates3

    # [EN] 2026-10-02 (user request: hydrogens added in the Builder must
    # belong to their residues): keep chains/residues in the pDynamo system
    # too -- System.FromConnectivity() (below) builds NO sequence, so a PDB
    # protein edited in the Builder used to come back with every residue
    # gone. Tried first; the old path stays as the fallback.
    system = _build_sequenced_from_snapshot ( snapshot )
    if system is None:
        connectivity = Connectivity ( )
        # [EN] 2026-09-25: the Builder's own formal_charge goes into the real
        # pDynamo Atom (formalCharge, read by DYFF's charged patterns and by
        # anything that sums the total charge).
        for ( symbol, name, formal_charge, where ) in snapshot["atoms"]:
            connectivity.AddNode ( Atom.WithOptions ( atomicNumber = PeriodicTable.AtomicNumber ( symbol ), label = name,
                                                       formalCharge = formal_charge ) )
        bond_type_by_order = { 1: BondType.Single, 2: BondType.Double, 3: BondType.Triple }
        for ( i, j, order, is_aromatic ) in snapshot["bonds"]:
            connectivity.AddEdge ( Bond.WithNodes ( connectivity.nodes[i], connectivity.nodes[j], isAromatic = is_aromatic,
                                                     type = bond_type_by_order.get ( order, BondType.Single ) ) )
        ConvertInputConnectivity ( connectivity, { } )
        system = System.FromConnectivity ( connectivity = connectivity )
        system.label = snapshot["label"]
        coordinates3 = Coordinates3.WithExtent ( len ( snapshot["atoms"] ) )
        for k, xyz in enumerate ( snapshot["coordinates"] ):
            coordinates3[k, 0], coordinates3[k, 1], coordinates3[k, 2] = float ( xyz[0] ), float ( xyz[1] ), float ( xyz[2] )
        system.coordinates3 = coordinates3
    _apply_cell ( snapshot["cell"], system )
    return system


def _apply_builder_cell ( vismol_object, system ):
    """ [EN] 2026-10-03: a Builder object with a periodic cell (set by the
    Solvate tool, atom_ops.set_builder_cell()) gets it as the pDynamo
    system's symmetry. """
    _apply_cell ( getattr ( vismol_object, "builder_cell", None ), system )


def _apply_cell ( cell, system ):
    """ cubic when a = b = c with right angles, orthorhombic for other
        right-angled boxes, triclinic otherwise. """
    if not cell:
        return
    try:
        from pScientific.Symmetry import ( PeriodicBoundaryConditions, CrystalSystemCubic,
                                           CrystalSystemOrthorhombic, CrystalSystemTriclinic )
        a, b, c, alpha, beta, gamma = cell
        right = all ( abs ( angle - 90.0 ) < 1e-6 for angle in ( alpha, beta, gamma ) )
        if right and abs ( a - b ) < 1e-6 and abs ( a - c ) < 1e-6:
            crystal, parameters = CrystalSystemCubic ( ), { "a": a }
        elif right:
            crystal, parameters = CrystalSystemOrthorhombic ( ), { "a": a, "b": b, "c": c }
        else:
            crystal, parameters = CrystalSystemTriclinic ( ), { "a": a, "b": b, "c": c,
                                                                "alpha": alpha, "beta": beta, "gamma": gamma }
        system.symmetry           = PeriodicBoundaryConditions.WithCrystalSystem ( crystal )
        system.symmetryParameters = system.symmetry.MakeSymmetryParameters ( **parameters )
    except Exception as exc:
        dprint ( "empty_object._apply_cell: could not set the cell ({})".format ( exc ) )


def _build_sequenced_pdynamo_system ( vismol_object, label = None ):
    """ The sequenced System of vismol_object, or None (see
        _build_sequenced_from_snapshot()). """
    system = _build_sequenced_from_snapshot ( _pdynamo_snapshot ( vismol_object, label ) )
    return system


def _build_sequenced_from_snapshot ( snapshot ):
    """ [EN] Same system as _build_system_from_snapshot() but WITH a sequence
    (entity = vismol chain name, component = "RESN.resi", atom label = atom
    name), built through Sequence.FromAtomPaths() + System.FromSequence() --
    the path util/hydrogen_builder.py already uses. Returns None (caller
    falls back to the sequence-less build) when:
      - any atom has no chain/residue, or a residue has two atoms with the
        same name (pDynamo atom paths must be unique);
      - a residue's atoms are not contiguous in atom_id order, or a chain
        is split -- pDynamo would group them and REORDER the atoms, and
        the rest of the Builder/EasyHybrid assumes pDynamo atom k is
        vismol atom k;
      - pDynamo raises for any other reason. """
    try:
        from pMolecule             import Atom, BondType, ConvertInputConnectivity, System
        from pMolecule.Sequence    import Sequence
        from pScientific           import PeriodicTable
        from pScientific.Geometry3 import Coordinates3

        if not snapshot["atoms"]:
            return None
        paths, atoms_list = [ ], [ ]
        seen_paths = set ( )
        seen_residues, seen_chains = set ( ), set ( )
        last_residue_key, last_chain = None, None
        for ( symbol, name, formal_charge, where ) in snapshot["atoms"]:
            if where is None:
                return None
            chain_name, residue_name, residue_index = where
            residue_key = where
            if residue_key != last_residue_key:
                if residue_key in seen_residues:
                    return None                       # residue split: would be reordered
                if chain_name != last_chain:
                    if chain_name in seen_chains:
                        return None                   # chain split: would be reordered
                    seen_chains.add ( chain_name )
                    last_chain = chain_name
                seen_residues.add ( residue_key )
                last_residue_key = residue_key
            path = "{}:{}.{}:{}".format ( chain_name, residue_name, residue_index, name )
            if path in seen_paths:
                return None
            seen_paths.add ( path )
            paths.append ( path )
            atoms_list.append ( Atom.WithOptions ( atomicNumber = PeriodicTable.AtomicNumber ( symbol ),
                                                   label = name, formalCharge = formal_charge ) )
        bond_type_by_order = { 1: BondType.Single, 2: BondType.Double, 3: BondType.Triple }
        # (i, j, type, isAromatic) tuples with INTEGER positions, not Bond
        # objects: pDynamo's Bond.FromIterable checks `bond.node1 in nodes` on
        # a plain list for every Bond object (O(atoms) per bond -- ~2 s of a
        # 6000-atom rebuild, i.e. most of a Builder Undo), while integer
        # positions are a direct nodes[i]. Positions == final atom order,
        # which is verified right after System.FromSequence() below.
        bonds = [ ( i, j, bond_type_by_order.get ( order, BondType.Single ), is_aromatic )
                  for ( i, j, order, is_aromatic ) in snapshot["bonds"] ]

        sequence = Sequence.FromAtomPaths ( paths, atoms = atoms_list )
        system   = System.FromSequence ( sequence, bonds = bonds )
        # same normalisation the sequence-less path applies before
        # System.FromConnectivity() (atom/bond aromaticity flags, ring
        # information) -- pDynamo's own typing relies on it
        ConvertInputConnectivity ( system.connectivity, { } )
        # pDynamo atom k must be vismol atom k (see the docstring)
        for k, atom in enumerate ( atoms_list ):
            if atom.index != k or system.atoms[k] is not atom:
                return None
        system.label = snapshot["label"]
        coordinates3 = Coordinates3.WithExtent ( len ( atoms_list ) )
        for k, xyz in enumerate ( snapshot["coordinates"] ):
            coordinates3[k, 0], coordinates3[k, 1], coordinates3[k, 2] = float ( xyz[0] ), float ( xyz[1] ), float ( xyz[2] )
        system.coordinates3 = coordinates3
        return system
    except Exception as exc:
        dprint ( "empty_object._build_sequenced_from_snapshot: falling back to a sequence-less "
                 "system ({})".format ( exc ) )
        return None


# [EN] 2026-10-17: background rebuild of a big Builder object's pDynamo system.
BACKGROUND_SYNC_MIN_ATOMS = 1500
BACKGROUND_SYNC_DELAY_MS  = 300   # the worker starts after this pause: clicks and the redraw don't share the GIL with it
_BACKGROUND_JOBS = { }          # id ( vismol_object ) -> _BackgroundSync


class _BackgroundSync:
    """ Rebuilds one Builder object's pDynamo system in a worker thread.
        request(): snapshot now (GTK thread) and build it; a request made while
        a build runs is built next, so bursts of edits collapse into the
        newest one. The newest result is installed in an idle callback. Until
        then p_session.psystem has a pending flush for the entry: any reader
        waits for the build (flush()), so nothing reads a stale system. """

    def __init__ ( self, vismol_object, main, e_id ):
        import threading
        self.obj, self.main, self.e_id = vismol_object, main, e_id
        self.lock    = threading.Lock ( )
        self.gen     = 0             # newest request
        self.latest  = None          # ( gen, snapshot ) of the newest request
        self.queued  = None          # ( gen, snapshot ) not taken by the worker yet
        self.result  = None          # ( gen, system or None, error )
        self.applied = 0             # newest gen installed (or given up)
        self.thread  = None
        self.timer   = None          # GLib source that starts the worker

    def request ( self ):
        from gi.repository import GLib
        snapshot = _pdynamo_snapshot ( self.obj, label = self.obj.name )
        with self.lock:
            self.gen += 1
            self.latest = self.queued = ( self.gen, snapshot )
        psystem = self.main.p_session.psystem
        if hasattr ( psystem, "set_pending" ):
            psystem.set_pending ( self.e_id, self.flush )
        # . (re)start the pause: a burst of edits builds once, after the last one
        if self.timer is not None:
            GLib.source_remove ( self.timer )
        self.timer = GLib.timeout_add ( BACKGROUND_SYNC_DELAY_MS, self._start_worker )

    def _start_worker ( self ):
        import threading
        self.timer = None
        with self.lock:
            if self.thread is not None or self.queued is None:
                return False                         # running (takes the queued job next) or nothing to do
            self.thread = threading.Thread ( target = self._work, daemon = True )
        self.thread.start ( )
        return False

    def _stop_timer ( self ):
        if self.timer is not None:
            from gi.repository import GLib
            GLib.source_remove ( self.timer )
            self.timer = None

    def cancel ( self ):
        """ A synchronous rebuild replaces whatever is pending. """
        self._stop_timer ( )
        with self.lock:
            self.applied = self.gen
            self.queued = None
        psystem = self.main.p_session.psystem
        if hasattr ( psystem, "clear_pending" ):
            psystem.clear_pending ( self.e_id )

    def _work ( self ):
        from gi.repository import GLib
        while True:
            with self.lock:
                job, self.queued = self.queued, None
                if job is None:
                    self.thread = None
                    return
            gen, snapshot = job
            try:
                system, error = _build_system_from_snapshot ( snapshot ), None
            except Exception as exc:
                system, error = None, exc
            with self.lock:
                self.result = ( gen, system, error )
            GLib.idle_add ( self._install_idle )

    def _install_idle ( self ):
        try:
            self._install ( )
        except Exception as exc:
            dprint ( "empty_object._BackgroundSync: install failed ({})".format ( exc ) )
        return False

    def _install ( self ):
        with self.lock:
            result, newest = self.result, self.gen
        if result is None:
            return
        gen, system, error = result
        if gen != newest or gen <= self.applied:
            return                                   # a newer build is coming, or this one is done
        self.applied = gen
        p_session = self.main.p_session
        p_session.psystem.clear_pending ( self.e_id )
        obj = self.obj
        if error is not None or system is None:
            dprint ( "WARNING empty_object.sync_pdynamo_system (background): could not rebuild the pDynamo "
                     "system of '{}' -- Builder editing continues unaffected. Error: {}".format ( obj.name, error ) )
            return
        if ( obj not in obj.vm_session.vm_objects_dic.values ( ) or getattr ( obj, "e_id", None ) != self.e_id
                or p_session.psystem.peek ( self.e_id ) is None ):
            return
        _swap_rebuilt_system_into_psystem ( p_session, self.e_id, system )
        self.main.main_treeview.refresh ( )

    def flush ( self ):
        """ Called by p_session.psystem when the entry is read: the newest
            system installed now (waiting for the worker, or building it). """
        self._stop_timer ( )
        with self.lock:
            self.queued = None                       # built right here if the worker has not taken it
        thread = self.thread
        if thread is not None:
            thread.join ( )
        with self.lock:
            result, latest, newest = self.result, self.latest, self.gen
        if self.applied >= newest:
            return
        if result is None or result[0] != newest:
            try:
                system, error = _build_system_from_snapshot ( latest[1] ), None
            except Exception as exc:
                system, error = None, exc
            with self.lock:
                self.result = ( newest, system, error )
        self._install ( )


def _background_sync_job ( vismol_object, main ):
    """ The background job for this Builder object, or None when the object
        is synced synchronously (small, not the Builder target, not
        registered yet, or EASYHYBRID_BUILDER_SYNC=sync). """
    import os
    job = _BACKGROUND_JOBS.get ( id ( vismol_object ) )
    if job is not None and job.obj is not vismol_object:
        job = None
    e_id = getattr ( vismol_object, "e_id", None )
    psystem = main.p_session.psystem
    eligible = ( os.environ.get ( "EASYHYBRID_BUILDER_SYNC", "" ).lower ( ) != "sync"
                 and hasattr ( psystem, "peek" )
                 and vismol_object is getattr ( vismol_object.vm_session, "builder_target_object", None )
                 and e_id is not None and psystem.peek ( e_id ) is not None
                 and len ( vismol_object.atoms ) >= BACKGROUND_SYNC_MIN_ATOMS )
    if not eligible:
        if job is not None:
            job.cancel ( )
            _BACKGROUND_JOBS.pop ( id ( vismol_object ), None )
        return None
    if job is None or job.e_id != e_id:
        if job is not None:
            job.cancel ( )
        job = _BACKGROUND_JOBS[id ( vismol_object )] = _BackgroundSync ( vismol_object, main, e_id )
    return job


def sync_pdynamo_system ( vismol_object ):
    """ [EN] (Re)builds vismol_object's LINKED pDynamo System from
    scratch (_build_pdynamo_system_from_vismol_object() above) and (re)
    registers it, so the two never drift apart -- see this module's own
    docstring for the full list of when this gets called and why NOT
    from move_atom()/every mouse-motion event.

    First call for a given vismol_object (no e_id yet): goes through
    p_session.add_new_system_to_psession() -- the SAME registration
    method real, file-loaded systems use, so the new system gets every
    piece of EasyHybrid bookkeeping (e_working_folder, e_selections,
    e_charges_backup, ...) that method sets up, not just the bare
    pDynamo object -- then links vismol_object.e_id to the new system's
    e_id and flips is_builder_only to False (it now genuinely isn't
    builder-only anymore).

    Every call after that (vismol_object.e_id already set): rebuilds the
    System object itself from scratch (structure changed -- new/removed
    atoms or bonds, or a moved position), but reuses the SAME e_id slot
    in p_session.psystem instead of registering a new one -- re-running
    add_new_system_to_psession() every time would both leak a fresh e_id
    per edit (a new treeview row every time, instead of updating the
    existing one) AND needlessly redo one-time setup (working folder,
    colour palette, ...) that should stay stable for this object's whole
    lifetime. Any 'e_'-prefixed attribute already present on the OLD
    system object (all that one-time EasyHybrid bookkeeping) is copied
    onto the freshly-rebuilt one first, so it survives the swap.

    Finally, calls main_treeview.refresh() -- a full treeview rebuild
    (clears and re-adds every system + every vobject, see that method's
    own code) rather than trying to surgically patch just this one row;
    simpler and safe to call this often for a Builder-sized session.

    Silently does nothing (prints a warning instead of raising) if
    anything in the pDynamo/eSession call chain fails -- the Builder's
    OWN vismol_object state is the source of truth either way, so a
    pDynamo-sync failure shouldn't take down editing itself. """
    # [EN] 2026-10-04 BUG FIX: never (re)create a pDynamo system for an
    # object that is no longer in the session (deleted from the treeview
    # while a Builder tool still held it) -- that made orphan systems.
    if vismol_object not in vismol_object.vm_session.vm_objects_dic.values ( ):
        dprint ( "empty_object.sync_pdynamo_system: '{}' is no longer in the session -- not synced.".format (
                 getattr ( vismol_object, "name", "?" ) ) )
        return
    try:
        main       = vismol_object.vm_session.main
        p_session  = main.p_session

        # [EN] 2026-10-17: big Builder objects are rebuilt in a worker thread
        # (_BackgroundSync); the edit returns at once
        job = _background_sync_job ( vismol_object, main )
        if job is not None:
            job.request ( )
            return

        new_system = _build_pdynamo_system_from_vismol_object ( vismol_object, label = vismol_object.name )

        existing_e_id = getattr ( vismol_object, "e_id", None )

        # [EN] 2026-09-25 BUG FIX (user's own report: "quando eu deleto o
        # sistema com o builder aberto, nao consigo mais criar um sistema
        # novo" -- reproduced live): `existing_e_id in p_session.psystem`
        # only checks the KEY exists, not that a real System actually
        # occupies it -- but p_session.psystem[e_id] can legitimately be
        # None (pDynamoSession's own startup convention, AND main_window.
        # delete_system()'s own "no systems left" fallback -- see
        # [[project_treeview_refresh_none_system_bug]]). If the Builder's
        # OWN target object was still open and being edited when its
        # underlying system got deleted out from under it (e_id stays
        # set on the orphaned vismol_object), the NEXT edit's own sync
        # would find `existing_e_id in psystem` True (the stale key is
        # still there, just pointing at None) and take the SWAP branch
        # below instead of this one -- _swap_rebuilt_system_into_psystem()
        # only copies bookkeeping from an EXISTING system (a no-op when
        # that slot is None), so the newly-built System ends up missing
        # EVERY one-time 'e_'-prefixed attribute add_new_system_to_
        # psession() normally sets up (e_color_palette, e_working_folder,
        # e_selections, ...) -- the very next main_treeview.refresh()
        # (called a few lines below, for EVERY system in psystem, not
        # just this one) then crashes on THAT system's own missing
        # e_color_palette while rendering its row, silently (caught by
        # this function's own broad except) but PERMANENTLY, since the
        # broken system is never removed from psystem and every FUTURE
        # refresh() (for ANY reason, not just this object) hits the same
        # crash before it can render anything -- an empty treeview
        # forever after, even though brand-new systems keep registering
        # correctly in psystem itself. Checking the VALUE (not just the
        # key) correctly falls through to the full add_new_system_to_
        # psession() registration instead, which sets up every required
        # attribute from scratch, exactly as if this were a genuinely
        # new system (which, functionally, it now is).
        if existing_e_id is not None and p_session.psystem.get ( existing_e_id ) is not None:
            _swap_rebuilt_system_into_psystem ( p_session, existing_e_id, new_system )
        else:
            p_session.add_new_system_to_psession ( system = new_system, name = vismol_object.name )
            vismol_object.e_id            = new_system.e_id
            vismol_object.is_builder_only = False

        main.main_treeview.refresh ( )

    except Exception as exc:
        dprint ( "WARNING empty_object.sync_pdynamo_system: failed to (re)build/register the linked "
                "pDynamo system for '{}' -- Builder editing continues unaffected. Error: {}".format (
                getattr ( vismol_object, "name", "?" ), exc ) )


def _swap_rebuilt_system_into_psystem ( p_session, e_id, new_system ):
    """ [EN] Installs new_system into p_session.psystem[e_id], preserving
    every 'e_'-prefixed EasyHybrid bookkeeping attribute (e_working_folder,
    e_selections, e_qc_table, ...) already present on whatever system
    currently occupies that slot -- factored out of sync_pdynamo_system()
    (above) so finish_editing_existing_system() (below) can reuse the
    exact same "keep the same e_id, don't lose one-time setup" logic when
    folding an "Edit in Builder" session back onto the ORIGINAL system it
    started from, instead of duplicating it. """
    old_system = p_session.psystem.get ( e_id )
    if old_system is not None:
        for attr_name, attr_value in old_system.__dict__.items ( ):
            if attr_name.startswith ( "e_" ) and not hasattr ( new_system, attr_name ):
                setattr ( new_system, attr_name, attr_value )
    new_system.e_id = e_id
    p_session.psystem[e_id] = new_system


def hide_other_vobjects_for_builder ( vm_session, keep_vobject ):
    """ [EN] 2026-09-25, user's own explicit request: "quando um novo
    vobject e criado para edicao, temos que desativar a visualizacao de
    todos os outros objetos existentes" -- called from builder_sidebar.
    py's open_window() (covers all 3 real callers: a blank "New" canvas,
    and both "Edit in Builder" entry points via begin_editing_existing_
    system(), which sets vm_session.builder_target_object BEFORE calling
    sidebar.open_window()) right after the target object is resolved.

    Deliberately HIDES (sets `.active = False`, the SAME flag render()'s
    own top-level "for vm_object in vm_objects_dic.values(): if vm_
    object.active: ..." loop already checks) rather than discarding
    anything -- unlike the REMOVED confirm_and_discard_other_vobjects()
    (see [[project_builder_edit_nondestructive_fix]] -- deleting a
    system's other vobjects turned out to be destructive exactly when
    editing became its own new system), this touches nothing about any
    object's identity/data, only its OWN CURRENT visibility, and is
    fully reversible by restore_hidden_vobjects_for_builder() below.
    This is also what keeps the temp clone's own ORIGINAL vobject (a
    near-identical geometric duplicate, for an "Edit in Builder"
    session) from rendering doubled-up with the clone being edited.

    Remembers exactly which vobjects IT deactivated (vm_session.
    builder_hidden_vobjects) -- one that was ALREADY inactive for some
    unrelated reason before this ran is left alone and never gets
    reactivated by the restore step, avoiding a surprising side effect. """
    hidden = [ ]
    for vobject in vm_session.vm_objects_dic.values ( ):
        if vobject is keep_vobject:
            continue
        if getattr ( vobject, "active", False ):
            vobject.active = False
            hidden.append ( vobject )
    vm_session.builder_hidden_vobjects = hidden
    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.queue_draw ( )


def restore_hidden_vobjects_for_builder ( vm_session ):
    """ Counterpart to hide_other_vobjects_for_builder() above -- called
    from builder_sidebar.py's close_window(). Reactivates exactly the
    vobjects THIS Builder session deactivated (harmless no-op for any
    that no longer exist, e.g. deleted mid-session -- just sets an
    attribute on an orphaned Python object nothing renders anymore). """
    hidden = getattr ( vm_session, "builder_hidden_vobjects", None ) or [ ]
    for vobject in hidden:
        vobject.active = True
    vm_session.builder_hidden_vobjects = [ ]
    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.queue_draw ( )


def begin_editing_existing_system ( vm_session, e_id, vismol_object = None ):
    """ [EN] Entry point for "Edit in Builder" on an already-loaded system
    (treeview_menu.py's _menu_edit_in_builder()) -- see the plan file
    (Builder: edit an existing system + dihedral rotation tool) for the
    full design.

    `vismol_object`: 2026-09-24 addition -- the SPECIFIC VismolObject to
    edit, when the caller already knows it (both "Edit in Builder" entry
    points now show an explicit vobject picker, see each entry point's
    own dialog in builder_entry_dialog.py). When None (old callers, or a
    system that's already down to exactly one real vobject), falls back
    to the OLD inline lookup below -- fixed to exclude is_surface objects
    (it didn't before: a molecular-surface child object sharing the same
    e_id as its parent could get picked as the "original" to edit purely
    by dict-iteration-order/`.active` luck, since only `get_active_
    vobject()` ever excluded surfaces, and this function had its own
    separate, non-excluding copy of that same lookup).

    [EN] 2026-09-25, user's own explicit request ("nao devemos misturar
    os dois comportamentos... a operacao de 'criar novo sistema' deve
    ser nao destrutiva: o sistema de origem deve permanecer exatamente
    como estava"): this NO LONGER discards any OTHER vobject sharing
    `e_id` before editing (a real, now-fixed bug: any sibling used to be
    deleted unconditionally here, corrupting the original system even
    when the edit ultimately became its own separate new system, per
    finish_editing_existing_system()'s own atom-count-changed outcome --
    see that function's own updated docstring for the 2 possible
    outcomes). The original system, and every vobject in it, is left
    completely alone by this function -- only a fresh CLONE is ever
    touched during editing.

    Editing never touches the original system directly: this clones it
    first (p_session.clone_system() -- the SAME machinery the treeview's
    own "Clone System" menu item already uses, so the clone gets its own
    fresh e_id via add_new_system_to_psession()'s plain auto-increment
    counter, entirely independent of the original), bootstraps the
    clone's VismolObject into a Builder-editable state (manual_bonds/
    manual_bond_orders/manual_aromatic_bonds, via atom_ops.
    bootstrap_manual_bonds_from_existing() -- a normally-loaded object has
    none of these yet), and sets it as vm_session.builder_target_object so
    opening the Builder sidebar resumes editing THIS object instead of
    creating a blank one (see builder_sidebar.py's open_window(), which
    only calls create_empty_vismol_object() when builder_target_object is
    still None).

    If the Builder sidebar is already open (mid-edit of something else),
    it's closed first -- going through its own close_window(), so
    whatever was being edited there is correctly finalised/discarded
    first, matching this project's existing "one editable object at a
    time" rule (Section 1.1) rather than leaving two half-finished
    sessions tangled together.

    Returns the new temporary (clone) VismolObject, or None (with a
    status-bar message explaining why) if e_id doesn't resolve to a
    system with a usable vobject. """
    main      = vm_session.main
    p_session = main.p_session

    original_system = p_session.psystem.get ( e_id )
    if original_system is None:
        main.bottom_notebook.status_teeview_add_new_item (
            message = 'Edit in Builder: system {} no longer exists.'.format ( e_id ), system = None )
        return None

    original_vobject = vismol_object
    if original_vobject is None:
        for candidate in vm_session.vm_objects_dic.values ( ):
            if ( getattr ( candidate, "e_id", None ) == e_id and getattr ( candidate, "active", True )
                 and not getattr ( candidate, "is_surface", False ) ):
                original_vobject = candidate
                break
        if original_vobject is None:
            for candidate in vm_session.vm_objects_dic.values ( ):
                if getattr ( candidate, "e_id", None ) == e_id and not getattr ( candidate, "is_surface", False ):
                    original_vobject = candidate
                    break
    if original_vobject is None:
        main.bottom_notebook.status_teeview_add_new_item (
            message = 'Edit in Builder: no visible object for system {}.'.format ( e_id ), system = None )
        return None

    # [EN] User's own explicit request: warn before editing a system that
    # already has an MM (force field) or QC (potential) model attached --
    # editing REBUILDS the pDynamo System from scratch (see this
    # function's own docstring and finish_editing_existing_system()
    # below: _build_pdynamo_system_from_vismol_object() creates a brand
    # new System via System.FromConnectivity(), which carries NONE of the
    # old system.mmModel/mmState/qcModel/qcState/nbModel over -- only
    # 'e_'-prefixed EasyHybrid bookkeeping survives the swap, see
    # _swap_rebuilt_system_into_psystem()). Shown here (not duplicated at
    # each entry point) so BOTH the treeview's per-row "Edit in Builder"
    # AND the toolbar's generic "New/Edit" chooser (builder_entry_dialog.
    # py) get it automatically. A "no" answer aborts cleanly -- nothing
    # is cloned, nothing changes, same as any other early-return above.
    has_mm_model = getattr ( original_system, "mmModel", None ) is not None
    has_qc_model = getattr ( original_system, "qcModel", None ) is not None
    if has_mm_model or has_qc_model:
        simple_dialog = getattr ( main, "simple_dialog", None )
        if simple_dialog is not None:
            proceed = simple_dialog.question (
                    "System '{}' already has {} assigned.\n\n"
                    "Editing it in the Builder rebuilds its structure from scratch -- "
                    "this assignment will be LOST, and will need to be redone afterward.{}\n\n"
                    "Continue anyway?".format (
                            "\n\nIt was made with DYFF + residue library: closing the Builder offers to apply "
                            "the library again to the edited structure." if getattr ( original_vobject, "dyff_library_options", None ) else "",
                            original_system.label,
                            " / ".join ( filter ( None, [ "an MM force field" if has_mm_model else None,
                                                           "a QC potential" if has_qc_model else None ] ) ) ) )
            if not proceed:
                return None

    sidebar = getattr ( main, "builder_sidebar_window", None )
    if sidebar is not None and getattr ( sidebar, "visible", False ):
        sidebar.close_window ( )

    before_ids = set ( vm_session.vm_objects_dic.keys ( ) )

    p_session.set_psystem_coordinates_from_vobject ( vobject = original_vobject, system_id = e_id )
    p_session.clone_system ( e_id  = e_id,
                              vobject = original_vobject,
                              name  = "{}_edit".format ( original_system.label ),
                              tag   = 'BLD',
                              color = [ 0, 1, 1 ] )

    new_ids = set ( vm_session.vm_objects_dic.keys ( ) ) - before_ids
    if not new_ids:
        main.bottom_notebook.status_teeview_add_new_item (
            message = 'Edit in Builder: cloning system {} failed.'.format ( e_id ), system = None )
        return None
    temp_vobject = vm_session.vm_objects_dic[ next ( iter ( new_ids ) ) ]

    temp_vobject.builder_edit_source_e_id         = e_id
    temp_vobject.builder_edit_original_atom_count = len ( temp_vobject.atoms )
    # [EN] 2026-09-25, user's own explicit request ("nao devemos misturar
    # os dois comportamentos" -- see this module's own confirm_and_
    # discard_other_vobjects()/its removal, below): a DIRECT reference to
    # the SPECIFIC vobject the user actually picked to edit -- now that
    # siblings sharing the same e_id are no longer discarded (a system
    # can legitimately still have several), finish_editing_existing_
    # system()'s own fold-back path needs to fold back onto THIS ONE
    # specifically, not just "whichever vobject happens to match e_id
    # first" (which would be ambiguous/wrong the moment more than one
    # exists).
    temp_vobject.builder_edit_source_vobject      = original_vobject

    # [EN] BUG FIX (found by live user testing): bootstrapping from
    # temp_vobject's OWN .bonds (the clone's freshly, independently
    # RE-DETECTED topology -- see bootstrap_manual_bonds_from_existing()'s
    # own docstring for exactly why that can silently diverge from what
    # original_vobject actually displays) was dropping real bonds. Reads
    # from original_vobject instead -- the actual ground truth -- and
    # then rebuilds temp_vobject's OWN .bonds/index_bonds/representations
    # from that right away (not deferred to the first structural edit),
    # so what the Builder shows the instant it opens already matches the
    # original, not whatever the clone's own auto-detection produced.
    from gui.windows.builder.atom_ops import ( bootstrap_manual_bonds_from_existing, _reapply_manual_bonds,
                                                 _activate_new_sphere_representation,
                                                 _refresh_bond_dependent_representations )
    bootstrap_manual_bonds_from_existing ( temp_vobject, source_vobject = original_vobject, source_system = original_system )
    # . a system parametrized with DYFF + residue library (Prepare DYFF System): the
    #   options travel with the edit, so closing the Builder can re-apply the library
    if getattr ( original_vobject, "dyff_library_options", None ):
        import copy
        temp_vobject.dyff_library_options = copy.deepcopy ( original_vobject.dyff_library_options )

    temp_vobject.cov_radii_array = None
    temp_vobject.electronegativity_array = None
    temp_vobject.index_bonds = None
    temp_vobject.bonds = None
    temp_vobject.non_bonded_atoms = None
    _reapply_manual_bonds ( temp_vobject )

    temp_vobject.create_representation ( rep_type = "sticks" )
    temp_vobject.create_representation ( rep_type = "stick_spheres" )
    _activate_new_sphere_representation ( temp_vobject, "stick_spheres" )
    temp_vobject.create_representation ( rep_type = "nonbonded" )
    temp_vobject.core_representations["picking_dots"] = None
    temp_vobject.core_representations["picking_text"] = None

    # [EN] BUG FIX (found by live user testing -- a real crash,
    # "IndexError: index N is out of bounds", clicking to place an atom
    # right after opening "Edit in Builder"): temp_vobject is a clone of
    # a NORMALLY-loaded system (via clone_system()), which can carry
    # OTHER representations active by default (e.g. 'cartoon' for a
    # protein) that this Builder never keeps in sync -- deactivate them
    # right now, before the very first click, not just after the first
    # structural edit (atom_ops.py's own add_atom()/remove_atom()/etc.
    # already do this too, for every LATER edit -- see _refresh_bond_
    # dependent_representations()'s own extended docstring for the full
    # story). Harmless no-op for a rep_type this object doesn't have.
    _refresh_bond_dependent_representations ( temp_vobject )

    vm_session.builder_target_object = temp_vobject

    if sidebar is not None:
        sidebar.open_window ( )

    return temp_vobject


def finish_editing_existing_system ( vismol_object, carry_customizations = False ):
    """ [EN] Called from builder_sidebar.py's close_window() when the
    object being edited is a temp clone created by
    begin_editing_existing_system() above (recognised by its own
    builder_edit_source_e_id tag). Implements the user's own confirmed
    rule: if the ATOM COUNT never changed during editing, the edit is
    safe to fold back onto the ORIGINAL system in place (same e_id,
    preserving everything else in the app that keys data by atom index
    against that e_id -- selections, QC region, restraints, ...); if the
    count DID change, folding back would silently desync all of that, so
    the original is left completely untouched and the temp clone (which
    already has its own permanent e_id from clone_system()) simply
    becomes a normal, separate, permanent system.

    [EN] 2026-09-25, user's own explicit request -- this IS the "Caso A"
    (fold back onto the original -- "a edicao deve ocorrer no proprio
    sistema... nao deve ser criado um novo sistema") vs "Caso B" (keep as
    its own new system -- "o sistema original nao deve ser alterado")
    split from the user's own message, still decided the SAME way it
    already was (by whether the atom count changed -- folding back a
    changed atom count is not just a preference, it would genuinely
    desync every OTHER atom-indexed piece of app state, so this isn't a
    free choice to begin with). What actually changed this round: the
    ORIGINAL system's OTHER vobjects (if any) are no longer discarded
    before editing even starts (see begin_editing_existing_system()'s own
    updated docstring) -- so Case B now genuinely leaves the original
    system, EVERY vobject in it, and every property on it, exactly as it
    was, matching the user's own explicit Case B checklist. Case A's own
    fold-back logic below already only ever touched `original_vobject`
    (the ONE vobject actually being edited, via builder_edit_source_
    vobject below) and the shared System object at `source_e_id` --
    never anything ELSE the user might not have wanted touched.

    Either way, vm_session.builder_target_object is left for the caller
    to reset -- this function only deals with the two systems/vobjects
    involved. Returns the vobject that holds the result (the original
    after a fold-back, the clone otherwise; None when not an edit session).
    carry_customizations: also copy the Builder's DYFF customizations
    onto the original in a fold-back (see builder_sidebar.close_window()). """
    source_e_id = getattr ( vismol_object, "builder_edit_source_e_id", None )
    if source_e_id is None:
        return None

    vm_session = vismol_object.vm_session
    main       = vm_session.main
    p_session  = main.p_session

    # [EN] 2026-09-25: prefer the DIRECT reference to the exact vobject
    # begin_editing_existing_system() started from (see its own comment
    # on builder_edit_source_vobject) -- a system can have more than one
    # vobject sharing source_e_id now that siblings are no longer auto-
    # discarded, so "first match by e_id" alone would be ambiguous.
    # Falls back to that old lookup only if the direct reference is
    # missing (an older/edge-case object) or has itself since been
    # removed from the session by something else in the meantime.
    original_vobject = getattr ( vismol_object, "builder_edit_source_vobject", None )
    if original_vobject is not None and original_vobject not in vm_session.vm_objects_dic.values ( ):
        original_vobject = None
    if original_vobject is None:
        for candidate in vm_session.vm_objects_dic.values ( ):
            if getattr ( candidate, "e_id", None ) == source_e_id:
                original_vobject = candidate
                break

    original_count = getattr ( vismol_object, "builder_edit_original_atom_count", None )
    atom_count_changed = ( original_count is None ) or ( len ( vismol_object.atoms ) != original_count )

    # [EN] 2026-09-25: same fix as sync_pdynamo_system() above -- checks
    # the VALUE at source_e_id is a real System, not just that the key
    # exists (p_session.psystem[e_id] can legitimately be None). original_
    # vobject is not None already rules out the most likely way to hit
    # this (delete_system() removes the vobject too), but this stays
    # correct even if some OTHER path ever leaves a vobject behind with
    # no matching live system.
    if not atom_count_changed and original_vobject is not None and p_session.psystem.get ( source_e_id ) is not None:
        from gui.windows.builder.atom_ops import ( _reapply_manual_bonds, _activate_new_sphere_representation,
                                                     set_atom_element )

        n_atoms = len ( vismol_object.atoms )
        original_vobject.frames[0, :n_atoms, :] = vismol_object.frames[0, :n_atoms, :]
        original_vobject.mass_center = np.mean ( original_vobject.frames[0], axis = 0 )

        # [EN] Reuses set_atom_element() (not a plain `.symbol = ...`
        # assignment) so every element-derived per-atom attribute (color/
        # vdw_rad/cov_rad/ball_rad/electronegativity -- see that
        # function's own docstring) gets recomputed too, same as any
        # other in-Builder element change. recompute_bonds/
        # update_representation are both deferred to the single batch
        # rebuild below (this loop may touch many atoms).
        for i in range ( n_atoms ):
            edited_atom = vismol_object.atoms[i]
            if original_vobject.atoms[i].symbol != edited_atom.symbol:
                set_atom_element ( original_vobject, i, edited_atom.symbol, name = edited_atom.name,
                                    recompute_bonds = False, update_representation = False )

        original_vobject.manual_bonds          = set ( getattr ( vismol_object, "manual_bonds", None ) or set ( ) )
        original_vobject.manual_bond_orders    = dict ( getattr ( vismol_object, "manual_bond_orders", None ) or { } )
        original_vobject.manual_aromatic_bonds = set ( getattr ( vismol_object, "manual_aromatic_bonds", None ) or set ( ) )

        original_vobject.cov_radii_array = None
        original_vobject.electronegativity_array = None
        original_vobject.index_bonds = None
        original_vobject.bonds = None
        original_vobject.non_bonded_atoms = None
        _reapply_manual_bonds ( original_vobject )

        original_vobject.create_representation ( rep_type = "sticks" )
        original_vobject.create_representation ( rep_type = "stick_spheres" )
        _activate_new_sphere_representation ( original_vobject, "stick_spheres" )
        original_vobject.create_representation ( rep_type = "nonbonded" )
        original_vobject.core_representations["picking_dots"] = None
        original_vobject.core_representations["picking_text"] = None

        # [EN] 2026-10-04: "incorporate with DYFF" on closing the Builder --
        # formal charges, type overrides, edited DYFF parameters and partial
        # charges go with the fold-back (before the rebuild: the pDynamo
        # atoms take their formal charges from the vobject)
        if carry_customizations:
            from gui.windows.builder.dyff_parameters import carry_builder_customizations
            carry_builder_customizations ( vismol_object, original_vobject )

        new_system = _build_pdynamo_system_from_vismol_object ( original_vobject, label = original_vobject.name )
        _swap_rebuilt_system_into_psystem ( p_session, source_e_id, new_system )
        final_vobject = original_vobject

        # The temp clone was only ever a scratch workspace for this
        # session -- now that its edits are folded back onto the
        # original, remove it completely (same bookkeeping
        # discard_builder_object_if_unused() does for a never-promoted
        # blank-canvas object, plus popping its OWN real psystem entry,
        # which that function never had to do).
        temp_e_id = getattr ( vismol_object, "e_id", None )
        if temp_e_id is not None:
            p_session.psystem.pop ( temp_e_id, None )
            # [EN] 2026-10-04 BUG FIX: the clone was the ACTIVE system while it
            # was edited -- leaving active_id on its removed e_id made the
            # treeview/statusbar fail (KeyError) right after the fold-back
            if getattr ( p_session, "active_id", None ) == temp_e_id:
                p_session.active_id = source_e_id
        index = getattr ( vismol_object, "index", None )
        name  = getattr ( vismol_object, "name", None )
        if index is not None:
            vm_session.vm_objects_dic.pop ( index, None )
        if name is not None:
            vm_session.vobject_names.pop ( name, None )

    else:
        # Atom count changed (or the original vanished in the meantime) --
        # leave the original untouched and keep the temp clone exactly as
        # a normal, separate, permanent system from now on.
        vismol_object.builder_edit_source_e_id         = None
        vismol_object.builder_edit_original_atom_count = None
        final_vobject = vismol_object

    main.main_treeview.refresh ( )
    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.queue_draw ( )
    return final_vobject
