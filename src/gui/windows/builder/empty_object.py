#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
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
#      bonds/positions -- called once here, at creation (see
#      create_empty_vismol_object() below), and again by atom_ops.
#      adjust_hydrogens() every time it runs (which itself already runs
#      after every structural edit -- place/replace an atom, a bond-drag,
#      cycling a bond's order, deleting an atom/bond -- see each of
#      those call sites' own comments in click_mode.py/vismol_glcore.py).
#      Same "always fully rebuild rather than incrementally patch"
#      philosophy already used throughout atom_ops.py for the
#      VismolObject's own bonds/topology (_reapply_manual_bonds()) --
#      NOT wired into move_atom() (called on every single mouse-motion
#      event during a live drag), which would rebuild a whole pDynamo
#      System dozens of times a second for no benefit.
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
#      this fix. Now that sync_pdynamo_system() gives the object a REAL
#      e_id from the moment it's created, that -1/root-level fallback
#      path is only ever exercised in the (should no longer happen, but
#      harmless if it somehow does) case where system creation itself
#      fails -- see sync_pdynamo_system()'s own try/except.
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

    if show_molecule:
        vismol_object.create_representation ( rep_type = "lines" )
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
    atoms are added. """
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

    sync_pdynamo_system ( vismol_object )

    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.queue_draw ( )

    return vismol_object


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

    Returns the new System (not yet registered with any eSession/
    treeview -- see sync_pdynamo_system(), which calls this and then
    handles registration). """
    from pMolecule             import Atom, Bond, BondType, Connectivity, ConvertInputConnectivity, System
    from pScientific            import PeriodicTable
    from pScientific.Geometry3  import Coordinates3

    connectivity = Connectivity ( )

    for atom_id in sorted ( vismol_object.atoms.keys ( ) ):
        atom = vismol_object.atoms[atom_id]
        atomic_number = PeriodicTable.AtomicNumber ( atom.symbol )
        connectivity.AddNode ( Atom.WithOptions ( atomicNumber = atomic_number, label = atom.name ) )

    bond_type_by_order = { 1: BondType.Single, 2: BondType.Double, 3: BondType.Triple }
    manual_bonds       = getattr ( vismol_object, "manual_bonds", None ) or set ( )
    manual_bond_orders = getattr ( vismol_object, "manual_bond_orders", None ) or { }

    for ( i, j ) in manual_bonds:
        order     = manual_bond_orders.get ( ( i, j ), 1 )
        bond_type = bond_type_by_order.get ( order, BondType.Single )
        connectivity.AddEdge ( Bond.WithNodes ( connectivity.nodes[i], connectivity.nodes[j],
                                                 isAromatic = False, type = bond_type ) )

    ConvertInputConnectivity ( connectivity, { } )

    system       = System.FromConnectivity ( connectivity = connectivity )
    system.label = label if label else vismol_object.name

    n_atoms      = len ( vismol_object.atoms )
    coordinates3 = Coordinates3.WithExtent ( n_atoms )
    for atom_id in range ( n_atoms ):
        pos = vismol_object.frames[0, atom_id]
        coordinates3[atom_id, 0] = float ( pos[0] )
        coordinates3[atom_id, 1] = float ( pos[1] )
        coordinates3[atom_id, 2] = float ( pos[2] )
    system.coordinates3 = coordinates3

    return system


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
    try:
        main       = vismol_object.vm_session.main
        p_session  = main.p_session

        new_system = _build_pdynamo_system_from_vismol_object ( vismol_object, label = vismol_object.name )

        existing_e_id = getattr ( vismol_object, "e_id", None )

        if existing_e_id is not None and existing_e_id in p_session.psystem:
            old_system = p_session.psystem[existing_e_id]
            for attr_name, attr_value in old_system.__dict__.items ( ):
                if attr_name.startswith ( "e_" ) and not hasattr ( new_system, attr_name ):
                    setattr ( new_system, attr_name, attr_value )
            new_system.e_id = existing_e_id
            p_session.psystem[existing_e_id] = new_system
        else:
            p_session.add_new_system_to_psession ( system = new_system, name = vismol_object.name )
            vismol_object.e_id            = new_system.e_id
            vismol_object.is_builder_only = False

        main.main_treeview.refresh ( )

    except Exception as exc:
        dprint ( "WARNING empty_object.sync_pdynamo_system: failed to (re)build/register the linked "
                "pDynamo system for '{}' -- Builder editing continues unaffected. Error: {}".format (
                getattr ( vismol_object, "name", "?" ), exc ) )
