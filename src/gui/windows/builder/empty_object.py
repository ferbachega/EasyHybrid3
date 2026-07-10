#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Molecule Builder -- empty object creation
#
#  Description:
#      First building block of the "Builder" tool (draw molecules from
#      scratch). Creates a VismolObject that is NOT backed by a pDynamo
#      system/psystem entry -- deliberately kept separate from the
#      p_session/QC-MM machinery for now (explicit decision: the builder
#      starts as a pure visual sketchpad; promoting a drawn molecule into
#      a real pDynamo system is a separate, later step).
#
#      Because of that choice, this module does NOT use
#      eSession._add_vismol_object() (which unconditionally looks up
#      self.main.p_session.psystem[vismol_object.e_id], which does not
#      exist for a builder-only object). register_builder_object()
#      below is a deliberately minimal, parallel registration path: it
#      does only the bookkeeping needed for the object to render
#      correctly and show up in session-level object lookups
#      (self.vm_session.vm_objects_dic, used directly by
#      vm_glcore.render(), the terminal's `list`/`show`/`select`
#      commands, etc.).
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
#      this fix. See main_treeview.py's own comments at the fix site for
#      the rest of the story (the -1 sentinel used for the treestore's
#      strictly-int e_id column, root-level placement instead of nesting
#      under a system node, and the matching guard added to
#      on_treeview_mouse_button_release_event() so right-clicking a
#      Builder row doesn't also raise inside treeview_menu.open_menu()).
#      STILL NOT done: a Builder-specific right-click context menu
#      (rename, delete object, promote to a real pDynamo system, ...) --
#      right-click on a Builder row is a silent no-op for now.
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

    if getattr ( vm_session, "vm_glcore", None ) is not None:
        vm_session.vm_glcore.queue_draw ( )

    return vismol_object
