#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  builder_sidebar.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  Builder editing sidebar -- small utility window with the controls for
#  everything built into gui/windows/builder/{atom_ops,click_mode}.py so
#  far: turning "placemode" (vm_session.builder_atom_mode) on/off,
#  picking the active tool (Add / Delete -- vm_session.builder_tool),
#  picking the element new atoms are placed as (vm_session.
#  builder_atom_symbol), and Undo (atom_ops.undo()).
#
#  [EN] DESIGN CHOICE: opening this window IS what puts EasyHybrid into
#  Builder editing mode (matching the request that this panel's
#  PRESENCE tracks whether placemode is active, rather than being an
#  always-visible panel with its own separate open/closed state) -- if
#  there's no vm_session.builder_target_object yet, one is created
#  automatically (empty_object.create_empty_vismol_object(), the exact
#  same helper the terminal's `new` command already uses). Closing the
#  window turns editing back off. The "Editing: ON/OFF" toggle INSIDE
#  the sidebar is a separate, lighter-weight PAUSE: it stops atom-
#  placement/deletion clicks from doing anything without closing the
#  window or losing your current tool/element selection -- handy for
#  just rotating/inspecting the molecule for a moment mid-edit.
#
#  Mirrors the exact open_window()/close_window() pattern already used
#  by setup/easyhybrid_terminal.py's TerminalWindow (same Gtk.Builder +
#  add_from_file + connect_signals(self) flow), for consistency with the
#  rest of this codebase's auxiliary windows.
# ============================================================================
import os
import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import Gtk

from gui.windows.builder.empty_object import create_empty_vismol_object
from gui.windows.builder.empty_object import discard_builder_object_if_unused
from gui.windows.builder.empty_object import finish_editing_existing_system
from gui.windows.builder.atom_ops     import undo as atom_ops_undo
from gui.windows.builder.atom_ops     import clean_up_structure
from gui.windows.builder.atom_ops     import optimize_geometry_dyff


class BuilderSidebarWindow ( ):
    """ Builder editing sidebar. """

    def __init__ ( self, main = None ):
        """ Class initialiser """
        self.main       = main
        self.vm_session = main.vm_session
        self.visible    = False

    def open_window ( self ):
        """ Opens the sidebar -- and, per this module's own design note
        above, turns Builder editing mode ON, creating a fresh empty
        target object first if none exists yet. """
        if self.visible:
            self.window.present ( )
            return

        self.builder = Gtk.Builder ( )
        self.builder.add_from_file ( os.path.join ( self.main.home, 'src/gui/windows/builder/builder_sidebar.glade' ) )
        self.builder.connect_signals ( self )

        self.window = self.builder.get_object ( 'window' )
        self.window.set_title ( 'Builder' )
        self.window.set_keep_above ( True )

        self.placemode_toggle  = self.builder.get_object ( 'placemode_toggle' )
        self.tool_add_radio    = self.builder.get_object ( 'tool_add_radio' )
        self.tool_delete_radio = self.builder.get_object ( 'tool_delete_radio' )
        self.tool_move_radio   = self.builder.get_object ( 'tool_move_radio' )
        self.tool_bond_order_radio = self.builder.get_object ( 'tool_bond_order_radio' )
        self.tool_attach_fragment_radio = self.builder.get_object ( 'tool_attach_fragment_radio' )
        self.tool_add_structure_radio   = self.builder.get_object ( 'tool_add_structure_radio' )
        self.tool_rotate_dihedral_radio = self.builder.get_object ( 'tool_rotate_dihedral_radio' )
        self.dihedral_angle_scale = self.builder.get_object ( 'dihedral_angle_scale' )
        self.fragments_button = self.builder.get_object ( 'fragments_button' )
        self.structures_button = self.builder.get_object ( 'structures_button' )
        self.element_C_radio   = self.builder.get_object ( 'element_C_radio' )
        self.element_N_radio   = self.builder.get_object ( 'element_N_radio' )
        self.element_O_radio   = self.builder.get_object ( 'element_O_radio' )
        self.element_H_radio   = self.builder.get_object ( 'element_H_radio' )
        self.element_other_radio = self.builder.get_object ( 'element_other_radio' )
        self._other_element_symbol = None
        self.bond_order_single_radio   = self.builder.get_object ( 'bond_order_single_radio' )
        self.bond_order_double_radio   = self.builder.get_object ( 'bond_order_double_radio' )
        self.bond_order_triple_radio   = self.builder.get_object ( 'bond_order_triple_radio' )
        self.bond_order_aromatic_radio = self.builder.get_object ( 'bond_order_aromatic_radio' )
        self.set_bond_order_picking_button = self.builder.get_object ( 'set_bond_order_picking_button' )
        self.transform_selection_button = self.builder.get_object ( 'transform_selection_button' )
        self.undo_button       = self.builder.get_object ( 'undo_button' )
        self.clean_up_button   = self.builder.get_object ( 'clean_up_button' )
        self.clean_up_adjust_hydrogens_checkbutton = self.builder.get_object ( 'clean_up_adjust_hydrogens_checkbutton' )
        self.optimize_dyff_button = self.builder.get_object ( 'optimize_dyff_button' )

        if getattr ( self.vm_session, "builder_target_object", None ) is None:
            vismol_object = create_empty_vismol_object ( self.vm_session, name = "builder_molecule" )
            self.vm_session.builder_target_object = vismol_object

        self.vm_session.builder_atom_mode          = True
        self.vm_session.builder_tool               = "add"
        self.vm_session.builder_atom_symbol        = "C"
        self.vm_session.builder_bond_order         = 1
        self.vm_session.builder_bond_aromatic      = False
        self.vm_session.builder_bond_pick_first_atom = None
        self.vm_session.builder_dihedral_pick_atoms  = [ ]
        self.vm_session.builder_dihedral_axis        = None
        self.vm_session.builder_dihedral_subgroup    = None
        self.vm_session.builder_atom_types_labeled_atoms = [ ]

        # reflete o estado inicial nos widgets sem disparar os handlers
        # de "toggled" de novo (eles ja rodariam com os valores certos,
        # mas sinalizar aqui deixa explicito o que esta acontecendo)
        self.placemode_toggle.set_active ( True )
        self.placemode_toggle.set_label ( "Editing: ON" )
        self.tool_add_radio.set_active ( True )
        self.element_C_radio.set_active ( True )
        self.bond_order_single_radio.set_active ( True )

        self.window.show_all ( )
        self.visible = True

    def close_window ( self, *args ):
        """ Closes the sidebar and turns Builder editing mode back OFF.
        Accepts *args so it can be connected directly to BOTH the
        "Close" button's `clicked` signal (button) -> and the window's
        own `delete-event` (widget, event) -- see builder_sidebar.glade,
        both point at on_close_button_clicked() below, which just calls
        this. Returns True so, when triggered via delete-event (the
        window's own X button), GTK doesn't ALSO try to run its own
        default destroy handling on top of the explicit self.window.
        destroy() already done here.

        [EN] Also discards the Builder's target object if no structural
        edit ever happened during this session (see
        discard_builder_object_if_unused()'s own docstring) -- and, ONLY
        in that case, clears vm_session.builder_target_object back to
        None too, so the *next* open_window() creates a genuinely fresh
        object instead of "resuming" a dangling reference to one that
        was just removed from the session. When a real edit DID happen,
        builder_target_object is deliberately left set, matching this
        module's existing "closing the sidebar pauses editing, reopening
        it resumes the same molecule" design. """
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is not None and getattr ( target_object, "is_builder_only", False ):
            discard_builder_object_if_unused ( target_object )
            self.vm_session.builder_target_object = None
        elif target_object is not None and getattr ( target_object, "builder_edit_source_e_id", None ) is not None:
            # [EN] "Edit in Builder" session on an already-loaded system
            # (see empty_object.begin_editing_existing_system()) -- unlike
            # the blank-canvas case above, this object already has a real,
            # committed e_id (from clone_system()), so there's no "was it
            # ever promoted" question; finish_editing_existing_system()
            # decides whether to fold it back onto the original (atom
            # count unchanged) or keep it as its own new system (atom
            # count changed) -- either way this Builder session is over,
            # so builder_target_object is reset to None unconditionally,
            # not left set for a "reopen resumes" pause like the
            # blank-canvas case.
            finish_editing_existing_system ( target_object )
            self.vm_session.builder_target_object = None

        self.vm_session.builder_atom_mode = False
        self.vm_session.builder_bond_pick_first_atom = None
        self.vm_session.builder_dihedral_pick_atoms  = [ ]
        self.vm_session.builder_dihedral_axis        = None
        self.vm_session.builder_dihedral_subgroup    = None
        # [EN] The Atom Types window's own 3D-view index-label overlay
        # (see atom_types_window.py's _refresh_rows()) points at atoms
        # belonging to the object THIS Builder session was editing --
        # stale/dangling once that session ends, regardless of whether
        # the Atom Types window itself is still open or already closed.
        self.vm_session.builder_atom_types_labeled_atoms = [ ]
        if getattr ( self, "window", None ) is not None:
            self.window.destroy ( )
        self.visible = False
        return True

    # ------------------------------------------------------------------
    #  Signal handlers (referenced by name in builder_sidebar.glade)
    # ------------------------------------------------------------------

    def on_placemode_toggled ( self, button ):
        """ The lighter-weight PAUSE toggle -- see this module's design
        note at the top for how this differs from actually closing the
        window. """
        is_active = button.get_active ( )
        self.vm_session.builder_atom_mode = is_active
        button.set_label ( "Editing: ON" if is_active else "Editing: OFF" )

    def on_tool_changed ( self, button ):
        """ [EN] GtkRadioButton fires "toggled" for BOTH the button that
        just became active AND the one that just became inactive -- only
        act on the one reporting active=True, otherwise this runs TWICE
        per click, the second time with the (now wrong) previous tool. """
        if not button.get_active ( ):
            return
        if button is self.tool_add_radio:
            self.vm_session.builder_tool = "add"
        elif button is self.tool_delete_radio:
            self.vm_session.builder_tool = "delete"
        elif button is self.tool_move_radio:
            self.vm_session.builder_tool = "move"
        elif button is self.tool_bond_order_radio:
            self.vm_session.builder_tool = "bond_order"
        elif button is self.tool_attach_fragment_radio:
            self.vm_session.builder_tool = "attach_fragment"
        elif button is self.tool_add_structure_radio:
            self.vm_session.builder_tool = "add_structure"
        elif button is self.tool_rotate_dihedral_radio:
            self.vm_session.builder_tool = "rotate_dihedral"

        # [EN] Leaving the "Rotate Dihedral" tool (or never having been on
        # it) drops any in-progress 4-atom pick and disarms the slider --
        # same stale-state reasoning as the Bond Order block below, just
        # for this tool's own armed state.
        if self.vm_session.builder_tool != "rotate_dihedral":
            self.vm_session.builder_dihedral_pick_atoms = [ ]
            self.vm_session.builder_dihedral_axis       = None
            self.vm_session.builder_dihedral_subgroup   = None
            self.disarm_dihedral_slider ( )

        # [EN] Leaving the "Bond Order" tool (or never having been on it)
        # drops any pending "first atom armed, waiting for its partner"
        # click_mode.handle_click_to_set_bond_order() may have left behind
        # -- otherwise switching to Add/Delete/Move mid-pick and back
        # could silently resume bonding to an atom the user has long since
        # forgotten they clicked (and, worse, one that clean-up/undo might
        # have since removed -- handle_click_to_set_bond_order() already
        # guards against a STALE reference, but not against a stale
        # SURVIVING one the user simply doesn't remember picking).
        if self.vm_session.builder_tool != "bond_order":
            self.vm_session.builder_bond_pick_first_atom = None

    def on_element_changed ( self, button ):
        """ Same "only act on the newly-active one" reasoning as
        on_tool_changed() above.

        [EN] "Other..." (user's own request -- a button to pick ANY
        element, not just the 4 quick-pick ones here) opens the full
        periodic table picker (periodic_table_dialog.choose_element())
        the moment it becomes the active radio, rather than sitting
        there as its own separate, silent "mode" -- picking a symbol
        sets builder_atom_symbol AND relabels this radio to show what
        was picked (so it stays visible which element is actually
        armed); cancelling with nothing EVER picked yet on this radio
        falls back to Carbon (the same default the sidebar opens with)
        instead of leaving an ambiguous "Other..." selected with no
        real symbol behind it -- cancelling a RE-pick (something was
        already chosen here before) just keeps that previous choice
        active, unchanged. """
        if not button.get_active ( ):
            return

        if button is self.element_other_radio:
            from gui.windows.builder.periodic_table_dialog import choose_element
            symbol = choose_element ( self.main )
            if symbol is not None:
                self._other_element_symbol = symbol
                self.element_other_radio.set_label ( symbol )
                self.vm_session.builder_atom_symbol = symbol
            elif self._other_element_symbol is None:
                self.element_C_radio.set_active ( True )   # nothing chosen yet -- fall back to the default
            return

        symbol_by_button = {
            self.element_C_radio: "C",
            self.element_N_radio: "N",
            self.element_O_radio: "O",
            self.element_H_radio: "H",
        }
        symbol = symbol_by_button.get ( button )
        if symbol is not None:
            self.vm_session.builder_atom_symbol = symbol

    def on_bond_order_changed ( self, button ):
        """ Same "only act on the newly-active one" reasoning as
        on_tool_changed()/on_element_changed() above. Sets the order AND
        aromatic flag together -- this selection applies to bonds
        CREATED from now on (drag-to-bond, 'b' key), to an EXISTING bond
        via Ctrl+click on it (click_mode.apply_selected_bond_order()),
        and to an existing bond via the "Set Bond Order (pk1/pk2)"
        button below (on_set_bond_order_picking_button_clicked()). """
        if not button.get_active ( ):
            return
        if button is self.bond_order_single_radio:
            self.vm_session.builder_bond_order    = 1
            self.vm_session.builder_bond_aromatic = False
        elif button is self.bond_order_double_radio:
            self.vm_session.builder_bond_order    = 2
            self.vm_session.builder_bond_aromatic = False
        elif button is self.bond_order_triple_radio:
            self.vm_session.builder_bond_order    = 3
            self.vm_session.builder_bond_aromatic = False
        elif button is self.bond_order_aromatic_radio:
            self.vm_session.builder_bond_order    = 1
            self.vm_session.builder_bond_aromatic = True

    def on_set_bond_order_picking_button_clicked ( self, button ):
        """ Applies the "Bond order" selection above to the EXISTING bond
        between the 2 atoms currently picked via the measurement picking
        tool (pk1/pk2) -- see click_mode.handle_set_bond_order_picking()
        for the full behaviour (requires pk1/pk2 to already be bonded,
        updates both atoms' hydrogens afterwards). Shows the resulting
        status string (success or why it couldn't be done, e.g. "pick 2
        atoms first") on the main window's status bar, same as how
        terminal commands report their own results. """
        from gui.windows.builder.click_mode import handle_set_bond_order_picking
        msg = handle_set_bond_order_picking ( self.vm_session )
        if getattr ( self.main, "statusbar_main", None ) is not None:
            self.main.statusbar_main.push ( 1, msg )

    def on_fragments_button_clicked ( self, button ):
        """ Opens the Fragment Library browser (fragment_library_window.
        FragmentLibraryWindow, instantiated once on main_window -- see
        main_window.py's own __init__) -- picking a fragment there sets
        vm_session.builder_selected_fragment and switches builder_tool to
        "attach_fragment" (see that window's own _load_and_select()), so
        the next click on a hydrogen atom here attaches it. """
        self.main.fragment_library_window.open_window ( )

    def on_structures_button_clicked ( self, button ):
        """ Opens the Structure Library browser (structure_library_window.
        StructureLibraryWindow, instantiated once on main_window -- see
        main_window.py's own __init__), same pattern as on_fragments_
        button_clicked() above -- picking a structure there sets vm_session.
        builder_selected_structure and switches builder_tool to
        "add_structure" (see that window's own _load_and_select()), so the
        next click ANYWHERE here places it (no hydrogen reference needed,
        unlike Fragments -- see click_mode.handle_click_to_add_structure()'s
        own docstring). """
        self.main.structure_library_window.open_window ( )

    def on_atom_types_button_clicked ( self, button ):
        """ Opens the Atom Types window (atom_types_window.
        AtomTypesWindow, instantiated once on main_window -- same
        pattern as fragment_library_window above) on whichever atoms are
        currently in the app's own VIEWING selection (shift+click/
        shift+drag -- a separate, pre-existing mechanism from any
        Builder click tool, confirmed to already work alongside them
        with no conflict) -- filtered to atoms belonging to THIS
        Builder's own target object, since a viewing selection isn't
        scoped to any one object and could include atoms from something
        else entirely on screen. """
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is None:
            return

        selection = self.vm_session.selections.get ( self.vm_session.current_selection )
        selected_atoms = getattr ( selection, "selected_atoms", None ) or set ( )
        atom_ids = { atom.atom_id for atom in selected_atoms if atom.vm_object is target_object }

        if not atom_ids:
            if getattr ( self.main, "statusbar_main", None ) is not None:
                self.main.statusbar_main.push ( 1, "Atom Types: select some atoms first (shift+click or shift+drag)." )
            return

        self.main.atom_types_window.open_window ( target_object, atom_ids )

    def on_transform_selection_button_clicked ( self, button ):
        """ Opens the Transform Selection window (transform_selection_
        window.TransformSelectionWindow, instantiated once on main_window)
        on whichever atoms are currently in the app's own VIEWING
        selection -- same filtering-to-the-Builder's-own-target-object
        reasoning as on_atom_types_button_clicked() just above. """
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is None:
            return

        selection = self.vm_session.selections.get ( self.vm_session.current_selection )
        selected_atoms = getattr ( selection, "selected_atoms", None ) or set ( )
        atom_ids = { atom.atom_id for atom in selected_atoms if atom.vm_object is target_object }

        if not atom_ids:
            if getattr ( self.main, "statusbar_main", None ) is not None:
                self.main.statusbar_main.push ( 1, "Transform Selection: select some atoms first (shift+click or shift+drag)." )
            return

        self.main.transform_selection_window.open_window ( target_object, atom_ids )

    def arm_dihedral_slider ( self, angle_deg ):
        """ Called by click_mode.handle_click_to_pick_dihedral_atom() (via
        vm_session.main.builder_sidebar_window, the same cross-window-sync
        pattern fragment_library_window.py already uses for THIS sidebar's
        Tool radio) once a full 4-atom dihedral pick succeeds: sets the
        slider to the dihedral's CURRENT measured angle (so dragging away
        from it is relative to where it already is, not some arbitrary
        0) and enables it. The suppress flag stops this programmatic
        set_value() from itself triggering on_dihedral_slider_value_
        changed() below as if the user had dragged it. """
        if self.dihedral_angle_scale is None:
            return
        self._suppress_dihedral_slider_signal = True
        self.dihedral_angle_scale.set_value ( angle_deg )
        self._suppress_dihedral_slider_signal = False
        self.dihedral_angle_scale.set_sensitive ( True )

    def disarm_dihedral_slider ( self ):
        """ Counterpart to arm_dihedral_slider() -- disables the slider
        again (tool switched away, pick cancelled, sidebar closed). Does
        NOT reset its displayed value -- there's nothing meaningful to
        reset it TO once disarmed, and leaving the last angle visible is
        harmless since it's not readable as "live" while disabled. """
        if self.dihedral_angle_scale is None:
            return
        self.dihedral_angle_scale.set_sensitive ( False )

    def on_dihedral_slider_value_changed ( self, scale ):
        """ Rotates the currently-armed dihedral (vm_session.
        builder_dihedral_axis/builder_dihedral_subgroup, set by click_mode.
        handle_click_to_pick_dihedral_atom()) so its CURRENT angle
        (re-measured fresh every call via util.geometric_analysis.
        get_dihedral() -- never trusted/accumulated, so this always means
        "set the dihedral to exactly this value") becomes the slider's new
        value. No pDynamo sync here -- see on_dihedral_slider_released()
        below for why that only happens once, at drag end. """
        if getattr ( self, "_suppress_dihedral_slider_signal", False ):
            return

        target_object = getattr ( self.vm_session, "builder_target_object", None )
        axis          = getattr ( self.vm_session, "builder_dihedral_axis", None )
        subgroup      = getattr ( self.vm_session, "builder_dihedral_subgroup", None )
        picked        = getattr ( self.vm_session, "builder_dihedral_pick_atoms", None )
        if target_object is None or axis is None or subgroup is None or not picked or len ( picked ) < 4:
            return

        from util.geometric_analysis import get_dihedral
        atom_ids = [ a.atom_id for a in picked ]
        try:
            current_angle = get_dihedral ( target_object, *atom_ids )
        except ValueError:
            return

        delta_deg = scale.get_value ( ) - current_angle
        while delta_deg > 180.0:
            delta_deg -= 360.0
        while delta_deg <= -180.0:
            delta_deg += 360.0
        if abs ( delta_deg ) < 1e-9:
            return

        import math
        from gui.windows.builder.atom_ops import rotate_atoms_around_bond
        rotate_atoms_around_bond ( target_object, axis[0], axis[1], subgroup, math.radians ( delta_deg ) )

    def on_dihedral_slider_released ( self, widget, event ):
        """ Syncs the linked pDynamo system exactly once, at drag END --
        NOT on every value-changed event during the drag, matching the
        project's own "don't rebuild a whole pDynamo System per mouse-
        motion event" rule (move_atom() follows the same pattern; see
        empty_object.py's module docstring). Returns False so GTK's own
        default GtkScale button-release handling still runs. """
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is not None:
            from gui.windows.builder.empty_object import sync_pdynamo_system
            sync_pdynamo_system ( target_object )
        return False

    def on_undo_button_clicked ( self, button ):
        """ Calls atom_ops.undo() on the CURRENT target object -- see
        that function's own docstring (and push_undo_snapshot()'s) for
        exactly which actions this can step back through (place/replace
        an atom, a whole drag-to-create-a-bonded-atom gesture, deleting
        an atom or bond, changing a bond's order). """
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is None:
            return
        atom_ops_undo ( target_object )

    def on_clean_up_button_clicked ( self, button ):
        """ Runs atom_ops.clean_up_structure() on the WHOLE current
        target object (atom_ids=None -- see that function's own
        docstring for the difference between whole-molecule and
        localised modes; this button always uses whole-molecule, per
        explicit request). Pushes an undo snapshot first, same as every
        other Builder mutation, since this can move every atom in the
        object at once.

        adjust_hydrogen_count comes straight from the sidebar's own
        checkbox (user's own request) -- when unticked, Clean Up only
        relaxes heavy-atom positions/angles/planarity and leaves
        whichever hydrogens already exist exactly where they are,
        neither repositioned nor added/removed (see clean_up_structure()'s
        own docstring for the full reasoning). """
        from gui.windows.builder.atom_ops import push_undo_snapshot
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is None:
            return
        adjust_hydrogen_count = ( self.clean_up_adjust_hydrogens_checkbutton is None
                                   or self.clean_up_adjust_hydrogens_checkbutton.get_active ( ) )
        push_undo_snapshot ( target_object )
        clean_up_structure ( target_object, atom_ids = None, adjust_hydrogen_count = adjust_hydrogen_count )

        from gui.windows.builder.empty_object import sync_pdynamo_system
        sync_pdynamo_system ( target_object )

    def on_optimize_dyff_button_clicked ( self, button ):
        """ User's own request: "otimizar a geometria com o DYFF, vamos
        criar um sistema provisorio no background, otimizar por 100 passos
        e atualizar as coordenadas, depois podemos descartar o sistema
        provisorio" -- runs atom_ops.optimize_geometry_dyff() (a REAL
        force-field minimisation on a disposable scratch pDynamo system,
        see that function's own docstring) on the WHOLE current target
        object. Pushes an undo snapshot first, same as Clean Up, since
        this can move every atom in the object at once. """
        from gui.windows.builder.atom_ops import push_undo_snapshot
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is None:
            return
        push_undo_snapshot ( target_object )
        ok, message = optimize_geometry_dyff ( target_object )
        if ok:
            self.main.statusbar_main.push ( 1, message )
            from gui.windows.builder.empty_object import sync_pdynamo_system
            sync_pdynamo_system ( target_object )
        else:
            self.main.simple_dialog.error ( msg = message )

    def on_close_button_clicked ( self, *args ):
        return self.close_window ( )
