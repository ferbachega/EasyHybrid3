#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  transform_selection_window.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  "Transform Selection" window -- user's own explicit request: "vamos
#  fazer uma ferramenta para transladar e rotacionar (x, y, z) uma selecao
#  de atomos (usando o viewing selection." Live 6-slider (translate X/Y/Z +
#  rotate X/Y/Z) manipulation of the app's own existing "viewing selection"
#  mechanism (shift+click / shift+drag -- the SAME mechanism [[project_
#  atom_types_window]] already reuses, confirmed independent of any active
#  Builder click tool), scoped to whichever selected atoms belong to the
#  current Builder target object -- same filtering builder_sidebar.py's own
#  on_atom_types_button_clicked() already does.
#
#  User's own explicit design choices (asked via AskUserQuestion before
#  building): LIVE sliders (not numeric-entry-plus-Apply), rotating around
#  the selection's own CENTROID (not the object's origin).
#
#  Mirrors atom_types_window.py's open_window()/close_window() pattern,
#  and builder_sidebar.py's own dihedral-slider live-update pattern
#  (arm/disarm + suppress-signal-during-programmatic-set_value + sync only
#  on button-release) -- see this module's own handlers below for exactly
#  how that pattern is adapted from "one absolute angle, re-measured every
#  time" (dihedral) to "6 independent nudge sliders, delta-from-last-value,
#  snapped back to 0 on release" (translate/rotate has no natural absolute
#  value to re-measure the way a 4-atom dihedral does).
# ============================================================================
import os
import math
import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import Gtk

from gui.windows.builder.atom_ops import (
    translate_atoms, rotate_atoms_around_point, selection_centroid, push_undo_snapshot,
)

_AXIS_VECTORS = { "x": ( 1.0, 0.0, 0.0 ), "y": ( 0.0, 1.0, 0.0 ), "z": ( 0.0, 0.0, 1.0 ) }


class TransformSelectionWindow ( ):
    """ Transform Selection (translate/rotate) window. """

    def __init__ ( self, main = None ):
        """ Class initialiser """
        self.main       = main
        self.vm_session = main.vm_session
        self.visible    = False
        self.target_object = None
        self.atom_ids      = [ ]
        # [EN] "delta from last value" state for the 6 sliders -- see
        # on_slider_released()'s own docstring for why this (not a
        # dihedral-style "re-measure an absolute value") is the right
        # convention here.
        self._last_slider_values = { "tx": 0.0, "ty": 0.0, "tz": 0.0, "rx": 0.0, "ry": 0.0, "rz": 0.0 }
        self._suppress_signal    = False
        self._undo_pushed_this_gesture = False

    def open_window ( self, target_object, atom_ids ):
        """ (Re)opens the window scoped to `atom_ids` (a list/set of
        atom_id ints, ALREADY filtered by the caller to atoms belonging to
        `target_object` -- see builder_sidebar.py's own on_transform_
        selection_button_clicked()). Re-captures from scratch on every
        call, same "re-scan rather than trust stale state" philosophy as
        atom_types_window.py's own open_window(). """
        self.target_object = target_object
        self.atom_ids      = sorted ( atom_ids )
        self._reset_sliders ( )

        if not self.visible:
            self.builder = Gtk.Builder ( )
            self.builder.add_from_file ( os.path.join ( self.main.home, 'src/gui/windows/builder/transform_selection_window.glade' ) )
            self.builder.connect_signals ( self )

            self.window       = self.builder.get_object ( 'window' )
            self.header_label = self.builder.get_object ( 'header_label' )
            self.status_label = self.builder.get_object ( 'status_label' )
            self.scales = {
                "tx": self.builder.get_object ( 'translate_x_scale' ),
                "ty": self.builder.get_object ( 'translate_y_scale' ),
                "tz": self.builder.get_object ( 'translate_z_scale' ),
                "rx": self.builder.get_object ( 'rotate_x_scale' ),
                "ry": self.builder.get_object ( 'rotate_y_scale' ),
                "rz": self.builder.get_object ( 'rotate_z_scale' ),
            }

            self.window.show_all ( )
            self.visible = True
        else:
            self.window.present ( )

        for scale in self.scales.values ( ):
            scale.set_sensitive ( len ( self.atom_ids ) > 0 )
        self.header_label.set_text ( "{} -- {} atom(s) selected".format (
                self.target_object.name, len ( self.atom_ids ) ) )
        self.status_label.set_text ( "" )

    def _reset_sliders ( self ):
        """ Snaps every slider back to 0 (suppressing the value-changed
        signal this itself would otherwise trigger -- same guard convention
        builder_sidebar.py's own arm_dihedral_slider() uses) and resets the
        delta-tracking state to match. Called both when the window (re)opens
        and after every drag-release (on_slider_released() below). """
        self._last_slider_values = { key: 0.0 for key in self._last_slider_values }
        self._undo_pushed_this_gesture = False
        scales = getattr ( self, "scales", None )
        if not scales:
            return
        self._suppress_signal = True
        for scale in scales.values ( ):
            scale.set_value ( 0.0 )
        self._suppress_signal = False

    def _live_atom_ids ( self ):
        """ Filters self.atom_ids down to atoms still present on target_
        object -- same defensive "an earlier atom_id may have since been
        removed by Undo/the Delete tool" guard every other Builder click
        handler with a captured atom-id list already uses (e.g. click_mode.
        handle_click_to_pick_dihedral_atom()'s own stale-reference check). """
        if self.target_object is None:
            return [ ]
        return [ aid for aid in self.atom_ids if aid in self.target_object.atoms ]

    def _apply_delta ( self, key, new_value ):
        """ Shared by all 6 on_*_changed() handlers below: computes the
        delta from this slider's own LAST value, applies it (translate:
        plain per-axis offset; rotate: around the selection's own centroid,
        recomputed FRESH every call -- see atom_ops.selection_centroid()'s
        own docstring for why that matters across interleaved gestures),
        and remembers new_value as the new "last value" for the next
        delta. Pushes ONE undo snapshot for the whole gesture (the FIRST
        real delta applied since the last reset), matching this project's
        established "one snapshot per user gesture, not per slider tick"
        convention (e.g. the dihedral tool's own push_undo_snapshot() at
        arm-time, or clean_up_structure()'s at button-click time). """
        if self._suppress_signal:
            return
        atom_ids = self._live_atom_ids ( )
        if not atom_ids:
            return

        delta = new_value - self._last_slider_values[key]
        self._last_slider_values[key] = new_value
        if abs ( delta ) < 1e-9:
            return

        if not self._undo_pushed_this_gesture:
            push_undo_snapshot ( self.target_object )
            self._undo_pushed_this_gesture = True

        if key == "tx":
            translate_atoms ( self.target_object, atom_ids, delta, 0.0, 0.0 )
        elif key == "ty":
            translate_atoms ( self.target_object, atom_ids, 0.0, delta, 0.0 )
        elif key == "tz":
            translate_atoms ( self.target_object, atom_ids, 0.0, 0.0, delta )
        else:
            axis  = _AXIS_VECTORS[ key[1] ]   # "rx"/"ry"/"rz" -> "x"/"y"/"z"
            pivot = selection_centroid ( self.target_object, atom_ids )
            rotate_atoms_around_point ( self.target_object, atom_ids, pivot, axis, math.radians ( delta ) )

    # ------------------------------------------------------------------
    #  Signal handlers (referenced by name in transform_selection_window.glade)
    # ------------------------------------------------------------------

    def on_translate_x_changed ( self, scale ):
        self._apply_delta ( "tx", scale.get_value ( ) )

    def on_translate_y_changed ( self, scale ):
        self._apply_delta ( "ty", scale.get_value ( ) )

    def on_translate_z_changed ( self, scale ):
        self._apply_delta ( "tz", scale.get_value ( ) )

    def on_rotate_x_changed ( self, scale ):
        self._apply_delta ( "rx", scale.get_value ( ) )

    def on_rotate_y_changed ( self, scale ):
        self._apply_delta ( "ry", scale.get_value ( ) )

    def on_rotate_z_changed ( self, scale ):
        self._apply_delta ( "rz", scale.get_value ( ) )

    def on_slider_released ( self, widget, event ):
        """ Wired to EVERY one of the 6 sliders' own button-release-event
        (same widget, same handler -- only one of them will actually have
        just been dragged). Syncs the linked pDynamo system exactly once,
        at gesture END -- NOT on every value-changed event during the drag
        (matching this project's own "don't rebuild a whole pDynamo System
        per mouse-motion event" rule, see empty_object.py's module
        docstring) -- then snaps every slider back to 0 so the next drag,
        on ANY axis, starts from a clean baseline instead of accumulating
        toward each slider's own +/-5 A / +/-180 degree range limit.
        Returns False so GTK's own default GtkScale button-release
        handling still runs. """
        if self.target_object is not None:
            from gui.windows.builder.empty_object import sync_pdynamo_system
            sync_pdynamo_system ( self.target_object )
        self._reset_sliders ( )
        return False

    def on_refresh_button_clicked ( self, button ):
        """ Re-captures the app's CURRENT viewing selection (in case it
        changed since this window opened), filtered to this Builder's own
        target object -- same filtering on_transform_selection_button_
        clicked() does the first time. """
        target_object = getattr ( self.vm_session, "builder_target_object", None )
        if target_object is None:
            self.status_label.set_text ( "No Builder target object." )
            return
        selection = self.vm_session.selections.get ( self.vm_session.current_selection )
        selected_atoms = getattr ( selection, "selected_atoms", None ) or set ( )
        atom_ids = { atom.atom_id for atom in selected_atoms if atom.vm_object is target_object }
        self.open_window ( target_object, atom_ids )
        self.status_label.set_text ( "Refreshed." )

    def on_close_button_clicked ( self, *args ):
        """ Accepts *args so it can be connected directly to BOTH the
        "Close" button's `clicked` signal (button) -> and the window's own
        `delete-event` (widget, event) -- same convention every other
        self-contained Builder window in this project already uses. """
        if getattr ( self, "window", None ) is not None:
            self.window.destroy ( )
        self.visible = False
        return True
