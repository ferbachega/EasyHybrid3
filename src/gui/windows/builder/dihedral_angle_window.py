#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  dihedral_angle_window.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  "Dihedral Angle" window -- user's own explicit request, 2026-09-24: "a
#  edicao do diedro deve ser feita em uma janela a parte, como acontece em
#  'Transform Selections'." Extracted from builder_sidebar.glade/.py's own
#  inline "Dihedral angle" label + GtkScale (part of the same request to
#  slim down/narrow the Builder sidebar) into its own small window,
#  mirroring transform_selection_window.py's own open_window()/
#  close_window() pattern.
#
#  Opened/updated by builder_sidebar.py's arm_dihedral_slider() the moment
#  click_mode.handle_click_to_pick_dihedral_atom() completes a 4-atom
#  dihedral pick (armed with the dihedral's own CURRENT measured angle, so
#  dragging away from it is relative to where it already is, not some
#  arbitrary 0) -- same cross-window-sync mechanism already used
#  everywhere else in this Builder (fragment_library_window.py syncing the
#  sidebar's own Tool radio, etc.). Closed via disarm_dihedral_slider()
#  whenever the pick is cancelled/rejected, the tool is switched away from
#  "Dihedral", or the Builder sidebar itself closes.
#
#  All the actual rotation logic (re-measuring the dihedral's current
#  angle fresh every slider move via util.geometric_analysis.get_dihedral(),
#  never accumulating a delta -- unlike Transform Selection's own 6
#  "nudge, snap to 0" sliders, which have no natural absolute value to
#  re-measure the way a 4-atom dihedral does) is moved here UNCHANGED from
#  builder_sidebar.py's own former on_dihedral_slider_value_changed()/
#  on_dihedral_slider_released().
# ============================================================================
import os
import math
import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import Gtk


class DihedralAngleWindow ( ):
    """ Dihedral Angle window. """

    def __init__ ( self, main = None ):
        """ Class initialiser """
        self.main       = main
        self.vm_session = main.vm_session
        self.visible    = False

    def open_window ( self, angle_deg ):
        """ (Re)opens the window, setting the slider to the dihedral's own
        CURRENT measured angle (`angle_deg`, computed by the caller --
        click_mode.handle_click_to_pick_dihedral_atom() -- right after a
        4-atom pick completes). The suppress flag stops this programmatic
        set_value() from itself triggering on_dihedral_slider_value_
        changed() below as if the user had dragged it. """
        if not self.visible:
            self.builder = Gtk.Builder ( )
            self.builder.add_from_file ( os.path.join ( self.main.home, 'src/gui/windows/builder/dihedral_angle_window.glade' ) )
            self.builder.connect_signals ( self )

            self.window = self.builder.get_object ( 'window' )
            self.dihedral_angle_scale = self.builder.get_object ( 'dihedral_angle_scale' )

            self.window.show_all ( )
            self.visible = True
        else:
            self.window.present ( )

        self._suppress_dihedral_slider_signal = True
        self.dihedral_angle_scale.set_value ( angle_deg )
        self._suppress_dihedral_slider_signal = False

    def close_window ( self, *args ):
        """ Accepts *args so it can be connected directly to BOTH the
        "Close" button's `clicked` signal (button) -> and the window's own
        `delete-event` (widget, event) -- same convention every other
        self-contained Builder window in this project already uses.
        Does NOT disarm vm_session.builder_dihedral_axis/_subgroup/
        _pick_atoms itself -- that stays builder_sidebar.py's own
        on_tool_changed()'s job (switching tools, or this window closing
        because the sidebar itself closed, both already clear that state
        from there), matching how this window never owned that state to
        begin with, only displayed/drove it. """
        if getattr ( self, "window", None ) is not None:
            self.window.destroy ( )
        self.visible = False
        return True

    # ------------------------------------------------------------------
    #  Signal handlers (referenced by name in dihedral_angle_window.glade)
    # ------------------------------------------------------------------

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

    def on_close_button_clicked ( self, *args ):
        return self.close_window ( )
