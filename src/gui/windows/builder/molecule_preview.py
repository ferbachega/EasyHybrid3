#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  molecule_preview.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  Small, self-contained 3D preview -- user's own request: "vamos criar uma
#  instancia do vismol onde, anexado a janela de fragmento, o usuario pode
#  ver uma previa dos fragmentos, com somente a opcao de zoom e rotacao."
#
#  A SEPARATE, independent VismolSession + its own GL context/widget --
#  VismolSession/VismolGTKWidget/VismolGLCore are NOT app-wide singletons:
#  vismol_session.py's own __init__ already builds a fresh VismolGTKWidget/
#  VismolGLCore pair for EVERY VismolSession instance, and this codebase's
#  own (dead) old Builder prototype (src/graphics_engine/src/builder/
#  builder_window.py) already does exactly this -- instantiates a second,
#  completely independent VismolSession. MoleculePreview below is embedded
#  as a plain GTK widget inside a GtkBox placeholder (fragment_library_
#  window.glade / structure_library_window.glade's own "preview_box"), same
#  "self.container.pack_start(self.vm_session.vm_widget, ...)" pattern
#  gui/main/main_window.py itself uses for the app's MAIN 3D view.
#
#  Deliberately NOT the app's EasyHybridSession subclass -- no pDynamo
#  system, no e_id, no treeview row, no Builder editing state: this is
#  pure, disposable, READ-ONLY display of one small molecule at a time.
#  Reuses atom_ops.add_structure_at_position() UNCHANGED (see [[project_
#  builder_structure_library]]) to actually build the displayed
#  VismolObject -- it already does exactly "batch-insert atoms + bonds
#  (with orders/aromaticity) + one representation rebuild", the same thing
#  this preview needs, and neither it nor anything it calls (add_atom(),
#  _reapply_manual_bonds(), _refresh_bond_dependent_representations(),
#  _activate_new_sphere_representation()) ever touches vismol_object.
#  vm_session.main -- confirmed by grep, the ONLY atom_ops.py function that
#  does is optimize_geometry_dyff(), never called here -- so it works
#  unmodified against a bare VismolSession with no `.main` at all.
# ============================================================================
import numpy as np
import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import Gdk

from vismol.core.vismol_session import VismolSession
from vismol.core.vismol_object import VismolObject

from gui.windows.builder.atom_ops import add_structure_at_position


def _install_zoom_and_rotate_only_filter ( vm_glcore ):
    """ [EN] "somente zoom e rotacao": VismolGLCore.mouse_pressed() already
    decides rotate (plain left-drag) / zoom (right-drag, or the SEPARATE
    mouse_scroll() -- untouched here either way) / pan (middle-drag) purely
    from which button/modifier was down at press time -- there is no single
    flag to just turn pan off. Wraps vm_glcore's OWN mouse_pressed method,
    on this ONE instance only (a plain Python attribute override, not a
    change to the shared VismolGLCore class -- the main app's own 3D view,
    and any other VismolSession, is completely unaffected), so that a
    middle-button press (pan) or a shift+left-button press (the multi-
    select box -- see vismol_glcore.py's own mouse_pressed(), "if self.
    shift and not self.ctrl") never reaches the original method at all;
    every other button/modifier combination (plain left = rotate,
    right-drag = zoom, scroll wheel = zoom) is forwarded completely
    unchanged.

    [EN] Why wrap vm_glcore.mouse_pressed rather than the WIDGET's own
    mouse_pressed(widget, event) forwarder (vismol_gtkwidget.py): the
    widget-level method is already bound and connected to the "button-
    press-event" GTK signal INSIDE VismolGTKWidget.__init__ -- reassigning
    widget.mouse_pressed afterwards would NOT change what that already-
    connected signal invokes (GTK signals capture the bound-method
    reference at connect() time). The WIDGET's own mouse_pressed() body,
    by contrast, looks up self.vm_glcore.mouse_pressed FRESH on every
    single event (plain attribute access, not a signal), so overriding it
    on the vm_glcore INSTANCE here takes effect immediately and for every
    future click -- verified live (see this feature's own memory file /
    test harness) that blocked gestures leave mouse_pan/show_selection_box
    False while the base mouse_rotate/mouse_zoom paths still engage
    normally. """
    original_mouse_pressed = vm_glcore.mouse_pressed

    def _filtered_mouse_pressed ( button_number, mouse_x, mouse_y, button_state = None ):
        button = int ( button_number )
        if button == 2:
            return   # middle-drag pan -- swallowed
        if button == 1 and getattr ( vm_glcore, "shift", False ):
            return   # shift+left multi-select box -- swallowed
        return original_mouse_pressed ( button_number, mouse_x, mouse_y, button_state )

    vm_glcore.mouse_pressed = _filtered_mouse_pressed


class MoleculePreview ( ):
    """ Owns one small, disposable VismolSession. `.widget` is the plain
    GTK widget to pack into a placeholder GtkBox in a window's own .glade
    file (e.g. `preview_box.pack_start(preview.widget, True, True, 0)`).
    `.show(name, atoms, bonds)` replaces whatever it is currently showing;
    `.clear()` empties the scene (e.g. when the owning window closes). """

    def __init__ ( self ):
        self.vm_session = VismolSession ( toolkit = "Gtk_3.0" )
        self.widget      = self.vm_session.vm_widget
        _install_zoom_and_rotate_only_filter ( self.vm_session.vm_glcore )

        # [EN] REAL BUG FOUND AND FIXED, 2026-09-24 -- user's own report:
        # "os sticks estao muito finos, quase invisiveis... representacao
        # de ligacoes duplas e triplas [nao aparece]". Root cause: this
        # bare VismolSession's own vismol.core.vismol_config.VismolConfig
        # defaults `sticks_radius` to 0.010 -- 16x thinner than the main
        # app's own gui.config.VismolConfig default (0.16, confirmed by
        # reading both files directly). SAME root cause for BOTH symptoms,
        # not two separate bugs: SticksRepresentation.draw_representation()
        # (representations.py) sets the double/triple-bond PARALLEL-
        # CYLINDER separation as `radius * 1.0` (a GL uniform, `u_sep_loc`)
        # -- i.e. the multi-bond offset is directly proportional to this
        # SAME radius, so at 0.010 the two/three cylinders of a double/
        # triple bond sit only ~0.01 A apart, visually indistinguishable
        # from a single stick even though the underlying bond-order data
        # (already correct -- see atom_ops.add_structure_at_position()/
        # attach_fragment_at_hydrogen(), which this preview reuses
        # unchanged) was fine all along; `multiple_bonds` itself already
        # defaults to True in this same bare config, so nothing there
        # needed changing either. Fixed by overriding just THIS instance's
        # own gl_parameters (not vismol_config.py's shared module-level
        # default, which might be intentional for other bare-session use
        # cases) to match the main app's own proven-visible value.
        self.vm_session.vm_config.gl_parameters['sticks_radius'] = 0.16

    def show ( self, name, atoms, bonds ):
        """ atoms: [(symbol, x, y, z), ...]. bonds: [(i, j, order), ...]
        (fragment_library.load_fragment()'s own shape -- no aromaticity)
        OR [(i, j, order, is_aromatic), ...] (structure_library.
        load_structure()'s shape) -- either is accepted directly, entries
        without an explicit 4th element are treated as non-aromatic. """
        self.clear ( )

        normalized_bonds = [ ( bond[0], bond[1], bond[2], ( bond[3] if len ( bond ) > 3 else False ) )
                              for bond in bonds ]

        vismol_object = VismolObject ( vismol_session = self.vm_session, index = 0, name = name, active = True )
        vismol_object.frames = np.zeros ( ( 1, 0, 3 ), dtype = np.float32 )

        if atoms:
            add_structure_at_position ( vismol_object, { "name": name, "atoms": atoms, "bonds": normalized_bonds },
                                         0.0, 0.0, 0.0 )

        self.vm_session.vm_objects_dic[0] = vismol_object

        vm_glcore = self.vm_session.vm_glcore
        if atoms:
            vm_glcore.center_on_coordinates ( vismol_object, vismol_object.mass_center )
        vm_glcore.queue_draw ( )

    def clear ( self ):
        """ Empties the preview scene -- resets atom_dic_id/atom_id_counter
        too (not just vm_objects_dic), so repeatedly browsing through many
        fragments/structures in one session doesn't leak an ever-growing
        dict of stale Atom references from every previously-shown preview
        (harmless individually, but this window can realistically be
        clicked through dozens of times while someone browses a library). """
        self.vm_session.vm_objects_dic = { }
        self.vm_session.atom_dic_id    = { }
        self.vm_session.atom_id_counter = np.uint32 ( 0 )
        if getattr ( self.vm_session, "vm_glcore", None ) is not None:
            self.vm_session.vm_glcore.queue_draw ( )
