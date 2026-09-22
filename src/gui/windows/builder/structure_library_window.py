#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  structure_library_window.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  Structure Library browser -- small utility window listing every ".mol2"
#  structure found under gui/windows/builder/structures/<category>/ (see
#  structure_library.discover_structure_files()). Picking one (double-click,
#  or select + "Use") loads it (structure_library.load_structure()), stores
#  the result in vm_session.builder_selected_structure, and switches the
#  Builder to the "Add Structure" tool (vm_session.builder_tool =
#  "add_structure") so the very next click anywhere places it -- see
#  click_mode.handle_click_to_add_structure().
#
#  Mirrors fragment_library_window.py's exact open_window()/close_window()
#  self-contained Glade+Python window pattern (own separate module/.glade
#  file rather than parameterising that one -- see structure_library.py's
#  own top-of-file note for why the two features are genuinely distinct,
#  not just a config toggle on the same one).
# ============================================================================
import os
import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import Gtk

from gui.windows.builder.structure_library import discover_structure_files, load_structure, StructureError
from gui.windows.builder.molecule_preview import MoleculePreview


class StructureLibraryWindow ( ):
    """ Structure Library browser. """

    def __init__ ( self, main = None ):
        """ Class initialiser """
        self.main       = main
        self.vm_session = main.vm_session
        self.visible    = False

    def open_window ( self ):
        """ Opens the browser, (re)populating the tree from whatever
        structure files currently exist on disk -- re-scanning on every
        open (rather than caching), same reasoning as fragment_library_
        window.py's own open_window(): a structure dropped into the folder
        while EasyHybrid is already running shows up the next time this
        window is opened, no restart needed. """
        if self.visible:
            self.window.present ( )
            return

        self.builder = Gtk.Builder ( )
        self.builder.add_from_file ( os.path.join ( self.main.home, 'src/gui/windows/builder/structure_library_window.glade' ) )
        self.builder.connect_signals ( self )

        self.window        = self.builder.get_object ( 'window' )
        self.treestore      = self.builder.get_object ( 'structure_treestore' )
        self.treeview        = self.builder.get_object ( 'structure_treeview' )
        self.status_label   = self.builder.get_object ( 'status_label' )
        self.preview_box    = self.builder.get_object ( 'preview_box' )

        # [EN] User's own request -- same small, disposable 3D preview as
        # fragment_library_window.py's own open_window() (see molecule_
        # preview.py's own top-of-file note) -- a fresh MoleculePreview
        # every time this window opens.
        self.preview = MoleculePreview ( )
        self.preview_box.pack_start ( self.preview.widget, True, True, 0 )
        self.preview.widget.show ( )

        self._populate_tree ( )

        self.window.show_all ( )
        self.visible = True

    def _populate_tree ( self ):
        """ Fills self.treestore from discover_structure_files() -- one
        top-level (expanded) row per category folder, one child row per
        ".mol2" file inside it. Category rows carry an empty string in
        the (hidden) path column, which on_structure_row_activated()/
        on_use_button_clicked() both check to ignore a category-row
        "selection" instead of trying to load a folder as a structure. """
        self.treestore.clear ( )
        for category_name, entries in discover_structure_files ( ):
            category_iter = self.treestore.append ( None, [ category_name, "" ] )
            for display_name, file_path in entries:
                self.treestore.append ( category_iter, [ display_name, file_path ] )
        self.treeview.expand_all ( )

    def _load_and_select ( self, file_path ):
        """ Shared by on_structure_row_activated()/on_use_button_clicked():
        loads file_path via structure_library.load_structure(), and on
        success stores it as the Builder's currently-selected structure
        and switches to the "Add Structure" tool. On failure (a malformed
        structure file -- see StructureError's own callers in structure_
        library.py for the specific, actionable messages), shows the error
        in this window's own status label instead of silently doing
        nothing. """
        try:
            structure = load_structure ( file_path )
        except StructureError as e:
            self.status_label.set_text ( "Error: {}".format ( e ) )
            return

        self.vm_session.builder_selected_structure = structure
        self.vm_session.builder_tool = "add_structure"

        # [EN] Keeps the Builder sidebar's own Tool radio group visually in
        # sync, if that window happens to be open -- same "vm_session.
        # builder_tool is the actual source of truth either way" reasoning
        # as fragment_library_window.py's own _load_and_select().
        sidebar = getattr ( self.main, "builder_sidebar_window", None )
        if sidebar is not None and getattr ( sidebar, "visible", False ):
            radio = getattr ( sidebar, "tool_add_structure_radio", None )
            if radio is not None:
                radio.set_active ( True )

        self.status_label.set_text ( "Selected: {} ({} atoms) -- click anywhere in the Builder to place it.".format (
                structure["name"], len ( structure["atoms"] ) ) )

    # ------------------------------------------------------------------
    #  Signal handlers (referenced by name in structure_library_window.glade)
    # ------------------------------------------------------------------

    def on_selection_changed ( self, selection ):
        """ Live preview while just BROWSING the list -- same reasoning as
        fragment_library_window.py's own on_selection_changed(). Never
        touches vm_session.builder_selected_structure/builder_tool (pure
        visual feedback); a malformed structure file or a category row
        just clears the preview, the status label stays the place a real
        Use/double-click attempt reports StructureError details. """
        model, tree_iter = selection.get_selected ( )
        if tree_iter is None:
            self.preview.clear ( )
            return
        file_path = model[tree_iter][1]
        if not file_path:
            self.preview.clear ( )
            return
        try:
            structure = load_structure ( file_path )
        except StructureError:
            self.preview.clear ( )
            return
        self.preview.show ( structure["name"], structure["atoms"], structure["bonds"] )

    def on_structure_row_activated ( self, treeview, path, column ):
        """ Double-click (or Enter on a selected row). """
        file_path = self.treestore[path][1]
        if file_path:
            self._load_and_select ( file_path )

    def on_use_button_clicked ( self, button ):
        """ Same as double-clicking a row, for whichever row is
        currently selected (if any). """
        selection = self.treeview.get_selection ( )
        model, tree_iter = selection.get_selected ( )
        if tree_iter is None:
            self.status_label.set_text ( "Select a structure first." )
            return
        file_path = model[tree_iter][1]
        if not file_path:
            self.status_label.set_text ( "Select a structure, not a category." )
            return
        self._load_and_select ( file_path )

    def on_close_button_clicked ( self, *args ):
        """ Accepts *args so it can be connected directly to BOTH the
        "Close" button's `clicked` signal (button) -> and the window's
        own `delete-event` (widget, event), same convention fragment_
        library_window.py's own close handler uses. Deliberately does NOT
        clear vm_session.builder_selected_structure/builder_tool --
        closing this browser is just putting it away, not un-selecting
        whatever structure the user already chose (they may still want to
        place several copies of it before picking a different one). """
        if getattr ( self, "window", None ) is not None:
            self.window.destroy ( )   # cascades to the packed preview widget/GL context too
        self.preview  = None
        self.visible = False
        return True
