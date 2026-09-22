#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  fragment_library_window.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  Fragment Library browser -- small utility window listing every ".mol2"
#  fragment found under gui/windows/builder/fragments/<category>/ (see
#  fragment_library.discover_fragment_files()). Picking one (double-click,
#  or select + "Use") loads it (fragment_library.load_fragment()), stores
#  the result in vm_session.builder_selected_fragment, and switches the
#  Builder to the "Attach Fragment" tool (vm_session.builder_tool =
#  "attach_fragment") so the very next click on a hydrogen atom attaches it
#  -- see click_mode.handle_click_to_attach_fragment().
#
#  Mirrors the exact open_window()/close_window() self-contained Glade+
#  Python window pattern already used by builder_sidebar.py.
# ============================================================================
import os
import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import Gtk

from gui.windows.builder.fragment_library import discover_fragment_files, load_fragment, FragmentError
from gui.windows.builder.molecule_preview import MoleculePreview


class FragmentLibraryWindow ( ):
    """ Fragment Library browser. """

    def __init__ ( self, main = None ):
        """ Class initialiser """
        self.main       = main
        self.vm_session = main.vm_session
        self.visible    = False

    def open_window ( self ):
        """ Opens the browser, (re)populating the tree from whatever
        fragment files currently exist on disk -- re-scanning on every
        open (rather than caching) means a fragment dropped into the
        folder structure while EasyHybrid is already running shows up
        the next time this window is opened, no restart needed, matching
        the "easy to customize" goal of the whole feature. """
        if self.visible:
            self.window.present ( )
            return

        self.builder = Gtk.Builder ( )
        self.builder.add_from_file ( os.path.join ( self.main.home, 'src/gui/windows/builder/fragment_library_window.glade' ) )
        self.builder.connect_signals ( self )

        self.window        = self.builder.get_object ( 'window' )
        self.treestore      = self.builder.get_object ( 'fragment_treestore' )
        self.treeview        = self.builder.get_object ( 'fragment_treeview' )
        self.status_label   = self.builder.get_object ( 'status_label' )
        self.preview_box    = self.builder.get_object ( 'preview_box' )

        # [EN] User's own request -- a small, disposable 3D preview (see
        # molecule_preview.py's own top-of-file note), packed into the
        # sized GtkBox placeholder the .glade file already reserves for
        # it. A fresh MoleculePreview every time this window OPENS (not
        # kept alive while closed) -- cheap (a handful of atoms at most)
        # and avoids keeping a second GL context alive for the whole
        # app session just because this browser was opened once.
        self.preview = MoleculePreview ( )
        self.preview_box.pack_start ( self.preview.widget, True, True, 0 )
        self.preview.widget.show ( )

        self._populate_tree ( )

        self.window.show_all ( )
        self.visible = True

    def _populate_tree ( self ):
        """ Fills self.treestore from discover_fragment_files() -- one
        top-level (expanded) row per category folder, one child row per
        ".mol2" file inside it. Category rows carry an empty string in
        the (hidden) path column, which on_fragment_row_activated()/
        on_use_button_clicked() both check to ignore a category-row
        "selection" instead of trying to load a folder as a fragment. """
        self.treestore.clear ( )
        for category_name, entries in discover_fragment_files ( ):
            category_iter = self.treestore.append ( None, [ category_name, "" ] )
            for display_name, file_path in entries:
                self.treestore.append ( category_iter, [ display_name, file_path ] )
        self.treeview.expand_all ( )

    def _load_and_select ( self, file_path ):
        """ Shared by on_fragment_row_activated()/on_use_button_clicked():
        loads file_path via fragment_library.load_fragment(), and on
        success stores it as the Builder's currently-selected fragment
        and switches to the "Attach Fragment" tool. On failure (a
        malformed fragment file -- see FragmentError's own callers in
        fragment_library.py for the specific, actionable messages),
        shows the error in this window's own status label instead of
        silently doing nothing, so a mistake while authoring a NEW
        fragment file is easy to diagnose. """
        try:
            fragment = load_fragment ( file_path )
        except FragmentError as e:
            self.status_label.set_text ( "Error: {}".format ( e ) )
            return

        self.vm_session.builder_selected_fragment = fragment
        self.vm_session.builder_tool = "attach_fragment"

        # [EN] Keeps the Builder sidebar's own Tool radio group visually in
        # sync, if that window happens to be open -- vm_session.builder_tool
        # is the actual source of truth the click-dispatch logic reads (see
        # vismol_glcore.py's render() hook), so the tool works correctly
        # either way; this is purely so the sidebar doesn't show "Add"
        # highlighted while the Builder is actually in fragment-attach mode.
        sidebar = getattr ( self.main, "builder_sidebar_window", None )
        if sidebar is not None and getattr ( sidebar, "visible", False ):
            radio = getattr ( sidebar, "tool_attach_fragment_radio", None )
            if radio is not None:
                radio.set_active ( True )

        self.status_label.set_text ( "Selected: {} -- click a hydrogen atom in the Builder to attach it.".format ( fragment["name"] ) )

    # ------------------------------------------------------------------
    #  Signal handlers (referenced by name in fragment_library_window.glade)
    # ------------------------------------------------------------------

    def on_selection_changed ( self, selection ):
        """ User's own request -- live preview while just BROWSING the
        list (single click / arrow keys), not only after committing to a
        choice via double-click/Use. Deliberately does NOT touch
        vm_session.builder_selected_fragment/builder_tool -- this is pure
        visual feedback, exactly like fragment_library.load_fragment()
        being called here is only ever used to feed the preview, never to
        select the fragment for actual attachment (that stays _load_
        and_select()'s job alone). A malformed fragment file, or a
        category row (empty path), just clears the preview instead of
        showing an error here -- the status label (from a real Use/
        double-click attempt) is still the place that reports FragmentError
        details. """
        model, tree_iter = selection.get_selected ( )
        if tree_iter is None:
            self.preview.clear ( )
            return
        file_path = model[tree_iter][1]
        if not file_path:
            self.preview.clear ( )
            return
        try:
            fragment = load_fragment ( file_path )
        except FragmentError:
            self.preview.clear ( )
            return
        self.preview.show ( fragment["name"], fragment["atoms"], fragment["bonds"] )

    def on_fragment_row_activated ( self, treeview, path, column ):
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
            self.status_label.set_text ( "Select a fragment first." )
            return
        file_path = model[tree_iter][1]
        if not file_path:
            self.status_label.set_text ( "Select a fragment, not a category." )
            return
        self._load_and_select ( file_path )

    def on_close_button_clicked ( self, *args ):
        """ Accepts *args so it can be connected directly to BOTH the
        "Close" button's `clicked` signal (button) -> and the window's
        own `delete-event` (widget, event) -- same convention builder_
        sidebar.py's own close handler already uses. Deliberately does
        NOT clear vm_session.builder_selected_fragment/builder_tool --
        closing this browser is just putting it away, not un-selecting
        whatever fragment the user already chose (they may still want to
        attach several copies of it before picking a different one). """
        if getattr ( self, "window", None ) is not None:
            self.window.destroy ( )   # cascades to the packed preview widget/GL context too
        self.preview  = None
        self.visible = False
        return True
