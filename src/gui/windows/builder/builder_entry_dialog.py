#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  builder_entry_dialog.py
#
#  Copyright 2022-2026 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
# ============================================================================
#  Entry point for opening the Builder from a GENERIC trigger -- the
#  toolbar button (main_window.py's own button-press handler, "_show_cell"
#  toolbutton) that, unlike the treeview's own per-row "Edit in Builder"
#  menu item (treeview_menu.py's _menu_edit_in_builder()), doesn't already
#  know which molecule (if any) the user wants to work on.
#
#  Per the user's own explicit request: ask "New system" vs "Edit an
#  existing system" FIRST, rather than always defaulting to a blank
#  canvas -- choose_and_open_builder() below is the whole of that. It's a
#  plain, throwaway Gtk.Dialog built directly in Python (no separate
#  .glade file -- matches how this codebase already builds other simple,
#  one-off choice dialogs, e.g. main_window.py's own Gtk.ColorChooserDialog/
#  Gtk.FileChooserDialog usage), reusing the EXISTING SystemComboBox widget
#  (gui/widgets/custom_widgets_mod.py) to list currently-loaded systems --
#  no new ListStore/model wiring needed.
#
#  Deliberately does NOT duplicate the "this system already has a force
#  field/potential -- it'll be lost" warning here: that lives entirely in
#  empty_object.begin_editing_existing_system() (this module's own
#  caller for the "edit existing" branch below), so it fires no matter
#  which entry point (this dialog, or the treeview's row menu) led there --
#  single source of truth.
# ============================================================================
import gi
gi.require_version ( "Gtk", "3.0" )
from gi.repository import Gtk


def choose_and_open_builder ( main ):
    """ Shows the New/Edit chooser and acts on it. If no system is
    currently loaded at all, there's nothing to choose between -- skips
    the dialog entirely and just opens a blank canvas, same as this
    toolbar button's old, unconditional behaviour. """
    vm_session = main.vm_session
    p_session  = main.p_session

    existing_ids = [ e_id for e_id, system in p_session.psystem.items ( ) if system is not None ]

    if not existing_ids:
        main.builder_sidebar_window.open_window ( )
        return

    dialog = Gtk.Dialog ( title = "Open Builder", transient_for = main.window, modal = True )
    dialog.add_buttons ( Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OK, Gtk.ResponseType.OK )
    dialog.set_default_response ( Gtk.ResponseType.OK )

    content = dialog.get_content_area ( )
    content.set_spacing ( 8 )
    content.set_border_width ( 10 )

    new_radio  = Gtk.RadioButton.new_with_label ( None, "New (blank) system" )
    edit_radio = Gtk.RadioButton.new_with_label_from_widget ( new_radio, "Edit an existing system" )
    content.pack_start ( new_radio,  False, False, 0 )
    content.pack_start ( edit_radio, False, False, 0 )

    from gui.widgets.custom_widgets import SystemComboBox
    system_combo = SystemComboBox ( main = main )
    system_combo.set_active ( 0 )
    system_combo.set_sensitive ( False )
    content.pack_start ( system_combo, False, False, 0 )

    def _on_mode_toggled ( button ):
        system_combo.set_sensitive ( edit_radio.get_active ( ) )
    new_radio.connect  ( "toggled", _on_mode_toggled )
    edit_radio.connect ( "toggled", _on_mode_toggled )

    dialog.show_all ( )
    response = dialog.run ( )

    edit_existing = edit_radio.get_active ( )
    e_id          = system_combo.get_system_id ( ) if edit_existing else None
    dialog.destroy ( )

    if response != Gtk.ResponseType.OK:
        return

    if edit_existing:
        if e_id is None:
            main.bottom_notebook.status_teeview_add_new_item (
                    message = 'Edit in Builder: no system selected.', system = None )
            return
        from gui.windows.builder.empty_object import begin_editing_existing_system
        begin_editing_existing_system ( vm_session, e_id )
    else:
        main.builder_sidebar_window.open_window ( )
