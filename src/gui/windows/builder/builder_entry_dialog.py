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


def _close_any_open_builder_session ( main ):
    """ [EN] 2026-09-25 BUG FIX (user's own report: "se crio um sistema,
    e depois tendo criar outro do zero, nao funciona" -- reproduced live:
    clicking this toolbar button's "New (blank) system" option again
    WHILE the sidebar is still open from a PREVIOUS session doesn't start
    anything fresh at all -- BuilderSidebarWindow.open_window()'s own
    first line is `if self.visible: self.window.present(); return`, so it
    just re-raises the SAME still-open window, silently resuming the
    SAME vismol_object/pDynamo system -- any further drawing just adds
    MORE atoms to that one system instead of creating a new treeview row,
    exactly matching the reported symptom (molecule visibly grows in the
    3D view, but no new row ever appears).

    Mirrors empty_object.begin_editing_existing_system()'s OWN existing
    "if the Builder sidebar is already open (mid-edit of something
    else), it's closed first" step (same "one editable object at a time"
    rule) -- that path already got this right; choose_and_open_builder()'s
    OWN two "New" branches did not, until now. Closing first runs the
    real close_window() logic (finalise/promote or discard, whichever
    applies), so nothing from the previous session is silently lost or
    left half-finished -- it just properly ENDS before a new one begins. """
    sidebar = getattr ( main, "builder_sidebar_window", None )
    if sidebar is not None and getattr ( sidebar, "visible", False ):
        sidebar.close_window ( )


def choose_and_open_builder ( main ):
    """ Shows the New/Edit chooser and acts on it. If no system is
    currently loaded at all, there's nothing to choose between -- skips
    the dialog entirely and just opens a blank canvas, same as this
    toolbar button's old, unconditional behaviour. """
    vm_session = main.vm_session
    p_session  = main.p_session

    existing_ids = [ e_id for e_id, system in p_session.psystem.items ( ) if system is not None ]

    if not existing_ids:
        _close_any_open_builder_session ( main )
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

    # [EN] 2026-09-24, user's own explicit request: "o usuario tem que
    # escolher qual vobject tera como referencia" -- a System can have
    # more than one real VismolObject (several docking poses sharing one
    # e_id, for instance), so picking a System alone was ambiguous (see
    # empty_object.begin_editing_existing_system()'s own updated
    # docstring for the bug this used to cause). ALWAYS shown (even for a
    # single-candidate system), per the user's own explicit choice,
    # rather than only appearing when there's a real choice to make.
    label_vobject = Gtk.Label ( label = "Vobject to edit:" )
    label_vobject.set_xalign ( 0 )
    content.pack_start ( label_vobject, False, False, 0 )
    vobject_combo = Gtk.ComboBoxText ( )
    vobject_combo.set_sensitive ( False )
    content.pack_start ( vobject_combo, False, False, 0 )

    def _populate_vobject_combo ( ):
        vobject_combo.remove_all ( )
        e_id = system_combo.get_system_id ( )
        if e_id is None:
            return
        vobjects = main.get_vobjects_for_system ( e_id )
        for vobject in vobjects:
            vobject_combo.append ( str ( vobject.index ), vobject.name )
        if vobjects:
            vobject_combo.set_active ( 0 )

    def _on_mode_toggled ( button ):
        is_edit = edit_radio.get_active ( )
        system_combo.set_sensitive ( is_edit )
        vobject_combo.set_sensitive ( is_edit )
        if is_edit:
            _populate_vobject_combo ( )
    new_radio.connect  ( "toggled", _on_mode_toggled )
    edit_radio.connect ( "toggled", _on_mode_toggled )
    system_combo.connect ( "changed", lambda w: _populate_vobject_combo ( ) if edit_radio.get_active ( ) else None )

    dialog.show_all ( )
    response = dialog.run ( )

    edit_existing = edit_radio.get_active ( )
    e_id          = system_combo.get_system_id ( ) if edit_existing else None
    vobject_index = vobject_combo.get_active_id ( ) if edit_existing else None
    dialog.destroy ( )

    if response != Gtk.ResponseType.OK:
        return

    if edit_existing:
        if e_id is None or vobject_index is None:
            main.bottom_notebook.status_teeview_add_new_item (
                    message = 'Edit in Builder: no system/vobject selected.', system = None )
            return
        vismol_object = vm_session.vm_objects_dic.get ( int ( vobject_index ) )
        if vismol_object is None:
            main.bottom_notebook.status_teeview_add_new_item (
                    message = 'Edit in Builder: selected vobject no longer exists.', system = None )
            return
        # [EN] 2026-09-25 BUG FIX (user's own report: editing a VObject
        # while choosing to keep it as a SEPARATE new system was removing
        # sibling VObjects/properties from the ORIGINAL system -- see
        # empty_object.py's own note on confirm_and_discard_other_
        # vobjects()'s removal): no longer discards this system's OTHER
        # vobjects before editing -- the original system must stay
        # EXACTLY as it was, whichever of the 2 outcomes (fold back onto
        # the original, or keep the edit as its own new system) this
        # session ends up with. See empty_object.begin_editing_existing_
        # system()'s own builder_edit_source_vobject for how the fold-
        # back path still finds the RIGHT vobject to fold onto now that
        # siblings can legitimately still be there.
        from gui.windows.builder.empty_object import begin_editing_existing_system
        begin_editing_existing_system ( vm_session, e_id, vismol_object = vismol_object )
    else:
        _close_any_open_builder_session ( main )
        main.builder_sidebar_window.open_window ( )


def choose_vobject_and_edit ( main, e_id ):
    """ [EN] 2026-09-24 addition -- vobject picker for entry points that
    already know WHICH system to edit (the treeview's own per-system-row
    "Edit in Builder", treeview_menu.py's _menu_edit_in_builder()), so
    only the vobject choice is still ambiguous. ALWAYS shows the picker
    dialog, even when the system currently has just one candidate vobject
    -- the user's own explicit preference (same reasoning as choose_and_
    open_builder()'s own vobject combo above: consistency over skipping
    steps in the common case), rather than silently auto-picking when
    there's "nothing to choose". """
    vm_session = main.vm_session
    vobjects   = main.get_vobjects_for_system ( e_id )
    if not vobjects:
        main.bottom_notebook.status_teeview_add_new_item (
                message = 'Edit in Builder: no editable object for this system.', system = None )
        return

    dialog = Gtk.Dialog ( title = "Edit in Builder", transient_for = main.window, modal = True )
    dialog.add_buttons ( Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OK, Gtk.ResponseType.OK )
    dialog.set_default_response ( Gtk.ResponseType.OK )

    content = dialog.get_content_area ( )
    content.set_spacing ( 8 )
    content.set_border_width ( 10 )

    label = Gtk.Label ( label = "Vobject to edit:" )
    label.set_xalign ( 0 )
    content.pack_start ( label, False, False, 0 )

    combo = Gtk.ComboBoxText ( )
    for vobject in vobjects:
        combo.append ( str ( vobject.index ), vobject.name )
    combo.set_active ( 0 )
    content.pack_start ( combo, False, False, 0 )

    dialog.show_all ( )
    response      = dialog.run ( )
    vobject_index = combo.get_active_id ( )
    dialog.destroy ( )

    if response != Gtk.ResponseType.OK or vobject_index is None:
        return

    vismol_object = vm_session.vm_objects_dic.get ( int ( vobject_index ) )
    if vismol_object is None:
        main.bottom_notebook.status_teeview_add_new_item (
                message = 'Edit in Builder: selected vobject no longer exists.', system = None )
        return

    # [EN] 2026-09-25: no longer discards this system's other vobjects --
    # see the sibling comment in choose_and_open_builder() above.
    from gui.windows.builder.empty_object import begin_editing_existing_system
    begin_editing_existing_system ( vm_session, e_id, vismol_object = vismol_object )
