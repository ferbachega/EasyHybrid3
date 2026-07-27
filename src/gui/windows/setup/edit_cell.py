#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  edit_cell.py
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#  Maintainer:
#      Fernando Bachega <ferbachega@gmail.com> or <easyhybrid3@gmail.com>
#
#  Description:
#      "System > Cell and Symmetry > Edit Cell" window.
#
#      Lets the user inspect/edit the crystal system and cell parameters
#      (a, b, c, alpha, beta, gamma) of a pDynamo system. If the system
#      already has a symmetry defined, the window is pre-filled with its
#      current values. Since NPT simulations change the cell volume over
#      time, the cell parameters can differ between the several loaded
#      objects/frames of a same system -- an "Object / frame" combobox lets
#      the user inspect the cell of each one of them individually.
#
import os

import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox


# . Default parameters used when a system does not have a cell yet and the
#   user turns periodicity on for the first time.
_DEFAULT_CELL = { 'a' : 20.0, 'b' : 20.0, 'c' : 20.0, 'alpha' : 90.0, 'beta' : 90.0, 'gamma' : 90.0 }


class EditCellWindow:
    """ "Edit Cell" window: edits the crystal system/cell parameters of a
    pDynamo system, one object/system (and, within it, one coordinate
    frame) at a time. """

    def __init__(self, main = None):
        """ Class initialiser """
        self.main       = main
        self.vm_session = main.vm_session
        self.Visible    = False
        self.home       = main.home
        self.p_session  = main.p_session

        # . Guard used while the entries are being filled programmatically,
        #   so that the "changed" signals used to enforce crystal-system
        #   constraints don't recurse/overwrite user input.
        self._updating  = False

    #-------------------------------------------------------------------------
    #  W I N D O W   L I F E C Y C L E
    #-------------------------------------------------------------------------
    def open_window (self):
        """ Function doc """
        if self.Visible == False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/setup/edit_cell.glade'))
            self.builder.connect_signals(self)

            self.window = self.builder.get_object('window')
            self.window.set_title('Edit Cell')

            # -------------------- widget shortcuts --------------------
            self.checkbox_periodic  = self.builder.get_object('checkbox_periodic')
            self.box_cell_content   = self.builder.get_object('box_cell_content')
            self.label_info         = self.builder.get_object('label_info')

            self.entry_a     = self.builder.get_object('entry_a')
            self.entry_b     = self.builder.get_object('entry_b')
            self.entry_c     = self.builder.get_object('entry_c')
            self.entry_alpha = self.builder.get_object('entry_alpha')
            self.entry_beta  = self.builder.get_object('entry_beta')
            self.entry_gamma = self.builder.get_object('entry_gamma')

            '''-------------------- systems combobox --------------------'''
            self.box_system = self.builder.get_object('box_system')
            self.combobox_systems = SystemComboBox(self.main)
            self.combobox_systems.connect("changed", self.on_combobox_systemsbox_changed)
            self.box_system.pack_start(self.combobox_systems, False, False, 0)
            '''------------------------------------------------------------'''

            '''-------------------- object/frame combobox --------------------'''
            self.box_coordinates = self.builder.get_object('box_coordinates')
            self.coordinates_combobox = CoordinatesComboBox()
            self.coordinates_combobox.connect("changed", self.on_combobox_coordinates_changed)
            self.box_coordinates.pack_start(self.coordinates_combobox, False, False, 0)
            '''------------------------------------------------------------'''

            '''-------------------- crystal system combobox --------------------'''
            self.box_crystal_system = self.builder.get_object('box_crystal_system')
            self.combobox_crystal_system = Gtk.ComboBoxText()
            for name in self.p_session.get_crystal_system_names():
                self.combobox_crystal_system.append(name, name)
            self.combobox_crystal_system.connect("changed", self.on_combobox_crystal_system_changed)
            self.box_crystal_system.pack_start(self.combobox_crystal_system, False, False, 0)
            '''------------------------------------------------------------'''

            self.window.show_all()
            self.window.connect('destroy', self.close_window)

            # . Start showing the currently active system, if any.
            if self.p_session.psystem:
                self.combobox_systems.set_active_system(e_id = self.p_session.active_id)

            self.Visible = True

    def close_window (self, button = None, data = None):
        """ Function doc """
        self.window.destroy()
        self.Visible = False

    def update (self):
        """ Called by main.uptade_interface_windows_and_dialogs() (e.g. after
        the treeview changes) so the window keeps showing up-to-date data
        for the currently selected system, if it happens to be open. """
        if not self.Visible:
            return
        system_id = self.combobox_systems.get_system_id()
        if system_id is not None and system_id in self.p_session.psystem:
            self.refresh_from_system(system_id)

    #-------------------------------------------------------------------------
    #  S Y S T E M   /   O B J E C T   S E L E C T I O N
    #-------------------------------------------------------------------------
    def on_combobox_systemsbox_changed (self, widget):
        """ Called when the user picks a different system. Refreshes the
        object/frame combobox and shows that system's current cell (if any)."""
        system_id = self.combobox_systems.get_system_id()

        if system_id is not None:
            self.coordinates_combobox.set_model(self.main.vobject_liststore_dict[system_id])
            size = len(list(self.main.vobject_liststore_dict[system_id]))
            self.coordinates_combobox.set_active(size - 1)

            self.refresh_from_system(system_id)

    def on_combobox_coordinates_changed (self, widget):
        """
        Called when the user picks a different object/frame of the current
        system. NPT simulations change the cell volume over time, so each
        loaded frame/object may have its own cell size -- when the selected
        frame has its own displayed cell (vobject.cell_parameters), we show
        those values instead of (only) the ones stored in the live pDynamo
        system, so the user can inspect how the cell evolved.
        """
        if self._updating:
            return

        vobject_id = self.coordinates_combobox.get_vobject_id()
        if vobject_id is None:
            return

        vobject = self.main.vm_session.vm_objects_dic.get(vobject_id, None)
        if vobject is not None and getattr(vobject, 'cell_parameters', None):
            cell = vobject.cell_parameters
            self._set_entries(cell['a'], cell['b'], cell['c'], cell['alpha'], cell['beta'], cell['gamma'])
            self.label_info.set_text('Showing the cell parameters stored for this object/frame.')

    #-------------------------------------------------------------------------
    #  R E A D I N G  /  W R I T I N G   T H E   p D y n a m o   S Y S T E M
    #-------------------------------------------------------------------------
    def refresh_from_system (self, system_id):
        """ Loads the crystal system/cell parameters currently defined for
        "system_id" (if any) into the window. """
        has_symmetry, crystal_system_name, cell = self.p_session.get_cell_and_symmetry(system_e_id = system_id)

        self._updating = True
        try:
            if has_symmetry:
                self.checkbox_periodic.set_active(True)
                self.combobox_crystal_system.set_active_id(crystal_system_name)
                self._set_entries(*cell)
                self.label_info.set_text('Showing the current cell parameters of this system.')
            else:
                self.checkbox_periodic.set_active(False)
                self.combobox_crystal_system.set_active_id('Cubic')
                self._set_entries(_DEFAULT_CELL['a'], _DEFAULT_CELL['b'], _DEFAULT_CELL['c'],
                                   _DEFAULT_CELL['alpha'], _DEFAULT_CELL['beta'], _DEFAULT_CELL['gamma'])
                self.label_info.set_text('This system has no cell parameters yet. '
                                          'Check "Periodic system" to define one.')
        finally:
            self._updating = False

        self.box_cell_content.set_sensitive(self.checkbox_periodic.get_active())
        self.update_editable_fields()

    def _set_entries (self, a, b, c, alpha, beta, gamma):
        """ Function doc """
        self.entry_a.set_text('{:.4f}'.format(a))
        self.entry_b.set_text('{:.4f}'.format(b))
        self.entry_c.set_text('{:.4f}'.format(c))
        self.entry_alpha.set_text('{:.4f}'.format(alpha))
        self.entry_beta.set_text('{:.4f}'.format(beta))
        self.entry_gamma.set_text('{:.4f}'.format(gamma))

    def _get_entries (self):
        """ Function doc """
        try:
            a     = float(self.entry_a.get_text())
            b     = float(self.entry_b.get_text())
            c     = float(self.entry_c.get_text())
            alpha = float(self.entry_alpha.get_text())
            beta  = float(self.entry_beta.get_text())
            gamma = float(self.entry_gamma.get_text())
            return a, b, c, alpha, beta, gamma
        except ValueError:
            return None

    #-------------------------------------------------------------------------
    #  C R Y S T A L   S Y S T E M   C O N S T R A I N T S
    #-------------------------------------------------------------------------
    def update_editable_fields (self):
        """ Enables only the cell entries that are independent for the
        currently selected crystal system; the others are derived
        automatically and shown disabled. """
        crystal_system_name = self.combobox_crystal_system.get_active_id()
        if crystal_system_name is None:
            return

        free_parameters = self.p_session.get_crystal_system_free_parameters(crystal_system_name)

        self.entry_a.set_sensitive('a' in free_parameters)
        self.entry_b.set_sensitive('b' in free_parameters)
        self.entry_c.set_sensitive('c' in free_parameters)
        self.entry_alpha.set_sensitive('alpha' in free_parameters)
        self.entry_beta.set_sensitive('beta' in free_parameters)
        self.entry_gamma.set_sensitive('gamma' in free_parameters)

        self.enforce_constraints()

    def enforce_constraints (self):
        """ Copies the current values of the free parameters into the
        parameters that are dependent for the currently selected crystal
        system (e.g. b = c = a and all angles = 90 for a Cubic cell), so
        that the entries always show a self-consistent cell. """
        crystal_system_name = self.combobox_crystal_system.get_active_id()
        values = self._get_entries()
        if crystal_system_name is None or values is None:
            return

        a, b, c, alpha, beta, gamma = values

        if crystal_system_name == 'Cubic':
            b, c                = a, a
            alpha, beta, gamma  = 90.0, 90.0, 90.0
        elif crystal_system_name == 'Tetragonal':
            b                   = a
            alpha, beta, gamma  = 90.0, 90.0, 90.0
        elif crystal_system_name == 'Orthorhombic':
            alpha, beta, gamma  = 90.0, 90.0, 90.0
        elif crystal_system_name == 'Monoclinic':
            alpha, gamma        = 90.0, 90.0
        elif crystal_system_name == 'Hexagonal':
            b                   = a
            alpha, beta, gamma  = 90.0, 90.0, 120.0
        elif crystal_system_name == 'Rhombohedral':
            b, c                = a, a
            beta, gamma         = alpha, alpha
        # Triclinic: all six parameters are independent -- nothing to force.

        was_updating    = self._updating
        self._updating  = True
        try:
            self._set_entries(a, b, c, alpha, beta, gamma)
        finally:
            self._updating = was_updating

    def on_entry_a_changed (self, widget):
        """ Function doc """
        if self._updating:
            return
        self.enforce_constraints()

    def on_entry_alpha_changed (self, widget):
        """ Function doc """
        if self._updating:
            return
        self.enforce_constraints()

    #-------------------------------------------------------------------------
    #  S I G N A L S
    #-------------------------------------------------------------------------
    def on_combobox_crystal_system_changed (self, widget):
        """ Function doc """
        if self._updating:
            return
        self.update_editable_fields()

    def on_checkbox_periodic_toggled (self, widget):
        """ Function doc """
        self.box_cell_content.set_sensitive(widget.get_active())

    def on_button_cancel_clicked (self, widget):
        """ Function doc """
        self.close_window()

    def on_button_apply_clicked (self, widget):
        """ Function doc """
        self.apply_changes()

    def on_button_ok_clicked (self, widget):
        """ Function doc """
        if self.apply_changes():
            self.close_window()

    #-------------------------------------------------------------------------
    #  A P P L Y
    #-------------------------------------------------------------------------
    def apply_changes (self):
        """ Writes the values currently shown in the window back into the
        selected pDynamo system (creating its symmetry if necessary) and
        refreshes the 3D-viewport cell of the currently selected
        object/frame to match. Returns True on success. """
        system_id = self.combobox_systems.get_system_id()
        if system_id is None:
            self.main.simple_dialog.info(msg = 'Please select a system first.')
            return False

        if not self.checkbox_periodic.get_active():
            # . User does not want a periodic cell for this system:
            #   nothing else to do here (removing an existing symmetry is
            #   not supported by pDynamo, so we simply leave it untouched).
            return True

        crystal_system_name = self.combobox_crystal_system.get_active_id()
        values = self._get_entries()

        if crystal_system_name is None or values is None:
            self.main.simple_dialog.info(msg = 'Please provide valid numeric cell parameters.')
            return False

        a, b, c, alpha, beta, gamma = values

        ok = self.p_session.set_cell_and_symmetry(system_e_id         = system_id,
                                                    crystal_system_name = crystal_system_name,
                                                    a                   = a,
                                                    b                   = b,
                                                    c                   = c,
                                                    alpha               = alpha,
                                                    beta                = beta,
                                                    gamma               = gamma)

        if not ok:
            self.main.simple_dialog.info(msg = 'Could not set the cell parameters for this system.')
            return False

        # . Refresh the 3D-viewport cell box of the currently selected
        #   object/frame, so every object can keep showing its own cell.
        vobject_id = self.coordinates_combobox.get_vobject_id()
        if vobject_id is not None:
            vobject = self.main.vm_session.vm_objects_dic.get(vobject_id, None)
            if vobject is not None:
                vobject.set_cell(a, b, c, alpha, beta, gamma, color = [0.7, 0.7, 0.2])
                self.vm_session.show_cell(vobject)

        self.main.bottom_notebook.status_teeview_add_new_item(
            message = 'Cell parameters updated ({}) for system {}'.format(
                crystal_system_name, self.p_session.psystem[system_id].label))

        self.label_info.set_text('Cell parameters applied.')
        return True
