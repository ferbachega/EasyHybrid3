#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Prepare Packmol System window
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
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
#      "Run Packmol" window -- a generic multi-component packing tool:
#      the user builds a list of components (a structure file on disk +
#      copy count + geometric constraint: inside/outside a box or
#      sphere, or fixed at a given position/orientation), and this
#      window writes the corresponding Packmol .inp file and runs the
#      real `packmol` executable as a subprocess (see
#      util/packmol_runner.py). Deliberately scoped to just produce the
#      packed PDB -- building a topology for it (tleap) is a separate,
#      already-existing step ("Prepare AMBER System (tLeap)"), not
#      something this window chains into automatically.
#
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, GLib

import os

from gui.widgets.custom_widgets import FolderChooserButton
from gui.windows.setup.windows_and_dialogs.dialogs.packmol_component_dialog import AddPackmolComponentDialog
from util import packmol_runner


class PreparePackmolWindow:
    """ "Run Packmol" window. """

    def __init__(self, main=None):
        """ Class initialiser """
        self.main       = main
        self.vm_session = main.vm_session
        self.p_session  = main.p_session
        self.home       = main.home
        self.Visible    = False

        # . The component list's single source of truth -- see
        #   packmol_runner.build_packmol_input()'s docstring for each
        #   entry's shape. liststore_components always mirrors this list
        #   1:1 by row order (rebuilt in full by _refresh_components()
        #   after every mutation, rather than edited in place).
        self.components = []

        self.current_process   = None   # subprocess.Popen or None
        self.current_log_path  = None
        self.current_output_path = None
        self._log_read_offset  = 0
        self._poll_timeout_id  = None

    #-------------------------------------------------------------------------
    #  W I N D O W   L I F E C Y C L E
    #-------------------------------------------------------------------------
    def open_window(self):
        """ Function doc """
        if self.Visible:
            self.window.present()
            return

        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/setup/prepare_packmol_window.glade'))
        self.builder.connect_signals(self)

        self.window = self.builder.get_object('window')
        self.window.set_title('Run Packmol')

        # -------------------- widget shortcuts --------------------
        self.label_info             = self.builder.get_object('label_info')
        self.entry_system_name      = self.builder.get_object('entry_system_name')
        self.entry_packmol_command  = self.builder.get_object('entry_packmol_command')
        self.spinbtn_tolerance      = self.builder.get_object('spinbtn_tolerance')
        self.checkbox_use_seed      = self.builder.get_object('checkbox_use_seed')
        self.spinbtn_seed           = self.builder.get_object('spinbtn_seed')
        self.treeview_components    = self.builder.get_object('treeview_components')
        self.liststore_components   = self.builder.get_object('liststore_components')
        self.textview_log           = self.builder.get_object('textview_log')
        self.button_run             = self.builder.get_object('button_run')

        # -------------------- working-folder chooser --------------------
        self.box_work_folder = self.builder.get_object('box_work_folder')
        self.work_folder_chooser = FolderChooserButton(self.main, 'folder', self.home)
        self.box_work_folder.pack_start(self.work_folder_chooser.btn, False, False, 0)
        scratch = os.environ.get('PDYNAMO3_SCRATCH')
        if scratch and os.path.isdir(scratch):
            self.work_folder_chooser.set_folder(folder=scratch)

        # -------------------- components popup menu --------------------
        # Same "edit-existing-row actions live in a right-click popup,
        # only 'Add' stays a visible button" pattern as
        # prepare_namd_run.py's pipeline treeview.
        self._components_popup_menu = Gtk.Menu()
        for label, handler in (
            ('Edit...', self.on_edit_component_clicked),
            ('Duplicate', self.on_duplicate_component_clicked),
            ('Move up', self.on_move_component_up_clicked),
            ('Move down', self.on_move_component_down_clicked),
            ('Remove', self.on_remove_component_clicked),
            (None, None),
            ('Clear all', self.on_clear_components_clicked),
        ):
            if label is None:
                self._components_popup_menu.append(Gtk.SeparatorMenuItem())
                continue
            menu_item = Gtk.MenuItem(label=label)
            menu_item.connect('activate', handler)
            self._components_popup_menu.append(menu_item)
        self._components_popup_menu.show_all()

        # -------------------- packmol executable --------------------
        packmol_command = self.main.vm_session.vm_config.gl_parameters.get('packmol_command')
        if not packmol_command or not os.path.isfile(packmol_command):
            packmol_command = packmol_runner.find_packmol_executable()
        self.entry_packmol_command.set_text(packmol_command or '')
        if not packmol_command:
            self.label_info.set_text(
                'packmol was not found automatically (checked PATH) -- '
                'please enter the full path to the executable below.')

        self.current_log_path = None
        self._log_read_offset = 0

        self._refresh_components()

        self.window.show_all()
        self.window.connect('destroy', self.close_window)
        self.Visible = True

    def close_window(self, button=None, data=None):
        """ Function doc """
        if not self.Visible:
            return
        self.window.destroy()
        self.Visible = False

    #-------------------------------------------------------------------------
    #  C O M P O N E N T   L I S T
    #-------------------------------------------------------------------------
    def _constraint_display(self, component):
        label = packmol_runner.CONSTRAINT_TYPES[component['constraint_type']]['label']
        params = ' '.join('{:.2f}'.format(p) for p in component['params'])
        return '{} ({})'.format(label, params)

    def _refresh_components(self):
        self.liststore_components.clear()
        for component in self.components:
            self.liststore_components.append([
                os.path.basename(component['structure_path']),
                component['number'],
                self._constraint_display(component),
            ])

    def _selected_component_index(self):
        model, treeiter = self.treeview_components.get_selection().get_selected()
        if treeiter is None:
            return None
        return model.get_path(treeiter).get_indices()[0]

    def _on_components_treeview_button_press(self, widget, event):
        """ Right-click context menu for the component list -- selects
            the row under the pointer first (so a right-click on a row
            the user hadn't already left-clicked still acts on THAT row),
            then pops the menu up there; a right-click on empty space
            leaves the selection alone (Edit/Duplicate/Move/Remove act on
            whatever was already selected, if anything).
        """
        if event.button != 3:
            return False
        path_info = widget.get_path_at_pos(int(event.x), int(event.y))
        if path_info is not None:
            path = path_info[0]
            widget.get_selection().select_path(path)
        self._components_popup_menu.popup(None, None, None, None, event.button, event.time)
        return True

    def on_button_add_component_clicked(self, widget):
        dialog = AddPackmolComponentDialog(main=self.main)
        if dialog.ok:
            self.components.append(dialog.component)
            self._refresh_components()

    def on_edit_component_clicked(self, widget):
        index = self._selected_component_index()
        if index is None:
            return
        dialog = AddPackmolComponentDialog(main=self.main, component=self.components[index])
        if dialog.ok:
            self.components[index] = dialog.component
            self._refresh_components()

    def on_duplicate_component_clicked(self, widget):
        index = self._selected_component_index()
        if index is None:
            return
        self.components.insert(index + 1, dict(self.components[index]))
        self._refresh_components()

    def on_remove_component_clicked(self, widget):
        index = self._selected_component_index()
        if index is None:
            return
        del self.components[index]
        self._refresh_components()

    def on_move_component_up_clicked(self, widget):
        index = self._selected_component_index()
        if index is None or index == 0:
            return
        self.components[index - 1], self.components[index] = self.components[index], self.components[index - 1]
        self._refresh_components()
        self.treeview_components.get_selection().select_path(index - 1)

    def on_move_component_down_clicked(self, widget):
        index = self._selected_component_index()
        if index is None or index >= len(self.components) - 1:
            return
        self.components[index + 1], self.components[index] = self.components[index], self.components[index + 1]
        self._refresh_components()
        self.treeview_components.get_selection().select_path(index + 1)

    def on_clear_components_clicked(self, widget):
        self.components = []
        self._refresh_components()

    #-------------------------------------------------------------------------
    #  S I G N A L S :  C H E C K B O X E S
    #-------------------------------------------------------------------------
    def on_checkbox_use_seed_toggled(self, widget):
        self.spinbtn_seed.set_sensitive(widget.get_active())

    #-------------------------------------------------------------------------
    #  S I G N A L S :  B U T T O N S
    #-------------------------------------------------------------------------
    def on_button_cancel_clicked(self, widget):
        self.close_window()

    #-------------------------------------------------------------------------
    #  R U N
    #-------------------------------------------------------------------------
    def on_button_run_clicked(self, widget):
        if not self.components:
            self.main.simple_dialog.info(msg='Please add at least one component to pack.')
            return

        packmol_command = self.entry_packmol_command.get_text().strip()
        if not packmol_command or not (os.path.isfile(packmol_command) and os.access(packmol_command, os.X_OK)):
            self.main.simple_dialog.info(msg='Please provide a valid path to the packmol executable.')
            return
        self.main.vm_session.vm_config.gl_parameters['packmol_command'] = packmol_command

        work_folder = self.work_folder_chooser.get_folder()
        if not work_folder:
            self.main.simple_dialog.info(msg='Please choose a working folder.')
            return
        os.makedirs(work_folder, exist_ok=True)

        system_name = self.entry_system_name.get_text().strip() or 'packmol_system'
        input_path  = os.path.join(work_folder, system_name + '.inp')
        output_path = os.path.join(work_folder, system_name + '.pdb')
        log_path    = os.path.join(work_folder, system_name + '.log')

        tolerance = self.spinbtn_tolerance.get_value()
        seed = int(self.spinbtn_seed.get_value()) if self.checkbox_use_seed.get_active() else None

        packmol_runner.build_packmol_input(
            input_path, output_path, self.components, tolerance=tolerance, seed=seed)

        self.current_output_path = output_path
        self.current_log_path = log_path
        self._log_read_offset = 0
        buf = self.textview_log.get_buffer()
        buf.set_text('')

        self.label_info.set_text('Running Packmol...')
        self._set_controls_sensitive(False)

        self.current_process = packmol_runner.run_packmol(input_path, work_folder, packmol_command, log_path)
        if self._poll_timeout_id is None:
            self._poll_timeout_id = GLib.timeout_add(500, self._poll_job)

    def _set_controls_sensitive(self, sensitive):
        self.button_run.set_sensitive(sensitive)

    def _append_log(self, text):
        buf = self.textview_log.get_buffer()
        buf.insert(buf.get_end_iter(), text)
        self.textview_log.scroll_to_iter(buf.get_end_iter(), 0.0, False, 0.0, 0.0)

    def _tail_log(self):
        if not self.current_log_path or not os.path.isfile(self.current_log_path):
            return
        with open(self.current_log_path, 'r', errors='replace') as f:
            f.seek(self._log_read_offset)
            new_text = f.read()
            self._log_read_offset = f.tell()
        if new_text:
            self._append_log(new_text)

    def _poll_job(self):
        """ GLib.timeout_add callback (every 500ms) -- same non-blocking
            polling pattern as prepare_namd_run.py's _poll_job()/
            namd_runner.run_namd(). Returns True to keep polling, False
            to stop.
        """
        self._tail_log()

        if self.current_process is None:
            return False

        returncode = self.current_process.poll()
        if returncode is None:
            return True   # still running

        self.current_process = None
        self._poll_timeout_id = None
        self._set_controls_sensitive(True)

        ok = ((returncode == 0)
              and packmol_runner.check_packmol_log_success(self.current_log_path)
              and os.path.isfile(self.current_output_path))
        if not ok:
            self.label_info.set_text('Packmol failed -- see the log above.')
            return False

        self.label_info.set_text('Done! Packed structure written to "{}".'.format(self.current_output_path))
        return False
