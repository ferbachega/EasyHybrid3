#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Prepare Ligand (Antechamber) window
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
#      "Prepare Ligand (Antechamber)" window -- runs AmberTools'
#      antechamber (atom typing + partial charges) and parmchk2
#      (missing bonded parameters) on a ligand/small-molecule file,
#      producing a .mol2 + .frcmod pair that can be handed straight to
#      the "Prepare AMBER System (tLeap)" window's Additional Files
#      list (see button_add_to_tleap below).
#

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

import os
import traceback

from gui.widgets.custom_widgets import FolderChooserButton
from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from util import antechamber_runner


class PrepareLigandAntechamberWindow:
    """ "Prepare Ligand (Antechamber)" window. """

    def __init__(self, main=None):
        """ Class initialiser """
        self.main       = main
        self.vm_session = main.vm_session
        self.p_session  = main.p_session
        self.home       = main.home
        self.Visible    = False

    #-------------------------------------------------------------------------
    #  W I N D O W   L I F E C Y C L E
    #-------------------------------------------------------------------------
    def open_window(self):
        """ Function doc """
        if self.Visible:
            self.window.present()
            return

        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/setup/prepare_ligand_antechamber_window.glade'))
        self.builder.connect_signals(self)

        self.window = self.builder.get_object('window')
        self.window.set_title('Prepare Ligand (Antechamber)')

        # -------------------- widget shortcuts --------------------
        self.label_info               = self.builder.get_object('label_info')
        self.entry_residue_name       = self.builder.get_object('entry_residue_name')
        self.spinbtn_net_charge       = self.builder.get_object('spinbtn_net_charge')
        self.spinbtn_multiplicity     = self.builder.get_object('spinbtn_multiplicity')
        self.combo_charge_method      = self.builder.get_object('combo_charge_method')
        self.combo_atom_type          = self.builder.get_object('combo_atom_type')
        self.entry_antechamber_command = self.builder.get_object('entry_antechamber_command')
        self.entry_parmchk2_command    = self.builder.get_object('entry_parmchk2_command')
        self.button_add_to_tleap      = self.builder.get_object('button_add_to_tleap')
        self.radio_input_file         = self.builder.get_object('radio_input_file')
        self.radio_input_vobject      = self.builder.get_object('radio_input_vobject')
        self.grid_input_file          = self.builder.get_object('grid_input_file')
        self.grid_input_vobject       = self.builder.get_object('grid_input_vobject')

        # -------------------- input file chooser --------------------
        self.box_input_file = self.builder.get_object('box_input_file')
        self.input_file_chooser = FolderChooserButton(self.main, 'file', None)
        self.box_input_file.pack_start(self.input_file_chooser.btn, False, False, 0)
        self.input_file_chooser.folder = None
        self.input_file_chooser.label.set_text('No file selected')

        # -------------------- input from a loaded system/object --------------------
        self.box_system = self.builder.get_object('box_system')
        self.combobox_systems = SystemComboBox(self.main)
        self.combobox_systems.connect("changed", self.on_combobox_systemsbox_changed)
        self.box_system.pack_start(self.combobox_systems, False, False, 0)

        self.box_coordinates = self.builder.get_object('box_coordinates')
        self.coordinates_combobox = CoordinatesComboBox()
        self.box_coordinates.pack_start(self.coordinates_combobox, False, False, 0)

        if self.p_session.psystem:
            self.combobox_systems.set_active_system(e_id=self.p_session.active_id)

        # -------------------- working-folder chooser --------------------
        self.box_work_folder = self.builder.get_object('box_work_folder')
        self.work_folder_chooser = FolderChooserButton(self.main, 'folder', self.home)
        self.box_work_folder.pack_start(self.work_folder_chooser.btn, False, False, 0)
        scratch = os.environ.get('PDYNAMO3_SCRATCH')
        if scratch and os.path.isdir(scratch):
            self.work_folder_chooser.set_folder(folder=scratch)

        # -------------------- comboboxes --------------------
        self._populate_combo(self.combo_charge_method, antechamber_runner.CHARGE_METHODS)
        self._populate_combo(self.combo_atom_type, antechamber_runner.ATOM_TYPES, preferred='gaff2')

        # -------------------- executables --------------------
        gl_parameters = self.main.vm_session.vm_config.gl_parameters
        antechamber_command = gl_parameters.get('antechamber_command')
        if not antechamber_command or not os.path.isfile(antechamber_command):
            antechamber_command = antechamber_runner.find_antechamber_executable()
        self.entry_antechamber_command.set_text(antechamber_command or '')

        parmchk2_command = gl_parameters.get('parmchk2_command')
        if not parmchk2_command or not os.path.isfile(parmchk2_command):
            parmchk2_command = antechamber_runner.find_parmchk2_executable()
        self.entry_parmchk2_command.set_text(parmchk2_command or '')

        if not antechamber_command or not parmchk2_command:
            self.label_info.set_text(
                'antechamber/parmchk2 were not found automatically (checked PATH and $AMBERHOME/bin) -- '
                'please enter the full paths below.')

        # -------------------- last successful run (for "Add to tLeap") --------------------
        self.last_mol2 = None
        self.last_frcmod = None
        self.button_add_to_tleap.hide()

        self.window.show_all()
        self.button_add_to_tleap.hide()  # show_all() would re-show it -- keep hidden until a run succeeds
        self.window.connect('destroy', self.close_window)
        self.Visible = True

    def close_window(self, button=None, data=None):
        """ Function doc """
        if not self.Visible:
            return
        self.window.destroy()
        self.Visible = False

    #-------------------------------------------------------------------------
    #  H E L P E R S
    #-------------------------------------------------------------------------
    def _populate_combo(self, combo, items, preferred=None):
        """ Fills a GtkComboBoxText with `items`, selecting `preferred`
            if present, otherwise the first item. """
        combo.remove_all()
        active_index = 0
        for index, item in enumerate(items):
            combo.append_text(item)
            if preferred is not None and item == preferred:
                active_index = index
        if items:
            combo.set_active(active_index)

    #-------------------------------------------------------------------------
    #  S I G N A L S
    #-------------------------------------------------------------------------
    def on_button_cancel_clicked(self, widget):
        """ Function doc """
        self.close_window()

    def on_radio_input_source_toggled(self, widget):
        """ Function doc """
        use_vobject = self.radio_input_vobject.get_active()
        self.grid_input_file.set_sensitive(not use_vobject)
        self.grid_input_vobject.set_sensitive(use_vobject)

    def on_combobox_systemsbox_changed(self, widget):
        """ Function doc """
        system_id = self.combobox_systems.get_system_id()
        if system_id is not None:
            self.coordinates_combobox.set_model(self.main.vobject_liststore_dict[system_id])
            size = len(list(self.main.vobject_liststore_dict[system_id]))
            self.coordinates_combobox.set_active(size - 1)

    def on_button_add_to_tleap_clicked(self, widget):
        """ Adds the last successfully generated .mol2/.frcmod to the
            "Prepare AMBER System (tLeap)" window's own Additional Files
            list, opening that window first if it isn't already.
        """
        if not (self.last_mol2 and self.last_frcmod):
            return

        tleap_window = self.main.prepare_amber_system_window
        tleap_window.open_window()
        for path in (self.last_mol2, self.last_frcmod):
            ext = os.path.splitext(path)[1].lower().lstrip('.')
            tleap_window.extra_files.append(path)
            tleap_window.liststore_extra_files.append([os.path.basename(path), ext])

        self.label_info.set_text('Added to "Prepare AMBER System" -- switch to that window to continue.')

    #-------------------------------------------------------------------------
    #  R U N
    #-------------------------------------------------------------------------
    def on_button_run_clicked(self, widget):
        """ Function doc """
        antechamber_command = self.entry_antechamber_command.get_text().strip()
        if not antechamber_command or not (os.path.isfile(antechamber_command) and os.access(antechamber_command, os.X_OK)):
            self.main.simple_dialog.info(msg='Please provide a valid path to the antechamber executable.')
            return

        parmchk2_command = self.entry_parmchk2_command.get_text().strip()
        if not parmchk2_command or not (os.path.isfile(parmchk2_command) and os.access(parmchk2_command, os.X_OK)):
            self.main.simple_dialog.info(msg='Please provide a valid path to the parmchk2 executable.')
            return

        # Persisted the same way as tleap_command -- see vm_config.gl_parameters.
        gl_parameters = self.main.vm_session.vm_config.gl_parameters
        gl_parameters['antechamber_command'] = antechamber_command
        gl_parameters['parmchk2_command'] = parmchk2_command

        work_folder = self.work_folder_chooser.get_folder()
        if not work_folder:
            self.main.simple_dialog.info(msg='Please choose a working folder.')
            return

        residue_name = self.entry_residue_name.get_text().strip() or 'LIG'

        # -------------------- resolve the input structure --------------------
        input_format = None
        if self.radio_input_vobject.get_active():
            # . From a loaded system/object -- export it to a PDB (same
            #   helper already used by the tLeap window) and feed that
            #   to antechamber. Meant for a small, already-isolated
            #   ligand object, not a whole protein/complex.
            vobject_id = self.coordinates_combobox.get_vobject_id()
            if vobject_id is None:
                self.main.simple_dialog.info(msg='Please select a system/object.')
                return
            vobject = self.main.vm_session.vm_objects_dic.get(vobject_id)
            if vobject is None:
                self.main.simple_dialog.info(msg='Please select a system/object.')
                return

            input_path = os.path.join(work_folder, residue_name + '_input.pdb')
            try:
                antechamber_runner.write_ligand_pdb(vobject=vobject, residue_name=residue_name, output_path=input_path)
            except Exception as error:
                traceback.print_exc()
                self.main.simple_dialog.info(msg='Could not export the structure to PDB:\n{}'.format(error))
                return
            input_format = 'pdb'
        else:
            input_path = self.input_file_chooser.get_folder()
            if not input_path or not os.path.isfile(input_path):
                self.main.simple_dialog.info(msg='Please select a ligand file.')
                return
        atom_type = self.combo_atom_type.get_active_text()
        charge_method = self.combo_charge_method.get_active_text()

        self.button_add_to_tleap.hide()
        self.last_mol2 = None
        self.last_frcmod = None

        self.label_info.set_text('Running antechamber (this can take a while for AM1-BCC charges)...')
        while Gtk.events_pending():
            Gtk.main_iteration()

        output_mol2 = os.path.join(work_folder, residue_name + '.mol2')
        result = antechamber_runner.run_antechamber(
            input_path=input_path,
            output_mol2_path=output_mol2,
            charge_method=charge_method,
            net_charge=int(self.spinbtn_net_charge.get_value()),
            multiplicity=int(self.spinbtn_multiplicity.get_value()),
            residue_name=residue_name,
            atom_type=atom_type,
            antechamber_command=antechamber_command,
            workdir=work_folder,
            input_format=input_format,
        )

        if not result['ok']:
            self.label_info.set_text('antechamber failed -- see the error dialog for details.')
            self.main.simple_dialog.error_details(
                parent=self.window,
                msg='antechamber did not produce a valid parametrized .mol2 for "{}".'.format(
                    os.path.basename(input_path)),
                details=result['stdout'] + result['stderr'],
                title='Antechamber Error',
            )
            return

        self.label_info.set_text('Running parmchk2...')
        while Gtk.events_pending():
            Gtk.main_iteration()

        output_frcmod = os.path.join(work_folder, residue_name + '.frcmod')
        result2 = antechamber_runner.run_parmchk2(
            mol2_path=result['mol2'],
            output_frcmod_path=output_frcmod,
            atom_type=atom_type,
            parmchk2_command=parmchk2_command,
            workdir=work_folder,
        )

        if not result2['ok']:
            self.label_info.set_text('parmchk2 failed -- see the error dialog for details.')
            self.main.simple_dialog.error_details(
                parent=self.window,
                msg='antechamber succeeded ("{}"), but parmchk2 could not generate the .frcmod.'.format(
                    os.path.basename(result['mol2'])),
                details=result2['stdout'] + result2['stderr'],
                title='parmchk2 Error',
            )
            return

        self.last_mol2 = result['mol2']
        self.last_frcmod = result2['frcmod']
        self.button_add_to_tleap.show()
        self.label_info.set_text('Done! "{}" and "{}" were generated.'.format(
            os.path.basename(self.last_mol2), os.path.basename(self.last_frcmod)))
