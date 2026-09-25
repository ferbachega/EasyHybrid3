#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: "Run AutoDock Vina" window
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
#      GUI port of ~/programs/adtVinaFR: prepares a receptor and a
#      ligand (each either a folder of PDB files, matching the original
#      standalone script, or a loaded EasyHybrid System/Object -- one
#      frame or a whole range exported as an ensemble of PDB snapshots)
#      into PDBQT via OpenBabel, then cross-docks every receptor PDBQT
#      against every ligand PDBQT with AutoDock Vina. See
#      util/vina_runner.py for the actual OpenBabel/Vina wrapping
#      (GTK-free, independently testable) -- this file only wires that
#      up to the window's widgets and runs it off the GTK main thread.
#
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, GLib

import os
import shutil
import subprocess
import threading
import traceback

import numpy as np

from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.widgets.custom_widgets import FolderChooserButton
from pdynamo.pDynamo2EasyHybrid.helpers import export_special_PDB
from util import vina_runner
from util import pdbqt_parser
from util import vismol_pdbqt_builder
from util.vismol_pdbqt_builder import DockingResult


class PrepareVinaDockingWindow:

    def __init__(self, main=None):
        """ Class initialiser """
        self.main = main
        self.home = main.home
        self.p_session = main.p_session
        self.vm_session = main.vm_session
        self.Visible = False

        self._receptor_pdbqt_list = []
        self._ligand_pdbqt_list = []
        self._results = []
        self._cancel_event = None
        self._receptor_vobject = None
        self._receptor_id_to_frame = None
        self._docking_results = []
        self._docking_box_vobject = None

    # -------------------------------------------------------------------
    #  Window setup
    # -------------------------------------------------------------------

    def open_window(self):
        """ Function doc """
        if self.Visible:
            self.window.present()
            return

        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(
            self.home, 'src/gui/windows/setup/prepare_vina_docking_window.glade'))
        self.builder.connect_signals(self)

        self.window = self.builder.get_object('window')

        # . Output & execution
        self.box_output_folder = self.builder.get_object('box_output_folder')
        self.output_folder_chooser = FolderChooserButton(self.main, 'folder', self.home)
        self.box_output_folder.pack_start(self.output_folder_chooser.btn, False, False, 0)
        self.entry_system_name = self.builder.get_object('entry_system_name')
        self.entry_obabel_command = self.builder.get_object('entry_obabel_command')
        self.entry_vina_command = self.builder.get_object('entry_vina_command')

        # . Receptor tab
        self.combo_receptor_source = self.builder.get_object('combo_receptor_source')
        self.box_receptor_single_file = self.builder.get_object('box_receptor_single_file')
        self.box_receptor_folder = self.builder.get_object('box_receptor_folder')
        self.box_receptor_object = self.builder.get_object('box_receptor_object')
        self.box_receptor_single_pdb = self.builder.get_object('box_receptor_single_pdb')
        self.receptor_single_file_chooser = FolderChooserButton(self.main, 'file', self.home)
        self.box_receptor_single_pdb.pack_start(self.receptor_single_file_chooser.btn, False, False, 0)
        self.box_receptor_pdb_folder = self.builder.get_object('box_receptor_pdb_folder')
        self.receptor_folder_chooser = FolderChooserButton(self.main, 'folder', self.home)
        self.box_receptor_pdb_folder.pack_start(self.receptor_folder_chooser.btn, False, False, 0)
        self.box_receptor_system = self.builder.get_object('box_receptor_system')
        self.combobox_receptor_systems = SystemComboBox(self.main)
        self.combobox_receptor_systems.connect("changed", self._on_receptor_system_changed)
        self.box_receptor_system.pack_start(self.combobox_receptor_systems, True, True, 0)
        self.box_receptor_coordinates = self.builder.get_object('box_receptor_coordinates')
        self.coordinates_combobox_receptor = CoordinatesComboBox()
        self.box_receptor_coordinates.pack_start(self.coordinates_combobox_receptor, True, True, 0)
        self.chk_receptor_all_frames = self.builder.get_object('chk_receptor_all_frames')
        self.spinbtn_receptor_frame_init = self.builder.get_object('spinbtn_receptor_frame_init')
        self.spinbtn_receptor_frame_last = self.builder.get_object('spinbtn_receptor_frame_last')
        self.spinbtn_receptor_step = self.builder.get_object('spinbtn_receptor_step')
        self.label_receptor_status = self.builder.get_object('label_receptor_status')
        self.button_load_receptor_vismol = self.builder.get_object('button_load_receptor_vismol')
        self.spinbtn_receptor_goto_frame = self.builder.get_object('spinbtn_receptor_goto_frame')

        # . Ligand tab
        self.combo_ligand_source = self.builder.get_object('combo_ligand_source')
        self.box_ligand_folder = self.builder.get_object('box_ligand_folder')
        self.box_ligand_object = self.builder.get_object('box_ligand_object')
        self.box_ligand_pdb_folder = self.builder.get_object('box_ligand_pdb_folder')
        self.ligand_folder_chooser = FolderChooserButton(self.main, 'folder', self.home)
        self.box_ligand_pdb_folder.pack_start(self.ligand_folder_chooser.btn, False, False, 0)
        self.box_ligand_system = self.builder.get_object('box_ligand_system')
        self.combobox_ligand_systems = SystemComboBox(self.main)
        self.combobox_ligand_systems.connect("changed", self._on_ligand_system_changed)
        self.box_ligand_system.pack_start(self.combobox_ligand_systems, True, True, 0)
        self.box_ligand_coordinates = self.builder.get_object('box_ligand_coordinates')
        self.coordinates_combobox_ligand = CoordinatesComboBox()
        self.box_ligand_coordinates.pack_start(self.coordinates_combobox_ligand, True, True, 0)
        self.chk_ligand_all_frames = self.builder.get_object('chk_ligand_all_frames')
        self.spinbtn_ligand_frame_init = self.builder.get_object('spinbtn_ligand_frame_init')
        self.spinbtn_ligand_frame_last = self.builder.get_object('spinbtn_ligand_frame_last')
        self.spinbtn_ligand_step = self.builder.get_object('spinbtn_ligand_step')
        self.spinbtn_ligand_ph = self.builder.get_object('spinbtn_ligand_ph')
        self.label_ligand_status = self.builder.get_object('label_ligand_status')
        self.liststore_ligands = self.builder.get_object('liststore_ligand_list')
        self.treeview_ligands = self.builder.get_object('treeview_ligand_list')

        # . Docking tab
        self.spinbtn_center_x = self.builder.get_object('spinbtn_center_x')
        self.spinbtn_center_y = self.builder.get_object('spinbtn_center_y')
        self.spinbtn_center_z = self.builder.get_object('spinbtn_center_z')
        self.spinbtn_size_x = self.builder.get_object('spinbtn_size_x')
        self.spinbtn_size_y = self.builder.get_object('spinbtn_size_y')
        self.spinbtn_size_z = self.builder.get_object('spinbtn_size_z')
        self.spinbtn_box_padding = self.builder.get_object('spinbtn_box_padding')
        self.spinbtn_exhaustiveness = self.builder.get_object('spinbtn_exhaustiveness')
        self.spinbtn_num_modes = self.builder.get_object('spinbtn_num_modes')
        self.spinbtn_energy_range = self.builder.get_object('spinbtn_energy_range')
        self.chk_cpu_limit = self.builder.get_object('chk_cpu_limit')
        self.spinbtn_cpu = self.builder.get_object('spinbtn_cpu')
        self.chk_seed = self.builder.get_object('chk_seed')
        self.spinbtn_seed = self.builder.get_object('spinbtn_seed')
        self.label_docking_status = self.builder.get_object('label_docking_status')
        self.textview_log = self.builder.get_object('textview_log')
        self.button_run = self.builder.get_object('button_run')
        self.button_stop = self.builder.get_object('button_stop')

        # . Results tab
        self.liststore_results = self.builder.get_object('liststore_results')
        self.treeview_results = self.builder.get_object('treeview_results')

        self._sync_executable_paths()

        if self.p_session.psystem:
            self.combobox_receptor_systems.set_active_system(e_id=self.p_session.active_id)
            self.combobox_ligand_systems.set_active_system(e_id=self.p_session.active_id)

        # . Initial visibility sync -- show_all() below would otherwise
        #   force every conditionally-shown widget visible regardless of
        #   the combo/checkbox state it should be following (same
        #   pattern every other window in this codebase uses).
        self.on_combo_receptor_source_changed(self.combo_receptor_source)
        self.on_combo_ligand_source_changed(self.combo_ligand_source)
        self.on_chk_cpu_limit_toggled(self.chk_cpu_limit)
        self.on_chk_seed_toggled(self.chk_seed)

        self.window.show_all()
        self.on_combo_receptor_source_changed(self.combo_receptor_source)
        self.on_combo_ligand_source_changed(self.combo_ligand_source)
        self.on_chk_cpu_limit_toggled(self.chk_cpu_limit)
        self.on_chk_seed_toggled(self.chk_seed)
        self._refresh_ligand_liststore()
        self.Visible = True

    def close_window(self, button=None, data=None):
        """ Function doc """
        if not self.Visible:
            return
        self._unpin_docking_results()
        self.window.destroy()
        self.Visible = False

    def _sync_executable_paths(self):
        """ Same persisted-preference pattern namd_runner.py's
            'namd_command' uses: a previously-typed/validated path
            persists in gl_parameters across sessions; otherwise fall
            back to a fresh PATH search.
        """
        gl_parameters = self.vm_session.vm_config.gl_parameters
        obabel_command = gl_parameters.get('obabel_command')
        if not obabel_command or not os.path.isfile(obabel_command):
            obabel_command = vina_runner.find_obabel_executable()
        self.entry_obabel_command.set_text(obabel_command or '')

        vina_command = gl_parameters.get('vina_command')
        if not vina_command or not os.path.isfile(vina_command):
            vina_command = vina_runner.find_vina_executable()
        self.entry_vina_command.set_text(vina_command or '')

    # -------------------------------------------------------------------
    #  Source-type / sensitivity toggles
    # -------------------------------------------------------------------

    def on_combo_receptor_source_changed(self, widget):
        source = widget.get_active_id()
        self.box_receptor_single_file.set_visible(source == 'file')
        self.box_receptor_folder.set_visible(source == 'folder')
        self.box_receptor_object.set_visible(source == 'vobject')

    def on_combo_ligand_source_changed(self, widget):
        source = widget.get_active_id()
        self.box_ligand_folder.set_visible(source == 'folder')
        self.box_ligand_object.set_visible(source == 'vobject')

    def on_chk_receptor_all_frames_toggled(self, widget):
        active = widget.get_active()
        self.spinbtn_receptor_frame_init.set_sensitive(active)
        self.spinbtn_receptor_frame_last.set_sensitive(active)
        self.spinbtn_receptor_step.set_sensitive(active)

    def on_chk_ligand_all_frames_toggled(self, widget):
        active = widget.get_active()
        self.spinbtn_ligand_frame_init.set_sensitive(active)
        self.spinbtn_ligand_frame_last.set_sensitive(active)
        self.spinbtn_ligand_step.set_sensitive(active)

    def on_chk_cpu_limit_toggled(self, widget):
        self.spinbtn_cpu.set_sensitive(widget.get_active())

    def on_chk_seed_toggled(self, widget):
        self.spinbtn_seed.set_sensitive(widget.get_active())

    def _on_receptor_system_changed(self, widget):
        system_id = self.combobox_receptor_systems.get_system_id()
        if system_id is None:
            return
        self.coordinates_combobox_receptor.set_model(self.main.vobject_liststore_dict[system_id])
        size = len(list(self.main.vobject_liststore_dict[system_id]))
        self.coordinates_combobox_receptor.set_active(size - 1)

    def _on_ligand_system_changed(self, widget):
        system_id = self.combobox_ligand_systems.get_system_id()
        if system_id is None:
            return
        self.coordinates_combobox_ligand.set_model(self.main.vobject_liststore_dict[system_id])
        size = len(list(self.main.vobject_liststore_dict[system_id]))
        self.coordinates_combobox_ligand.set_active(size - 1)

    # -------------------------------------------------------------------
    #  Executable browse
    # -------------------------------------------------------------------

    def _browse_for_executable(self, entry):
        dialog = Gtk.FileChooserDialog(
            title="Select executable", transient_for=self.window,
            action=Gtk.FileChooserAction.OPEN)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                            Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            entry.set_text(dialog.get_filename())
        dialog.destroy()

    def on_button_browse_obabel_clicked(self, widget):
        self._browse_for_executable(self.entry_obabel_command)

    def on_button_browse_vina_clicked(self, widget):
        self._browse_for_executable(self.entry_vina_command)

    # -------------------------------------------------------------------
    #  Frame-range resolution / PDB materialisation (shared by receptor/ligand)
    # -------------------------------------------------------------------

    def _resolve_frame_indices(self, vobject, chk_all_frames, spin_init, spin_last, spin_step):
        if not chk_all_frames.get_active():
            return [self.vm_session.frame]
        f_init = int(spin_init.get_value())
        f_last = int(spin_last.get_value())
        if f_last == -1:
            f_last = vobject.frames.shape[0]
        step = max(1, int(spin_step.get_value()))
        return list(range(f_init, f_last, step))

    def _export_vobject_frames_to_pdb_folder(self, vobject, frame_indices, folder, basename):
        """ Materialises `frame_indices` of `vobject` as one PDB per
            frame in `folder` -- once on disk, this is indistinguishable
            from a user-supplied "folder of PDB files", so the SAME
            vina_runner.prepare_*_ensemble() functions handle both input
            sources unmodified.
        """
        os.makedirs(folder, exist_ok=True)
        for frame in frame_indices:
            path = os.path.join(folder, '{}_{:04d}.pdb'.format(basename, frame))
            export_special_PDB(vobject, frame, path)
        return folder

    def _get_work_folder(self):
        folder = self.output_folder_chooser.get_folder()
        if not folder:
            self.main.simple_dialog.info(msg='Please choose a working folder first.')
            return None
        return folder

    def _get_system_name(self):
        return self.entry_system_name.get_text().strip() or 'vina_run'

    # -------------------------------------------------------------------
    #  Prepare receptor / ligand
    # -------------------------------------------------------------------

    def _resolve_obabel_command(self):
        obabel_command = self.entry_obabel_command.get_text().strip()
        if not obabel_command or not (os.path.isfile(obabel_command) and os.access(obabel_command, os.X_OK)):
            self.main.simple_dialog.info(
                msg='Could not find the OpenBabel ("obabel") executable. Please set it manually.')
            return None
        self.vm_session.vm_config.gl_parameters['obabel_command'] = obabel_command
        return obabel_command

    def _resolve_vina_command(self):
        vina_command = self.entry_vina_command.get_text().strip()
        if not vina_command or not (os.path.isfile(vina_command) and os.access(vina_command, os.X_OK)):
            self.main.simple_dialog.info(
                msg='Could not find the AutoDock Vina ("vina") executable. Please set it manually.')
            return None
        self.vm_session.vm_config.gl_parameters['vina_command'] = vina_command
        return vina_command

    def on_button_prepare_receptor_clicked(self, widget):
        obabel_command = self._resolve_obabel_command()
        work_folder = self._get_work_folder()
        if obabel_command is None or work_folder is None:
            return

        system_name = self._get_system_name()
        source = self.combo_receptor_source.get_active_id()
        if source == 'file':
            pdb_file = self.receptor_single_file_chooser.get_folder()
            if not pdb_file:
                self.main.simple_dialog.info(msg='Please choose a receptor PDB file.')
                return
            pdb_folder = os.path.join(work_folder, system_name + '_receptor_PDB')
            os.makedirs(pdb_folder, exist_ok=True)
            shutil.copy(pdb_file, os.path.join(pdb_folder, os.path.basename(pdb_file)))
        elif source == 'folder':
            pdb_folder = self.receptor_folder_chooser.get_folder()
            if not pdb_folder:
                self.main.simple_dialog.info(msg='Please choose a folder with receptor PDB file(s).')
                return
        else:
            vobject_id = self.coordinates_combobox_receptor.get_vobject_id()
            vobject = self.vm_session.vm_objects_dic.get(vobject_id) if vobject_id is not None else None
            if vobject is None:
                self.main.simple_dialog.info(msg='Please select a System and an Object for the receptor.')
                return
            frame_indices = self._resolve_frame_indices(
                vobject, self.chk_receptor_all_frames,
                self.spinbtn_receptor_frame_init, self.spinbtn_receptor_frame_last, self.spinbtn_receptor_step)
            pdb_folder = self._export_vobject_frames_to_pdb_folder(
                vobject, frame_indices, os.path.join(work_folder, system_name + '_receptor_PDB'), 'receptor')

        pdbqt_folder = os.path.join(work_folder, system_name + '_receptor_PDBQT')
        try:
            pdbqt_list, failures = vina_runner.prepare_receptor_ensemble(
                pdb_folder, pdbqt_folder, obabel_bin=obabel_command)
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg='Could not prepare the receptor(s):\n{}'.format(error))
            return

        self._receptor_pdbqt_list = pdbqt_list
        if failures:
            text = '{} receptor(s) ready, {} failed:\n'.format(len(pdbqt_list), len(failures))
            text += '\n'.join('{}: {}'.format(os.path.basename(f), msg.splitlines()[-1] if msg else '')
                               for f, msg in failures)
            self.label_receptor_status.set_text(text)
        else:
            self.label_receptor_status.set_text('{} receptor(s) ready (PDBQT).'.format(len(pdbqt_list)))

    def on_button_load_receptor_vismol_clicked(self, widget):
        """ Loads self._receptor_pdbqt_list directly into VisMol -- a
            single receptor becomes a 1-frame object; more than one
            becomes ONE object with a frame per receptor. No temporary
            PDB file, no pDynamo ImportSystem() call reading a file --
            see vismol_pdbqt_builder.py's own module docstring for how
            this builds a minimal, in-memory pDynamo System instead.
        """
        if not self._receptor_pdbqt_list:
            self.main.simple_dialog.info(msg='Prepare the receptor(s) first (Receptor tab).')
            return
        try:
            if len(self._receptor_pdbqt_list) > 1:
                vm_object = vismol_pdbqt_builder.load_receptor_ensemble_pdbqt(
                    self.p_session, pdbqt_parser, self._receptor_pdbqt_list)
            else:
                vm_object = vismol_pdbqt_builder.load_receptor_pdbqt(
                    self.p_session, pdbqt_parser, self._receptor_pdbqt_list[0])
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg='Could not load the receptor into VisMol:\n{}'.format(error))
            return
        self._receptor_vobject = vm_object
        self._receptor_id_to_frame = None
        self.label_receptor_status.set_text(
            self.label_receptor_status.get_text() + '\nLoaded into VisMol as "{}" ({} frame(s)).'.format(
                vm_object.name, vm_object.frames.shape[0]))

    def on_button_goto_receptor_frame_clicked(self, widget):
        """ Lets the user pin the receptor to one specific ensemble
            frame directly (independent of any docked ligand pose's own
            pin) -- e.g. to inspect receptor_id N on its own before/
            without importing any docking result. Selecting a pose in
            the Results tab afterward may pin the receptor again to
            that pose's own associated frame, same as any other pin.
        """
        if self._receptor_vobject is None:
            self.main.simple_dialog.info(msg='Load the receptor into VisMol first.')
            return
        frame = int(self.spinbtn_receptor_goto_frame.get_value())
        self._receptor_vobject.pin_frame(frame)
        self.vm_session.vm_glcore.updated_coords = True
        self.vm_session.vm_widget.queue_draw()

    def on_button_prepare_ligand_clicked(self, widget):
        obabel_command = self._resolve_obabel_command()
        work_folder = self._get_work_folder()
        if obabel_command is None or work_folder is None:
            return

        system_name = self._get_system_name()
        source = self.combo_ligand_source.get_active_id()
        if source == 'folder':
            pdb_folder = self.ligand_folder_chooser.get_folder()
            if not pdb_folder:
                self.main.simple_dialog.info(msg='Please choose a folder with ligand file(s) (.pdb/.mol/.mol2/.sdf/.xyz).')
                return
        else:
            vobject_id = self.coordinates_combobox_ligand.get_vobject_id()
            vobject = self.vm_session.vm_objects_dic.get(vobject_id) if vobject_id is not None else None
            if vobject is None:
                self.main.simple_dialog.info(msg='Please select a System and an Object for the ligand.')
                return
            frame_indices = self._resolve_frame_indices(
                vobject, self.chk_ligand_all_frames,
                self.spinbtn_ligand_frame_init, self.spinbtn_ligand_frame_last, self.spinbtn_ligand_step)
            pdb_folder = self._export_vobject_frames_to_pdb_folder(
                vobject, frame_indices, os.path.join(work_folder, system_name + '_ligand_PDB'), 'ligand')

        pdbqt_folder = os.path.join(work_folder, system_name + '_ligand_PDBQT')
        ph = self.spinbtn_ligand_ph.get_value()
        try:
            pdbqt_list, failures = vina_runner.prepare_ligand_ensemble(
                pdb_folder, pdbqt_folder, ph=ph, obabel_bin=obabel_command)
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg='Could not prepare the ligand(s):\n{}'.format(error))
            return

        added = self._add_ligands_to_list(pdbqt_list)
        if failures:
            text = '{} ligand(s) added ({} total in list), {} failed:\n'.format(
                added, len(self._ligand_pdbqt_list), len(failures))
            text += '\n'.join('{}: {}'.format(os.path.basename(f), msg.splitlines()[-1] if msg else '')
                               for f, msg in failures)
            self.label_ligand_status.set_text(text)
        else:
            self.label_ligand_status.set_text(
                '{} ligand(s) added ({} total in list).'.format(added, len(self._ligand_pdbqt_list)))

    def _add_ligands_to_list(self, paths):
        """ Appends `paths` to self._ligand_pdbqt_list, skipping any
            path already present (so re-preparing the same folder, or
            re-importing the same PDBQT, doesn't duplicate an entry),
            then refreshes the treeview. Returns how many were
            actually newly added.
        """
        existing = set(self._ligand_pdbqt_list)
        new_paths = [p for p in paths if p not in existing]
        self._ligand_pdbqt_list.extend(new_paths)
        self._refresh_ligand_liststore()
        return len(new_paths)

    def _refresh_ligand_liststore(self):
        self.liststore_ligands.clear()
        for path in self._ligand_pdbqt_list:
            self.liststore_ligands.append([os.path.basename(path), path])

    def on_button_import_ligand_pdbqt_clicked(self, widget):
        """ Adds one or more already-prepared ligand .pdbqt files
            directly to self._ligand_pdbqt_list -- no OpenBabel
            conversion, they're used exactly as selected (unlike
            "Prepare Ligand(s)", which always converts from a source
            format via vina_runner.prepare_ligand_ensemble()).
        """
        dialog = Gtk.FileChooserDialog(
            title='Select ligand PDBQT file(s)', parent=self.window, action=Gtk.FileChooserAction.OPEN)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
        dialog.set_select_multiple(True)
        file_filter = Gtk.FileFilter()
        file_filter.set_name('PDBQT files')
        file_filter.add_pattern('*.pdbqt')
        dialog.add_filter(file_filter)
        response = dialog.run()
        paths = dialog.get_filenames() if response == Gtk.ResponseType.OK else []
        dialog.destroy()
        if not paths:
            return
        added = self._add_ligands_to_list(paths)
        self.label_ligand_status.set_text(
            '{} ligand(s) imported ({} total in list).'.format(added, len(self._ligand_pdbqt_list)))

    def on_button_remove_ligand_clicked(self, widget):
        model, tree_paths = self.treeview_ligands.get_selection().get_selected_rows()
        paths_to_remove = {model[tp][1] for tp in tree_paths}
        if not paths_to_remove:
            return
        self._ligand_pdbqt_list = [p for p in self._ligand_pdbqt_list if p not in paths_to_remove]
        self._refresh_ligand_liststore()
        self.label_ligand_status.set_text(
            '{} ligand(s) removed ({} left in list).'.format(len(paths_to_remove), len(self._ligand_pdbqt_list)))

    def on_button_clear_ligand_list_clicked(self, widget):
        self._ligand_pdbqt_list = []
        self._refresh_ligand_liststore()
        self.label_ligand_status.set_text('Ligand list cleared.')

    # -------------------------------------------------------------------
    #  Docking box from the current 3D-view selection
    # -------------------------------------------------------------------

    def on_button_box_from_selection_clicked(self, widget):
        selected_list, _residue_dict, vobject = self.vm_session.build_index_list_from_atom_selection(
            return_vobject=True)
        if not selected_list:
            self.main.simple_dialog.info(
                msg='No atoms are currently selected in the 3D view. Select the binding-pocket '
                    'atoms/residues first, then click this button.')
            return
        coords = np.array([vobject.atoms[i].coords() for i in selected_list], dtype=np.float64)
        padding = self.spinbtn_box_padding.get_value()
        box = vina_runner.compute_box_from_coordinates(coords, padding=padding)
        self.spinbtn_center_x.set_value(box['center_x'])
        self.spinbtn_center_y.set_value(box['center_y'])
        self.spinbtn_center_z.set_value(box['center_z'])
        self.spinbtn_size_x.set_value(box['size_x'])
        self.spinbtn_size_y.set_value(box['size_y'])
        self.spinbtn_size_z.set_value(box['size_z'])

    def on_button_center_on_centroid_clicked(self, widget):
        """ Sets ONLY Center (not Size) to the plain centroid (average
            position) of the current 3D-view selection -- see
            vina_runner.compute_centroid()'s docstring for how this
            differs from on_button_box_from_selection_clicked() above,
            which centres on the selection's bounding-box midpoint
            instead. Size is left exactly as it was.
        """
        selected_list, _residue_dict, vobject = self.vm_session.build_index_list_from_atom_selection(
            return_vobject=True)
        if not selected_list:
            self.main.simple_dialog.info(
                msg='No atoms are currently selected in the 3D view. Select the atoms/residues '
                    'to centre on first, then click this button.')
            return
        coords = np.array([vobject.atoms[i].coords() for i in selected_list], dtype=np.float64)
        centroid = vina_runner.compute_centroid(coords)
        self.spinbtn_center_x.set_value(centroid['center_x'])
        self.spinbtn_center_y.set_value(centroid['center_y'])
        self.spinbtn_center_z.set_value(centroid['center_z'])

    def on_button_show_docking_box_clicked(self, widget):
        """ Draws (or, on a later click, moves/resizes in place) the
            current Center/Size as a wireframe box in the 3D view,
            using VisMol's own unit-cell representation -- see
            vismol_pdbqt_builder.create_or_update_docking_box()'s own
            docstring for why a small dummy-atom VObject hosts it.
        """
        box = self._get_box()
        name = '{}_docking_box'.format(self._get_system_name())
        try:
            self._docking_box_vobject = vismol_pdbqt_builder.create_or_update_docking_box(
                self.p_session, self.vm_session, pdbqt_parser, self._docking_box_vobject, box, name)
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg='Could not draw the docking box:\n{}'.format(error))

    def on_button_hide_docking_box_clicked(self, widget):
        if self._docking_box_vobject is not None:
            self.vm_session.hide_cell(self._docking_box_vobject)

    # -------------------------------------------------------------------
    #  Run / Stop
    # -------------------------------------------------------------------

    def _get_box(self):
        return {
            'center_x': self.spinbtn_center_x.get_value(), 'center_y': self.spinbtn_center_y.get_value(),
            'center_z': self.spinbtn_center_z.get_value(),
            'size_x': self.spinbtn_size_x.get_value(), 'size_y': self.spinbtn_size_y.get_value(),
            'size_z': self.spinbtn_size_z.get_value(),
        }

    def _get_params(self):
        params = {
            'exhaustiveness': int(self.spinbtn_exhaustiveness.get_value()),
            'num_modes': int(self.spinbtn_num_modes.get_value()),
            'energy_range': self.spinbtn_energy_range.get_value(),
        }
        if self.chk_cpu_limit.get_active():
            params['cpu'] = int(self.spinbtn_cpu.get_value())
        if self.chk_seed.get_active():
            params['seed'] = int(self.spinbtn_seed.get_value())
        return params

    def _set_controls_sensitive(self, sensitive):
        self.button_run.set_sensitive(sensitive)
        self.button_stop.set_sensitive(not sensitive)

    def on_button_run_clicked(self, widget):
        if not self._receptor_pdbqt_list:
            self.main.simple_dialog.info(msg='Prepare the receptor(s) first (Receptor tab).')
            return
        if not self._ligand_pdbqt_list:
            self.main.simple_dialog.info(msg='Prepare the ligand(s) first (Ligand tab).')
            return
        vina_command = self._resolve_vina_command()
        if vina_command is None:
            return
        work_folder = self._get_work_folder()
        if work_folder is None:
            return

        output_folder = os.path.join(work_folder, self._get_system_name() + '_OUTPUTS')
        self._output_folder = output_folder
        box = self._get_box()
        params = self._get_params()

        buf = self.textview_log.get_buffer()
        buf.set_text('')
        self._unpin_docking_results()
        self.liststore_results.clear()
        self._results = []
        self.label_docking_status.set_text('Starting...')
        self._set_controls_sensitive(False)

        self._cancel_event = threading.Event()
        thread = threading.Thread(
            target=self._run_docking_thread,
            args=(vina_command, box, params, output_folder, self._cancel_event))
        thread.start()

    def on_button_stop_clicked(self, widget):
        if self._cancel_event is not None:
            self._cancel_event.set()
            self.label_docking_status.set_text('Stopping (finishing the current job)...')

    def _run_docking_thread(self, vina_command, box, params, output_folder, cancel_event):
        def progress(event, info):
            GLib.idle_add(self._on_progress, event, dict(info))

        try:
            rows = vina_runner.run_vina_docking(
                vina_command, self._receptor_pdbqt_list, self._ligand_pdbqt_list,
                output_folder, box, params=params,
                progress_callback=progress, cancel_event=cancel_event)
        except Exception as error:
            traceback.print_exc()
            GLib.idle_add(self._on_docking_error, str(error))
            return
        GLib.idle_add(self._on_docking_finished, rows)

    def _append_log(self, text):
        buf = self.textview_log.get_buffer()
        buf.insert(buf.get_end_iter(), text)
        self.textview_log.scroll_to_iter(buf.get_end_iter(), 0.0, False, 0.0, 0.0)

    def _on_progress(self, event, info):
        if event == 'job_start':
            self.label_docking_status.set_text('Docking {}/{}: {} vs {}...'.format(
                info['index'], info['total'], info['ligand'], info['receptor']))
            self._append_log('[{}/{}] {} vs {} ... '.format(
                info['index'], info['total'], info['ligand'], info['receptor']))
        elif event == 'job_done':
            best = min((r['affinity'] for r in info['rows']), default=None)
            self._append_log('done (best affinity {} kcal/mol)\n'.format(best))
            for row in info['rows']:
                self.liststore_results.append([
                    row['receptor'], row['ligand'], row['mode'],
                    row['affinity'], row['rmsd_lb'], row['rmsd_ub'],
                    -1, -1, -1, -1])  # vobject_id, frame_index, receptor_vobject_id,
                                      # receptor_frame_index -- set once this pose is imported
        elif event == 'job_failed':
            self._append_log('FAILED: {}\n'.format(info['error']))
        elif event == 'cancelled':
            self._append_log('Cancelled by user.\n')
        return False

    def _on_docking_finished(self, rows):
        self._results = rows
        self._set_controls_sensitive(True)
        self.label_docking_status.set_text(
            'Finished -- {} result row(s). See the Results tab. Every job\'s command/config '
            'was also saved to {} (run_all.sh reruns this whole batch outside EasyHybrid).'.format(
                len(rows), self._output_folder))
        return False

    def _on_docking_error(self, message):
        self._set_controls_sensitive(True)
        self.label_docking_status.set_text('Failed -- see log.')
        self._append_log('ERROR: {}\n'.format(message))
        return False

    # -------------------------------------------------------------------
    #  Results
    # -------------------------------------------------------------------

    def on_button_import_previous_results_clicked(self, widget):
        """ Points at a previous "..._OUTPUTS" folder (from this tool's
            own run_vina_docking(), which always writes an
            affinity_logs.txt there) and ADDS every row it can still
            find a real docked PDBQT for to the table below, alongside
            whatever is already loaded -- no docking is re-run.
        """
        dialog = Gtk.FileChooserDialog(
            title='Select a previous docking OUTPUTS folder', parent=self.window,
            action=Gtk.FileChooserAction.SELECT_FOLDER)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
        response = dialog.run()
        folder = dialog.get_filename() if response == Gtk.ResponseType.OK else None
        dialog.destroy()
        if not folder:
            return

        try:
            rows, missing_pairs = vina_runner.parse_affinity_log(folder)
        except FileNotFoundError:
            self.main.simple_dialog.info(
                msg='No affinity_logs.txt found in:\n{}\n\n'
                    'This must be an "..._OUTPUTS" folder generated by this tool.'.format(folder))
            return
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg='Could not read previous results:\n{}'.format(error))
            return

        for row in rows:
            self._results.append(row)
            self.liststore_results.append([
                row['receptor'], row['ligand'], row['mode'], row['affinity'],
                row['rmsd_lb'], row['rmsd_ub'], -1, -1, -1, -1])

        status = '{} previous result row(s) added ({} total now).'.format(len(rows), len(self._results))
        if missing_pairs:
            status += ' {} pair(s) skipped (docked PDBQT no longer found).'.format(len(missing_pairs))
        self.label_docking_status.set_text(status)

    def on_button_clear_results_clicked(self, widget):
        self._unpin_docking_results()
        self.liststore_results.clear()
        self._results = []

    def _unpin_docking_results(self):
        """ Releases every pin_frame() set by
            on_treeview_results_selection_changed() -- called before
            starting a new run, when results are cleared, and when the
            window closes, so a pinned object never outlives the
            docking session that pinned it (per design: pins are scoped
            to this window's own lifetime, not a persistent user choice
            like the trajectory player's own frame).
        """
        if self._receptor_vobject is not None:
            self._receptor_vobject.unpin_frame()
        for result in self._docking_results:
            if result.ligand_object is not None:
                result.ligand_object.unpin_frame()
        self._docking_results = []

    def _resolve_receptor_frame(self, receptor_id):
        """ Maps a docking result row's own 'receptor' id (the receptor
            PDBQT's basename, e.g. from self._results) to a frame index
            in self._receptor_vobject -- the concrete pose_receptor_map
            this tool needs, built lazily from the receptor VObject's
            OWN metadata (see vismol_pdbqt_builder.load_receptor_pdbqt/
            load_receptor_ensemble_pdbqt) rather than kept as a second,
            separately-maintained copy.

            Returns -1 if the receptor was never loaded into VisMol
            (self._receptor_vobject is None), or 0 for a single (non-
            ensemble) receptor -- there is no ambiguity there, it has
            exactly one frame regardless of how many receptor_ids the
            results table itself lists.
        """
        if self._receptor_vobject is None:
            return -1
        frames_meta = self._receptor_vobject.metadata.get('frames')
        if frames_meta is None:
            return 0
        if self._receptor_id_to_frame is None:
            self._receptor_id_to_frame = {v['receptor_id']: k for k, v in frames_meta.items()}
        return self._receptor_id_to_frame.get(receptor_id, -1)

    def on_button_import_pose_clicked(self, widget):
        """ Imports EVERY pose currently in the Results table, GROUPED
            BY LIGAND -- one separate VObject per distinct ligand name,
            each with one frame per that ligand's own pose (frame 0 =
            its own best affinity). Populates every liststore_results
            row's hidden vobject_id/frame_index/receptor_vobject_id/
            receptor_frame_index columns (6-9) so
            on_treeview_results_selection_changed() can jump the 3D
            view to the right ligand+receptor frame combination just by
            selecting a row afterward.

            An EARLIER version of this combined EVERY row in the whole
            table into ONE object regardless of ligand name -- correct
            back when the only way several differently-named "ligand"
            rows could appear together was an ensemble of conformers of
            the SAME molecule (several PDBQT files, one chemical
            entity). That assumption broke once ligands could be
            queued/imported independently (see the ligand-list feature
            above): a batch can now genuinely contain several UNRELATED
            ligands, and force-combining them into one object made no
            chemical sense -- confirmed by the user directly. Grouping
            by ligand name restores per-ligand-ensemble combining
            (several conformers sharing one name still merge into one
            trajectory, exactly as before) while giving every distinct
            ligand its own object.

            All poses within ONE group must share the same atom count
            to be combined into that group's object (they do whenever
            they really are conformers/receptor variants of the SAME
            ligand) -- a mismatch WITHIN a group (unexpected -- it
            shares a name, but not a topology) stops combining further
            poses into THAT group only, reported afterward, without
            affecting any other group's own object.

            No temporary PDB file is written anywhere in this method --
            every pose's atoms/coordinates come directly from parsing
            its own docked PDBQT's MODEL block (pdbqt_parser.py), and
            each group's first (topology-defining) pose is built into a
            VObject via vismol_pdbqt_builder.load_pdbqt_as_vobject() (an
            in-memory pDynamo System, see that module's docstring for
            why a real System is still needed even with no file
            involved).
        """
        if not self._results:
            self.main.simple_dialog.info(msg='No results to import yet -- run the docking first.')
            return

        # . Group row indices by ligand name, preserving first-seen
        #   order (purely cosmetic -- only affects the order objects
        #   are created in / listed in the status text).
        groups = {}
        group_order = []
        for i, row in enumerate(self._results):
            key = row['ligand']
            if key not in groups:
                groups[key] = []
                group_order.append(key)
            groups[key].append(i)

        # . Cache: one parse per distinct docked_pdbqt file, since many
        #   rows (different modes of the same job) share one file.
        models_by_file = {}
        def get_model(row):
            """ None (not an exception) if `row`'s own mode has no
                matching MODEL in its docked_pdbqt -- confirmed real:
                AutoDock Vina 1.2.3 can print a mode in its results
                table that it never actually writes to --out (a
                severely clashing pose, confirmed via an absurd
                +142 kcal/mol affinity, silently dropped during Vina's
                own write-out step) -- vina_runner.run_single_docking_job()
                now filters those out of self._results BEFORE they ever
                reach this table, but this stays defensive (None,
                handled as a skip) rather than crashing via an
                uncaught StopIteration, in case any other engine/path
                ever produces the same kind of log/file mismatch.
            """
            path = row['docked_pdbqt']
            if path not in models_by_file:
                models_by_file[path] = pdbqt_parser.parse_pdbqt_models(path)
            models = models_by_file[path]
            return next((m for m in models if m['model'] == row['mode']), None)

        all_skipped = []
        imported_summaries = []

        for ligand_name in group_order:
            group_indices = sorted(groups[ligand_name], key=lambda i: self._results[i]['affinity'])

            best_row = self._results[group_indices[0]]
            object_name = '{}_poses'.format(ligand_name)
            first_model = get_model(best_row)
            if first_model is None:
                all_skipped.append(
                    '{}: best-affinity pose (mode {}) has no matching MODEL in its PDBQT -- '
                    'nothing imported for this ligand.'.format(ligand_name, best_row['mode']))
                continue
            try:
                # . dedupe_policy='rename' -- a ligand PDBQT can
                #   legitimately have several genuinely distinct atoms
                #   sharing one generic name (confirmed real: a 5-atom
                #   ligand with 3 plain "C" atoms), so collisions here
                #   must be disambiguated by renaming, never dropped
                #   (dropping is the RECEPTOR policy, the default --
                #   see vismol_pdbqt_builder.py's docstrings).
                vobject = vismol_pdbqt_builder.load_pdbqt_as_vobject(
                    self.p_session, pdbqt_parser, first_model['atoms'], object_name,
                    tag='vina_pose', dedupe_policy='rename')
            except Exception as error:
                traceback.print_exc()
                all_skipped.append('{}: could not import ({}).'.format(ligand_name, error))
                continue

            # . rename_duplicate_atom_names() (inside load_pdbqt_as_vobject())
            #   works on a COPY, so first_model['atoms'] itself is
            #   still the pre-dedupe topology here -- fine,
            #   pdbqt_topology_names_match() below only cares about
            #   relative consistency between poses, not about matching
            #   vobject.atoms' (possibly renamed) names.
            reference_topology = first_model['atoms']
            imported_indices = [group_indices[0]]
            for pose_i in group_indices[1:]:
                row = self._results[pose_i]
                model = get_model(row)
                if model is None:
                    all_skipped.append('{} vs {}, mode {} (no MODEL for this mode in its PDBQT)'.format(
                        row['ligand'], row['receptor'], row['mode']))
                    continue
                this_topology = model['atoms']
                if len(this_topology) != len(reference_topology):
                    all_skipped.append(
                        '{} vs {}, mode {}: atom count mismatch within the "{}" group ({} vs {}) -- '
                        'stopped combining further poses into this object; {} frame(s) kept.'.format(
                            row['ligand'], row['receptor'], row['mode'], ligand_name,
                            len(reference_topology), len(this_topology), len(imported_indices)))
                    break
                if not pdbqt_parser.pdbqt_topology_names_match(reference_topology, this_topology):
                    # . Same atom count, different order -- almost
                    #   always the differing-conformer-file case
                    #   (different atom ordering across independently-
                    #   prepared ligand files), not a different
                    #   molecule. Skip only THIS pose and keep trying
                    #   the rest, rather than aborting the whole group.
                    all_skipped.append('{} vs {}, mode {} (atom order mismatch)'.format(
                        row['ligand'], row['receptor'], row['mode']))
                    continue
                coords = pdbqt_parser.get_coords_array(model)
                vismol_pdbqt_builder.append_coord_frame(vobject, coords)
                imported_indices.append(pose_i)

            # . Representations (bonds, QC highlighting, ...) are built
            #   from frame 0 at System-creation time above --
            #   reapplying them now that every frame is in place avoids
            #   the dynamic_bonds-vs-frames length mismatch a QC/
            #   dynamic representation built too early can cause (see
            #   [[project_vismol_dynamic_bonds_bug]]).
            self.p_session._apply_QC_representation_to_vobject(vismol_object=vobject)

            vobject.metadata['type'] = 'ligand_poses'
            vobject.metadata['poses'] = {}
            new_results = []
            for frame_index, i in enumerate(imported_indices):
                row = self._results[i]
                receptor_frame = self._resolve_receptor_frame(row['receptor'])
                vobject.metadata['poses'][frame_index] = {
                    'pose_rank': frame_index, 'score': row['affinity'],
                    'receptor_id': row['receptor'], 'ligand_id': row['ligand'],
                    'docked_pdbqt': row['docked_pdbqt'],
                }
                self.liststore_results[i][6] = vobject.index
                self.liststore_results[i][7] = frame_index
                receptor_vobject_id = self._receptor_vobject.index if (
                    self._receptor_vobject is not None and receptor_frame >= 0) else -1
                self.liststore_results[i][8] = receptor_vobject_id
                self.liststore_results[i][9] = receptor_frame
                if receptor_vobject_id >= 0:
                    receptor_text = '{} [frame {}]'.format(row['receptor'], receptor_frame)
                    self.liststore_results[i][0] = receptor_text
                new_results.append(DockingResult(
                    ligand_object=vobject, ligand_frame=frame_index,
                    receptor_object=self._receptor_vobject if receptor_vobject_id >= 0 else None,
                    receptor_frame=receptor_frame if receptor_vobject_id >= 0 else -1,
                    pose_rank=frame_index, score=row['affinity'], receptor_id=row['receptor']))
            self._docking_results.extend(new_results)

            imported_summaries.append('{}: {} pose(s)'.format(object_name, len(imported_indices)))

        if all_skipped:
            self.main.simple_dialog.info(
                msg='{} issue(s) while importing:\n\n{}'.format(len(all_skipped), '\n'.join(all_skipped)))

        if imported_summaries:
            status_text = '{} object(s) imported -- {}.'.format(len(imported_summaries), '; '.join(imported_summaries))
        else:
            status_text = 'Nothing was imported -- see dialog for details.'
        status_text += ' Select a row to jump to its frame.'
        self.label_docking_status.set_text(status_text)

    def on_treeview_results_selection_changed(self, selection):
        """ Once a ligand's poses have been imported (see
            on_button_import_pose_clicked()), selecting any of its rows
            pins the ligand object to its own pose frame AND, if the
            receptor was also loaded into VisMol, pins the receptor
            object to ITS OWN matching frame -- these are independent
            VismolObjects with independent frame numbers (see
            VismolObject.pin_frame(), which overrides only the object
            it's called on; every other, unrelated object in the scene
            keeps following the session-global vm_session.frame
            normally). Never assumes ligand_frame == receptor_frame.
        """
        model, treeiter = selection.get_selected()
        if treeiter is None:
            return
        vobject_id, frame_index, r_vobject_id, r_frame_index = (
            model[treeiter][6], model[treeiter][7], model[treeiter][8], model[treeiter][9])
        if vobject_id >= 0 and frame_index >= 0 and vobject_id in self.vm_session.vm_objects_dic:
            self.vm_session.vm_objects_dic[vobject_id].pin_frame(frame_index)
        if r_vobject_id >= 0 and r_frame_index >= 0 and r_vobject_id in self.vm_session.vm_objects_dic:
            self.vm_session.vm_objects_dic[r_vobject_id].pin_frame(r_frame_index)
        # . Same redraw-trigger pair vm_session.set_frame() itself uses
        #   (eSession.py) -- pin_frame() only sets an attribute, so
        #   nothing would repaint the glArea without this.
        self.vm_session.vm_glcore.updated_coords = True
        self.vm_session.vm_widget.queue_draw()
