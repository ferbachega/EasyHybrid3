#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Run SMD window
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
#      "Run SMD" window -- Steered Molecular Dynamics (NAMD's own
#      ug/node48.html) as a standalone sibling of "Run NAMD"
#      (prepare_namd_run.py), NOT a mode bolted onto that window, per
#      the user's own request. Deliberately reuses as much of that
#      window's proven machinery as possible instead of reinventing it:
#
#        - The ENTIRE "Input" tab (AMBER/CHARMM topology pickers,
#          working folder, system name, namd executable/n_procs) is
#          copied verbatim from prepare_namd_run_window.glade -- same
#          widget ids, so every handler below that touches it
#          (on_input_type_radio_toggled, on_amber_coords_source_changed,
#          _validate_common_inputs, _validate_input_files, ...) is
#          copied close to verbatim from prepare_namd_run.py too.
#        - namd_runner.build_namd_config() (extended this session with
#          `smd`/`smd_*`/`fixed_atoms_file` kwargs specifically for this
#          window), compute_cell_from_coordinates(),
#          namd_runner.run_namd(), check_namd_log_success(), and
#          process_manager_window.py's psutil process-tree helpers are
#          all reused directly -- no parallel implementation of any of
#          that.
#        - Only ever runs ONE job (no pipeline concept at all -- SMD is
#          inherently a single continuous pull), so _start_process/
#          _poll_job/_tail_log/_append_log/Stop are simplified,
#          single-job versions of prepare_namd_run.py's own methods.
#
#      What's actually NEW here: the "Steered MD" tab. The pulled group
#      (SMDFile) and the optional fixed-atoms anchor (fixedAtomsFile)
#      are both built from EasyHybrid's own EXISTING named selections
#      (`p_session.psystem[id].e_selections`, the same dict the
#      Selections window manages) via namd_runner.flag_pdb_occupancy()
#      -- no new atom-picking UI was built for this, on purpose (the
#      user's own call: build the selection first with EasyHybrid's
#      existing tools, then just pick it by name here). The pulling
#      direction (SMDDir) can EITHER be typed by hand OR filled by
#      clicking "Pick 2 atoms", which reads vismol's own picking
#      mechanism (vm_session.picking_selections, the same one already
#      used for on-screen distance measurement) and computes the vector
#      between the last 2 atoms clicked.
#

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk, GLib, Pango

import os
import re
import time
import traceback

from gui.widgets.custom_widgets import FolderChooserButton
from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.windows.setup.process_manager_window import _terminate_process_tree, _force_kill_pids, _pid_alive
from pdynamo.pDynamo2EasyHybrid.helpers import export_special_PDB
from util import namd_runner


# Same per-line highlight rules as prepare_namd_run.py's own Log tab --
# see namd_runner.check_namd_log_success() for the same markers used
# there to decide pass/fail.
_LOG_ERROR_RE = re.compile(r'FATAL|ERROR:|Error!', re.IGNORECASE)
_LOG_WARNING_RE = re.compile(r'Warning', re.IGNORECASE)
_LOG_SUCCESS_RE = re.compile(r'WRITING COORDINATES|End of program|WallClock:')

_NO_FIXED_ATOMS = '(none)'


class PrepareSMDWindow:
    """ "Run SMD" window (Steered Molecular Dynamics). """

    def __init__(self, main=None):
        """ Class initialiser """
        self.main       = main
        self.vm_session = main.vm_session
        self.p_session  = main.p_session
        self.home       = main.home
        self.Visible    = False

        self.current_process   = None   # subprocess.Popen or None
        self.current_log_path  = None
        self._log_read_offset  = 0
        self._poll_timeout_id  = None
        self._input_files      = None   # dict from _validate_input_files()
        self._work_folder      = None
        self._system_name      = None

    #-------------------------------------------------------------------------
    #  W I N D O W   L I F E C Y C L E
    #-------------------------------------------------------------------------
    def open_window(self):
        """ Function doc """
        if self.Visible:
            self.window.present()
            return

        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/setup/prepare_smd_window.glade'))
        self.builder.connect_signals(self)

        self.window = self.builder.get_object('window')

        # -------------------- Input tab (copied from "Run NAMD") --------------------
        self.radio_input_amber  = self.builder.get_object('radio_input_amber')
        self.radio_input_charmm = self.builder.get_object('radio_input_charmm')
        self.box_amber_input    = self.builder.get_object('box_amber_input')
        self.box_charmm_input   = self.builder.get_object('box_charmm_input')
        self.entry_prmtop_path = self.builder.get_object('entry_prmtop_path')
        self.entry_inpcrd_path = self.builder.get_object('entry_inpcrd_path')
        self.combo_amber_coords_source = self.builder.get_object('combo_amber_coords_source')
        self.label_amber_inp = self.builder.get_object('label_amber_inp')
        self.button_browse_inpcrd = self.builder.get_object('button_browse_inpcrd')
        self.label_system = self.builder.get_object('label_system')
        self.label_object = self.builder.get_object('label_object')
        self.box_amber_coords_system = self.builder.get_object('box_amber_coords_system')
        self.box_amber_coords_object = self.builder.get_object('box_amber_coords_object')
        self.entry_charmm_psf_path = self.builder.get_object('entry_charmm_psf_path')
        self.entry_charmm_coordinates_path = self.builder.get_object('entry_charmm_coordinates_path')
        self.combo_charmm_coords_source = self.builder.get_object('combo_charmm_coords_source')
        self.box_charmm_coords_file = self.builder.get_object('box_charmm_coords_file')
        self.box_charmm_coords_vobject = self.builder.get_object('box_charmm_coords_vobject')
        self.treeview_charmm_parameters = self.builder.get_object('treeview_charmm_parameters')
        self.liststore_charmm_parameters = self.builder.get_object('liststore_charmm_parameters')
        self.checkbox_charmm_xplor = self.builder.get_object('checkbox_charmm_xplor')
        self.entry_namd_command = self.builder.get_object('entry_namd_command')
        self.spinbtn_n_procs   = self.builder.get_object('spinbtn_n_procs')
        self.entry_system_name = self.builder.get_object('entry_system_name')
        self.label_info        = self.builder.get_object('label_info')
        self.textview_log      = self.builder.get_object('textview_log')
        log_buffer = self.textview_log.get_buffer()
        self._log_tag_error = log_buffer.create_tag('namd_log_error', foreground='#e01b24', weight=Pango.Weight.BOLD)
        self._log_tag_warning = log_buffer.create_tag('namd_log_warning', foreground='#c47a00')
        self._log_tag_success = log_buffer.create_tag('namd_log_success', foreground='#26a269', weight=Pango.Weight.BOLD)
        self.button_stop        = self.builder.get_object('button_stop')

        self.amber_coords_system_combo = SystemComboBox(self.main)
        self.amber_coords_system_combo.connect("changed", self._on_amber_coords_system_changed)
        self.box_amber_coords_system.pack_start(self.amber_coords_system_combo, False, False, 0)
        self.amber_coords_object_combo = CoordinatesComboBox()
        self.box_amber_coords_object.pack_start(self.amber_coords_object_combo, False, False, 0)

        self.charmm_coords_system_combo = SystemComboBox(self.main)
        self.charmm_coords_system_combo.connect("changed", self._on_charmm_coords_system_changed)
        self.builder.get_object('box_charmm_coords_system').pack_start(self.charmm_coords_system_combo, False, False, 0)
        self.charmm_coords_object_combo = CoordinatesComboBox()
        self.builder.get_object('box_charmm_coords_object').pack_start(self.charmm_coords_object_combo, False, False, 0)

        self.box_work_folder = self.builder.get_object('box_work_folder')
        self.work_folder_chooser = FolderChooserButton(self.main, 'folder', self.home)
        self.box_work_folder.pack_start(self.work_folder_chooser.btn, False, False, 0)
        scratch = os.environ.get('PDYNAMO3_SCRATCH')
        if scratch and os.path.isdir(scratch):
            self.work_folder_chooser.set_folder(folder=scratch)

        namd_command = self.main.vm_session.vm_config.gl_parameters.get('namd_command')
        if not namd_command or not os.path.isfile(namd_command):
            namd_command = namd_runner.find_namd_executable()
        self.entry_namd_command.set_text(namd_command or '')
        if not namd_command:
            self.label_info.set_text(
                'namd3/namd2 was not found automatically on PATH -- please enter the full path below.')

        self.charmm_parameters = []

        # -------------------- Steered MD tab --------------------
        self.spinbtn_smd_temperature = self.builder.get_object('spinbtn_smd_temperature')
        self.spinbtn_smd_timestep    = self.builder.get_object('spinbtn_smd_timestep')
        self.spinbtn_smd_cutoff      = self.builder.get_object('spinbtn_smd_cutoff')
        self.combo_smd_watermodel    = self.builder.get_object('combo_smd_watermodel')
        self.radio_smd_ensemble_nvt  = self.builder.get_object('radio_smd_ensemble_nvt')
        self.radio_smd_ensemble_npt  = self.builder.get_object('radio_smd_ensemble_npt')
        self.spinbtn_smd_run_steps   = self.builder.get_object('spinbtn_smd_run_steps')

        self.box_smd_selection_system = self.builder.get_object('box_smd_selection_system')
        self.box_smd_selection_object = self.builder.get_object('box_smd_selection_object')
        self.combo_smd_pulled_group = self.builder.get_object('combo_smd_pulled_group')
        self.combo_smd_fixed_atoms  = self.builder.get_object('combo_smd_fixed_atoms')

        self.spinbtn_smd_k = self.builder.get_object('spinbtn_smd_k')
        self.spinbtn_smd_k2 = self.builder.get_object('spinbtn_smd_k2')
        self.spinbtn_smd_vel = self.builder.get_object('spinbtn_smd_vel')
        self.spinbtn_smd_output_freq = self.builder.get_object('spinbtn_smd_output_freq')
        self.spinbtn_smd_dir_x = self.builder.get_object('spinbtn_smd_dir_x')
        self.spinbtn_smd_dir_y = self.builder.get_object('spinbtn_smd_dir_y')
        self.spinbtn_smd_dir_z = self.builder.get_object('spinbtn_smd_dir_z')
        self.button_run_smd = self.builder.get_object('button_run_smd')

        self.smd_selection_system_combo = SystemComboBox(self.main)
        self.smd_selection_system_combo.connect("changed", self._on_smd_selection_system_changed)
        self.box_smd_selection_system.pack_start(self.smd_selection_system_combo, False, False, 0)
        self.smd_selection_object_combo = CoordinatesComboBox()
        self.box_smd_selection_object.pack_start(self.smd_selection_object_combo, False, False, 0)

        self.window.show_all()
        # . show_all() forces every hidden widget momentarily visible --
        #   same re-sync-after-the-fact pattern prepare_namd_run.py's
        #   own open_window() already documents/relies on.
        self.on_input_type_radio_toggled(self.radio_input_amber)
        self.on_amber_coords_source_changed(self.combo_amber_coords_source)
        self.box_charmm_coords_file.set_visible(self.combo_charmm_coords_source.get_active_id() != 'vobject')
        self.box_charmm_coords_vobject.set_visible(self.combo_charmm_coords_source.get_active_id() == 'vobject')
        self._refresh_smd_selection_combos()
        self.window.reshow_with_initial_size()
        self.window.connect('destroy', self.close_window)
        self.Visible = True

    def close_window(self, button=None, data=None):
        """ Function doc """
        if not self.Visible:
            return
        self.window.destroy()
        self.Visible = False

    def on_button_cancel_clicked(self, widget):
        self.close_window()

    #-------------------------------------------------------------------------
    #  H E L P E R S  ( c o p i e d   f r o m   " R u n   N A M D " )
    #-------------------------------------------------------------------------
    def _browse_file(self, entry, title, patterns, pattern_name):
        dialog = Gtk.FileChooserDialog(
            title=title, parent=self.window, action=Gtk.FileChooserAction.OPEN,
        )
        dialog.add_buttons(
            Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
            Gtk.STOCK_OPEN, Gtk.ResponseType.OK,
        )
        if patterns:
            file_filter = Gtk.FileFilter()
            file_filter.set_name(pattern_name)
            for pattern in patterns:
                file_filter.add_pattern(pattern)
            dialog.add_filter(file_filter)
        if dialog.run() == Gtk.ResponseType.OK:
            entry.set_text(dialog.get_filename())
        dialog.destroy()

    def _append_log(self, text):
        buf = self.textview_log.get_buffer()
        for line in text.splitlines(keepends=True):
            if _LOG_ERROR_RE.search(line):
                tag = self._log_tag_error
            elif _LOG_WARNING_RE.search(line):
                tag = self._log_tag_warning
            elif _LOG_SUCCESS_RE.search(line):
                tag = self._log_tag_success
            else:
                tag = None
            end_iter = buf.get_end_iter()
            if tag is not None:
                buf.insert_with_tags(end_iter, line, tag)
            else:
                buf.insert(end_iter, line)
        self.textview_log.scroll_to_iter(buf.get_end_iter(), 0.0, False, 0.0, 0.0)

    def _job_running(self):
        return self.current_process is not None

    def _set_controls_sensitive(self, sensitive):
        self.button_run_smd.set_sensitive(sensitive)
        self.button_stop.set_sensitive(not sensitive)

    #-------------------------------------------------------------------------
    #  S I G N A L S :  A M B E R   I N P U T   ( c o p i e d )
    #-------------------------------------------------------------------------
    def on_button_browse_prmtop_clicked(self, widget):
        self._browse_file(self.entry_prmtop_path, "Select AMBER topology",
                           ["*.prmtop", "*.top", "*.parm7"], "AMBER topology (*.prmtop, *.top, *.parm7)")

    def on_entry_prmtop_path_changed(self, widget):
        # . No manual-periodic-cell UI in this window (unlike "Run
        #   NAMD"'s Setup tab) -- nothing to prefill from the prmtop's
        #   own BOX_DIMENSIONS, so this is a deliberate no-op. Kept as a
        #   real handler (not removed from the glade) purely because
        #   the copied Input tab fragment already wires this signal.
        pass

    def on_button_browse_inpcrd_clicked(self, widget):
        self._browse_file(self.entry_inpcrd_path, "Select AMBER coordinates",
                           ["*.inpcrd", "*.crd", "*.rst7"], "AMBER coordinates (*.inpcrd, *.crd, *.rst7)")

    #-------------------------------------------------------------------------
    #  S I G N A L S :  C H A R M M   I N P U T   ( c o p i e d )
    #-------------------------------------------------------------------------
    def on_input_type_radio_toggled(self, widget):
        if not widget.get_active():
            return   # fires for both the button losing AND gaining active state
        self.box_amber_input.set_visible(self.radio_input_amber.get_active())
        self.box_charmm_input.set_visible(self.radio_input_charmm.get_active())

    def on_amber_coords_source_changed(self, widget):
        is_vobject = widget.get_active_id() == 'vobject'
        self.label_amber_inp.set_visible(not is_vobject)
        self.entry_inpcrd_path.set_visible(not is_vobject)
        self.button_browse_inpcrd.set_visible(not is_vobject)

        self.label_system.set_visible(is_vobject)
        self.label_object.set_visible(is_vobject)
        self.box_amber_coords_system.set_visible(is_vobject)
        self.box_amber_coords_object.set_visible(is_vobject)

    def on_charmm_coords_source_changed(self, widget):
        is_vobject = widget.get_active_id() == 'vobject'
        self.box_charmm_coords_file.set_visible(not is_vobject)
        self.box_charmm_coords_vobject.set_visible(is_vobject)

    def _on_amber_coords_system_changed(self, widget):
        system_id = self.amber_coords_system_combo.get_system_id()
        if system_id is not None:
            self.amber_coords_object_combo.set_model(self.main.vobject_liststore_dict[system_id])
            size = len(list(self.main.vobject_liststore_dict[system_id]))
            self.amber_coords_object_combo.set_active(size - 1)

    def _on_charmm_coords_system_changed(self, widget):
        system_id = self.charmm_coords_system_combo.get_system_id()
        if system_id is not None:
            self.charmm_coords_object_combo.set_model(self.main.vobject_liststore_dict[system_id])
            size = len(list(self.main.vobject_liststore_dict[system_id]))
            self.charmm_coords_object_combo.set_active(size - 1)

    def on_button_browse_charmm_psf_clicked(self, widget):
        self._browse_file(self.entry_charmm_psf_path, "Select CHARMM psf", ["*.psf"], "CHARMM topology (*.psf)")

    def on_button_browse_charmm_coordinates_clicked(self, widget):
        self._browse_file(self.entry_charmm_coordinates_path, "Select starting coordinates",
                           ["*.pdb"], "Coordinates (*.pdb)")

    def on_button_add_charmm_parameter_clicked(self, widget):
        dialog = Gtk.FileChooserDialog(
            title="Add CHARMM parameter file(s) (.prm/.par/.str/.xplor)",
            parent=self.window, action=Gtk.FileChooserAction.OPEN,
        )
        dialog.set_select_multiple(True)
        dialog.add_buttons(
            Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
            Gtk.STOCK_OPEN, Gtk.ResponseType.OK,
        )
        file_filter = Gtk.FileFilter()
        file_filter.set_name("CHARMM parameter files (*.prm, *.par, *.str, *.xplor, *.rtf)")
        for pattern in ("*.prm", "*.par", "*.str", "*.xplor", "*.rtf"):
            file_filter.add_pattern(pattern)
        dialog.add_filter(file_filter)

        if dialog.run() == Gtk.ResponseType.OK:
            for path in dialog.get_filenames():
                self.charmm_parameters.append(path)
                self.liststore_charmm_parameters.append([os.path.basename(path)])
        dialog.destroy()

    def on_button_clear_charmm_parameters_clicked(self, widget):
        self.charmm_parameters = []
        self.liststore_charmm_parameters.clear()

    #-------------------------------------------------------------------------
    #  C O M M O N   V A L I D A T I O N   ( c o p i e d )
    #-------------------------------------------------------------------------
    def _validate_common_inputs(self):
        if self._job_running():
            self.main.simple_dialog.info(msg='A job is already running -- stop it first.')
            return None

        namd_command = self.entry_namd_command.get_text().strip()
        if not namd_command or not (os.path.isfile(namd_command) and os.access(namd_command, os.X_OK)):
            self.main.simple_dialog.info(msg='Please provide a valid path to the namd3/namd2 executable.')
            return None
        self.main.vm_session.vm_config.gl_parameters['namd_command'] = namd_command

        work_folder = self.work_folder_chooser.get_folder()
        if not work_folder:
            self.main.simple_dialog.info(msg='Please choose a working folder.')
            return None

        system_name = self.entry_system_name.get_text().strip() or 'smd_run'

        return namd_command, work_folder, system_name

    def _get_input_type(self):
        return 'charmm' if self.radio_input_charmm.get_active() else 'amber'

    def _get_selected_vobject(self, system_combo, object_combo):
        vobject_id = object_combo.get_vobject_id()
        if vobject_id is None:
            return None
        return self.main.vm_session.vm_objects_dic.get(vobject_id)

    def _write_amber_crd_from_vobject(self, vobject, output_path):
        positions = [vobject.atoms[index].coords(-1) for index in vobject.atoms.keys()]
        namd_runner.write_amber_crd(positions, output_path, title=vobject.name)

    def _validate_input_files(self, work_folder):
        if self._get_input_type() == 'charmm':
            psf = self.entry_charmm_psf_path.get_text().strip()
            parameters = list(self.charmm_parameters)
            if not psf or not os.path.isfile(psf):
                self.main.simple_dialog.info(msg='Please select a valid CHARMM psf file.')
                return None
            if not parameters:
                self.main.simple_dialog.info(msg='Please add at least one CHARMM parameter file.')
                return None

            if self.combo_charmm_coords_source.get_active_id() == 'vobject':
                vobject = self._get_selected_vobject(self.charmm_coords_system_combo, self.charmm_coords_object_combo)
                if vobject is None:
                    self.main.simple_dialog.info(msg='Please select a system/object for the coordinates.')
                    return None
                coordinates = os.path.join(work_folder, '_smd_coords_from_vobject.pdb')
                try:
                    export_special_PDB(vobject=vobject, frame=-1, output=coordinates)
                except Exception as error:
                    traceback.print_exc()
                    self.main.simple_dialog.info(msg='Could not export coordinates from the selected object:\n{}'.format(error))
                    return None
            else:
                coordinates = self.entry_charmm_coordinates_path.get_text().strip()
                if not coordinates or not os.path.isfile(coordinates):
                    self.main.simple_dialog.info(msg='Please select a valid coordinates (pdb) file.')
                    return None

            return {
                'type': 'charmm', 'psf': psf, 'parameters': parameters,
                'coordinates': coordinates, 'xplor': self.checkbox_charmm_xplor.get_active(),
            }

        prmtop = self.entry_prmtop_path.get_text().strip()
        if not prmtop or not os.path.isfile(prmtop):
            self.main.simple_dialog.info(msg='Please select a valid AMBER prmtop file.')
            return None

        if self.combo_amber_coords_source.get_active_id() == 'vobject':
            vobject = self._get_selected_vobject(self.amber_coords_system_combo, self.amber_coords_object_combo)
            if vobject is None:
                self.main.simple_dialog.info(msg='Please select a system/object for the coordinates.')
                return None
            inpcrd = os.path.join(work_folder, '_smd_coords_from_vobject.crd')
            try:
                self._write_amber_crd_from_vobject(vobject, inpcrd)
            except Exception as error:
                traceback.print_exc()
                self.main.simple_dialog.info(msg='Could not export coordinates from the selected object:\n{}'.format(error))
                return None
        else:
            inpcrd = self.entry_inpcrd_path.get_text().strip()
            if not inpcrd or not os.path.isfile(inpcrd):
                self.main.simple_dialog.info(msg='Please select a valid AMBER inpcrd file.')
                return None

        return {'type': 'amber', 'prmtop': prmtop, 'inpcrd': inpcrd}

    #-------------------------------------------------------------------------
    #  S I G N A L S :  S T E E R E D   M D   T A B
    #-------------------------------------------------------------------------
    def _on_smd_selection_system_changed(self, widget):
        """ Mirrors _on_amber_coords_system_changed() -- populates the
            "Object / frame" combo for whichever system was just picked
            AND refreshes the pulled/fixed-group combos from that
            system's own p_session.psystem[id].e_selections (a plain
            {name: [atom indexes]} dict the Selections window already
            manages -- see this module's own docstring for why no new
            atom-picking UI was built here instead).
        """
        system_id = self.smd_selection_system_combo.get_system_id()
        if system_id is not None:
            self.smd_selection_object_combo.set_model(self.main.vobject_liststore_dict[system_id])
            size = len(list(self.main.vobject_liststore_dict[system_id]))
            self.smd_selection_object_combo.set_active(size - 1)
        self._refresh_smd_selection_combos()

    def _refresh_smd_selection_combos(self):
        self.combo_smd_pulled_group.remove_all()
        self.combo_smd_fixed_atoms.remove_all()
        self.combo_smd_fixed_atoms.append_text(_NO_FIXED_ATOMS)
        self.combo_smd_fixed_atoms.set_active(0)

        system_id = self.smd_selection_system_combo.get_system_id()
        if system_id is None:
            return
        system = self.p_session.psystem.get(system_id)
        if system is None:
            return
        for name in system.e_selections.keys():
            self.combo_smd_pulled_group.append_text(name)
            self.combo_smd_fixed_atoms.append_text(name)
        if len(system.e_selections):
            self.combo_smd_pulled_group.set_active(0)

    def on_button_smd_pick_direction_clicked(self, widget):
        """ Reads the last 2 atoms clicked (picked) in the 3D view --
            vismol's own picking_selections_list, the same mechanism
            already used for on-screen distance measurement -- and
            fills the X/Y/Z fields with the vector between them (second
            minus first). Does NOT normalize -- NAMD's own SMDDir
            doesn't require it either (see build_namd_config()'s smd_dir
            docstring).
        """
        picking = self.vm_session.picking_selections
        atom1 = picking.picking_selections_list[0]
        atom2 = picking.picking_selections_list[1]
        if atom1 is None or atom2 is None:
            self.main.simple_dialog.info(msg='Pick (click) 2 atoms in the 3D view first, then try again.')
            return
        xyz1 = atom1.coords()
        xyz2 = atom2.coords()
        self.spinbtn_smd_dir_x.set_value(xyz2[0] - xyz1[0])
        self.spinbtn_smd_dir_y.set_value(xyz2[1] - xyz1[1])
        self.spinbtn_smd_dir_z.set_value(xyz2[2] - xyz1[2])

    #-------------------------------------------------------------------------
    #  R U N
    #-------------------------------------------------------------------------
    def on_button_run_smd_clicked(self, widget):
        """ Function doc """
        validated = self._validate_common_inputs()
        if validated is None:
            return
        namd_command, work_folder, system_name = validated

        input_files = self._validate_input_files(work_folder)
        if input_files is None:
            return
        is_charmm = input_files['type'] == 'charmm'

        selection_vobject = self._get_selected_vobject(
            self.smd_selection_system_combo, self.smd_selection_object_combo)
        if selection_vobject is None:
            self.main.simple_dialog.info(msg='Please select a system/object to read selections from.')
            return
        selection_system_id = self.smd_selection_system_combo.get_system_id()
        e_selections = self.p_session.psystem[selection_system_id].e_selections

        pulled_name = self.combo_smd_pulled_group.get_active_text()
        if not pulled_name:
            self.main.simple_dialog.info(msg='Please select a pulled group (build one first with the Selections tool).')
            return
        pulled_indexes = e_selections[pulled_name]

        fixed_name = self.combo_smd_fixed_atoms.get_active_text()
        fixed_indexes = None if (not fixed_name or fixed_name == _NO_FIXED_ATOMS) else e_selections[fixed_name]

        dx = self.spinbtn_smd_dir_x.get_value()
        dy = self.spinbtn_smd_dir_y.get_value()
        dz = self.spinbtn_smd_dir_z.get_value()
        if dx == 0.0 and dy == 0.0 and dz == 0.0:
            self.main.simple_dialog.info(msg='The pulling direction is (0, 0, 0) -- set it by hand or "Pick 2 atoms" first.')
            return

        try:
            full_pdb = os.path.join(work_folder, system_name + '_smd_full.pdb')
            export_special_PDB(vobject=selection_vobject, frame=-1, output=full_pdb)

            smd_pdb = os.path.join(work_folder, system_name + '_smd_pulled.pdb')
            namd_runner.flag_pdb_occupancy(full_pdb, pulled_indexes, smd_pdb)

            fixed_pdb = None
            if fixed_indexes is not None:
                fixed_pdb = os.path.join(work_folder, system_name + '_smd_fixed.pdb')
                namd_runner.flag_pdb_occupancy(full_pdb, fixed_indexes, fixed_pdb)
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg='Could not build the SMD/fixed-atoms PDB(s):\n{}'.format(error))
            return

        try:
            if is_charmm:
                cell = namd_runner.compute_cell_from_charmm(input_files['coordinates'])
            else:
                cell = namd_runner.compute_cell_from_coordinates(input_files['prmtop'], input_files['inpcrd'])
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg='Could not compute the periodic cell from the input coordinates:\n{}'.format(error))
            return

        input_kwargs = (
            {
                'charmm_psf': input_files['psf'], 'charmm_parameters': input_files['parameters'],
                'charmm_xplor': input_files['xplor'], 'coordinates': input_files['coordinates'],
            } if is_charmm else
            {'parmfile': input_files['prmtop'], 'ambercoor': input_files['inpcrd']}
        )

        timestep = self.spinbtn_smd_timestep.get_value()
        # . SMDVel is Å/TIMESTEP in NAMD's own config, not the more
        #   intuitive Å/ps this window's own field shows -- see
        #   build_namd_config()'s own smd_vel docstring for why the
        #   conversion happens HERE, not inside that function.
        smd_vel_per_timestep = self.spinbtn_smd_vel.get_value() * (timestep / 1000.0)

        config_text = namd_runner.build_namd_config(
            output_basename=system_name, **input_kwargs,
            cell=cell,
            temperature=self.spinbtn_smd_temperature.get_value(),
            timestep=timestep,
            cutoff=self.spinbtn_smd_cutoff.get_value(),
            watermodel=self.combo_smd_watermodel.get_active_id() or 'tip3',
            run_steps=int(self.spinbtn_smd_run_steps.get_value()),
            ensemble='NPT' if self.radio_smd_ensemble_npt.get_active() else 'NVT',
            langevin=True,
            fixed_atoms_file=fixed_pdb,
            smd=True,
            smd_file=smd_pdb,
            smd_k=self.spinbtn_smd_k.get_value(),
            smd_k2=self.spinbtn_smd_k2.get_value(),
            smd_vel=smd_vel_per_timestep,
            smd_dir=(dx, dy, dz),
            smd_output_freq=int(self.spinbtn_smd_output_freq.get_value()),
        )

        config_path = os.path.join(work_folder, system_name + '.namd')
        with open(config_path, 'w') as f:
            f.write(config_text)
        log_path = os.path.join(work_folder, system_name + '.log')

        self._input_files = input_files
        self._work_folder = work_folder
        self._system_name = system_name

        self._start_process(config_path, work_folder, namd_command, log_path,
                             int(self.spinbtn_n_procs.get_value()))
        self.label_info.set_text('Running SMD ("{}")...'.format(system_name))

    #-------------------------------------------------------------------------
    #  P R O C E S S   L I F E C Y C L E   ( s i n g l e - j o b   o n l y )
    #-------------------------------------------------------------------------
    def _start_process(self, config_path, work_folder, namd_command, log_path, n_procs):
        self.current_process = namd_runner.run_namd(config_path, work_folder, namd_command, log_path, n_procs=n_procs)
        self.current_log_path = log_path
        self._log_read_offset = 0
        buf = self.textview_log.get_buffer()
        buf.set_text('')
        self._set_controls_sensitive(False)
        if self._poll_timeout_id is None:
            self._poll_timeout_id = GLib.timeout_add(500, self._poll_job)

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
            polling pattern process_manager_window.py/prepare_namd_run.py
            already use. Returns True to keep polling, False to stop.
        """
        self._tail_log()

        if self.current_process is None:
            return False

        returncode = self.current_process.poll()
        if returncode is None:
            return True   # still running

        ok = (returncode == 0) and namd_runner.check_namd_log_success(self.current_log_path)
        self.current_process = None
        self._poll_timeout_id = None
        self._set_controls_sensitive(True)

        if not ok:
            self.label_info.set_text('SMD run failed -- see the log above.')
            return False

        self.label_info.set_text('SMD run finished -- importing the result back into EasyHybrid...')
        dcd_path = os.path.join(self._work_folder, self._system_name + '.dcd')
        self._import_result_back(self._input_files, [dcd_path], self._system_name)
        return False

    #-------------------------------------------------------------------------
    #  I M P O R T   B A C K   ( c o p i e d ,   s i n g l e - d c d )
    #-------------------------------------------------------------------------
    def _import_result_back(self, input_files, dcd_paths, system_name):
        try:
            if input_files['type'] == 'charmm':
                self.p_session.load_a_new_pDynamo_system_from_dict(
                    input_files={
                        'charmm_psf': input_files['psf'], 'charmm_par': input_files['parameters'],
                        'coordinates': input_files['coordinates'],
                    },
                    system_type=1, name=system_name, tag='SMD',
                    working_folder=self._work_folder,
                )
            else:
                top_path, crd_path = namd_runner.ensure_amber_top_crd(
                    input_files['prmtop'], input_files['inpcrd'], self._work_folder)
                self.p_session.load_a_new_pDynamo_system_from_dict(
                    input_files={'amber_prmtop': top_path, 'coordinates': crd_path},
                    system_type=0, name=system_name, tag='SMD',
                    working_folder=self._work_folder,
                )
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg='SMD finished, but the system could not be imported back into '
                                              'EasyHybrid:\n{}'.format(error))
            self.label_info.set_text('SMD finished, but importing the result failed -- see the message above.')
            return

        existing_dcd_paths = [p for p in dcd_paths if os.path.isfile(p)]
        if not existing_dcd_paths:
            self.label_info.set_text(
                'SMD finished and "{}" was added, but no trajectory (.dcd) frame was written '
                '(run too short for the output frequency) -- nothing to import.'.format(system_name))
            return

        system_id = list(self.p_session.psystem.keys())[-1]
        vobject = None
        try:
            for dcd_path in existing_dcd_paths:
                parameters = {
                    'system_id': system_id, 'data_path': dcd_path, 'data_type': 'dcd',
                    'new_vobj_name': None if vobject else system_name + '_traj',
                    'vobject_id': vobject.index if vobject else None, 'vobject': vobject,
                    'logfile': None, 'first': None, 'last': None, 'stride': None,
                }
                self.p_session.import_data(parameters)
                vobject = parameters['vobject']
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg='The system was imported, but the trajectory could not be loaded:\n{}'.format(error))
            self.label_info.set_text('SMD finished and "{}" was added, but loading the trajectory failed.'.format(system_name))
            return

        self.label_info.set_text('Done! "{}" and its trajectory were added to the treeview.'.format(system_name))

    #-------------------------------------------------------------------------
    #  S T O P   ( c o p i e d )
    #-------------------------------------------------------------------------
    def on_button_stop_clicked(self, widget):
        if self.current_process is None:
            return
        if not self.main.simple_dialog.question('Stop the current SMD run?'):
            return

        if self._poll_timeout_id is not None:
            GLib.source_remove(self._poll_timeout_id)
            self._poll_timeout_id = None

        pids_to_watch = _terminate_process_tree(self.current_process.pid)
        self.label_info.set_text('Stopping...')
        GLib.timeout_add(200, self._check_stop_progress, pids_to_watch, time.time())

    def _check_stop_progress(self, pids_to_watch, started_at):
        GRACE_PERIOD_SECONDS = 5.0

        still_alive = any(_pid_alive(pid) for pid in pids_to_watch)
        if still_alive and (time.time() - started_at) < GRACE_PERIOD_SECONDS:
            return True

        if still_alive:
            _force_kill_pids(pids_to_watch)

        if self.current_process is not None:
            try:
                self.current_process.wait(timeout=1)
            except Exception:
                pass
        self.current_process = None

        self._set_controls_sensitive(True)
        self.label_info.set_text('Stopped by user.')
        return False
