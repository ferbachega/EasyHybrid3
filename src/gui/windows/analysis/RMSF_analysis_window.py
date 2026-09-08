#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: RMSF (root-mean-square fluctuation) analysis window
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
#      "RMSF" analysis window -- per-atom root-mean-square fluctuation
#      over a trajectory, relative to the mean position of each atom over
#      the analysed frames. The Object picked via System/Object (same
#      SystemComboBox/CoordinatesComboBox pair used throughout the rest
#      of the app) is analysed IN FULL by default; "Limit to current
#      viewer selection" narrows that down to whatever is currently
#      selected in the 3D view instead (via
#      build_index_list_from_atom_selection(), the same mechanism
#      RMSD_analysis_window.py uses as its ONLY input) -- cross-checked
#      against the chosen Object so a selection made on some other,
#      unrelated Object can't silently be applied here.
#      The actual math lives in util/md_analysis.py (GTK-free, reusable).
#
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

import os

from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.widgets.custom_widgets import FolderChooserButton
from gui.windows.setup.windows_and_dialogs import TextWindow
from util import md_analysis

# . Standard PDB/CHARMM/AMBER protein backbone atom names. Checked
#   together with residue.is_protein (vismol's own residue classifier,
#   see graphics_engine/src/vismol/model/residue.py) rather than by name
#   alone -- a bare name check would also match e.g. a water oxygen
#   ("O") or a generic ligand atom named "C", silently pulling
#   non-backbone atoms into what is supposed to be a backbone-only set.
_BACKBONE_ATOM_NAMES = {"N", "CA", "C", "O"}


class RMSFAnalysisWindow:

    def __init__(self, main=None):
        """ Class initialiser """
        self.main       = main
        self.home       = main.home
        self.p_session  = self.main.p_session
        self.vm_session = self.main.vm_session
        self.Visible    = False

    def open_window(self):
        """ Function doc """
        if self.Visible:
            self.window.present()
            return

        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/analysis/RMSF_analysis.glade'))
        self.builder.connect_signals(self)

        self.window = self.builder.get_object('window')
        self.window.set_title('RMSF Analysis')

        self.builder.get_object('btn_close').connect('clicked', self.close_window)
        self.builder.get_object('chk_frames_from').connect('toggled', self.on_frame_range_checkbox)
        self.builder.get_object('chk_step_size').connect('toggled', self.on_step_size_checkbox)

        # . expand=True/fill=True (rather than this codebase's more usual
        #   False/False) so the combobox itself stretches to fill
        #   box_system/box_coordinates's full width -- those two GtkBoxes
        #   are already set hexpand=True in the .glade, which only grows
        #   the BOX; the child still needs its own pack_start flags to
        #   actually fill that grown space instead of staying shrink-wrapped.
        self.box_system = self.builder.get_object('box_system')
        self.combobox_systems = SystemComboBox(self.main)
        self.combobox_systems.connect("changed", self.on_combobox_systems_changed)
        self.box_system.pack_start(self.combobox_systems, True, True, 0)

        self.box_coordinates = self.builder.get_object('box_coordinates')
        self.coordinates_combobox = CoordinatesComboBox()
        self.box_coordinates.pack_start(self.coordinates_combobox, True, True, 0)

        self.grid_save_output = self.builder.get_object('grid_save_output')
        self.entry_output_filename = self.builder.get_object('entry_output_filename')
        self.builder.get_object('chk_save_file').connect('toggled', self.on_save_file_checkbox)

        self.box_output_folder = self.builder.get_object('box_output_folder')
        self.output_folder_chooser = FolderChooserButton(self.main, 'folder', self.home)
        self.box_output_folder.pack_start(self.output_folder_chooser.btn, False, False, 0)

        if self.p_session.psystem:
            self.combobox_systems.set_active_system(e_id=self.p_session.active_id)

        self.window.show_all()
        self.Visible = True

    def on_save_file_checkbox(self, widget):
        self.grid_save_output.set_sensitive(widget.get_active())

    def on_combobox_systems_changed(self, widget):
        system_id = self.combobox_systems.get_system_id()
        if system_id is None:
            return
        # . Default the output folder to this system's own working
        #   folder (system.e_working_folder, set at import time -- same
        #   attribute PES_analysis_window.py's own export already
        #   defaults to) each time the System changes, rather than only
        #   once at window-open time.
        system = self.p_session.psystem.get(system_id)
        if system is not None and getattr(system, 'e_working_folder', None):
            self.output_folder_chooser.set_folder(folder=system.e_working_folder)
        self.coordinates_combobox.set_model(self.main.vobject_liststore_dict[system_id])
        size = len(list(self.main.vobject_liststore_dict[system_id]))
        self.coordinates_combobox.set_active(size - 1)

    def close_window(self, button=None, data=None):
        """ Function doc """
        if not self.Visible:
            return
        self.window.destroy()
        self.Visible = False

    def on_step_size_checkbox(self, widget):
        self.builder.get_object('entry_step_size').set_sensitive(widget.get_active())

    def on_frame_range_checkbox(self, widget):
        active = widget.get_active()
        self.builder.get_object('entry_frame_init').set_sensitive(active)
        self.builder.get_object('entry_frame_last').set_sensitive(active)

    def _save_text_to_file(self, text, default_filename):
        """ Writes `text` directly to the folder/file name set in the
            "Save Output" frame (Folder: defaults to the chosen System's
            own e_working_folder, see on_combobox_systems_changed() --
            File name: defaults to the static text set in the .glade,
            "rmsf.txt") -- no "Save as" popup, per this window's own
            convention: the destination is a persistent, visible field
            in the window, not a one-off dialog.
        """
        folder = self.output_folder_chooser.get_folder()
        if not folder:
            self.main.simple_dialog.info(msg='Please choose an output folder.')
            return False

        filename = self.entry_output_filename.get_text().strip() or default_filename
        if not filename.lower().endswith('.txt'):
            filename += '.txt'
        filepath = os.path.join(folder, filename)

        try:
            with open(filepath, 'w') as f:
                f.write(text)
        except Exception as error:
            self.main.simple_dialog.info(msg='Could not save the file:\n{}'.format(error))
            return False
        return True

    def _filter_backbone(self, selection_list, vismol_object):
        filtered = []
        for atom_index in selection_list:
            atom = vismol_object.atoms[atom_index]
            residue = atom.residue
            if residue is not None and residue.is_protein and atom.name in _BACKBONE_ATOM_NAMES:
                filtered.append(atom_index)
        return filtered

    def on_btn_run_analysis(self, widget, event=None):
        """ Function doc """
        vobject_id = self.coordinates_combobox.get_vobject_id()
        if vobject_id is None:
            self.main.simple_dialog.info(msg='Please select a System and an Object (trajectory) first.')
            return False
        vismol_object = self.main.vm_session.vm_objects_dic.get(vobject_id)
        if vismol_object is None:
            self.main.simple_dialog.info(msg='Please select a valid Object (trajectory).')
            return False

        if self.builder.get_object('chk_limit_to_selection').get_active():
            selected_list, residue_dict, selected_vobject = self.vm_session.build_index_list_from_atom_selection(return_vobject=True)
            if selected_vobject is not vismol_object or not selected_list:
                self.main.simple_dialog.info(
                    msg='"Limit to current viewer selection" is checked, but there is no valid '
                        'selection on the chosen Object. Select atoms in the 3D view first, or '
                        'uncheck it to analyse the whole Object.')
                return False
            selection_list = selected_list
        else:
            # . Default: the whole Object -- vismol_object.atoms is keyed
            #   0-based, same indexing as vismol_object.frames' second
            #   axis (see the comment on the frame slice below).
            selection_list = list(vismol_object.atoms.keys())

        if self.builder.get_object('chk_backbone_only').get_active():
            selection_list = self._filter_backbone(selection_list, vismol_object)
            if not selection_list:
                self.main.simple_dialog.info(
                    msg='No backbone atoms (N, CA, C, O of a protein residue) found in the current selection.')
                return False

        f_init = 0
        f_last = -1
        step_size = 1
        if self.builder.get_object('chk_frames_from').get_active():
            f_init = int(self.builder.get_object('entry_frame_init').get_text())
            f_last = int(self.builder.get_object('entry_frame_last').get_text())
        if self.builder.get_object('chk_step_size').get_active():
            step_size = int(self.builder.get_object('entry_step_size').get_text())

        if f_last == -1:
            f_last = vismol_object.frames.shape[0]

        # . selection_list indexes are 0-based, same ordering as
        #   vismol_object.frames' second axis (confirmed by
        #   RMSD_analysis_window.py's own identical use of it).
        frames = vismol_object.frames[f_init:f_last:step_size][:, selection_list, :]
        if frames.shape[0] < 2:
            self.main.simple_dialog.info(msg='RMSF needs at least 2 frames in the selected range.')
            return False

        rmsf = md_analysis.compute_rmsf(frames)

        group_by_residue = self.builder.get_object('radio_per_residue').get_active()

        if group_by_residue:
            # . Grouped by the residue OBJECT itself (not by name/number
            #   alone), so residues sharing a number across different
            #   chains are never merged together.
            residue_order  = []
            residue_values = {}
            residue_labels = {}
            for position, atom_index in enumerate(selection_list):
                residue = vismol_object.atoms[atom_index].residue
                key = id(residue)
                if key not in residue_values:
                    residue_values[key] = []
                    residue_labels[key] = (residue.name, residue.index)
                    residue_order.append(key)
                residue_values[key].append(rmsf[position])

            x_values = [residue_labels[key][1] for key in residue_order]
            y_values = [sum(residue_values[key]) / len(residue_values[key]) for key in residue_order]

            text = 'residue\tn_atoms\tRMSF\n'
            for key in residue_order:
                resname, resindex = residue_labels[key]
                values = residue_values[key]
                text += '{}{}\t{}\t{:.4f}\n'.format(resname, resindex, len(values), sum(values) / len(values))
            plot_title = 'RMSF per residue'
        else:
            x_values = list(selection_list)
            y_values = list(rmsf)

            text = 'index\tsymbol\tRMSF\n'
            for position, atom_index in enumerate(selection_list):
                atom = vismol_object.atoms[atom_index]
                text += '{}\t{}\t{:.4f}\n'.format(atom_index, atom.symbol, rmsf[position])
            plot_title = 'RMSF per atom'

        textwindow = TextWindow(text)

        if self.builder.get_object('chk_save_file').get_active():
            default_name = 'rmsf_per_residue.txt' if group_by_residue else 'rmsf_per_atom.txt'
            self._save_text_to_file(text, default_name)

        if self.builder.get_object('chk_plot').get_active():
            from util.easyplot import XYPlot
            self.plot = XYPlot()
            self.plot.add(X=x_values, Y=y_values,
                           symbol=None, sym_color=[1, 1, 1], sym_fill=False,
                           line='solid', line_color=[0, 0, 0], energy_label=None)
            self.plot.Ymin_list = [0]

            window = Gtk.Window()
            window.set_title(plot_title)
            window.add(self.plot)
            window.set_default_size(700, 420)
            window.show_all()
