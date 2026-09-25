#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Ramachandran (phi/psi) analysis window
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
#      Ramachandran plot window -- picks a System/Object (same
#      SystemComboBox/CoordinatesComboBox pair used throughout the rest
#      of the app), computes backbone phi/psi for every protein residue
#      via util/ramachandran.py (GTK-free), and renders it with
#      util/rama_plot.py's RamachandranPlot widget, embedded directly in
#      this window (box_plot) rather than a separate popup, so the
#      "residue under the mouse" hover feedback has somewhere to land
#      (label_status).
#
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

import os

from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.widgets.custom_widgets import FolderChooserButton
from gui.windows.setup.windows_and_dialogs import TextWindow
from util import ramachandran
from util.rama_plot import RamachandranPlot, CATEGORY_LABELS


class RamachandranAnalysisWindow:

    def __init__(self, main=None):
        """ Class initialiser """
        self.main       = main
        self.home       = main.home
        self.p_session  = self.main.p_session
        self.vm_session = self.main.vm_session
        self.Visible    = False
        self._last_results = []

    def open_window(self):
        """ Function doc """
        if self.Visible:
            self.window.present()
            return

        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/analysis/ramachandran_analysis_window.glade'))
        self.builder.connect_signals(self)

        self.window = self.builder.get_object('window')
        self.window.set_title('Ramachandran Plot')

        self.label_status = self.builder.get_object('label_status')
        self._default_status_text = self.label_status.get_text()

        self.builder.get_object('btn_close').connect('clicked', self.close_window)

        self.box_system = self.builder.get_object('box_system')
        self.combobox_systems = SystemComboBox(self.main)
        self.combobox_systems.connect("changed", self.on_combobox_systems_changed)
        self.box_system.pack_start(self.combobox_systems, True, True, 0)

        self.box_coordinates = self.builder.get_object('box_coordinates')
        self.coordinates_combobox = CoordinatesComboBox()
        self.box_coordinates.pack_start(self.coordinates_combobox, True, True, 0)

        self.grid_save_output = self.builder.get_object('grid_save_output')
        self.entry_output_filename = self.builder.get_object('entry_output_filename')

        self.box_output_folder = self.builder.get_object('box_output_folder')
        self.output_folder_chooser = FolderChooserButton(self.main, 'folder', self.home)
        self.box_output_folder.pack_start(self.output_folder_chooser.btn, False, False, 0)

        self.plot = RamachandranPlot()
        self.plot.on_hover_callback = self._on_plot_hover
        self.box_plot = self.builder.get_object('box_plot')
        self.box_plot.pack_start(self.plot, True, True, 0)

        if self.p_session.psystem:
            self.combobox_systems.set_active_system(e_id=self.p_session.active_id)

        self.window.show_all()
        self.Visible = True

    def on_save_file_checkbox(self, widget):
        self.grid_save_output.set_sensitive(widget.get_active())

    def on_all_frames_checkbox(self, widget):
        self.builder.get_object('frame_range').set_sensitive(widget.get_active())

    def on_combobox_systems_changed(self, widget):
        system_id = self.combobox_systems.get_system_id()
        if system_id is None:
            return
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

    def _on_plot_hover(self, point):
        if point is None:
            self.label_status.set_text(self._default_status_text)
            return
        frame_text = '' if point['frame'] is None else ' (frame {})'.format(point['frame'])
        self.label_status.set_text(
            'Resíduo {}{} — cadeia {} — φ = {:.1f}°, ψ = {:.1f}°{} [{}]'.format(
                point['resname'], point['resnum'], point['chain'] or '-',
                point['phi'], point['psi'], frame_text,
                CATEGORY_LABELS.get(point['category'], point['category'])))

    def _save_text_to_file(self, text, default_filename):
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

    def _residues_with_selected_atoms(self, vismol_object):
        """ Returns a set of (chain_name, resnum) for every residue that
            has at least one currently-selected atom, restricted to
            `vismol_object` -- used to filter (not re-compute) results
            when "Limit to current viewer selection" is checked: phi/psi
            still needs each residue's full backbone context (its own
            and its neighbours' N/CA/C), so the selection can only
            narrow down which already-computed points get PLOTTED.
        """
        selected_list, _residue_dict, selected_vobject = \
            self.vm_session.build_index_list_from_atom_selection(return_vobject=True)
        if selected_vobject is not vismol_object or not selected_list:
            return None
        keys = set()
        for atom_index in selected_list:
            residue = vismol_object.atoms[atom_index].residue
            if residue is not None:
                keys.add((residue.chain.name if residue.chain is not None else None, residue.index))
        return keys

    def on_btn_run_analysis(self, widget, event=None):
        """ Function doc """
        vobject_id = self.coordinates_combobox.get_vobject_id()
        if vobject_id is None:
            self.main.simple_dialog.info(msg='Please select a System and an Object first.')
            return False
        vismol_object = self.main.vm_session.vm_objects_dic.get(vobject_id)
        if vismol_object is None:
            self.main.simple_dialog.info(msg='Please select a valid Object.')
            return False

        if self.builder.get_object('chk_all_frames').get_active():
            f_init = int(self.builder.get_object('entry_frame_init').get_text())
            f_last = int(self.builder.get_object('entry_frame_last').get_text())
            if f_last == -1:
                f_last = vismol_object.frames.shape[0]
            step_size = 1
            if self.builder.get_object('chk_step_size').get_active():
                step_size = int(self.builder.get_object('entry_step_size').get_text())
            frame_indices = list(range(f_init, f_last, step_size))
            if not frame_indices:
                self.main.simple_dialog.info(msg='The selected frame range is empty.')
                return False
        else:
            frame_indices = [self.vm_session.frame]

        results = ramachandran.compute_phi_psi(vismol_object, frame_indices=frame_indices)

        if self.builder.get_object('chk_limit_to_selection').get_active():
            allowed_keys = self._residues_with_selected_atoms(vismol_object)
            if not allowed_keys:
                self.main.simple_dialog.info(
                    msg='"Limit to current viewer selection" is checked, but there is no valid '
                        'selection on the chosen Object. Select atoms in the 3D view first, or '
                        'uncheck it to plot every residue.')
                return False
            results = [r for r in results if (r['chain'], r['resnum']) in allowed_keys]

        if not results:
            self.main.simple_dialog.info(
                msg='No phi/psi angles could be computed for this Object. This tool needs '
                    'protein residues (backbone N, CA, C atoms present) with an unbroken '
                    'peptide bond to at least one neighbour -- a purely non-protein system, '
                    'or a selection limited to isolated/terminal residues, produces no points.')
            return False

        self._last_results = results
        self.plot.set_data(results)

        n_residues = len(set((r['chain'], r['resnum']) for r in results))
        n_frames = len(set(r['frame'] for r in results))
        self.label_status.set_text(
            '{} pontos plotados ({} resíduos, {} frame{}). Passe o mouse sobre um ponto '
            'para identificar o resíduo.'.format(
                len(results), n_residues, n_frames, 's' if n_frames != 1 else ''))

        if self.builder.get_object('chk_save_file').get_active():
            text = self._results_to_text(results)
            self._save_text_to_file(text, 'ramachandran.txt')

        return False

    def _results_to_text(self, results):
        text = 'chain\tresidue\tresnum\tcategory\tframe\tphi\tpsi\n'
        for r in results:
            text += '{}\t{}\t{}\t{}\t{}\t{:.4f}\t{:.4f}\n'.format(
                r['chain'], r['resname'], r['resnum'], r['category'],
                '' if r['frame'] is None else r['frame'], r['phi'], r['psi'])
        return text

    def on_btn_export_png_clicked(self, widget):
        """ Same "Save As" pattern as PES_analysis_window.py's own PNG
            export: redraws the plot via Cairo at 4x the on-screen
            resolution (see util/easyplot/export_utils.py), it is not a
            screenshot of the widget.
        """
        if not self._last_results:
            self.main.simple_dialog.info(msg='Run the analysis first.')
            return

        from util.easyplot import export_plot_to_png

        dialog = Gtk.FileChooserDialog(
            title="Export figure as PNG",
            transient_for=self.window,
            action=Gtk.FileChooserAction.SAVE,
        )
        dialog.add_buttons(
            Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
            Gtk.STOCK_SAVE, Gtk.ResponseType.OK,
        )
        dialog.set_do_overwrite_confirmation(True)
        dialog.set_current_name('ramachandran.png')

        filter_png = Gtk.FileFilter()
        filter_png.set_name("PNG image (*.png)")
        filter_png.add_pattern("*.png")
        dialog.add_filter(filter_png)

        response = dialog.run()
        filepath = None
        if response == Gtk.ResponseType.OK:
            filepath = dialog.get_filename()
            if not filepath.lower().endswith('.png'):
                filepath += '.png'
        dialog.destroy()

        if filepath is None:
            return

        try:
            export_plot_to_png(self.plot, filepath, scale=4)
        except Exception as e:
            self.main.simple_dialog.info(msg='Could not export the plot:\n{}'.format(e))
