#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: Radial distribution function g(r) analysis window
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
#      "RDF Analysis" window -- the pair radial distribution function g(r)
#      (plus running coordination number) between two named selections
#      over a trajectory's periodic cell. Reuses the exact same
#      named-selection mechanism as the SMD window's pulled/fixed-group
#      pickers (system.e_selections, populated via Selection -> Send to
#      Selection List) instead of building a new atom picker, and the
#      same per-frame orthorhombic cell (vismol_object.cell_coordinates)
#      reimaging_trajectory.py already reads. The actual math lives in
#      util/md_analysis.py (GTK-free, reusable, NumPy-only).
#
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

import os

from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.windows.setup.windows_and_dialogs import TextWindow
from util import md_analysis

_ALL_ATOMS = "(all atoms)"


class RDFAnalysisWindow:

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
        self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/analysis/RDF_analysis.glade'))
        self.builder.connect_signals(self)

        self.window = self.builder.get_object('window')

        self.label_info        = self.builder.get_object('label_info')
        self.combo_group_a     = self.builder.get_object('combo_group_a')
        self.combo_group_b     = self.builder.get_object('combo_group_b')
        self.spinbtn_r_max     = self.builder.get_object('spinbtn_r_max')
        self.spinbtn_bin_width = self.builder.get_object('spinbtn_bin_width')

        self.builder.get_object('btn_close').connect('clicked', self.close_window)
        self.builder.get_object('chk_frames_from').connect('toggled', self.on_frame_range_checkbox)
        self.builder.get_object('chk_step_size').connect('toggled', self.on_step_size_checkbox)

        self.box_system = self.builder.get_object('box_system')
        self.combobox_systems = SystemComboBox(self.main)
        self.combobox_systems.connect("changed", self.on_combobox_systems_changed)
        self.box_system.pack_start(self.combobox_systems, False, False, 0)

        self.box_coordinates = self.builder.get_object('box_coordinates')
        self.coordinates_combobox = CoordinatesComboBox()
        self.coordinates_combobox.connect("changed", self.on_combobox_coordinates_changed)
        self.box_coordinates.pack_start(self.coordinates_combobox, False, False, 0)

        if self.p_session.psystem:
            self.combobox_systems.set_active_system(e_id=self.p_session.active_id)

        self.window.show_all()
        self.Visible = True

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

    def on_combobox_systems_changed(self, widget):
        system_id = self.combobox_systems.get_system_id()
        if system_id is None:
            return
        self.coordinates_combobox.set_model(self.main.vobject_liststore_dict[system_id])
        size = len(list(self.main.vobject_liststore_dict[system_id]))
        self.coordinates_combobox.set_active(size - 1)

    def on_combobox_coordinates_changed(self, widget):
        self._refresh_group_combos()

    def _refresh_group_combos(self):
        self.combo_group_a.remove_all()
        self.combo_group_b.remove_all()
        self.combo_group_a.append_text(_ALL_ATOMS)
        self.combo_group_b.append_text(_ALL_ATOMS)
        self.combo_group_a.set_active(0)
        self.combo_group_b.set_active(0)

        system_id = self.coordinates_combobox.get_system_id()
        if system_id is None:
            return
        system = self.p_session.psystem.get(system_id)
        if system is None:
            return
        for name in system.e_selections.keys():
            self.combo_group_a.append_text(name)
            self.combo_group_b.append_text(name)

    def _export_plot_png(self, parent, plot, default_name):
        dialog = Gtk.FileChooserDialog(title='Export g(r) plot', transient_for=parent, action=Gtk.FileChooserAction.SAVE)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        dialog.set_do_overwrite_confirmation(True)
        dialog.set_current_name(default_name)
        if dialog.run() == Gtk.ResponseType.OK:
            path = dialog.get_filename()
            if not path.lower().endswith('.png'):
                path += '.png'
            plot.export_png(path)
        dialog.destroy()

    def _resolve_group(self, combo, system_id, n_atoms):
        name = combo.get_active_text()
        if name is None or name == _ALL_ATOMS:
            return list(range(n_atoms))
        return list(self.p_session.psystem[system_id].e_selections[name])

    def on_btn_run_analysis(self, widget, event=None):
        """ Function doc """
        vobject_id = self.coordinates_combobox.get_vobject_id()
        system_id = self.coordinates_combobox.get_system_id()
        if vobject_id is None or system_id is None:
            self.main.simple_dialog.info(msg='Please select a System and an Object (trajectory) first.')
            return

        vismol_object = self.main.vm_session.vm_objects_dic.get(vobject_id)
        if vismol_object is None:
            self.main.simple_dialog.info(msg='Please select a valid Object (trajectory).')
            return

        n_atoms = len(vismol_object.atoms)
        group_a = self._resolve_group(self.combo_group_a, system_id, n_atoms)
        group_b = self._resolve_group(self.combo_group_b, system_id, n_atoms)
        if not group_a or not group_b:
            self.main.simple_dialog.info(msg='Both Group A and Group B must contain at least one atom.')
            return
        same_selection = (self.combo_group_a.get_active_text() == self.combo_group_b.get_active_text())

        if vismol_object.cell_coordinates is None:
            self.main.simple_dialog.info(
                msg='This Object has no periodic cell information -- g(r) needs a periodic '
                    'box to normalise by density (see util/md_analysis.py). '
                    'Import/build the trajectory with its cell, or set one via "Reimaging (Wrapping)" first.')
            return

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

        frame_indices = range(f_init, f_last, step_size)
        n_cell_frames = vismol_object.cell_coordinates.shape[0]

        coords_a = [vismol_object.frames[f][group_a] for f in frame_indices]
        coords_b = [vismol_object.frames[f][group_b] for f in frame_indices]
        # . cell_coordinates[f][-1] is the box's far corner (a, b, c) for an
        #   orthorhombic cell built from the origin -- same convention
        #   reimaging_trajectory.py already relies on. Falls back to the
        #   single stored cell (index 0) for a trajectory whose box was
        #   only ever set once (no per-frame XST/XSC data).
        boxes = [vismol_object.cell_coordinates[min(f, n_cell_frames - 1)][-1] for f in frame_indices]

        if len(coords_a) < 1:
            self.main.simple_dialog.info(msg='No frames in the selected range.')
            return

        r_max = self.spinbtn_r_max.get_value()
        bin_width = self.spinbtn_bin_width.get_value()

        r, g_r, coordination = md_analysis.compute_rdf(
            coords_a, coords_b, boxes, same_selection=same_selection,
            r_max=r_max, bin_width=bin_width)

        text = 'r (Å)\tg(r)\tN(r)\n'
        for i in range(len(r)):
            text += '{:.4f}\t{:.6f}\t{:.4f}\n'.format(r[i], g_r[i], coordination[i])
        textwindow = TextWindow(text)

        if self.builder.get_object('chk_plot').get_active():
            # . 2026-10-21: own g(r) plot (util/rdf_plot.py): g(r) from 0, N(r) on
            #   the right axis, axis titles, hover readout, PNG export
            from util.rdf_plot import RDFPlot
            name_a = self.combo_group_a.get_active_text() or _ALL_ATOMS
            name_b = self.combo_group_b.get_active_text() or _ALL_ATOMS
            self.plot = RDFPlot()
            self.plot.set_data(r, g_r, coordination, title='{} – {}  ({} frames)'.format(name_a, name_b, len(coords_a)))

            window = Gtk.Window()
            window.set_title('g(r)')
            box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=4, margin=6)
            box.pack_start(self.plot, True, True, 0)
            export = Gtk.Button(label='Export PNG…')
            export.set_halign(Gtk.Align.END)
            plot = self.plot
            default_name = 'rdf_{}_{}.png'.format(name_a, name_b).replace(' ', '_').replace('(', '').replace(')', '')
            export.connect('clicked', lambda b: self._export_plot_png(window, plot, default_name))
            box.pack_start(export, False, False, 0)
            window.add(box)
            window.set_default_size(720, 460)
            window.show_all()

        self.label_info.set_text(
            'Done -- {} frames, {} atoms in A, {} in B.'.format(len(coords_a), len(group_a), len(group_b)))
