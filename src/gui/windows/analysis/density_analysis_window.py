#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Box density analysis window
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
#      "Box Density Analysis" window -- mass density (g/cm^3) of the
#      periodic box along a trajectory, over a user-defined frame range
#      and step. Total mass comes from the pDynamo system's atoms; the
#      volume of each frame comes from the per-frame cell stored in
#      vismol_object.cell_coordinates (same source RDF_analysis_window.py
#      and reimaging_trajectory.py use). Writes a log (per-frame volume
#      and density, plus mean/std/min/max at the end) and plots the
#      density with EasyPlot. The math lives in util/md_analysis.py.
#
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

import os
import time
import numpy as np

from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from util import md_analysis


class DensityAnalysisWindow:

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
        self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/analysis/density_analysis.glade'))
        self.builder.connect_signals(self)

        self.window = self.builder.get_object('window')

        self.label_info     = self.builder.get_object('label_info')
        self.entry_log_file = self.builder.get_object('entry_log_file')

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
        """ Suggests a log file name in the system's working folder,
            named after the selected Object. """
        system_id  = self.coordinates_combobox.get_system_id()
        vobject_id = self.coordinates_combobox.get_vobject_id()
        if system_id is None or vobject_id is None:
            return
        system = self.p_session.psystem.get(system_id)
        vismol_object = self.vm_session.vm_objects_dic.get(vobject_id)
        if system is None or vismol_object is None:
            return

        folder = getattr(system, 'e_working_folder', None) or os.getcwd()
        name = ''.join(c if (c.isalnum() or c in '-_.') else '_' for c in str(vismol_object.name))
        self.entry_log_file.set_text(os.path.join(folder, 'density_{}.log'.format(name)))

    def on_btn_browse_log(self, widget):
        dialog = Gtk.FileChooserDialog(title='Save density log as...',
                                       parent=self.window,
                                       action=Gtk.FileChooserAction.SAVE)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                           Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        dialog.set_do_overwrite_confirmation(True)

        current = self.entry_log_file.get_text().strip()
        if current:
            folder = os.path.dirname(current)
            if os.path.isdir(folder):
                dialog.set_current_folder(folder)
            dialog.set_current_name(os.path.basename(current))
        else:
            dialog.set_current_name('density.log')

        if dialog.run() == Gtk.ResponseType.OK:
            self.entry_log_file.set_text(dialog.get_filename())
        dialog.destroy()

    def _get_frame_range(self, n_frames):
        """ Returns the list of frame indices chosen by the user (last
            frame inclusive, -1 = last frame of the trajectory), or raises
            ValueError with a message for the user. """
        f_init = 0
        f_last = n_frames - 1
        step_size = 1
        try:
            if self.builder.get_object('chk_frames_from').get_active():
                f_init = int(self.builder.get_object('entry_frame_init').get_text())
                f_last = int(self.builder.get_object('entry_frame_last').get_text())
            if self.builder.get_object('chk_step_size').get_active():
                step_size = int(self.builder.get_object('entry_step_size').get_text())
        except ValueError:
            raise ValueError('Frame range and step size must be integer numbers.')

        if f_last == -1:
            f_last = n_frames - 1
        if step_size < 1:
            raise ValueError('Step size must be at least 1.')
        if not (0 <= f_init < n_frames) or not (0 <= f_last < n_frames):
            raise ValueError('Frame range must be within 0 and {} (this Object has {} frames).'.format(
                n_frames - 1, n_frames))
        if f_init > f_last:
            raise ValueError('The first frame ({}) is after the last frame ({}).'.format(f_init, f_last))

        return list(range(f_init, f_last + 1, step_size))

    def on_btn_run_analysis(self, widget, event=None):
        """ Function doc """
        vobject_id = self.coordinates_combobox.get_vobject_id()
        system_id  = self.coordinates_combobox.get_system_id()
        if vobject_id is None or system_id is None:
            self.main.simple_dialog.info(msg='Please select a System and an Object (trajectory) first.')
            return

        vismol_object = self.vm_session.vm_objects_dic.get(vobject_id)
        system = self.p_session.psystem.get(system_id)
        if vismol_object is None or system is None:
            self.main.simple_dialog.info(msg='Please select a valid Object (trajectory).')
            return

        if vismol_object.cell_coordinates is None:
            self.main.simple_dialog.info(
                msg='This Object has no periodic cell information, so the box volume is unknown. '
                    'Import the trajectory with its cell (e.g. a DCD from an NPT run of a periodic system), '
                    'or set one via "Reimaging (Wrapping)" first.')
            return

        n_atoms = len(system.atoms)
        if n_atoms != len(vismol_object.atoms):
            self.main.simple_dialog.info(
                msg='The Object has {} atoms but the System has {}; cannot compute the box mass.'.format(
                    len(vismol_object.atoms), n_atoms))
            return

        n_frames = vismol_object.frames.shape[0]
        try:
            frame_indices = self._get_frame_range(n_frames)
        except ValueError as error:
            self.main.simple_dialog.info(msg=str(error))
            return

        log_file = self.entry_log_file.get_text().strip()
        if not log_file:
            self.main.simple_dialog.info(msg='Please choose a log file.')
            return

        # . Same fallback as RDF_analysis_window.py: a trajectory whose
        #   box was only stored once (or fewer times than it has frames)
        #   reuses the last stored cell -- in that case the density is
        #   (partly) constant, and the log says so.
        n_cell_frames = vismol_object.cell_coordinates.shape[0]
        cells = [vismol_object.cell_coordinates[min(f, n_cell_frames - 1)] for f in frame_indices]

        # . pDynamo's Atom.mass is the element's standard atomic weight
        #   (PeriodicTable, looked up by atomicNumber), not a topology
        #   mass; an atom with an undefined element raises AttributeError.
        try:
            total_mass = float(sum(atom.mass for atom in system.atoms))
        except AttributeError:
            self.main.simple_dialog.info(
                msg='At least one atom of this System has no defined element, so its mass is unknown.')
            return
        try:
            volumes, densities = md_analysis.compute_box_density(total_mass, cells)
        except ValueError as error:
            self.main.simple_dialog.info(msg=str(error))
            return

        mean = float(np.mean(densities))
        std  = float(np.std(densities))

        warning = None
        if n_cell_frames != n_frames:
            warning = ('The Object has {} frames but only {} stored cell(s); frames beyond the last '
                       'stored cell reuse it (constant volume).'.format(n_frames, n_cell_frames))

        try:
            self._write_log(log_file, system, vismol_object, total_mass, frame_indices,
                            volumes, densities, warning)
        except OSError as error:
            self.main.simple_dialog.info(msg='Could not write the log file:\n{}'.format(error))
            return

        if self.builder.get_object('chk_plot').get_active():
            self._plot(frame_indices, densities, mean, vismol_object.name)

        info = 'Done: {} frames. Mean density = {:.4f} ± {:.4f} g/cm³\nLog: {}'.format(
            len(frame_indices), mean, std, log_file)
        if warning:
            info += '\nWarning: ' + warning
        self.label_info.set_text(info)

    def _write_log(self, log_file, system, vismol_object, total_mass, frame_indices,
                   volumes, densities, warning=None):
        mean = np.mean(densities)
        lines = []
        lines.append('#' + '-' * 60)
        lines.append('#  EasyHybrid - Box Density Analysis')
        lines.append('#' + '-' * 60)
        lines.append('#  Date        : {}'.format(time.strftime('%Y-%m-%d %H:%M:%S')))
        lines.append('#  System      : {}'.format(getattr(system, 'label', '')))
        lines.append('#  Object      : {}'.format(vismol_object.name))
        lines.append('#  Atoms       : {}'.format(len(system.atoms)))
        lines.append('#  Total mass  : {:.4f} amu (g/mol)'.format(total_mass))
        lines.append('#  Frames      : {} to {} (step {}), {} frames analysed'.format(
            frame_indices[0], frame_indices[-1],
            frame_indices[1] - frame_indices[0] if len(frame_indices) > 1 else 1,
            len(frame_indices)))
        if warning:
            lines.append('#  WARNING     : {}'.format(warning))
        lines.append('#' + '-' * 60)
        lines.append('#{:>9s} {:>16s} {:>16s}'.format('frame', 'volume(A^3)', 'density(g/cm^3)'))
        for f, v, d in zip(frame_indices, volumes, densities):
            lines.append('{:>10d} {:>16.4f} {:>16.6f}'.format(f, v, d))
        lines.append('#' + '-' * 60)
        lines.append('#  Average volume  : {:.4f} A^3'.format(np.mean(volumes)))
        lines.append('#  Minimum density : {:.6f} g/cm^3'.format(np.min(densities)))
        lines.append('#  Maximum density : {:.6f} g/cm^3'.format(np.max(densities)))
        lines.append('#  Std. deviation  : {:.6f} g/cm^3'.format(np.std(densities)))
        lines.append('#  AVERAGE DENSITY : {:.6f} g/cm^3'.format(mean))
        lines.append('#' + '-' * 60)

        with open(log_file, 'w') as fh:
            fh.write('\n'.join(lines) + '\n')

    def _plot(self, frame_indices, densities, mean, name):
        from util.easyplot import XYPlot
        self.plot = XYPlot()
        X = list(frame_indices)
        Y = [float(d) for d in densities]
        self.plot.add(X=X, Y=Y,
                      symbol=None, sym_color=[1, 1, 1], sym_fill=False,
                      line='solid', line_color=[0, 0, 0], energy_label=None)
        # . Average as a horizontal red line across the analysed range.
        self.plot.add(X=[X[0], X[-1]], Y=[mean, mean],
                      symbol=None, sym_color=[1, 1, 1], sym_fill=False,
                      line='solid', line_color=[0.85, 0.1, 0.1], energy_label=None)

        # . XYPlot.define_xy_limits() floors/ceils the Y axis to integers
        #   with an extra -1/+1 margin (made for energies), which flattens
        #   density fluctuations of ~0.01 g/cm^3 into a straight line.
        #   Tighten the Y axis around the data instead (add() is the only
        #   caller of define_xy_limits(), so this survives redraws).
        y_lo, y_hi = min(Y), max(Y)
        margin = 0.1 * (y_hi - y_lo) if y_hi > y_lo else max(0.005, 0.01 * abs(y_hi))
        self.plot.Ymin = y_lo - margin
        self.plot.Ymax = y_hi + margin
        self.plot.deltaY = self.plot.Ymax - self.plot.Ymin

        window = Gtk.Window()
        window.set_title('Box density (g/cm³): {}  |  mean = {:.4f}'.format(name, mean))
        window.add(self.plot)
        window.set_default_size(700, 420)
        window.show_all()
