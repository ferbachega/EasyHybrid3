#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Search Transition State (Baker) window
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
#      "Search Transition State" window -- Baker's eigenvector-following
#      saddle-point search (pSimulation.GeometryOptimization.
#      BakerSaddleOptimize_SystemGeometry, backed by pScientific.
#      ObjectiveFunctionIterators.BakerOptimizer). Modeled directly on
#      geometry_optimization_window.py: same CoordinatesComboBox +
#      SaveTrajectoryBox + p_session.run_simulation() dispatch, so this
#      job gets the exact same Process Manager / job-history / result-
#      import treatment every other "Geometry_Optimization" job already
#      gets -- the backend (p_methods/geometry_optimization.py) was
#      extended with a "BakerSaddle" optimizer branch rather than
#      inventing a new simulation_type, specifically so nothing about
#      that existing machinery needed to change.
#
import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

import os

from gui.widgets.custom_widgets import SaveTrajectoryBox
from gui.widgets.custom_widgets import CoordinatesComboBox

HOME = os.environ.get('HOME')


class TransitionStateSearchWindow:
    """ "Search Transition State" (Baker) window. """

    def __init__(self, main=None):
        """ Class initialiser """
        self.main      = main
        self.home      = main.home
        self.p_session = main.p_session
        self.Visible   = False
        self.sym_tag   = 'ts_search'

    def open_window(self):
        """ Function doc """
        if self.Visible:
            self.window.present()
            return

        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/simulation/transition_state_search_window.glade'))
        self.builder.connect_signals(self)

        self.window = self.builder.get_object('transition_state_search_window')
        self.window.set_title('Search Transition State (Baker)')
        self.window.set_keep_above(True)

        # -------------------- starting coordinates --------------------
        self.box_coordinates = self.builder.get_object('box_coordinates')
        self.combobox_starting_coordinates = CoordinatesComboBox()
        self.combobox_starting_coordinates.connect("changed", self.on_name_combo_changed)
        self.box_coordinates.pack_start(self.combobox_starting_coordinates, False, False, 0)
        self._starting_coordinates_model_update(init=True)

        # -------------------- Baker-specific widgets --------------------
        self.spinbtn_follow_mode        = self.builder.get_object('spinbtn_follow_mode')
        self.combobox_hessian_updating  = self.builder.get_object('combobox_hessian_updating')
        self.spinbtn_analytic_frequency = self.builder.get_object('spinbtn_analytic_frequency')

        # -------------------- trajectory saving --------------------
        self.save_trajectory_box = SaveTrajectoryBox(parent=self.main, home=self.home)
        self.builder.get_object('ts_search_parm_box').pack_end(self.save_trajectory_box.box, True, True, 0)

        # -------------------- working folder default --------------------
        if self.main.p_session.psystem[self.main.p_session.active_id]:
            folder = self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder
            if folder is None:
                folder = HOME
            self.update_working_folder_chooser(folder=folder)

            output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
            self.save_trajectory_box.set_filename(output_name)

        self.window.show_all()
        self.Visible = True

    def close_window(self, button=None, data=None):
        """ Function doc """
        self.window.destroy()
        self.Visible = False

    def on_name_combo_changed(self, widget):
        """ Function doc """
        pass

    def run_dialog(self):
        """ Function doc """
        dialog = Gtk.MessageDialog(
            transient_for=self.main.window,
            message_type=Gtk.MessageType.ERROR,
            buttons=Gtk.ButtonsType.OK,
            text="Folder not found",
        )
        dialog.format_secondary_text(
            "The folder you have selected does not appear to be valid. Please select a different folder or create a new one."
        )
        dialog.run()
        dialog.destroy()

    def run_opt(self, button):
        """ Runs a Baker saddle-point search based on the parameters set in the GUI. """
        simParameters = {
            "simulation_type"     : "Geometry_Optimization",
            "optimizer"           : "BakerSaddle",
            "trajectory_name"     : None,
            "dialog"              : False,
            "folder"              : os.getcwd(),
            "maximumIterations"   : 50,
            "logFrequency"        : 1,
            "save_frequency"      : 10,
            "rmsGradientTolerance": 0.1,
            "followMode"          : 0,
            "analyticFrequency"   : 0,
            "hessianUpdatingOption": "BOFILL",
            "save_format"         : None,
            "save_traj"           : False,
            "save_pdb"            : False,
        }

        # . Starting coordinates.
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            model = self.combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]
            vobject = self.main.vm_session.vm_objects_dic[vobject_id]
            self.main.p_session.set_psystem_coordinates_from_vobject(vobject)
        simParameters['obj1_key6'] = vobject.key6

        # . Common parameters.
        simParameters["logFrequency"]         = int(self.builder.get_object('entry_log_frequence').get_text())
        simParameters["maximumIterations"]    = int(self.builder.get_object('entry_max_int').get_text())
        simParameters["rmsGradientTolerance"] = float(self.builder.get_object('entry_rmsd_tol').get_text())
        simParameters["vobject_name"]         = self.save_trajectory_box.builder.get_object('entry_trajectory_name').get_text()

        # . Baker-specific parameters.
        simParameters["followMode"]            = int(self.spinbtn_follow_mode.get_value())
        simParameters["analyticFrequency"]     = int(self.spinbtn_analytic_frequency.get_value())
        # . get_active_id(), not get_active_text() -- the combobox's
        #   display text is the descriptive "MS (Murtagh-Sargent)" label,
        #   while the id (set via the .glade <item id="..."> attribute)
        #   is the short code SymmetricMatrix.Update()'s own
        #   `option` argument actually expects ("BFGS"/"BOFILL"/"MS"/
        #   "POWELL"); passing the display text through would silently
        #   fail the backend's option.upper() lookup and fall back to
        #   BFGS regardless of what was actually selected.
        simParameters["hessianUpdatingOption"] = self.combobox_hessian_updating.get_active_id()

        # . Trajectory saving (if enabled).
        if self.save_trajectory_box.builder.get_object('checkbox_save_traj').get_active():
            simParameters["save_traj"]       = True
            simParameters["dialog"]          = True
            simParameters["folder"]          = self.save_trajectory_box.folder_chooser_button.get_folder()
            simParameters["trajectory_name"] = self.save_trajectory_box.builder.get_object('entry_trajectory_name').get_text()
            simParameters["save_frequency"]  = int(self.save_trajectory_box.builder.get_object('entry_trajectory_frequency').get_text())
            simParameters["trajectory_name"] += ".ptGeo"

            saveFormat = self.save_trajectory_box.builder.get_object('combobox_format').get_active()
            if   saveFormat == 0: simParameters["save_format"] = ".ptGeo"
            elif saveFormat == 1: simParameters["save_format"] = ".mdcrd"
            elif saveFormat == 2: simParameters["save_format"] = ".dcd"
            elif saveFormat == 3: simParameters["save_format"] = ".xyz"

            self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder = simParameters["folder"]

        if not os.path.exists(simParameters['folder']):
            self.run_dialog()
            return None

        self.main.p_session.run_simulation(parameters=simParameters)

        self.window.destroy()
        self.Visible = False

    def _starting_coordinates_model_update(self, init=False):
        """ Function doc """
        if self.Visible or init:
            e_id = self.main.p_session.active_id
            self.combobox_starting_coordinates.set_model(self.main.vobject_liststore_dict[e_id])
            size = len(self.main.vobject_liststore_dict[e_id])
            self.combobox_starting_coordinates.set_active(size - 1)

    def update(self, parameters=None):
        """ Function doc """
        if self.Visible:
            self._starting_coordinates_model_update()
            self.update_working_folder_chooser()
            if self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
                self.save_trajectory_box.set_filename(output_name)

    def update_working_folder_chooser(self, folder=None):
        """ Function doc """
        if folder:
            self.save_trajectory_box.set_folder(folder=folder)
        else:
            folder = self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder
            if folder:
                self.save_trajectory_box.set_folder(folder=folder)

    def restore_the_parameters_to_the_window(self, parameters):
        """ Restore exported parameters into the GTK window (Process Manager "reopen job"). """
        self.combobox_starting_coordinates.set_active(parameters['cb1_active'])

        self.builder.get_object('entry_log_frequence').set_text(str(parameters['logFrequency']))
        self.builder.get_object('entry_max_int').set_text(str(parameters['maximumIterations']))
        self.builder.get_object('entry_rmsd_tol').set_text(str(parameters['rmsGradientTolerance']))

        self.spinbtn_follow_mode.set_value(parameters.get('followMode', 0))
        self.spinbtn_analytic_frequency.set_value(parameters.get('analyticFrequency', 0))

        # . set_active_id(), matching the id-based read in run_opt() --
        #   see that method's own comment for why the text/id distinction
        #   matters here.
        self.combobox_hessian_updating.set_active_id(parameters.get('hessianUpdatingOption', 'BOFILL'))

        if parameters.get('save_traj'):
            checkbox = self.save_trajectory_box.builder.get_object('checkbox_save_traj')
            checkbox.set_active(True)
            self.save_trajectory_box.folder_chooser_button.set_folder(parameters.get('folder', ""))
            trajectory_name = parameters.get('trajectory_name')
            if trajectory_name:
                self.save_trajectory_box.builder.get_object('entry_trajectory_name').set_text(trajectory_name)
            entry_freq = self.save_trajectory_box.builder.get_object('entry_trajectory_frequency')
            entry_freq.set_text(str(parameters.get('save_frequency', "")))
            self.save_trajectory_box.builder.get_object('combobox_format').set_active(0)
