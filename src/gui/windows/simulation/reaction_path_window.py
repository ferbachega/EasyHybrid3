#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: IRC / Reaction Path window
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
#      "IRC / Reaction Path" window -- steepest-descent reaction path
#      starting from a first-order saddle point (pSimulation.
#      SteepestDescentReactionPath.SteepestDescentPath_SystemGeometry,
#      fromSaddle=True always here). Modeled directly on
#      transition_state_search_window.py: same CoordinatesComboBox +
#      SaveTrajectoryBox + p_session.run_simulation() dispatch, so this
#      job gets the exact same Process Manager / job-history / result-
#      import treatment every other "Geometry_Optimization" job already
#      gets -- the backend (p_methods/geometry_optimization.py) was
#      extended with a "ReactionPath" optimizer branch rather than
#      inventing a new simulation_type.
#
import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

import os

from gui.widgets.custom_widgets import SaveTrajectoryBox
from gui.widgets.custom_widgets import CoordinatesComboBox

HOME = os.environ.get('HOME')


class ReactionPathWindow:
    """ "IRC / Reaction Path" window. """

    def __init__(self, main=None):
        """ Class initialiser """
        self.main      = main
        self.home      = main.home
        self.p_session = main.p_session
        self.Visible   = False
        self.sym_tag   = 'reaction_path'

    def open_window(self):
        """ Function doc """
        if self.Visible:
            self.window.present()
            return

        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/simulation/reaction_path_window.glade'))
        self.builder.connect_signals(self)

        self.window = self.builder.get_object('reaction_path_window')
        self.window.set_title('IRC / Reaction Path')
        self.window.set_keep_above(True)

        # -------------------- starting coordinates --------------------
        self.box_coordinates = self.builder.get_object('box_coordinates')
        self.combobox_starting_coordinates = CoordinatesComboBox()
        self.combobox_starting_coordinates.connect("changed", self.on_name_combo_changed)
        self.box_coordinates.pack_start(self.combobox_starting_coordinates, False, False, 0)
        self._starting_coordinates_model_update(init=True)

        # -------------------- reaction-path-specific widgets --------------------
        self.entry_function_step     = self.builder.get_object('entry_function_step')
        self.entry_path_step         = self.builder.get_object('entry_path_step')
        self.checkbox_mass_weighting = self.builder.get_object('checkbox_mass_weighting')

        # -------------------- trajectory saving --------------------
        # [EN] Mandatory here, unlike every other Geometry_Optimization
        # window: SteepestDescentPath_SystemGeometry takes a FIXED number
        # of FIXED-length steps per branch with no convergence check at
        # all (see ObjectiveFunctionIterator.Continue() in pDynamo3 --
        # pure iteration count, nothing else) and its per-step direction
        # is the gradient normalized to a unit vector before being scaled
        # -- once a branch is near a stationary point the "true" gradient
        # is ~0 and this normalization amplifies whatever numerical noise
        # is left, so the walk continues as an effectively random kick of
        # size pathStep every remaining iteration. The FINAL system
        # coordinates (what every other optimizer in this app presents as
        # "the result") are therefore NOT a chemically meaningful
        # end-point for this tool -- only the saved trajectory, inspected
        # frame by frame (with the RMS Gradient column already logged by
        # SteepestDescentPathFinder.LogIteration), is. Locking the
        # checkbox on prevents the single most misleading way to use this
        # window.
        self.save_trajectory_box = SaveTrajectoryBox(parent=self.main, home=self.home)
        self.builder.get_object('rp_parm_box').pack_end(self.save_trajectory_box.box, True, True, 0)
        self.save_trajectory_box.set_active(True)
        self.save_trajectory_box.builder.get_object('checkbox_save_traj').set_sensitive(False)

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
        """ Runs a steepest-descent reaction path based on the parameters set in the GUI. """
        simParameters = {
            "simulation_type"     : "Geometry_Optimization",
            "optimizer"           : "ReactionPath",
            "trajectory_name"     : None,
            "dialog"              : False,
            "folder"              : os.getcwd(),
            "maximumIterations"   : 50,
            "logFrequency"        : 1,
            "save_frequency"      : 1,
            "fromSaddle"          : True,
            "functionStep"        : 0.2,
            "pathStep"            : 0.025,
            "useMassWeighting"    : False,
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
        simParameters["logFrequency"]      = int(self.builder.get_object('entry_log_frequence').get_text())
        simParameters["maximumIterations"] = int(self.builder.get_object('entry_max_int').get_text())
        simParameters["vobject_name"]      = self.save_trajectory_box.builder.get_object('entry_trajectory_name').get_text()

        # . Reaction-path-specific parameters.
        simParameters["functionStep"]     = float(self.entry_function_step.get_text())
        simParameters["pathStep"]         = float(self.entry_path_step.get_text())
        simParameters["useMassWeighting"] = self.checkbox_mass_weighting.get_active()

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

        self.entry_function_step.set_text(str(parameters.get('functionStep', 0.2)))
        self.entry_path_step.set_text(str(parameters.get('pathStep', 0.025)))
        self.checkbox_mass_weighting.set_active(parameters.get('useMassWeighting', False))

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
