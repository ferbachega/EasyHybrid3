#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Selection utilities for pDynamo systems
#
#  Copyright 2022-2025 Fernando Bachega
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
#      Provides functions for selecting atoms and residues in pDynamo systems
#      to facilitate QM/MM partitioning and molecular simulations.
#
import gi
from pprint import pprint
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
#from GTKGUI.gtkWidgets.filechooser import FileChooser
from gui.widgets.custom_widgets import FolderChooserButton
from gui.widgets.custom_widgets import SaveTrajectoryBox
from gui.widgets.custom_widgets import CoordinatesComboBox

import gc
import os

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')


#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Selection utilities for pDynamo systems
#
#  Author: Fernando Bachega
#  Maintainer: <ferbachega@gmail.com> or <easyhybrid3@gmail.com>
#

import os
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

# Custom GUI widgets
from gui.widgets.custom_widgets import (
    FolderChooserButton,
    SaveTrajectoryBox,
    CoordinatesComboBox,
)

# Environment variables
VISMOL_HOME = os.environ.get("VISMOL_HOME")
HOME = os.environ.get("HOME")


class MolecularDynamicsWindow:
    """GTK-based window for configuring and launching Molecular Dynamics simulations."""

    def __init__(self, main=None) -> None:
        """Initialize the Molecular Dynamics setup window."""
        self.main = main
        self.p_session = main.p_session
        self.home = main.home

        # GUI state
        self.window: Gtk.Window | None = None
        self.Visible: bool = False

        # GUI data models
        self.residue_liststore = Gtk.ListStore(str, str, str)
        self.job_liststore = Gtk.ListStore(str, str, str)

        # Available integrators
        self.md_intergators_dict: dict[int, str] = {
            0: "Velocity Verlet Dynamics",
            1: "Leap Frog Dynamics",
            2: "Langevin Dynamics",
        }

        # Available temperature scaling options
        self.temp_scale_options: dict[int, str] = {
            0: "constant",
            1: "linear",
            2: "exponential",
        }

        self.sym_tag: str = "mol_dyn"
        self.joblist: list = []
        self.last_parameters: dict | None = None  # stores last run parameters

    def open_window(self) -> None:
        """Open the Molecular Dynamics setup window."""
        if self.Visible:
            self.window.present()
            return

        # Load GUI from Glade file
        self.builder = Gtk.Builder()
        self.builder.add_from_file(
            os.path.join(self.home, "src/gui/windows/simulation/molecular_dynamics_window.glade")
        )

        self.window = self.builder.get_object("setup_md_window")
        self.window.set_title("Molecular Dynamics Setup")
        self.window.set_keep_above(True)
        self.window.connect("destroy", self.close_window)

        # Hide unused widgets initially
        for widget_name in [
            "check_pressure_control",
            "entry_pressure",
            "entry_temp_coupling",
            "label_temp_coupling",
            "label_pressure_coupling",
            "entry_pressure_coupling",
            "label_collision_freq",
            "entry_collision_frequency",
        ]:
            self.builder.get_object(widget_name).hide()

        # Setup coordinates combobox
        self.box_coordinates = self.builder.get_object("box_coordinates")
        self.combobox_starting_coordinates = CoordinatesComboBox()
        self.box_coordinates.pack_start(self.combobox_starting_coordinates, False, False, 0)
        self._starting_coordinates_model_update(init=True)
        size = len(self.combobox_starting_coordinates.get_model())
        self.combobox_starting_coordinates.set_active(size - 1)

        # Setup integrator combobox
        self.md_methods_liststore = Gtk.ListStore(str)
        for method in self.md_intergators_dict.values():
            self.md_methods_liststore.append([method])
        self.methods_combo = self.builder.get_object("comobobox_md_integrator")
        self.methods_combo.set_model(self.md_methods_liststore)
        renderer_text = Gtk.CellRendererText()
        self.methods_combo.pack_start(renderer_text, True)
        self.methods_combo.add_attribute(renderer_text, "text", 0)
        self.methods_combo.connect("changed", self.combobox_change)

        # Setup temperature scaling combobox
        self.temp_scale_option_store = Gtk.ListStore(str)
        for option in self.temp_scale_options.values():
            self.temp_scale_option_store.append([option])
        self.temp_scale_options_combo = self.builder.get_object("combobox_temperature_scale_options")
        self.temp_scale_options_combo.set_model(self.temp_scale_option_store)
        renderer_text = Gtk.CellRendererText()
        self.temp_scale_options_combo.pack_start(renderer_text, True)
        self.temp_scale_options_combo.add_attribute(renderer_text, "text", 0)
        self.temp_scale_options_combo.connect("changed", self.combobox_change)
        self.temp_scale_options_combo.set_active(0)

        # Setup trajectory saving box
        self.save_trajectory_box = SaveTrajectoryBox(parent=self.main, home=self.home)
        self.builder.get_object("folder_chooser_box").pack_end(
            self.save_trajectory_box.box, True, True, 0
        )
        try:
            working_folder = self.p_session.psystem[self.p_session.active_id].e_working_folder
            self.save_trajectory_box.set_folder(working_folder)
            output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
            self.save_trajectory_box.set_filename(output_name)
        except Exception:
            self.save_trajectory_box.set_folder(HOME)
        self.save_trajectory_box.set_active()

        # Setup buttons
        self.builder.get_object("button_run").connect("clicked", self.run)
        self.builder.get_object("button_cancel").connect("clicked", self.close_window)

        # Finalize window
        self.window.show_all()
        self.methods_combo.set_active(0)
        self.Visible = True

    def close_window(self, *args) -> None:
        """Close the Molecular Dynamics setup window."""
        self.window.destroy()
        self.Visible = False

    def _starting_coordinates_model_update(self, init: bool = False) -> None:
        """Update the coordinates combobox model based on the active system."""
        e_id = self.main.p_session.active_id
        model = self.main.vobject_liststore_dict[e_id]
        self.combobox_starting_coordinates.set_model(model)
        self.combobox_starting_coordinates.set_active(len(model) - 1)

    def run(self, button: Gtk.Button) -> None:
        """Collect user inputs and launch the Molecular Dynamics simulation."""
        
        # MD configuration
        MD_method = {0: "Verlet", 1: "LeapFrog", 2: "Langevin"}
        temp_scale_options = {0: "constant", 1: "linear", 2: "exponential"}
        
        
        # Extract parameters from GUI
        integrator_id = self.builder.get_object("comobobox_md_integrator").get_active()
        number_of_steps = int(self.builder.get_object("entry_number_of_steps").get_text())
        temp_scale_id = self.builder.get_object("combobox_temperature_scale_options").get_active()
        temp_start = float(self.builder.get_object("entry_temp_start").get_text())
        
        if temp_scale_id ==1 or  temp_scale_id ==2:
            temp_end = float(self.builder.get_object("entry_temp_end").get_text())
        else:
            temp_end = None
        temp_scale_factor = int(self.builder.get_object("entry_temp_scale_factor").get_text())
        time_step = float(self.builder.get_object("entry_time_step").get_text())
        log_frequence = int(self.builder.get_object("entry_log_frequency").get_text())
        random_seed = int(self.builder.get_object("entry_random_seed").get_text())
        collision_frequency = float(self.builder.get_object("entry_collision_frequency").get_text())
        temp_coupling = float(self.builder.get_object("entry_temp_coupling").get_text())

        pressure_control = self.builder.get_object("check_pressure_control").get_active()
        pressure = float(self.builder.get_object("entry_pressure").get_text())
        pressure_coupling = float(self.builder.get_object("entry_pressure_coupling").get_text())

        if temp_scale_id == 0:  # "constant" → no end temperature
            temp_end = None

        # Update system coordinates from selected vobject
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            model = self.combobox_starting_coordinates.get_model()
            _, vobject_id = model[tree_iter][:2]
            vobject = self.main.vm_session.vm_objects_dic[vobject_id]
            self.main.p_session.set_psystem_coordinates_from_vobject(vobject)



        parameters: dict = {
            "folder": HOME,
            "integrator": MD_method[integrator_id],
            "logFrequency": log_frequence,
            "seed": random_seed,
            "normal_deviate_generator": None,
            "steps": number_of_steps,
            "timeStep": time_step,
            "trajectories": None,
            # Temperature scaling
            "temperatureScaleFrequency": temp_scale_factor,
            "temperatureScaleOption": temp_scale_options[temp_scale_id],
            "temperatureStart": temp_start,
            "temperatureStop": temp_end,
            # LeapFrog/Pressure-related
            "pressure": pressure,
            "pressureControl": pressure_control,
            "pressureCoupling": pressure_coupling,
            # Temperature coupling
            "temperature": temp_start,
            "temperatureControl": True,
            "temperatureCoupling": temp_coupling,
            # Langevin-specific
            "collisionFrequency": collision_frequency,
            "simulation_type": "Molecular_Dynamics",
            "restraints": None,
        }

        # Trajectory saving
        if self.save_trajectory_box.get_active():
            parameters.update(
                {
                    "trajectory_name": self.save_trajectory_box.get_filename(),
                    "folder": self.save_trajectory_box.get_folder(),
                    "trajectory_format": self.save_trajectory_box.get_format(),
                    "trajectory_frequency": self.save_trajectory_box.get_trajectory_frequency(),
                }
            )
        else:
            parameters["trajectory_name"] = None
        
        parameters['obj1_key6'] = vobject.key6
        
        # Save last run parameters for rerun
        self.last_parameters = parameters.copy()

        # Run the simulation
        self.p_session.run_simulation(parameters)

        # Close the window
        self.window.destroy()
        self.Visible = False

    def combobox_change(self, widget: Gtk.ComboBox | None = None) -> None:
        """Handle changes in integrator or temperature scaling comboboxes."""
        if widget == self.builder.get_object("comobobox_md_integrator"):
            active = widget.get_active()
            # Reset visibility of relevant fields
            if active == 0:  # Velocity Verlet
                self._toggle_widgets(show=[], hide=[
                    "check_pressure_control", "entry_pressure",
                    "entry_temp_coupling", "label_temp_coupling",
                    "label_pressure_coupling", "entry_pressure_coupling",
                    "label_collision_freq", "entry_collision_frequency",
                ])
                self.temp_scale_options_combo.set_sensitive(True)

            elif active == 1:  # Leap Frog
                self._toggle_widgets(show=[
                    "check_pressure_control", "entry_pressure",
                    "label_pressure_coupling", "entry_pressure_coupling",
                    "entry_temp_coupling", "label_temp_coupling",
                ], hide=[
                    "label_collision_freq", "entry_collision_frequency",
                ])
                self.temp_scale_options_combo.set_active(0)
                self.temp_scale_options_combo.set_sensitive(False)

            elif active == 2:  # Langevin
                self._toggle_widgets(show=[
                    "label_collision_freq", "entry_collision_frequency",
                ], hide=[
                    "check_pressure_control", "entry_pressure",
                    "entry_temp_coupling", "label_temp_coupling",
                    "label_pressure_coupling", "entry_pressure_coupling",
                ])
                self.temp_scale_options_combo.set_active(0)
                self.temp_scale_options_combo.set_sensitive(False)

        elif widget == self.builder.get_object("combobox_temperature_scale_options"):
            temp_scale_id = widget.get_active()
            sensitive = temp_scale_id != 0
            self.builder.get_object("entry_temp_end").set_sensitive(sensitive)
            self.builder.get_object("label_temp_end").set_sensitive(sensitive)

    def _toggle_widgets(self, show: list[str] | None = None, hide: list[str] | None = None) -> None:
        """Utility function to show or hide multiple widgets."""
        show = show or []
        hide = hide or []
        for name in show:
            self.builder.get_object(name).show()
        for name in hide:
            self.builder.get_object(name).hide()

    def update_working_folder_chooser(self, folder: str | None = None) -> None:
        """Update the working folder shown in the SaveTrajectoryBox."""
        if folder:
            self.save_trajectory_box.set_folder(folder=folder)
        else:
            folder = self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder
            if folder:
                self.save_trajectory_box.set_folder(folder=folder)

    def update(self, parameters: dict | None = None) -> None:
        """Refresh window contents when system state changes."""
        
        if self.Visible:
            self._starting_coordinates_model_update()
            self.update_working_folder_chooser()
            if self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
                self.save_trajectory_box.set_filename(output_name)

    def restore_the_parameters_to_the_window(self, parameters: dict | None = None) -> None:
        """
        Reapply previously used simulation parameters back into the GUI widgets.
        If parameters are not provided, use the last stored ones from `run`.
        """
        pprint(parameters)
        
        if parameters is None:
            parameters = self.last_parameters

        if not parameters:
            print("No previous parameters available to rerun.")
            return
        
        
        ##--------------------------------------------------------------
        self.combobox_starting_coordinates.set_active(parameters['cb1_active'])
        #--------------------------------------------------------------
        
        # Refill widgets with stored values
        integrator_map = {"Verlet": 0, "LeapFrog": 1, "Langevin": 2}
        temp_scale_map = {"constant": 0, "linear": 1, "exponential": 2}

        self.methods_combo.set_active(integrator_map.get(parameters["integrator"], 0))
        self.builder.get_object("entry_number_of_steps").set_text(str(parameters["steps"]))
        self.temp_scale_options_combo.set_active(
            temp_scale_map.get(parameters["temperatureScaleOption"], 0)
        )
        self.builder.get_object("entry_temp_start").set_text(str(parameters["temperatureStart"]))
        self.builder.get_object("entry_temp_end").set_text(
            str(parameters["temperatureStop"]) if parameters["temperatureStop"] else ""
        )
        self.builder.get_object("entry_temp_scale_factor").set_text(
            str(parameters["temperatureScaleFrequency"])
        )
        self.builder.get_object("entry_time_step").set_text(str(parameters["timeStep"]))
        self.builder.get_object("entry_log_frequency").set_text(str(parameters["logFrequency"]))
        self.builder.get_object("entry_random_seed").set_text(str(parameters["seed"]))
        self.builder.get_object("entry_collision_frequency").set_text(
            str(parameters["collisionFrequency"])
        )
        self.builder.get_object("entry_temp_coupling").set_text(
            str(parameters["temperatureCoupling"])
        )

        # Pressure-related
        self.builder.get_object("check_pressure_control").set_active(
            parameters["pressureControl"]
        )
        self.builder.get_object("entry_pressure").set_text(str(parameters["pressure"]))
        self.builder.get_object("entry_pressure_coupling").set_text(
            str(parameters["pressureCoupling"])
        )

        # Trajectory saving
        if parameters.get("trajectory_name"):
            self.save_trajectory_box.set_folder(parameters["folder"])
            self.save_trajectory_box.set_filename(parameters["trajectory_name"])
            self.save_trajectory_box.set_format(parameters["trajectory_format"])
            self.save_trajectory_box.set_trajectory_frequency(parameters["trajectory_frequency"])
            self.save_trajectory_box.set_active(True)
        else:
            self.save_trajectory_box.set_active(False)











class MolecularDynamicsWindow_OLD():

    def __init__(self, main = None):
        """ Class initialiser """
        self.windon_md_main      = main
        
        self.main                = main
        self.p_session           = main.p_session
        
        self.home                = main.home
        self.Visible             =  False        
        self.residue_liststore   = Gtk.ListStore(str, str, str)
        self.job_liststore       = Gtk.ListStore(str, str, str)
        
        
        self.md_intergators_dict = {
                                   0:"Velocity Verlet Dynamics", 
                                   1:"Leap Frog Dynamics",
                                   2:"Langevin Dynamics"
                                   }
                       
        self.temp_scale_options = {0: "constant", 
                                   1: "linear",
                                   2: "exponential"}
        
        self.sym_tag = 'mol_dyn'
        self.joblist = []

    def open_window (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'src/gui/windows/simulation/molecular_dynamics_window.glade'))

            self.window = self.builder.get_object('setup_md_window')
            self.window.set_title('Molecular Dynamics Setup')
            self.window.set_keep_above(True)
            self.window.connect("destroy", self.close_window)
        
                
            #------------------------------------------------------------------------------------------
            self.builder.get_object('check_pressure_control').hide()
            self.builder.get_object('entry_pressure').hide()
            self.builder.get_object('entry_pressure').hide()

            self.builder.get_object('entry_temp_coupling').hide()
            self.builder.get_object('label_temp_coupling').hide()

            self.builder.get_object('label_pressure_coupling').hide()
            self.builder.get_object('entry_pressure_coupling').hide()

            self.builder.get_object('label_collision_freq').hide()
            self.builder.get_object('entry_collision_frequency').hide()
            #------------------------------------------------------------------------------------------


            
            #
            #'''--------------------------------------------------------------------------------------------'''
            # combobox
            #self.starting_coords_liststore = self.main.vobject_liststore_dict[self.main.p_session.active_id]
            #self.combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
            #self.combobox_starting_coordinates.set_model(self.starting_coords_liststore)
            #renderer_text = Gtk.CellRendererText()
            #self.combobox_starting_coordinates.pack_start(renderer_text, True)
            #self.combobox_starting_coordinates.add_attribute(renderer_text, "text", 0)
            #
            #self.starting_coords_liststore = self.main.vobject_liststore_dict[self.main.p_session.active_id]
            self.box_coordinates = self.builder.get_object('box_coordinates')
            self.combobox_starting_coordinates = CoordinatesComboBox() #self.builder.get_object('coordinates_combobox')
            self.box_coordinates.pack_start(self.combobox_starting_coordinates, False, False, 0)
            self._starting_coordinates_model_update(init = True)
                       
            
            size = len(self.combobox_starting_coordinates.get_model())
            self.combobox_starting_coordinates.set_active(size-1)
            #'''--------------------------------------------------------------------------------------------'''

            
            
            #'''--------------------------------------------------------------------------------------------'''
            self.md_methods_liststore = Gtk.ListStore(str)
            for key, method in self.md_intergators_dict.items():
                self.md_methods_liststore.append([method])
            
            self.methods_combo = self.builder.get_object('comobobox_md_integrator')
            self.methods_combo.set_model(self.md_methods_liststore)
            #
            renderer_text = Gtk.CellRendererText()
            self.methods_combo.pack_start(renderer_text, True)
            self.methods_combo.add_attribute(renderer_text, "text", 0)
            #'''--------------------------------------------------------------------------------------------'''
            self.methods_combo.connect("changed", self.combobox_change)
            
            #'''--------------------------------------------------------------------------------------------'''


            
            #'''--------------------------------------------------------------------------------------------'''
            self.temp_scale_option_store = Gtk.ListStore(str)
            temp_scale_options = ["constant", "linear","exponential"]
            #temp_scale_options
            for key,temp_scale_option in self.temp_scale_options.items():
                self.temp_scale_option_store.append([temp_scale_option])
            
            self.temp_scale_options_combo = self.builder.get_object('combobox_temperature_scale_options')
            self.temp_scale_options_combo.set_model(self.temp_scale_option_store)
            renderer_text = Gtk.CellRendererText()
            self.temp_scale_options_combo.pack_start(renderer_text, True)
            self.temp_scale_options_combo.add_attribute(renderer_text, "text", 0)
            #'''--------------------------------------------------------------------------------------------'''
            
            
            #'''--------------------------------------------------------------------------------------------'''
            self.temp_scale_options_combo.connect("changed", self.combobox_change)
            self.temp_scale_options_combo.set_active(0)
            #'''--------------------------------------------------------------------------------------------'''
            
            
            #'''--------------------------------------------------------------------------------------------'''
            self.save_trajectory_box = SaveTrajectoryBox(parent = self.main, home = self.home)
            self.builder.get_object('folder_chooser_box').pack_end(self.save_trajectory_box.box, True, True, 0)
            
            try:
                self.save_trajectory_box.set_folder(self.p_session.psystem[self.p_session.active_id].e_working_folder)
                output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
                self.save_trajectory_box.set_filename (output_name )
            except:
                self.save_trajectory_box.set_folder(HOME)
            
            '''
            if self.p_session.psystem[self.p_session.active_id].e_working_folder:
                self.save_trajectory_box.set_folder(self.p_session.psystem[self.p_session.active_id].e_working_folder)
            
            else:
                self.save_trajectory_box.set_folder(HOME)
              
            
            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
                self.save_trajectory_box.set_filename (output_name )
            
            else:
                pass
            #'''    
            self.save_trajectory_box.set_active ()
            #'''--------------------------------------------------------------------------------------------'''


            #-------------------------------------------------------------------------------------------
            button_send_job_to_list = self.builder.get_object('button_run')
            button_send_job_to_list.connect('clicked', self.run)
            self.builder.get_object('button_cancel').connect('clicked', self.close_window)
            #-------------------------------------------------------------------------------------------
            
            
            self.window.show_all()
            self.methods_combo.set_active(0)
            
            self.Visible  = True

        else:
            self.window.present()
            
    def close_window (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False


    def _starting_coordinates_model_update (self, init = False):
        """ Function doc """
        #------------------------------------------------------------------------------------
        '''The combobox accesses, according to the id of the active system, 
        listostore of the dictionary object_list state_dict'''
        if self.Visible:

            e_id = self.main.p_session.active_id 
            self.combobox_starting_coordinates.set_model(self.main.vobject_liststore_dict[e_id])
            #------------------------------------------------------------------------------------
            size = len(self.main.vobject_liststore_dict[e_id])
            self.combobox_starting_coordinates.set_active(size-1)
            #------------------------------------------------------------------------------------
        else:
            if init:
                e_id = self.main.p_session.active_id 
                self.combobox_starting_coordinates.set_model(self.main.vobject_liststore_dict[e_id])
                #------------------------------------------------------------------------------------
                size = len(self.main.vobject_liststore_dict[e_id])
                self.combobox_starting_coordinates.set_active(size-1)
                #------------------------------------------------------------------------------------
            else:
                pass

    
    def run (self, button):
        '''
        "folder"                    : os.getcwd()             ,
        'integrator'                : MD_method[integrator_id], # verlet / leapfrog /langevin
        'logFrequency'              : log_frequence           ,
        'seed'                      : random_seed             ,
        'normal_deviate_generator'    : None                    ,
        'steps'                     : number_of_steps         ,
        'timeStep'                  : time_step               ,
        'trajectories'              : None                    ,
        
        #VelocityVerletDynamics
        'temperatureScaleFrequency' : temp_scale_factor                ,
        'temperatureScaleOption'    : temp_scale_options[temp_scale_id],  # "linear" , "constant" ,
        'temperatureStart'          : temp_start                       ,
        'temperature_stop'           : temp_end                         ,
        
        #LeapFrogDynamics
        'pressure'                  : pressure           ,  #  LeapFrogDynamics 
        'temperatureControl'           : pressure_control   ,  #  LeapFrogDynamics 
        'temperatureCoupling'          : pressure_coupling  ,  #  LeapFrogDynamics 
                                                         
        'temperature'               : temp_start         ,               
        'temperatureControl'        : True               ,  # True / False LeapFrogDynamics / LangevinDynamics
        'temperatureCoupling'       : temp_coupling      ,  #  LeapFrogDynamics / LangevinDynamics
        
        #LangevinDynamics
        'collisionFrequency'        : collision_frequency,
        
        
        #Mandatory keys in self.parameters:
        #    "MD_method"	 	 : string containing the integrator algorithm name
        #    "protocol" 		 : string indicating if is a normal run or for heating
        #    "nsteps"   		 : Number of steps to be taken in the simulation
        #    "trajectory_name":
        #Optinal  :
        #    "temperature" 			  : float with the simulation temperature. If not passed we assume 300.15K as default.
        #    "coll_freq"  			  : integer with the colision frequency. Generally set for Langevin integrator. 
        #    "pressure"   			  : float with the simulation pressure. If not passed we assume 1.0bar as default.
        #    "pressure_coupling"		  : boolean indicating if is to control the simulation pressure.
        #    'temperatureScaleOption': string with the type of temperature scaling. Default is 'linear' ( relevant for "heating" protocol)
        #    "temperature_scale"		  : float with the  temperature scaling step. Default is 10K  ( relevant for "heating" protocol)
        #    "start_temperatue"		  : float with the start temperature for heating protocol
        #    "timeStep"   			  : float indicating the size of integration time step. 0.001 ps is taken as default.					
        #    "sampling_factor"		  : integer indicating in which frequency to save/collect structure/data. default 0.
        #    "seed"					  : integer indicating the seed for rumdomness of the simulations.
        #    "logFrequency"     	  : integer indicating the frequency of the screen log output.
        
        '''		

        
        integrator_id       = self.builder.get_object('comobobox_md_integrator').get_active()
        number_of_steps     = int(self.builder.get_object('entry_number_of_steps').get_text())
        temp_scale_id       = self.builder.get_object('combobox_temperature_scale_options').get_active()
        temp_start          = float(self.builder.get_object('entry_temp_start').get_text())
        temp_end            = float(self.builder.get_object('entry_temp_end').get_text())
        temp_scale_factor   = int(self.builder.get_object('entry_temp_scale_factor').get_text())
        time_step           = float(self.builder.get_object('entry_time_step').get_text())
        log_frequence       = int(self.builder.get_object('entry_log_frequency').get_text())
        random_seed         = int(self.builder.get_object('entry_random_seed').get_text())
        collision_frequency = float(self.builder.get_object('entry_collision_frequency').get_text())
        
        
        temp_coupling       =  float(self.builder.get_object('entry_temp_coupling').get_text())
        
        if self.builder.get_object('check_pressure_control').get_active():
            pressure_control = True
        else:
            pressure_control = False
       
        pressure          = float(self.builder.get_object('entry_pressure').get_text())
        pressure_coupling = float(self.builder.get_object('entry_pressure_coupling').get_text())
        
        if temp_scale_id == 0:
            temp_end = None
        else:
            pass
        
        '''---------------------------------------------------------------------------------'''
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]
            vobject = self.main.vm_session.vm_objects_dic[vobject_id]
            
            '''This function imports the coordinates of a vobject into the dynamo system in memory.''' 
            self.main.p_session.set_psystem_coordinates_from_vobject(vobject)
        '''---------------------------------------------------------------------------------'''

        
        
        MD_method = {
                     0 : "Verlet"   ,
                     1 : "LeapFrog" ,
                     2 : "Langevin" ,
                     }
        
        temp_scale_options = {
                              0 : 'constant'    ,
                              1 : 'linear'      ,
                              2 : 'exponential' , 
                             }
        
        
        temp_scale_id
        
        parameters = {
                    
                    "folder"                    : HOME                    ,
                    'integrator'                : MD_method[integrator_id], # verlet / leapfrog /langevin
                    'logFrequency'             : log_frequence           ,
                    'seed'                      : random_seed             ,
                    'normal_deviate_generator'    : None                    ,
                    'steps'                     : number_of_steps         ,
                    'timeStep'                 : time_step               ,
                    'trajectories'              : None                    ,
                    
                    #VelocityVerletDynamics
                    'temperatureScaleFrequency' : temp_scale_factor                ,
                    'temperatureScaleOption'    : temp_scale_options[temp_scale_id],  # "linear" , "constant" ,
                    'temperatureStart'          : temp_start                       ,
                    'temperatureStop'           : temp_end                         ,
                    
                    #LeapFrogDynamics
                    'pressure'                  : pressure           ,  #  LeapFrogDynamics 
                    'pressureControl'           : pressure_control   ,  #  LeapFrogDynamics 
                    'pressureCoupling'          : pressure_coupling  ,  #  LeapFrogDynamics 
                                                                     
                    'temperature'               : temp_start         ,               
                    'temperatureControl'        : True               ,  # True / False LeapFrogDynamics / LangevinDynamics
                    'temperatureCoupling'       : temp_coupling      ,  #  LeapFrogDynamics / LangevinDynamics
                    
                    #LangevinDynamics
                    'collisionFrequency'        : collision_frequency,
                    
                    }
        
        
        
 
        parameters["simulation_type"] = "Molecular_Dynamics"
        
        if self.save_trajectory_box.get_active():
            parameters['trajectory_name'  ]    = self.save_trajectory_box.get_filename() 
            parameters['folder']               = self.save_trajectory_box.get_folder()  
            parameters['trajectory_format']    = self.save_trajectory_box.get_format()  
            parameters['trajectory_frequency'] = self.save_trajectory_box.get_trajectory_frequency()  
        else:
            parameters['trajectory_name'  ]    = None
        
        # restraints
        parameters['restraints'] = None
        
        
        #pprint(parameters)
        self.p_session.run_simulation(parameters)
        self.window.destroy()
        self.Visible    =  False
        
    def update_working_folder_chooser (self, folder = None):
        """ Function doc """
        if folder:
            #print('update_working_folder_chooser')
            self.save_trajectory_box.set_folder(folder = folder)
        else:
            self.save_trajectory_box.set_folder(folder = HOME)

    def combobox_change (self, widget = None):
        """ Function doc """


        
        
        if widget  == self.builder.get_object('comobobox_md_integrator'):
            
            #velocity verlet molecular dynamics
            if 0 == widget.get_active():
                self.builder.get_object('check_pressure_control').hide()
                self.builder.get_object('entry_pressure').hide()
                self.builder.get_object('entry_pressure').hide()
                
                
                self.builder.get_object('entry_temp_coupling').hide()
                self.builder.get_object('label_temp_coupling').hide()
                
                self.builder.get_object('label_pressure_coupling').hide()
                self.builder.get_object('entry_pressure_coupling').hide()
                
                self.builder.get_object('label_collision_freq').hide()
                self.builder.get_object('entry_collision_frequency').hide()
                
                self.temp_scale_options_combo.set_sensitive(True)
                
            #Leap Frog molecular dynamics
            elif 1 == widget.get_active():
                self.builder.get_object('check_pressure_control').show()
                self.builder.get_object('entry_pressure').show()
                
                self.builder.get_object('label_pressure_coupling').show()
                self.builder.get_object('entry_pressure_coupling').show()
                
                
                self.builder.get_object('entry_temp_coupling').show()
                self.builder.get_object('label_temp_coupling').show()
                
                
                self.builder.get_object('label_collision_freq').hide()
                self.builder.get_object('entry_collision_frequency').hide()
                
                self.temp_scale_options_combo.set_active(0)
                self.temp_scale_options_combo.set_sensitive(False)
                
            #Langevin molecular dynamics
            elif 2 == widget.get_active():
                self.builder.get_object('check_pressure_control').hide()
                self.builder.get_object('entry_pressure').hide()
                
                self.builder.get_object('entry_temp_coupling').hide()
                self.builder.get_object('label_temp_coupling').hide()
                
                self.builder.get_object('label_pressure_coupling').hide()
                self.builder.get_object('entry_pressure_coupling').hide()
                
                self.builder.get_object('label_collision_freq').show()
                self.builder.get_object('entry_collision_frequency').show()
            
                self.temp_scale_options_combo.set_active(0)
                self.temp_scale_options_combo.set_sensitive(False)
            else:
                pass

        
        
        
        elif widget == self.builder.get_object('combobox_temperature_scale_options'):
            temp_scale_id = widget.get_active()
            
            if temp_scale_id == 0:
                self.builder.get_object('entry_temp_end').set_sensitive(False)
                self.builder.get_object('label_temp_end').set_sensitive(False)
            else:
                self.builder.get_object('entry_temp_end').set_sensitive(True)
                self.builder.get_object('label_temp_end').set_sensitive(True)

        else:
            pass

    def update_working_folder_chooser (self, folder = None):
        """ Function doc """
        if folder:
            #print('update_working_folder_chooser')
            self.save_trajectory_box.set_folder(folder = folder)
        else:
            
            folder = self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder
            if folder:
                self.save_trajectory_box.set_folder(folder = folder)
            else:
                pass


    def update (self, parameters = None):
        """ Function doc """
        
        if self.Visible:
            self._starting_coordinates_model_update()
            self.update_working_folder_chooser()
        
        
        
            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
                self.save_trajectory_box.set_filename (output_name )
            else:
                pass
