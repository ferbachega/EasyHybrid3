#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  easyhybrid_pDynamo_selection.py
#  
#  Copyright 2022 Fernando <fernando@winter>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
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
#  
import gi
import threading
import multiprocessing

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
import gc
import os
from gui.widgets.custom_widgets import FolderChooserButton
from gui.widgets.custom_widgets import SaveTrajectoryBox

from gui.widgets.custom_widgets import CoordinatesComboBox

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')

##====================================================================================
class GeometryOptimization(Gtk.Window):
    """ Class doc """
    
    def __init__(self, main = None ):
        """ Class initialiser """
        self.main     = main
        self.home     = main.home 
        self.p_session= main.p_session 
        self.Visible  =  False        
        #self.residue_liststore = Gtk.ListStore(str, str, str)
        self.sym_tag  = 'geo_opt'
        self.opt_methods = { 
                            0 : 'ConjugatedGradient',
                            1 : 'SteepestDescent'   ,
                            2 : 'LFBGS'             ,
                            3 : 'QuasiNewton'       ,
                            4 : 'FIRE'              ,
                             }

    
    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = self.main.builder #Gtk.Builder()
            
            #self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/geometry_optimization_window.glade'))
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'src/gui/windows/simulation/geometry_optimization_window.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('geometry_optimization_window')
            self.window.set_title('Geometry Optmization Window')
            self.window.set_keep_above(True)
            
            
            '''--------------------------------------------------------------------------------------------'''
            self.method_store = Gtk.ListStore(str)
            methods = [
                        'Conjugated-Gradient',
                        'Steepest-Descent'   ,
                        'L-FBGS'             ,
                        'Quasi-Newton'       ,
                        'FIRE'               ,

                ]
            for method in methods:
                self.method_store.append([method])
                #print (method)
            self.methods_combo = self.builder.get_object('combobox_geo_opt')
            self.methods_combo.set_model(self.method_store)
            self.methods_combo.connect("changed", self.on_name_combo_changed)
            self.methods_combo.set_model(self.method_store)
            
            renderer_text = Gtk.CellRendererText()
            self.methods_combo.pack_start(renderer_text, True)
            self.methods_combo.add_attribute(renderer_text, "text", 0)
            self.methods_combo.set_active(0)            
            
            # - - - - - - - - - - - - - Starting Coordinates ComboBox - - - - - - - - - - - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            #----------------------------------------------------------------------------------------------
            self.box_coordinates = self.builder.get_object('box_coordinates')
            self.combobox_starting_coordinates = CoordinatesComboBox() #self.builder.get_object('coordinates_combobox')
            self.combobox_starting_coordinates.connect("changed", self.on_name_combo_changed)
            self.box_coordinates.pack_start(self.combobox_starting_coordinates, False, False, 0)
            self._starting_coordinates_model_update(init = True)
            #----------------------------------------------------------------------------------------------
            
            


            self.save_trajectory_box = SaveTrajectoryBox(parent = self.main, home = self.home)
            self.builder.get_object('geo_opt_parm_box').pack_end(self.save_trajectory_box.box, True, True, 0)
                        
            # updating data 
            
            #------------------------------------------------------------------------------------------------
            #print(self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder)
            if self.main.p_session.psystem[self.main.p_session.active_id]:
                if self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder == None:
                    folder = HOME
                else:
                    folder = self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder
                self.update_working_folder_chooser(folder = folder)            
            #------------------------------------------------------------------------------------------------

            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
                self.save_trajectory_box.set_filename (output_name )
            else:
                pass

            self.window.show_all()
            self.Visible  = True   
            #self.window.set_size_request(900, 900)
        else:
            self.window.present()
    
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False

    def run_opt(self, button):
        """
        Run a geometry optimization based on the parameters set in the GUI.
        """

        # ---------------------------------------------------------------
        # Default simulation parameters
        simParameters = {
            "simulation_type"     : "Geometry_Optimization",
            "trajectory_name"     : None,
            "dialog"              : False,
            "folder"              : os.getcwd(),
            "optimizer"           : "ConjugatedGradient",
            "maximumIterations"   : 600,
            "logFrequency"        : 10,
            "save_frequency"      : 10,
            "rmsGradientTolerance": 0.1,
            "save_format"         : None,
            "save_traj"           : False,
            "save_pdb"            : False,
        }

        # ---------------------------------------------------------------
        # Get starting coordinates from the combobox
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            # Select vismol object from combobox model
            model = self.combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]
            vobject = self.main.vm_session.vm_objects_dic[vobject_id]

            # Import coordinates of the vismol object into the pDynamo system
            self.main.p_session.get_coordinates_from_vobject_to_pDynamo_system(vobject)

        # ---------------------------------------------------------------
        # Update parameters from GUI entries
        simParameters["optimizer"]           = self.opt_methods[
            self.builder.get_object('combobox_geo_opt').get_active()
        ]
        simParameters["logFrequency"]        = int(
            self.builder.get_object('entry_log_frequence').get_text()
        )
        simParameters["maximumIterations"]   = int(
            self.builder.get_object('entry_max_int').get_text()
        )
        simParameters["rmsGradientTolerance"] = float(
            self.builder.get_object('entry_rmsd_tol').get_text()
        )
        simParameters["vobject_name"]        = self.save_trajectory_box.builder.get_object(
            'entry_trajectory_name'
        ).get_text()

        # ---------------------------------------------------------------
        # Configure trajectory saving (if enabled)
        if self.save_trajectory_box.builder.get_object('checkbox_save_traj').get_active():
            simParameters["save_traj"]       = True
            simParameters["dialog"]          = True
            simParameters["folder"]          = self.save_trajectory_box.folder_chooser_button.get_folder()
            simParameters["trajectory_name"] = self.save_trajectory_box.builder.get_object(
                'entry_trajectory_name'
            ).get_text()
            simParameters["save_frequency"]  = int(
                self.save_trajectory_box.builder.get_object('entry_trajectory_frequency').get_text()
            )

            # Append file extension to trajectory name
            simParameters["trajectory_name"] += ".ptGeo"

            # Select save format
            saveFormat = self.save_trajectory_box.builder.get_object('combobox_format').get_active()
            if   saveFormat == 0: simParameters["save_format"] = ".ptGeo"
            elif saveFormat == 1: simParameters["save_format"] = ".mdcrd"
            elif saveFormat == 2: simParameters["save_format"] = ".dcd"
            elif saveFormat == 3: simParameters["save_format"] = ".xyz"

            # Update working folder in the active system
            self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder = simParameters["folder"]

        # ---------------------------------------------------------------
        # Check if the output folder exists
        if not os.path.exists(simParameters['folder']):
            self.run_dialog()
            return None

        # ---------------------------------------------------------------
        # Run simulation
        self.main.p_session.run_simulation(parameters=simParameters)

        # Close window
        self.window.destroy()
        self.Visible = False

    
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

    
    def on_name_combo_changed(self, widget):
        """ Function doc """
        #print('eba - apagar')
    #=================================================================================
    
    def run_dialog (self):
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


    def update (self, parameters = None):
        """ Function doc """
        self._starting_coordinates_model_update()
        if self.Visible:
            self.update_working_folder_chooser()
            

            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
                self.save_trajectory_box.set_filename (output_name )
            else:
                pass


    def update_working_folder_chooser (self, folder = None):
        """ Function doc """
        #folder = self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder
        if folder:
            #print('update_working_folder_chooser')
            self.save_trajectory_box.set_folder(folder = folder)
        else:
            
            folder = self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder
            if folder:
                self.save_trajectory_box.set_folder(folder = folder)
            else:
                pass
            #self.save_trajectory_box.set_folder(folder = folder)
            #self.save_trajectory_box.set_folder(folder = self.main.p_session.systems[self.main.p_session.active_id]['working_folder'])
    
    def restore_the_parameters_to_the_window(self, parameters):
        """
        Restore exported parameters into the GTK window.
        """

        # ---------------------------------------------------------------
        # Select optimization method
        e_id = parameters['system']

        # Invert the opt_methods dictionary to get tag -> key
        self.opt_tags = {value: key for key, value in self.opt_methods.items()}

        # Set the optimizer method in the combo box
        optimizer_key = self.opt_tags.get(parameters['optimizer'])
        if optimizer_key is not None:
            self.methods_combo.set_active(optimizer_key)

        # ---------------------------------------------------------------
        # Configure basic parameters
        self.builder.get_object('entry_log_frequence').set_text(
            str(parameters['logFrequency'])
        )
        self.builder.get_object('entry_max_int').set_text(
            str(parameters['maximumIterations'])
        )
        self.builder.get_object('entry_rmsd_tol').set_text(
            str(parameters['rmsGradientTolerance'])
        )

        # ---------------------------------------------------------------
        # Configure trajectory settings (if enabled)
        if parameters.get('save_traj'):
            checkbox = self.save_trajectory_box.builder.get_object('checkbox_save_traj')
            checkbox.set_active(True)

            # Set folder
            self.save_trajectory_box.folder_chooser_button.set_folder(
                parameters.get('folder', "")
            )

            # Set trajectory name (if provided)
            trajectory_name = parameters.get('trajectory_name')
            if trajectory_name:
                entry_name = self.save_trajectory_box.builder.get_object('entry_trajectory_name')
                entry_name.set_text(trajectory_name)

            # Set trajectory frequency
            entry_freq = self.save_trajectory_box.builder.get_object('entry_trajectory_frequency')
            entry_freq.set_text(str(parameters.get('save_frequency', "")))

            # Set trajectory format (default: first option)
            combobox_format = self.save_trajectory_box.builder.get_object('combobox_format')
            combobox_format.set_active(0)

        # ---------------------------------------------------------------
        # Configure starting coordinates
        #model = self.combobox_starting_coordinates.get_model()
        #name, vobject_id = model[tree_iter][:2]
        #vobject = self.main.vm_session.vm_objects_dic[vobject_id]
