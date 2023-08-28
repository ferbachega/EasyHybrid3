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
from pprint import pprint
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
#from GTKGUI.gtkWidgets.filechooser import FileChooser
from gui.gtk_widgets import FolderChooserButton
from gui.gtk_widgets import SaveTrajectoryBox
from gui.gtk_widgets import CoordinatesComboBox

import gc
import os

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')



class MolecularDynamicsWindow():

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

    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'gui/windows/molecular_dynamics_window_new.glade'))

            self.window = self.builder.get_object('setup_md_window')
            self.window.set_title('Molecular Dynamics Setup')
            self.window.set_keep_above(True)
            self.window.connect("destroy", self.CloseWindow)
        
                
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
            
            if self.p_session.psystem[self.p_session.active_id].e_working_folder:
                self.save_trajectory_box.set_folder(self.p_session.psystem[self.p_session.active_id].e_working_folder)
            
            else:
                self.save_trajectory_box.set_folder(HOME)
                
            
            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
                self.save_trajectory_box.set_filename (output_name )
            
            else:
                pass
                
            self.save_trajectory_box.set_active ()
            #'''--------------------------------------------------------------------------------------------'''


            #-------------------------------------------------------------------------------------------
            button_send_job_to_list = self.builder.get_object('button_run')
            button_send_job_to_list.connect('clicked', self.run)
            self.builder.get_object('button_cancel').connect('clicked', self.CloseWindow)
            #-------------------------------------------------------------------------------------------
            
            
            self.window.show_all()
            self.methods_combo.set_active(0)
            
            self.Visible  = True

        else:
            self.window.present()
            
    def CloseWindow (self, button, data  = None):
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
        'normalDeviateGenerator'    : None                    ,
        'steps'                     : number_of_steps         ,
        'timeStep'                  : time_step               ,
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
        #    "temperature_scale_option": string with the type of temperature scaling. Default is 'linear' ( relevant for "heating" protocol)
        #    "temperature_scale"		  : float with the  temperature scaling step. Default is 10K  ( relevant for "heating" protocol)
        #    "start_temperatue"		  : float with the start temperature for heating protocol
        #    "timeStep"   			  : float indicating the size of integration time step. 0.001 ps is taken as default.					
        #    "sampling_factor"		  : integer indicating in which frequency to save/collect structure/data. default 0.
        #    "seed"					  : integer indicating the seed for rumdomness of the simulations.
        #    "log_frequency"     	  : integer indicating the frequency of the screen log output.
        
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
            self.main.p_session.get_coordinates_from_vobject_to_pDynamo_system(vobject)
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
                    'logFrequency'              : log_frequence           ,
                    'seed'                      : random_seed             ,
                    'normalDeviateGenerator'    : None                    ,
                    'steps'                     : number_of_steps         ,
                    'timeStep'                  : time_step               ,
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
        self._starting_coordinates_model_update()
        if self.Visible:
            self.update_working_folder_chooser()
        
        
        
            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
                self.save_trajectory_box.set_filename (output_name )
            else:
                pass
