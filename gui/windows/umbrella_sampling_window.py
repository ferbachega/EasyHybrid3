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

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
import gc
import os
from gui.gtk_widgets import FolderChooserButton
from gui.gtk_widgets import SaveTrajectoryBox

from gui.windows.PES_scan_window            import compute_sigma_a1_a3 
from gui.windows.PES_scan_window            import get_dihedral 
from gui.windows.PES_scan_window            import get_angle 
from gui.windows.PES_scan_window            import get_distance 

from gui.windows.PES_scan_window            import atomic_dic 
from gui.windows.PES_scan_window            import texto_d1 
from gui.windows.PES_scan_window            import texto_d2d1 


HOME        = os.environ.get('HOME')



##====================================================================================
class UmbrellaSamplingWindow(Gtk.Window):
    """ Class doc """
    
    def __init__(self, main = None ):
        """ Class initialiser """
        self.main       = main
        self.home       = main.home
        self.p_session  = main.p_session 
        self.vm_session = main.vm_session 
        self.Visible    =  False        
        
        self.input_types  = {
                            0 : 'From Coordinates (sequential)',
                            1 : 'From Trajectory (parallel)'   ,
                           }
        
        
        self.reaction_coord_types = { 
                            0: "simple distance", 
                            1: "multiple distance", 
                            2: 'dihedral'
                            
                            }
        
        self.opt_methods = { 
                            0 : 'ConjugatedGradient',
                            1 : 'SteepestDescent'   ,
                            2 : 'LFBGS'             ,
                            3 : 'QuasiNewton'       ,
                            4 : 'FIRE'              ,
                             }

    
        self.md_methods = {
                           0:"Velocity Verlet Dynamics", 
                           1:"Leap Frog Dynamics",
                           2:"Langevin Dynamics"
                            }
    
    
    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = self.main.builder #Gtk.Builder()
            
            #self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/geometry_optimization_window.glade'))
            self.builder = Gtk.Builder()                
            self.builder.add_from_file(os.path.join (self.home,'gui/windows/umbrella_sampling_window.glade') )
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('umbrella_samppling_window')
            self.window.set_title('Umbrella Samppling Window')
            self.window.set_keep_above(True)
            self.window.connect("destroy", self.CloseWindow)
            

            '''--------------------------------------------------------------------------------------------'''
            self.input_type_combo = self.builder.get_object('combobox_input_type')
            self.input_type_store = Gtk.ListStore(str)

            for key, input_type in self.input_types.items():
                self.input_type_store.append([input_type])
            
            renderer_text = Gtk.CellRendererText()
            
            self.input_type_combo.set_model(self.input_type_store)
            self.input_type_combo.pack_start(renderer_text, True)
            self.input_type_combo.add_attribute(renderer_text, "text", 0)            
            self.input_type_combo.connect("changed", self.on_combobox_inputtype)
            '''--------------------------------------------------------------------------------------------'''
            
            # - - - - - - - - - - - - - Starting Coordinates ComboBox - - - - - - - - - - - - - - - - -
            self.combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
            #---------------------------------------------------------------------------------------------
            self._starting_coordinates_model_update(init = True)
            
            
            self.combobox_starting_coordinates.connect("changed", self.on_name_combo_changed)
            renderer_text = Gtk.CellRendererText()
            self.combobox_starting_coordinates.pack_start(renderer_text, True)
            self.combobox_starting_coordinates.add_attribute(renderer_text, "text", 0)
            #----------------------------------------------------------------------------------------------
            
            
            # - - - - - - - - - - - - - Folder Chooser Button - - - - - - - - - - - - - - - - -
            #----------------------------------------------------------------------------------------------
            self.folder_chooser_button = FolderChooserButton(main =  self.main, sel_type = 'folder', home =  self.home)
            self.builder.get_object('folder_chooser_box').pack_start(self.folder_chooser_button.btn, True, True, 0)
            system_id      = self.p_session.active_id
            #----------------------------------------------------------------------------------------------
            
            
            # - - - - - - - - - - - - - - - - Spin Button Number Of CPUS - - - - - - - - - - - - - - - - -
            #----------------------------------------------------------------------------------------------
            self.spinbutton = self.builder.get_object('ncpus_spinbutton')

            # set the range of values that the spinbutton can take
            self.spinbutton.set_range(1, 1000)

            # set the increment step when the user clicks the up/down buttons
            self.spinbutton.set_increments(1, 10)

            # set the number of decimal places to display
            self.spinbutton.set_digits(0)

            # set the value of the spinbutton
            self.spinbutton.set_value(1)
            #----------------------------------------------------------------------------------------------
            
            
            # - - - - - - - - - - - - - Reaction Coordinates 1 ComboBox - - - - - - - - - - - - - - - - -
            #----------------------------------------------------------------------------------------------
            self.combobox_reaction_coord1 = self.builder.get_object('combobox_reaction_coord1')
            self.reaction_coord1_store = Gtk.ListStore(str)

            for key, reaction_coord_types in self.reaction_coord_types.items():
                self.reaction_coord1_store.append([reaction_coord_types])
            
            self.combobox_reaction_coord1.set_model(self.reaction_coord1_store)
            renderer_text = Gtk.CellRendererText()
            self.combobox_reaction_coord1.pack_start(renderer_text, True)
            self.combobox_reaction_coord1.add_attribute(renderer_text, "text", 0)
            self.combobox_reaction_coord1.connect("changed", self.change_cb_coordType1)

            #----------------------------------------------------------------------------------------------

        
            # - - - - - - - - - - - - - Reaction Coordinates 2 ComboBox - - - - - - - - - - - - - - - - -
            #----------------------------------------------------------------------------------------------
            self.combobox_reaction_coord2 = self.builder.get_object('combobox_reaction_coord2')
            self.reaction_coord2_store = Gtk.ListStore(str)

            for key, reaction_coord_types in self.reaction_coord_types.items():
                self.reaction_coord2_store.append([reaction_coord_types])


            self.combobox_reaction_coord2.set_model(self.reaction_coord2_store)
            renderer_text = Gtk.CellRendererText()
            self.combobox_reaction_coord2.pack_start(renderer_text, True)
            self.combobox_reaction_coord2.add_attribute(renderer_text, "text", 0)
            self.combobox_reaction_coord2.connect("changed", self.change_cb_coordType2)

            #----------------------------------------------------------------------------------------------

        
        
            # - - - - - - - - - - - - - -  Geometry Optimization ComboBox - - - - - - - - - - - - - - - - -
            #----------------------------------------------------------------------------------------------
            self.opt_methods_combo = self.builder.get_object('combobox_geo_opt')
            self.opt_method_store = Gtk.ListStore(str)

            for key, method in self.opt_methods.items():
                self.opt_method_store.append([method])
            
            self.opt_methods_combo.set_model(self.opt_method_store)
            renderer_text = Gtk.CellRendererText()
            self.opt_methods_combo.pack_start(renderer_text, True)
            self.opt_methods_combo.add_attribute(renderer_text, "text", 0)
            #----------------------------------------------------------------------------------------------
            
            
            
            # - - - - - - - - - - - - - - - - Molecular Dynamics ComboBox - - - - - - - - - - - - - - - - -
            #----------------------------------------------------------------------------------------------
            self.md_integrators_combobox = self.builder.get_object('md_integrator_combobox')
            self.md_integrators_store    = Gtk.ListStore(str)

            for key, integrators in self.md_methods.items():
                self.md_integrators_store.append([integrators])
            
            self.md_integrators_combobox.set_model(self.md_integrators_store)
            renderer_text = Gtk.CellRendererText()
            self.md_integrators_combobox.pack_start(renderer_text, True)
            self.md_integrators_combobox.add_attribute(renderer_text, "text", 0)
            self.md_integrators_combobox.connect("changed", self.on_md_integrator_combobox)
            #----------------------------------------------------------------------------------------------
            
            
            #---------------------------------------------------------------------------------------------
            self.builder.get_object('button_import_picking_selection_coord1').connect("clicked", self.import_picking_selection_data)
            self.builder.get_object('button_import_picking_selection_coord2').connect("clicked", self.import_picking_selection_data)
            
            self.builder.get_object('check_mass_restraints1').connect("clicked", self.toggle_mass_restraint1)
            self.builder.get_object('check_mass_restraints2').connect("clicked", self.toggle_mass_restraint2)
            
            self.builder.get_object('checkbox_reaction_coordinate2').connect("clicked", self.on_checkbox_reaction_coordinate2)
            self.builder.get_object('checkbox_geometry_optimization').connect("clicked", self.on_checkbox_geometry_optimization)
            
            self.builder.get_object('box_reaction_coordinate2').set_sensitive(False)
            self.builder.get_object('frame_geometry_optimization').set_sensitive(False)
            
            self.builder.get_object('button_run').connect("clicked", self.run)
            self.builder.get_object('button_cancel').connect('clicked', self.CloseWindow)
            #---------------------------------------------------------------------------------------------
            
            
            
            #self._starting_coordinates_model_update()
            self.window.show_all()
            
            self.input_type_combo.set_active(0)
            self.combobox_reaction_coord1.set_active(0)
            self.combobox_reaction_coord2.set_active(0)
            self.opt_methods_combo.set_active(0)
            self.md_integrators_combobox.set_active(0)

            self.Visible  = True   

        else:
            self.window.present()
            
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False


    def change_cb_coordType1 (self, combobox):
        """ Function doc """
        
        _type = self.combobox_reaction_coord1.get_active()        
        
        if _type == 0:
            self.builder.get_object('label_atom3_index_coord1').hide()
            self.builder.get_object('entry_atom3_index_coord1').hide()
            self.builder.get_object('label_atom3_name_coord1').hide()
            self.builder.get_object('entry_atom3_name_coord1').hide()
            
            self.builder.get_object('label_atom4_index_coord1').hide()
            self.builder.get_object('entry_atom4_index_coord1').hide()
            self.builder.get_object('label_atom4_name_coord1').hide()
            self.builder.get_object('entry_atom4_name_coord1').hide()

            self.builder.get_object('check_mass_restraints1').set_sensitive(False)

        if _type == 1:
            self.builder.get_object('label_atom3_index_coord1').show()
            self.builder.get_object('entry_atom3_index_coord1').show()
            self.builder.get_object('label_atom3_name_coord1') .show()
            self.builder.get_object('entry_atom3_name_coord1') .show()
            
            self.builder.get_object('label_atom4_index_coord1').hide()
            self.builder.get_object('entry_atom4_index_coord1').hide()
            self.builder.get_object('label_atom4_name_coord1') .hide()
            self.builder.get_object('entry_atom4_name_coord1') .hide()
            self.builder.get_object('check_mass_restraints1').set_sensitive(True)

        if _type == 2:
            self.builder.get_object('label_atom3_index_coord1').show()
            self.builder.get_object('entry_atom3_index_coord1').show()
            self.builder.get_object('label_atom3_name_coord1') .show()
            self.builder.get_object('entry_atom3_name_coord1') .show()
            
            self.builder.get_object('label_atom4_index_coord1').show()
            self.builder.get_object('entry_atom4_index_coord1').show()
            self.builder.get_object('label_atom4_name_coord1') .show()
            self.builder.get_object('entry_atom4_name_coord1') .show()
            self.builder.get_object('check_mass_restraints1').set_sensitive(False)

                    
    def change_cb_coordType2 (self, combobox):
        """ Function doc """
        
        _type = self.combobox_reaction_coord2.get_active()        
        
        if _type == 0:
            self.builder.get_object('label_atom3_index_coord2').hide()
            self.builder.get_object('entry_atom3_index_coord2').hide()
            self.builder.get_object('label_atom3_name_coord2').hide()
            self.builder.get_object('entry_atom3_name_coord2').hide()
            
            self.builder.get_object('label_atom4_index_coord2').hide()
            self.builder.get_object('entry_atom4_index_coord2').hide()
            self.builder.get_object('label_atom4_name_coord2').hide()
            self.builder.get_object('entry_atom4_name_coord2').hide()

            self.builder.get_object('check_mass_restraints1').set_sensitive(False)

        if _type == 1:
            self.builder.get_object('label_atom3_index_coord2').show()
            self.builder.get_object('entry_atom3_index_coord2').show()
            self.builder.get_object('label_atom3_name_coord2') .show()
            self.builder.get_object('entry_atom3_name_coord2') .show()
            
            self.builder.get_object('label_atom4_index_coord2').hide()
            self.builder.get_object('entry_atom4_index_coord2').hide()
            self.builder.get_object('label_atom4_name_coord2') .hide()
            self.builder.get_object('entry_atom4_name_coord2') .hide()
            self.builder.get_object('check_mass_restraints1').set_sensitive(True)

        if _type == 2:
            self.builder.get_object('label_atom3_index_coord2').show()
            self.builder.get_object('entry_atom3_index_coord2').show()
            self.builder.get_object('label_atom3_name_coord2') .show()
            self.builder.get_object('entry_atom3_name_coord2') .show()
            
            self.builder.get_object('label_atom4_index_coord2').show()
            self.builder.get_object('entry_atom4_index_coord2').show()
            self.builder.get_object('label_atom4_name_coord2') .show()
            self.builder.get_object('entry_atom4_name_coord2') .show()
            self.builder.get_object('check_mass_restraints2').set_sensitive(False)


    def on_combobox_inputtype (self, combobox):
        """ Function doc """
        _type =  combobox.get_active()  
        
        if _type == 0:
            self.builder.get_object('label_input_trajectory').hide()
            self.builder.get_object('folder_chooser_box').hide()
            self.builder.get_object('label_number_of_cpus').hide()
            self.builder.get_object('ncpus_spinbutton').hide()
            
            self.builder.get_object('label_starting_coordinates').show()
            self.builder.get_object('combobox_starting_coordinates').show()
        
        if _type == 1:
            self.builder.get_object('label_input_trajectory').show()
            self.builder.get_object('folder_chooser_box').show()
            self.builder.get_object('label_number_of_cpus').show()
            self.builder.get_object('ncpus_spinbutton').show()
            
            self.builder.get_object('label_starting_coordinates').hide()
            self.builder.get_object('combobox_starting_coordinates').hide()


    def on_md_integrator_combobox (self, widget = None):
        """ Function doc 
        
            0  :  "Velocity Verlet Dynamics", 
            1  :  "Leap Frog Dynamics",
            2  :  "Langevin Dynamics"
        
        """        
        if widget  == self.builder.get_object('md_integrator_combobox'):
            
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
                
                #self.temp_scale_options_combo.set_sensitive(True)
                
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
                
                #self.temp_scale_options_combo.set_active(0)
                #self.temp_scale_options_combo.set_sensitive(False)
                
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
            
                #self.temp_scale_options_combo.set_active(0)
                #self.temp_scale_options_combo.set_sensitive(False)
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


    def run (self, button):
        """ Function doc """
        
        '''this combobox has the reference to the starting coordinates of a simulation'''
        
        parameters = {
        
                      "simulation_type":"Umbralla_Sampling",
                      "trajectory_name": None                  , 
                      "dialog"         : False                 , 
                      
                      "folder"         :os.getcwd()            , 
                      
                      
                      
                      
                      #  - - - - - - - REACTION COORDINATES - - - - - - - - -
                      
                      "second_coordinate": False                           ,
                      "ATOMS_RC1":None                                     ,
                      "ATOMS_RC2":None                                     ,
                      
                      "ATOMS_RC1_NAMES":None                               ,
                      "ATOMS_RC2_NAMES":None                               ,                      
                      
                      "nsteps_RC1":0                                       ,
                      "nsteps_RC2":0                                       ,
                      
                      "force_constant_1":4000.0                            ,
                      "force_constant_2":4000.0                            ,
                      
                      "dincre_RC1":0.1                                     ,
                      "dincre_RC2":0.1                                     ,
                      "dminimum_RC1":0.0                                   ,
                      "dminimum_RC2":0.0                                   ,
                      
                      "sigma_pk1pk3_rc1":1.0                               ,
                      "sigma_pk3pk1_rc1":-1.0                              ,
                      "sigma_pk1pk3_rc2":1.0                               ,
                      "sigma_pk3pk1_rc2":-1.0                              ,
                      
                      "rc_type_1":"simple_distance"                        ,
                      "rc_type_2":"simple_distance"                        ,
                      
                      "MC_RC1":False                                       ,
                      "MC_RC2":False                                       ,
                      
                      
                      
                      
                      # - - - - - GEOMETRY OPT - - - - - - - -
                      "optimizer"      :"ConjugatedGradient"   ,
                      
                      "maxIterations"  : 600                    ,
                      "log_frequency"  : 10                     ,
                      "save_frequency" : 10                     ,
                      "rmsGradient"    : 0.1                    ,
                      #"save_format"    :None                   ,
                      #"save_traj"      :False                  ,
                      #"save_pdb"       :False                  }
                      
                      
                      # - - - - MOLECULAR DYNAMICS - - - - - -
                      'seed'                      : None         ,
                      'normalDeviateGenerator'    : None         ,
                      'steps'                     : None         ,
                      'timeStep'                  : None         ,
                      'trajectories'              : None         ,
                      
                      #VelocityVerletDynamics
                      'temperatureScaleFrequency' : None         ,
                      'temperatureScaleOption'    : None         ,  # "linear" , "constant" ,
                      'temperatureStart'          : None         ,
                      'temperatureStop'           : None         ,
                      
                      #LeapFrogDynamics
                      'pressure'                  : None         ,  #  LeapFrogDynamics 
                      'pressureControl'           : None         ,  #  LeapFrogDynamics 
                      'pressureCoupling'          : None         ,  #  LeapFrogDynamics 
                                                                       
                      'temperature'               : None         ,               
                      'temperatureControl'        : None         ,  # True / False LeapFrogDynamics / LangevinDynamics
                      'temperatureCoupling'       : None         ,  #  LeapFrogDynamics / LangevinDynamics
                      
                      #LangevinDynamics
                      'collisionFrequency'        : None         ,
        
        
        }
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        simParameters={ "simulation_type":"Geometry_Optimization",
                        "trajectory_name": None                  , 
                        "dialog"         : False                 , 
                        "folder"         :os.getcwd()            , 
                        "optimizer"      :"ConjugatedGradient"   ,
                        "maxIterations"  :600                    ,
                        "log_frequency"  :10                     ,
                        "save_frequency" :10                     ,
                        "rmsGradient"    :0.1                    ,
                        "save_format"    :None                   ,
                        "save_traj"      :False                  ,
                        "save_pdb"       :False                  }
        #----------------------------------------------------------------------------------
        #combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]
            vobject = self.main.vm_session.vm_objects_dic[vobject_id]
            
            '''This function imports the coordinates of a vobject into the dynamo system in memory.''' 
            #print('vobject:', vobject.name, len(vobject.frames) )
            self.main.p_session.get_coordinates_from_vobject_to_pDynamo_system(vobject)
                
        simParameters["optimizer"]      = self.opt_methods[self.builder.get_object('combobox_geo_opt').get_active()]
        simParameters["log_frequency"]  = int  ( self.builder.get_object('entry_log_frequence').get_text())
        simParameters["maxIterations"]  = int  ( self.builder.get_object('entry_max_int').get_text() )
        simParameters["rmsGradient"]    = float( self.builder.get_object('entry_rmsd_tol').get_text() )
        simParameters["vobject_name"]   = self.save_trajectory_box.builder.get_object('entry_trajectory_name').get_text()
        
        #------------------------------------------------------------------------------------
        if self.save_trajectory_box.builder.get_object('checkbox_save_traj').get_active():
            simParameters["save_traj"]       = True
            simParameters["dialog"]          = True
            simParameters["folder"]          = self.save_trajectory_box.folder_chooser_button.get_folder()
            simParameters["trajectory_name"] = self.save_trajectory_box.builder.get_object('entry_trajectory_name').get_text()
            saveFormat                       = self.save_trajectory_box.builder.get_object('combobox_format').get_active()
            simParameters["save_frequency"]  = int( self.save_trajectory_box.builder.get_object('entry_trajectory_frequency').get_text() ) 
            simParameters["trajectory_name"] = simParameters["trajectory_name"] + ".ptGeo"
            if   saveFormat == 0: simParameters["save_format"] = ".ptGeo"
            elif saveFormat == 1: simParameters["save_format"] = ".mdcrd"
            elif saveFormat == 2: simParameters["save_format"] = ".dcd"
            elif saveFormat == 3: simParameters["save_format"] = ".xyz"
            self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder = simParameters["folder"] 
        #-------------------------------------------------------------------------------------    
        #print(simParameters)
        
        isExist = os.path.exists(simParameters['folder'])
        if isExist:
            pass
        else:
            self.run_dialog()
            return None
        
        #self.main.refresh_main_statusbar(message = 'Running geometry optimization...')
        
        self.main.p_session.run_simulation( parameters = simParameters )
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

    
    def on_name_combo_changed(self, widget):
        """ Function doc """
    
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

    def update_working_folder_chooser (self, folder = None):
        """ Function doc """
        if folder:
            #print('update_working_folder_chooser')
            self.save_trajectory_box.set_folder(folder = folder)
        else:
            pass

   
    def on_checkbox_reaction_coordinate2 (self, widget):
        """ Function doc """
        if widget.get_active():
            self.builder.get_object('box_reaction_coordinate2').set_sensitive(True)
        else:
            self.builder.get_object('box_reaction_coordinate2').set_sensitive(False)


    def on_checkbox_geometry_optimization (self, widget):
        """ Function doc """
        if widget.get_active():
            self.builder.get_object('frame_geometry_optimization').set_sensitive(True)
        else:
            self.builder.get_object('frame_geometry_optimization').set_sensitive(False)

   


    def toggle_mass_restraint1 (self, widget):
        """ Function doc """
        self.refresh_dmininum(coord1 =  True)
    
    def toggle_mass_restraint2 (self, widget):
        """ Function doc """
        self.refresh_dmininum(coord2 =  True)

    def import_picking_selection_data (self, widget):
        """  
                   R                    R
                    \                  /
                     A1--A2  . . . . A3
                    /                  \ 
                   R                    R
                     ^   ^            ^
                     |   |            |
                    pk1-pk2  . . . . pk3
                       d1       d2	

                q1 =  1 / (mpk1 + mpk3)  =  [ mpk1 * r (pk3_pk2)  -   mpk3 * r (pk1_pk2) ]
        """       
        atom1 = self.vm_session.picking_selections.picking_selections_list[0]
        atom2 = self.vm_session.picking_selections.picking_selections_list[1]
        atom3 = self.vm_session.picking_selections.picking_selections_list[2]
        atom4 = self.vm_session.picking_selections.picking_selections_list[3]
        if atom1:
            self.vobject = atom1.vm_object
        else:
            return None
            
        if widget == self.builder.get_object('button_import_picking_selection_coord1'):
            if atom1:
                self.builder.get_object('entry_atom1_index_coord1').set_text(str(atom1.index-1) )
                self.builder.get_object('entry_atom1_name_coord1' ).set_text(str(atom1.name) )
            else: print('use picking selection to chose the central atom')            
            #-------
            if atom2:
                self.builder.get_object('entry_atom2_index_coord1').set_text(str(atom2.index-1) )
                self.builder.get_object('entry_atom2_name_coord1' ).set_text(str(atom2.name) )
            else: print('use picking selection to chose the central atom')
            #-------
            if atom3:
                self.builder.get_object('entry_atom3_index_coord1').set_text(str(atom3.index-1) )
                self.builder.get_object('entry_atom3_name_coord1' ).set_text(str(atom3.name) )
            else: print('use picking selection to chose the central atom')
             #-------
            if atom4:
                self.builder.get_object('entry_atom4_index_coord1').set_text(str(atom4.index-1) )
                self.builder.get_object('entry_atom4_name_coord1' ).set_text(str(atom4.name) )
            else: print('use picking selection to chose the central atom')
            
            self.refresh_dmininum( coord1 =  True, coord2 =  False)
            
        else:
            if atom1:
                self.builder.get_object('entry_atom1_index_coord2').set_text(str(atom1.index-1) )
                self.builder.get_object('entry_atom1_name_coord2' ).set_text(str(atom1.name) )
            else: print('use picking selection to chose the central atom')
            #-------
            if atom2:
                self.builder.get_object('entry_atom2_index_coord2').set_text(str(atom2.index-1) )
                self.builder.get_object('entry_atom2_name_coord2' ).set_text(str(atom2.name) )
            else: print('use picking selection to chose the central atom')            
            #-------
            if atom3:
                self.builder.get_object('entry_atom3_index_coord2').set_text(str(atom3.index-1) )
                self.builder.get_object('entry_atom3_name_coord2' ).set_text(str(atom3.name) )
            else: print('use picking selection to chose the central atom')           
            #-------
            if atom4:
                self.builder.get_object('entry_atom4_index_coord2').set_text(str(atom4.index-1) )
                self.builder.get_object('entry_atom4_name_coord2' ).set_text(str(atom4.name) )
            else: print('use picking selection to chose the central atom')
    
            self.refresh_dmininum( coord1 =  False, coord2 =  True )

    def refresh_dmininum (self, coord1 =  False, coord2 = False):
        """ Function doc """
        print('coord1',coord1 )
        print('coord2',coord2 )
        if coord1:
            _type = self.combobox_reaction_coord1.get_active()
            print('_type', _type)
            if _type == 0:
                index1 = int(self.builder.get_object('entry_atom1_index_coord1').get_text() )
                index2 = int(self.builder.get_object('entry_atom2_index_coord1').get_text() )

                dist1 = get_distance(self.vobject, index1, index2 )
                self.builder.get_object('entry_dmin_coord1').set_text(str(dist1))
            
            elif _type == 1:
                index1 = int(self.builder.get_object('entry_atom1_index_coord1').get_text() )
                index2 = int(self.builder.get_object('entry_atom2_index_coord1').get_text() )
                index3 = int(self.builder.get_object('entry_atom3_index_coord1').get_text() )
                
                dist1 = get_distance(self.vobject, index1, index2 )
                dist2 = get_distance(self.vobject, index2, index3 )
                
                if self.builder.get_object('check_mass_restraints1').get_active():
                    self.sigma_pk1_pk3, self.sigma_pk3_pk1  = compute_sigma_a1_a3(self.vobject, index1, index3)
                    #print('distance a1 - a2:', dist1 - dist2)
                    DMINIMUM =  (self.sigma_pk1_pk3 * dist1) -(self.sigma_pk3_pk1 * dist2*-1)
                    self.builder.get_object('entry_dmin_coord1').set_text(str(DMINIMUM))
                else:
                    DMINIMUM =  dist1- dist2
                    self.builder.get_object('entry_dmin_coord1').set_text(str(DMINIMUM))
        else:
            pass    
              
        if coord2:
            print('herehere')
            _type = self.combobox_reaction_coord2.get_active()
            try:
                if _type == 0:
                    index1 = int(self.builder.get_object('entry_atom1_index_coord2').get_text() )
                    index2 = int(self.builder.get_object('entry_atom2_index_coord2').get_text() )

                    dist1 = get_distance(self.vobject, index1, index2 )
                    self.builder.get_object('check_mass_restraints2').set_text(str(dist1))
                if _type == 1:
                    index1 = int(self.builder.get_object('entry_atom1_index_coord2').get_text() )
                    index2 = int(self.builder.get_object('entry_atom2_index_coord2').get_text() )
                    index3 = int(self.builder.get_object('entry_atom3_index_coord2').get_text() )
                    
                    dist1 = get_distance(self.vobject, index1, index2 )
                    dist2 = get_distance(self.vobject, index2, index3 )
                    
                    if self.builder.get_object('check_mass_restraints2').get_active():
                        self.sigma_pk1_pk3_rc2, self.sigma_pk3_pk1_rc2  = compute_sigma_a1_a3(self.vobject, index1, index3)
                        #print('distance a1 - a2:', dist1 - dist2)
                        DMINIMUM =  (self.sigma_pk1_pk3_rc2 * dist1) -(self.sigma_pk3_pk1_rc2 * dist2*-1)
                        self.builder.get_object('entry_dmin_coord2').set_text(str(DMINIMUM))
                    else:
                        DMINIMUM =  dist1- dist2
                        self.builder.get_object('entry_dmin_coord2').set_text(str(DMINIMUM))
            except:
                print(texto_d1)
                print(texto_d2d1)
        else:
            pass
    
