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
from gui.widgets.custom_widgets import FolderChooserButton
from gui.widgets.custom_widgets import SaveTrajectoryBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.widgets.custom_widgets import ReactionCoordinateBox
from gui.windows.setup.windows_and_dialogs import ExportScriptDialog

from gui.windows.simulation.PES_scan_window            import compute_sigma_a1_a3 
from util.geometric_analysis            import get_dihedral 
from util.geometric_analysis            import get_angle 
from util.geometric_analysis            import get_distance 

#from gui.windows.PES_scan_window            import texto_d1 
#from gui.windows.PES_scan_window            import texto_d2d1 
#from util.periodic_table import atomic_dic 
from pprint import pprint

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
        self.sym_tag    = "umb_sam"
        self.input_types  = {
                            0 : 'From Coordinates (sequential)',
                            1 : 'From Trajectory (parallel)'   ,
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
            self.builder.add_from_file(os.path.join (self.home,'src/gui/windows/simulation/umbrella_sampling_window.glade') )
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
            
            

            #----------------------------------------------------------------------------------------------
            # - - - - - - - - - - - - - Starting Coordinates ComboBox - - - - - - - - - - - - - - - - -
            self.box_coordinates = self.builder.get_object('box_coordinates')
            self.combobox_starting_coordinates = CoordinatesComboBox() #self.builder.get_object('coordinates_combobox')
            #self.combobox_starting_coordinates.connect("changed", self.on_name_combo_changed)
            self.box_coordinates.pack_start(self.combobox_starting_coordinates, False, False, 0)
            self._starting_coordinates_model_update(init = True)
            #----------------------------------------------------------------------------------------------
            
            
            
            #----------------------------------------------------------------------------------------------
            # - - - - - - - - - - - - - Folder Chooser Button - - - - - - - - - - - - - - - - -
            #----------------------------------------------------------------------------------------------
            self.folder_chooser_button = FolderChooserButton(main =  self.main, sel_type = 'folder', home =  self.home)
            self.builder.get_object('folder_chooser_box').pack_start(self.folder_chooser_button.btn, True, True, 0)
            system_id      = self.p_session.active_id
            #----------------------------------------------------------------------------------------------
            #----------------------------------------------------------------------------------------------
            self.folder_chooser_button2 = FolderChooserButton(main =  self.main, sel_type = 'folder', home =  self.home)
            self.builder.get_object('folder_chooser_box2').pack_start(self.folder_chooser_button2.btn, True, True, 0)
            #----------------------------------------------------------------------------------------------
            
        
        
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
            
            
            
            #----------------------------------------------------------------------------------------------
            # - - - - - - - - - - - - - Reaction Coordinates 1 ComboBox - - - - - - - - - - - - - - - - -
            #----------------------------------------------------------------------------------------------
            self.RC_box1 = ReactionCoordinateBox(main = self.main, mode = 1)
            self.builder.get_object('rc1_aligment').add(self.RC_box1)
            self.RC_box2 = ReactionCoordinateBox(main = self.main, mode = 1)
            self.builder.get_object('rc2_aligment').add(self.RC_box2)
            #----------------------------------------------------------------------------------------------

        
        
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
            self.builder.get_object('checkbox_reaction_coordinate2').connect("clicked", self.on_checkbox_reaction_coordinate2)
            self.builder.get_object('checkbox_geometry_optimization').connect("clicked", self.on_checkbox_geometry_optimization)
            
            self.RC_box2.set_sensitive(False)
            self.builder.get_object('frame_geometry_optimization').set_sensitive(False)
            
            self.builder.get_object('button_run').connect("clicked", self.run)
            self.builder.get_object('button_cancel').connect('clicked', self.CloseWindow)
            self.builder.get_object('button_export').connect('clicked', self.on_btn_export)
            #---------------------------------------------------------------------------------------------
            
            
            
            #self._starting_coordinates_model_update()
            self.window.show_all()
            
            self.input_type_combo.set_active(0)
            self.RC_box1.set_rc_type(0)
            self.RC_box2.set_rc_type(0)
            self.opt_methods_combo.set_active(0)
            self.md_integrators_combobox.set_active(0)
            self.update_working_folder_chooser ( )
            
            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
                self.builder.get_object('entry_traj_name').set_text(output_name)
            else:
                pass
            
            
            self.Visible  = True   

        else:
            self.window.present()
            

    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False


    def on_combobox_inputtype (self, combobox):
        """ Function doc """
        _type =  combobox.get_active()  
        
        if _type == 0:
            self.builder.get_object('label_input_trajectory').hide()
            self.builder.get_object('folder_chooser_box').hide()
            self.builder.get_object('label_number_of_cpus').hide()
            self.builder.get_object('ncpus_spinbutton').hide()
            
            self.builder.get_object('label_starting_coordinates').show()
            self.combobox_starting_coordinates.show()
        
            self.RC_box1.entry_step_size .set_sensitive(True)
            self.RC_box1.entry_nsteps    .set_sensitive(True)
            self.RC_box1.entry_dmin_coord.set_sensitive(True)
            self.RC_box2.entry_step_size .set_sensitive(True)
            self.RC_box2.entry_nsteps    .set_sensitive(True)
            self.RC_box2.entry_dmin_coord.set_sensitive(True)
        
            self.RC_box1.label_step_size     .set_sensitive(True)
            self.RC_box1.label_nsteps        .set_sensitive(True)        
            #self.RC_box1.label_force_constant.set_sensitive(True)        
            self.RC_box1.label_dmin          .set_sensitive(True)        
            
            self.RC_box2.label_step_size     .set_sensitive(True)
            self.RC_box2.label_nsteps        .set_sensitive(True)        
            #self.RC_box2.label_force_constant.set_sensitive(True)        
            self.RC_box2.label_dmin          .set_sensitive(True)        
        
        
        
        if _type == 1:
            self.builder.get_object('label_input_trajectory').show()
            self.builder.get_object('folder_chooser_box').show()
            self.builder.get_object('label_number_of_cpus').show()
            self.builder.get_object('ncpus_spinbutton').show()
            
            self.builder.get_object('label_starting_coordinates').hide()
            self.combobox_starting_coordinates.hide()
            
            self.RC_box1.entry_step_size .set_sensitive(False)
            self.RC_box1.entry_nsteps    .set_sensitive(False)
            self.RC_box1.entry_dmin_coord.set_sensitive(False)
            self.RC_box2.entry_step_size .set_sensitive(False)
            self.RC_box2.entry_nsteps    .set_sensitive(False)
            self.RC_box2.entry_dmin_coord.set_sensitive(False)


            self.RC_box1.label_step_size.set_sensitive(False)
            self.RC_box1.label_nsteps   .set_sensitive(False)        
            self.RC_box1.label_dmin     .set_sensitive(False)        
            
            self.RC_box2.label_step_size.set_sensitive(False)
            self.RC_box2.label_nsteps   .set_sensitive(False)        
            self.RC_box2.label_dmin     .set_sensitive(False)


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
                self.builder.get_object('entry_traj_name').set_text(output_name)
            else:
                pass


    def update_working_folder_chooser (self, folder = None):
        """ Function doc """
        if folder:
            #print('update_working_folder_chooser')
            self.folder_chooser_button.set_folder(folder = folder)
        else:
            
            try:
                folder = self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder
                if folder:
                    self.folder_chooser_button2.set_folder(folder = folder)
                else:
                    pass
            except:
                self.folder_chooser_button2.set_folder(folder = self.home )
   

    def on_checkbox_reaction_coordinate2 (self, widget):
        """ Function doc """
        if widget.get_active():
            self.RC_box2.set_sensitive(True)
        else:
            self.RC_box2.set_sensitive(False) 


    def on_checkbox_geometry_optimization (self, widget):
        """ Function doc """
        if widget.get_active():
            self.builder.get_object('frame_geometry_optimization').set_sensitive(True)
        else:
            self.builder.get_object('frame_geometry_optimization').set_sensitive(False)

   
    def get_input_setup_box_info (self, parameters):
        """ Function doc """
        parameters['input_type'] = self.builder.get_object('combobox_input_type').get_active()
        
        if parameters['input_type'] == 0:
            vobject_id = self.combobox_starting_coordinates.get_vobject_id()
            vobject = self.main.vm_session.vm_objects_dic[vobject_id]
            self.main.p_session.get_coordinates_from_vobject_to_pDynamo_system(vobject)
            parameters['vobj'] = vobject.name
            parameters['source_folder'] = None
            
        elif parameters['input_type'] == 1:
            parameters['vobj'] = None
            parameters['source_folder'] = self.folder_chooser_button.get_folder ()
            
        
        return parameters
    
        
    def get_parameters (self, ):
        
        """ Function doc """
        
        '''this combobox has the reference to the starting coordinates of a simulation'''
        parameters = {"simulation_type" : "Umbrella_Sampling"}

        parameters = self.get_input_setup_box_info ( parameters)
        parameters["folder"] = self.folder_chooser_button2.get_folder()        


        #vobject_id = self.combobox_starting_coordinates.get_vobject_id()
        #vobject = self.main.vm_session.vm_objects_dic[vobject_id]
        #self.main.p_session.get_coordinates_from_vobject_to_pDynamo_system(vobject)

        '''
        
        - If a trajectory folder is given then parallelization is allowed 
            - NmaxThreads >= 1
        
        - If a vobject is provided and a second reaction coordinate is active 
          (2D), then parallelization is allowed
            - NmaxThreads >= 1
        
        - Else: sequetial
        
        '''
        parameters["RC1"] = self.RC_box1.get_rc_data()
        if self.builder.get_object('checkbox_reaction_coordinate2').get_active():
            parameters["RC2"] = self.RC_box2.get_rc_data()
            parameters["NmaxThreads"] =  int(self.builder.get_object('ncpus_spinbutton').get_value())
        else:
            parameters["RC2"] = None
            
            if parameters['input_type'] == 1:
                parameters["NmaxThreads"] = int(self.builder.get_object('ncpus_spinbutton').get_value())
            else:
                parameters["NmaxThreads"] = 1
            
        
        
        #----------------------------------------------------------------------
        #                            GEO OPT
        #----------------------------------------------------------------------
        if self.builder.get_object('checkbox_geometry_optimization').get_active():
            parameters['OPT_parm'] = self._define_OPT_parameters ()
        else:
            parameters['OPT_parm'] = None
        #----------------------------------------------------------------------
        
        
        #----------------------------------------------------------------------
        #                          M. DYNAMICS
        #----------------------------------------------------------------------        
        parameters['MD_parm'] = self._define_MD_parameters ( )
        #----------------------------------------------------------------------
        
        
        parameters['traj_folder_name'] = self.builder.get_object('entry_traj_name').get_text()
        
        return parameters
    



        
    def run (self, button):
        """ Function doc """
        
        parameters = self.get_parameters()
        
        #pprint(parameters)

        isExist = os.path.exists(parameters['folder'])
        if isExist:
            pass
        else:
            self.run_dialog()
            return None
        
        
        self.main.p_session.run_simulation( parameters = parameters )
        self.window.destroy()
        self.Visible    =  False


    def on_btn_export (self, widget):
        """ Function doc """

        #header += '\nparameters = '
        parameters = self.get_parameters()
        try:
            parameters.pop('vobject_name')
        except:
            pass
        #pprint.pprint(parameters)

        export_dialog =  ExportScriptDialog(self.main, parameters = parameters)
        

        


    def _define_OPT_parameters (self):
        parameters = {}
        opt_id       = self.builder.get_object('combobox_geo_opt').get_active()
        parameters['optimizer']        =  self.opt_methods[opt_id]
        parameters["maxIterations"]    =  int  (self.builder.get_object('entry_max_int') .get_text())
        parameters["rmsGradient"]      =  float(self.builder.get_object('entry_rmsd_tol').get_text())
        return parameters
    
    
    def _define_MD_parameters (self):
        """ Function doc """
        integrator_id       = self.builder.get_object('md_integrator_combobox').get_active()
        number_of_steps_eq  = int(self.builder.get_object('entry_number_of_steps_eq').get_text())
        number_of_steps_dc  = int(self.builder.get_object('entry_number_of_steps_dc').get_text())
        
        
        temp_start          = float(self.builder.get_object('entry_temp_start').get_text())
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
        

        
        MD_method = {
                     0 : "Verlet"   ,
                     1 : "LeapFrog" ,
                     2 : "Langevin" ,
                     }
        
        

        parameters = {
                    "simulation_type"           : "Umbrella_Sampling",
                    "folder"                    : HOME                    ,
                    'integrator'                : MD_method[integrator_id], # verlet / leapfrog /langevin
                    'logFrequency'              : log_frequence           ,
                    'seed'                      : random_seed             ,
                    'normalDeviateGenerator'    : None                    ,
                    'steps_eq'                  : number_of_steps_eq      ,
                    'steps_dc'                  : number_of_steps_dc      ,
                    'timeStep'                  : time_step               ,
                    #'trajectories'              : None                    ,
                    'trajectory_frequency'      : int(self.builder.get_object('entry_traj_frequency').get_text()), 
                    'trajectory_frequency_dc_ptRes'   : int(self.builder.get_object('entry_traj_frequency_dc1').get_text()), 
                    'trajectory_frequency_dc_ptGeo'   : int(self.builder.get_object('entry_traj_frequency_dc2').get_text()), 
                    
                    #VelocityVerletDynamics
                    'temperatureScaleFrequency' : temp_scale_factor                ,
                    'temperatureStart'          : temp_start                       ,
                    
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
        return parameters
