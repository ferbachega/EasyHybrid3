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
import copy
from gui.widgets.custom_widgets import FolderChooserButton
from gui.widgets.custom_widgets import SaveTrajectoryBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.windows.setup.windows_and_dialogs import ExportScriptDialog

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')

##====================================================================================
class ChainOfStatesOptWindow(Gtk.Window):
    """ Class doc """
    
    def __init__(self, main = None ):
        """ Class initialiser """
        self.main     = main
        self.home     = main.home 
        self.p_session= main.p_session 
        self.Visible  =  False        
        #self.residue_liststore = Gtk.ListStore(str, str, str)
        self.sym_tag  = 'chain_of_state_opt'
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
            self.builder.add_from_file(os.path.join(self.home,'src/gui/windows/simulation/chain_of_states_opt_window.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('chain_of_states_opt_window')
            self.window.set_title('Chain of State Optmization Window')
            #self.window.set_default_size(500, 600)
            self.window.set_keep_above(True)
            #'''--------------------------------------------------------------------------------------------'''



            # - - - - - - - - - - - - - Starting Coordinates ComboBox 1 - - - - - - - - - - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            #self.combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates1')
            ##---------------------------------------------------------------------------------------------
            #
            #self.combobox_starting_coordinates.connect("changed", self.on_name_combo_changed)
            #renderer_text = Gtk.CellRendererText()
            #self.combobox_starting_coordinates.pack_start(renderer_text, True)
            #self.combobox_starting_coordinates.add_attribute(renderer_text, "text", 0)
            #----------------------------------------------------------------------------------------------
            
            #----------------------------------------------------------------------------------------------
            self.box_coordinates1 = self.builder.get_object('box_coordinates1')
            self.combobox_starting_coordinates = CoordinatesComboBox() #self.builder.get_object('coordinates_combobox')
            self.box_coordinates1.pack_start(self.combobox_starting_coordinates, False, False, 0)
            #----------------------------------------------------------------------------------------------
            
            

            # - - - - - - - - - - - - - Starting Coordinates ComboBox 2 - - - - - - - - - - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            #self.combobox_starting_coordinates2 = self.builder.get_object('combobox_starting_coordinates2')
            ##---------------------------------------------------------------------------------------------
            #
            #self.combobox_starting_coordinates2.connect("changed", self.on_name_combo_changed)
            #renderer_text = Gtk.CellRendererText()
            #self.combobox_starting_coordinates2.pack_start(renderer_text, True)
            #self.combobox_starting_coordinates2.add_attribute(renderer_text, "text", 0)
            ##----------------------------------------------------------------------------------------------
            
            #----------------------------------------------------------------------------------------------
            self.box_coordinates2 = self.builder.get_object('box_coordinates2')
            self.combobox_starting_coordinates2 = CoordinatesComboBox() #self.builder.get_object('coordinates_combobox')
            self.box_coordinates2.pack_start(self.combobox_starting_coordinates2, False, False, 0)
            #----------------------------------------------------------------------------------------------
            
            
            
            
            
            
            self._starting_coordinates_model_update(init = True)

            
            
            
            '''--------------------------------------------------------------------------------------------'''     
            self.folder_chooser_button = FolderChooserButton(main =  self.main, sel_type = 'folder', home =  self.home)
            self.builder.get_object('folder_chooser_box').pack_start(self.folder_chooser_button.btn, True, True, 0)
            system_id      = self.p_session.active_id
            
            #------------------------------------------------------------------------------------------------
            if self.main.p_session.psystem[self.p_session.active_id]:
                if self.main.p_session.psystem[self.p_session.active_id].e_working_folder == None:
                    folder = HOME
                else:
                    folder = self.main.p_session.psystem[self.p_session.active_id].e_working_folder
                self.folder_chooser_button.set_folder(folder = folder)
            #------------------------------------------------------------------------------------------------

            #------------------------------------------------------------------------------------------------
            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
                self.builder.get_object('traj_name').set_text(output_name)            
            else:
                pass
            #------------------------------------------------------------------------------------------------
            
            self.builder.get_object('button_cancel').connect('clicked', self.CloseWindow)
            self.builder.get_object('button_run').connect('clicked', self.run)
            self.builder.get_object('button_export').connect('clicked', self.on_btn_export)
            self.window.set_size_request(-1, -1)
            self.window.show_all()
            
            self.builder.get_object('label_poolfactory').hide()
            self.builder.get_object('entry_poolfactory').hide()
            
            
            self.Visible  = True   
            

        else:
            self.window.present()
       
            
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False


    def get_parameters (self):
        """ Function doc """
        '''this combobox has the reference to the starting coordinates of a simulation'''
        
        parameters={       "simulation_type"                            :"Nudged_Elastic_Band"  ,
                           "number_of_structures"                       : 11                    ,
                           "logFrequency"                               : 1                     ,
                           "maximumIterations"                          : 600                   ,
                           "rmsGradientTolerance"                       : 0.1                   ,
                           "trajectory_name"                            : 'new_trajectory'       ,
                           "folder"                                     : os.getcwd()            ,
                           
                           #pdynamo default parameters
                           "fixedTerminalImages"                        : True                                ,
                           "forceOneSingleImageOptimization"            : False                               ,
                           "forceSingleImageOptimizations"              : False                               ,
                           "forceSplineRedistributionCheckPerIteration" : False                               ,
                           "freezeRMSGradientTolerance"                 :   0.0                               ,
                           "optimizer"                                  : None                                ,
                           "poolFactory"                                : None                                ,
                           "rmsGradientToleranceScale"                  :   0.25                              ,
                           "splineRedistributionTolerance"              :   1.5                               ,
                           "springForceConstant"                        : 500.0                               ,
                           "useSplineRedistribution"                    : False                               }
        
        
        
        '''        
        #parameters={    #"simulation_type"                 :"Nudged_Elastic_Band"  ,
        #                #"number_of_structures"            : 11                    , 
        #                #"logFrequency"                   : 1                     , 
        #                #"maximumIterations"               : 600                   ,
        #                #"rmsGradientTolerance"            : 0.1                   ,
        #                #"spring_force_constant"           : 500                   ,
        #                #"spline_redistribution_tolerance" : 1.5                   ,
        #                
        #                #"fixed_terminal_images"                           : True  ,
        #                #"force_spline_redistribution_check_per_iteration" : False ,           
        #                
        #                #"trajectory_name"        : 'new_trajectory'       ,
        #                #"folder"                 : os.getcwd()            , 
        #             }
        
        

        #----------------------------------------------------------------------------------
        #parameters['number_of_structures'           ] = int  (self.builder.get_object('entry_mun_of_structures').get_text()              )
        #parameters["logFrequency"                  ] = int  (self.builder.get_object('entry_log_frequency').get_text()                  )
        #parameters["maximumIterations"              ] = int  (self.builder.get_object('entry_max_int').get_text()                        )
        #parameters["rmsGradientTolerance"           ] = float(self.builder.get_object('entry_rmsd_tol').get_text()                       )
        #parameters["spring_force_constant"          ] = int  (self.builder.get_object('entry_spring_force_constant').get_text()          )
        #parameters["spline_redistribution_tolerance"] = float(self.builder.get_object('entry_spline_redistribution_tolerance').get_text())
        #
        #parameters["fixed_terminal_images"                          ] = self.builder.get_object('check_fixed_terminal_images').get_active()
        #parameters["force_spline_redistribution_check_per_iteration"] = self.builder.get_object('check_force_spline_redistribution_check_per_iteration').get_active()
        #
        #parameters["trajectory_name"] = self.builder.get_object('traj_name').get_text()
        #parameters["folder"         ] = self.folder_chooser_button.get_folder()
        '''
        parameters['number_of_structures'           ] = int  (self.builder.get_object('entry_mun_of_structures').get_text()              )  
        parameters["logFrequency"                  ] = int  (self.builder.get_object('entry_log_frequency').get_text()                  )  
        parameters["maximumIterations"              ] = int  (self.builder.get_object('entry_max_int').get_text()                        )  
        parameters["rmsGradientTolerance"           ] = float(self.builder.get_object('entry_rmsd_tol').get_text()                       )
        
        parameters["fixedTerminalImages"            ] = self.builder.get_object('check_fixed_terminal_images').get_active()
        
        parameters["forceOneSingleImageOptimization"           ] = self.builder.get_object('check_force_one_single_image_optimization').get_active() # True/Fase
        parameters["forceSingleImageOptimizations"             ] = self.builder.get_object('check_force_single_image_optimizations').get_active() # True/False
        parameters["forceSplineRedistributionCheckPerIteration"] = self.builder.get_object('check_force_spline_redistribution_check_per_iteration').get_active()# True/False
        
        parameters["freezeRMSGradientTolerance"     ] = float(self.builder.get_object('entry_freeze_RMS_gradient_tolerance').get_text() )# #0.0
        
        parameters["optimizer"                      ] =  None
        parameters["poolFactory"                    ] =  int(self.builder.get_object('entry_poolfactory').get_text())
        parameters["rmsGradientToleranceScale"      ] =  float(self.builder.get_object('entry_optimizer_tolerance_scaling').get_text())#0.25
        
        parameters["splineRedistributionTolerance"  ] = float(self.builder.get_object('entry_spline_redistribution_tolerance').get_text())
        parameters["springForceConstant"            ] = float(self.builder.get_object('entry_spring_force_constant').get_text()          )
        parameters["useSplineRedistribution"        ] = self.builder.get_object('check_use_spline_redistribution').get_active() # True/False
        

        parameters["trajectory_name"] = self.builder.get_object('traj_name').get_text()
        parameters["folder"         ] = self.folder_chooser_button.get_folder()
        
        
        '''---------------------------------------------------------------------------------'''
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]
            vobject = self.main.vm_session.vm_objects_dic[vobject_id]
            
            '''This function imports the coordinates of a vobject into the dynamo system in memory.''' 
            self.main.p_session.get_coordinates_from_vobject_to_pDynamo_system(vobject)
            parameters["reac_coordinates"] = copy.deepcopy(self.main.p_session.psystem[self.main.p_session.active_id].coordinates3)
            #print (list(parameters["reac_coordinates"]))
        '''---------------------------------------------------------------------------------'''
        
        
        '''---------------------------------------------------------------------------------'''
        tree_iter = self.combobox_starting_coordinates2.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.combobox_starting_coordinates2.get_model()
            name, vobject_id = model[tree_iter][:2]
            vobject = self.main.vm_session.vm_objects_dic[vobject_id]
            
            '''This function imports the coordinates of a vobject into the dynamo system in memory.''' 
            self.main.p_session.get_coordinates_from_vobject_to_pDynamo_system(vobject)
            parameters["prod_coordinates"] = copy.deepcopy(self.main.p_session.psystem[self.main.p_session.active_id].coordinates3)
            #print (list(parameters["prod_coordinates"]))
        '''---------------------------------------------------------------------------------'''
        return parameters

    def run(self, button):
        """ Function doc """
        parameters = self.get_parameters()
        parameters['external_coords'] = False
        
        isExist = os.path.exists(parameters['folder'])
        if isExist:
            pass
        else:
            self.run_dialog()
            return None
        
        #self.main.refresh_main_statusbar(message = 'Running geometry optimization...')
        
        self.main.p_session.run_simulation( parameters = parameters )
        self.window.destroy()
        self.Visible    =  False

    
    def on_btn_export (self, parameters):
        """ Function doc """
        parameters = self.get_parameters()
        parameters['external_coords'] = True
        try:
            parameters.pop('vobject_name')
        except:
            pass
        #pprint.pprint(parameters)
        export_dialog =  ExportScriptDialog(self.main, parameters = parameters)
        
        
    def _starting_coordinates_model_update (self, init = False):
        """ Function doc """
        #------------------------------------------------------------------------------------
        '''The combobox accesses, according to the id of the active system, 
        listostore of the dictionary object_list state_dict'''
        if self.Visible:

            e_id = self.main.p_session.active_id 
            self.combobox_starting_coordinates.set_model(self.main.vobject_liststore_dict[e_id])
            self.combobox_starting_coordinates2.set_model(self.main.vobject_liststore_dict[e_id])
            #------------------------------------------------------------------------------------
            size = len(self.main.vobject_liststore_dict[e_id])
            self.combobox_starting_coordinates.set_active(size-1)
            #------------------------------------------------------------------------------------
        else:
            if init:
                e_id = self.main.p_session.active_id 
                self.combobox_starting_coordinates.set_model(self.main.vobject_liststore_dict[e_id])
                self.combobox_starting_coordinates2.set_model(self.main.vobject_liststore_dict[e_id])
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
        
        
            #------------------------------------------------------------------------------------------------
            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
                self.builder.get_object('traj_name').set_text(output_name)            
            else:
                pass
            #------------------------------------------------------------------------------------------------

    def update_working_folder_chooser (self, folder = None):
        """ Function doc """
        if folder:
            #print('update_working_folder_chooser')
            self.folder_chooser_button.set_folder(folder = folder)
        else:
            
            folder = self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder
            if folder:
                self.folder_chooser_button.set_folder(folder = folder)
            else:
                pass
   
#=====================================================================================
   
