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

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')

##====================================================================================
class NormalModesWindow(Gtk.Window):
    """ Class doc """
    
    def __init__(self, main = None ):
        """ Class initialiser """
        self.main     = main
        self.home     = main.home 
        self.p_session= main.p_session 
        self.Visible  =  False        
        #self.residue_liststore = Gtk.ListStore(str, str, str)
        self.sym_tag  = 'normal_modes'
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
            self.builder.add_from_file(os.path.join(self.home,'gui/windows/normal_modes_window.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('normal_modes_window')
            self.window.set_title('Normal Modes Window')
            self.window.set_keep_above(True)
            
            
            '''--------------------------------------------------------------------------------------------'''
            
            # - - - - - - - - - - - - - Starting Coordinates ComboBox - - - - - - - - - - - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            self.combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
            #---------------------------------------------------------------------------------------------
            self._starting_coordinates_model_update(init = True)
            
            self.combobox_starting_coordinates.connect("changed", self.on_name_combo_changed)
            renderer_text = Gtk.CellRendererText()
            self.combobox_starting_coordinates.pack_start(renderer_text, True)
            self.combobox_starting_coordinates.add_attribute(renderer_text, "text", 0)
            #----------------------------------------------------------------------------------------------
            
            


            self.save_trajectory_box = SaveTrajectoryBox(parent = self.main, home = self.home)

            
            self.builder.get_object('parameters_box').pack_end(self.save_trajectory_box.box, True, True, 0)
                        
            # updating data 
            
            #------------------------------------------------------------------------------------------------
            if self.main.p_session.psystem[self.main.p_session.active_id]:
                if self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder == None:
                    folder = HOME
                else:
                    folder = self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder
                self.update_working_folder_chooser(folder = folder)            
            #------------------------------------------------------------------------------------------------

            if  self.p_session.psystem[self.p_session.active_id]:
                system = self.p_session.psystem[self.p_session.active_id]
                self.save_trajectory_box.set_filename (str(system.e_step_counter)+'-'+system.e_tag +'_'+ self.sym_tag )
            else:
                pass


            '''
            if 'tag' in self.main.p_session.psystem[self.main.p_session.active_id].keys():
                pass
            else:
                self.main.p_session.systems[self.main.p_session.active_id]['tag'] = 'molsys'
            
            tag  = self.main.p_session.systems[self.main.p_session.active_id]['tag']
            step = str(self.main.p_session.systems[self.main.p_session.active_id]['step_counter'])
            tag  = step+'_'+tag+'_geo_opt'  
            self.save_trajectory_box.builder.get_object('entry_trajectory_name').set_text(tag)
            '''
            #self._starting_coordinates_model_update()
            self.window.show_all()
            self.save_trajectory_box.builder.get_object('checkbox_save_traj').set_active(True)
            self.save_trajectory_box.builder.get_object('combobox_format').hide()
            self.save_trajectory_box.builder.get_object('label_format').hide()
            self.save_trajectory_box.builder.get_object('entry_trajectory_frequency').hide()
            self.save_trajectory_box.builder.get_object('label_trajectory_frequency').hide()
            
            self.Visible  = True   

    
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False


    def run(self, button):
        """ Function doc """
        
        '''this combobox has the reference to the starting coordinates of a simulation'''
        parameters={ "simulation_type"   :"Normal_Modes"  ,
                        "trajectory_name": None           , 
                        "dialog"         : False          , 
                        "folder"         : None           , 
                        "cycles"         : None           ,
                        "frames"         : None           ,
                        "temperature"    : None           ,
                        #"save_frequency" : None           ,
                        #"rmsGradient"    : None           ,
                        #"save_format"    : None           ,
                        "save_traj"      : False          ,
                        #"save_pdb"       : False          }
                        }
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
                
        parameters["cycles"]       = int( self.builder.get_object('entry_cycles').get_text())
        parameters["frames"]       = int( self.builder.get_object('entry_frames_per_cycles').get_text())
        parameters["temperature"]  = int( self.builder.get_object('entry_temperature').get_text() )
        
        #parameters["rmsGradient"]    = float( self.builder.get_object('entry_rmsd_tol').get_text() )
        parameters["vobject_name"]   = self.save_trajectory_box.builder.get_object('entry_trajectory_name').get_text()
        
        #------------------------------------------------------------------------------------
        if self.save_trajectory_box.builder.get_object('checkbox_save_traj').get_active():
            parameters["save_traj"]       = True
            parameters["dialog"]          = True
            parameters["folder"]          = self.save_trajectory_box.folder_chooser_button.get_folder()
            parameters["trajectory_name"] = self.save_trajectory_box.builder.get_object('entry_trajectory_name').get_text()

            #self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder = parameters["folder"] 
        #-------------------------------------------------------------------------------------    
        #print(parameters)
        
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

    def update_working_folder_chooser (self, folder = None):
        """ Function doc """
        if folder:
            #print('update_working_folder_chooser')
            self.save_trajectory_box.set_folder(folder = folder)
        else:
            pass
            #self.save_trajectory_box.set_folder(folder = self.main.p_session.systems[self.main.p_session.active_id]['working_folder'])
   
#=====================================================================================
   
