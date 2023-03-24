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

from gui.gtk_widgets import CoordinatesComboBox
from gui.gtk_widgets import SystemComboBox

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')

##====================================================================================
class WHAMWindow(Gtk.Window):
    """ Class doc """
    
    def __init__(self, main = None ):
        """ Class initialiser """
        self.main     = main
        self.home     = main.home 
        self.p_session= main.p_session 
        self.Visible  =  False        


    
    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = self.main.builder #Gtk.Builder()
            
            #self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/geometry_optimization_window.glade'))
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'gui/windows/WHAM_analysis_window.glade'))
            self.builder.connect_signals(self)

            self.window = self.builder.get_object('window')
            self.window.set_title('WHAM Analysis')
            self.window.set_keep_above(True)
            
            
            
            # - - - - - - - systems combobox - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            self.box = self.builder.get_object('box_system')
            self.combobox_systems = SystemComboBox(self.main)
            self.combobox_systems.connect("changed", self.on_combobox_systemsbox_changed)
            '''--------------------------------------------------------------------------------------------'''
            self.box.pack_start(self.combobox_systems, False, False, 0)
            '''--------------------------------------------------------------------------------------------'''

            
            '''--------------------------------------------------------------------------------------------'''
            self.method_store = Gtk.ListStore(str)
            methods = [
                        '1D',
                        '2D',
                        ]
            for method in methods:
                self.method_store.append([method])
                #print (method)
            self.methods_combo = self.builder.get_object('PMF_combobox')
            self.methods_combo.set_model(self.method_store)
            #self.methods_combo.connect("changed", self.on_name_combo_changed)
            self.methods_combo.set_model(self.method_store)
            
            renderer_text = Gtk.CellRendererText()
            self.methods_combo.pack_start(renderer_text, True)
            self.methods_combo.add_attribute(renderer_text, "text", 0)
            self.methods_combo.set_active(0)            
            
            self.builder.get_object('button_add_trajectory_blocks').connect("clicked", self.button_add_trajectories)
            
            self.window.connect("destroy", self.CloseWindow)
            self.window.show_all()
            self.Visible  = True   
        
        else:
            self.window.present()
    
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False


    def run (self, button):
        """ Function doc """
        
        
        e_id = self.combobox_systems.get_system_id()
        system = self.p_session.psystem[e_id]
        
        
        _type     = self.methods_combo.get_active()

        bins      = self.builder.get_object('entry_bins'     ).get_text()
        frequency = self.builder.get_object('entry_frequency').get_text()
        max_int   = self.builder.get_object('entry_max_int'  ).get_text()
        RMS_grad  = self.builder.get_object('entry_RMS_grad' ).get_text()
        temp      = self.builder.get_object('entry_temp'     ).get_text()
        
        
        
        #fileNames = glob.glob ( os.path.join ( outPath, "bAla_phi_*_psi_*.ptRes" ) )
        fileNames = None
        fileNames.sort ( )
        
        parameters = {
                      'system'               : system    ,
                      'fileNames'            : fileNames ,
                      'type'                 : _type     ,
                      'logFrequency'         : frequency ,
                      'maximumIterations'    : max_int   ,
                      'rmsGradientTolerance' : RMS_grad  ,
                      'temperature'          : temp       
                      }


        
        if _type ==1:
            parameters['bins'] = [bins, bins] 

        else:
            parameters['bins'] = [bins] 

        
        print(parameters)
        
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

    
    def button_add_trajectories (self, button):
        """ Function doc """
        dialog = Gtk.FileChooserDialog("Select Folders", None,
                                       Gtk.FileChooserAction.SELECT_FOLDER,
                                       (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                                        Gtk.STOCK_OPEN, Gtk.ResponseType.OK))
        dialog.set_select_multiple(True)
        
        parameters = {}
        
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            # Get the list of selected folder paths
            folder_paths = dialog.get_filenames()
            print(folder_paths)
            
            
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
        if self.Visible:
            self.update_working_folder_chooser()
            
        
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

#=====================================================================================
    def on_combobox_systemsbox_changed(self, widget):
        """ Function doc """
        system_id = self.combobox_systems.get_system_id()
       
            
        #self.update_window ( selections = False, restraints = True)
