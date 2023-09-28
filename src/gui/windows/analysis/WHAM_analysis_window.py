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
from gui.widgets.custom_widgets import SystemComboBox

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
        self.liststore= Gtk.ListStore(bool, str)

    
    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = self.main.builder #Gtk.Builder()
            
            #self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/geometry_optimization_window.glade'))
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'src/gui/windows/analysis/WHAM_analysis_window.glade'))
            self.builder.connect_signals(self)

            self.window = self.builder.get_object('window')
            self.window.set_title('WHAM Analysis')
            #self.window.set_keep_above(True)
            
            
            
            
            self.treeview = self.builder.get_object('treeview_trajectory_blocks')
            self.treeview.set_model(self.liststore)
            

            renderer_radio = Gtk.CellRendererToggle()
            renderer_radio.connect("toggled", self.on_cell_active_toggled)
            column_radio = Gtk.TreeViewColumn("", renderer_radio, active=0 )
            self.treeview.append_column(column_radio)
                
            
            
            
            renderer_pixbuf = Gtk.CellRendererPixbuf()
            renderer_text = Gtk.CellRendererText()

            column_text = Gtk.TreeViewColumn("Block")#, renderer_text, text=2)
            column_text.pack_start(renderer_text, True)
            column_text.add_attribute(renderer_text, "text", 1)
            self.treeview.append_column(column_text)
            
            
            
            
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
            
            
            
            
            '''--------------------------------------------------------------------------------------------'''     
            self.folder_chooser_button = FolderChooserButton(main =  self.main, sel_type = 'folder', home =  self.home)
            self.builder.get_object('folder_chooser_box').pack_start(self.folder_chooser_button.btn, True, True, 0)
            
            
            
            
            
            
            self.builder.get_object('button_add_trajectory_blocks').connect("clicked", self.button_add_trajectories)
            
            self.button_cancel = self.builder.get_object('button_cancel')
            self.button_cancel.connect("clicked", self.clear_treeview)
            self.window.connect("destroy", self.CloseWindow)
            self.window.show_all()
            
            self.builder.get_object('entry_bins_RC2').hide()
            self.builder.get_object('label_bins_RC2').hide()

            
            self.Visible  = True   
        
        else:
            self.window.present()
    
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False


    def run (self, button):
        """ Function doc """
        
        try:
            e_id   = self.combobox_systems.get_system_id()
            system = self.p_session.psystem[e_id]
        except:
            system = None
        
        _type     = self.methods_combo.get_active()

        bins      = self.builder.get_object('entry_bins'     ).get_text()
        frequency = self.builder.get_object('entry_frequency').get_text()
        max_int   = self.builder.get_object('entry_max_int'  ).get_text()
        RMS_grad  = self.builder.get_object('entry_RMS_grad' ).get_text()
        temp      = self.builder.get_object('entry_temp'     ).get_text()
        
        
        
        #fileNames = glob.glob ( os.path.join ( outPath, "bAla_phi_*_psi_*.ptRes" ) )
        fileNames = []
        for item in self.liststore:
            if item[0]:
                fileNames.append(item[1])
            else:
                pass
                
        fileNames.sort ( )
        
        parameters = {'analysis_type'        : 'wham'    , 
                      'system'               : system    ,
                      'fileNames'            : fileNames ,
                      'type'                 : _type     ,
                      'logFrequency'         : int(frequency) ,
                      'maximumIterations'    : int(max_int  ) ,
                      'rmsGradientTolerance' : float(RMS_grad ) ,
                      'temperature'          : int(temp     )  
                      }


        
        if _type ==1:
            bins_RC1 = self.builder.get_object('entry_bins'     ).get_text()
            bins_RC2 = self.builder.get_object('entry_bins_RC2').get_text()
            parameters['bins'] = [int(bins_RC1), int(bins_RC2)] 

        else:
            parameters['bins'] = [int(bins)] 
        
        parameters["folder"]  = self.folder_chooser_button.get_folder()
        
        
        #----------------------------------------------------------------------
        if parameters["folder"]:
            pass
        else:
            print('folder: ', parameters["folder"] )
        #----------------------------------------------------------------------
        
        #----------------------------------------------------------------------
        parameters["logfile"] = self.builder.get_object('log_name').get_text() 
        #----------------------------------------------------------------------
        
        
        #print(parameters)
        TrueFalse = self.p_session.run_analysis(parameters)
        
        if TrueFalse:
            self.window.destroy()
            self.Visible    =  False

        else:
            pass
            self.window.present()

    def on_cell_active_toggled (self, widget, path):
        """ Function doc """
        self.liststore[path][0] = not self.liststore[path][0]
    
    def clear_treeview (self, widget):
        """ Function doc """
        #print('clear_treeview')
        self.liststore.clear()
    
    
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
            #print(folder_paths)
            
            liststore = self.builder.get_object('liststore1')
            
            for data in folder_paths:
                self.liststore.append([True, data])

        dialog.destroy()
            
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
        #self._starting_coordinates_model_update()
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
    def on_combobox_PMF_type(self, widget):
        _type     = self.methods_combo.get_active()
        
        if _type ==1:
            self.builder.get_object('entry_bins_RC2').show()
            self.builder.get_object('label_bins_RC2').show()
            
            #parameters['bins'] = [int(bins), int(bins)] 

        else:
            self.builder.get_object('entry_bins_RC2').hide()
            self.builder.get_object('label_bins_RC2').hide()
        
        #self.builder.get_object('label_bins_RC2').hide()
        
    
    def on_combobox_systemsbox_changed(self, widget):
        """ Function doc """
        system_id = self.combobox_systems.get_system_id()
        
        working_folder = self.p_session.psystem[system_id].e_working_folder
        self.folder_chooser_button.set_folder(working_folder)    
        #self.update_window ( selections = False, restraints = True)
