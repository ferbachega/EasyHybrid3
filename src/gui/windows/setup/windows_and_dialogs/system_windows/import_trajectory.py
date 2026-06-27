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

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Pango
from gi.repository import Gdk
from gi.repository import GdkPixbuf
import cairo


#---------------------------------------
from pBabel                    import*                                     
from pCore                     import*  
#---------------------------------------
from pMolecule                 import*                              
from pMolecule.MMModel         import*
from pMolecule.NBModel         import*                                     
from pMolecule.QCModel         import*
#---------------------------------------
from pScientific               import*                                     
from pScientific.Arrays        import*                                     
from pScientific.Geometry3     import*                                     
from pScientific.RandomNumbers import*                                     
from pScientific.Statistics    import*
from pScientific.Symmetry      import*
#---------------------------------------                              
from pSimulation               import*
#---------------------------------------

#import Pickle
from gui.widgets.custom_widgets import FolderChooserButton
from gui.widgets.custom_widgets import VismolTrajectoryFrame
from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.widgets.custom_widgets import get_colorful_square_pixel_buffer
from gui.widgets.custom_widgets import ReactionCoordinateBox

#from gui.widgets.custom_widgets import get_distance
from util.geometric_analysis            import get_distance 
from util.geometric_analysis            import get_dihedral 
from util.geometric_analysis            import get_angle 

from pdynamo.p_methods import LogFile


import util.orca_qc_keywords as orca_keys

from util.file_parser import get_file_type  
from util.file_parser import read_MOL2  
import pprint
import numpy as np
import gc
import os

import traceback

# --- cross-module imports added during refactor ---
from gui.windows.setup.windows_and_dialogs.dialogs.simple_dialog import SimpleDialog


VISMOL_HOME      = os.environ.get('VISMOL_HOME')
HOME             = os.environ.get('HOME')
PDYNAMO3_SCRATCH = os.environ.get('PDYNAMO3_SCRATCH')


class ImportTrajectoryWindow:
    """ Class doc """
    def __init__(self, main = None):
        """ Class initialiser """
        self.main                = main
        self.p_session           = main.p_session
        self.home                = main.home
        
        self.Visible             =  False        
        self.starting_coords_liststore = Gtk.ListStore(str, int)
        
        self.data_type_dict = {
                        0:'pklfile'    , #  single file
                        1:'pklfolder'  , #  trajectory
                        2:'pklfolder2D', #  2d trajectory  
                        3:'pdbfile'    ,
                        4:'pdbfolder'  ,
                        5:'dcd',
                        6:'crd',
                        7:'xyz',
                        8:'mol2',
                        9:'netcdf',
                       10:'log_file',
                       11:'charges'  }
        
        self.folder_type_list = ['pklfolder', 'pklfolder2D', 'pdbfolder']

    def open_window (self, sys_selected = None):
        """ Function doc """
        if self.Visible  ==  False:
            '''--------------------------------------------------------------------------------------------'''
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'src/gui/windows/setup/import_trajectory_window.glade'))
            self.builder.connect_signals(self)
            self.window = self.builder.get_object('import_trajectory_window')
            self.window.set_title('Import Data Window')
            self.window.set_keep_above(True)
            '''--------------------------------------------------------------------------------------------'''
            self.box = self.builder.get_object('box_system')
            #------------------------------------------------------------------------------------
            
            self.combobox_pdynamo_system = SystemComboBox (self.main)
            self.box.pack_start(self.combobox_pdynamo_system, False, False, 0)
            if sys_selected:
                self.combobox_pdynamo_system.set_active_system(sys_selected)
            else:
                self.combobox_pdynamo_system.set_active(0)
            #self.combobox_pdynamo_system = self.builder.get_object('combobox_pdynamo_system')
            #renderer_text = Gtk.CellRendererText()
            #self.combobox_pdynamo_system.add_attribute(renderer_text, "text", 0)
            #
            #
            #
            ##------------------------------------------------------------------------------------
            #"""This block of instructions is used to organize the system combobox. 
            #The menu function provides the system ID, but the position in the combobox 
            #list must be calculated dynamically, as the index of items will not match 
            #the system ID when a system has been deleted."""
            ##self._define_system_liststore()
            ##self.system_liststore = Gtk.ListStore(str, int)
            ##names      = [ ]
            #
            #counter    = 0
            #set_active = 0
            #for key , system in self.p_session.psystem.items():
            #    if system:
            #        #name = system.label
            #        #self.system_liststore.append([name, int(key)])
            #        if sys_selected == int(key):
            #            set_active = counter
            #            counter += 1
            #        else:
            #            counter += 1
            #    else:
            #        pass
            ##------------------------------------------------------------------------------------
            #
            #
            ##self.combobox_pdynamo_system.set_model(self.system_liststore)
            #self.combobox_pdynamo_system.set_model(self.main.system_liststore)
            #
            #if sys_selected:
            #    self.combobox_pdynamo_system.set_active(set_active)
            #else:
            #    self.combobox_pdynamo_system.set_active(set_active)
            ##------------------------------------------------------------------------------------
            
            
            
            
            
            #'''--------------------------------------------------------------------------------------------------
            self.folder_chooser_button = FolderChooserButton(main =  self.main)
            self.builder.get_object('folder_chooser_box').pack_start(self.folder_chooser_button.btn, True, True, 0)
            #'''--------------------------------------------------------------------------------------------------
            
            
            
            #------------------------------------------------------------------------------------
            self.combobox_starting_coordinates = self.builder.get_object('vobjects_combobox')
            #self.combobox_starting_coordinates.set_model(self.starting_coords_liststore)
            
            renderer_text = Gtk.CellRendererText()
            self.combobox_starting_coordinates.pack_start(renderer_text, True)
            self.combobox_starting_coordinates.add_attribute(renderer_text, "text", 0)
            
            size = len(self.starting_coords_liststore)
            self.combobox_starting_coordinates.set_active(size-1)
            #------------------------------------------------------------------------------------
            
            
            #------------------------------------------------------------------------------------
            self.builder.get_object('vobjects_combobox').set_sensitive(False)
            #------------------------------------------------------------------------------------

            #'''--------------------------------------------------------------------------------------------'''
            self.combox = self.builder.get_object('combobox_coordinate_type')
            self.combox.connect("changed", self.on_combobox_coordinate_type)
    
            self.folder_chooser_button.btn.connect('clicked', self.update_logfile_chooser_btn)
            self.on_combobox_pdynamo_system(None)
            self.combox.set_active(0)
            self.window.show_all()
            self.Visible  = True
        
        else:
            self.window.present()
    
    def close_window (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False
    
    def update_logfile_chooser_btn (self, widget):
        """ Function doc """
        #print(self.folder_chooser_button.folder)
        trajfolder = self.folder_chooser_button.get_folder()
        #print(trajfolder)
        
        
        isdir = os.path.isdir(trajfolder)
        if isdir:
            files = os.listdir(trajfolder)
        
        
            logfiles = []
            for file_name in files:
                if file_name.endswith("log"):
                    logfiles.append(file_name)


            #basename  = os.path.basename(trajfolder)
            #logfile   = basename[:-5]+'log'

            #logfile   = os.path.join(trajfolder, logfile)
            #if exists(logfiles[0]):
            if len(logfiles):
                fullpath = os.path.join(trajfolder, logfiles[0])
                self.builder.get_object('file_chooser_btn_logfile').set_filename(fullpath)
        
    def on_combobox_pdynamo_system (self, widget):
        """ Function doc """
        tree_iter = self.combobox_pdynamo_system.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.combobox_pdynamo_system.get_model()
            name, sys_id = model[tree_iter][:2]
            #print ('name/ system_id: ', name, sys_id)
        else:
            return False

        self.combobox_starting_coordinates = self.builder.get_object('vobjects_combobox')
        self.combobox_starting_coordinates.set_model(self.main.vobject_liststore_dict[sys_id])
        self.folder_chooser_button.folder = self.main.p_session.psystem[sys_id].e_working_folder

    def on_combobox_coordinate_type (self, widget):
        """ Function doc """
        data_type = self.builder.get_object('combobox_coordinate_type').get_active() 

        
        
        if data_type == 10:
            self.builder.get_object('frame_stride_box').set_sensitive(False) 
            self.builder.get_object('folder_chooser_box').set_sensitive(False) 
            self.builder.get_object('entry_create_a_new_vobj').set_sensitive(False)
            self.builder.get_object('vobjects_combobox').set_sensitive(True)

        
        elif data_type == 11:
            self.builder.get_object('frame_stride_box').set_sensitive(False) 
            self.builder.get_object('folder_chooser_box').set_sensitive(True) 
            self.builder.get_object('entry_create_a_new_vobj').set_sensitive(False)
            self.builder.get_object('vobjects_combobox').set_sensitive(False)
        
        elif data_type == 0:
            self.builder.get_object('frame_stride_box').set_sensitive(False) 
            self.builder.get_object('folder_chooser_box').set_sensitive(True) 
            self.builder.get_object('entry_create_a_new_vobj').set_sensitive(True)
            self.builder.get_object('vobjects_combobox').set_sensitive(True)

        elif data_type in [1,3,5,6,7,8,9]:
            self.builder.get_object('frame_stride_box').set_sensitive(True) 
            self.builder.get_object('folder_chooser_box').set_sensitive(True) 
            self.builder.get_object('entry_create_a_new_vobj').set_sensitive(True)
            self.builder.get_object('vobjects_combobox').set_sensitive(True)

        elif data_type in [2,4]:
            self.builder.get_object('frame_stride_box').set_sensitive(False) 
            self.builder.get_object('folder_chooser_box').set_sensitive(True) 
            self.builder.get_object('entry_create_a_new_vobj').set_sensitive(True)
            self.builder.get_object('vobjects_combobox').set_sensitive(False)
        else:
            self.builder.get_object('frame_stride_box').set_sensitive(True) 
            self.builder.get_object('folder_chooser_box').set_sensitive(True) 
            self.builder.get_object('entry_create_a_new_vobj').set_sensitive(True)
            self.builder.get_object('vobjects_combobox').set_sensitive(True)        
        
        #print ('\ndata_type: ', data_type, '\ndata_type_dict: ',self.data_type_dict[data_type])
        data_type = self.data_type_dict[data_type]
        
        if  data_type in self.folder_type_list:
            self.folder_chooser_button.sel_type = 'folder'
        else:
            self.folder_chooser_button.sel_type = 'file'        

    def on_vobject_combo_changed (self, widget):
        '''this combobox has the reference to the starting coordinates of a simulation'''
        #combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]

    def on_name_combo_changed (self, widget):
        """ Function doc """
        traj_type = self.builder.get_object('combobox_coordinate_type').get_active() 
        #print (traj_type, self.data_type_dict[traj_type])
        traj_type = self.data_type_dict[traj_type]
        
        if  traj_type in self.folder_type_list:
            self.folder_chooser_button.sel_type = 'folder'
        else:
            self.folder_chooser_button.sel_type = 'file'

    def on_radio_button_change (self, widget):
        """ Function doc """
        if self.builder.get_object('radiobutton_import_as_new_object').get_active():
            self.builder.get_object('vobjects_combobox').set_sensitive(False)
            self.builder.get_object('entry_create_a_new_vobj').set_sensitive(True)

        if self.builder.get_object('radiobutton_append_to_a_vobject').get_active():
            self.builder.get_object('entry_create_a_new_vobj').set_sensitive(False)
            self.builder.get_object('vobjects_combobox').set_sensitive(True)

    def import_data (self, button):
        """ Function doc """
        
        parameters = {
                      'system_id'    : None,
                                     
                      'data_path'    : None,
                      'data_type'    : 0   ,
                      
                      'new_vobj_name': None, 
                      'vobject_id'   : None,
                      'vobject'      : None,
                      
                      'isAppend'     : False,
                      
                      'logfile'      : None,
                      
                      'first'        : None,
                      'last'         : None,
                      'stride'       : None,
                     }
        
        #entry_name     = None
        # - - - system - - - 
        tree_iter = self.combobox_pdynamo_system.get_active()
        model = self.combobox_pdynamo_system.get_model()
        name, sys_id = model[tree_iter][:2]
        parameters['system_id'] = sys_id
        
        # data, a file or a folder - many files
        parameters['data_path'] = self.folder_chooser_button.get_folder()
        
        # A log file - optional
        parameters['logfile']   = self.builder.get_object('file_chooser_btn_logfile').get_filename()
        
        
        data_type = self.builder.get_object('combobox_coordinate_type').get_active()        
        parameters['data_type'] = self.data_type_dict[data_type]
        
        try:
            parameters['first']  = int(self.builder.get_object('entry_first').get_text()  )
            parameters['last']   = int(self.builder.get_object('entry_last').get_text()   )
            parameters['stride'] = int(self.builder.get_object('entry_stride').get_text() )
        except:
            pass
        
        #----------------------------------------------------------------------------------------------
        '''Here, we determine whether it is necessary to create a new object 
        or append the coordinates to an existing one.'''
        if self.builder.get_object('radiobutton_import_as_new_object').get_active():
            parameters['new_vobj_name'] = self.builder.get_object('entry_create_a_new_vobj').get_text()
        else:
            #-----------------------------------------------------------------------------
            tree_iter = self.combobox_starting_coordinates.get_active_iter()
            if tree_iter is not None:
                
                '''selecting the vismol object from the content that is in the combobox '''
                model = self.combobox_starting_coordinates.get_model()
                name, vobject_id = model[tree_iter][:2]
                vobject = self.main.vm_session.vm_objects_dic[vobject_id]
                
                parameters['vobject_id'] = vobject_id
                parameters['vobject']    = vobject
                parameters['isAppend']   = True
                
            #-------------------------------------------------------------------------------------------

        #print('\n parameters: ', parameters)
        try:
            self.main.p_session.import_data ( parameters ) 
        except Exception as e:
            tb = traceback.format_exc()
            self.main.bottom_notebook.status_teeview_add_new_item(
                message=f"Error loading data: {e}", system=None)
                #message=f"Error loading data: {e}\n{tb}", system=None)
            
            simpledialog = SimpleDialog(self.main )
            simpledialog.error(f"Error loading data: {e}")
            
            print(f"Error loading data: {e}\n{tb}")
           
            #system = None
        
    def update (self):
        """ Function doc """
        #print('VismolGoToAtomWindow2 update')
        pass
