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


class ExportScriptDialog:
    """ 
    A GTK-based dialog for exporting simulation scripts. 
    Allows the user to choose a filename, folder, and automatically 
    generates a Python script for running simulations with pDynamo. 
    """
    
    def __init__(self, main=None, parameters=None):
        """ 
        Initialize the dialog. 
        - main: reference to the main application object. 
        - parameters: dictionary with simulation configuration parameters. 
        """
        self.main       = main
        self.home       = main.home                  # Home directory of the main application
        self.parameters = parameters
        
        # Load the GUI layout from a Glade file
        self.builder = Gtk.Builder()
        self.builder.add_from_file(
            os.path.join(self.home, 'src/gui/windows/setup/export_script_window2.glade')
        )
        
        # Connect signals for window and button actions
        self.builder.get_object('dialog').connect('destroy', self.close_window)
        self.builder.get_object('button_cancel').connect('clicked', self.close_window)
        self.builder.get_object('button_export').connect('clicked', self.on_button_export_clicked)
        
        # Set up the folder chooser with the default path from parameters
        self.folder_chooser = self.builder.get_object('btn_folder_chooser')
        self.folder_chooser.set_current_folder(self.parameters['folder'])
         
        # Define dialog size and run it (modal behavior)
        self.builder.get_object('dialog').set_default_size(400, 100)
        self.builder.get_object('dialog').run()
    

    def on_button_export_clicked(self, widget):
        """ 
        Event handler for the "Export" button. 
        Collects input values, exports the system, and generates a Python script 
        that reproduces the simulation setup. 
        """
        
        # Path to internal pDynamo methods
        pmethod_path = os.path.join(self.main.home, 'src/pdynamo')
        
        # Get filename and folder from user input
        filename     = self.builder.get_object('entry_filename').get_text()
        outputfolder = self.builder.get_object('btn_folder_chooser').get_current_folder()
        
        # Export parameters for the system
        exp_parameters = { 
            'vobject_id' : 0,
            'folder'     : outputfolder,
            'filename'   : filename,
            'format'     : 0,
            'last'       : -1,
            'system_id'  : self.main.p_session.active_id
        }
        # Export the system as a .pkl file
        self.main.p_session.export_system(exp_parameters, False)
        fullfilename = os.path.join(outputfolder, filename + '.pkl')

        # Special case: Nudged Elastic Band simulation requires coordinate export
        if self.parameters['simulation_type'] == 'Nudged_Elastic_Band':
            Pickle(os.path.join(outputfolder, filename + 'reac_coordinates.pkl'),
                   self.parameters['reac_coordinates'])
            Pickle(os.path.join(outputfolder, filename + 'prod_coordinates.pkl'),
                   self.parameters['prod_coordinates'])
            self.parameters['reac_coordinates'] = os.path.join(outputfolder, filename + 'reac_coordinates.pkl')
            self.parameters['prod_coordinates'] = os.path.join(outputfolder, filename + 'prod_coordinates.pkl')
            
        # Build the Python script header with imports and system loading
        header = ''' 
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
import sys

# This section should be edited
#.Enter the path to the folder ../EasyHybid/pDynamo here
sys.path.insert(1, '{}')
import p_methods  as pMethods

# Importing system
system = ImportSystem ('{}')
# system.qcState.DeterminePaths('/home/usr/scratch') # Uncomment for ORCA / DFTB+ / xTB

# Summary of the imported system
system.Summary()
'''.format(pmethod_path, fullfilename)
        
        # Append parameters dictionary in pretty format
        header += '\nparameters = '
        pp = pprint.PrettyPrinter(indent=4)
        formatted_data = pp.pformat(self.parameters)
        header += formatted_data
        header += "\nparameters['system'] = system"
        
        # Add simulation object based on the type
        if self.parameters['simulation_type'] == 'Umbrella_Sampling':
            header += '\n\nsimObj =  pMethods.UmbrellaSampling()'
        
        elif self.parameters['simulation_type'] == 'Relaxed_Surface_Scan':
            header += '\n\nsimObj =  pMethods.RelaxedSurfaceScan()'
        
        elif self.parameters['simulation_type'] == 'Nudged_Elastic_Band':
            header += '\n\nsimObj =  pMethods.ChainOfStatesOptimizePath()'
        
        # Add execution line
        header += '\nsimObj.run(parameters)' 
        
        # Save the generated script as a .py file
        fullfilename = os.path.join(outputfolder, filename + '.py')
        with open(fullfilename, 'w') as file:
            file.write(header)
        
        # Close the dialog
        self.builder.get_object('dialog').destroy()
        self.Visible = False
        

    def close_window(self, button, data=None):
        """ 
        Close the dialog window and mark it as not visible. 
        """
        self.builder.get_object('dialog').destroy()
        self.Visible = False


class ExportDataWindow:
    """ 
    A GTK-based window for exporting molecular system data.  
    The user can choose the system, file format, trajectory frames, and output folder.  
    Supports single-file and multi-file export, including 2D trajectories. 
    """

    def __init__(self, main=None):
        """ 
        Initialize the ExportDataWindow. 
        
        Parameters:
        - main: reference to the main application object. 
        """
        self.main    = main
        self.home    = main.home
        self.Visible = False     # Tracks whether the window is already open
        self.starting_coords_liststore = Gtk.ListStore(str, int)  # Stores available coordinates
        self.is_single_frame  = True   # Default: exporting a single frame
        self.is_2D_trajectory = False  # Default: not a 2D trajectory


    def close_window(self, button, data=None):
        """ 
        Close the export window and update visibility state. 
        """
        self.window.destroy()
        self.Visible = False


    def open_window(self, sys_selected=None):
        """ 
        Open the Export Data window.  
        Creates and initializes all GUI components if the window is not already open. 
        
        Parameters:
        - sys_selected: optional system ID to preselect in the system combobox. 
        """
        if not self.Visible:
            # Load GUI from Glade file
            self.builder = Gtk.Builder()
            self.builder.add_from_file(
                os.path.join(self.home, 'src/gui/windows/setup/export_system_window.glade')
            )
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('export_data_window')
            self.window.set_title('Export Data Window')
            
            # ------------------ File formats combobox ------------------ #
            self.format_store = Gtk.ListStore(str)
            self.formats = {
                0: 'pkl - pDynamo system',
                1: 'yaml - pDynamo system',
                2: 'pkl - pDynamo coordinates',
                3: 'pdb - Protein Data Bank',
                4: 'xyz - cartesian coordinates',
                5: 'mol',
                6: 'mol2',
                7: 'crd'
            }
            for _, _format in self.formats.items():
                self.format_store.append([_format])
            
            self.combobox_fileformat = self.builder.get_object('combobox_fileformat')
            self.combobox_fileformat.set_model(self.format_store)
            renderer_text = Gtk.CellRendererText()
            self.combobox_fileformat.pack_start(renderer_text, True)
            self.combobox_fileformat.add_attribute(renderer_text, "text", 0)
            # ----------------------------------------------------------- #
            
            
            # ------------------ System combobox ------------------ #
            self.box = self.builder.get_object('box_system')
            self.psystem_combo = SystemComboBox(self.main)
            self.box.pack_start(self.psystem_combo, False, False, 0)
            self.psystem_combo.connect("changed", self.combobox_pdynamo_system)
            if sys_selected:
                self.psystem_combo.set_active_system(sys_selected)
            else:
                self.psystem_combo.set_active(0)
            # ----------------------------------------------------- #
            

            sys_selected = self.psystem_combo.get_system_id()
            
            # ------------------ Folder chooser ------------------ #
            self.folder_chooser_button = FolderChooserButton(main=self.main)
            self.builder.get_object('folder_chooser_box').pack_start(
                self.folder_chooser_button.btn, True, True, 0
            )
            # ----------------------------------------------------- #
            

            # ------------------ Starting coordinates combobox ------------------ #
            self.box_coordinates = self.builder.get_object('box_coordinates')
            try:
                # Populate coordinates from the selected system
                self.combobox_starting_coordinates = CoordinatesComboBox(
                    coordinates_liststore=self.main.vobject_liststore_dict[sys_selected]
                )
            except:
                self.combobox_starting_coordinates = CoordinatesComboBox(coordinates_liststore=-1)
            
            self.box_coordinates.pack_start(self.combobox_starting_coordinates, False, False, 0)
            size = len(self.main.vobject_liststore_dict[sys_selected])
            
            self.combobox_starting_coordinates.connect("changed", self.on_vobject_combo_changed)
            self.combobox_starting_coordinates.set_active(size-1)
            # ------------------------------------------------------------------- #
            

            self.checkbox_keep_window = self.builder.get_object('checkbox_keep_window_open')
            
            # Set default folder to system’s working folder if available
            folder = self.main.p_session.psystem[sys_selected].e_working_folder
            if folder:
                self.folder_chooser_button.set_folder(folder=folder)

            # Show window
            self.window.show_all()
            self.is_single_file = True

            # Disable multiple-frame entries by default
            self.builder.get_object('entry_first').set_sensitive(False)
            self.builder.get_object('label_first').set_sensitive(False)
            self.builder.get_object('entry_stride').set_sensitive(False)
            self.builder.get_object('label_stride').set_sensitive(False)
            self.builder.get_object('entry_first2').set_sensitive(False)
            self.builder.get_object('entry_stride2').set_sensitive(False)
            
            self.combobox_fileformat.set_active(0)
            self.Visible = True
        else:
            # If window already open, bring it to front
            self.window.present()


    def combobox_pdynamo_system(self, widget):
        """ 
        Event handler when system combobox changes.  
        Updates coordinates combobox and folder chooser accordingly. 
        """
        if self.Visible:
            sys_id = widget.get_system_id()
            self.combobox_starting_coordinates.set_model(
                self.main.vobject_liststore_dict[sys_id]
            )
            size = len(self.main.vobject_liststore_dict[sys_id])
            self.combobox_starting_coordinates.set_active(size-1)

            folder = self.main.p_session.psystem[sys_id].e_working_folder
            if folder:
                self.folder_chooser_button.set_folder(folder=folder)


    def on_vobject_combo_changed(self, widget):
        """ 
        Event handler for changing starting coordinates.  
        Detects if the selected trajectory is 2D and updates UI elements accordingly. 
        """
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            model = self.combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]
            vismol_object = self.main.vm_session.vm_objects_dic[vobject_id]
        
            if getattr(vismol_object, 'idx_2D_xy', False):
                # Show UI for 2D trajectory export
                self.builder.get_object('label_y_rc').show()
                self.builder.get_object('entry_first2').show()
                self.builder.get_object('entry_last2').show()
                self.builder.get_object('entry_stride2').show()
                self.is_2D_trajectory = True
            else:
                # Hide 2D-related inputs
                self.builder.get_object('label_y_rc').hide()
                self.builder.get_object('entry_first2').hide()
                self.builder.get_object('entry_last2').hide()
                self.builder.get_object('entry_stride2').hide()
                self.is_2D_trajectory = False
        self._fileformat_update()


    def on_combobox_fileformat(self, widget):
        """ 
        Event handler for changing file format combobox.  
        Updates UI based on file format selection. 
        """
        self._fileformat_update()


    def _fileformat_update(self):
        """ 
        Internal function to update UI depending on file format  
        and number of frames in the selected object. 
        """
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            model = self.combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]
        
        if len(self.main.vm_session.vm_objects_dic[vobject_id].frames) > 1:
            self.is_single_frame = True
            if self.combobox_fileformat.get_active() in [0, 1]:
                self.builder.get_object('radiobutton_singlefile').set_active(True)
                self.builder.get_object('radiobutton_multiplefile').set_sensitive(False)
            else:
                self.builder.get_object('radiobutton_multiplefile').set_sensitive(True)
        else:
            # Single-frame object → force single-file export
            self.builder.get_object('radiobutton_singlefile').set_active(True)
            self.builder.get_object('radiobutton_multiplefile').set_sensitive(False)
            self.is_single_frame = True


    def on_name_combo_changed(self, widget):
        """ 
        Event handler for name-type combobox.  
        Changes folder chooser type between file and folder selection. 
        """
        if widget.get_active() == 0:
            self.folder_chooser_button.sel_type = 'folder'
        else:
            self.folder_chooser_button.sel_type = 'file'


    def on_radio_button_change(self, widget):
        """ 
        Event handler for switching between single-file and multiple-file export.  
        Enables/disables UI fields accordingly. 
        """
        if self.builder.get_object('radiobutton_singlefile').get_active():
            self.is_single_file = True
            # Disable multiple-frame options
            self.builder.get_object('entry_first').set_sensitive(False)
            self.builder.get_object('label_first').set_sensitive(False)
            self.builder.get_object('entry_stride').set_sensitive(False)
            self.builder.get_object('label_stride').set_sensitive(False)
            self.builder.get_object('entry_first2').set_sensitive(False)
            self.builder.get_object('entry_stride2').set_sensitive(False)
        else:
            self.is_single_file = False
            # Enable multiple-frame options
            self.builder.get_object('entry_first').set_sensitive(True)
            self.builder.get_object('label_first').set_sensitive(True)
            self.builder.get_object('entry_stride').set_sensitive(True)
            self.builder.get_object('label_stride').set_sensitive(True)
            self.builder.get_object('entry_first2').set_sensitive(True)
            self.builder.get_object('entry_stride2').set_sensitive(True)


    def export_data(self, button):
        """ 
        Collects all user-selected parameters and exports system data.  
        Supports different formats, single/multiple-file export, and 2D trajectories. 
        """
        # Default parameters
        parameters = {
            'system_id': None,
            'vobject_id': None, 
            'format': None,
            'is_single_file': True,
            'is_2D_trajectory': False,
            'first': 0,
            'last': -1,
            'stride': 1,
            'first2': 0,
            'last2': -1,
            'stride2': 1,
            'export_QC_atoms_only': False,
            'filename': 'exported_system',
            'folder': None,
        }
        
        # ----------------- Get selected starting coordinates ----------------- #
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            model = self.combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]
        parameters['vobject_id'] = vobject_id
        # --------------------------------------------------------------------- #
        
        # ----------------- File format ----------------- #
        _format = self.combobox_fileformat.get_active()
        parameters['format'] = _format
        # ------------------------------------------------ #
        
        # ----------------- Output folder and filename ----------------- #
        folder = self.folder_chooser_button.get_folder()
        filename = self.builder.get_object('entry_filename').get_text()
        parameters['folder'] = folder
        parameters['filename'] = filename
        # --------------------------------------------------------------- #
        
        # ----------------- Selected system ----------------- #
        tree_iter = self.psystem_combo.get_active_iter()
        if tree_iter is not None:
            model = self.psystem_combo.get_model()
            name, system_id = model[tree_iter][:2]
        parameters['system_id'] = system_id
        # --------------------------------------------------- #
        
        # ----------------- Single or multiple file ----------------- #
        if self.builder.get_object('radiobutton_singlefile').get_active():
            parameters['is_single_file'] = True
        else:
            parameters['is_single_file'] = False
        # ----------------------------------------------------------- #
        
        parameters['is_2D_trajectory'] = self.is_2D_trajectory
        
        # ----------------- Frame ranges ----------------- #
        parameters['first']  = int(self.builder.get_object('entry_first').get_text())
        parameters['last']   = int(self.builder.get_object('entry_last').get_text())
        parameters['stride'] = int(self.builder.get_object('entry_stride').get_text())
        parameters['first2'] = int(self.builder.get_object('entry_first2').get_text())
        parameters['last2']  = int(self.builder.get_object('entry_last2').get_text())
        parameters['stride2'] = int(self.builder.get_object('entry_stride2').get_text())
        # ------------------------------------------------ #
        
        parameters['export_QC_atoms_only'] = self.builder.get_object('checkbox_export_QC_atoms_only').get_active()
        
        # Attach the actual system object
        parameters['system'] = self.main.p_session.psystem[parameters['system_id']]
        
        # ----------------- File extension handling ----------------- #
        format_dict = {0:'pkl', 1:'yaml', 2:'pkl', 3:'pdb', 4:'xyz', 5:'mol', 6:'mol2', 7:'crd'}
        if parameters['is_single_file']:
            _format = '.' + format_dict[parameters['format']]
        else:
            _format = '.ptGeo'
        # ------------------------------------------------------------ #
        
        #print(parameters)  # Debug
        # Call export function from session
        try:
            self.main.p_session.export_system(parameters)
            
            system = self.main.p_session.psystem[parameters['system_id']]
            path = os.path.join(parameters['folder'],parameters['filename']+'.'+format_dict[parameters['format']])
            self.main.bottom_notebook.status_teeview_add_new_item(message = 'Data exported successfully: {}'.format(path), system = system)
        
        except Exception as e:
            error_str = str(e)  # converte a mensagem de erro para string
            print("Error:", error_str)
            self.main.bottom_notebook.status_teeview_add_new_item(message = 'Error: Could not export data.', system = None)
            simpledialog = SimpleDialog(self.main)
            simpledialog.error("Error: Could not export data.")
            
        #'''
        
        if self.checkbox_keep_window.get_active():
            pass
        else:
            self.window.destroy()
            self.Visible = False
        
        
        
        
        
        
        
        # Optional: status update (commented out in your code)
        '''
        try:
            self.main.p_session.export_system(parameters)
            self.main.bottom_notebook.status_teeview_add_new_item(
                message=':  {}  saved'.format(os.path.join(
                    parameters['folder'], 
                    parameters['filename']+_format)
                ), 
                system=parameters['system']
            )
        except:
            print('Failed when trying to export system data: ', parameters['system'].label)
        
        if self.checkbox_keep_window.get_active():
            pass
        else:
            self.window.destroy()
            self.Visible = False
        '''
        

    def update(self):
        """ 
        Placeholder function for future updates to the window. 
        """
        pass
