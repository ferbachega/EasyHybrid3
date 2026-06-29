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


class ImportANewSystemWindow(Gtk.Window):
    """
    Window for importing a new molecular system into EasyHybrid.
    
    This class provides a GTK-based interface that allows the user to 
    load molecular systems prepared in different formats (AMBER, CHARMM, 
    OPLS, DYFF, pDynamo, etc.) and configure their working folder.
    """
    def __init__(self, main = None, home = None):
        """
        Initialize the ImportANewSystemWindow class.

        Parameters
        ----------
        main : object
            Reference to the main EasyHybrid instance.
        home : str
            Path to the EasyHybrid home directory.
        """
        
        self.easyhybrid_main = main
        self.home = main.home
        self.Visible = False
        self.vm_session = self.easyhybrid_main.vm_session

        # Store for loaded residue data (filename, type, number of atoms)
        self.residue_liststore = Gtk.ListStore(str, str, str)

        # Default working folder
        #self.folder = self.vm_session.vm_config.gl_parameters["workspace_path"]
        
        # Descriptions for different system types
        self.charmm_txt = '''When using the traditional CHARMM/PSF format for system preparation, it is necessary to have three types of files: parameter files (in formats such as prm or par), topology files in the form of psf, and coordinate files (which can be in various formats including chm, crd, pdb, xyz, etc. It's worth noting that multiple parameter files may be required depending on the system being simulated.'''
        self.amber_txt = '''In order to properly run simulations utilizing the AMBER force field, it is necessary to have two types of files: topologies (in the form of either top or prmtop files) and coordinates (which can be in the form of crd, pdb, xyz, or other similar file types).'''
        self.OPLS_txt = '''When preparing systems natively in pDynamo utilizing the OPLS force field, it is necessary to have two components: a folder containing the OPLS parameters, and a topology and coordinate file (in formats such as pdb, mol2, or mol).'''
        self.DYFF_txt = '''The DYFF force field is a generic force field specifically designed for use within the pDynamo program. When utilizing DYFF to prepare systems natively in pDynamo, it is essential to have two components: a folder containing the necessary force field parameters, and a coordinate file in formats such as pdb, mol2, or mol.'''
        self.gmx_txt = '''Systems prepared natively in GROMACS using the CHARMM(or AMBER) force field require: A parameter/topology (top) file and coordinate file (pdb, mol2, mol, ...) .  '''
        self.xyz_pdb_text = '''If you load a coordinate file in this manner, you won't be able to perform any molecular mechanics simulations.'''
        self.pkl_text = '''If you load a coordinate file in this manner, you won't be able to perform any molecular mechanics simulations.'''
    
    
        self.color_pallet = {
            0: (124 / 255, 252 / 255, 0 / 255, 255 / 255),   # Lawn green
            1: (238 / 255, 130 / 255, 238 / 255, 255 / 255), # Violet
            2: (255 / 255, 255 / 255, 0 / 255, 255 / 255),   # Yellow
            3: (135 / 255, 206 / 255, 250 / 255, 255 / 255), # Light sky blue
            4: (255 / 255, 99 / 255, 71 / 255, 255 / 255),   # Tomato
            5: (0 / 255, 250 / 255, 154 / 255, 255 / 255),   # Medium spring green
            6: (123 / 255, 104 / 255, 238 / 255, 255 / 255), # Medium slate blue
            7: (238 / 255, 232 / 255, 170 / 255, 255 / 255), # Pale goldenrod
            8: (224 / 255, 255 / 255, 255 / 255, 255 / 255), # Light cyan
            9: (219 / 255, 112 / 255, 147 / 255, 255 / 255), # Pale violet red
            10: (255 / 255, 215 / 255, 0 / 255, 255 / 255),  # Gold
        }
        
        # Counter used to assign colors cyclically
        self.color_counter = 0
    
    def open_window(self):
        """Open and configure the 'Import New System' window."""

        if not self.Visible:
            # Load the Glade UI definition file
            self.builder = Gtk.Builder()
            self.builder.add_from_file(
                os.path.join(self.home, 'src/gui/windows/setup/import_system_window_new.glade')
            )
            self.builder.connect_signals(self)

            # Retrieve and configure the main window
            self.window = self.builder.get_object('ImportNewSystemWindow')
            self.window.set_border_width(10)
            self.window.set_default_size(500, 370)

            # --------------------------------------------------------------------------------------------
            # Create a ListStore to hold system types
            self.system_type_store = Gtk.ListStore(str, int)
            system_types = [
                ["AMBER", 0],
                ["CHARMM", 1],
                ["OPLS", 2],
                ['DYFF', 5],
                ["pdynamo files (*.pkl, *.yaml)", 3],
                ["other (*.pdb, *.xyz, *.mol2)", 4],
            ]
            for system_type in system_types:
                self.system_type_store.append(system_type)

            # Create a combo box to select the system type
            self.system_types_combo = Gtk.ComboBox.new_with_model(self.system_type_store)
            self.system_types_combo.set_tooltip_text(self.amber_txt)
            self.box_combo = self.builder.get_object('box')
            self.box_combo.pack_start(self.system_types_combo, True, True, 0)

            # Hide unused label (legacy support)
            self.builder.get_object('gtk_label_fftype').hide()

            # Connect callback for when the system type changes
            self.system_types_combo.connect("changed", self.on_name_combo_changed)
            self.system_types_combo.set_model(self.system_type_store)

            # Renderer to display system type names in the combo box
            renderer_text = Gtk.CellRendererText()
            self.system_types_combo.pack_start(renderer_text, True)
            self.system_types_combo.add_attribute(renderer_text, "text", 0)
            # --------------------------------------------------------------------------------------------

            # Configure the treeview for listing loaded files
            self.treeview = self.builder.get_object('gtktreeview_import_system')
            for i, column_title in enumerate(['file', "type", "number of atoms"]):
                renderer = Gtk.CellRendererText()
                column = Gtk.TreeViewColumn(column_title, renderer, text=i)
                self.treeview.append_column(column)

            # --------------------------------------------------------------------------------------------
            # System name entry field
            self.entry_system_name = self.builder.get_object('entry_system_name')
            self.entry_system_name.connect('changed', self.on_entry_system_name_change)

            # Show the window and connect signals
            self.window.show_all()
            self.builder.connect_signals(self)

            # Hide OPLS folder chooser by default
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()

            # Mark window as visible
            self.Visible = True

            # Initialize dictionary to store input files
            self.files = {
                'amber_prmtop': None,
                'charmm_par': [],
                'charmm_psf': None,
                'charmm_extra': None,
                'prm_folder': [],
                'coordinates': None,
            }
            # Set default system type to AMBER
            self.system_types_combo.set_active(0)

            # Connect window and button signals
            self.window.connect('destroy', self.close_window)
            self.builder.get_object('button_load_files').connect('clicked', self.on_button_load_files_clicked)
            self.builder.get_object('button_remove_files').connect('clicked', self.on_button_delete_files_clicked)
            self.builder.get_object('button_cancel').connect('clicked', self.close_window)
            self.builder.get_object('import_import_system').connect('clicked', self.on_button_import_system_clicked)

            # Configure the system color button with the current color from the palette
            color = Gdk.RGBA(
                self.color_pallet[self.color_counter][0],
                self.color_pallet[self.color_counter][1],
                self.color_pallet[self.color_counter][2],
                self.color_pallet[self.color_counter][3],
            )
            self.color_button = self.builder.get_object('button_color')
            self.color_button.set_rgba(color)
            
            
            tag = self.vm_session.gen_random_tag_string()
            self.builder.get_object('entry_system_tag').set_text(tag)
            
            
            # --------------------------------------------------------------------------------------------
            # Configure working folder widgets
            self.entry_working_folder = self.builder.get_object('entry_working_folder')
            self.btn_choose_folder = self.builder.get_object('btn_choose_folder')
            self.set_working_folder_path()

            # Checkbox to toggle creation of a new working folder
            self.cb_create_folder_change = self.builder.get_object('cb_create_folder')
            self.cb_create_folder_change.connect('toggled', self.on_cb_create_folder_change)
            
            self.btn_choose_folder = self.builder.get_object('btn_choose_folder')
            self.btn_choose_folder.connect('clicked', self.on_btn_choose_folder)
            # --------------------------------------------------------------------------------------------

        else:
            # If already visible, just bring the window to the front
            self.window.present()
    
    def close_window (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False

    def on_entry_system_name_change (self, widget):
        """ Function doc """
        name = self.entry_system_name.get_text()
        
        #tag  = name.replace(' ','')
        #size = len(tag)
        #
        #if size > 15:
        #    tag = tag[:15]
        #else:
        #    pass
        #
        #self.builder.get_object('entry_system_tag').set_text(tag)
        self.on_entry_widget_change()
        
    def on_name_combo_changed(self, widget):
        """ Function doc """
        #fftype = self.system_types_combo.get_active()
        
        active_iter = self.system_types_combo.get_active_iter()
        if active_iter:
            model = self.system_types_combo.get_model()
            fftype = model[active_iter][1]
            #print(model[active_iter][0])
        
        #print (fftype)
        self.files    = {
            'amber_prmtop': None,
            'charmm_par'  : [],
            'charmm_psf'  : None,
            'charmm_extra': None, 
            'prm_folder' : [],
            'coordinates' : None,
            'charges'     : None,
            }
        self.residue_liststore = Gtk.ListStore(str, str, str)
        self.treeview.set_model(self.residue_liststore)
            
        if fftype == 0: #AMBER
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()
            #self.builder.get_object('gtk_label_fftype').set_text(self.amber_txt)
            self.system_types_combo.set_tooltip_text(self.amber_txt)
            
        if fftype == 1: # "CHARMM":
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()
            #self.builder.get_object('gtk_label_fftype').set_text(self.charmm_txt)
            self.system_types_combo.set_tooltip_text(self.charmm_txt)
            
        if fftype == 10: #"GROMACS":
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()
            #self.builder.get_object('gtk_label_fftype').set_text(self.gmx_txt)
            self.system_types_combo.set_tooltip_text(self.OPLS_txt)
            
        if fftype == 2:#"OPLS":
            self.builder.get_object('gtkbox_OPLS_folderchooser').show()
            #self.builder.get_object('gtk_label_fftype').set_text(self.OPLS_txt)
            
            '''Eventually, the user may have another set of parameters, but 
            by default, the pDynamo3 parameter directories are first searched.'''
            path = os.environ.get('PDYNAMO3_PARAMETERS')
            path = os.path.join(path,'forceFields/opls/protein')
            self.builder.get_object('OPLS_folderchooserbutton').set_filename(path)
            
            
        if fftype == 3: #"pDynamo files(*.pkl,*.yaml)":
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()
            
        if fftype == 4: #"Other(*.pdb,*.xyz,*.mol2...)":
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()
            self.system_types_combo.set_tooltip_text(self.xyz_pdb_text)

        if fftype == 5:#"DYFF":
            self.builder.get_object('gtkbox_OPLS_folderchooser').show()
            #self.builder.get_object('gtk_label_fftype').set_text(self.DYFF_txt)
            #try:
            '''Eventually, the user may have another set of parameters, but 
            by default, the pDynamo3 parameter directories are first searched.'''
            path = os.environ.get('PDYNAMO3_PARAMETERS')            
            path = os.path.join(path,'forceFields/dyff/dyff-1.0')
            self.builder.get_object('OPLS_folderchooserbutton').set_filename(path)
            self.system_types_combo.set_tooltip_text(self.DYFF_txt)
                
    def filetype_parser(self, filein, systemtype):
        filetype = get_file_type(filein)
        #print (filetype, systemtype)
        if filetype in ['top', 'prmtop', 'TOP', 'PRMTOP', ]:
            if systemtype == 0:
                self.files['amber_prmtop'] = filein
                return 'amber parameters/topologies'
                
   
        
        elif filetype in ['par', 'prmtop', 'prm', 'PAR', 'PRM', 'str', 'rtf']:
            if systemtype == 1:
                self.files['charmm_par'].append(filein)
                return 'charmm parameters'
        
        
        
        
        elif filetype in ['psf', 'PSF','psfx', 'PSFX']:
            if systemtype == 1:
                self.files['charmm_psf'] = filein
                return 'charmm topologies'
        
        
        
        elif filetype in ['pdb', 'PDB','mol','MOL','mol2','MOL2', 'xyz', 'XYZ', 'crd', 'inpcrd', 'chm', 'arc']:
            #if systemtype == 1:
            self.files['coordinates'] = filein
            if filetype == 'mol2':
                atoms, bonds = read_MOL2(filein)
                self.files['charges'] = atoms['charges']
                pass
            return 'coordinates'
        
        
        
        elif filetype in ['pkl', 'PKL']:
            #if systemtype == 1:
            self.files['coordinates'] = filein
            return 'pDynamo coordinates'
        else:
            return 'unknow'

    def on_button_delete_files_clicked (self, button):
        """ Function doc """
        self.residue_liststore = Gtk.ListStore(str, str, str)
        self.treeview.set_model(self.residue_liststore)
    
    def on_button_load_files_clicked (self, button):
        """ Function doc """
        files = self.easyhybrid_main.filechooser.open(select_multiple = True)
        #print(files)

        for _file in files:
            #for res in self.VObj.chains[chain].residues:
                ##print(res.resi, res.resn, chain,  len(res.atoms) ) 
            
            active_iter = self.system_types_combo.get_active_iter()
            if active_iter:
                model = self.system_types_combo.get_model()
                #print(model[active_iter][1]) # access the value of the active item
                systemtype = model[active_iter][1]
            
            filetype = self.filetype_parser( _file, systemtype)
            self.residue_liststore.append(list([_file, filetype, 'unk' ]))
        self.treeview.set_model(self.residue_liststore)
        self.files['prm_folder'] =  self.builder.get_object('OPLS_folderchooserbutton').get_filename()
            
    def on_button_import_system_clicked (self, button):
        #print('ok_button_import_a_new_system')
        #systemtype = self.system_types_combo.get_active()
        active_iter = self.system_types_combo.get_active_iter()
        if active_iter:
            model = self.system_types_combo.get_model()
            #print(model[active_iter][1]) # access the value of the active item
            systemtype = model[active_iter][1]
        else:
            return None
        self.files['prm_folder'] =  self.builder.get_object('OPLS_folderchooserbutton').get_filename()
        #name =  self.builder.get_object('entry_system_name').get_text()
        name =  self.entry_system_name.get_text()
        
        tag  =  self.builder.get_object('entry_system_tag').get_text()
        #color  =  self.builder.get_object('button_color').get_color()
        color  =  self.builder.get_object('button_color').get_rgba()
        red   = color.red 
        green = color.green 
        blue  = color.blue  
        
        
        
        #print(self.files, systemtype, color, [red, green,blue ])

        #'''
        
        #                  Working Folder
        #wfolder  = self.folder_chooser_button.get_folder()
        if self.builder.get_object('cb_create_folder').get_active():
            wfolder = self.builder.get_object('entry_working_folder').get_text()
            
            if not os.path.exists(wfolder):
                try:
                    os.makedirs(wfolder)
                except:
                    print('Failed to create the working directory {}'.format(wfolder))
        else:
            wfolder = None
        
        '''
        self.easyhybrid_main.p_session.load_a_new_pDynamo_system_from_dict(input_files    = self.files, 
                                                                           system_type    = systemtype, 
                                                                           name           = name      ,
                                                                           tag            = tag       ,
                                                                           color          = [red, green, blue],
                                                                           working_folder = wfolder)
    
        self.close_window(button, data  = None)
        #'''
        
        #'''
        try:
            self.easyhybrid_main.p_session.load_a_new_pDynamo_system_from_dict(input_files    = self.files, 
                                                                               system_type    = systemtype, 
                                                                               name           = name      ,
                                                                               tag            = tag       ,
                                                                               color          = [red, green, blue],
                                                                               working_folder = wfolder)
        
            self.close_window(button, data  = None)
            if self.color_counter > 10:
                self.color_counter == 0 
            else:
                self.color_counter += 1
        except Exception as e:
            error_str = str(e)  # converte a mensagem de erro para string
            print("Error:", error_str)
            self.easyhybrid_main.bottom_notebook.status_teeview_add_new_item(message = 'Could not import the system.', system = None)
            simpledialog = SimpleDialog(self.easyhybrid_main)
            simpledialog.error("Could not import the system.\n\n{}".format(error_str))

    def on_entry_widget_change (self, widget = None):
        """ Function doc """
        #system_name = widget.get_text()
        self.set_working_folder_path()
    
    def set_working_folder_path (self, path = None):
        """ Function doc """
        system_name = self.entry_system_name.get_text()
        
        folder = self.vm_session.vm_config.gl_parameters["workspace_path"]
        path  = os.path.join(folder, system_name)
        self.entry_working_folder.set_text(path)
 
    def on_cb_create_folder_change (self, widget):
        """ Function doc """
        if self.builder.get_object('cb_create_folder').get_active():
            #self.builder.get_object('cb_create_folder').get_active(False):
            self.builder.get_object('entry_working_folder').set_sensitive(True)
            self.builder.get_object('btn_choose_folder').set_sensitive(True)
            
        else:
            #self.builder.get_object('cb_create_folder').get_active(True):
            #self.builder.get_object('entry_working_folder').get_active(True):
            self.builder.get_object('entry_working_folder').set_sensitive(False)
            self.builder.get_object('btn_choose_folder').set_sensitive(False)

    def on_btn_choose_folder (self, widget):
        """ Function doc """
        dialog = Gtk.FileChooserDialog(
            title="Please choose a folder",
            parent=self.window,
            action=Gtk.FileChooserAction.SELECT_FOLDER
        )
        
        dialog.add_buttons(
            Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
            Gtk.STOCK_OPEN, Gtk.ResponseType.OK
        )
        path = self.vm_session.vm_config.gl_parameters['startup_path']
        dialog.set_current_folder(path)

        response = dialog.run()

        if response == Gtk.ResponseType.OK:
            self.folder = dialog.get_filename()
            self.set_working_folder_path()

        dialog.destroy()

    def update (self):
        """ Function doc """
        #print('VismolGoToAtomWindow2 update')
        pass
