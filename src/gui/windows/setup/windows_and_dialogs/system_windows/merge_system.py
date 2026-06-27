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


VISMOL_HOME      = os.environ.get('VISMOL_HOME')
HOME             = os.environ.get('HOME')
PDYNAMO3_SCRATCH = os.environ.get('PDYNAMO3_SCRATCH')


class MergeSystemWindow(Gtk.Window):
    """ 
    A GTK3 window class to facilitate merging of two molecular or simulation systems.
    This window allows the user to select two systems, choose coordinate sets, and
    merge them with specific metadata like name, tag, and color.
    """
    
    def __init__(self, main=None):
        """
        Initialize the MergeSystemWindow instance.

        Parameters:
        main : object
            Reference to the main application object, providing access to the home
            directory, session, and other application-wide resources.
        """
        self.main = main
        self.home = main.home  # Home directory for UI files and resources
        self.p_session = main.p_session  # The session responsible for performing merges
        self.Visible = False  # Tracks whether the window is currently visible
        # ListStore holds system selections: first column for selection state, second for system name
        self.liststore = Gtk.ListStore(bool, str)
        self.selected_system_id = None  # Optional: preselected system for convenience
    
    def open_window (self, system_id = None):
        """
        Open and display the MergeSystem window. If already open, bring it to the front.

        Parameters:
        system_id : int or None
            Optional system ID to preselect in the first system combo box.
        """
        if self.Visible  ==  False:
            self.builder = self.main.builder #Gtk.Builder()
            
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'src/gui/windows/setup/merge_system.glade'))
            self.builder.connect_signals(self)

            self.window = self.builder.get_object('window')
            self.window.set_title('Merge System')
            self.window.set_keep_above(True)
            
            
            # - - - - - - - systems combobox - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            self.box1 = self.builder.get_object('box_system1')
            self.combobox_systems1 = SystemComboBox(self.main)
            self.combobox_systems1.index = 1
            '''--------------------------------------------------------------------------------------------'''
            self.box1.pack_start(self.combobox_systems1, False, False, 0)
            '''--------------------------------------------------------------------------------------------'''
            if self.selected_system_id:
                # Preselect a system if provided
                self.combobox_systems1.set_active_system (e_id = self.selected_system_id)
            # Connect the combobox change signal to update corresponding coordinates
            self.combobox_systems1.connect("changed", self.on_combobox_systems_changed)        
            
            
            # - - - - - - - systems combobox - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            self.box2 = self.builder.get_object('box_system2')
            self.combobox_systems2 = SystemComboBox(self.main)
            self.combobox_systems2.index = 2
            '''--------------------------------------------------------------------------------------------'''
            self.box2.pack_start(self.combobox_systems2, False, False, 0)
            '''--------------------------------------------------------------------------------------------'''
            self.combobox_systems2.connect("changed", self.on_combobox_systems_changed)        

            #------------------------------------------------------------------#
            self.box_coordinates1 = self.builder.get_object('box_coordinates1')
            self.coordinates_combobox1 = CoordinatesComboBox(self.main.vobject_liststore_dict[self.selected_system_id])            
            self.box_coordinates1.pack_start(self.coordinates_combobox1, False, False, 0)
            #------------------------------------------------------------------#
            self.coordinates_combobox1.index = 1

            #------------------------------------------------------------------#
            self.box_coordinates2 = self.builder.get_object('box_coordinates2')
            self.coordinates_combobox2 = CoordinatesComboBox(self.main.vobject_liststore_dict[self.selected_system_id])            
            self.box_coordinates2.pack_start(self.coordinates_combobox2, False, False, 0)
            #------------------------------------------------------------------#
            self.coordinates_combobox2.index2 = 2

            self.button_ok     = self.builder.get_object('button_merge')
            self.button_ok.connect("clicked", self.merge)

            self.button_cancel = self.builder.get_object('button_cancel')
            self.button_cancel.connect("clicked", self.close_window)

            # Display the fully constructed window
            self.window.show_all()
            self.Visible  = True   
        
        else:
            self.window.present()
    
    def close_window (self, button, data  = None):
        """
        Close and destroy the MergeSystem window, marking it as not visible.

        Parameters:
        button : Gtk.Button
            The button that triggered this closure.
        data : any
            Optional extra data (not used here).
        """
        self.window.destroy()
        self.Visible    =  False

    def on_combobox_systems_changed(self, widget):
        """
        Callback function to update coordinate comboboxes when the selected system changes.

        Parameters:
        widget : SystemComboBox
            The system combobox that triggered the event.
        """
        e_id = widget.get_system_id()
        if widget.index == 1:
            # Update coordinates combobox for system 1
            self.coordinates_combobox1.set_model(self.main.vobject_liststore_dict[e_id])
            self.coordinates_combobox1.set_active_vobject(-1)
        elif widget.index == 2:
            # Update coordinates combobox for system 2
            self.coordinates_combobox2.set_model(self.main.vobject_liststore_dict[e_id])
            self.coordinates_combobox2.set_active_vobject(-1)

    def merge(self, widget):
        """
        Execute the merge operation using selected systems, coordinates, and user-defined parameters.

        Parameters:
        widget : Gtk.Button
            The button that triggered the merge.
        """
        system1_e_id = self.combobox_systems1.get_system_id()
        system2_e_id = self.combobox_systems2.get_system_id()
        
        vobject1_id = self.coordinates_combobox1.get_vobject_id()
        vobject2_id = self.coordinates_combobox2.get_vobject_id()
        
        # Retrieve metadata from UI entries
        name = self.builder.get_object('entry_name').get_text()
        tag = self.builder.get_object('entry_tag').get_text()
        color = self.builder.get_object('button_color').get_rgba()
        
        # Extract RGB components
        red = color.red
        green = color.green
        blue = color.blue

        # Call the session's merge function with all parameters
        self.p_session.merge_system(
            e_id1=system1_e_id,
            e_id2=system2_e_id,
            vob_id1=vobject1_id,
            vob_id2=vobject2_id,
            name=name,
            tag=tag,
            color=[red, green, blue]
        )
        # Close the window after merging
        self.Visible = False
        self.window.destroy()
