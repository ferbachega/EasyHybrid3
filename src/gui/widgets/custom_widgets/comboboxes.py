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
from gi.repository import Gtk
from gi.repository import Gdk
from gi.repository import GdkPixbuf
import os

import threading
import time

from util.geometric_analysis            import get_distance 
from util.geometric_analysis            import get_dihedral 
from util.geometric_analysis            import get_angle 
#from util.periodic_table import atomic_dic 

from util.sequence_plot import GtkSequenceViewer
from pprint import pprint

import numpy as np

class SystemComboBox(Gtk.ComboBox):
    """
    A custom GTK3 ComboBox widget for selecting a system from the application's
    available systems. Supports optional linkage to a coordinates combobox.
    """
    
    def __init__(self, main=None, coord_combobox=False):
        """
        Initialize the SystemComboBox instance.

        Parameters:
        main : object
            Reference to the main application object, which provides access
            to the global system liststore and other resources.
        coord_combobox : bool
            Flag indicating whether this combobox is linked to a coordinates
            combobox for updating available coordinate sets when the system changes.
        """
        # Initialize base Gtk.ComboBox
        Gtk.ComboBox.__init__(self)
        
        self.main = main  # Reference to main application object
        
        # Set the model for the combobox using the main application's system liststore
        self.system_liststore = self.main.system_liststore
        self.set_model(self.system_liststore)
        
        # -------------------- RENDERERS --------------------
        # Add a pixbuf renderer to display system icons (e.g., molecule images)
        renderer_pixbuf = Gtk.CellRendererPixbuf()
        self.pack_start(renderer_pixbuf, True)
        self.add_attribute(renderer_pixbuf, "pixbuf", 2)  # Column 2 holds icon data
        
        # Add a text renderer to display the system's name
        renderer_text2 = Gtk.CellRendererText()
        self.pack_start(renderer_text2, True)
        self.add_attribute(renderer_text2, "text", 0)  # Column 0 holds system names
        
        # -------------------- SIGNALS --------------------
        # Connect the "changed" signal to update linked coordinate comboboxes if needed
        self.connect("changed", self.on_change)
        
        # Store whether this combobox should update a coordinates combobox
        self.coord_combobox = coord_combobox
    
    def get_system_id(self, widget = None):
        _, system_id, pixbuf = self._get_system()
        return system_id
    
    def get_system_name(self, widget = None):
        name, system_id, pixbuf = self._get_system()
        return name
        
    def _get_system(self, widget = None):
        """ Function doc """
        index = self.get_active()
        
        if index == -1:
            '''_id = -1 means no item inside the combobox'''
            return None, None, None
        
        else:    
            _, system_id, pixbuf = self.system_liststore[index]
            #print(_, system_id, pixbuf)
            return _, system_id, pixbuf

    def set_active_system (self, e_id = 0):
        """ Function doc """
        
        counter    = 0
        set_active = 0
        for key , system in self.main.p_session.psystem.items():
            if system:
                if e_id == int(key):
                    set_active = counter
                    counter += 1
                else:
                    counter += 1
            else:
                pass
        
        self.set_active(set_active)

    def on_change (self, widget):
        """ Function doc """
        #print('AQUI', self.coord_combobox)
        
        if self.coord_combobox:
            system_id = self.get_system_id()
            system  = self.main.p_session.get_system(system_id)
            
            if system_id is not None:
                self.coord_combobox.set_model(self.main.vobject_liststore_dict[system_id])
                #self.refresh_selection_liststore (system_id)            
                size  =  len(list(self.main.vobject_liststore_dict[system_id]))
                self.coord_combobox.set_active(size-1)
        else:
            pass


class CoordinatesComboBox(Gtk.ComboBox):
    """
    A custom GTK3 ComboBox widget for selecting a coordinate set (vobject) associated
    with a molecular or simulation system. Supports optional pixbuf icons for visual
    representation.
    """
    
    def __init__(self, coordinates_liststore=None, system_combobox=None, pixbuf=True):
        """
        Initialize the CoordinatesComboBox instance.

        Parameters:
        coordinates_liststore : Gtk.ListStore
            The ListStore containing coordinate sets for the selected system.
            Typically each row holds metadata about a coordinate set.
        system_combobox : SystemComboBox or None
            Optional reference to the system combobox this coordinates combobox is linked to.
            Can be used to update the coordinates dynamically when the system changes.
        pixbuf : bool
            Flag to determine whether to display icons in the combobox using a CellRendererPixbuf.
        """
        # Initialize base Gtk.ComboBox
        Gtk.ComboBox.__init__(self)
        
        # Set the ListStore containing coordinate sets as the model for this combobox
        self.coordinates_liststore = coordinates_liststore
        self.set_model(self.coordinates_liststore)
        
        # -------------------- RENDERERS --------------------
        if pixbuf:
            # Add an optional pixbuf renderer to display a small image for each coordinate set
            renderer_pixbuf = Gtk.CellRendererPixbuf()
            self.pack_start(renderer_pixbuf, True)
            self.add_attribute(renderer_pixbuf, "pixbuf", 3)  # Column 3 holds icon data
        
        # Add a text renderer to display the coordinate set name
        renderer_text = Gtk.CellRendererText()
        self.pack_start(renderer_text, True)
        self.add_attribute(renderer_text, "text", 0)  # Column 0 holds coordinate names

    def get_vobject(self):
        """
        Retrieve the currently selected vobject (coordinate object) from the combobox.

        Returns:
        tuple
            Information about the selected coordinate set, typically including:
            - name: The display name of the vobject
            - key1: Primary identifier (e.g., system ID)
            - key2: Secondary identifier (e.g., coordinate set ID)
            - pixbuf: Optional icon associated with the vobject
        """
        # Fetch all coordinate info from the selected row
        name, key1, key2, pixbuf = self._get_coordinates_info()
        
        # Placeholder for linking to actual vobject in VM session
        # self.vobject = self.main.vm_session.vm_objects_dic[vobject_index]
        
    def get_vobject_id (self):
        """ Returns the id of the vobject (id number of the vobject - 
        generated in eSession.py / self._add_vismol_object() ) """
        name, key1,  key2, pixbuf  = self._get_coordinates_info ()
        return key1

    
    def get_system_id (self):
        """ Returns the id of the vobject (id number of the vobject - 
        generated in eSession.py / self._add_vismol_object() ) """
        name, key1,  key2, pixbuf  = self._get_coordinates_info ()
        return key2
        
        
    def _get_coordinates_info (self):
        """ Function doc """
        cb_id =  self.get_active()
        model = self.get_model()
        name, key1,  key2, pixbuf = model[cb_id]
        #print (name, key1,  key2)
        return name, key1,  key2, pixbuf


    def set_active_vobject (self, pos = 0):
        """ Function doc """
        #self.coordinates_liststore[cb_id]
        if pos ==-1:
            n = len(self.get_model())
            #print(776, n)
            self.set_active(n-1)
        else:
            self.set_active(pos)

    
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
