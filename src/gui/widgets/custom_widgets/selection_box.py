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

class VismolSelectionTypeBox(Gtk.Box):
    """ Class doc """
    
    def __init__ (self, vm_session = None):
        """ Class initialiser """
        Gtk.Box.__init__(self)
        #self.set_orientation(Gtk.Orientation.VERTICAL)
        #self.set_orientation(Gtk.Orientation.HORIZONTAL)
        #self.set_spacing(5)
        self.box           = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
        self.vm_session = vm_session
        #combobox
        
        self.Vismol_selection_modes_ListStore = Gtk.ListStore(str)
        data = ['atom'    , 
                'residue' ,
                'chain'   , 
                'molecule',
                'c-alpha' ,
                'solvent' ,
                #'protein' , 
                #'C alpha' ,
                #'solvent' ,
                #'atom name',
                #'element',
                ]
        for i in data:
            #print (i)
            self.Vismol_selection_modes_ListStore.append([i])
            
        self.combobox_selection_type = Gtk.ComboBox.new_with_model(self.Vismol_selection_modes_ListStore)
        
        
        self.combobox_selection_type.set_model(self.Vismol_selection_modes_ListStore)
        
        self.renderer_text = Gtk.CellRendererText()
        self.combobox_selection_type.pack_start(self.renderer_text, True)
        self.combobox_selection_type.add_attribute(self.renderer_text, "text", 0)
        
        self.combobox_selection_type.connect('changed', self.on_combobox_selection_type)
        #
        
        #labels        
        self.label_selecting_by = Gtk.Label('selecting by: ')
        
        #toggle_button
        self.toggle_button_selecting_mode = Gtk.ToggleButton('Viewing')
        self.toggle_button_selecting_mode.connect('clicked', self.on_toggle_button_selecting_mode)
        
        
        #Atom name entry box
        self.entry_atom_names = None#Gtk.Entry()
        self.entry_elements   = None
        
        # Packing 
        self.box.pack_start(self.toggle_button_selecting_mode, False, False, 0)
        self.box.pack_start(self.label_selecting_by          , False, False, 0)
        self.box.pack_start(self.combobox_selection_type     , False, False, 0)
        
        self.combobox_selection_type.set_active(1)
        self.box.show_all()
    
    def on_entry_data_change (self, widget):
        """ Function doc """
        
        string = widget.get_text()
        
        if widget == self.entry_atom_names:
            self.vm_session.selections[self.vm_session.current_selection].selected_atom_names_list = []
            selected_atom_names_list = self.vm_session.selections[self.vm_session.current_selection].selected_atom_names_list 
        
            keys = string.split('+')
            for key in keys:
                selected_atom_names_list.append(key.strip())
            
            print ('entry_atom_names', self.vm_session.selections[self.vm_session.current_selection].selected_atom_names_list)
        
        
        
        elif widget == self.entry_elements:
            self.vm_session.selections[self.vm_session.current_selection].selected_element_list = []
            selected_element_list = self.vm_session.selections[self.vm_session.current_selection].selected_element_list  
            
            keys = string.split('+')
            for key in keys:
                selected_element_list.append(key.strip())
            
            print ('entry_elements', self.vm_session.selections[self.vm_session.current_selection].selected_element_list)

        
        else: 
            pass
               
    def on_combobox_selection_type (self, combobox):
        """ Function doc """
        self.active = combobox.get_active()
        #print('on_combobox_selection_type', self.active)
        if self.active == 0:
            self.vm_session.viewing_selection_mode(sel_type = 'atom')
            #self.show_or_hide_entries (name = 'atom_names', show = False)
            #self.show_or_hide_entries (name = 'element', show = False)

        elif self.active == 1:
            self.vm_session.viewing_selection_mode(sel_type = 'residue')
            #self.show_or_hide_entries (name = 'atom_names', show = False)
            #self.show_or_hide_entries (name = 'element', show = False)
        
        elif self.active == 2:
            self.vm_session.viewing_selection_mode(sel_type = 'chain')
            #self.show_or_hide_entries (name = 'atom_names', show = False)
            #self.show_or_hide_entries (name = 'element', show = False)
        
        elif self.active == 3:
            self.vm_session.viewing_selection_mode(sel_type = 'molecule')
            #self.show_or_hide_entries (name = 'atom_names', show = False)
            #self.show_or_hide_entries (name = 'element', show = False)
        
        elif self.active == 4:
            self.vm_session.viewing_selection_mode(sel_type = 'c-alpha')
        
        elif self.active == 5:
            self.vm_session.viewing_selection_mode(sel_type = 'solvent')
            #self.show_or_hide_entries (name = 'atom_names', show = False)
            #self.show_or_hide_entries (name = 'element', show = False)
        
        #elif self.active == 5:
        #    self.vm_session.viewing_selection_mode(sel_type = 'solvent')
        #    #self.show_or_hide_entries (name = 'atom_names', show = False)
        #    #self.show_or_hide_entries (name = 'element', show = False)
        #
        #elif self.active == 6:
        #    self.vm_session.viewing_selection_mode(sel_type = 'atom name')
        #    #self.show_or_hide_entries (name = 'atom_names', show = True)
        #    #self.show_or_hide_entries (name = 'element', show = False)
        #
        #elif self.active == 7:
        #    self.vm_session.viewing_selection_mode(sel_type = 'element')
        #    #self.show_or_hide_entries (name = 'atom_names', show = False)
        #    #self.show_or_hide_entries (name = 'element', show = True)
        #        
        else:pass
        
    def change_sel_type_in_combobox (self, sel_type):
        """ Function doc """
        #print('change_sel_type_in_combobox', sel_type)
        if sel_type == 'atom':
            self.combobox_selection_type.set_active(0)
        
        elif sel_type == 'residue':
            self.combobox_selection_type.set_active(1)
        
        elif sel_type == 'chain':
            self.combobox_selection_type.set_active(2)
        
        elif sel_type == 'molecule':
            self.combobox_selection_type.set_active(3)
        
        elif sel_type == 'C alpha':
            self.combobox_selection_type.set_active(4)
        
        elif sel_type == 'solvent':
            self.combobox_selection_type.set_active(5)
        
        elif sel_type == 'atom name':
            self.combobox_selection_type.set_active(6)
        
        elif sel_type == 'element':
            self.combobox_selection_type.set_active(7)
        
        else: pass
        
    def on_toggle_button_selecting_mode (self, button):
        """ Function doc """
        if button.get_active():
            state = "on"
            
            self.vm_session.picking_selection_mode = True
            button.set_label('Picking')
            #print(self.combobox_selection_type.get_active())
            self.vm_session._selection_function_set (None)
            self.vm_session.vm_glcore.queue_draw()
            
            #self.vm_session._picking_selection_mode = True
            
            self.combobox_selection_type.set_sensitive(False)
            self.label_selecting_by.set_sensitive(False)
            
            circle = Gdk.Cursor(Gdk.CursorType.CROSSHAIR)
            button.get_window().set_cursor(circle)
            
        else:
            state = "off"
            self.vm_session.picking_selection_mode = False
            button.set_label('Viewing')
            self.vm_session.vm_glcore.queue_draw()
            
            self.combobox_selection_type.set_sensitive(True)
            self.label_selecting_by.set_sensitive(True)
            
            
            circle = Gdk.Cursor(Gdk.CursorType.ARROW)
            button.get_window().set_cursor(circle)
            
    def change_toggle_button_selecting_mode_status (self, status = False):
        """ Function doc """
        self.toggle_button_selecting_mode.set_active(status)

    def update (self):
        """ Function doc """
        print('VismolSelectionTypeBox update')
