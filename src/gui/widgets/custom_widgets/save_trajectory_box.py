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

# --- imports entre modulos adicionados na refatoracao ---
from gui.widgets.custom_widgets.file_choosers import FolderChooserButton

class SaveTrajectoryBox:
    """ Class doc """
    
    def __init__ (self, parent, home):
        """ Class initialiser """
        
        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(home,'src/gui/widgets/trajectory_box.glade'))
        self.builder.connect_signals(self)
        
        self.box = self.builder.get_object('trajectory_box')
        
        self.folder_chooser_button = FolderChooserButton(main =  parent, home = home)
        self.builder.get_object('folder_chooser_box').pack_start(self.folder_chooser_button.btn, True, True, 0)
        self.folder_chooser_button.btn.set_sensitive(False)
                
        '''--------------------------------------------------------------------------------------------'''
        self.format_store = Gtk.ListStore(str)
        formats = [
                    "pDynamo / pkl" ,
                    #"amber / crd"   ,
                    #"charmm / dcd"  ,
                    #"xyz"           ,
            ]
        for format in formats:
            self.format_store.append([format])
            #print (format)
        self.formats_combo = self.builder.get_object('combobox_format')
        self.formats_combo.set_model(self.format_store)
        #self.formats_combo.connect("changed", self.on_name_combo_changed)
        self.formats_combo.set_model(self.format_store)
        
        renderer_text = Gtk.CellRendererText()
        self.formats_combo.pack_start(renderer_text, True)
        self.formats_combo.add_attribute(renderer_text, "text", 0)
        '''--------------------------------------------------------------------------------------------'''
        self.formats_combo.set_active(0)
        self.set_folder(None)
        #simParameters["trajectory_name"] = self.save_trajectory_box.builder.get_object('entry_trajectory_name').get_text()
    #====================================================================================
    def on_toggle_save_checkbox (self, widget):
        """ Function doc """
        if self.builder.get_object('checkbox_save_traj').get_active():
            self.builder.get_object('entry_trajectory_name').set_sensitive(True)
            
            self.builder.get_object('label_working_folder').set_sensitive(True)
            #self.builder.get_object('file_chooser_working_folder').set_sensitive(True)
            self.builder.get_object('label_format').set_sensitive(True)
            self.builder.get_object('combobox_format').set_sensitive(True)
            self.builder.get_object('label_trajectory_frequency').set_sensitive(True)
            self.builder.get_object('entry_trajectory_frequency').set_sensitive(True)
            self.folder_chooser_button.btn.set_sensitive(True)

        else:
            self.builder.get_object('entry_trajectory_name').set_sensitive(False)

            self.builder.get_object('label_working_folder').set_sensitive(False)
            #self.builder.get_object('file_chooser_working_folder').set_sensitive(False)
            self.builder.get_object('label_format').set_sensitive(False)
            self.builder.get_object('combobox_format').set_sensitive(False)
            self.builder.get_object('label_trajectory_frequency').set_sensitive(False)
            self.builder.get_object('entry_trajectory_frequency').set_sensitive(False)
            self.folder_chooser_button.btn.set_sensitive(False)
       
    
    def get_trajectory_frequency (self):
        """ Function doc """
        return int(self.builder.get_object('entry_trajectory_frequency').get_text())
    
    def set_trajectory_frequency (self, freq = 10):
        """ Function doc """
        self.builder.get_object('entry_trajectory_frequency').set_text(str(freq))
    
    def get_format (self):
        return  self.formats_combo.get_active()
    
    def set_format (self, _format):
        return  self.formats_combo.set_active(_format)
        
    def get_active (self):
        """ Function doc  """
        if self.builder.get_object('checkbox_save_traj').get_active():
            return True
        else:
            return False
    
    def set_active (self, active = True):
        """ Function doc  """
        self.builder.get_object('checkbox_save_traj').set_active(active) 

    #====================================================================================
    def get_folder (self):
        """ Function doc """
        return self.folder_chooser_button.get_folder()
    #====================================================================================    
    def set_folder (self, folder):
        """ Function doc """
        self.folder_chooser_button.set_folder(folder)
    
    def get_filename (self):
        """ Function doc """
        return self.builder.get_object('entry_trajectory_name').get_text()
    
    def set_filename (self, filename = 'filename'):
        """ Function doc """
        return self.builder.get_object('entry_trajectory_name').set_text(filename)
