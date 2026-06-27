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
from gui.widgets.custom_widgets.comboboxes import CoordinatesComboBox

class SequenceViewerBox(Gtk.Box):
    """ Class doc """
    
    def __init__ (self, main = None, mode = 0 ):
        """ Class initialiser """
        Gtk.Box.__init__(self)
        
        self.main = main
        self.seqview   = GtkSequenceViewer(main = main)
        #-----------------------------------------------------------------------
        self.seqview.add_events(Gdk.EventMask.POINTER_MOTION_MASK)
        self.seqview.add_events(Gdk.EventMask.BUTTON_PRESS_MASK)
        self.seqview.add_events(Gdk.EventMask.BUTTON_RELEASE_MASK)

        self.seqview.connect("motion-notify-event" , self.seqview.on_motion)
        self.seqview.connect("button_press_event"  , self.seqview.button_press)
        self.seqview.connect("button-release-event", self.seqview.button_release)
        #-----------------------------------------------------------------------

        #self.add()
        
        #--------------------------------------------------------------#
        self.box_horizontal1 = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 10)
        
        self.label1  = Gtk.Label()
        self.label1.set_text('Coordinates:')
        self.box_horizontal1.pack_start(self.label1, False, False, 0)
        self.coordinates_combobox = CoordinatesComboBox(self.main.vobject_liststore_dict[self.main.p_session.active_id])
        
        self.box_horizontal1.pack_start(self.coordinates_combobox, False, False, 0)
        #--------------------------------------------------------------#




        #------------------------------------------------------------------#
        #                  CHAIN combobox and Label
        #------------------------------------------------------------------#
        self.label2  = Gtk.Label()
        self.label2.set_text('Chain:')
        self.box_horizontal1.pack_start(self.label2, False, False, 0)

        self.liststore_chains = Gtk.ListStore(str)
        self.combobox_chains = Gtk.ComboBox.new_with_model(self.liststore_chains)
        self.combobox_chains.connect("changed", self.on_combobox_chains_changed)
        renderer_text = Gtk.CellRendererText() 
        self.combobox_chains.pack_start(renderer_text, True)
        self.combobox_chains.add_attribute(renderer_text, "text", 0)
        #vbox.pack_start(self.combobox_chains, False, False, True)
        self.box_horizontal1.pack_start(self.combobox_chains, False, False, 0)
        
        #Gtk.show_all()



    def on_combobox_chains_changed (self, widget):
        """ Function doc """
        
        tree_iter = self.combobox_chains.get_active_iter()
        if tree_iter is not None:
            model = self.combobox_chains.get_model()
            chain = model[tree_iter][0]
            #print("Selected: country=%s" % country)
        
        self.current_filter_chain = chain
        #print("%s Chain selected!" % self.current_filter_chain)
        # we update the filter, which updates in turn the view
        self.chain_filter.refilter()
