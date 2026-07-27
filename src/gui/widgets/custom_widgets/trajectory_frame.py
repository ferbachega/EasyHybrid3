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

class VismolTrajectoryFrame(Gtk.Frame):
    """ Class doc """
    
    def __init__ (self, vm_session = None):
        """ Class initialiser """
        self.vm_session = vm_session 
        self.frame_forward = True
        self.frame      = Gtk.Frame()
        #self.frame.set_shadow_type(Gtk.SHADOW_IN)
        self.frame.set_border_width(4)
        
        self.box        = Gtk.Box(orientation = Gtk.Orientation.VERTICAL, spacing = 6) 
        
        self.box.set_margin_top    (3)
        self.box.set_margin_bottom (3)
        self.box.set_margin_left   (3)
        self.box.set_margin_right  (3)
        
        self.value      = 0
        self.upper      = 100
        self.scale      = Gtk.Scale()
        self.running    = False
        #self.adjustment = Gtk.Adjustment(self.value, 1, 1, 0, 1, 0)
        self.adjustment     = Gtk.Adjustment(value         = self.value,
                                             lower         = 0,
                                             upper         = 100,
                                             step_increment= 1,
                                             page_increment= 1,
                                             page_size     = 1)
        
        
        
        self.scale.set_adjustment ( self.adjustment)
        self.scale.set_digits(0)
        #self.scale.connect("change_value", self.on_scaler_frame_change_change_value)
        self.scale.connect("value_changed", self.on_scaler_frame_value_changed)
        self.box.pack_start(self.scale, True, True, 0)
        
        self.vbox =  Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
        
        self.functions = [self.reverse, self.stop, self.play, self.forward ]
        c = 0
        for label in ['<<','#', '>','>>']:
            button = Gtk.Button(label)
            self.vbox.pack_start(button, True, True, 0)
            
            if self.functions[c]:
                button.connect("clicked", self.functions[c])
            c += 1 
            
        self.box.pack_start(self.vbox, True, True, 0)
        
        
        #----------------------------------------------------------------------------
        '''
        self.label2 =  Gtk.Label('Object:')
        self.combobox_vobjects = Gtk.ComboBox()
        #self.combobox_vobjects = Gtk.ComboBox.new_with_model(self.vm_session.Vismol_Objects_ListStore)
        
        self.combobox_vobjects.connect("changed", self.on_combobox_vobjects_changed)
        self.renderer_text = Gtk.CellRendererText()
        self.combobox_vobjects.pack_start(self.renderer_text, True)
        self.combobox_vobjects.add_attribute(self.renderer_text, "text", 1)
        '''
        #self.box.pack_start(self.combobox, True, True, 0)
        #----------------------------------------------------------------------------
        
        #----------------------------------------------------------------------------
        self.vbox2 = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
        self.label = Gtk.Label('FPS:')
        self.label_step = Gtk.Label('Step:')
        
        #self.entry = Gtk.Entry()
        #self.entry.set_text(str(25))
        
        self.fps_adjustment = Gtk.Adjustment(value          = 24 , 
                                             upper          = 10000, 
                                             step_increment = 1  , 
                                             page_increment = 10 )
        
        self.step_adjustment = Gtk.Adjustment(value          = 1 , 
                                             upper          = 10000, 
                                             step_increment = 1  , 
                                             page_increment = 10 )
        
        
        #self.fps_adjustment = Gtk.Adjustment(value         = float, 
        #                                     lower         = float, 
        #                                     upper         = float, 
        #                                     step_increment= float, 
        #                                     page_increment= float, 
        #                                     page_size     = float)
        self.entry = Gtk.SpinButton()
        self.entry.set_adjustment ( self.fps_adjustment)
        
        self.entry_step = Gtk.SpinButton()
        self.entry_step.set_adjustment ( self.step_adjustment)
        
        #self.vbox2.pack_start(self.label2, True, True, 0)
        #self.vbox2.pack_start(self.combobox_vobjects, True, True, 0)
        self.vbox2.pack_start(self.label, True, True, 0)
        self.vbox2.pack_start(self.entry, True, True, 0)
        
        self.vbox2.pack_start(self.label_step, True, True, 0)
        self.vbox2.pack_start(self.entry_step, True, True, 0)
        
        self.box.pack_start(self.vbox2, True, True, 0)
        #----------------------------------------------------------------------------

    def play (self, button):
        # Create and start a new thread when the button is clicked
        if self.running:
            pass
        
        else:
            #thread = threading.Thread(target=self.long_task)
            thread = threading.Thread(target=self.long_task_boomerang)
            thread.start()
            
            self.running = True
        
    def long_task_boomerang(self):
        self.stop_thread = False
        i = 0
        
        while self.stop_thread == False:
            if self.stop_thread:
                return
            
            if self.frame_forward:
            
                value = self.forward(step =  self.entry_step.get_value()) 
                #value = self.reverse(None) 

                if value >= self.upper-1:
                    #value = 0
                    self.frame_forward = False
                    self.scale.set_value(int(value))
                    self.vm_session.set_frame(int(value))
                #time.sleep(0.01)
                time.sleep(1/self.entry.get_value())
            
            else:
                value = self.reverse(step =  self.entry_step.get_value()) 
                if value <= 0:
                    self.frame_forward = True
                    self.scale.set_value(int(value))
                    self.vm_session.set_frame(int(value))
                time.sleep(1/self.entry.get_value())
    
    def long_task(self):
        # This flag is used to stop the thread
        self.stop_thread = False
        i = 0
        while self.stop_thread == False:
            if self.stop_thread:
                return
            value = self.forward(None) 
            #value = self.reverse(None) 

            if value >= self.upper-1:
                value = 0
                self.scale.set_value(int(value))
                self.vm_session.set_frame(int(value))
            #time.sleep(0.01)
            time.sleep(1/self.entry.get_value())
            
    def stop (self, button):
        """ Function doc """
        self.stop_thread = True
        self.running = False    

    def forward (self, button = None, step = 1):
        """ Function doc """
        value =  int(self.scale.get_value())
        value = value+step
        self.scale.set_value(int(value))
        self.vm_session.set_frame(int(value))
        #print(value)
        return value
        
        
    def reverse (self, button =  None, step = 1):
        """ Function doc """
        value = int(self.scale.get_value())
        #print(value)
        if value == 0:
            pass
        else:
            value = value-step
            if value < 0:
                value = 0
            else: 
                pass
                
        self.vm_session.set_frame(int(value))
        self.scale.set_value(value)
        return value
        #print(value)
    
    def get_box (self):
        """ Function doc """
        #self.add(self.box)
        return self.box
        #return self.frame
        
    def on_combobox_vobjects_changed (self, widget):
        """ Function doc """
        #print('\n\n',widget)
        #print('\n\n',widget.get_active())
        
        cb_index = widget.get_active()
        if cb_index in self.vm_session.vobjects_dic:
            self.VObj = self.vm_session.vobjects_dic[widget.get_active()]
            #self.VObj = self.vm_session.vobjects[widget.get_active()]
            number_of_frames = len(self.VObj.frames)
            self.scale.set_range(0, int(number_of_frames))
            self.scale.set_value(self.vm_session.get_frame())
        else:
            pass

    def on_scaler_frame_value_changed (self, hscale, text= None,  data=None):
        """ Function doc """
        value = hscale.get_value()
        pos   = hscale.get_value_pos ()
        self.vm_session.set_frame(int(value)) 
        self.scale.set_value(value)
        #print(value, pos)

    def update (self):
        """ Function doc """
        #print('VismolTrajectoryFrame update')
        #for index , vobject in self.vm_session.vobjects_dic.items():
        #last_obj = len(self.vm_session.vobjects) -1
        last_obj = len(self.vm_session.vobjects_dic.items()) -1
        self.combobox_vobjects.set_active(last_obj)
    
    def change_range (self, upper = 100):
        """ Function doc """
        #print('upper =', upper)
        self.adjustment     = Gtk.Adjustment(value         = self.value,
                                             lower         = 0,
                                             upper         = upper,
                                             step_increment= 1,
                                             page_increment= 1,
                                             page_size     = 1)
        self.scale.set_adjustment ( self.adjustment)
        self.scale.set_digits(0)
        self.upper = upper
