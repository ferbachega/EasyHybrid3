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
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk
import cairo
import random
import numpy as np
import sys

try:
    from util.colormaps import COLOR_MAPS  
except:
    from colormaps import COLOR_MAPS
from pprint import pprint

class Canvas(Gtk.DrawingArea):
    """ Class doc """
    
    def __init__ (self, bg_color = [1,1,1], cmap =  'jet' ):
        """ Class initialiser """
        super().__init__( )
        #Gtk.DrawingArea.__init__(main)
        
        self.color_map = COLOR_MAPS[cmap]
        
        self.bg_color = None
        self.set_bg_color (bg_color[0], bg_color[1], bg_color[2])
        
        self.set_events(Gdk.EventMask.ALL_EVENTS_MASK)
        self.connect("draw", self.on_draw)

        
    def set_bg_color (self, r = 1 , g = 1 , b = 1):
        """ Function doc """
        self.bg_color = [r, g, b]
    
    
    def _make_bg (self, cr):
        """ Function doc """
        cr.set_source_rgb( self.bg_color[0],self.bg_color[1], self.bg_color[2])
        cr.paint()
    
    
    def get_i_and_j_from_click_event (self, event):
        """ Function doc """
        (x, y) = int(event.x), int(event.y)
        
        x_on_plot = int((x-self.bx)/self.factor_x)
        y_on_plot = int((y-self.by)/self.factor_y)
        
        return x_on_plot, y_on_plot, x, y

    
    def check_clicked_points(self, x_on_plot, y_on_plot):
        """ Function doc """
        
        if int(x_on_plot) < 0 or int(x_on_plot) >= self.size_x :
            print('canvas')
            self.points = []
            self.selected_dot = None
            return False
        elif int(y_on_plot) < 0 or int(y_on_plot) >= self.size_y :
            print('canvas')
            self.points = []
            self.selected_dot = None #set no None the redot at IJ matrix
            return False
        
        else:
            return True
            
            self.points.append((x_on_plot, y_on_plot))
        
        print("Mouse clicker at:",  x, y, int(x_on_plot), int(y_on_plot), 
                                    (x-self.bx)/self.factor_x, #/(self.x_final-self.bx), 
                                    (y-self.by)/self.factor_y) #/(self.y_final-self.by) ) #, x - self.x_final , y - self.y_final , self.x_final ,  self.y_final ,   (x-self.bx) /self.factor_x ,    (y-self.by)/self.factor_y )
        print(self.points)

    
    def on_mouse_button_press(self, widget, event):
        #(x, y) = int(event.x), int(event.y)
        ##x, y = device.get_position(widget)
        #print("Mouse moved to:", x, y)
        (x, y) = int(event.x), int(event.y)
        
        x_on_plot = int((x-self.bx)/self.factor_x)
        y_on_plot = int((y-self.by)/self.factor_y)
        

        
        if int(x_on_plot) < 0 or int(x_on_plot) >= self.size_x :
            #print('canvas2222')
            self.points = []
            self.selected_dot = None
        
        elif int(y_on_plot) < 0 or int(y_on_plot) >= self.size_y :
            #print('canvas1111')
            self.points = []
            self.selected_dot = None
        
        else:
            self.points.append((x_on_plot, y_on_plot))
        
        print("Mouse clicker at:",  x, y, int(x_on_plot), int(y_on_plot), 
                                    (x-self.bx)/self.factor_x, #/(self.x_final-self.bx), 
                                    (y-self.by)/self.factor_y) #/(self.y_final-self.by) ) #, x - self.x_final , y - self.y_final , self.x_final ,  self.y_final ,   (x-self.bx) /self.factor_x ,    (y-self.by)/self.factor_y )
        print(self.points)
        self.queue_draw()

    
    def on_motion(self, widget, event):
        (x, y) = int(event.x), int(event.y)
        x_on_plot = int((x-self.bx)/self.factor_x)
        y_on_plot = int((y-self.by)/self.factor_y)
        #print("Mouse moved to:", x, y, int(x_on_plot), int(y_on_plot))

    def get_pixel_rgb(self, x, y):
        # Converte os dados da superfície para um array de bytes
        buf = self.surface.get_data()

        # Cairo usa um formato ARGB32, onde cada pixel é representado por 4 bytes (A, R, G, B)
        array = np.frombuffer(buf, dtype=np.uint8)
        array = array.reshape((self.height, self.width, 4))  # Reshape para uma matriz 2D com 4 valores por pixel

        # Pegar o valor do pixel em x, y
        pixel = array[y, x]

        # O formato é ARGB, então extraímos apenas R, G e B
        a, r, g, b = pixel
        return (r, g, b)
