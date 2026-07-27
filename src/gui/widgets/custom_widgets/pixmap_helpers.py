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

def getColouredPixmap( r, g, b, a=255 ):
    """ Given components, return a colour swatch pixmap """
    CHANNEL_BITS=8
    WIDTH=10
    HEIGHT=10
    swatch = GdkPixbuf.Pixbuf.new( GdkPixbuf.Colorspace.RGB, True, CHANNEL_BITS, WIDTH, HEIGHT ) 
    swatch.fill( (r<<24) | (g<<16) | (b<<8) | a ) # RGBA
    return swatch


def get_colorful_square_pixel_buffer (system = None,  atomic_symbol = 'C'):
    """ Function doc """
    if system is not None:
        color        =  system.e_color_palette[atomic_symbol]
        res_color    = [int(color[0]*255),int(color[1]*255),int(color[2]*255)] 
        pixelbuffer  =  getColouredPixmap( res_color[0], res_color[1], res_color[2] )
    
    else:
        res_color    = [255,255,255] 
        pixelbuffer  =  getColouredPixmap( res_color[0], res_color[1], res_color[2] )
    return pixelbuffer
