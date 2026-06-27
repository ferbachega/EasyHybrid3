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


class TextWindow:
    """ Class doc """
    
    def __init__ (self, text = 'No text', title = None):
        """ Class initialiser """
        self.window = Gtk.Window(title=title)
        self.window.set_default_size(1100, 600)
        
        #if title:
        #    self.window.title = title
        
        
        self.textview = Gtk.TextView()
        self.textbuffer = self.textview.get_buffer()
        self.textbuffer.set_text(text)
        
        
        # Create a Pango font description with the desired font family and size
        fontdesc = Pango.FontDescription()
        fontdesc.set_family("Monospace")
        fontdesc.set_size(12 * Pango.SCALE)  # 12 point size
        
        # Apply the font description to the text view
        self.textview.modify_font(fontdesc)
        
        # Set the text color to black
        style = self.textview.get_style_context()
        style.add_class("text-black")
        
        
        scrolledwindow = Gtk.ScrolledWindow()
        scrolledwindow.set_hexpand(True)
        scrolledwindow.set_vexpand(True)
        scrolledwindow.add(self.textview)
        
        self.window.add(scrolledwindow)
        self.window.show_all()


class InfoWindow:
    """ 
    Create a text window. Currently used to display system 
    information or log output. 
    """
    
    def __init__ (self, system = None, text = None):
        """ Class initialiser """
        
        if text:
            pass
        else:
            text = ''
        
        
        if system:
            log  = LogFile(system)
            path = log.path
            with open(path, "r") as f:
                header = f.read()
        else:
            header = '' 
        
        text = header+text
        textwindow = TextWindow(text)
