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


class TrajectoryPlayerWindow:
    """ Class doc """
    
    def __init__ (self, main = None):
        """ Class initialiser """
        self.main = main
        self.vm_session = main.vm_session
        self.p_session  = main.p_session
            
        self.Visible = False
        self.upper = 1
        
        
    def open_window (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.window = Gtk.Window()
            self.window.connect('destroy-event', self.close_window)
            self.window.connect('delete-event', self.close_window)
            self.vm_traj_obj = VismolTrajectoryFrame(self.vm_session)
            self.vm_traj_obj.change_range(upper = self.upper)
            self.vm_traj_obj.scale.set_value(self.vm_session.frame)
            self.window.add(self.vm_traj_obj.get_box())
            
            self.window.set_default_size(300, 50)  
            self.window.set_title('EasyHybrid Player')  
            self.window.set_keep_above(True)
            
            self.window.show_all()
            self.Visible  = True
    
        else:
            self.window.present()
        

    def close_window (self, button, data  = None):
        """ Function doc """
        self.vm_traj_obj.stop(None)

        self.window.destroy()
        self.Visible    =  False
        self.main.trajectory_player_button.set_active(False)


    def change_range (self, upper = 100):
        """ Function doc """
        if self.Visible:
            print('upper =', upper)
            self.vm_traj_obj.change_range (upper = upper)
