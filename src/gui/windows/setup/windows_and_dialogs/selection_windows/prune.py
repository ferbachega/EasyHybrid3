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


class EasyHybridDialogPrune:
    """ Class doc """

    def __init__ (self, main         = None    , 
                        num_of_atoms = 0       ,  
                        name         = 'Unknow', 
                        tag          = 'UNK'   ,
                        e_id         = 1       , 
                        _type        = 0       ):
        
        """ Class initialiser """
        self.main = main
        self.home = main.home
        self.p_session = main.p_session
        self.vm_session = main.vm_session
        self.vobject_id = None 
        
        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(self.home,'src/gui/windows/setup/windows_and_dialogs.glade'))
        #self.builder.connect_signals(self)
        
        self.dialog       = self.builder.get_object('dialog_prune')
        
        self.label_msg    = self.builder.get_object('label_msg')
        
        
        if _type == 1: #means it's cloning
            self.dialog.set_title('Cloning Window')
            self.label_msg.set_text('The cloning process will generate a new system')
            self.builder.get_object('entry_name').set_text(name + '_cloned')
            self.builder.get_object('prune_dialog_button_prune').set_label("Clone it!")
            self.coordinates_combobox = CoordinatesComboBox(self.main.vobject_liststore_dict[e_id])
        
        else:          #means it's pruning 
            self.dialog.set_title('Pruning Window')
            self.builder.get_object('entry_name').set_text(name + '_pruned')
            self.coordinates_combobox = CoordinatesComboBox(self.main.vobject_liststore_dict[self.p_session.active_id]) 
        
        self.builder.get_object('entry_tag').set_text(tag)
        self.builder.get_object('entry_number_of_atoms').set_text(str(num_of_atoms))

        self.builder.get_object('prune_dialog_button_prune').connect("clicked", self.on_click_button_prune)
        self.builder.get_object('prune_dialog_button_cancel').connect("clicked", self.on_click_button_cancel)
        
        
        #           - - - - - - - coordinates combobox - - - - - - -
        '''--------------------------------------------------------------------------------------------'''
        self.box2 = self.builder.get_object('box_coordinates')
        self.coordinates_combobox.set_active_vobject ( pos = -1)
        self.box2.pack_start(self.coordinates_combobox, False, False, 0)
        self.box2.show_all()
        '''--------------------------------------------------------------------------------------------'''

        #print(self.coordinates_combobox)
        
        self.prune = False
        self.name  = None
        self.tag   = None
        answer = self.dialog.run()
        #print ('answer', answer)
    
    def on_click_button_prune (self, widget):
        """ Function doc """
        num_of_atoms = self.builder.get_object('entry_number_of_atoms').get_text( )
        self.name    = self.builder.get_object('entry_name').get_text( )
        self.tag     = self.builder.get_object('entry_tag').get_text( )
        
        
        self.vobject_id  = self.coordinates_combobox.get_vobject_id()
        color            = self.builder.get_object('button_color').get_rgba()
        #red   = color.red 
        #green = color.green 
        #blue  = color.blue 
        self.color = [color.red, color.green, color.blue]

        
        self.prune   = True

        self.dialog.destroy()

    def on_click_button_cancel (self, widget):
        """ Function doc """
        self.dialog.destroy()
        self.prune   = False

    def update (self):
        """ Function doc """
        #print('VismolGoToAtomWindow2 update')
        pass
