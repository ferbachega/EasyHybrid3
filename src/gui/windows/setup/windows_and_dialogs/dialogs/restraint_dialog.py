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


class AddHarmonicRestraintDialog:
    """ 
    A GTK-based dialog for adding or editing a harmonic restraint between two atoms.  
    The dialog allows the user to specify atom indices, names, target distance, 
    and force constant for the restraint. 
    """
    
    def __init__(self, main=None, 
                       atom1=None, 
                       atom2=None, 
                       distance=0.0, 
                       force=4000,
                       system_id=None,
                       edit=False):
        """ 
        Initialize the dialog. 
        
        Parameters:
        - main: reference to the main application object. 
        - atom1, atom2: atom objects or indices defining the restraint. 
        - distance: target distance for the harmonic restraint (default 0.0). 
        - force: force constant for the restraint (default 4000). 
        - system_id: (optional) ID of the system in use. 
        - edit: if True, load the dialog with pre-filled values for editing. 
        """
        self.main = main
        self.home = main.home
        
        # Load GUI layout from Glade file
        self.builder = Gtk.Builder()
        self.builder.add_from_file(
            os.path.join(self.home, 'src/gui/windows/setup/add_harmonic_restraint_dialog.glade')
        )
        self.builder.connect_signals(self)
        self.builder.get_object('dialog').connect('destroy', self.close_window)
        
        # If editing an existing restraint, populate fields with provided values
        if edit:
            self.builder.get_object('button_add').set_label('Ok')  # Change button label
            self.builder.get_object('entry_atom1_index_coord1').set_text(str(atom1))
            self.builder.get_object('entry_atom2_index_coord1').set_text(str(atom2))
            # Optional: could also use atom names if available
            self.builder.get_object('entry_dmin_coord1').set_text(str(distance))
            self.builder.get_object('entry_FORCE_coord1').set_text(str(force))
        else:
            # If adding a new restraint, populate with default atom data
            self.builder.get_object('entry_atom1_index_coord1').set_text(str(atom1.index-1))
            self.builder.get_object('entry_atom2_index_coord1').set_text(str(atom2.index-1))
            self.builder.get_object('entry_atom1_name_coord1').set_text(atom1.name)
            self.builder.get_object('entry_atom2_name_coord1').set_text(atom2.name)
            self.builder.get_object('entry_dmin_coord1').set_text(str(distance))
            self.builder.get_object('entry_FORCE_coord1').set_text(str(force))

        # Connect buttons to their callbacks
        self.builder.get_object('button_cancel').connect('clicked', self.close_window)
        self.builder.get_object('button_add').connect('clicked', self.on_button_ok_clicked)

        # Variables to store the dialog results
        self.dist  = None
        self.force = None
        self.ok    = False   # Will be set True if user confirms the dialog
        
        # Run the dialog (blocks until user interacts)
        self.builder.get_object('dialog').run()
    

    def on_button_ok_clicked(self, widget):
        """ 
        Event handler for the "Add"/"Ok" button.  
        Collects distance and force values entered by the user,  
        sets the confirmation flag, and closes the dialog. 
        """
        self.dist  = self.builder.get_object('entry_dmin_coord1').get_text()
        self.force = self.builder.get_object('entry_FORCE_coord1').get_text()
        self.ok    = True
        self.builder.get_object('dialog').destroy()
        
    def close_window(self, button, data=None):
        """ 
        Event handler for cancel or dialog close.  
        Closes the dialog without applying changes. 
        """
        self.builder.get_object('dialog').destroy()
        self.Visible = False
