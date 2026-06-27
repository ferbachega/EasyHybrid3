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


class SetupSPARROWWindow:
    """ Class initialiser """
    
    def __init__ (self, main, setup_QC_model_window):
    
        self.main_session     = main
        self.vm_session       = main.vm_session
        self.home             = main.home
        self.Visible          = False        
        self.vismol_object    = None 
        self.window           = None
        
        self.hamiltonians = { 
                        0:'AM1',
                        1:'PM3',       # Use the MNDO-PM3 Hamiltonian
                        2:'PM6',       # Use the PM6 Hamiltonian
                       # 3:'PM6-D3',    # Use the PM6 Hamiltonian with Grimme's corrections for dispersion
                       # 4:'PM6-DH+',   # Use the PM6 Hamiltonian with corrections for dispersion and hydrogen-bonding
                       # 5:'PM6-DH2',   # Use the PM6 Hamiltonian with corrections for dispersion and hydrogen-bonding
                       # 6:'PM6-DH2X',  # Use PM6 with corrections for dispersion and hydrogen and halogen bonding
                       # 7:'PM6-D3H4',  # Use PM6 with Rezac and Hobza's D3H4 correction
                       # 8:'PM6-D3H4X', # Use PM6 with Brahmkshatriya, et al.'s D3H4X correction
                       # 9:'PM6-ORG',   # Use the PM6-ORG Hamiltonian
                       #10:'PM7',       # Use the PM7 Hamiltonian
                       #11:'PM7-TS',    # Use the PM7-TS Hamiltonian (only for barrier heights)
                       #12:'MNDO',      # Use the MNDO hamiltonian
                       #13:'MNDOD',     # Use the MNDO-d hamiltonian
                       }
        
        self.parameters = {
                          'method'          : 'PM6',
                          #'gfn'          : 1  ,
                          #'parallel'     : 1  ,
                          #'acc'          : 1.0,
                          #'iterations'   : 300,
                          #'fermi_temp'   : 300.0,
                          'add_keywords' : '',
                          'scratch'      : None #os.path.join(scratch,'XTBScratch')
                          }
        
        self.setup_QC_model_window = setup_QC_model_window

    def open_window (self, vismol_object = None):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/setup/setup_qc_model_window.glade'))
            
            self.window = self.builder.get_object('window_setup_sparrow')
            self.window.set_keep_above(True)
            self.window.set_default_size(450, 200)
            
             
            self.button_ok         = self.builder.get_object('btn_sparrow_ok') 
            self.button_cancel     = self.builder.get_object('btn_sparrow_cancel') 
            self.button_ok.connect("clicked", self.on_button_ok)
            self.button_cancel.connect("clicked", self.close_window)
            
            self.box_sparrow = self.builder.get_object('box_sparrow_type')
            
            #-----------------------------------------------------------
            self.sparrow_cbox = Gtk.ComboBoxText()
            self.sparrow_cbox.set_entry_text_column(0)
            #self.gfn_cbox.connect("changed", self.on_combobox_changed)

            # Add items to the ComboBox
            for key in self.hamiltonians.keys():
                self.sparrow_cbox.append_text(self.hamiltonians[key])
            self.box_sparrow.pack_start(self.sparrow_cbox, False, False, 0)

            
            ##.Interface Widgets

            self.entry_keywords          = self.builder.get_object('entry_sparrow_keywords')
            self.entry_scratch           = self.builder.get_object('entry_sparrow_scratch')
            self.spinbtn_parallel        = self.builder.get_object('spinbtn_parallel1')
            self.check_box_denout_oldens = self.builder.get_object('check_box_denout_oldens')
            
            self.check_box_NOGPU      = self.builder.get_object('check_box_NOGPU')
            self.check_box_NOGPU.connect("toggled", self.check_box_nogpu_toggled)
            
            #.Interface Show All
            self.window.connect("destroy", self.close_window)
            self.window.show_all()
            
            
            if self.parameters['scratch'] is not None:
                pass
            else:
                if os.path.isdir(self.vm_session.vm_config.gl_parameters["tmp_path"]):
                    self.parameters['scratch'] = self.vm_session.vm_config.gl_parameters["tmp_path"]
                else:
                    self.parameters['scratch'] = PDYNAMO3_SCRATCH

            self.entry_scratch       .set_text(str(self.parameters['scratch']))
            
            self.sparrow_cbox.connect("changed", self.cbox_changed)
            self.sparrow_cbox.set_active(0)
            
            self.spinbtn_parallel.connect("value-changed", self.spinbtn_change)
            
            #self.refresh_orca_parameters (None)
            self.Visible  =  True
    
    def update_keywords (self):
        """ Function doc """
        self.parameters['method'] = self.sparrow_cbox.get_active()

        method = self.hamiltonians[self.parameters['method']]
        threads = int(self.spinbtn_parallel.get_value ())
        
        if self.check_box_NOGPU.get_active():
            nogpu = 'NOGPU'
        else:
            nogpu = ''
        text = '{} 1SCF GRADIENTS AUX THREADS={} {}'.format(method, threads, nogpu)
        self.entry_keywords.set_text(text)

    def spinbtn_change (self, widget, value = None):
        """ Function doc """
        self.update_keywords()

    def check_box_nogpu_toggled (self, widget):
        """ Function doc """
        self.update_keywords()
    
    def cbox_changed (self, widget):
        """ Function doc """
        self.update_keywords()

        
    def on_button_ok (self, widget):
        """ Function doc """
        self.parameters['denout_oldens'] = self.check_box_denout_oldens.get_active()
        self.parameters['method'] = self.entry_keywords.get_text()
        self.parameters['keywords'] = self.entry_keywords.get_text()
        self.parameters['scratch']  = self.entry_scratch .get_text()
        
        print(self.parameters)
        self.close_window (widget,None)

    def close_window (self, button, data  = None):
        """ Function doc """
        if self.window:
            self.window.destroy()
        self.Visible    =  False
