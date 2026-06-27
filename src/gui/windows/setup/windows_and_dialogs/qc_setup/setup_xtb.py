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


class SetupXTBWindow:
    """ Class doc """
    
    def __init__ (self, main, setup_QC_model_window):
        """ Class initialiser """
        self.main_session     = main
        self.vm_session       = main.vm_session
        self.home             = main.home
        self.Visible          = False        
        self.vismol_object    = None 
        self.window           = None
        
        

        self.parameters = {
                          'gfn'          : 1  ,
                          'parallel'     : 1  ,
                          'acc'          : 1.0,
                          'iterations'   : 300,
                          'fermi_temp'   : 300.0,
                          'add_keywords' : '',
                          'scratch'      : None#os.path.join(scratch,'XTBScratch')
                          }
        
        self.setup_QC_model_window = setup_QC_model_window

    def open_window (self, vismol_object = None):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/setup/setup_qc_model_window.glade'))
            #self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('window_setup_xtb')
            self.window.set_keep_above(True)
            self.window.set_default_size(450, 200)
            
            self.button_ok         = self.builder.get_object('btn_xtb_ok') 
            self.button_cancel     = self.builder.get_object('btn_xtb_cancel') 
            self.button_ok.connect("clicked", self.on_button_ok)
            self.button_cancel.connect("clicked", self.close_window)
            
            self.box_gfn = self.builder.get_object('box_gfn_type')
            
            #-----------------------------------------------------------
            self.gfn_cbox = Gtk.ComboBoxText()
            self.gfn_cbox.set_entry_text_column(0)
            #self.gfn_cbox.connect("changed", self.on_combobox_changed)

            # Add items to the ComboBox
            self.gfn_cbox.append_text("GFN1-xTB")
            self.gfn_cbox.append_text("GFN2-xTB")
            self.gfn_cbox.append_text("g-xTB")
            #self.gfn_cbox.append_text("Option 3")
            self.box_gfn.pack_start(self.gfn_cbox, False, False, 0)
            
            #.Interface Widgets
            self.spinbtn_parallel     = self.builder.get_object('spinbtn_parallel')
            self.entry_xtb_acc        = self.builder.get_object('entry_xtb_acc')
            self.entry_xtb_iterations = self.builder.get_object('entry_xtb_iterations')
            
            self.entry_xtb_fermi_temp = self.builder.get_object('entry_xtb_fermi_temp')
            self.entry_keywords       = self.builder.get_object('entry_keywords')
            self.entry_scratch        = self.builder.get_object('entry_scratch')
            
            
            #.Interface Show All
            self.window.connect("destroy", self.close_window)
            self.window.show_all()
            
            
            #.Assigning the previously adjusted parameters to the respective widgets.            
            self.gfn_cbox.set_active(self.parameters['gfn']-1)
            self.spinbtn_parallel.set_value(self.parameters['parallel'])
            
            if self.parameters['scratch'] is not None:
                pass
            else:
                if os.path.isdir(self.vm_session.vm_config.gl_parameters["tmp_path"]):
                    self.parameters['scratch'] = self.vm_session.vm_config.gl_parameters["tmp_path"]
                else:
                    self.parameters['scratch'] = PDYNAMO3_SCRATCH
                
            
            
            self.entry_xtb_acc       .set_text(str(self.parameters['acc'         ]))    
            self.entry_xtb_iterations.set_text(str(self.parameters['iterations'  ]))
            self.entry_xtb_fermi_temp.set_text(str(self.parameters['fermi_temp'  ]))
            self.entry_keywords      .set_text(str(self.parameters['add_keywords']))
            self.entry_scratch       .set_text(str(self.parameters['scratch']))
            
            
            
            
            
            
            
            #self.refresh_orca_parameters (None)
            self.Visible  =  True
            
    def on_button_ok (self, widget):
        """ Function doc """
        self.parameters['gfn']      = (self.gfn_cbox        .get_active() + 1)
        self.parameters['parallel'] = int(self.spinbtn_parallel.get_value ())
        
        self.parameters['acc'         ] = float(self.entry_xtb_acc       .get_text())
        self.parameters['iterations'  ] = int(self.entry_xtb_iterations.get_text())
        self.parameters['fermi_temp'  ] = float(self.entry_xtb_fermi_temp.get_text())
        self.parameters['add_keywords'] = self.entry_keywords      .get_text()
        self.parameters['scratch']      = self.entry_scratch       .get_text()
        
        print(self.parameters)
        print('on_button_ok')
        self.close_window (widget,None)
        
    def close_window (self, button, data  = None):
        """ Function doc """
        if self.window:
            self.window.destroy()
        self.Visible    =  False
