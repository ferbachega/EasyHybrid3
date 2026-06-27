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


class SetupDFTBplusWindow:
    """ Class doc """
    
    def __init__ (self, main, setup_QC_model_window):
        """ Class initialiser """
        self.main_session     = main
        self.home             = main.home
        self.Visible          = False        
        self.vismol_object    = None 
        self.window           = None
        self.setup_QC_model_window = setup_QC_model_window
        
        try:
            self.skf_folder = os.path.join(os.environ.get('PDYNAMO3_HOME'),'examples/dftbPlus/data/skf')
        except:
            self.skf_folder = os.environ.get('PDYNAMO3_HOME')
        
        try:
            self.scratch_folder = os.environ.get('PDYNAMO3_SCRATCH')
        except:
            self.skf_folder = os.environ.get('PDYNAMO3_HOME')
        
        
        self.Hubbard_Derivs_parameters = {
                                        'Br': -0.0573 ,
                                        'C' : -0.1492 ,
                                        'Ca': -0.0340 ,
                                        'Cl': -0.0697 ,
                                        'F' : -0.1623 ,
                                        'H' : -0.1857 ,
                                        'I' : -0.0433 ,
                                        'K' : -0.0339 ,
                                        'Mg': -0.02   ,
                                        'N' : -0.1535 ,
                                        'Na': -0.0454 ,
                                        'O' : -0.1575 ,
                                        'P' : -0.14   ,
                                        'S' : -0.11   ,
                                        'Zn': -0.03   ,
                                        }
        
        #self.parameters =

        #self.skf_folder           = None
        #self.scratch_folder       = None
        self.use_scc              = True
        self.delete_job_files     = True
        self.random_scratch       = False
                                  
        
        
        
        #self.ThirdOrderFull       = False
        #self.zeta                 = 4.00
        #self.HubbardDerivs        = None

        self.fermiTemperature     = 300
        self.gaussianBlurWidth    = 0.0
        self.maximumSCCIterations = 300
        self.sccTolerance         = 1.0e-8
        self.text_extended_input  = None




    def open_window (self, vismol_object = None):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/setup/setup_qc_model_window.glade'))
            #self.builder.connect_signals(self)
            
            
            #.Interface Widgets
            self.window = self.builder.get_object('window_setup_dftb')
            self.window.set_keep_above(True)
            self.window.set_default_size(450, 200)
            
            self.button_ok         = self.builder.get_object('dftb_button_ok') 
            self.button_cancel     = self.builder.get_object('dftb_button_cancel') 
            self.button_ok.connect("clicked", self.on_button_ok)
            self.button_cancel.connect("clicked", self.close_window)
            
            
            self.skf_folder_chooser = self.builder.get_object('folder_chooser_skf_files')            
            self.entry_dftb_scratch_folder = self.builder.get_object('entry_dftb_scratch_folder')
            
            self.checkbox_use_scc          = self.builder.get_object('checkbox_use_scc')
            self.checkbox_delete_job_files = self.builder.get_object('checkbox_delete_job_files')
            self.checkbox_random_scratch   = self.builder.get_object('checkbox_random_scratch')
            
            self.chk_extended_input   = self.builder.get_object('chk_extended_input')
            self.chk_extended_input.connect("toggled", self.on_chk_extended_input)
            
            
            
            
            self.entry_fermiTemperature      = self.builder.get_object('entry_fermiTemperature')       
            self.entry_gaussianBlurWidth     = self.builder.get_object('entry_gaussianBlurWidth')   
            self.entry_maximumSCCIterations  = self.builder.get_object('entry_maximumSCCIterations')
            self.entry_sccTolerance          = self.builder.get_object('entry_sccTolerance')        
            
            
            #. Interface Sholl All
            self.window.connect("destroy", self.close_window)
            self.window.show_all()
            self.on_chk_extended_input(None)
            
            #.Assigning the previously adjusted parameters to the respective widgets.            

            self.skf_folder_chooser.set_filename(self.skf_folder)
            self.entry_dftb_scratch_folder.set_text(self.scratch_folder)

            self.checkbox_use_scc         .set_active(self.use_scc         )
            self.checkbox_delete_job_files.set_active(self.delete_job_files)
            self.checkbox_random_scratch  .set_active(self.random_scratch  )


            #self.checkbox_ThirdOrderFull.set_active(self.ThirdOrderFull)
            #self.builder.get_object('entry_zeta').set_text(str(self.zeta))
            
            self.entry_fermiTemperature    .set_text(str(self.fermiTemperature    ))     
            self.entry_gaussianBlurWidth   .set_text(str(self.gaussianBlurWidth   ))    
            self.entry_maximumSCCIterations.set_text(str(self.maximumSCCIterations)) 
            self.entry_sccTolerance        .set_text(str(self.sccTolerance        ))         
  
    
            #self.refresh_orca_parameters (None)
            self.Visible  =  True
    
    

    def close_window (self, button, data  = None):
        """ Function doc """
        if self.window:
            self.window.destroy()
        self.Visible    =  False

    
    def on_chk_extended_input (self, widget):
        """ Function doc """
        
        if self.chk_extended_input.get_active():
            self.builder.get_object('scroll_text_viewer').show( )
        else:
            self.builder.get_object('scroll_text_viewer').hide( )
        
        '''
        self.builder.get_object('entry_zeta').set_text(str(4.00))
        
        if self.checkbox_ThirdOrderFull.get_active():
            self.builder.get_object('label_zeta').set_sensitive(True)
            self.builder.get_object('entry_zeta').set_sensitive(True)
            self.builder.get_object('btn_HubbardDerivs').set_sensitive(True)
        else:
            self.builder.get_object('label_zeta').set_sensitive(False)
            self.builder.get_object('entry_zeta').set_sensitive(False)
            self.builder.get_object('btn_HubbardDerivs').set_sensitive(False)
        '''
        
        #self.skf_folder       = self.skf_folder_chooser.get_filename()
        #self.scratch_folder   = self.entry_dftb_scratch_folder.get_text()
        #self.use_scc          = self.checkbox_use_scc         .get_active()
        #self.delete_job_files = self.checkbox_delete_job_files.get_active()
        #self.random_scratch   = self.checkbox_random_scratch  .get_active()
        #self.close_window (None, None)
    
    def on_button_ok (self, button):
        """ Function doc """
        self.skf_folder       = self.skf_folder_chooser.get_filename()
        self.scratch_folder   = self.entry_dftb_scratch_folder.get_text()
        self.use_scc          = self.checkbox_use_scc         .get_active()
        self.delete_job_files = self.checkbox_delete_job_files.get_active()
        self.random_scratch   = self.checkbox_random_scratch  .get_active()
        
        #if self.checkbox_ThirdOrderFull.get_active():
        #    self.ThirdOrderFull = True
        #    self.zeta = float(self.builder.get_object('entry_zeta').get_text())
        #    self.HubbardDerivs = self.Hubbard_Derivs_parameters
        #else:
        #    self.ThirdOrderFull = False
        #    self.zeta           = 0.0
        #    self.HubbardDerivs  = None

        self.fermiTemperature     = float(self.builder.get_object('entry_fermiTemperature').get_text()     )
        self.gaussianBlurWidth    = float(self.builder.get_object('entry_gaussianBlurWidth').get_text()    )
        self.maximumSCCIterations =   int(self.builder.get_object('entry_maximumSCCIterations').get_text() )
        self.sccTolerance         = float(self.builder.get_object('entry_sccTolerance').get_text()         )

        
        if self.chk_extended_input.get_active():
            textbuffer = self.builder.get_object('text_viewer').get_buffer ()
            self.text_extended_input = textbuffer.get_text (textbuffer.get_start_iter(), textbuffer.get_end_iter(), True)        
        else:
            self.text_extended_input = None
        
        self.close_window (None, None)
    
        
    
    
    def on_skf_folder_chooser_changed (self, widget):
        """ Function doc """
        pass
