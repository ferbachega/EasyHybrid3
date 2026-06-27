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

# --- cross-module imports added during refactor ---
from gui.windows.setup.windows_and_dialogs.qc_setup.setup_orca import SetupORCAWindow
from gui.windows.setup.windows_and_dialogs.qc_setup.setup_dftbplus import SetupDFTBplusWindow
from gui.windows.setup.windows_and_dialogs.qc_setup.setup_xtb import SetupXTBWindow
from gui.windows.setup.windows_and_dialogs.qc_setup.setup_mopac import SetupMOPACWindow


VISMOL_HOME      = os.environ.get('VISMOL_HOME')
HOME             = os.environ.get('HOME')
PDYNAMO3_SCRATCH = os.environ.get('PDYNAMO3_SCRATCH')


class EasyHybridSetupQCModelWindow:
    """ Class doc """
    
    def __init__(self, main = None):
        """ Class initialiser """
        self.main_session     = main
        self.home             = main.home
        self.Visible          = False        
        self.vismol_object    = None 
        self.window           = None
        self.methods_liststore = Gtk.ListStore(str, str, str)
        
        self.method_id    = 0    # 0 for am1, 1 for pm3  and...
        self.charge       = 0
        self.multiplicity = 1
        self.restricted   = True
        
        self.adjustment_charge = Gtk.Adjustment(value=0,
                                           lower=-100,
                                           upper=100,
                                           step_increment=1,
                                           page_increment=1,
                                           page_size=0)
        
        self.adjustment_multiplicity = Gtk.Adjustment(value=1,
                                           lower=1,
                                           upper=100,
                                           step_increment=1,
                                           page_increment=1,
                                           page_size=0)
        
        self.QC_engines = {
                          0:'pDynamo',
                          1:'ORCA'   ,
                          2:'DFTB+'  ,
                          3:'xTB'    ,  
                          #4:'MOPAC'  ,
                          #5:'sparrow',
                          }
        
        
        self.methods_id_dictionary = {
                                      0 : 'am1'             ,
                                      1 : 'am1dphot'        ,
                                      2 : 'pm3'             ,
                                      3 : 'pm6'             ,
                                      4 : 'rm1'             ,
                                      5 : 'mndo'            ,
                                      #6 : 'ab initio - ORCA',
                                      #7 : 'DFTB+'           ,
                                      #8 : 'xTB'             ,
                                      #9 : 'MOPAC'           ,
                                      }
        
        self.setup_orca_window = SetupORCAWindow(self.main_session, self)
        self.setup_dftb_window = SetupDFTBplusWindow(self.main_session, self)
        self.setup_xtb_window  = SetupXTBWindow(self.main_session, self)
        self.setup_mopac_window= SetupMOPACWindow(self.main_session, self)
        
        self.orca_options = ''
        self.orca_scratch = ''
        self.orca_random_scratch = False
        
    def refresh (self, option = 'all'):
        """ Function doc """
        self.update_number_of_qc_atoms()
       
    def open_window (self, vismol_object = None):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/setup/setup_qc_model_window.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('SetupQCModelWindow')
            self.window.set_keep_above(True)
            
            #self.window.set_default_size(400, 150)
            self.vismol_object = vismol_object
            
            if self.vismol_object is not None:
                system = self.main_session.p_session.psystem[self.vismol_object.e_id]
            else:
                system = self.main_session.p_session.psystem[self.main_session.p_session.active_id]
                
                
                
            self.window.set_title('Quantum Chemistry - '+ system.label)
            
            
            
            
            pixbuf = get_colorful_square_pixel_buffer (system)
            self.builder.get_object('system_sqr_image').set_from_pixbuf(pixbuf)
            
            
            
            
            '''-----------------------------------------QC ENGINES------------------------------------------'''
            self.qc_engine_store = Gtk.ListStore(str)
            engines = self.QC_engines.values()          
            for engine in engines:
                self.qc_engine_store.append([engine])
                #print (method_type)
            
            self.qc_engine_combobox = self.builder.get_object('QCEngine_combobox')
            
            self.qc_engine_combobox.connect("changed", self.on_qc_engine_combo_changed) #mudar
            self.qc_engine_combobox.set_model(self.qc_engine_store)
            renderer_text = Gtk.CellRendererText()
            self.qc_engine_combobox.pack_start(renderer_text, True)
            self.qc_engine_combobox.add_attribute(renderer_text, "text", 0)
            self.qc_engine_combobox.set_active(0)
            '''--------------------------------------------------------------------------------------------'''
            
            '''--------------------------------------------------------------------------------------------'''
            self.methods_type_store = Gtk.ListStore(str)
            methods_types = self.methods_id_dictionary.values()          
            for method_type in methods_types:
                self.methods_type_store.append([method_type])
                #print (method_type)
            
            self.methods_combo = self.builder.get_object('QCModel_methods_combobox')
            self.methods_combo.connect("changed", self.on_name_combo_changed)
            self.methods_combo.set_model(self.methods_type_store)
            renderer_text = Gtk.CellRendererText()
            self.methods_combo.pack_start(renderer_text, True)
            self.methods_combo.add_attribute(renderer_text, "text", 0)
            self.methods_combo.set_active(0)
            
            '''--------------------------------------------------------------------------------------------'''


            self.spinbutton_charge       = self.builder.get_object('spinbutton_charge'      )
            self.spinbutton_multiplicity = self.builder.get_object('spinbutton_multiplicity')
            self.spinbutton_charge      .set_adjustment(self.adjustment_charge)
            self.spinbutton_multiplicity.set_adjustment(self.adjustment_multiplicity)
            
            self.window.show_all()                                               
            self.builder.connect_signals(self)                                   
            self.methods_combo.set_active(self.method_id)
            

            

            
            
            self.button_ok = self.builder.get_object('button_ok') 
            self.button_cancel = self.builder.get_object('button_cancel') 
            
            self.button_ok.connect("clicked", self.on_button_ok)
            self.button_cancel.connect("clicked", self.close_window)
            self.window.connect("destroy", self.close_window)
            
            
            self.button_orca_setup = self.builder.get_object('button_setup_orca') 
            self.button_orca_setup.connect("clicked", self.on_button_setup_orca)
            self.button_orca_setup = self.builder.get_object('button_setup_dftb') 
            self.button_orca_setup.connect("clicked", self.on_button_setup_dftb)
            
            self.builder.get_object('button_setup_orca').hide()
            self.builder.get_object('button_setup_dftb').hide()
            
            
            ''' Updating the number of atoms '''
            self.update_number_of_qc_atoms ()
            
            self.Visible  =  True

        else:
            ''' Updating the number of atoms '''
            self.update_number_of_qc_atoms ()
            self.window.present()
    
    def update_number_of_qc_atoms (self):
        """ Function doc """
        self.entry_number_of_qc_atoms = self.builder.get_object('entry_number_of_qc_atoms')
        
        '''   Estiamting the QC charge  '''
        '''----------------------------------------------------------------------------------------------'''
        system  = self.main_session.p_session.psystem[self.main_session.p_session.active_id]
        estimated_charge = 0.0
        
        if len(system.e_qc_table) > 0:
            
            is_mmState   = getattr(system, 'mmState', False) 
            
            if is_mmState:
                for index in system.e_qc_table:
                    try:
                        estimated_charge += system.mmState.charges[index]
                    except:
                        print('System object has no attribute mmState - pure QC system')
                        estimated_charge += 0
            else:
                estimated_charge = 0
            
            estimated_charge = int(round(estimated_charge))
            self.spinbutton_charge.set_value (estimated_charge)
            
            number_of_qc_atoms = len(system.e_qc_table)
            self.entry_number_of_qc_atoms.set_text(str(number_of_qc_atoms))
        
        else:
            
            if system.mmModel:
                for charge in system.mmState.charges:
                    estimated_charge += charge
                
                estimated_charge = int(round(estimated_charge))
                self.spinbutton_charge.set_value (estimated_charge)
            else:
                print('mmModel not found, pure QC system.')
            number_of_qc_atoms = len(system.atoms)
            self.entry_number_of_qc_atoms.set_text(str(number_of_qc_atoms)+ ' (all)')
        '''----------------------------------------------------------------------------------------------'''

    def close_window (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False

        if self.setup_orca_window.Visible:
            self.setup_orca_window.close_window ( None,  None)
        
        if self.setup_dftb_window.Visible:
            self.setup_dftb_window.close_window ( None,  None)
        
        if self.setup_xtb_window.Visible:
            self.setup_dftb_window.close_window ( None,  None)

    def on_spian_button_change (self, widget):
        """ Function doc """
        self.charge       = self.spinbutton_charge.get_value_as_int()
        self.multiplicity = self.spinbutton_multiplicity.get_value_as_int()
        
    def on_qc_engine_combo_changed (self, combobox):
        self.qc_engine_id = self.builder.get_object('QCEngine_combobox').get_active()
        print(self.qc_engine_id)
        
        if self.qc_engine_id ==0:            
            self.builder.get_object('button_setup_orca').hide()
            self.builder.get_object('button_setup_dftb').hide()
            self.builder.get_object('expander_DIISSCF_converger').show()
            
            self.builder.get_object('label_method').show()
            self.builder.get_object('QCModel_methods_combobox').show()
            
        elif self.qc_engine_id == 1:
            self.builder.get_object('button_setup_orca').show()
            self.builder.get_object('button_setup_dftb').hide()
            self.builder.get_object('expander_DIISSCF_converger').hide()

            self.builder.get_object('label_method').hide()
            self.builder.get_object('QCModel_methods_combobox').hide()

            self.setup_orca_window.open_window()
        
        elif self.qc_engine_id == 2:
            self.builder.get_object('button_setup_dftb').set_label('DFTB+ Setup')
            self.builder.get_object('button_setup_dftb').show()
            
            self.builder.get_object('button_setup_orca').hide()
            self.builder.get_object('expander_DIISSCF_converger').hide()
            self.builder.get_object('label_method').hide()
            self.builder.get_object('QCModel_methods_combobox').hide()
            self.setup_dftb_window.open_window()
        
        elif self.qc_engine_id == 3:
            self.builder.get_object('button_setup_dftb').set_label('xTB Setup')
            self.builder.get_object('button_setup_dftb').show()
            
            self.builder.get_object('button_setup_orca').hide()
            self.builder.get_object('expander_DIISSCF_converger').hide()
            self.builder.get_object('label_method').hide()
            self.builder.get_object('QCModel_methods_combobox').hide()
            self.setup_xtb_window.open_window()

        elif self.qc_engine_id == 4:
            self.builder.get_object('button_setup_dftb').set_label('MOPAC Setup')
            self.builder.get_object('button_setup_dftb').show()
            
            self.builder.get_object('button_setup_orca').hide()
            self.builder.get_object('expander_DIISSCF_converger').hide()
            self.builder.get_object('label_method').hide()
            self.builder.get_object('QCModel_methods_combobox').hide()
            self.setup_mopac_window.open_window()
        else:
            pass

    def on_name_combo_changed (self, combobox):
        """ Function doc """
        pass
        #
        self.method_id = self.builder.get_object('QCModel_methods_combobox').get_active()
        print(self.method_id)
        if self.method_id in [0,1,2,3,4,5]:            
            self.builder.get_object('button_setup_orca').hide()
            self.builder.get_object('button_setup_dftb').hide()
            self.builder.get_object('expander_DIISSCF_converger').show()
            
        elif self.method_id == 6:
            self.builder.get_object('button_setup_orca').show()
            self.builder.get_object('button_setup_dftb').hide()
            self.builder.get_object('expander_DIISSCF_converger').hide()
            self.setup_orca_window.open_window()
        
        elif self.method_id == 7:
            self.builder.get_object('button_setup_dftb').set_label('DFTB+ Setup')
            self.builder.get_object('button_setup_dftb').show()
            self.builder.get_object('button_setup_orca').hide()
            self.builder.get_object('expander_DIISSCF_converger').hide()
            self.setup_dftb_window.open_window()
        
        elif self.method_id == 8:
            self.builder.get_object('button_setup_dftb').set_label('xTB Setup')
            self.builder.get_object('button_setup_dftb').show()
            self.builder.get_object('button_setup_orca').hide()
            self.builder.get_object('expander_DIISSCF_converger').hide()
            self.setup_xtb_window.open_window()
        
        elif self.method_id == 9:
            #self.builder.get_object('button_setup_dftb').set_label('xTB Setup')
            #self.builder.get_object('button_setup_dftb').show()
            #self.builder.get_object('button_setup_orca').hide()
            #self.builder.get_object('expander_DIISSCF_converger').hide()
            print('here')
            self.setup_mopac_window.open_window()
        else:
            pass

    def on_button_ok (self, button):
        """ Function doc """
        
        if self.builder.get_object('radio_button_restricted').get_active():
            #print("%s is active" % (self.builder.get_object('radio_button_restricted').get_label()))
            self.restricted = True
        else:
            #print("%s is not active" % (self.builder.get_object('radio_button_restricted').get_label()))
            self.restricted = False
        

        
        parameters = {
                    'qcengine'         : self.QC_engines[self.qc_engine_id]             ,
                    'qcmodel'          : 'mndo'           ,
                    'charge'           : self.spinbutton_charge.get_value_as_int()      ,
                    'multiplicity'     : self.spinbutton_multiplicity.get_value_as_int(),
                    
                    'method'           : self.methods_id_dictionary[self.method_id]   ,
                    
                    'isSpinRestricted' : self.restricted  ,
                     }

        if self.qc_engine_id == 1: # ORCA
            parameters['orca_options'  ]   = self.orca_options
            parameters['orca_scratch'  ]   = self.orca_scratch
            parameters['random_scratch'  ] = self.orca_random_scratch
        
        elif self.qc_engine_id == 2: # DFTB+
            parameters["deleteJobFiles"      ] = self.setup_dftb_window.delete_job_files # True/False 
            parameters["extendedInput"       ] = self.setup_dftb_window.text_extended_input #None
            #parameters['method'              ] = 'DFTB+'
            parameters["fermiTemperature"    ] = self.setup_dftb_window.fermiTemperature
            parameters["gaussianBlurWidth"   ] = self.setup_dftb_window.gaussianBlurWidth
            #parameters["hamiltonian"         ] = "DFTB"
            parameters["maximumSCCIterations"] = self.setup_dftb_window.maximumSCCIterations
            parameters["randomScratch"       ] = self.setup_dftb_window.random_scratch
            parameters["sccTolerance"        ] = self.setup_dftb_window.sccTolerance
            parameters["scratch"             ] = self.setup_dftb_window.scratch_folder
            parameters["skfPath"             ] = self.setup_dftb_window.skf_folder 
            parameters["useSCC"              ] = self.setup_dftb_window.use_scc
        
        elif self.qc_engine_id == 3: # xTB
            parameters['gfn'       ] = self.setup_xtb_window.parameters['gfn']  
            parameters['parallel'  ] = self.setup_xtb_window.parameters['parallel'] 
            parameters['acc'       ] = self.setup_xtb_window.parameters['acc']
            parameters['iterations'] = self.setup_xtb_window.parameters['iterations']
            parameters['fermi_temp'] = self.setup_xtb_window.parameters['fermi_temp'  ] 
            parameters['keywords'  ] = self.setup_xtb_window.parameters['add_keywords'] 
            parameters['scratch'   ] = self.setup_xtb_window.parameters['scratch'] 
        
        elif self.qc_engine_id == 4:
            #parameters['method'    ]  = 'MOPAC'
            parameters['hamiltonian'] = self.setup_mopac_window.parameters['method']  
            parameters['keywords'  ]  = self.setup_mopac_window.parameters['add_keywords'] 
            parameters['scratch'   ]  = self.setup_mopac_window.parameters['scratch'] 


        parameters['energyTolerance'  ] = float(self.builder.get_object('entry_energyTolerance').get_text())
        parameters['densityTolerance' ] = float(self.builder.get_object('entry_densityTolerance').get_text())
        parameters['maximumIterations'] = int(self.builder.get_object('entry_maximumIterations').get_text())

        print(parameters)
        
        self.main_session.p_session.define_a_new_QCModel(system        = None,  
                                                         parameters    = parameters, 
                                                         vismol_object = self.vismol_object)
        #self.main_session.update_gui_widgets ()
        self.window.destroy()
        self.Visible    =  False

    def on_button_setup_orca (self, button):
        """ Function doc """
        self.setup_orca_window.open_window()
    
    def on_button_setup_dftb (self, button):
        """ Function doc """
        if self.method_id == 7:
            self.setup_dftb_window.open_window()
        
        elif self.method_id == 8:
            print(self.method_id, self.setup_xtb_window)
            self.setup_xtb_window.open_window() 
        
        else:
            pass
        
    def update (self):
        """ Function doc """
        pass
