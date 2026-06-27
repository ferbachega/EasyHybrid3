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


class SetupORCAWindow:
    """ Class doc """
    
    def __init__ (self, main, setup_QC_model_window):
        """ Class initialiser """
        self.main_session     = main
        self.home             = main.home
        self.Visible          = False        
        self.vismol_object    = None 
        
        self.setup_QC_model_window = setup_QC_model_window
        self.restricted   = True
        
        self.adjustment_nCPU = Gtk.Adjustment(value=1,
                                              lower=1,
                                              upper=1000,
                                              step_increment=1,
                                              page_increment=1,
                                              page_size=0)
                                           
                                           
        self.orca_methods_dict = {}
        self.restrited_label = ''
        
        #---------------------------------------------------------------------
        # ORCA
        # this liststore contains the types of methods (HF, local and gradient corrections, ...)
        self.liststore_orca_type_of_method = Gtk.ListStore(str)
        for key, _type in orca_keys.type_of_method.items():
            self.liststore_orca_type_of_method.append([_type])

            # this liststore contains the methods (HF, B1LYP, B3LYP ...)
            self.orca_methods_dict[key] = Gtk.ListStore(str)
            
            for method in orca_keys.methods[key]:
                self.orca_methods_dict[key].append([method])
        #---------------------------------------------------------------------
        
        
        #---------------------------------------------------------------------
        #excited states
        self.excited_states_liststore_dict = {}
        self.excited_states_liststore_dict['semiemp'] = Gtk.ListStore(str)
        for excited_state in orca_keys.list_semiemp_excited.values():
            self.excited_states_liststore_dict['semiemp'].append([excited_state])
        
        self.excited_states_liststore_dict['HF'] = Gtk.ListStore(str)
        for excited_state in orca_keys.list_HF_excited.values():
            self.excited_states_liststore_dict['HF'].append([excited_state])
        
        
        self.excited_states_liststore_dict['DFT'] = Gtk.ListStore(str)
        for excited_state in orca_keys.list_DFT_excited.values():
            self.excited_states_liststore_dict['DFT'].append([excited_state])
        
        #---------------------------------------------------------------------
        
        
        
        #---------------------------------------------------------------------
        self.liststore_orca_scf_convergence = Gtk.ListStore(str)
        for scf_convergence in orca_keys.list_scf_convergence.values():
            self.liststore_orca_scf_convergence.append([scf_convergence])
        #---------------------------------------------------------------------
        
        
        #---------------------------------------------------------------------
        #basis set and basis set group
        self.orca_type_basis_dict = {}                       # contain  several liststores
        self.liststore_orca_type_basis = Gtk.ListStore(str)  # 
        
        for key , type_basis in  orca_keys.basis_set_group_dict.items():

            self.liststore_orca_type_basis.append([type_basis])
            self.orca_type_basis_dict[key] =  Gtk.ListStore(str)
            for basis in orca_keys.type_basis_dict[key]:
                self.orca_type_basis_dict[key].append([basis])   
        #---------------------------------------------------------------------
        '''
        This is a dictionary that stores combobox 
        information in case you need to regenerate the window.
        '''
        self.random_scratch = False
        self.comboboxes = {}


    
    def open_window (self, vismol_object = None):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/setup/setup_qc_model_window.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('window_setup_orca')
            self.window.set_keep_above(True)

            #-------------------------------------------------------------------------------
            self.orca_type_of_method_combo = self.builder.get_object('combobox_orca_type')
            self.orca_type_of_method_combo.connect("changed", self.on_combobox_orca_type_changed)
            self.orca_type_of_method_combo.set_model(self.liststore_orca_type_of_method)
            renderer_text = Gtk.CellRendererText()
            self.orca_type_of_method_combo.pack_start(renderer_text, True)
            self.orca_type_of_method_combo.add_attribute(renderer_text, "text", 0)
            
            #-------------------------------------------------------------------------------
            
            #-------------------------------------------------------------------------------
            self.orca_methods_combo = self.builder.get_object('combobox_orca_method')
            self.orca_methods_combo.connect("changed", self.on_combobox_orca_method_changed)
            self.orca_methods_combo.set_model( self.orca_methods_dict[0] )
            renderer_text = Gtk.CellRendererText()
            self.orca_methods_combo.pack_start(renderer_text, True)
            self.orca_methods_combo.add_attribute(renderer_text, "text", 0)
            #
            #-------------------------------------------------------------------------------
            
            
            #-------------------------------------------------------------------------------
            self.combobox_orca_excited_states = self.builder.get_object('combobox_orca_excited_states')
            #self.combobox_orca_excited_states.connect("changed", self.on_combobox_orca_method_changed)
            self.combobox_orca_excited_states.set_model( self.excited_states_liststore_dict['HF'] )
            renderer_text = Gtk.CellRendererText()
            self.combobox_orca_excited_states.pack_start(renderer_text, True)
            self.combobox_orca_excited_states.add_attribute(renderer_text, "text", 0)
            #
            #-------------------------------------------------------------------------------
            
            
            #-------------------------------------------------------------------------------
            self.orca_scf_conv_combo = self.builder.get_object('combobox_orca_scf_convergence')
            #self.orca_scf_conv_combo.connect("changed", self.on_combobox_orca_method_changed)
            self.orca_scf_conv_combo.set_model( self.liststore_orca_scf_convergence )
            renderer_text = Gtk.CellRendererText()
            self.orca_scf_conv_combo.pack_start(renderer_text, True)
            self.orca_scf_conv_combo.add_attribute(renderer_text, "text", 0)
            #
            #-------------------------------------------------------------------------------
            
            
            
            #-------------------------------------------------------------------------------
            self.orca_basis_group_combo = self.builder.get_object('combobox_orca_basis_set_group')
            self.orca_basis_group_combo.connect("changed", self.on_orca_basis_group_combo_changed)
            self.orca_basis_group_combo.set_model( self.liststore_orca_type_basis )
            renderer_text = Gtk.CellRendererText()
            self.orca_basis_group_combo.pack_start(renderer_text, True)
            self.orca_basis_group_combo.add_attribute(renderer_text, "text", 0)
            #-------------------------------------------------------------------------------
            
            
            #-------------------------------------------------------------------------------
            self.orca_basis_set_combo = self.builder.get_object('combobox_orca_basis_set')
            #self.orca_basis_set_combo.connect("changed", self.on_orca_basis_set_combo_changed)
            self.orca_basis_set_combo.set_model( self.orca_type_basis_dict[0] )
            renderer_text = Gtk.CellRendererText()
            self.orca_basis_set_combo.pack_start(renderer_text, True)
            self.orca_basis_set_combo.add_attribute(renderer_text, "text", 0)
            #-------------------------------------------------------------------------------            
            
            
            #-------------------------------------------------------------------------------            
            self.spinbutton_orca_ncpus = self.builder.get_object('spin_button_orca_ncpus')
            self.spinbutton_orca_ncpus.set_adjustment(self.adjustment_nCPU)
            #-------------------------------------------------------------------------------            
            
            
            #-------------------------------------------------------------------------------            
            try:
                scratch = os.environ.get('PDYNAMO3_SCRATCH')
                self.builder.get_object('entry_orca_scratch_folder').set_text(scratch)
            except:
                scratch  = os.environ.get('HOME')
                self.builder.get_object('entry_orca_scratch_folder').set_text(scratch)
                
            #-------------------------------------------------------------------------------            
            
            
            
            
            
            '''-----------------------------------------------------------------------------------'''
            
            self.combobox_list = [self.orca_type_of_method_combo,
                                  self.orca_methods_combo,
                                  self.combobox_orca_excited_states,
                                  self.orca_scf_conv_combo,
                                  self.orca_basis_group_combo,
                                  self.orca_basis_set_combo,
                                 ]
            
            self.window.show_all()
            self.orca_type_of_method_combo.set_active(0)
            self.orca_methods_combo.set_active(0)
            
            self.orca_scf_conv_combo.set_active(0)
            self.orca_basis_group_combo.set_active(0)
            
            
            self.spinbutton_orca_ncpus.connect('value-changed', self.refresh_orca_parameters)

            self.orca_type_of_method_combo.connect('changed', self.refresh_orca_parameters)
            self.orca_methods_combo.connect('changed', self.refresh_orca_parameters)
            self.combobox_orca_excited_states.connect('changed', self.refresh_orca_parameters)
            self.orca_scf_conv_combo.connect('changed', self.refresh_orca_parameters)
            self.orca_basis_group_combo.connect('changed', self.refresh_orca_parameters)
            self.orca_basis_set_combo.connect('changed', self.refresh_orca_parameters)


            self.checkbox_ranbom_scratch = self.builder.get_object('checkbox_ranbom_scratch') 

            self.button_ok         = self.builder.get_object('orca_button_ok ') 
            self.button_cancel     = self.builder.get_object('orca_button_cancel ') 
            
            self.button_ok.connect("clicked", self.on_button_ok)
            #self.button_ok.connect("clicked", self.on_button_ok2)
            self.button_cancel.connect("clicked", self.close_window)
            
            self.window.connect("destroy", self.close_window)
            self.refresh_orca_parameters (None)
            self.Visible  =  True
            
    def on_orca_basis_group_combo_changed (self, widget):
        """ Function doc """
        _id = widget.get_active()
        self.orca_basis_set_combo.set_model( self.orca_type_basis_dict[_id])
        self.orca_basis_set_combo.set_active( 0 )        
    
    def on_combobox_orca_type_changed (self, widget):
        """ Function doc """
        _id = widget.get_active()
        self.orca_methods_combo.set_model( self.orca_methods_dict[_id] )
        self.orca_methods_combo.set_active( 0 )
        
        self.restrited = self.setup_QC_model_window.builder.get_object('radio_button_restricted').get_active()
        
        if _id in [0]:
            self.combobox_orca_excited_states.set_model( self.excited_states_liststore_dict['HF'] )
            if self.restrited:
                self.restrited_label = 'RHF'
            else:
                self.restrited_label = 'UHF'
                
                
        elif _id in [1,2,3,4]:
            self.combobox_orca_excited_states.set_model( self.excited_states_liststore_dict['DFT'] )
            if self.restrited:
                self.restrited_label = 'RKS'
            else:
                self.restrited_label = 'UKS'
        
        elif _id == 7:
            self.combobox_orca_excited_states.set_model( self.excited_states_liststore_dict['semiemp'] )
            if self.restrited:
                self.restrited_label = 'RHF'
            else:
                self.restrited_label = 'UHF'
        
        else:
            self.combobox_orca_excited_states.set_model( None )
            if self.restrited:
                self.restrited_label = 'RHF'
            else:
                self.restrited_label = 'UHF'
            
            
        self.combobox_orca_excited_states.set_active(0)
        
    def on_combobox_orca_method_changed (self, widget):
        """ Function doc """
        
    def build_ORCA_widgets(self):
        """ Function doc """
        self.builder.get_object('combobox_orca_type')
                         
    def close_window (self, button, data  = None):
        """ Function doc """
        if self.window:
            self.window.destroy()
        self.Visible    =  False
    
    def on_button_ok (self, button):
        """ Function doc """
        #print('on_button_ok')
        textbuffer = self.builder.get_object('textview_orca').get_buffer ()
        text = textbuffer.get_text (textbuffer.get_start_iter(), textbuffer.get_end_iter(), True)
        
        self.setup_QC_model_window.orca_options = text
        self.setup_QC_model_window.orca_scratch = self.builder.get_object('entry_orca_scratch_folder').get_text()
        
        self.setup_QC_model_window.orca_random_scratch = self.builder.get_object('checkbox_orca_random_scratch').get_active()
        
        print (text)
        print ('checkbox_random_scratch', self.setup_QC_model_window.orca_random_scratch, self.builder.get_object('checkbox_orca_random_scratch').get_active())
        self.close_window (None)
    
    def on_button_ok2 (self, button):
        """ Function doc """
        #print('on_button_ok2')
        pass

    def refresh_orca_parameters (self, widget):
        """ Function doc """
        #self.orca_type_of_method_combo,
        
        _iter = self.orca_methods_combo.get_active_iter()
        model = self.orca_methods_combo.get_model()
        method = model[_iter][0].split()
        method = method[0]

        
        _iter =self.orca_basis_set_combo.get_active_iter()
        model =self.orca_basis_set_combo.get_model()
        basis = model[_iter][0].split()
        basis = basis[0]
        
        
        _iter =self.orca_scf_conv_combo.get_active_iter()
        model =self.orca_scf_conv_combo.get_model()
        scf_conv = model[_iter][0].split()
        scf_conv = scf_conv[0]        
        if scf_conv == 'Default':
            scf_conv = ''
        
        
        _iter =self.combobox_orca_excited_states.get_active_iter()
        model =self.combobox_orca_excited_states.get_model()
        e_states = model[_iter][0].split()
        e_states = e_states[0]        
        
        if e_states == 'TD-DFT':
            e_states = '''
%tddft
     nroots 8 # the number of excited states to be calculated.
     maxdim 30 # the maximum dimension of the expansion space in the Davidson procedure.
 end #tddft'''
        else:
            e_states = ''
        
        
        #-----------------------------------------------------------------------
        ncpu = self.spinbutton_orca_ncpus.get_value_as_int()
        if ncpu > 1:
            ncpu = '\n%pal nprocs {} \n end'.format(ncpu)
        else:
            ncpu = ''
        
        #-----------------------------------------------------------------------
        
        
        
        
        textbuffer = self.builder.get_object('textview_orca').get_buffer ()
        text = textbuffer.get_text (textbuffer.get_start_iter(), textbuffer.get_end_iter(), True)

        new_text = '{} {} {} {} {} {}'.format(self.restrited_label,  method, basis, scf_conv, e_states, ncpu )
        textbuffer.set_text(new_text)
