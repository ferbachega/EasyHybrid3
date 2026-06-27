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


class SolvateSystemWindow:
    """ Class doc """
    
    def __init__(self, main = None):
        """ Class initialiser """
        self.main = main
        self.vm_session      = main.vm_session
        self.Visible         = False        
        self.home            = main.home
        self.p_session       = main.p_session
    
    def open_window (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/setup/solvate_system_window.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('window')
            self.window.set_title('Solvate System Window')
            #self.window.set_keep_above(True)
            
            
            '''--------------------------------------------------------------------------------------------'''
            self.box = self.builder.get_object('box_system')
            self.combobox_systems = SystemComboBox(self.main )
            self.combobox_systems.connect("changed", self.on_combobox_systemsbox_changed)
            self.box.pack_start(self.combobox_systems, False, False, 0)
            '''--------------------------------------------------------------------------------------------'''

            # - - - - - - - coordinates combobox - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            self.box_coordinates = self.builder.get_object('box_coordinates')
            self.coordinates_combobox = CoordinatesComboBox() #self.builder.get_object('coordinates_combobox')
            self.box_coordinates.pack_start(self.coordinates_combobox, False, False, 0)
            '''--------------------------------------------------------------------------------------------'''
            
            
            '''--------------------------------------------------------------------------------------------'''
            self.box_cation  = self.builder.get_object('box_cation')
            self.cation_filechooser =  FolderChooserButton(self.main, 'file', self.home)
            self.box_cation.pack_start(self.cation_filechooser.btn, False, False, 0)
            
            self.box_anion  = self.builder.get_object('box_anion')
            self.anion_filechooser =  FolderChooserButton(self.main, 'file', self.home)
            self.box_anion.pack_start(self.anion_filechooser.btn, False, False, 0)

            #self.box_solvent_box     =  self.builder.get_object('box_solvent_box')
            #self.solvent_filechooser =  FolderChooserButton(self.main, 'file', self.home)
            #self.box_solvent_box.pack_start(self.solvent_filechooser.btn, False, False, 0)
            '''--------------------------------------------------------------------------------------------'''
            
            self.box_solvent_box2         = self.builder.get_object('box_solvent_box2')
            self.combobox_solvent_system  = SystemComboBox(self.main )
            self.combobox_solvent_system.connect("changed", self.on_combobox_solvent_system_changed)
            self.box_solvent_box2.pack_start(self.combobox_solvent_system, False, False, 0)
            
            
            
            self.window.show_all()                                               
            #self.solvent_filechooser.btn.hide()
            
            
            self.window.connect('destroy', self.close_window)                                               
            #self.combobox_systems.set_active(0)
            self.visible    =  True
    
    def close_window (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False

    def on_combobox_systemsbox_changed(self, widget):
        """ Function doc """
        system_id = self.combobox_systems.get_system_id()
       
        if system_id is not None:
            self.coordinates_combobox.set_model(self.main.vobject_liststore_dict[system_id])
            #self.refresh_selection_liststore (system_id)            
            size  =  len(list(self.main.vobject_liststore_dict[system_id]))
            self.coordinates_combobox.set_active(size-1)
            
    def on_combobox_solvent_system_changed (self, widget):
        """ Function doc """
        print('on_combobox_solvent_system_changed')
        solvent_box_id = self.combobox_solvent_system.get_system_id()
        solvent_box = self.p_session.psystem[solvent_box_id]
        if solvent_box.symmetry:
            x = solvent_box.symmetryParameters.a
            y = solvent_box.symmetryParameters.b
            z = solvent_box.symmetryParameters.c

            #alpha  = solvent_box.symmetryParameters.alpha                                                                                                              │
            #beta   = solvent_box.symmetryParameters.beta                                                                                                               │
            #gamma  = solvent_box.symmetryParameters.gamma
            
            self.builder.get_object('entry_a').set_text(str(x))
            self.builder.get_object('entry_b').set_text(str(y))
            self.builder.get_object('entry_c').set_text(str(z))
        else:
            msg = 'The selected system has no cell parameters!'
            self.main.simple_dialog.info(msg = msg )



    def on_checkbox_add_ions_toggled (self, widget):
        """ Function doc """
        if self.builder.get_object('checkbox_add_ions').get_active():
            self.builder.get_object('frame_add_ions').set_sensitive(True)
        else:
            self.builder.get_object('frame_add_ions').set_sensitive(False)
    
    def on_checkbox_add_solvent_box_toggled (self, widget):
        """ Function doc """
        
        #if self.builder.get_object('checkbox_add_ions').get_active():
        #    self.builder.get_object('frame_add_ions').set_sensitive(True)
        #else:
        #    self.builder.get_object('frame_add_ions').set_sensitive(False)
        
        print('on_checkbox_add_solvent_box_toggled')

    def on_button_run_clicked (self, widget):
        """ Function doc """
        system_id      = self.combobox_systems.get_system_id()
        solvent_box_id = self.combobox_solvent_system.get_system_id()        
        parameters = {}
        
        parameters['XBox'] = float(self.builder.get_object('entry_a').get_text())
        parameters['YBox'] = float(self.builder.get_object('entry_b').get_text())
        parameters['ZBox'] = float(self.builder.get_object('entry_c').get_text())
        
        if self.builder.get_object('checkbox_add_ions').get_active():
            #---------------------------- I o n s ----------------------------------------
            parameters['NPositive'] = int(self.builder.get_object('spinbtn_cations').get_value())
            if parameters['NPositive'] == 0:
                parameters['cation']    = None
            else:
                parameters['cation']    = self.cation_filechooser.get_folder ()
            
            
            parameters['NNegative'] = int(self.builder.get_object('spinbtn_anions').get_value())
            if parameters['NNegative'] == 0:
                parameters['anion']     = None
            else:
                parameters['anion']     = self.anion_filechooser.get_folder ()
            #------------------------------------------------------------------------------
        else:
            parameters['NPositive'] = 0
            parameters['NNegative'] = 0
            parameters['cation']    = None
            parameters['anion']     = None
        
        
        #------------------------------------------------------------------------------
        if self.builder.get_object('checkbox_reorient').get_active():
            parameters['reorient']  = True
        else:
            parameters['reorient']  = False
        #------------------------------------------------------------------------------

        
        #------------------------------------------------------------------------------
        #if self.builder.get_object('checkbox_solvente_box_from_file').get_active():
        #    parameters['solvent'] = None
        #else:
        parameters['solvent'] = self.p_session.psystem[solvent_box_id]
        #------------------------------------------------------------------------------
        
        #print(parameters)
        self.p_session.solvate_system (e_id = system_id, parameters = parameters)


class MakeSolventBoxWindow:
    """ Class doc """
    
    def __init__(self, main = None):
        """ Class initialiser """
        self.main = main
        self.vm_session      = main.vm_session
        self.Visible         = False        
        self.home            = main.home
        self.p_session       = main.p_session
        
    def open_window (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/setup/make_solvent_box_window.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('window')
            self.window.set_title('Make Solvent Box Window')
            self.window.set_keep_above(True)
            
            # - - - - - - - systems combobox - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            self.box = self.builder.get_object('box_system')
            self.combobox_systems = SystemComboBox(self.main )
            self.combobox_systems.connect("changed", self.on_combobox_systemsbox_changed)
            '''--------------------------------------------------------------------------------------------'''
            self.box.pack_start(self.combobox_systems, False, False, 0)



            # - - - - - - - coordinates combobox - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            self.box_coordinates = self.builder.get_object('box_coordinates')
            self.coordinates_combobox = CoordinatesComboBox() #self.builder.get_object('coordinates_combobox')
            self.box_coordinates.pack_start(self.coordinates_combobox, False, False, 0)
            '''--------------------------------------------------------------------------------------------'''
        
            
            self.btn_run = self.builder.get_object('button_run')
            self.btn_run.connect('clicked', self.run)
            system  = self.main.p_session.get_system()
            self.combobox_systems.set_active_iter(system.e_liststore_iter)
            
            self.window.show_all()                                               
            self.window.connect('destroy', self.close_window)                                               
            #self.combobox_systems.set_active(0)
            self.visible    =  True
            '''--------------------------------------------------------------------------------------------'''
        
        #print(idnum, text )
    def close_window (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False
        
    def on_combobox_systemsbox_changed(self, widget):
        """ Function doc """
        system_id = self.combobox_systems.get_system_id()
       
        if system_id is not None:
            self.coordinates_combobox.set_model(self.main.vobject_liststore_dict[system_id])
            #self.refresh_selection_liststore (system_id)            
            size  =  len(list(self.main.vobject_liststore_dict[system_id]))
            self.coordinates_combobox.set_active(size-1)
            
            #self.update_window ( selections = False, restraints = True)
        

    def run (self, widget):
        """ Function doc """
        
        parameters = {}
        
        parameters['_Density'] = int(self.builder.get_object('entry_density').get_text())
        parameters['_Steps']   = int(self.builder.get_object('entry_number_of_steps').get_text())
        parameters['_XBox']    = int(self.builder.get_object('entry_size_X').get_text())
        parameters['_YBox']    = int(self.builder.get_object('entry_size_Y').get_text())
        parameters['_ZBox']    = int(self.builder.get_object('entry_size_Z').get_text())
        parameters['_Refine']  = True
        
        system_id = self.combobox_systems.get_system_id()
        '''selecting the vismol object from the content that is in the combobox '''
        vobject_id = self.coordinates_combobox.get_vobject_id()
        vobject    = self.main.vm_session.vm_objects_dic[vobject_id]
        
        '''This function imports the coordinates of a vobject into the dynamo system in memory.''' 
        self.main.p_session.set_psystem_coordinates_from_vobject(vobject   = vobject, 
                                                                           system_id = system_id )
        
        parameters['molecule'] = self.main.p_session.psystem[system_id]
        self.p_session.make_solvent_box(parameters)
        self.close_window(None)
