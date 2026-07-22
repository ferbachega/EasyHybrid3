#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Selection utilities for pDynamo systems
#
#  Copyright 2022-2026 Fernando Bachega
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
    """
    "Solvate System" window.

    Adds counter-ions and/or superimposes a pre-built, periodic solvent box
    (e.g. water) around a solute system.
    """

    def __init__(self, main = None):
        """ Class initialiser """
        self.main = main
        self.vm_session      = main.vm_session
        self.Visible         = False        
        self.home            = main.home
        self.p_session       = main.p_session
    
    #-------------------------------------------------------------------------
    #  W I N D O W   L I F E C Y C L E
    #-------------------------------------------------------------------------
    def open_window (self):
        """ Function doc """
        if self.Visible == False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/setup/solvate_system_window.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('window')
            self.window.set_title('Solvate System')

            # -------------------- widget shortcuts --------------------
            self.label_info          = self.builder.get_object('label_info')
            self.frame_solvent_box   = self.builder.get_object('frame_solvent_box')
            self.frame_add_ions      = self.builder.get_object('frame_add_ions')
            self.checkbox_add_solvent_box = self.builder.get_object('checkbox_add_solvent_box')
            self.checkbox_add_ions        = self.builder.get_object('checkbox_add_ions')

            '''--------------------------------------------------------------------------------------------'''
            self.box = self.builder.get_object('box_system')
            self.combobox_systems = SystemComboBox(self.main )
            self.combobox_systems.connect("changed", self.on_combobox_systemsbox_changed)
            self.box.pack_start(self.combobox_systems, False, False, 0)
            '''--------------------------------------------------------------------------------------------'''

            # - - - - - - - coordinates (object/frame) combobox - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            self.box_coordinates = self.builder.get_object('box_coordinates')
            self.coordinates_combobox = CoordinatesComboBox()
            self.box_coordinates.pack_start(self.coordinates_combobox, False, False, 0)
            '''--------------------------------------------------------------------------------------------'''
            
            
            '''--------------------------------------------------------------------------------------------'''
            self.box_cation  = self.builder.get_object('box_cation')
            self.cation_filechooser =  FolderChooserButton(self.main, 'file', self.home)
            self.box_cation.pack_start(self.cation_filechooser.btn, False, False, 0)
            
            self.box_anion  = self.builder.get_object('box_anion')
            self.anion_filechooser =  FolderChooserButton(self.main, 'file', self.home)
            self.box_anion.pack_start(self.anion_filechooser.btn, False, False, 0)

            # . Start the ion file choosers empty instead of pointing at the
            #   EasyHybrid install directory (which is not a valid ion file
            #   and would otherwise crash a "Run" click if the user forgot
            #   to pick a file).
            self.cation_filechooser.folder = None
            self.cation_filechooser.label.set_text('No file selected')
            self.anion_filechooser.folder  = None
            self.anion_filechooser.label.set_text('No file selected')
            '''--------------------------------------------------------------------------------------------'''
            
            self.box_solvent_box2         = self.builder.get_object('box_solvent_box2')
            self.combobox_solvent_system  = SystemComboBox(self.main )
            self.combobox_solvent_system.connect("changed", self.on_combobox_solvent_system_changed)
            self.box_solvent_box2.pack_start(self.combobox_solvent_system, False, False, 0)
            
            
            self.window.show_all()

            self.window.connect('destroy', self.close_window)

            # . Pre-select the currently active system as the solute, and,
            #   if there is one, a periodic system as the solvent box.
            if self.p_session.psystem:
                self.combobox_systems.set_active_system(e_id = self.p_session.active_id)

            self._select_default_solvent_system()

            self.frame_add_ions.set_sensitive(self.checkbox_add_ions.get_active())
            self.frame_solvent_box.set_sensitive(self.checkbox_add_solvent_box.get_active())

            self.Visible = True
    
    def close_window (self, button = None, data = None):
        """ Function doc """
        if not self.Visible:
            return
        self.window.destroy()
        self.Visible = False

    #-------------------------------------------------------------------------
    #  H E L P E R S
    #-------------------------------------------------------------------------
    def _select_default_solvent_system (self):
        """ Pre-selects, in the solvent-box combobox, the first loaded
        system that actually looks like a periodic solvent box (i.e. has a
        symmetry/cell defined) and is not the solute itself. Leaves the
        combobox as-is (nothing selected) if no such system exists. """
        solute_id = self.combobox_systems.get_system_id()

        for index, row in enumerate(self.combobox_solvent_system.system_liststore):
            _, system_id, _ = row
            system = self.p_session.psystem.get(system_id, None)
            if system is not None and system.symmetry and system_id != solute_id:
                self.combobox_solvent_system.set_active(index)
                return

    #-------------------------------------------------------------------------
    #  S I G N A L S :  S Y S T E M   /   O B J E C T   S E L E C T I O N
    #-------------------------------------------------------------------------
    def on_combobox_systemsbox_changed(self, widget):
        """ Function doc """
        system_id = self.combobox_systems.get_system_id()
       
        if system_id is not None:
            self.coordinates_combobox.set_model(self.main.vobject_liststore_dict[system_id])
            #self.refresh_selection_liststore (system_id)            
            size  =  len(list(self.main.vobject_liststore_dict[system_id]))
            self.coordinates_combobox.set_active(size-1)
            
    def on_combobox_solvent_system_changed (self, widget):
        """ Auto-fills the ion-placement box (X, Y, Z) from the selected
        solvent-box system's own cell, since that is normally the box the
        user wants the counter-ions distributed in as well. """
        solvent_box_id = self.combobox_solvent_system.get_system_id()
        if solvent_box_id is None:
            return

        solvent_box = self.p_session.psystem[solvent_box_id]
        if solvent_box.symmetry:
            x = solvent_box.symmetryParameters.a
            y = solvent_box.symmetryParameters.b
            z = solvent_box.symmetryParameters.c

            self.builder.get_object('entry_a').set_text('{:.3f}'.format(x))
            self.builder.get_object('entry_b').set_text('{:.3f}'.format(y))
            self.builder.get_object('entry_c').set_text('{:.3f}'.format(z))

            self.label_info.set_text('Ion-placement box synced to the "{}" solvent box.'.format(solvent_box.label))
        else:
            msg = 'The selected system has no cell parameters, so it cannot be used as a solvent box!'
            self.main.simple_dialog.info(msg = msg )

    #-------------------------------------------------------------------------
    #  S I G N A L S :  C H E C K B O X E S   /   B U T T O N S
    #-------------------------------------------------------------------------
    def on_checkbox_add_ions_toggled (self, widget):
        """ Function doc """
        self.frame_add_ions.set_sensitive(widget.get_active())
    
    def on_checkbox_add_solvent_box_toggled (self, widget):
        """ Function doc """
        self.frame_solvent_box.set_sensitive(widget.get_active())
        if widget.get_active():
            self.label_info.set_text('A solvent box will be added around the (ionized) solute.')
        else:
            self.label_info.set_text('No solvent box will be added -- only counter-ions (if any).')

    def on_button_make_solvent_box_clicked (self, widget):
        """ Opens the "Make Solvent Box" tool, so the user can build a new
        pre-equilibrated solvent box without leaving this workflow. """
        self.main.make_solvent_box_window.open_window()

    def on_button_cancel_clicked (self, widget):
        """ Function doc """
        self.close_window()

    #-------------------------------------------------------------------------
    #  R U N
    #-------------------------------------------------------------------------
    def on_button_run_clicked (self, widget):
        """ Function doc """
        system_id = self.combobox_systems.get_system_id()

        if system_id is None:
            self.main.simple_dialog.info(msg = 'Please select a solute system.')
            return

        # . Make sure the coordinates in memory for the solute match the
        #   object/frame currently selected in the "Object / frame" combobox
        #   (this used to be a purely decorative combobox that had no effect
        #   on the actual solvation).
        vobject_id = self.coordinates_combobox.get_vobject_id()
        if vobject_id is not None:
            vobject = self.main.vm_session.vm_objects_dic.get(vobject_id, None)
            if vobject is not None:
                self.p_session.set_psystem_coordinates_from_vobject(vobject = vobject, system_id = system_id)

        parameters = {}

        parameters['add_solvent'] = self.checkbox_add_solvent_box.get_active()

        try:
            parameters['XBox'] = float(self.builder.get_object('entry_a').get_text())
            parameters['YBox'] = float(self.builder.get_object('entry_b').get_text())
            parameters['ZBox'] = float(self.builder.get_object('entry_c').get_text())
        except ValueError:
            self.main.simple_dialog.info(msg = 'Please provide valid numeric values for the ion-placement box (X, Y, Z).')
            return

        if self.checkbox_add_ions.get_active():
            #---------------------------- I o n s ----------------------------------------
            parameters['NPositive'] = int(self.builder.get_object('spinbtn_cations').get_value())
            if parameters['NPositive'] == 0:
                parameters['cation']    = None
            else:
                parameters['cation']    = self.cation_filechooser.get_folder ()
                if not parameters['cation']:
                    self.main.simple_dialog.info(msg = 'Please choose a cation file, or set the cation quantity to 0.')
                    return

            parameters['NNegative'] = int(self.builder.get_object('spinbtn_anions').get_value())
            if parameters['NNegative'] == 0:
                parameters['anion']     = None
            else:
                parameters['anion']     = self.anion_filechooser.get_folder ()
                if not parameters['anion']:
                    self.main.simple_dialog.info(msg = 'Please choose an anion file, or set the anion quantity to 0.')
                    return
            #------------------------------------------------------------------------------
        else:
            parameters['NPositive'] = 0
            parameters['NNegative'] = 0
            parameters['cation']    = None
            parameters['anion']     = None
        
        
        #------------------------------------------------------------------------------
        parameters['reorient'] = self.builder.get_object('checkbox_reorient').get_active()
        #------------------------------------------------------------------------------

        
        #------------------------------------------------------------------------------
        if parameters['add_solvent']:
            solvent_box_id = self.combobox_solvent_system.get_system_id()
            if solvent_box_id is None:
                self.main.simple_dialog.info(msg = 'Please select a pre-built solvent-box system, '
                                                    'or uncheck "Add Solvent Box" to only add ions.')
                return
            if solvent_box_id == system_id:
                self.main.simple_dialog.info(msg = 'The solvent box cannot be the same system as the solute.')
                return
            parameters['solvent'] = self.p_session.psystem[solvent_box_id]
        else:
            parameters['solvent'] = None
        #------------------------------------------------------------------------------
        
        self.label_info.set_text('Running...')
        try:
            self.p_session.solvate_system (e_id = system_id, parameters = parameters)
        except (ValueError, KeyError) as error:
            self.main.simple_dialog.info(msg = str(error))
            self.label_info.set_text('Could not solvate the system -- see the message above.')
            return
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg = 'Unexpected error while solvating the system:\n{}'.format(error))
            self.label_info.set_text('Could not solvate the system -- see the message above.')
            return

        self.label_info.set_text('Done! A new system was added to the treeview.')


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
        if self.Visible == False:
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

            if self.p_session.psystem:
                system  = self.main.p_session.get_system()
                self.combobox_systems.set_active_iter(system.e_liststore_iter)
            
            self.window.show_all()                                               
            self.window.connect('destroy', self.close_window)                                               
            #self.combobox_systems.set_active(0)
            self.Visible    =  True
            '''--------------------------------------------------------------------------------------------'''
        
        #print(idnum, text )
    def close_window (self, button = None, data  = None):
        """ Function doc """
        if not self.Visible:
            return
        self.window.destroy()
        self.Visible    =  False

    def on_button_cancel_clicked (self, widget):
        """ Function doc """
        self.close_window()
        
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
        
        system_id = self.combobox_systems.get_system_id()
        if system_id is None:
            self.main.simple_dialog.info(msg = 'Please select a system.')
            return

        try:
            parameters = {}

            parameters['_Density'] = int(self.builder.get_object('entry_density').get_text())
            parameters['_Steps']   = int(self.builder.get_object('entry_number_of_steps').get_text())
            parameters['_XBox']    = int(self.builder.get_object('entry_size_X').get_text())
            parameters['_YBox']    = int(self.builder.get_object('entry_size_Y').get_text())
            parameters['_ZBox']    = int(self.builder.get_object('entry_size_Z').get_text())
            parameters['_Refine']  = True
        except ValueError:
            self.main.simple_dialog.info(msg = 'Please provide valid numeric values for density, '
                                                'number of steps and box size (X, Y, Z).')
            return

        '''selecting the vismol object from the content that is in the combobox '''
        vobject_id = self.coordinates_combobox.get_vobject_id()
        vobject    = self.main.vm_session.vm_objects_dic[vobject_id]
        
        '''This function imports the coordinates of a vobject into the dynamo system in memory.''' 
        self.main.p_session.set_psystem_coordinates_from_vobject(vobject   = vobject, 
                                                                           system_id = system_id )
        
        parameters['molecule'] = self.main.p_session.psystem[system_id]

        try:
            self.p_session.make_solvent_box(parameters)
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg = 'Unexpected error while building the solvent box:\n{}'.format(error))
            return

        self.close_window(None)
