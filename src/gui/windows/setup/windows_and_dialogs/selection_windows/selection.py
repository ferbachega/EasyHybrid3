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


class EasyHybridSelectionWindow:
    """Window for managing and executing atom/residue/molecule selections.

    This class provides a GTK-based window that allows users to configure and
    run selection operations (Expand, Around, Complete, etc.) in the EasyHybrid
    environment. The user can choose the selection type, select by atom/residue/
    molecule/chain, and set a selection radius.
    """
    def __init__(self, main = None):
        """Initialize the selection window.

        Args:
            main (object, optional): Reference to the main EasyHybrid application
                instance. Used to access sessions, home directory, and utilities.
        """
        self.main = main
        self.vm_session      = main.vm_session
        self.Visible         = False        
        self.home            = main.home
        self.p_session       = main.p_session
        
        
        self.chain = ''
        self.resn  = ''
        self.resi  = ''
        self.atom  = ''

        self._type_dict={
                        0 : "Expand"   ,
                        1 : "Around"   ,
                        2 : "Complete" ,
                        3 : "ByComplement" ,
                        }
        
        self.select_by_dict ={
                             0 : "Atom"      ,
                             1 : "Residue"   ,
                             2 : "Molecule"  ,
                             4 : "Chain"    ,
                             }
        
        
    def open_window (self):
        """Open the selection window.

        Loads the GTK interface from a Glade file, builds the combo boxes for
        selection type and selection target, configures spin buttons, and shows
        the window. If already open, brings it to the foreground.
        """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/setup/modify_selection_window.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('modify_selection_window')
            self.window.set_title('Modify Selection Window')
            self.window.set_keep_above(True)
             
            
            self.box_combo_methods = self.builder.get_object('box_combo_methods')
            self.box_select_by     = self.builder.get_object('box_select_by')
            
            method_store = Gtk.ListStore(str)
            for key, item in self._type_dict.items():
                method_store.append([item])


            #------------------------------------------------------------------
            self.method_combo = Gtk.ComboBox.new_with_model(method_store)
            renderer_text = Gtk.CellRendererText()
            self.method_combo.pack_start(renderer_text, True)
            self.method_combo.add_attribute(renderer_text, "text", 0)            
            
            self.method_combo.set_entry_text_column(0)
            self.method_combo.set_active(0)
            self.box_combo_methods.pack_start(self.method_combo, False, False, 0)
            #------------------------------------------------------------------

            
            
            select_by_store = Gtk.ListStore(str)
            for key, item in self.select_by_dict.items():
                select_by_store.append([item])


            #------------------------------------------------------------------
            self.select_by_combo = Gtk.ComboBox.new_with_model( select_by_store)
            renderer_text = Gtk.CellRendererText()
            self.select_by_combo.pack_start(renderer_text, True)
            self.select_by_combo.add_attribute(renderer_text, "text", 0)            
            
            self.select_by_combo.set_entry_text_column(0)
            self.select_by_combo.set_active(1)
            self.box_select_by.pack_start(self.select_by_combo, False, False, 0)
            #------------------------------------------------------------------
            
            
            
            #------------------------------------------------------------------
            self.radius_spinbutton  = self.builder.get_object('radius_spinbutton' )
            #------------------------------------------------------------------
            self.radius_adjustment = Gtk.Adjustment(value          = 14 , 
                                                 upper          = 100, 
                                                 step_increment = 1  , 
                                                 page_increment = 10 )

            self.radius_spinbutton.set_adjustment ( self.radius_adjustment)
            #------------------------------------------------------------------


            self.window.show_all()
            self.Visible  = True
    
        else:
            self.window.present()
    

    def close_window (self, button, data  = None):
        """Close the selection window.

        Args:
            button (Gtk.Button): The button that triggered the event.
            data (Any, optional): Additional signal data.
        """
        self.window.destroy()
        self.Visible    =  False
    
    
    def run_selection (self, button):
        """Run the selection operation based on user input.

        Collects parameters from the UI (selection type, target, radius) and
        executes an advanced selection via `vm_session`.

        Args:
            button (Gtk.Button): The button that triggered the event.
        """
        #self.method_combo
        #self.select_by_combo

        _radius      =  self.radius_spinbutton.get_value ()
        _type        =  self.method_combo.get_active()
        _select_by   =  self.select_by_combo.get_active()
        

        
        _type      = self._type_dict    [_type]
        _select_by = self.select_by_dict[_select_by]



        TrueFalse, msg = self.vm_session.advanced_selection( selection        = None,
                                                                 _type        = _type ,
                                                                 selecting_by = _select_by,
                                                                 radius        = _radius,  
                                                                 grid_size     = _radius)
        if TrueFalse:
            print(msg)
            pass
            #self.main.simple_dialog.info(msg = msg )
        else:
            self.main.simple_dialog.error(msg = msg )


class PDynamoSelectionWindow:
    """ Class doc """
    def __init__(self, main = None):
        """ Class initialiser """
        self.main = main
        self.vm_session      = main.vm_session
        self.Visible         = False        
        self.home            = main.home
        self.p_session       = main.p_session
        
        
        self.chain = ''
        self.resn  = ''
        self.resi  = ''
        self.atom  = ''

 
    def open_window (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/setup/pDynamo_selection.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('pdynamo_selection_window')
            self.window.set_title('pDynamo Selection Window')
            self.window.set_keep_above(True)
            
            #self.chain_entry  = self.builder.get_object('chain_entry')
            #self.resn_entry   = self.builder.get_object('resn_entry' )
            #self.resi_entry   = self.builder.get_object('resi_entry' )
            #self.atom_entry   = self.builder.get_object('atom_entry' )
            self.builder.get_object('chain_entry').set_text(str(self.chain) )
            self.builder.get_object('resn_entry' ).set_text(str(self.resn ) )
            self.builder.get_object('resi_entry' ).set_text(str(self.resi ) )
            self.builder.get_object('atom_entry' ).set_text(str(self.atom ) )
            
            
            self.box_combo_methods = self.builder.get_object('box_combo_methods')
            
            #self.self.method_combobox   = self.builder.get_object('self.method_combobox'   )
            #self.radius_spinbutton = self.builder.get_object('radius_spinbutton' )
            #self.action_combobox   = self.builder.get_object('action_combobox'   )

            method_store = Gtk.ListStore(str)
            method_store.append(["ByComponent"])
            method_store.append(["Complement"])
            method_store.append(["ByAtom"])


            #------------------------------------------------------------------
            self.method_combo = Gtk.ComboBox.new_with_model(method_store)
            renderer_text = Gtk.CellRendererText()
            self.method_combo.pack_start(renderer_text, True)
            self.method_combo.add_attribute(renderer_text, "text", 0)            
            
            self.method_combo.set_entry_text_column(0)
            self.method_combo.set_active(0)
            self.box_combo_methods.pack_start(self.method_combo, False, False, 0)
            #------------------------------------------------------------------

            
            
            
            
            self.radius_spinbutton  = self.builder.get_object('radius_spinbutton' )
            #------------------------------------------------------------------
            self.fps_adjustment = Gtk.Adjustment(value          = 24 , 
                                                 upper          = 100, 
                                                 step_increment = 1  , 
                                                 page_increment = 10 )

            self.radius_spinbutton.set_adjustment ( self.fps_adjustment)
            #------------------------------------------------------------------




            '''
            self.box_combo_action = self.builder.get_object('box_combo_action' )
            
            action_store = Gtk.ListStore(str)
            action_store.append(["ByComponent"])
            action_store.append(["Complement"])
            action_store.append(["ByLinearPolymer"])
            #------------------------------------------------------------------
            self.action_combo = Gtk.ComboBox.new_with_model(action_store)
            renderer_text = Gtk.CellRendererText()
            self.action_combo.pack_start(renderer_text, True)
            self.action_combo.add_attribute(renderer_text, "text", 0)            
            
            self.action_combo.set_entry_text_column(0)
            self.action_combo.set_active(0)
            
            self.box_combo_action.pack_start(self.action_combo, False, False, 0)
            #------------------------------------------------------------------
            '''
                
                
                

            self.window.show_all()
            self.Visible  = True
    
        else:
            self.window.present()
    
    def import_data (self, button):
        """ Function doc """
        entry_name    = None
        idnum     = self.combobox_pdynamo_system.get_active()
        text      = self.combobox_pdynamo_system.get_active_text()
        
        #print(idnum, text )
    def close_window (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False
    
    
    def run_selection (self, button):
        """ Function doc """
        #print('run_selection', self.method_combo.get_active(), self.radius_spinbutton.get_value ())
        self.chain = self.builder.get_object('chain_entry').get_text()
        self.resn  = self.builder.get_object('resn_entry' ).get_text()
        self.resi  = self.builder.get_object('resi_entry' ).get_text()
        self.atom  = self.builder.get_object('atom_entry' ).get_text()
        
        #print(chain,resn,resi, atom)
        #'''
        #atom1 = self.vm_session.picking_selections.picking_selections_list[0]
        #print (atom1.chain, atom1.resn, atom1.resi, atom1.name)
        #print ("%s:%s.%s:%s" %(chain,resn,resi,atom))
        
        _centerAtom = "%s:%s.%s:%s" %(self.chain, 
                                      self.resn,
                                      self.resi,
                                      self.atom)
        _radius     =  self.radius_spinbutton.get_value ()
        _method     =  self.method_combo.get_active()
        
        #self.main.p_session.selections (_centerAtom, _radius, _method )
        indexes = self.main.p_session.p_selections (system_id = None, 
                                                  _centerAtom = _centerAtom, 
                                                      _radius = _radius, 
                                                      _method = _method)
        
        #print (indexes)
        self.main.vm_session.selections[self.vm_session.current_selection].selecting_by_indexes(self.atom1.vm_object, 
                                                           indexes, 
                                                           clear=True)
        if self.main.vm_session.selection_box_frame:
            self.main.vm_session.selection_box_frame.change_toggle_button_selecting_mode_status(False)
        else:
            self.main.vm_session._picking_selection_mode = False
        
        self.main.vm_session.selections[self.main.vm_session.current_selection].active = True
        self.main.vm_session.vm_glcore.queue_draw()
        
        
    def get_info_from_selection (self, button):
        """ Function doc """
        #chain = self.builder.get_object('chain_entry').get_text()
        #resn  = self.builder.get_object('resn_entry' ).get_text()
        #resi  = self.builder.get_object('resi_entry' ).get_text()
        #atom  = self.builder.get_object('atom_entry' ).get_text()
        
        atom1 = self.vm_session.picking_selections.picking_selections_list[0]
        self.atom1 = atom1
        if atom1:
            #atom1.chain, atom1.resn, atom1.resi, atom1.name
            
            if atom1.chain.name =='':
                chain = '*'
            else:
                chain = atom1.chain.name
            #residue = atom1.vm_object.residues[atom1.residue]
            self.builder.get_object('chain_entry').set_text(str(chain )      )
            self.builder.get_object('resn_entry' ).set_text(str(atom1.residue.name) )
            self.builder.get_object('resi_entry' ).set_text(str(atom1.residue.index) )
            self.builder.get_object('atom_entry' ).set_text(str(atom1.name) )
        else:
            print('use picking selection to chose the central atom')
            
    def update (self):
        """ Function doc """
        #print('VismolGoToAtomWindow2 update')
        pass
