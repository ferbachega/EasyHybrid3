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


class EasyHybridGoToAtomWindow(Gtk.Window):
    """Window for navigating through systems, residues, and atoms.

    This class provides a GTK3-based interface to select a molecular
    system, view chains, residues, and atoms, and perform visualization 
    actions such as centering and selecting.
    """
    def __init__(self, main = None, system_liststore = None):
        """Initialize the Go-To-Atom window.

        Args:
            main: Reference to the main application object.
            system_liststore: GTK ListStore containing the available systems.
        """
        
        self.main  = main
        self.vm_session    = self.main.vm_session
        self.p_session     = self.main.p_session
        self.system_liststore      = system_liststore
        self.coordinates_liststore = Gtk.ListStore(str, int)
        
        
        self.residue_liststore = Gtk.ListStore(GdkPixbuf.Pixbuf, int, str, str, int)
        self.atom_liststore    = Gtk.ListStore(GdkPixbuf.Pixbuf, int, str, str, float, int, )
        self.residue_filter    = False
        self.visible           = False
        
        self.shift = False
        
        # Dictionary for assigning colors to residues.
        self.residues_dictionary = {
                               'WAT': [165,42,42], 
                               'SOL': [165,42,42], 
                               'HOH': [165,42,42], 
                               
                               'CYS': [207,40,168], 
                               'ASP': [128,0,128], 
                               'SER': [0, 128, 0], 
                               'GLN': [0, 128, 0], 
                               'LYS': [255, 0, 0],
                               'ILE': [30, 144, 255], 
                               'PRO': [255, 255,0], 
                               'THR': [0, 128, 0], 
                               'PHE': [6,176,176], 
                               'ASN': [0, 128, 0], 
                               'GLY': [255, 165,0], 
                               'HIS': [6,176,176], 
                               
                               # amber
                               "HID": [6,176,176],
                               "HIE": [6,176,176],
                               "HIP": [6,176,176],
                               "ASH": [0, 128, 0],
                               "GLH": [0, 128, 0],
                               "CYX": [207,40,168],
                               
                               # charmm
                               "HSD": [6,176,176], 
                               "HSE": [6,176,176], 
                               "HSP": [6,176,176], 
                               
                               
                               'LEU': [30, 144, 255], 
                               'ARG': [255, 0, 0], 
                               'TRP': [6,176,176], 
                               'ALA': [30, 144, 255], 
                               'VAL': [30, 144, 255], 
                               'GLU': [128,0,128], 
                               'TYR': [6,176,176], 
                               'MET': [30, 144, 255]}

    def open_window (self):
        """Open the Go-To-Atom window.

        Builds and displays the interface elements such as system 
        combobox, residue and chain lists, and atom treeviews.
        
        Args:
            main: Reference to the main application object.
            system_liststore: GTK ListStore containing the available systems.
        
        """
       
        #'''
        if self.visible  ==  False:
            
            #self.vm_session.Vismol_Objects_ListStore
            
            #------------------------------------------------------------------#
            #                  SYSTEM combobox and Label
            #------------------------------------------------------------------#
            self.box_vertical    = Gtk.Box(orientation = Gtk.Orientation.VERTICAL,   spacing = 10)
            self.box_horizontal1 = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 10)
             

            
            self.label1  = Gtk.Label()
            self.label1.set_text('System:')
            self.box_horizontal1.pack_start(self.label1, False, False, 0)

            self.combobox_systems = SystemComboBox(self.main)
            self.combobox_systems.connect("changed", self.on_combobox_systems_changed)
            self.box_horizontal1.pack_start(self.combobox_systems, False, False, 0)
            #------------------------------------------------------------------#

            #------------------------------------------------------------------#
            #                  COORDINATES combobox and Label
            #------------------------------------------------------------------#
            self.label1  = Gtk.Label()
            self.label1.set_text('Coordinates:')
            self.box_horizontal1.pack_start(self.label1, False, False, 0)
            self.coordinates_combobox = CoordinatesComboBox(self.main.vobject_liststore_dict[self.main.p_session.active_id])
            
            self.box_horizontal1.pack_start(self.coordinates_combobox, False, False, 0)
            #------------------------------------------------------------------#
            

            #------------------------------------------------------------------#
            #                  CHAIN combobox and Label
            #------------------------------------------------------------------#
            self.box_horizontal2 = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
            
            
            self.label2  = Gtk.Label()
            self.label2.set_text('Chain:')
            self.box_horizontal2.pack_start(self.label2, False, False, 0)

            self.liststore_chains = Gtk.ListStore(str)
            
            self.combobox_chains = Gtk.ComboBox.new_with_model(self.liststore_chains)
            self.combobox_chains.connect("changed", self.on_combobox_chains_changed)
            renderer_text = Gtk.CellRendererText()
            self.combobox_chains.pack_start(renderer_text, True)
            self.combobox_chains.add_attribute(renderer_text, "text", 0)
            self.box_horizontal2.pack_start(self.combobox_chains, False, False, 0)
            

            #------------------------------------------------------------------#
            #                  RESIDUES combobox and Label
            #------------------------------------------------------------------#
            
            self.label3  = Gtk.Label()
            self.label3.set_text('Residue type:')
            self.box_horizontal2.pack_start(self.label3, False, False, 0)

            self.liststore_residues = Gtk.ListStore(str)
            
            self.combobox_residues = Gtk.ComboBox.new_with_model(self.liststore_residues)
            self.combobox_residues.connect("changed", self.on_combobox_residues_changed)
            renderer_text = Gtk.CellRendererText()
            self.combobox_residues.pack_start(renderer_text, True)
            self.combobox_residues.add_attribute(renderer_text, "text", 0)
            self.box_horizontal2.pack_start(self.combobox_residues, False, False, 0)
            
            
            #------------------------------------------------------------------#
            self.treeviewbox_horizontal = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
            #------------------------------------------------------------------------------------------
            
            
            #-----------------------------------------------------------------------------------------
            #                                      Chain filter
            #-----------------------------------------------------------------------------------------
            self.current_filter_chain = None
            # Creating the filter, feeding it with the liststore model
            self.chain_filter = self.residue_liststore.filter_new()
            # setting the filter function, note that we're not using the
            self.chain_filter.set_visible_func(self.chain_filter_func)
            #-----------------------------------------------------------------------------------------
            
            
            self.treeview = Gtk.TreeView(model = self.chain_filter)
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.MULTIPLE)
            
            self.treeview.connect("button-release-event", self.on_treeview_Objects_button_release_event)
            self.treeview.connect("row-activated", self.on_treeview_row_activated_event)
            self.treeview.connect("key-press-event",   self.key_pressed )
            self.treeview.connect("key-release-event", self.key_released)

            for i, column_title in enumerate(
                ['', "index", "Residue",  "Chain", 'size']
            ):
                if i == 0:
                    column = Gtk.TreeViewColumn(column_title, Gtk.CellRendererPixbuf(), pixbuf=0)
                    self.treeview.append_column(column)
                else:
                    renderer = Gtk.CellRendererText()
                    column = Gtk.TreeViewColumn(column_title, renderer, text=i)
                    self.treeview.append_column(column)

            
            self.current_filter_chain = None

            self.scrollable_treelist = Gtk.ScrolledWindow()
            self.scrollable_treelist.set_vexpand(True)
            self.scrollable_treelist.add(self.treeview)
            #------------------------------------------------------------------------------------------
            


            #------------------------------------------------------------------------------------------
            self.treeview_atom = Gtk.TreeView(model =self.atom_liststore)
            self.treeview_atom.connect("button-release-event", self.on_treeview_atom_button_release_event)
            self.treeview_atom.connect("row-activated", self.on_treeview_atom_row_activated_event)

            for i, column_title in enumerate(
                ['', "index", "name", "MM atom", 'MM charge']
            ):
                if i == 0:
                    column = Gtk.TreeViewColumn(column_title, Gtk.CellRendererPixbuf(), pixbuf=0)
                    self.treeview_atom.append_column(column)
                else:
                    renderer = Gtk.CellRendererText()
                    column = Gtk.TreeViewColumn(column_title, renderer, text=i)
                    self.treeview_atom.append_column(column)

            self.scrollable_treelist2 = Gtk.ScrolledWindow()
            self.scrollable_treelist2.set_vexpand(True)
            self.scrollable_treelist2.add(self.treeview_atom)
            #------------------------------------------------------------------------------------------
            
            
            self.box_vertical.pack_start(self.box_horizontal1, False, True, 0)
            self.box_vertical.pack_start(self.box_horizontal2, False, True, 0)
            self.treeviewbox_horizontal.pack_start(self.scrollable_treelist, True, True, 0)
            self.treeviewbox_horizontal.pack_start(self.scrollable_treelist2, True, True, 0)
            
            self.box_vertical.pack_start(self.treeviewbox_horizontal, False, True, 0)
            
            
            self.combobox_systems.set_active(0)
            self.window =  Gtk.Window()
            self.window.set_border_width(10)
            self.window.set_default_size(600, 600)  
            self.window.add(self.box_vertical)
            self.window.connect("destroy", self.close_window)
            self.window.set_title('Go to Atom Window') 
            self.window.show_all() 
                                                          
            #                                                                
            #self.builder.connect_signals(self)                                   
            
            self.visible  =  True
            #----------------------------------------------------------------
        else:
            self.window.present()
        #'''
    def close_window (self, button):
        """ Function doc """
        self.window.destroy()
        self.visible    =  False

    def getColouredPixmap( self, r, g, b, a=255 ):
        """ Given components, return a colour swatch pixmap """
        CHANNEL_BITS=8
        WIDTH=10
        HEIGHT=10
        swatch = GdkPixbuf.Pixbuf.new( GdkPixbuf.Colorspace.RGB, True, CHANNEL_BITS, WIDTH, HEIGHT ) 
        swatch.fill( (r<<24) | (g<<16) | (b<<8) | a ) # RGBA
        return swatch

    def refresh (self, option = 'all'):
        """ Function doc """
        self.update_window()
    
    def update_window (self, system_names = False, coordinates = True,  selections = False ):
        """ Function doc """

        if self.visible:
            
            _id = self.combobox_systems.get_active()
            if _id == -1:
                '''_id = -1 means no item inside the combobox'''
                return None
            else:    
                _, system_id, pixbuf = self.system_liststore[_id]
            
            if system_names:
                self.refresh_system_liststore ()
                self.combobox_systems.set_active(_id)
            
            if coordinates:
                self.refresh_coordinates_liststore ()
                
            
            #if selections:
            #    _, system_id = self.system_liststore[_id]
            #    self.refresh_selection_liststore(system_id)
        else:
            if system_names:
                self.refresh_system_liststore ()
            if coordinates:
                self.refresh_coordinates_liststore ()
            
    def refresh_coordinates_liststore(self, system_id = None):
        """
        Refreshes the coordinates ListStore associated with the currently selected system.

        This method updates the model of the coordinates_combobox with the ListStore corresponding 
        to the system currently selected in the systems combobox. It then resets the active selection 
        in the coordinates_combobox to a default state (-1), which typically means no selection.

        Parameters
        ----------
        system_id : optional
            Not used directly; the selected system ID is obtained from the combobox.
        """
        # Override any passed system_id with the ID obtained from the combobox.
        system_id = self.combobox_systems.get_system_id()

        # Update the coordinates combobox model with the ListStore for the selected system.
        # vobject_liststore_dict  maps system IDs to their associated ListStore objects.
        self.coordinates_combobox.set_model(self.main.vobject_liststore_dict[system_id])

        # Reset the active selection in the coordinates combobox. -1 typically indicates no selection.
        self.coordinates_combobox.set_active_vobject(-1)

    def refresh_system_liststore (self):
        """ Function doc """

    def on_combobox_residues_changed (self, widget):
        """
        Callback triggered when the selected residue in a combobox changes.

        This function retrieves the currently selected residue from the provided combobox widget, 
        updates an internal state variable `current_filter_residue`, and triggers a refiltering 
        of a proxy model (filter) if it exists. This pattern is common in GTK3 when using a filtered 
        TreeModel (Gtk.TreeModelFilter) to dynamically update the visible entries in a TreeView 
        based on selection criteria.

        Parameters
        ----------
        widget : Gtk.ComboBox
            The combobox that emitted the "changed" signal.
        """

        # Retrieve the iterator for the currently selected row.
        tree_iter = widget.get_active_iter()

        # Only proceed if a valid row is selected.
        if tree_iter is not None:
            # Obtain the ListStore (model) backing the combobox.
            model = widget.get_model()

            # Extract the first column's value of the selected row.
            residue = model[tree_iter][0]

            # Update the internal state variable to store the currently selected residue.
            self.current_filter_residue = residue

            # If a residue filter exists, trigger a refiltering of the model.
            # This will automatically update any TreeView or widget displaying filtered data.
            if self.residue_filter:
                self.residue_filter.refilter()

    def on_combobox_chains_changed (self, widget):
        """
        Callback triggered when the selected residue in a combobox changes.

        This function retrieves the currently selected residue from the provided combobox widget, 
        updates an internal state variable `current_filter_residue`, and triggers a refiltering 
        of a proxy model (filter) if it exists. This pattern is common in GTK3 when using a filtered 
        TreeModel (Gtk.TreeModelFilter) to dynamically update the visible entries in a TreeView 
        based on selection criteria.

        Parameters
        ----------
        widget : Gtk.ComboBox
            The combobox that emitted the "changed" signal.
        """
        # Retrieve the iterator for the currently selected row.
        tree_iter = self.combobox_chains.get_active_iter()
        if tree_iter is not None:
            model = self.combobox_chains.get_model()
            chain = model[tree_iter][0]
            #print("Selected: country=%s" % country)
        
        self.current_filter_chain = chain
        #print("%s Chain selected!" % self.current_filter_chain)
        # we update the filter, which updates in turn the view
        self.chain_filter.refilter()

    def _build_chain_liststore (self):
        """ Function doc """
        self.liststore_chains = Gtk.ListStore(str)
        self.liststore_chains.append(['all'])
        chains = self.vobject.chains.keys()

        #self.chain_liststore = Gtk.ListStore(str)

        for chain in chains:
            self.liststore_chains.append([chain])
        self.combobox_chains.set_model(self.liststore_chains)
        self.combobox_chains.set_active(0)
    
    def _build_res_liststore(self):
        self.residue_liststore.clear() #= Gtk.ListStore(GdkPixbuf.Pixbuf, int, str, str, int)
        for chain in self.vobject.chains:
            for index, resi in self.vobject.chains[chain].residues.items():
                #print(resi.index, resi.name, chain,  len(resi.atoms) ) 
                color  =  self.vobject.color_palette['C']
                res_color  = [int(color[0]*255),int(color[1]*255),int(color[2]*255)] 
                swatch = self.getColouredPixmap( res_color[0], res_color[1], res_color[2] )
                
                self.residue_liststore.append(list([swatch, index, resi.name , chain,  len(self.vobject.chains[chain].residues[index].atoms)]))

    def on_combobox_systems_changed (self, widget):
        """ Function doc """
        cb_id = widget.get_system_id()
        
        if cb_id is not None:
            
            self.update_window (coordinates = True)
            
            key =  self.coordinates_combobox.get_vobject_id()
            #_, key = self.coordinates_liststore[cb_id]
            
            self.vobject = self.vm_session.vm_objects_dic[key]
            self._build_chain_liststore()
            self._build_res_liststore()
            
            
            #-----------------------------------------------------------------------------------------
            #                                      Chain filter
            #-----------------------------------------------------------------------------------------
            self.current_filter_chain = None
            # Creating the filter, feeding it with the liststore model
            self.chain_filter = self.residue_liststore.filter_new()
            # setting the filter function, note that we're not using the
            self.chain_filter.set_visible_func(self.chain_filter_func)
            #-----------------------------------------------------------------------------------------
            
            

            #-----------------------------------------------------------------------------------------
            #                                      Residue combobox
            #-----------------------------------------------------------------------------------------
            self.liststore_residues = Gtk.ListStore(str)
            self.liststore_residues.append(['all'])
            
            resn_labels = {}
            for chain in self.vobject.chains.keys():
                for resi, residue in self.vobject.chains[chain].residues.items():
                    resn_labels[residue.name] = True
            
            #for resi, residue in self.vobject.residues.items():
            #    resn_labels[residue.name] = True
            
            for resn in resn_labels.keys():
                #print (resn)
                self.liststore_residues.append([resn])
            
            self.combobox_residues.set_model(self.liststore_residues)
            self.combobox_residues.set_active(0)
            
            #-----------------------------------------------------------------------------------------
            #                                      Residue filter
            #-----------------------------------------------------------------------------------------
            self.current_filter_residue = None
            # Creating the filter, feeding it with the liststore model
            self.residue_filter = self.chain_filter.filter_new()
            # setting the filter function, note that we're not using the
            self.residue_filter.set_visible_func(self.residue_filter_func)
            #-----------------------------------------------------------------------------------------        
            
            self.treeview.set_model(self.residue_filter)
            
    def on_treeview_atom_button_release_event(self, tree, event):
        if event.button == 2:
            selection     = tree.get_selection()
            model         = tree.get_model()
            (model, iter) = selection.get_selected()
            
            
            if iter != None:
                self.selectedID  = int(model.get_value(iter, 1))-1  # @+
                atom = self.vobject.atoms[self.selectedID]
                self.vm_session.vm_glcore.center_on_atom(atom)
       
    def on_treeview_atom_row_activated_event (self, tree, rowline , column):
        """ Function doc """
        selection   = tree.get_selection()
        model       = tree.get_model()
        data        = list(model[rowline])
        pickedID    = data[-1]
                
                
        key =  self.coordinates_combobox.get_vobject_id()
        self.vobject = self.vm_session.vm_objects_dic[key]
                
        atom_picked = self.vobject.atoms[pickedID]
                    
        self.vm_session.selections[self.vm_session.current_selection].selecting_by_atom([atom_picked],  True)
        self.vm_session.selections[self.vm_session.current_selection]._build_selected_atoms_coords_and_selected_objects_from_selected_atoms()
        self.vm_session.vm_glcore.queue_draw()

    def on_treeview_row_activated_event(self, tree, rowline , column ):
        #print (A,B,C)
        selection     = tree.get_selection()
        model         = tree.get_model()
        
        #print(model)
        #print(rowline, list(model[rowline]))
        
        data  = list(model[rowline])
        self.selectedID  = int(data[1])  # @+
        self.selectedObj = str(data[2])
        self.selectedChn = str(data[3])
        
        key =  self.coordinates_combobox.get_vobject_id()
        self.vobject = self.vm_session.vm_objects_dic[key]
        
        res = self.vobject.chains[self.selectedChn].residues[self.selectedID]
        frame = self.vm_session.get_frame ()
        res.get_center_of_mass(frame = frame)
        '''centering and selecting'''

        if self.shift:
            
            for chain in self.vobject.chains.keys():
                for resi, residue in self.vobject.chains[chain].residues.items():
                    if residue.name == res.name:
                        atom_keys = list(residue.atoms.values())
                        self.vm_session._selection_function_set({atom_keys[0]})
        else:
            atom_keys = list(res.atoms.values())
            self.vm_session._selection_function_set({atom_keys[0]})
            
        self.vm_session.vm_glcore.updated_coords = True
        self.vm_session.vm_glcore.queue_draw()
        
    def on_treeview_Objects_button_release_event(self, tree, event):
        #print ( tree, event)
        
        if event.button == 3:
            print ("button 3", event)

        if event.button == 2:
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.SINGLE)

            
            selection     = tree.get_selection()
            model         = tree.get_model()
            (model, iter) = selection.get_selected()

            if iter != None:
                self.selectedID  = int(model.get_value(iter, 1))  # @+
                self.selectedObj = str(model.get_value(iter, 2))
                self.selectedChn = str(model.get_value(iter, 3))
                res = self.vobject.chains[self.selectedChn].residues[self.selectedID]
                frame = self.vm_session.get_frame ()
                res.get_center_of_mass(frame = frame)
                
                
                charges         = list(self.p_session.psystem[self.vobject.e_id].mmState.charges)
                atomTypes       =      self.p_session.psystem[self.vobject.e_id].mmState.atomTypes
                atomTypeIndices = list(self.p_session.psystem[self.vobject.e_id].mmState.atomTypeIndices)
                
                self.vm_session.vm_glcore.center_on_coordinates(res.vm_object, res.mass_center)
        
                self.atom_liststore.clear()
                for atom in res.atoms.values():
                    color  =  self.vobject.color_palette[atom.symbol]
                    swatch = self.getColouredPixmap( int(color[0]*255), int(color[1]*255), int(color[2]*255) )
                    self.atom_liststore.append(list([swatch, int(atom.index), atom.name, atom.symbol ,  charges[atom.index-1] , int(atom.atom_id) ]))
            
            
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.MULTIPLE)

        
        if event.button == 1:
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.SINGLE)

            selection     = tree.get_selection()
            model         = tree.get_model()
            (model, iter) = selection.get_selected()
            
            
            if iter != None:
                self.selectedID  = int(model.get_value(iter, 1))  # @+
                self.selectedObj = str(model.get_value(iter, 2))
                self.selectedChn = str(model.get_value(iter, 3))
                res = self.vobject.chains[self.selectedChn].residues[self.selectedID]
                
                
                self.atom_liststore.clear()
                
                '''when we don't have an MM (molecular mechanics) system, the "system" object doesn't 
                have the mmState attribute, so we need to get the information of atom names and indexes
                 in another way.'''
                #--------------------------------------------------------------------------------------
                try:
                    charges         = list(self.p_session.psystem[self.vobject.e_id].mmState.charges)
                    atomTypes       =      self.p_session.psystem[self.vobject.e_id].mmState.atomTypes
                    atomTypeIndices = list(self.p_session.psystem[self.vobject.e_id].mmState.atomTypeIndices)
                except:
                    charges         = [0.0]*len(self.p_session.psystem[self.vobject.e_id].atoms)
                    atomTypes       = []
                    atomTypeIndices = []
                    for atom in self.p_session.psystem[self.vobject.e_id].atoms.items:
                        atomTypes.append('-')
                        atomTypeIndices.append(atom.index)
                #--------------------------------------------------------------------------------------
                for atom in res.atoms.values():
                    color  =  self.vobject.color_palette[atom.symbol]
                    swatch = self.getColouredPixmap( int(color[0]*255), int(color[1]*255), int(color[2]*255) )
                    self.atom_liststore.append(list([swatch, int(atom.index), str(atom.name), str(atomTypes[atomTypeIndices[atom.index-1]]) ,  charges[atom.index-1] , int(atom.atom_id)]))

            self.treeview.get_selection().set_mode(Gtk.SelectionMode.MULTIPLE)

    def on_chk_renderer_toggled(self, cell, path, model):
        print('on_chk_renderer_toggled -> model[path][0]: ', model[path][0])

    def residue_filter_func(self, model, iter, data):
        """Tests if the language in the row is the one in the filter"""
        if (
            self.current_filter_residue is None
            or self.current_filter_residue == "all"
        ):
            return True
        else:
            return model[iter][2] == self.current_filter_residue

    def chain_filter_func(self, model, iter, data):
        """Tests if the language in the row is the one in the filter"""
        if (
            self.current_filter_chain is None
            or self.current_filter_chain == "all"
        ):
            return True
        else:
            return model[iter][3] == self.current_filter_chain

    def update (self):
        """ Function doc """
        #print('VismolGoToAtomWindow2 update')
        pass

    def key_pressed  (self, widget, event):
        key = Gdk.keyval_name(event.keyval)
        
        if key == 'Shift_R' or key == 'Shift_L':
            self.shift = True
        print(widget, event, Gdk.keyval_name(event.keyval), self.shift)
    
    def key_released (self, widget, event):
        key = Gdk.keyval_name(event.keyval)
        if key == 'Shift_R' or key == 'Shift_L':
            self.shift = False
        print(widget, event, Gdk.keyval_name(event.keyval), self.shift)
