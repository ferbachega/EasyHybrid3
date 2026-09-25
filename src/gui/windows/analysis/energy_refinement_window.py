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
from util.debug import dprint
import gi
import os
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk, Pango

from gui.widgets.custom_widgets import FolderChooserButton
#from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.widgets.custom_widgets import ReactionCoordinateBox
from gui.widgets.custom_widgets import AdvancedReactionCoordinateBox
from pdynamo.LogFileWriter import detect_reaction_coordinates_from_log
from pprint import pprint

class EnergyRefinementWindow():

    def __init__(self, main = None ):
        """ Class initialiser """
        self.main       = main
        self.home       = main.home
        self.p_session  = main.p_session
        self.vm_session = main.vm_session
        self.Visible    =  False        
        
        self.vobject_liststore   = Gtk.ListStore(str, int)
        self.data_liststore      = Gtk.ListStore(str, int)

        # Backing stores for the "Advanced mode" (weighted distance list) RC
        # boxes -- same pattern as PotentialEnergyScanWindow/
        # UmbrellaSamplingWindow. Columns: atom1 name, atom1 index,
        # atom2 name, atom2 index, weight, dist.
        self.rc_liststore1 = Gtk.ListStore(str, str, str, str, str, str)
        self.rc_liststore2 = Gtk.ListStore(str, str, str, str, str, str)
        self.rc_liststore1.connect("row-inserted", self.on_row_inserted)
        self.rc_liststore1.connect("row-deleted", self.on_row_deleted)
        self.rc_liststore2.connect("row-inserted", self.on_row_inserted)
        self.rc_liststore2.connect("row-deleted", self.on_row_deleted)

        self.input_types_liststore = Gtk.ListStore(str)
        self.input_types = {
                        0:'Vobject'    , #  single file
                        1:'pklfolder - pDynamo Trajectory'  , #  trajectory
                        2:'pklfolder2D - pDynamo 2D Trajectory', #  2d trajectory  
                      # 3:'pdbfile'    ,
                      # 4:'pdbfolder'  ,
                      # 5:'dcd',
                      # 6:'crd',
                      # 7:'xyz',
                      # 8:'mol2',
                      # 9:'netcdf',
                      #10:'log_file'  
                           }
        for key, input_type in self.input_types.items():
                self.input_types_liststore.append([input_type])

    def open_window (self, vobject = None):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'src/gui/windows/analysis/energy_refinement_window.glade'))
            self.builder.connect_signals(self)
            #self.vobject = vobject
            self.window = self.builder.get_object('window')
            self.window.set_title('Energy Refinement Window')
            self.window.connect('destroy', self.close_window)
            
            
            
            
            
            
            #---------------------------------------------------------------------------------
            #                           INPUT TYPE COMBOBOX
            #---------------------------------------------------------------------------------
            self.comobobox_input_type = self.builder.get_object('comobobox_input_type')
            
            
            self.comobobox_input_type.set_model(self.input_types_liststore)
            self.comobobox_input_type.connect("changed", self.on_input_types_changed)
            self.comobobox_input_type.set_model(self.input_types_liststore)
            
            renderer_text = Gtk.CellRendererText()
            self.comobobox_input_type.pack_start(renderer_text, True)
            self.comobobox_input_type.add_attribute(renderer_text, "text", 0)
            #---------------------------------------------------------------------------------
           
           
           
            
            #----------------------------------------------------------------------------------------------
            #                                Starting Coordinates ComboBox  
            #----------------------------------------------------------------------------------------------
            self.box_coordinates = self.builder.get_object('box_coordinates')
            self.combobox_starting_coordinates = CoordinatesComboBox() #self.builder.get_object('coordinates_combobox')
            self.combobox_starting_coordinates.connect("changed", self.on_coordinates_combobox_change)
            self.box_coordinates.pack_start(self.combobox_starting_coordinates, False, False, 0)
            self._starting_coordinates_model_update(init = True)
            #----------------------------------------------------------------------------------------------
            
            
            
            #----------------------------------------------------------------------------------------------
            #                                    Reaction Coordinate Boxes
            #'''--------------------------------------------------------------------------------------------'''
            # "simple" (fixed-shape: distance/multiple-distance/dihedral) and
            # "advanced" (arbitrary weighted distance list) share the same
            # Alignment placeholder; only one of each pair is visible at a
            # time (see on_advanced_mode_toggled) -- same pattern as
            # PotentialEnergyScanWindow/UmbrellaSamplingWindow.
            self.RC_box1 = ReactionCoordinateBox(self.main)
            self.RC_box1_adv = AdvancedReactionCoordinateBox(main = self.main, liststore = self.rc_liststore1)
            rc1_box = Gtk.Box(orientation = Gtk.Orientation.VERTICAL)
            rc1_box.pack_start(self.RC_box1, True, True, 0)
            rc1_box.pack_start(self.RC_box1_adv, True, True, 0)
            self.builder.get_object('rc1_aligment').add(rc1_box)

            self.RC_box2 = ReactionCoordinateBox(self.main)
            self.RC_box2_adv = AdvancedReactionCoordinateBox(main = self.main, liststore = self.rc_liststore2)
            rc2_box = Gtk.Box(orientation = Gtk.Orientation.VERTICAL)
            rc2_box.pack_start(self.RC_box2, True, True, 0)
            rc2_box.pack_start(self.RC_box2_adv, True, True, 0)
            self.builder.get_object('rc2_aligment').add(rc2_box)
            self.build_advanced_treeviews()
            #'''--------------------------------------------------------------------------------------------'''






            '''--------------------------------------------------------------------------------------------'''     
            self.folder_chooser_button = FolderChooserButton(main =  self.main, sel_type = 'folder', home =  self.home)
            self.builder.get_object('folder_chooser_box').pack_start(self.folder_chooser_button.btn, True, True, 0)
            system_id      = self.p_session.active_id
            
            #------------------------------------------------------------------------------------------------
            if self.main.p_session.psystem[self.p_session.active_id]:
                if self.main.p_session.psystem[self.p_session.active_id].e_working_folder == None:
                    # NOTE: 'HOME' was an undefined name here (NameError bug) - it should
                    # refer to the window's own home path, exactly like self.home is used
                    # everywhere else in this file.
                    folder = self.home
                else:
                    folder = self.main.p_session.psystem[self.p_session.active_id].e_working_folder
                self.folder_chooser_button.set_folder(folder = folder)
            #------------------------------------------------------------------------------------------------
            
            
            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system('energy_ref')
                self.builder.get_object('entry_logfile').set_text(output_name)  
            

            self.builder.get_object('checkbtn_advanced_mode').connect('toggled', self.on_advanced_mode_toggled)
            self.builder.get_object('filechooser_file_folder').connect('file-set', self.on_trajectory_data_folder_selected)

            self.window.show_all()

            self.RC_box1.set_hide_scan_parameters()
            self.RC_box2.set_hide_scan_parameters()
            self.RC_box1_adv.set_hide_scan_parameters()
            self.RC_box2_adv.set_hide_scan_parameters()

            self.RC_box1.set_rc_type(0)
            self.RC_box2.set_rc_type(0)

            self.comobobox_input_type.set_active(0)

            self.combobox_starting_coordinates.set_active(0)
            self.builder.get_object('button_cancel').connect('clicked', self.close_window)
            self.builder.get_object('button_run').connect('clicked', self.on_button_run_clicked)

            self.on_advanced_mode_toggled(None)

            self.Visible  = True
    
        else:
            pass


    def _starting_coordinates_model_update (self, init = False):
        """ Function doc """
        #------------------------------------------------------------------------------------
        '''The combobox accesses, according to the id of the active system, 
        listostore of the dictionary object_list state_dict'''
        if self.Visible:

            e_id = self.main.p_session.active_id 
            self.combobox_starting_coordinates.set_model(self.main.vobject_liststore_dict[e_id])
            #------------------------------------------------------------------------------------
            size = len(self.main.vobject_liststore_dict[e_id])
            self.combobox_starting_coordinates.set_active(size-1)
            #------------------------------------------------------------------------------------
        else:
            if init:
                e_id = self.main.p_session.active_id 
                self.combobox_starting_coordinates.set_model(self.main.vobject_liststore_dict[e_id])
                #------------------------------------------------------------------------------------
                size = len(self.main.vobject_liststore_dict[e_id])
                self.combobox_starting_coordinates.set_active(size-1)
                #------------------------------------------------------------------------------------
            else:
                pass


    def on_input_types_changed(self, widget):
        """ Function doc """
        _id = self.comobobox_input_type.get_active()
        #print(_id)
        
        if _id == 0:
            self.builder.get_object('label_coordinates').show()
            self.combobox_starting_coordinates.show()
            
            self.builder.get_object('label_file_folder').hide()
            self.builder.get_object('filechooser_file_folder').hide()
            self.builder.get_object('box_reaction_coordinate2').set_sensitive(False)
        
        elif _id in [1,2]:
            self.builder.get_object('label_coordinates').hide()
            self.combobox_starting_coordinates.hide()
            
            self.builder.get_object('label_file_folder').show()
            self.builder.get_object('filechooser_file_folder').show()
            
            if _id ==1:
                self.builder.get_object('box_reaction_coordinate2').set_sensitive(False)
            elif _id ==2:
                self.builder.get_object('box_reaction_coordinate2').set_sensitive(True)
            else:
                pass

    def on_advanced_mode_toggled (self, widget):
        """ Switches both RC boxes between "simple" (fixed-shape) and
        "advanced" (weighted distance list) mode -- same pattern as
        PotentialEnergyScanWindow/UmbrellaSamplingWindow. """
        advanced = self.builder.get_object('checkbtn_advanced_mode').get_active()

        self.RC_box1.set_visible(not advanced)
        self.RC_box1_adv.set_visible(advanced)
        self.RC_box2.set_visible(not advanced)
        self.RC_box2_adv.set_visible(advanced)

        # [EN] BUG FIX: the window used to be fixed-size (resizable=False in
        # the glade), so the advanced boxes' weighted-distance treeviews
        # (6 columns each, x2 for RC1/RC2) had no room to lay out and the
        # whole window looked squeezed/too narrow -- reported by the user.
        # The window is now resizable (see the .glade), and this widens it
        # to a size that comfortably fits both treeviews the first time
        # advanced mode is turned on, instead of making the user manually
        # drag it wider every time.
        if advanced:
            height = self.window.get_allocated_height()
            self.window.resize(1150, height)

    def build_advanced_treeviews (self):
        """ Builds the two weighted-distance-list treeviews backing the
        "advanced" RC boxes (see AdvancedReactionCoordinateBox.get_rc_data
        for the row shape) -- ported from PotentialEnergyScanWindow. """
        self.treeview1 = Gtk.TreeView(model = self.rc_liststore1)
        self.treeview2 = Gtk.TreeView(model = self.rc_liststore2)

        columns = {
            "atm1": 0,
            "idx1": 1,
            "atm2": 2,
            "idx2": 3,
            "weight": 4,
            "dist": 5,
        }
        for title, i in columns.items():
            renderer = Gtk.CellRendererText()
            renderer.set_property("editable", True)
            renderer.connect("edited", self.on_cell_edited1, i)
            self.treeview1.append_column(Gtk.TreeViewColumn(title, renderer, text=i))

        for title, i in columns.items():
            renderer = Gtk.CellRendererText()
            renderer.set_property("editable", True)
            renderer.connect("edited", self.on_cell_edited2, i)
            self.treeview2.append_column(Gtk.TreeViewColumn(title, renderer, text=i))

        self.treeview1.connect("key-press-event", self.on_key_press)
        self.treeview2.connect("key-press-event", self.on_key_press)

        self.RC_box1_adv.scrolledbox.add(self.treeview1)
        self.RC_box2_adv.scrolledbox.add(self.treeview2)
        self.RC_box1_adv.treeview = self.treeview1
        self.RC_box2_adv.treeview = self.treeview2

    def on_key_press (self, widget, event):
        """ Deletes the selected row of an advanced RC treeview on Delete. """
        if event.keyval == Gdk.KEY_Delete:
            selection = widget.get_selection()
            model, treeiter = selection.get_selected()
            if treeiter is not None:
                model.remove(treeiter)
            return True
        return False

    def on_row_inserted (self, model, path, iter):
        """ Function doc """
        if model is self.rc_liststore1:
            self.RC_box1_adv.refresh_dmininum()
        else:
            self.RC_box2_adv.refresh_dmininum()

    def on_row_deleted (self, model, path):
        """ Function doc """
        if model is self.rc_liststore1:
            self.RC_box1_adv.refresh_dmininum()
        else:
            self.RC_box2_adv.refresh_dmininum()

    def on_cell_edited1 (self, widget, path, new_text, column_index):
        """ Function doc """
        self.rc_liststore1[path][column_index] = new_text

    def on_cell_edited2 (self, widget, path, new_text, column_index):
        """ Function doc """
        self.rc_liststore2[path][column_index] = new_text

    def on_trajectory_data_folder_selected (self, widget):
        """ Auto-detects the reaction coordinate(s) used to produce the
        trajectory in a freshly-picked pklfolder/pklfolder2D data folder,
        by parsing its own 'output.log' -- same mechanism as
        UmbrellaSamplingWindow.on_trajectory_folder_selected (see
        LogFileWriter.detect_reaction_coordinates_from_log for the
        recognized formats: RelaxedSurfaceScan/AdvancedRelaxedSurfaceScan's
        PES scans, or this window's own EnergyRefinement runs) -- so the
        user doesn't have to redefine the same atoms that already define
        this trajectory's reaction coordinate.

        Connected to filechooser_file_folder's own "file-set" signal (a
        native GtkFileChooserButton, unlike the other folder pickers in
        this window which go through the custom FolderChooserButton
        class), so the callback receives the WIDGET itself rather than a
        path string, unlike UmbrellaSamplingWindow's own version of this.

        Silently does nothing if output.log doesn't exist or isn't a
        recognizable EasyHybrid scan log -- this is a convenience, not a
        requirement; the user can still enter the RC by hand exactly as
        before.
        """
        folder = widget.get_filename()
        if not folder:
            return

        detected = detect_reaction_coordinates_from_log(os.path.join(folder, 'output.log'))
        if not detected:
            return

        rc1 = detected.get('RC1')
        rc2 = detected.get('RC2')
        if not (rc1 or rc2):
            return

        advanced = (rc1 or rc2)['rc_type'] == 'advanced'
        self.builder.get_object('checkbtn_advanced_mode').set_active(advanced)
        self.on_advanced_mode_toggled(None)

        if rc1:
            if advanced:
                self._apply_detected_advanced_rc(self.rc_liststore1, rc1)
            else:
                self._apply_detected_simple_rc(self.RC_box1, rc1)

        if rc2:
            # RC2 only actually gets used downstream when input_type is
            # pklfolder2D (see on_button_run_clicked / is_2D_xy) -- switch
            # to it automatically so the auto-filled RC2 isn't silently
            # ignored because the combobox was still left at pklfolder(1D).
            self.comobobox_input_type.set_active(2)
            if advanced:
                self._apply_detected_advanced_rc(self.rc_liststore2, rc2)
            else:
                self._apply_detected_simple_rc(self.RC_box2, rc2)

    def _apply_detected_simple_rc (self, rc_box, data):
        """ Fills ONLY a ReactionCoordinateBox's coordinate-type combobox
        and atom index/name entries from a
        detect_reaction_coordinates_from_log() result -- same approach as
        UmbrellaSamplingWindow._apply_detected_simple_rc. Doesn't touch
        force_constant/nsteps/dincre/dminimum at all, which is moot here
        anyway since RC_box1/RC_box2's scan-parameter grid is already
        hidden via set_hide_scan_parameters() (see open_window). """
        type_map = {
            'simple_distance': 0,
            'multiple_distance': 1,
            'multiple_distance*4atoms': 2,
            'dihedral': 3,
        }
        rc_box.combobox_reaction_coord1.set_active(type_map.get(data['rc_type'], 0))
        rc_box.change_cb_coordType1(rc_box.combobox_reaction_coord1)

        for i, (atom_index, atom_name) in enumerate(zip(data['ATOMS'], data['ATOM_NAMES']), start=1):
            entry_index = rc_box.builder.get_object('entry_atom{}_index_coord1'.format(i))
            entry_name  = rc_box.builder.get_object('entry_atom{}_name_coord1'.format(i))
            if entry_index:
                entry_index.set_text(str(atom_index))
            if entry_name:
                entry_name.set_text(str(atom_name))

    def _apply_detected_advanced_rc (self, liststore, data):
        """ Fills ONLY an AdvancedReactionCoordinateBox's weighted-distance
        treeview (via its backing liststore) from a
        detect_reaction_coordinates_from_log() result -- same approach as
        UmbrellaSamplingWindow._apply_detected_advanced_rc. """
        liststore.clear()
        for row in data['RC']:
            liststore.append(list(row))

    def _get_selected_vobject_info (self):
        """ Resolves the vismol object currently selected in
        combobox_starting_coordinates, along with its idx_2D_xy attribute
        (if any). Returns (vobject, idx_2D_xy) or (None, None) if nothing
        is selected.

        This used to be duplicated (identically) in on_coordinates_combobox_change
        and on_button_run_clicked - now both call this single implementation.
        """
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is None:
            return None, None

        model = self.combobox_starting_coordinates.get_model()
        name, vobject_id = model[tree_iter][:2]

        vobject   = self.main.vm_session.vm_objects_dic[vobject_id]
        idx_2D_xy = getattr(vobject, 'idx_2D_xy', None)

        return vobject, idx_2D_xy

    def on_coordinates_combobox_change (self, widget):
        """ Function doc """
        vobject, idx_2D_xy = self._get_selected_vobject_info()

        if vobject is not None:
            dprint("is_2D_xy" , bool(idx_2D_xy))
            dprint("idx_2D_xy", idx_2D_xy)

            self.builder.get_object('box_reaction_coordinate2').set_sensitive(bool(idx_2D_xy))








    def close_window (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False


    def on_button_run_clicked (self, button):
        """ Function doc """
        parameters = {"simulation_type":"Energy_Refinement"                ,
                      "dialog":True                                        ,                      
                      "system"     : self.p_session.psystem[self.p_session.active_id],
                      "system_name": self.p_session.psystem[self.p_session.active_id].label,
                      "initial_coordinates": None                          ,                       
                      "traj_type":'pklfolder'                              ,
                      'ignore_mm_charges': False                           ,
                      "NmaxThreads":1                                      ,
                      "show":False                                         }

        
        parameters['folder'] = self.folder_chooser_button.get_folder()
        
        input_type = self.comobobox_input_type.get_active()
        dprint('_type: ',input_type)
        #----------------------------------------------------------------------
        if input_type == 0:
            parameters["traj_type"] = 'vobject'

            vobject, idx_2D_xy = self._get_selected_vobject_info()
            if vobject is not None:
                parameters["is_2D_xy"]  = bool(idx_2D_xy)
                parameters["idx_2D_xy"] = idx_2D_xy if idx_2D_xy else False

                parameters["trajectory"] = vobject.frames
                parameters["filename"] = self.builder.get_object('entry_logfile').get_text()



        elif input_type in [1,2]:
            if input_type == 1:
                parameters["traj_type"]  = 'pklfolder'
                parameters["is_2D_xy"] = False
            else:
                parameters["traj_type"]  = 'pklfolder2D'
                parameters["is_2D_xy"] = True
            
            parameters["data_path"] = self.builder.get_object('filechooser_file_folder').get_filename()
            parameters["filename"] = self.builder.get_object('entry_logfile').get_text()
            files = os.listdir( parameters['data_path'])
            pkl_files = []
            
            for _file in files:
                # Check if the file is a text file
                if _file.endswith('.pkl'):
                    pkl_files.append(_file)

            dprint ('pDynamo pkl folder:' , parameters['traj_type'])
            dprint ('Number of pkl files:', len(pkl_files))
            parameters["trajectory"] = pkl_files
        
        else:
            pass
        #----------------------------------------------------------------------
        
        if self.builder.get_object('check_box_MM_chrg_to_zero').get_active():
            parameters["ignore_mm_charges"] = True
        else:
            pass
        
        # "Advanced mode" (weighted distance list) -- same pattern as
        # PotentialEnergyScanWindow/UmbrellaSamplingWindow, but since Energy
        # Refinement only ever MEASURES the reaction coordinate (never runs
        # a restrained scan with it), the advanced box's own get_rc_data()
        # dict (RC/nsteps/dminimum/force_constant/dincre, no "rc_type" key)
        # is tagged with rc_type='advanced' here so compute_reaction_coordinate()
        # in p_methods/energy.py can dispatch on it exactly like the other
        # three (simple_distance/multiple_distance/multiple_distance*4atoms)
        # RC shapes.
        advanced = self.builder.get_object('checkbtn_advanced_mode').get_active()
        rc_box1 = self.RC_box1_adv if advanced else self.RC_box1

        parameters["RC1"] = rc_box1.get_rc_data()
        if advanced:
            parameters["RC1"]["rc_type"] = "advanced"

        #if self.builder.get_object('label_check_button_reaction_coordinate2').get_active():
        #if parameters["traj_type"]  == 'pklfolder2D':

        if parameters["is_2D_xy"]:
            rc_box2 = self.RC_box2_adv if advanced else self.RC_box2
            parameters["RC2"] = rc_box2.get_rc_data()
            if advanced:
                parameters["RC2"]["rc_type"] = "advanced"
            parameters["NmaxThreads"] =  int(self.builder.get_object('n_CPUs_spinbutton').get_value())
            #parameters["traj_type"]   = 'pklfolder2D'

        else:
            parameters["RC2"] = None
            parameters["NmaxThreads"] = 1
            #parameters["traj_type"]   = 'pklfolder'
        
        #----------------------------------------------------------------------
        
        #pprint (parameters)
        self.p_session.run_simulation( parameters = parameters )


    def update (self):
        """ Function doc """
        #print('VismolGoToAtomWindow2 update')
        self._starting_coordinates_model_update()
        if self.Visible:
            self.update_working_folder_chooser()
            
            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system('energy_ref')
                self.builder.get_object('entry_logfile').set_text(output_name)  

    def update_working_folder_chooser (self, folder = None):
        """ Function doc """
        if folder:
            #print('update_working_folder_chooser')
            self.folder_chooser_button.set_folder(folder = folder)
        else:
            
            folder = self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder
            if folder:
                self.folder_chooser_button.set_folder(folder = folder)
            else:
                pass

