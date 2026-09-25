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

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
import gc
import os
from gui.widgets.custom_widgets import FolderChooserButton
from gui.widgets.custom_widgets import SaveTrajectoryBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.widgets.custom_widgets import ReactionCoordinateBox
from gui.widgets.custom_widgets import AdvancedReactionCoordinateBox

from gui.windows.setup.windows_and_dialogs import ExportScriptDialog

from pdynamo.LogFileWriter import detect_reaction_coordinates_from_log

from gui.windows.simulation.PES_scan_window            import compute_sigma_a1_a3
from util.geometric_analysis            import get_dihedral 
from util.geometric_analysis            import get_angle 
from util.geometric_analysis            import get_distance 

#from gui.windows.PES_scan_window            import texto_d1 
#from gui.windows.PES_scan_window            import texto_d2d1 
#from util.periodic_table import atomic_dic 
from pprint import pprint

HOME        = os.environ.get('HOME')



##====================================================================================
class UmbrellaSamplingWindow(Gtk.Window):
    """ Class doc """
    
    def __init__(self, main = None ):
        """ Class initialiser """
        self.main       = main
        self.home       = main.home
        self.p_session  = main.p_session 
        self.vm_session = main.vm_session 
        self.Visible    =  False        
        self.sym_tag    = "umb_sam"
        self.input_types  = {
                            0 : 'From Coordinates (sequential)',
                            1 : 'From Trajectory (parallel)'   ,
                           }
        
        

        self.opt_methods = { 
                            0 : 'ConjugatedGradient',
                            1 : 'SteepestDescent'   ,
                            2 : 'LFBGS'             ,
                            3 : 'QuasiNewton'       ,
                            4 : 'FIRE'              ,
                             }

    
        self.md_methods = {
                           0:"Velocity Verlet Dynamics", 
                           1:"Leap Frog Dynamics",
                           2:"Langevin Dynamics"
                            }
    
        self.rc_liststore1 = Gtk.ListStore(
            str,               # 0: atom1 name
            str,               # 1: atom1 number
            str,               # 2: atom2 name
            str,               # 3: atom2 number
            str,               # 4: weight
            str,               # 5: dist
            )
        self.rc_liststore2 = Gtk.ListStore(
            str,               # 0: atom1 name
            str,               # 1: atom1 number
            str,               # 2: atom2 name
            str,               # 3: atom2 number
            str,               # 4: weight
            str,               # 5: dist
            )
        
        self.rc_liststore1.connect("row-inserted", self.on_row_inserted)
        self.rc_liststore1.connect("row-deleted", self.on_row_deleted)
        self.rc_liststore1.connect("row-changed", self.on_row_changed)
        #self.rc_liststore1.connect("rows-reordered", self.on_rows_reordered)
        
        self.rc_liststore2.connect("row-inserted", self.on_row_inserted)
        self.rc_liststore2.connect("row-deleted", self.on_row_deleted)
        self.rc_liststore2.connect("row-changed", self.on_row_changed)
        #self.rc_liststore2.connect("rows-reordered", self.on_rows_reordered)
    
    
    def open_window (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = self.main.builder #Gtk.Builder()
            
            #self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/geometry_optimization_window.glade'))
            self.builder = Gtk.Builder()                
            self.builder.add_from_file(os.path.join (self.home,'src/gui/windows/simulation/umbrella_sampling_window.glade') )
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('umbrella_samppling_window')
            self.window.set_title('Umbrella Samppling Window')
            self.window.set_keep_above(True)
            self.window.connect("destroy", self.close_window)
            

            '''--------------------------------------------------------------------------------------------'''
            self.input_type_combo = self.builder.get_object('combobox_input_type')
            self.input_type_store = Gtk.ListStore(str)

            for key, input_type in self.input_types.items():
                self.input_type_store.append([input_type])
            
            renderer_text = Gtk.CellRendererText()
            
            self.input_type_combo.set_model(self.input_type_store)
            self.input_type_combo.pack_start(renderer_text, True)
            self.input_type_combo.add_attribute(renderer_text, "text", 0)            
            self.input_type_combo.connect("changed", self.on_combobox_inputtype)
            '''--------------------------------------------------------------------------------------------'''
            
            

            #----------------------------------------------------------------------------------------------
            # - - - - - - - - - - - - - Starting Coordinates ComboBox - - - - - - - - - - - - - - - - -
            self.box_coordinates = self.builder.get_object('box_coordinates')
            self.combobox_starting_coordinates = CoordinatesComboBox() #self.builder.get_object('coordinates_combobox')
            #self.combobox_starting_coordinates.connect("changed", self.on_name_combo_changed)
            self.box_coordinates.pack_start(self.combobox_starting_coordinates, False, False, 0)
            self._starting_coordinates_model_update(init = True)
            #----------------------------------------------------------------------------------------------
            
            
            
            #----------------------------------------------------------------------------------------------
            # - - - - - - - - - - - - - Folder Chooser Button - - - - - - - - - - - - - - - - -
            #----------------------------------------------------------------------------------------------
            self.folder_chooser_button = FolderChooserButton(main =  self.main, sel_type = 'folder', home =  self.home,
                                                              on_folder_selected = self.on_trajectory_folder_selected)
            self.builder.get_object('folder_chooser_box').pack_start(self.folder_chooser_button.btn, True, True, 0)
            system_id      = self.p_session.active_id
            #----------------------------------------------------------------------------------------------
            #----------------------------------------------------------------------------------------------
            self.folder_chooser_button2 = FolderChooserButton(main =  self.main, sel_type = 'folder', home =  self.home)
            self.builder.get_object('folder_chooser_box2').pack_start(self.folder_chooser_button2.btn, True, True, 0)
            #----------------------------------------------------------------------------------------------
            
        
        
            #----------------------------------------------------------------------------------------------
            # - - - - - - - - - - - - - - - - Spin Button Number Of CPUS - - - - - - - - - - - - - - - - -
            #----------------------------------------------------------------------------------------------
            self.spinbutton = self.builder.get_object('ncpus_spinbutton')
            # set the range of values that the spinbutton can take
            self.spinbutton.set_range(1, 1000)
            # set the increment step when the user clicks the up/down buttons
            self.spinbutton.set_increments(1, 10)
            # set the number of decimal places to display
            self.spinbutton.set_digits(0)
            # set the value of the spinbutton
            self.spinbutton.set_value(1)
            #----------------------------------------------------------------------------------------------
            
            
            
            #----------------------------------------------------------------------------------------------
            # - - - - - - - - - - - - - Reaction Coordinates 1 ComboBox - - - - - - - - - - - - - - - - -
            #----------------------------------------------------------------------------------------------
            # "Simple" (fixed-shape: distance/multiple-distance/dihedral) and
            # "advanced" (arbitrary weighted distance list) share the same
            # Alignment placeholder; only one of each pair is visible at a
            # time (see on_advanced_mode_toggled). Advanced is the default,
            # matching this window's own behavior before the simple box was
            # re-enabled (its construction used to be commented out below).
            self.RC_box1 = ReactionCoordinateBox(main = self.main, mode = 1)
            self.RC_box1_adv = AdvancedReactionCoordinateBox(main = self.main, liststore = self.rc_liststore1, mode = 0)
            rc1_box = Gtk.Box(orientation = Gtk.Orientation.VERTICAL)
            rc1_box.pack_start(self.RC_box1, True, True, 0)
            rc1_box.pack_start(self.RC_box1_adv, True, True, 0)
            self.builder.get_object("rc1_aligment").add(rc1_box)

            self.RC_box2 = ReactionCoordinateBox(main = self.main, mode = 1)
            self.RC_box2_adv = AdvancedReactionCoordinateBox(main = self.main, liststore = self.rc_liststore2, mode = 0)
            rc2_box = Gtk.Box(orientation = Gtk.Orientation.VERTICAL)
            rc2_box.pack_start(self.RC_box2, True, True, 0)
            rc2_box.pack_start(self.RC_box2_adv, True, True, 0)
            self.builder.get_object("rc2_aligment").add(rc2_box)
        
        
            #----------------------------------------------------------------------------------------------
            # - - - - - - - - - - - - - -  Geometry Optimization ComboBox - - - - - - - - - - - - - - - - -
            #----------------------------------------------------------------------------------------------
            self.opt_methods_combo = self.builder.get_object('combobox_geo_opt')
            self.opt_method_store = Gtk.ListStore(str)

            for key, method in self.opt_methods.items():
                self.opt_method_store.append([method])
            
            self.opt_methods_combo.set_model(self.opt_method_store)
            renderer_text = Gtk.CellRendererText()
            self.opt_methods_combo.pack_start(renderer_text, True)
            self.opt_methods_combo.add_attribute(renderer_text, "text", 0)
            #----------------------------------------------------------------------------------------------
            
            

            #----------------------------------------------------------------------------------------------
            # - - - - - - - - - - - - - - - - Molecular Dynamics ComboBox - - - - - - - - - - - - - - - - -
            #----------------------------------------------------------------------------------------------
            self.md_integrators_combobox = self.builder.get_object('md_integrator_combobox')
            self.md_integrators_store    = Gtk.ListStore(str)

            for key, integrators in self.md_methods.items():
                self.md_integrators_store.append([integrators])
            
            self.md_integrators_combobox.set_model(self.md_integrators_store)
            renderer_text = Gtk.CellRendererText()
            self.md_integrators_combobox.pack_start(renderer_text, True)
            self.md_integrators_combobox.add_attribute(renderer_text, "text", 0)
            self.md_integrators_combobox.connect("changed", self.on_md_integrator_combobox)
            #----------------------------------------------------------------------------------------------
            
            
            
            #---------------------------------------------------------------------------------------------
            self.builder.get_object('checkbox_reaction_coordinate2').connect("clicked", self.on_checkbox_reaction_coordinate2)
            self.builder.get_object('checkbox_geometry_optimization').connect("clicked", self.on_checkbox_geometry_optimization)
            self.builder.get_object('checkbtn_advanced_mode').connect("toggled", self.on_advanced_mode_toggled)

            self.RC_box2.set_sensitive(False)
            self.RC_box2_adv.set_sensitive(False)
            self.builder.get_object('frame_geometry_optimization').set_sensitive(False)

            self.builder.get_object('button_run').connect("clicked", self.run)
            self.builder.get_object('button_cancel').connect('clicked', self.close_window)
            self.builder.get_object('button_export').connect('clicked', self.on_btn_export)
            #---------------------------------------------------------------------------------------------



            #self._starting_coordinates_model_update()
            self.build_advanced_treeviews()
            self.window.show_all()

            self.input_type_combo.set_active(0)

            self.opt_methods_combo.set_active(0)
            self.md_integrators_combobox.set_active(0)
            self.update_working_folder_chooser ( )

            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
                self.builder.get_object('entry_traj_name').set_text(output_name)
            else:
                pass

            self.RC_box1.set_rc_mode(rc_mode=0)
            self.RC_box2.set_rc_mode(rc_mode=0)
            self.RC_box1_adv.set_rc_mode(rc_mode=0)
            self.RC_box2_adv.set_rc_mode(rc_mode=0)
            self.RC_box1.set_rc_type(0)
            self.RC_box2.set_rc_type(0)

            # Advanced mode is the default -- matches this window's own
            # behavior before the simple box was re-enabled.
            self.builder.get_object('checkbtn_advanced_mode').set_active(True)
            self.on_advanced_mode_toggled(None)
            self.Visible  = True

        else:
            self.window.present()
            
    def close_window (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False

    def build_advanced_treeviews (self):
        """ Function doc """

        self.treeview1 = Gtk.TreeView(model = self.rc_liststore1)
        self.treeview2 = Gtk.TreeView(model = self.rc_liststore2)

        columns = {"atm1"     : 0,
                   'idx1'    : 1,
                   "atm2"     : 2,
                   "idx2"    : 3,
                   "weight"    : 4,
                   "dist"    : 5,
                   }
        #treeview 2
        for title in columns.keys():
            renderer = Gtk.CellRendererText()
            renderer.set_property("editable", True)   # permite editar
            i = columns[title]
            renderer.connect("edited", self.on_cell_edited1, i)
            column = Gtk.TreeViewColumn(title, renderer, text=i)
            self.treeview1.append_column(column)

        #treeview 2
        for title in columns.keys():
            renderer = Gtk.CellRendererText()
            renderer.set_property("editable", True)   # permite editar
            i = columns[title]
            renderer.connect("edited", self.on_cell_edited2, i)
            column = Gtk.TreeViewColumn(title, renderer, text=i)
            self.treeview2.append_column(column)
            #column.set_fixed_width(20)
        # Enable context menu on right-click
        #self.treeview1.connect("button-press-event", self.on_button_press_event)
        #self.treeview2.connect("button-press-event", self.on_button_press_event)

        self.treeview1.connect("key-press-event", self.on_key_press)
        self.treeview2.connect("key-press-event", self.on_key_press)

        #self.scrolledbox
        self.RC_box1_adv.scrolledbox.add(self.treeview1)
        self.RC_box2_adv.scrolledbox.add(self.treeview2)
        self.RC_box1_adv.treeview = self.treeview1
        self.RC_box2_adv.treeview = self.treeview2

    def on_advanced_mode_toggled (self, widget) -> None:
        """Switches both RC boxes between "simple" (fixed-shape) and
        "advanced" (weighted distance list) mode."""
        advanced = self.builder.get_object('checkbtn_advanced_mode').get_active()

        self.RC_box1.set_visible(not advanced)
        self.RC_box1_adv.set_visible(advanced)
        self.RC_box2.set_visible(not advanced)
        self.RC_box2_adv.set_visible(advanced)

    def on_row_inserted(self, model, path, iter):
        if model is  self.rc_liststore1:
            #print("Linha inserida em:", path)
            self.RC_box1_adv.refresh_dmininum()
        else:
            self.RC_box2_adv.refresh_dmininum()

    def on_row_deleted(self, model, path):
        dprint("Row removed at:", path)
        if model is  self.rc_liststore1:
            #print("Linha inserida em:", path)
            self.RC_box1_adv.refresh_dmininum()
        else:
            self.RC_box2_adv.refresh_dmininum()
    
    def on_row_changed(self, model, path, iter):
        pass
        #print("Linha alterada em:", path)
        
    def on_cell_edited1 (self, widget, path, new_text, column_index):
        """
        widget        -> o CellRendererText
        path          -> índice da linha (ex: "0", "1", ...)
        new_text      -> texto digitado pelo usuário
        column_index  -> índice da coluna (0 ou 1) - enviado como user_data
        """
        self.rc_liststore1[path][column_index] = new_text
    
    def on_cell_edited2 (self, widget, path, new_text, column_index):
        """ Function doc """
        """
        path: string com o índice da linha ("0", "1", ...)
        new_text: novo valor digitado
        """
        self.rc_liststore2[path][column_index] = new_text

    def on_key_press(self, widget, event):
        # checking the Delete key
        if event.keyval == Gdk.KEY_Delete:
            selection = widget.get_selection()
            model, treeiter = selection.get_selected()

            if treeiter is not None:
                model.remove(treeiter)
            return True  # evita propagação do evento

        return False

    def _set_rc_input_fields_sensitive (self, sensitive) -> None:
        """Enables/disables the per-window step-size/nsteps/dmin fields on
        ALL FOUR RC boxes (simple + advanced, RC1 + RC2) -- both
        ReactionCoordinateBox and AdvancedReactionCoordinateBox expose the
        same entry_step_size/entry_nsteps/entry_dmin_coord/label_step_size
        attribute names, so this applies uniformly regardless of which
        pair is currently visible (see on_advanced_mode_toggled)."""
        for rc_box in (self.RC_box1, self.RC_box1_adv, self.RC_box2, self.RC_box2_adv):
            rc_box.entry_step_size .set_sensitive(sensitive)
            rc_box.entry_nsteps    .set_sensitive(sensitive)
            rc_box.entry_dmin_coord.set_sensitive(sensitive)
            rc_box.label_step_size .set_sensitive(sensitive)

    def on_combobox_inputtype (self, combobox):
        """ Function doc """
        _type =  combobox.get_active()

        if _type == 0:
            self.builder.get_object('label_input_trajectory').hide()
            self.builder.get_object('folder_chooser_box').hide()
            self.builder.get_object('label_number_of_cpus').hide()
            self.builder.get_object('ncpus_spinbutton').hide()

            self.builder.get_object('label_starting_coordinates').show()
            self.combobox_starting_coordinates.show()

            self._set_rc_input_fields_sensitive(True)

        if _type == 1:
            self.builder.get_object('label_input_trajectory').show()
            self.builder.get_object('folder_chooser_box').show()
            self.builder.get_object('label_number_of_cpus').show()
            self.builder.get_object('ncpus_spinbutton').show()

            self.builder.get_object('label_starting_coordinates').hide()
            self.combobox_starting_coordinates.hide()

            self._set_rc_input_fields_sensitive(False)



    def on_md_integrator_combobox (self, widget = None):
        """ Function doc 
        
            0  :  "Velocity Verlet Dynamics", 
            1  :  "Leap Frog Dynamics",
            2  :  "Langevin Dynamics"
        
        """        
        if widget  == self.builder.get_object('md_integrator_combobox'):
            
            #velocity verlet molecular dynamics
            if 0 == widget.get_active():
                self.builder.get_object('check_pressure_control').hide()
                self.builder.get_object('entry_pressure').hide()
                self.builder.get_object('entry_pressure').hide()
                
                
                self.builder.get_object('entry_temp_coupling').hide()
                self.builder.get_object('label_temp_coupling').hide()
                
                self.builder.get_object('label_pressure_coupling').hide()
                self.builder.get_object('entry_pressure_coupling').hide()
                
                self.builder.get_object('label_collision_freq').hide()
                self.builder.get_object('entry_collision_frequency').hide()
                
                #self.temp_scale_options_combo.set_sensitive(True)
                
            #Leap Frog molecular dynamics
            elif 1 == widget.get_active():
                self.builder.get_object('check_pressure_control').show()
                self.builder.get_object('entry_pressure').show()
                
                self.builder.get_object('label_pressure_coupling').show()
                self.builder.get_object('entry_pressure_coupling').show()
                
                
                self.builder.get_object('entry_temp_coupling').show()
                self.builder.get_object('label_temp_coupling').show()
                
                
                self.builder.get_object('label_collision_freq').hide()
                self.builder.get_object('entry_collision_frequency').hide()
                
                #self.temp_scale_options_combo.set_active(0)
                #self.temp_scale_options_combo.set_sensitive(False)
                
            #Langevin molecular dynamics
            elif 2 == widget.get_active():
                self.builder.get_object('check_pressure_control').hide()
                self.builder.get_object('entry_pressure').hide()
                
                self.builder.get_object('entry_temp_coupling').hide()
                self.builder.get_object('label_temp_coupling').hide()
                
                self.builder.get_object('label_pressure_coupling').hide()
                self.builder.get_object('entry_pressure_coupling').hide()
                
                self.builder.get_object('label_collision_freq').show()
                self.builder.get_object('entry_collision_frequency').show()
            
                #self.temp_scale_options_combo.set_active(0)
                #self.temp_scale_options_combo.set_sensitive(False)
            else:
                pass

        
        
        
        elif widget == self.builder.get_object('combobox_temperature_scale_options'):
            temp_scale_id = widget.get_active()
            
            if temp_scale_id == 0:
                self.builder.get_object('entry_temp_end').set_sensitive(False)
                self.builder.get_object('label_temp_end').set_sensitive(False)
            else:
                self.builder.get_object('entry_temp_end').set_sensitive(True)
                self.builder.get_object('label_temp_end').set_sensitive(True)

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

    
    def run_dialog (self):
        """ Function doc """
        dialog = Gtk.MessageDialog(
                    transient_for=self.main.window,
                    message_type=Gtk.MessageType.ERROR,
                    buttons=Gtk.ButtonsType.OK,
                    text="Folder not found",
                    )
        dialog.format_secondary_text(
            "The folder you have selected does not appear to be valid. Please select a different folder or create a new one."
        )
        dialog.run()
        dialog.destroy()

    
    def update (self, parameters = None):
        """ Function doc """
        
        if self.Visible:
            self._starting_coordinates_model_update()
            self.update_working_folder_chooser()
        
            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
                self.builder.get_object('entry_traj_name').set_text(output_name)
            else:
                pass


    def update_working_folder_chooser (self, folder = None):
        """ Function doc """
        if folder:
            #print('update_working_folder_chooser')
            self.folder_chooser_button.set_folder(folder = folder)
        else:
            
            try:
                folder = self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder
                if folder:
                    self.folder_chooser_button2.set_folder(folder = folder)
                else:
                    pass
            except:
                self.folder_chooser_button2.set_folder(folder = self.home )
   

    def on_checkbox_reaction_coordinate2 (self, widget):
        """ Function doc """
        active = widget.get_active()
        self.RC_box2.set_sensitive(active)
        self.RC_box2_adv.set_sensitive(active)

    def on_trajectory_folder_selected (self, folder):
        """ Auto-detects the reaction coordinate(s) used to produce the
        trajectory in a freshly-picked "From Trajectory" input folder, by
        parsing its own 'output.log' (written by whichever job created it
        -- see LogFileWriter.detect_reaction_coordinates_from_log for the
        recognized formats: RelaxedSurfaceScan/AdvancedRelaxedSurfaceScan's
        PES scans, or EnergyRefinement), and pre-fills RC1/RC2 with them so
        the user doesn't have to redefine the same atoms that already
        define this trajectory's reaction coordinate.

        Silently does nothing if output.log doesn't exist or isn't a
        recognizable EasyHybrid scan log -- this is a convenience, not a
        requirement; the user can still enter the RC by hand exactly as
        before.
        """
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
            checkbox2 = self.builder.get_object('checkbox_reaction_coordinate2')
            checkbox2.set_active(True)
            self.on_checkbox_reaction_coordinate2(checkbox2)
            if advanced:
                self._apply_detected_advanced_rc(self.rc_liststore2, rc2)
            else:
                self._apply_detected_simple_rc(self.RC_box2, rc2)

    def _apply_detected_simple_rc (self, rc_box, data):
        """ Fills ONLY a ReactionCoordinateBox's coordinate-type combobox
        and atom index/name entries from a detect_reaction_coordinates_from_log()
        result -- deliberately does NOT touch force_constant/nsteps/
        dincre/dminimum (unlike ReactionCoordinateBox.set_rc_data(), which
        would reset those to 0/0.0 since the detected dict never carries
        them -- they belong to the run that's being newly configured, not
        to the trajectory that's being imported as its input). """
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
        detect_reaction_coordinates_from_log() result -- same "clear and
        repopulate the liststore" approach PotentialEnergyScanWindow's own
        _restore_advanced_rc_data uses, minus the force_constant/nsteps/
        dincre/dminimum fields for the same reason as
        _apply_detected_simple_rc above. """
        liststore.clear()
        for row in data['RC']:
            liststore.append(list(row))

    def on_checkbox_geometry_optimization (self, widget):
        """ Function doc """
        if widget.get_active():
            self.builder.get_object('frame_geometry_optimization').set_sensitive(True)
        else:
            self.builder.get_object('frame_geometry_optimization').set_sensitive(False)

   
    def get_input_setup_box_info (self, parameters):
        """ Function doc """
        parameters['input_type'] = self.builder.get_object('combobox_input_type').get_active()
        
        if parameters['input_type'] == 0:
            vobject_id = self.combobox_starting_coordinates.get_vobject_id()
            vobject = self.main.vm_session.vm_objects_dic[vobject_id]
            self.main.p_session.set_psystem_coordinates_from_vobject(vobject)
            parameters['vobj'] = vobject.name
            parameters['source_folder'] = None
            
        elif parameters['input_type'] == 1:
            parameters['vobj'] = None
            parameters['source_folder'] = self.folder_chooser_button.get_folder ()
            
        
        return parameters
    
        
    def get_parameters (self, ):
        
        """ Function doc """
        
        '''this combobox has the reference to the starting coordinates of a simulation'''
        parameters = {"simulation_type" : "Umbrella_Sampling"}

        parameters = self.get_input_setup_box_info ( parameters)
        parameters["folder"] = self.folder_chooser_button2.get_folder()        


        #vobject_id = self.combobox_starting_coordinates.get_vobject_id()
        #vobject = self.main.vm_session.vm_objects_dic[vobject_id]
        #self.main.p_session.set_psystem_coordinates_from_vobject(vobject)

        '''
        
        - If a trajectory folder is given then parallelization is allowed 
            - NmaxThreads >= 1
        
        - If a vobject is provided and a second reaction coordinate is active 
          (2D), then parallelization is allowed
            - NmaxThreads >= 1
        
        - Else: sequetial
        
        '''
        advanced = self.builder.get_object('checkbtn_advanced_mode').get_active()
        rc_box1  = self.RC_box1_adv if advanced else self.RC_box1
        rc_box2  = self.RC_box2_adv if advanced else self.RC_box2

        parameters["RC1"] = rc_box1.get_rc_data()
        if self.builder.get_object('checkbox_reaction_coordinate2').get_active():
            parameters["RC2"] = rc_box2.get_rc_data()
            parameters["NmaxThreads"] =  int(self.builder.get_object('ncpus_spinbutton').get_value())
        else:
            parameters["RC2"] = None

            if parameters['input_type'] == 1:
                parameters["NmaxThreads"] = int(self.builder.get_object('ncpus_spinbutton').get_value())
            else:
                parameters["NmaxThreads"] = 1

        
        
        #----------------------------------------------------------------------
        #                            GEO OPT
        #----------------------------------------------------------------------
        if self.builder.get_object('checkbox_geometry_optimization').get_active():
            parameters['OPT_parm'] = self._define_OPT_parameters ()
        else:
            parameters['OPT_parm'] = None
        #----------------------------------------------------------------------
        
        
        #----------------------------------------------------------------------
        #                          M. DYNAMICS
        #----------------------------------------------------------------------        
        parameters['MD_parm'] = self._define_MD_parameters ( )
        #----------------------------------------------------------------------
        
        
        parameters['traj_folder_name'] = self.builder.get_object('entry_traj_name').get_text()
        
        return parameters
    
        
    def run (self, button):
        """ Function doc """
        
        parameters = self.get_parameters()
        
        #pprint(parameters)

        isExist = os.path.exists(parameters['folder'])
        if isExist:
            pass
        else:
            self.run_dialog()
            return None
        
        
        self.main.p_session.run_simulation( parameters = parameters )
        self.window.destroy()
        self.Visible    =  False


    def on_btn_export (self, widget):
        """ Function doc """

        #header += '\nparameters = '
        parameters = self.get_parameters()
        try:
            parameters.pop('vobject_name')
        except:
            pass
        #pprint.pprint(parameters)

        export_dialog =  ExportScriptDialog(self.main, parameters = parameters)
        

    def _define_OPT_parameters (self):
        parameters = {}
        opt_id       = self.builder.get_object('combobox_geo_opt').get_active()
        parameters['optimizer']        =  self.opt_methods[opt_id]
        parameters["maximumIterations"]    =  int  (self.builder.get_object('entry_max_int') .get_text())
        parameters["rmsGradientTolerance"]      =  float(self.builder.get_object('entry_rmsd_tol').get_text())
        return parameters
    
    
    def _define_MD_parameters (self):
        """ Function doc """
        integrator_id       = self.builder.get_object('md_integrator_combobox').get_active()
        number_of_steps_eq  = int(self.builder.get_object('entry_number_of_steps_eq').get_text())
        number_of_steps_dc  = int(self.builder.get_object('entry_number_of_steps_dc').get_text())
        
        
        temp_start          = float(self.builder.get_object('entry_temp_start').get_text())
        temp_scale_factor   = int(self.builder.get_object('entry_temp_scale_factor').get_text())
        time_step           = float(self.builder.get_object('entry_time_step').get_text())
        log_frequence       = int(self.builder.get_object('entry_log_frequency').get_text())
        random_seed         = int(self.builder.get_object('entry_random_seed').get_text())
        
        collision_frequency = None
        temp_coupling       = None #float(self.builder.get_object('entry_temp_coupling').get_text())
        
        if self.builder.get_object('check_pressure_control').get_active():
            pressure_control = True
        else:
            pressure_control = False
        
        #-------------------------------------------------------------------------------------------
        if integrator_id == 1: #LeapFrog
            temp_coupling     = float(self.builder.get_object('entry_temp_coupling').get_text())
            pressure          = float(self.builder.get_object('entry_pressure').get_text())
            pressure_coupling = float(self.builder.get_object('entry_pressure_coupling').get_text())
            temp_control      = True
        else:
            pressure          = None
            pressure_coupling = None
            temp_control      = True
        #-------------------------------------------------------------------------------------------
        
        if integrator_id == 2: #Langevin
            temp_coupling       = float(self.builder.get_object('entry_temp_coupling').get_text())
            collision_frequency = float(self.builder.get_object('entry_collision_frequency').get_text())
        else:
            collision_frequency = None
        
        MD_method = {
                     0 : "Verlet"   ,
                     1 : "LeapFrog" ,
                     2 : "Langevin" ,
                     }
        
        

        parameters = {
                    "simulation_type"           : "Umbrella_Sampling",
                    #"folder"                    : HOME                    ,
                    'integrator'                : MD_method[integrator_id], # verlet / leapfrog /langevin
                    'logFrequency'              : log_frequence           ,
                    'seed'                      : random_seed             ,
                    'normal_deviate_generator'  : None                    ,
                    'steps_eq'                  : number_of_steps_eq      ,
                    'steps_dc'                  : number_of_steps_dc      ,
                    'timeStep'                  : time_step               ,
                    #'trajectories'              : None                    ,
                    'trajectory_frequency'          : int(self.builder.get_object('entry_traj_frequency').get_text()), 
                    'trajectory_frequency_dc_ptRes' : int(self.builder.get_object('entry_traj_frequency_dc1').get_text()), 
                    'trajectory_frequency_dc_ptGeo' : int(self.builder.get_object('entry_traj_frequency_dc2').get_text()), 
                    
                    #VelocityVerletDynamics
                    'temperatureScaleFrequency' : temp_scale_factor                ,
                    'temperatureStart'          : temp_start                       ,
                    
                    #LeapFrogDynamics
                    'pressure'                  : pressure           ,  #  LeapFrogDynamics  1.0 default
                    'pressureControl'           : pressure_control   ,  #  LeapFrogDynamics  True / False
                    'pressureCoupling'          : pressure_coupling  ,  #  LeapFrogDynamics  2000
                                                                     
                    'temperature'               : temp_start         ,               
                    'temperatureControl'        : temp_control   ,  # True / False LeapFrogDynamics / LangevinDynamics
                    'temperatureCoupling'       : temp_coupling      ,  #  LeapFrogDynamics / LangevinDynamics
                    
                    #LangevinDynamics
                    'collisionFrequency'        : collision_frequency,  # Langevin 25
                    }
        return parameters

    def restore_the_parameters_to_the_window(self, parameters):
        """Update the GUI widgets with values from the parameters dictionary."""
        MD_ids = {
             "Verlet"   : 0 ,
             "LeapFrog" : 1 ,
             "Langevin" : 2 ,
             }
        
        
        pprint(parameters)
        if not parameters:
            return
        
        #--------------------------------------------------------------
        self.combobox_starting_coordinates.set_active(parameters['cb1_active'])
        #--------------------------------------------------------------    
        
        # ----------------------------
        # Input type
        # ----------------------------
        input_type = parameters.get('input_type', 0)
        self.input_type_combo.set_active(input_type)

        # Set starting coordinates or trajectory folder
        if input_type == 0 and parameters.get('vobj') is not None:
            # Set starting coordinates
            vobj_name = parameters['vobj']
            vobject_ids = self.main.vm_session.vm_objects_dic
            for vobj_id, vobj in vobject_ids.items():
                if vobj.name == vobj_name:
                    self.combobox_starting_coordinates.set_active(vobj_id)
                    break
            self.combobox_starting_coordinates.show()
            self.folder_chooser_button.btn.hide()
            self.spinbutton.hide()
        elif input_type == 1 and parameters.get('source_folder') is not None:
            # Set trajectory folder
            self.folder_chooser_button.set_folder(parameters['source_folder'])
            self.folder_chooser_button.btn.show()
            self.spinbutton.show()
            self.combobox_starting_coordinates.hide()

        # ----------------------------
        # Reaction Coordinates
        # ----------------------------
        # 'RC' present -> this job was set up in advanced (weighted
        # distance list) mode; 'rc_type' present -> simple (fixed-shape,
        # possibly dihedral) mode. Mirrors the same check
        # UmbrellaSampling.run() itself uses (p_methods/umbrella_sampling.py).
        rc1_advanced = bool(parameters.get('RC1')) and 'RC' in parameters['RC1']
        self.builder.get_object('checkbtn_advanced_mode').set_active(rc1_advanced)
        self.on_advanced_mode_toggled(None)

        if 'RC1' in parameters and parameters['RC1'] is not None:
            if rc1_advanced:
                self._restore_advanced_rc_data(self.RC_box1_adv, self.rc_liststore1, parameters['RC1'])
            else:
                self.RC_box1.set_rc_data(parameters['RC1'])

        if 'RC2' in parameters and parameters['RC2'] is not None:
            if rc1_advanced:
                self._restore_advanced_rc_data(self.RC_box2_adv, self.rc_liststore2, parameters['RC2'])
            else:
                self.RC_box2.set_rc_data(parameters['RC2'])
            self.builder.get_object('checkbox_reaction_coordinate2').set_active(True)
            self.RC_box2.set_sensitive(True)
            self.RC_box2_adv.set_sensitive(True)
        else:
            self.builder.get_object('checkbox_reaction_coordinate2').set_active(False)
            self.RC_box2.set_sensitive(False)
            self.RC_box2_adv.set_sensitive(False)

    def _restore_advanced_rc_data(self, rc_box, liststore, data) -> None:
        """Repopulates an AdvancedReactionCoordinateBox's weighted-distance
        table and summary fields from a previously saved RC1/RC2 dict (see
        AdvancedReactionCoordinateBox.get_rc_data for the shape). NOT the
        same as (the broken, dead) AdvancedReactionCoordinateBox.set_rc_data
        -- that method references combobox_reaction_coord1, which only
        exists on the SIMPLE ReactionCoordinateBox and would raise
        AttributeError if it were ever actually called."""
        liststore.clear()
        for row in data.get("RC", []):
            liststore.append(list(row))
        rc_box.builder.get_object("entry_dmin_coord1").set_text(str(data.get("dminimum", 0.0)))
        rc_box.builder.get_object("entry_nsteps1").set_text(str(data.get("nsteps", 0)))
        rc_box.builder.get_object("entry_FORCE_coord1").set_text(str(data.get("force_constant", 0.0)))
        step_size_entry = rc_box.builder.get_object("entry_step_size1")
        if step_size_entry:
            step_size_entry.set_text(str(data.get("dincre", 0.0)))

        # ----------------------------
        # Geometry Optimization
        # ----------------------------
        if parameters.get('OPT_parm') is not None:
            self.builder.get_object('checkbox_geometry_optimization').set_active(True)
            self.builder.get_object('frame_geometry_optimization').set_sensitive(True)
            opt_index = list(self.opt_methods.values()).index(parameters['OPT_parm']['optimizer'])
            self.opt_methods_combo.set_active(opt_index)
            self.builder.get_object('entry_max_int').set_text(str(parameters['OPT_parm']['maximumIterations']))
            self.builder.get_object('entry_rmsd_tol').set_text(str(parameters['OPT_parm']['rmsGradientTolerance']))
        else:
            self.builder.get_object('checkbox_geometry_optimization').set_active(False)
            self.builder.get_object('frame_geometry_optimization').set_sensitive(False)

        # ----------------------------
        # Molecular Dynamics
        # ----------------------------
        if parameters.get('MD_parm') is not None:
            md_parm = parameters['MD_parm']
            #integrator = md_parm['integrator']
            integrator_index = MD_ids[md_parm['integrator']]
            
            #integrator_index = list(md_parm['integrator'] for md_parm in md_parm.items() if md_parm).index(md_parm['integrator']) \
            #    if md_parm['integrator'] in ["Verlet", "LeapFrog", "Langevin"] else 0
            self.md_integrators_combobox.set_active(integrator_index)

            self.builder.get_object('entry_number_of_steps_eq').set_text(str(md_parm.get('steps_eq', 0)))
            self.builder.get_object('entry_number_of_steps_dc').set_text(str(md_parm.get('steps_dc', 0)))
            self.builder.get_object('entry_time_step').set_text(str(md_parm.get('timeStep', 0.0)))
            self.builder.get_object('entry_log_frequency').set_text(str(md_parm.get('logFrequency', 0)))
            self.builder.get_object('entry_random_seed').set_text(str(md_parm.get('seed', 0)))
            self.builder.get_object('entry_temp_start').set_text(str(md_parm.get('temperatureStart', 300)))
            self.builder.get_object('entry_temp_scale_factor').set_text(str(md_parm.get('temperatureScaleFrequency', 1)))
            self.builder.get_object('entry_collision_frequency').set_text(str(md_parm.get('collisionFrequency', 0.0)))
            self.builder.get_object('entry_pressure').set_text(str(md_parm.get('pressure', 1.0)))
            self.builder.get_object('entry_pressure_coupling').set_text(str(md_parm.get('pressureCoupling', 1.0)))
            self.builder.get_object('entry_temp_coupling').set_text(str(md_parm.get('temperatureCoupling', 1.0)))

        # ----------------------------
        # Trajectory folder name
        # ----------------------------
        if 'traj_folder_name' in parameters:
            self.builder.get_object('entry_traj_name').set_text(parameters['traj_folder_name'])

        # ----------------------------
        # Number of threads (CPUs)
        # ----------------------------
        if 'NmaxThreads' in parameters:
            self.spinbutton.set_value(parameters['NmaxThreads'])


class UmbrellaSamplingWindow_old(Gtk.Window):
    """ Class doc """
    
    def __init__(self, main = None ):
        """ Class initialiser """
        self.main       = main
        self.home       = main.home
        self.p_session  = main.p_session 
        self.vm_session = main.vm_session 
        self.Visible    =  False        
        self.sym_tag    = "umb_sam"
        self.input_types  = {
                            0 : 'From Coordinates (sequential)',
                            1 : 'From Trajectory (parallel)'   ,
                           }
        
        

        self.opt_methods = { 
                            0 : 'ConjugatedGradient',
                            1 : 'SteepestDescent'   ,
                            2 : 'LFBGS'             ,
                            3 : 'QuasiNewton'       ,
                            4 : 'FIRE'              ,
                             }

    
        self.md_methods = {
                           0:"Velocity Verlet Dynamics", 
                           1:"Leap Frog Dynamics",
                           2:"Langevin Dynamics"
                            }
    
    
    def open_window (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = self.main.builder #Gtk.Builder()
            
            #self.builder.add_from_file(os.path.join(VISMOL_HOME,'easyhybrid/gui/geometry_optimization_window.glade'))
            self.builder = Gtk.Builder()                
            self.builder.add_from_file(os.path.join (self.home,'src/gui/windows/simulation/umbrella_sampling_window.glade') )
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('umbrella_samppling_window')
            self.window.set_title('Umbrella Samppling Window')
            self.window.set_keep_above(True)
            self.window.connect("destroy", self.close_window)
            

            '''--------------------------------------------------------------------------------------------'''
            self.input_type_combo = self.builder.get_object('combobox_input_type')
            self.input_type_store = Gtk.ListStore(str)

            for key, input_type in self.input_types.items():
                self.input_type_store.append([input_type])
            
            renderer_text = Gtk.CellRendererText()
            
            self.input_type_combo.set_model(self.input_type_store)
            self.input_type_combo.pack_start(renderer_text, True)
            self.input_type_combo.add_attribute(renderer_text, "text", 0)            
            self.input_type_combo.connect("changed", self.on_combobox_inputtype)
            '''--------------------------------------------------------------------------------------------'''
            
            

            #----------------------------------------------------------------------------------------------
            # - - - - - - - - - - - - - Starting Coordinates ComboBox - - - - - - - - - - - - - - - - -
            self.box_coordinates = self.builder.get_object('box_coordinates')
            self.combobox_starting_coordinates = CoordinatesComboBox() #self.builder.get_object('coordinates_combobox')
            #self.combobox_starting_coordinates.connect("changed", self.on_name_combo_changed)
            self.box_coordinates.pack_start(self.combobox_starting_coordinates, False, False, 0)
            self._starting_coordinates_model_update(init = True)
            #----------------------------------------------------------------------------------------------
            
            
            
            #----------------------------------------------------------------------------------------------
            # - - - - - - - - - - - - - Folder Chooser Button - - - - - - - - - - - - - - - - -
            #----------------------------------------------------------------------------------------------
            self.folder_chooser_button = FolderChooserButton(main =  self.main, sel_type = 'folder', home =  self.home)
            self.builder.get_object('folder_chooser_box').pack_start(self.folder_chooser_button.btn, True, True, 0)
            system_id      = self.p_session.active_id
            #----------------------------------------------------------------------------------------------
            #----------------------------------------------------------------------------------------------
            self.folder_chooser_button2 = FolderChooserButton(main =  self.main, sel_type = 'folder', home =  self.home)
            self.builder.get_object('folder_chooser_box2').pack_start(self.folder_chooser_button2.btn, True, True, 0)
            #----------------------------------------------------------------------------------------------
            
        
        
            #----------------------------------------------------------------------------------------------
            # - - - - - - - - - - - - - - - - Spin Button Number Of CPUS - - - - - - - - - - - - - - - - -
            #----------------------------------------------------------------------------------------------
            self.spinbutton = self.builder.get_object('ncpus_spinbutton')
            # set the range of values that the spinbutton can take
            self.spinbutton.set_range(1, 1000)
            # set the increment step when the user clicks the up/down buttons
            self.spinbutton.set_increments(1, 10)
            # set the number of decimal places to display
            self.spinbutton.set_digits(0)
            # set the value of the spinbutton
            self.spinbutton.set_value(1)
            #----------------------------------------------------------------------------------------------
            
            
            
            #----------------------------------------------------------------------------------------------
            # - - - - - - - - - - - - - Reaction Coordinates 1 ComboBox - - - - - - - - - - - - - - - - -
            #----------------------------------------------------------------------------------------------
            self.RC_box1 = ReactionCoordinateBox(main = self.main, mode = 1)
            self.builder.get_object('rc1_aligment').add(self.RC_box1)
            self.RC_box2 = ReactionCoordinateBox(main = self.main, mode = 1)
            self.builder.get_object('rc2_aligment').add(self.RC_box2)
            #----------------------------------------------------------------------------------------------

        
        
            #----------------------------------------------------------------------------------------------
            # - - - - - - - - - - - - - -  Geometry Optimization ComboBox - - - - - - - - - - - - - - - - -
            #----------------------------------------------------------------------------------------------
            self.opt_methods_combo = self.builder.get_object('combobox_geo_opt')
            self.opt_method_store = Gtk.ListStore(str)

            for key, method in self.opt_methods.items():
                self.opt_method_store.append([method])
            
            self.opt_methods_combo.set_model(self.opt_method_store)
            renderer_text = Gtk.CellRendererText()
            self.opt_methods_combo.pack_start(renderer_text, True)
            self.opt_methods_combo.add_attribute(renderer_text, "text", 0)
            #----------------------------------------------------------------------------------------------
            
            

            #----------------------------------------------------------------------------------------------
            # - - - - - - - - - - - - - - - - Molecular Dynamics ComboBox - - - - - - - - - - - - - - - - -
            #----------------------------------------------------------------------------------------------
            self.md_integrators_combobox = self.builder.get_object('md_integrator_combobox')
            self.md_integrators_store    = Gtk.ListStore(str)

            for key, integrators in self.md_methods.items():
                self.md_integrators_store.append([integrators])
            
            self.md_integrators_combobox.set_model(self.md_integrators_store)
            renderer_text = Gtk.CellRendererText()
            self.md_integrators_combobox.pack_start(renderer_text, True)
            self.md_integrators_combobox.add_attribute(renderer_text, "text", 0)
            self.md_integrators_combobox.connect("changed", self.on_md_integrator_combobox)
            #----------------------------------------------------------------------------------------------
            
            
            
            #---------------------------------------------------------------------------------------------
            self.builder.get_object('checkbox_reaction_coordinate2').connect("clicked", self.on_checkbox_reaction_coordinate2)
            self.builder.get_object('checkbox_geometry_optimization').connect("clicked", self.on_checkbox_geometry_optimization)
            
            self.RC_box2.set_sensitive(False)
            self.builder.get_object('frame_geometry_optimization').set_sensitive(False)
            
            self.builder.get_object('button_run').connect("clicked", self.run)
            self.builder.get_object('button_cancel').connect('clicked', self.close_window)
            self.builder.get_object('button_export').connect('clicked', self.on_btn_export)
            #---------------------------------------------------------------------------------------------
            
            
            
            #self._starting_coordinates_model_update()
            self.window.show_all()
            
            self.input_type_combo.set_active(0)

            self.opt_methods_combo.set_active(0)
            self.md_integrators_combobox.set_active(0)
            self.update_working_folder_chooser ( )
            
            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
                self.builder.get_object('entry_traj_name').set_text(output_name)
            else:
                pass
            
            self.RC_box1.set_rc_type(0)
            self.RC_box2.set_rc_type(0)
            self.Visible  = True   

        else:
            self.window.present()
            

    def close_window (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False


    def on_combobox_inputtype (self, combobox):
        """ Function doc """
        _type =  combobox.get_active()  
        
        if _type == 0:
            self.builder.get_object('label_input_trajectory').hide()
            self.builder.get_object('folder_chooser_box').hide()
            self.builder.get_object('label_number_of_cpus').hide()
            self.builder.get_object('ncpus_spinbutton').hide()
            
            self.builder.get_object('label_starting_coordinates').show()
            self.combobox_starting_coordinates.show()
        
            self.RC_box1.entry_step_size .set_sensitive(True)
            self.RC_box1.entry_nsteps    .set_sensitive(True)
            self.RC_box1.entry_dmin_coord.set_sensitive(True)
            self.RC_box2.entry_step_size .set_sensitive(True)
            self.RC_box2.entry_nsteps    .set_sensitive(True)
            self.RC_box2.entry_dmin_coord.set_sensitive(True)
        
            self.RC_box1.label_step_size     .set_sensitive(True)
            #self.RC_box1.label_nsteps        .set_sensitive(True)        
            #self.RC_box1.label_force_constant.set_sensitive(True)        
            #self.RC_box1.label_dmin          .set_sensitive(True)        
            
            self.RC_box2.label_step_size     .set_sensitive(True)
            #self.RC_box2.label_nsteps        .set_sensitive(True)        
            #self.RC_box2.label_force_constant.set_sensitive(True)        
            #self.RC_box2.label_dmin          .set_sensitive(True)        
        
        
        
        if _type == 1:
            self.builder.get_object('label_input_trajectory').show()
            self.builder.get_object('folder_chooser_box').show()
            self.builder.get_object('label_number_of_cpus').show()
            self.builder.get_object('ncpus_spinbutton').show()
            
            self.builder.get_object('label_starting_coordinates').hide()
            self.combobox_starting_coordinates.hide()
            
            self.RC_box1.entry_step_size .set_sensitive(False)
            self.RC_box1.entry_nsteps    .set_sensitive(False)
            self.RC_box1.entry_dmin_coord.set_sensitive(False)
            self.RC_box2.entry_step_size .set_sensitive(False)
            self.RC_box2.entry_nsteps    .set_sensitive(False)
            self.RC_box2.entry_dmin_coord.set_sensitive(False)


            self.RC_box1.label_step_size.set_sensitive(False)
            #self.RC_box1.label_nsteps   .set_sensitive(False)        
            #self.RC_box1.label_dmin     .set_sensitive(False)        
            
            self.RC_box2.label_step_size.set_sensitive(False)
            #self.RC_box2.label_nsteps   .set_sensitive(False)        
            #self.RC_box2.label_dmin     .set_sensitive(False)


    def on_md_integrator_combobox (self, widget = None):
        """ Function doc 
        
            0  :  "Velocity Verlet Dynamics", 
            1  :  "Leap Frog Dynamics",
            2  :  "Langevin Dynamics"
        
        """        
        if widget  == self.builder.get_object('md_integrator_combobox'):
            
            #velocity verlet molecular dynamics
            if 0 == widget.get_active():
                self.builder.get_object('check_pressure_control').hide()
                self.builder.get_object('entry_pressure').hide()
                self.builder.get_object('entry_pressure').hide()
                
                
                self.builder.get_object('entry_temp_coupling').hide()
                self.builder.get_object('label_temp_coupling').hide()
                
                self.builder.get_object('label_pressure_coupling').hide()
                self.builder.get_object('entry_pressure_coupling').hide()
                
                self.builder.get_object('label_collision_freq').hide()
                self.builder.get_object('entry_collision_frequency').hide()
                
                #self.temp_scale_options_combo.set_sensitive(True)
                
            #Leap Frog molecular dynamics
            elif 1 == widget.get_active():
                self.builder.get_object('check_pressure_control').show()
                self.builder.get_object('entry_pressure').show()
                
                self.builder.get_object('label_pressure_coupling').show()
                self.builder.get_object('entry_pressure_coupling').show()
                
                
                self.builder.get_object('entry_temp_coupling').show()
                self.builder.get_object('label_temp_coupling').show()
                
                
                self.builder.get_object('label_collision_freq').hide()
                self.builder.get_object('entry_collision_frequency').hide()
                
                #self.temp_scale_options_combo.set_active(0)
                #self.temp_scale_options_combo.set_sensitive(False)
                
            #Langevin molecular dynamics
            elif 2 == widget.get_active():
                self.builder.get_object('check_pressure_control').hide()
                self.builder.get_object('entry_pressure').hide()
                
                self.builder.get_object('entry_temp_coupling').hide()
                self.builder.get_object('label_temp_coupling').hide()
                
                self.builder.get_object('label_pressure_coupling').hide()
                self.builder.get_object('entry_pressure_coupling').hide()
                
                self.builder.get_object('label_collision_freq').show()
                self.builder.get_object('entry_collision_frequency').show()
            
                #self.temp_scale_options_combo.set_active(0)
                #self.temp_scale_options_combo.set_sensitive(False)
            else:
                pass

        
        
        
        elif widget == self.builder.get_object('combobox_temperature_scale_options'):
            temp_scale_id = widget.get_active()
            
            if temp_scale_id == 0:
                self.builder.get_object('entry_temp_end').set_sensitive(False)
                self.builder.get_object('label_temp_end').set_sensitive(False)
            else:
                self.builder.get_object('entry_temp_end').set_sensitive(True)
                self.builder.get_object('label_temp_end').set_sensitive(True)

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

    
    def run_dialog (self):
        """ Function doc """
        dialog = Gtk.MessageDialog(
                    transient_for=self.main.window,
                    message_type=Gtk.MessageType.ERROR,
                    buttons=Gtk.ButtonsType.OK,
                    text="Folder not found",
                    )
        dialog.format_secondary_text(
            "The folder you have selected does not appear to be valid. Please select a different folder or create a new one."
        )
        dialog.run()
        dialog.destroy()

    
    def update (self, parameters = None):
        """ Function doc """
        
        if self.Visible:
            self._starting_coordinates_model_update()
            self.update_working_folder_chooser()
        
            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
                self.builder.get_object('entry_traj_name').set_text(output_name)
            else:
                pass


    def update_working_folder_chooser (self, folder = None):
        """ Function doc """
        if folder:
            #print('update_working_folder_chooser')
            self.folder_chooser_button.set_folder(folder = folder)
        else:
            
            try:
                folder = self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder
                if folder:
                    self.folder_chooser_button2.set_folder(folder = folder)
                else:
                    pass
            except:
                self.folder_chooser_button2.set_folder(folder = self.home )
   

    def on_checkbox_reaction_coordinate2 (self, widget):
        """ Function doc """
        if widget.get_active():
            self.RC_box2.set_sensitive(True)
        else:
            self.RC_box2.set_sensitive(False) 


    def on_checkbox_geometry_optimization (self, widget):
        """ Function doc """
        if widget.get_active():
            self.builder.get_object('frame_geometry_optimization').set_sensitive(True)
        else:
            self.builder.get_object('frame_geometry_optimization').set_sensitive(False)

   
    def get_input_setup_box_info (self, parameters):
        """ Function doc """
        parameters['input_type'] = self.builder.get_object('combobox_input_type').get_active()
        
        if parameters['input_type'] == 0:
            vobject_id = self.combobox_starting_coordinates.get_vobject_id()
            vobject = self.main.vm_session.vm_objects_dic[vobject_id]
            self.main.p_session.set_psystem_coordinates_from_vobject(vobject)
            parameters['vobj'] = vobject.name
            parameters['source_folder'] = None
            
        elif parameters['input_type'] == 1:
            parameters['vobj'] = None
            parameters['source_folder'] = self.folder_chooser_button.get_folder ()
            
        
        return parameters
    
        
    def get_parameters (self, ):
        
        """ Function doc """
        
        '''this combobox has the reference to the starting coordinates of a simulation'''
        parameters = {"simulation_type" : "Umbrella_Sampling"}

        parameters = self.get_input_setup_box_info ( parameters)
        parameters["folder"] = self.folder_chooser_button2.get_folder()        


        #vobject_id = self.combobox_starting_coordinates.get_vobject_id()
        #vobject = self.main.vm_session.vm_objects_dic[vobject_id]
        #self.main.p_session.set_psystem_coordinates_from_vobject(vobject)

        '''
        
        - If a trajectory folder is given then parallelization is allowed 
            - NmaxThreads >= 1
        
        - If a vobject is provided and a second reaction coordinate is active 
          (2D), then parallelization is allowed
            - NmaxThreads >= 1
        
        - Else: sequetial
        
        '''
        parameters["RC1"] = self.RC_box1.get_rc_data()
        if self.builder.get_object('checkbox_reaction_coordinate2').get_active():
            parameters["RC2"] = self.RC_box2.get_rc_data()
            parameters["NmaxThreads"] =  int(self.builder.get_object('ncpus_spinbutton').get_value())
        else:
            parameters["RC2"] = None
            
            if parameters['input_type'] == 1:
                parameters["NmaxThreads"] = int(self.builder.get_object('ncpus_spinbutton').get_value())
            else:
                parameters["NmaxThreads"] = 1
            
        
        
        #----------------------------------------------------------------------
        #                            GEO OPT
        #----------------------------------------------------------------------
        if self.builder.get_object('checkbox_geometry_optimization').get_active():
            parameters['OPT_parm'] = self._define_OPT_parameters ()
        else:
            parameters['OPT_parm'] = None
        #----------------------------------------------------------------------
        
        
        #----------------------------------------------------------------------
        #                          M. DYNAMICS
        #----------------------------------------------------------------------        
        parameters['MD_parm'] = self._define_MD_parameters ( )
        #----------------------------------------------------------------------
        
        
        parameters['traj_folder_name'] = self.builder.get_object('entry_traj_name').get_text()
        
        return parameters
    
        
    def run (self, button):
        """ Function doc """
        
        parameters = self.get_parameters()
        
        #pprint(parameters)

        isExist = os.path.exists(parameters['folder'])
        if isExist:
            pass
        else:
            self.run_dialog()
            return None
        
        
        self.main.p_session.run_simulation( parameters = parameters )
        self.window.destroy()
        self.Visible    =  False


    def on_btn_export (self, widget):
        """ Function doc """

        #header += '\nparameters = '
        parameters = self.get_parameters()
        try:
            parameters.pop('vobject_name')
        except:
            pass
        #pprint.pprint(parameters)

        export_dialog =  ExportScriptDialog(self.main, parameters = parameters)
        

    def _define_OPT_parameters (self):
        parameters = {}
        opt_id       = self.builder.get_object('combobox_geo_opt').get_active()
        parameters['optimizer']        =  self.opt_methods[opt_id]
        parameters["maximumIterations"]    =  int  (self.builder.get_object('entry_max_int') .get_text())
        parameters["rmsGradientTolerance"]      =  float(self.builder.get_object('entry_rmsd_tol').get_text())
        return parameters
    
    
    def _define_MD_parameters (self):
        """ Function doc """
        integrator_id       = self.builder.get_object('md_integrator_combobox').get_active()
        number_of_steps_eq  = int(self.builder.get_object('entry_number_of_steps_eq').get_text())
        number_of_steps_dc  = int(self.builder.get_object('entry_number_of_steps_dc').get_text())
        
        
        temp_start          = float(self.builder.get_object('entry_temp_start').get_text())
        temp_scale_factor   = int(self.builder.get_object('entry_temp_scale_factor').get_text())
        time_step           = float(self.builder.get_object('entry_time_step').get_text())
        log_frequence       = int(self.builder.get_object('entry_log_frequency').get_text())
        random_seed         = int(self.builder.get_object('entry_random_seed').get_text())
        
        collision_frequency = None
        temp_coupling       = None #float(self.builder.get_object('entry_temp_coupling').get_text())
        
        if self.builder.get_object('check_pressure_control').get_active():
            pressure_control = True
        else:
            pressure_control = False
        
        #-------------------------------------------------------------------------------------------
        if integrator_id == 1: #LeapFrog
            temp_coupling     = float(self.builder.get_object('entry_temp_coupling').get_text())
            pressure          = float(self.builder.get_object('entry_pressure').get_text())
            pressure_coupling = float(self.builder.get_object('entry_pressure_coupling').get_text())
            temp_control      = True
        else:
            pressure          = None
            pressure_coupling = None
            temp_control      = True
        #-------------------------------------------------------------------------------------------
        
        if integrator_id == 2: #Langevin
            temp_coupling       = float(self.builder.get_object('entry_temp_coupling').get_text())
            collision_frequency = float(self.builder.get_object('entry_collision_frequency').get_text())
        else:
            collision_frequency = None
        
        MD_method = {
                     0 : "Verlet"   ,
                     1 : "LeapFrog" ,
                     2 : "Langevin" ,
                     }
        
        

        parameters = {
                    "simulation_type"           : "Umbrella_Sampling",
                    #"folder"                    : HOME                    ,
                    'integrator'                : MD_method[integrator_id], # verlet / leapfrog /langevin
                    'logFrequency'              : log_frequence           ,
                    'seed'                      : random_seed             ,
                    'normal_deviate_generator'  : None                    ,
                    'steps_eq'                  : number_of_steps_eq      ,
                    'steps_dc'                  : number_of_steps_dc      ,
                    'timeStep'                  : time_step               ,
                    #'trajectories'              : None                    ,
                    'trajectory_frequency'          : int(self.builder.get_object('entry_traj_frequency').get_text()), 
                    'trajectory_frequency_dc_ptRes' : int(self.builder.get_object('entry_traj_frequency_dc1').get_text()), 
                    'trajectory_frequency_dc_ptGeo' : int(self.builder.get_object('entry_traj_frequency_dc2').get_text()), 
                    
                    #VelocityVerletDynamics
                    'temperatureScaleFrequency' : temp_scale_factor                ,
                    'temperatureStart'          : temp_start                       ,
                    
                    #LeapFrogDynamics
                    'pressure'                  : pressure           ,  #  LeapFrogDynamics  1.0 default
                    'pressureControl'           : pressure_control   ,  #  LeapFrogDynamics  True / False
                    'pressureCoupling'          : pressure_coupling  ,  #  LeapFrogDynamics  2000
                                                                     
                    'temperature'               : temp_start         ,               
                    'temperatureControl'        : temp_control   ,  # True / False LeapFrogDynamics / LangevinDynamics
                    'temperatureCoupling'       : temp_coupling      ,  #  LeapFrogDynamics / LangevinDynamics
                    
                    #LangevinDynamics
                    'collisionFrequency'        : collision_frequency,  # Langevin 25
                    }
        return parameters


    def restore_the_parameters_to_the_window(self, parameters):
        """Update the GUI widgets with values from the parameters dictionary."""
        MD_ids = {
             "Verlet"   : 0 ,
             "LeapFrog" : 1 ,
             "Langevin" : 2 ,
             }
        
        
        pprint(parameters)
        if not parameters:
            return
        
        #--------------------------------------------------------------
        self.combobox_starting_coordinates.set_active(parameters['cb1_active'])
        #--------------------------------------------------------------    
        
        # ----------------------------
        # Input type
        # ----------------------------
        input_type = parameters.get('input_type', 0)
        self.input_type_combo.set_active(input_type)

        # Set starting coordinates or trajectory folder
        if input_type == 0 and parameters.get('vobj') is not None:
            # Set starting coordinates
            vobj_name = parameters['vobj']
            vobject_ids = self.main.vm_session.vm_objects_dic
            for vobj_id, vobj in vobject_ids.items():
                if vobj.name == vobj_name:
                    self.combobox_starting_coordinates.set_active(vobj_id)
                    break
            self.combobox_starting_coordinates.show()
            self.folder_chooser_button.btn.hide()
            self.spinbutton.hide()
        elif input_type == 1 and parameters.get('source_folder') is not None:
            # Set trajectory folder
            self.folder_chooser_button.set_folder(parameters['source_folder'])
            self.folder_chooser_button.btn.show()
            self.spinbutton.show()
            self.combobox_starting_coordinates.hide()

        # ----------------------------
        # Reaction Coordinates
        # ----------------------------
        if 'RC1' in parameters and parameters['RC1'] is not None:
            self.RC_box1.set_rc_data(parameters['RC1'])

        if 'RC2' in parameters and parameters['RC2'] is not None:
            self.RC_box2.set_rc_data(parameters['RC2'])
            self.builder.get_object('checkbox_reaction_coordinate2').set_active(True)
            self.RC_box2.set_sensitive(True)
        else:
            self.builder.get_object('checkbox_reaction_coordinate2').set_active(False)
            self.RC_box2.set_sensitive(False)

        # ----------------------------
        # Geometry Optimization
        # ----------------------------
        if parameters.get('OPT_parm') is not None:
            self.builder.get_object('checkbox_geometry_optimization').set_active(True)
            self.builder.get_object('frame_geometry_optimization').set_sensitive(True)
            opt_index = list(self.opt_methods.values()).index(parameters['OPT_parm']['optimizer'])
            self.opt_methods_combo.set_active(opt_index)
            self.builder.get_object('entry_max_int').set_text(str(parameters['OPT_parm']['maximumIterations']))
            self.builder.get_object('entry_rmsd_tol').set_text(str(parameters['OPT_parm']['rmsGradientTolerance']))
        else:
            self.builder.get_object('checkbox_geometry_optimization').set_active(False)
            self.builder.get_object('frame_geometry_optimization').set_sensitive(False)

        # ----------------------------
        # Molecular Dynamics
        # ----------------------------
        if parameters.get('MD_parm') is not None:
            md_parm = parameters['MD_parm']
            #integrator = md_parm['integrator']
            integrator_index = MD_ids[md_parm['integrator']]
            
            #integrator_index = list(md_parm['integrator'] for md_parm in md_parm.items() if md_parm).index(md_parm['integrator']) \
            #    if md_parm['integrator'] in ["Verlet", "LeapFrog", "Langevin"] else 0
            self.md_integrators_combobox.set_active(integrator_index)

            self.builder.get_object('entry_number_of_steps_eq').set_text(str(md_parm.get('steps_eq', 0)))
            self.builder.get_object('entry_number_of_steps_dc').set_text(str(md_parm.get('steps_dc', 0)))
            self.builder.get_object('entry_time_step').set_text(str(md_parm.get('timeStep', 0.0)))
            self.builder.get_object('entry_log_frequency').set_text(str(md_parm.get('logFrequency', 0)))
            self.builder.get_object('entry_random_seed').set_text(str(md_parm.get('seed', 0)))
            self.builder.get_object('entry_temp_start').set_text(str(md_parm.get('temperatureStart', 300)))
            self.builder.get_object('entry_temp_scale_factor').set_text(str(md_parm.get('temperatureScaleFrequency', 1)))
            self.builder.get_object('entry_collision_frequency').set_text(str(md_parm.get('collisionFrequency', 0.0)))
            self.builder.get_object('entry_pressure').set_text(str(md_parm.get('pressure', 1.0)))
            self.builder.get_object('entry_pressure_coupling').set_text(str(md_parm.get('temperatureCoupling', 1.0)))
            self.builder.get_object('entry_temp_coupling').set_text(str(md_parm.get('temperatureCoupling', 1.0)))

        # ----------------------------
        # Trajectory folder name
        # ----------------------------
        if 'traj_folder_name' in parameters:
            self.builder.get_object('entry_traj_name').set_text(parameters['traj_folder_name'])

        # ----------------------------
        # Number of threads (CPUs)
        # ----------------------------
        if 'NmaxThreads' in parameters:
            self.spinbutton.set_value(parameters['NmaxThreads'])





















































