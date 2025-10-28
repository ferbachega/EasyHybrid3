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
import os
import pprint
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
#from GTKGUI.gtkWidgets.filechooser import FileChooser
#from easyhybrid.pDynamoMethods.pDynamo2Vismol import *
import gc
from gui.widgets.custom_widgets import FolderChooserButton
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.widgets.custom_widgets import ReactionCoordinateBox
from gui.widgets.custom_widgets import AdvancedReactionCoordinateBox
from gui.windows.setup.windows_and_dialogs import ExportScriptDialog

#from gui.windows.geometry_optimization_window import  FolderChooserButton
VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')




class AdvancedPotentialEnergyScanWindow:
    """GTK-based window for configuring and launching Potential Energy Scan (PES)."""

    def __init__(self, main=None) -> None:
        """Class initializer."""
        self.main = main
        self.home = main.home
        self.p_session = self.main.p_session
        self.vm_session = main.vm_session
        self.Visible: bool = False

        self.residue_liststore = Gtk.ListStore(str, str, str)
        self.opt_methods: dict[int, str] = {
            0: "ConjugatedGradient",
            1: "SteepestDescent",
            2: "LFBGS",
            3: "QuasiNewton",
            4: "FIRE",
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

        self.sym_tag: str = "PES_scan"
        self.last_parameters: dict | None = None  # store last used parameters

    def open_window(self) -> None:
        """Open the PES setup window."""
        if self.Visible:
            self.window.present()
            return

        self.builder = Gtk.Builder()
        self.builder.add_from_file(
            os.path.join(self.home, "src/gui/windows/simulation/PES_scan_window.glade")
        )
        self.builder.connect_signals(self)

        self.window = self.builder.get_object("pes_scan_window")
        self.window.set_title("Reaction Coordinate Scans")
        self.window.set_keep_above(True)

        # Reaction coordinate boxes
        self.RC_box1 = AdvancedReactionCoordinateBox(main = self.main, liststore = self.rc_liststore1)
        self.builder.get_object("rc1_aligment").add(self.RC_box1)
        self.RC_box2 = AdvancedReactionCoordinateBox(main = self.main, liststore = self.rc_liststore2)
        self.builder.get_object("rc2_aligment").add(self.RC_box2)
        #self.RC_box1.set_treeview_size(-1, 300)
        #self.RC_box2.set_treeview_size(-1, 300)        
        
        # Optimizer methods combobox
        self.method_store = Gtk.ListStore(str)
        methods = ["Conjugate Gradient", "FIRE", "L-BFGS", "Steepest Descent"]
        for method in methods:
            self.method_store.append([method])
        self.methods_combo = self.builder.get_object("combobox_methods")
        self.methods_combo.set_model(self.method_store)
        renderer_text = Gtk.CellRendererText()
        self.methods_combo.pack_start(renderer_text, True)
        self.methods_combo.add_attribute(renderer_text, "text", 0)
        self.methods_combo.set_active(0)

        # Coordinates combobox
        self.box_coordinates = self.builder.get_object("box_coordinates")
        self.combobox_starting_coordinates = CoordinatesComboBox()
        self.box_coordinates.pack_start(
            self.combobox_starting_coordinates, False, False, 0
        )
        self._starting_coordinates_model_update(init=True)
        size = len(self.main.vobject_liststore_dict[self.main.p_session.active_id])
        self.combobox_starting_coordinates.set_active(size - 1)

        # Folder chooser
        self.folder_chooser_button = FolderChooserButton(
            main=self.main, sel_type="folder", home=self.home
        )
        self.builder.get_object("folder_chooser_box").pack_start(
            self.folder_chooser_button.btn, True, True, 0
        )
        if self.main.p_session.psystem[self.p_session.active_id]:
            folder = (
                self.main.p_session.psystem[self.p_session.active_id].e_working_folder
                or os.environ.get("HOME")
            )
            self.folder_chooser_button.set_folder(folder=folder)

        # Trajectory name
        if self.p_session.psystem[self.p_session.active_id]:
            output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
            self.builder.get_object("traj_name").set_text(output_name)



        # Buttons
        self.builder.get_object("button_cancel").connect("clicked", self.close_window)
        self.builder.get_object("button_export").connect("clicked", self.on_btn_export)
        self.builder.get_object("checkbtn_TS-centered_mode").connect(
            "toggled", self.change_check_button_TS_centered_mode
        )

        self.change_check_button_reaction_coordinate(None)
        self.build_treeview()
        self.window.show_all()
        #self.builder.get_object("mass_restraints1").hide()
        
        # Configure RC boxes
        self.builder.get_object("checkbtn_TS-centered_mode").hide()
        self.RC_box1.set_rc_mode(rc_mode=0)
        self.RC_box2.set_rc_mode(rc_mode=0)
        #self.RC_box1.set_rc_type(0)
        #self.RC_box2.set_rc_type(0)
        self.Visible = True
        #self.builder.get_object("box_reaction_coordinate2").hide()

    def build_treeview (self):
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
        self.RC_box1.scrolledbox.add(self.treeview1)
        self.RC_box2.scrolledbox.add(self.treeview2)
        self.RC_box1.treeview = self.treeview1
        self.RC_box2.treeview = self.treeview2

    def on_key_press(self, widget, event):
        # checking the Delete key
        if event.keyval == Gdk.KEY_Delete:
            selection = widget.get_selection()
            model, treeiter = selection.get_selected()

            if treeiter is not None:
                model.remove(treeiter)
            return True  # evita propagação do evento

        return False

    def close_window(self, *args) -> None:
        """Close the PES window."""
        self.window.destroy()
        self.Visible = False

    def on_row_inserted(self, model, path, iter):
        if model is  self.rc_liststore1:
            print("Linha inserida em:", path)  
            self.RC_box1.refresh_dmininum()
        else:
            self.RC_box2.refresh_dmininum()

    def on_row_deleted(self, model, path):
        print("Linha removida em:", path)
        if model is  self.rc_liststore1:
            print("Linha inserida em:", path)  
            self.RC_box1.refresh_dmininum()
        else:
            self.RC_box2.refresh_dmininum()
    
    def on_row_changed(self, model, path, iter):
        print("Linha alterada em:", path)
        
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

    def on_button_press_event (self, widget):
        """ Function doc """
        print(widget)
    
    def change_check_button_TS_centered_mode(self, widget: Gtk.ToggleButton) -> None:
        """Enable/disable TS-centered mode."""
        active = self.builder.get_object("checkbtn_TS-centered_mode").get_active()
        mode = 1 if active else 0
        self.RC_box1.set_rc_mode(rc_mode=mode)
        #self.RC_box2.set_rc_mode(rc_mode=mode)

    def change_check_button_reaction_coordinate(self, widget):
        """Enable/disable second reaction coordinate."""
        active = self.builder.get_object("label_check_button_reaction_coordinate2").get_active()
        self.RC_box2.set_sensitive(active)
        self.builder.get_object("n_CPUs_spinbutton").set_sensitive(active)
        self.builder.get_object("n_CPUs_label").set_sensitive(active)

    def run_dialog(self, text = None, secondary_text = None):
        """Show error dialog."""
        if text is None:
             text="Folder not found"

        secondary_text = secondary_text or (
            "The folder you have selected does not appear to be valid. "
            "Please select a different folder or create a new one."
        )
        dialog = Gtk.MessageDialog(
            transient_for=self.main.window,
            message_type=Gtk.MessageType.ERROR,
            buttons=Gtk.ButtonsType.OK,
            text=text,
        )
        dialog.format_secondary_text(secondary_text)
        dialog.run()
        dialog.destroy()

    def run_scan(self, button: Gtk.Button) -> None:
        """Collect parameters and launch PES scan."""
        parameters = self.get_parameters()
        pprint.pprint(parameters)

        if not os.path.exists(parameters["folder"]):
            self.run_dialog()
            return
        
        traj_path = os.path.join(parameters["folder"], parameters["traj_folder_name"] + ".ptGeo")
        if os.path.exists(traj_path):
            self.run_dialog(
                text="Trajectory name not valid!",
                secondary_text="Please provide a new trajectory name as the specified folder already exists",
            )
            return
        
        # Save last run parameters
        self.last_parameters = parameters.copy()
        
        self.p_session.run_simulation(parameters=parameters)

    def get_parameters(self) -> dict:
        """Extract user-defined PES scan parameters from GUI."""
        
        #model = self.treeview1.get_model()
        #lines = []
        #for row in model:
        #    lines.append(row[:])
            #print(row[:]) 
        #print(lines)
        
        
        parameters: dict = {
                "simulation_type": "Advanced_Relaxed_Surface_Scan",
                "logFrequency"   : 50,
                }
        

        parameters["optimizer"] = self.opt_methods[self.methods_combo.get_active()]
        parameters["folder"] = self.folder_chooser_button.get_folder()
        parameters["maximumIterations"] = float(
            self.builder.get_object("entry_max_int").get_text()
        )
        
        parameters["rmsGradientTolerance"] = float(
            self.builder.get_object("entry_rmsd_tol").get_text()
        )
        
        parameters["traj_folder_name"] = self.builder.get_object("traj_name").get_text()
        parameters["vobject_name"] = self.builder.get_object("traj_name").get_text()
        
        vobject_id = self.combobox_starting_coordinates.get_vobject_id()
        vobject = self.main.vm_session.vm_objects_dic[vobject_id]
        parameters['obj1_key6'] = vobject.key6
        
        
        
        #
        self.main.p_session.set_psystem_coordinates_from_vobject(vobject)
        parameters["initial_coordinates"] = vobject.name
        #
        ts_mode = self.builder.get_object("checkbtn_TS-centered_mode").get_active()
        parameters["_is_ts_centered"] = ts_mode
        parameters["RC1"] = self.RC_box1.get_rc_data(ts_mode)
        
        if self.builder.get_object("label_check_button_reaction_coordinate2").get_active():
            parameters['RC2'] = self.RC_box2.get_rc_data(ts_mode)
            parameters["NmaxThreads"] = int(
                 self.builder.get_object("n_CPUs_spinbutton").get_value()
            )
            parameters["traj_type"] = "pklfolder2D"
        else:
            parameters["RC2"] = None
            parameters["NmaxThreads"] = 1
            parameters["traj_type"] = "pklfolder"
            
        
        print(parameters)
        #
        #if self.builder.get_object("label_check_button_reaction_coordinate2").get_active():
        #    parameters["RC2"] = self.RC_box2.get_rc_data(ts_mode)
        #    parameters["NmaxThreads"] = int(
        #        self.builder.get_object("n_CPUs_spinbutton").get_value()
        #    )
        #    parameters["traj_type"] = "pklfolder2D"
        #else:
        #    parameters["RC2"] = None
        #    parameters["NmaxThreads"] = 1
        #    parameters["traj_type"] = "pklfolder"
        #
        return  parameters

    def _starting_coordinates_model_update(self, init: bool = False) -> None:
        """Update coordinates combobox model."""
        e_id = self.main.p_session.active_id
        model = self.main.vobject_liststore_dict[e_id]
        self.combobox_starting_coordinates.set_model(model)
        if model:
            self.combobox_starting_coordinates.set_active(len(model) - 1)

    def update(self, parameters = None) -> None:
        """Refresh window contents when system state changes."""
        
        if self.Visible:
            self._starting_coordinates_model_update()
            self.update_working_folder_chooser()
            if self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system(self.sym_tag)
                self.builder.get_object("traj_name").set_text(output_name)

    def update_working_folder_chooser(self, folder = None) -> None:
        """Update folder chooser button to current working folder."""
        if folder:
            self.folder_chooser_button.set_folder(folder=folder)
        else:
            folder = self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder
            if folder:
                self.folder_chooser_button.set_folder(folder=folder)

    def on_btn_export(self, widget: Gtk.Button) -> None:
        """Export PES setup as script."""
        parameters = self.get_parameters()
        parameters.pop("vobject_name", None)
        ExportScriptDialog(self.main, parameters=parameters)

    def restore_the_parameters_to_the_window(self, parameters = None) -> None:
        """
        Reapply previously used scan parameters back into the GUI.
        If parameters are not provided, use the last stored ones.
        """
        pprint.pprint(parameters)
        if parameters is None:
            parameters = self.last_parameters

        if not parameters:
            print("No previous PES scan parameters available to rerun.")
            return
        
        #--------------------------------------------------------------
        self.combobox_starting_coordinates.set_active(parameters['cb1_active'])
        #--------------------------------------------------------------        
        
        # Optimizer
        opt_map = {v: k for k, v in self.opt_methods.items()}
        self.methods_combo.set_active(opt_map.get(parameters["optimizer"], 0))

        # Numeric fields
        self.builder.get_object("entry_max_int").set_text(
            str(parameters["maximumIterations"])
        )
        self.builder.get_object("entry_rmsd_tol").set_text(
            str(parameters["rmsGradientTolerance"])
        )
        self.builder.get_object("traj_name").set_text(parameters["traj_folder_name"])

        # Folder chooser
        self.folder_chooser_button.set_folder(parameters["folder"])

        # TS-centered mode
        self.builder.get_object("checkbtn_TS-centered_mode").set_active(
            parameters.get("_is_ts_centered", False)
        )

        # Reaction coordinates
        if parameters.get("RC1"):
            self.RC_box1.set_rc_data(parameters["RC1"])
        if parameters.get("RC2"):
            self.builder.get_object("label_check_button_reaction_coordinate2").set_active(True)
            self.RC_box2.set_rc_data(parameters["RC2"])
            self.builder.get_object("n_CPUs_spinbutton").set_value(
                parameters.get("NmaxThreads", 1)
            )
        else:
            self.builder.get_object("label_check_button_reaction_coordinate2").set_active(False)

        

def get_distance (vobject, index1, index2):
    """ Function doc """
    print( index1, index2)
    atom1 = vobject.atoms[index1]
    atom2 = vobject.atoms[index2]
    a1_coord = atom1.coords()
    a2_coord = atom2.coords()
    
    dx = a1_coord[0] - a2_coord[0]
    dy = a1_coord[1] - a2_coord[1]
    dz = a1_coord[2] - a2_coord[2]
    dist = (dx**2+dy**2+dz**2)**0.5
    print('distance a1 - a2:', dist)
    return dist
#===========================================================
def get_angle (vobject, index1, index2, index3):
    """ Function doc """
#===========================================================
def get_dihedral (vobject, index1, index2, index3, index4):
    """ Function doc """
#===========================================================
def compute_sigma_a1_a3 (vobject, index1, index3):

    """ example:
        pk1 ---> pk2 ---> pk3
         N  ---   H  ---  O	    
         
         where H is the moving atom
         calculation only includes N and O ! 
    """
    periodic_table = vobject.vm_session.periodic_table 
    atom1 = vobject.atoms[index1]
    atom3 = vobject.atoms[index3]
    
    symbol1 = atom1.symbol
    symbol3 = atom3.symbol    
    
    mass1 = periodic_table.get_atomic_mass(symbol1)
    mass3 = periodic_table.get_atomic_mass(symbol3)
    #mass1 = atomic_dic[symbol1][4]
    #mass3 = atomic_dic[symbol3][4]
    print(atom1.name, symbol1, mass1)
    print(atom3.name, symbol3, mass3)

    ##pk1_name
    ##pk3_name
    #mass1 = atomic_dic[pk1_name][4]
    #mass3 = atomic_dic[pk3_name][4]
    #
    sigma_pk1_pk3 =  mass1/(mass1+mass3)
    print ("sigma_pk1_pk3: ",sigma_pk1_pk3)
    #
    sigma_pk3_pk1 =  mass3/(mass1+mass3)
    sigma_pk3_pk1 = sigma_pk3_pk1*-1
    #
    print ("sigma_pk3_pk1: ", sigma_pk3_pk1)
    return(sigma_pk1_pk3, sigma_pk3_pk1)
