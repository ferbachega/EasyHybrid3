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
import os, sys, time
import gi 
import signal
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk, Pango
from gi.repository import GdkPixbuf


from gui.widgets.custom_widgets  import VismolSelectionTypeBox
from gui.widgets.custom_widgets  import FileChooser
from gui.widgets.custom_widgets  import get_colorful_square_pixel_buffer
from gui.widgets.custom_widgets  import ReactionCoordinateBox
from gui.widgets.custom_widgets  import SequenceViewerBox

from gui.windows.setup.windows_and_dialogs import ImportANewSystemWindow
from gui.windows.setup.windows_and_dialogs import EasyHybridDialogSetQCAtoms
from gui.windows.setup.windows_and_dialogs import EasyHybridSetupQCModelWindow
from gui.windows.setup.windows_and_dialogs import EasyHybridGoToAtomWindow
#from gui.windows.setup.windows_and_dialogs import PDynamoSelectionWindow
from gui.windows.setup.windows_and_dialogs import EasyHybridSelectionWindow
from gui.windows.setup.windows_and_dialogs import ExportDataWindow
from gui.windows.setup.windows_and_dialogs import EasyHybridDialogPrune
from gui.windows.setup.windows_and_dialogs import MakeSolventBoxWindow


from gui.windows.setup.windows_and_dialogs import ImportTrajectoryWindow
from gui.windows.setup.windows_and_dialogs import TrajectoryPlayerWindow
from gui.windows.setup.windows_and_dialogs import InfoWindow
from gui.windows.setup.windows_and_dialogs import MergeSystemWindow
from gui.windows.setup.windows_and_dialogs import SolvateSystemWindow
from gui.windows.setup.windows_and_dialogs import SimpleDialog
from gui.windows.setup.edit_frames_dialog import EditFrameDialog

from gui.windows.setup.easyhybrid_terminal    import TerminalWindow
from gui.windows.setup.selection_list_window  import *
from gui.windows.setup.setup_interface        import EasyHybridPreferencesWindow
from gui.windows.setup.process_manager_window import ProcessManagerWindow

from gui.windows.simulation.single_point_window          import SinglePointWindow
from gui.windows.simulation.geometry_optimization_window import GeometryOptimization
from gui.windows.simulation.PES_scan_window              import PotentialEnergyScanWindow 
from gui.windows.simulation.PES_advanced_scan_window     import AdvancedPotentialEnergyScanWindow 
from gui.windows.simulation.molecular_dynamics_window    import MolecularDynamicsWindow 
from gui.windows.simulation.umbrella_sampling_window     import UmbrellaSamplingWindow 
from gui.windows.simulation.chain_of_states_opt_window   import ChainOfStatesOptWindow 
from gui.windows.simulation.normal_modes_window          import NormalModesWindow 

from gui.windows.analysis.WHAM_analysis_window                    import WHAMWindow 
from gui.windows.analysis.normal_modes_analysis_window            import NormalModesAnalysisWindow 
from gui.windows.analysis.surface_analysis_window                 import SurfaceAnalysisWindow 
#from gui.windows.analysis.surface_list_window                     import SurfaceListWindow 
from gui.windows.analysis.energy_refinement_window                import EnergyRefinementWindow
from gui.windows.analysis.PES_analysis_window                     import PotentialEnergyAnalysisWindow
from gui.windows.analysis.distance_angle_dihedral_analysis_window import DistanceAngleDihedralAnalysisWindow
from gui.windows.analysis.RMSD_tool                               import RMSDToolWindow
from gui.windows.analysis.RMSD_analysis_window                    import RMSDAnalysisWindow #/home/fernando/programs/EasyHybrid3/src/gui/windows/analysis/RMSD_analysis_window.py
from gui.windows.analysis.align_trajectory                        import AlignTrajectoryWindow #/home/fernando/programs/EasyHybrid3/src/gui/windows/analysis/RMSD_analysis_window.py
from gui.windows.analysis.reimaging_trajectory                    import ReimagingTrajectoryWindow #/home/fernando/programs/EasyHybrid3/src/gui/windows/analysis/RMSD_analysis_window.py

from util.geometric_analysis import get_simple_distance
from util.sequence_plot import GtkSequenceViewer
from util.rama_plot import RamachandranWindow


from pdynamo.pDynamo2EasyHybrid import pDynamoSession
import numpy as np



from pCore                     import Align                                        , \
                                      Clone                                        , \
                                      logFile                                      , \
                                      Selection                                    , \
                                      TestScript_InputDataPath                     , \
                                      TestScript_OutputDataPath                    , \
                                      XHTMLLogFileWriter

class EasyHybridMainTreeView(Gtk.TreeView):
    """
    Main TreeView widget for managing systems and vismol objects
    in the EasyHybrid application.
    """

    def __init__(self):
        super().__init__()

        # -----------------------------------------------------------
        # Define TreeStore model
        # Each row stores information about systems and vismol objects
        self.treestore = Gtk.TreeStore(
            int,               # 0 system_e_id
            int,               # 1 vobject
            str,               # 2 name

            bool,              # 3 is the radio button visible?
            bool,              # 4 is the radio button active?

            bool,              # 5 is the toggle box visible?
            bool,              # 6 is the toggle box active?

            bool,              # 7 is the Frames visible?
            int,               # 8 number of frames

            GdkPixbuf.Pixbuf   # 9 color square
        )
        self.set_model(self.treestore)

        # List of Gtk.TreeIter objects for quick access
        self.tree_iters = []

        # Dictionary mapping "treeview keys" to TreeIters
        # Each vobject has a key to access its Gtk.TreeIter
        self.tree_iters_dict = {}

        # Counter to generate keys for tree_iters_dict
        self.tree_iters_counter = 0

        # Create TreeView columns and renderers
        self._create_treeview()

    # ---------------------------------------------------------------
    def add_new_system_to_treeview(self, system):
        """
        Add a new system entry to the TreeView.
        """
        sqr_color = get_colorful_square_pixel_buffer(system)

        # Deactivate all existing radio buttons
        for row in self.treestore:
            row[4] = False

        # Append new system row
        parent = self.treestore.append(
            None,
            [
                system.e_id,
                -1,
                f"{system.e_id} - {system.label}",
                True,   # radiobutton visible
                True,   # radiobutton active
                False,  # togglebox visible
                False,  # togglebox active
                False,  # frames visible
                0,      # number of frames
                sqr_color,
            ]
        )

        # Register TreeIter reference for this system
        e_id = system.e_id
        self.main.system_treeview_iters[e_id] = parent

        # Add entry to system ListStore (for global access)
        self.main.system_liststore_iters[e_id] = self.main.system_liststore.append(
            [f"{system.e_id} - {system.label}", system.e_id, sqr_color]
        )

        # Create a ListStore for vismol objects of this system
        self.main.vobject_liststore_dict[system.e_id] = Gtk.ListStore(str, int, int, sqr_color)

        # Update notebook text to reflect active system
        self.main.bottom_notebook.set_active_system_text_to_textbuffer()

    # ---------------------------------------------------------------
    def add_vismol_object_to_treeview(self, vismol_object, vobj_parent=False):
        """
        Add a vismol object (e.g., trajectory, surface, etc.) to the TreeView.
        """
        e_id = vismol_object.e_id
        system = self.main.p_session.psystem[e_id]
        sqr_color = get_colorful_square_pixel_buffer(system)

        print(system, vismol_object, e_id)

        # If a parent is specified (e.g., surface), attach there
        if vobj_parent:
            vismol_object.is_surface = True
            parent = vobj_parent
            sqr_color = None
        else:
            # Otherwise attach to the system node
            parent = self.main.system_treeview_iters[e_id]

        # Number of frames in vismol object
        size = len(vismol_object.frames)

        # Append vismol object row
        _iter = self.treestore.append(
            parent,
            [
                e_id,
                vismol_object.index,
                vismol_object.name,
                False,               # radiobutton visible
                False,               # radiobutton active
                True,                # togglebox visible
                vismol_object.active,  # togglebox active
                True,                # frames visible
                size,                # number of frames
                sqr_color,
            ]
        )

        # Track iter and dictionary reference
        self.tree_iters.append(_iter)
        self.tree_iters_dict[self.tree_iters_counter] = parent
        vismol_object.e_treeview_iter_parent_key = self.tree_iters_counter
        vismol_object.e_treeview_iter = _iter
        self.tree_iters_counter += 1

        # Expand last row (keeps UI responsive for newly added objects)
        for row in self.treestore:
            pass
        self.expand_row(row.path, True)

        # Refresh UI components linked to trajectories
        print('Refresh UI components linked to trajectories')
        self.refresh_trajectory_scalebar()

        # Update sequence viewer with new vismol object
        self.main.bottom_notebook.seqview.refresh_vobject_sequence_list(
            self.main.vm_session.vm_objects_dic
        )

    def _create_treeview (self):
        """ Function doc """
        #treeview = Gtk.TreeView(model=self.treestore)

        # column
        renderer_radio = Gtk.CellRendererToggle()
        renderer_radio.set_radio(True)
        renderer_radio.connect("toggled", self.on_cell_active_radio_toggled)
        column_radio = Gtk.TreeViewColumn("A", renderer_radio, active=4,  visible = 3)
        self.append_column(column_radio)

        # column
        renderer_pixbuf = Gtk.CellRendererPixbuf()
        renderer_text = Gtk.CellRendererText()
        renderer_text.set_property("ellipsize", Pango.EllipsizeMode.END)  # reticências no fim
        renderer_text.set_property("ellipsize-set", True)
        renderer_text.set_property("width", 100)
        #renderer_text.connect("edited", self.text_edited, model, column)
        
        #column_text = Gtk.TreeViewColumn("System", renderer_text, text=2)
        column_text = Gtk.TreeViewColumn("System")#, renderer_text, text=2)
        
        column_text.pack_start(renderer_pixbuf, False)
        column_text.add_attribute(renderer_pixbuf, "pixbuf", 9)
        
        column_text.pack_start(renderer_text, True)
        column_text.add_attribute(renderer_text, "text", 2)
        
        column_text.set_sizing(Gtk.TreeViewColumnSizing.AUTOSIZE)
        column_text.set_resizable(True)
        column_text.set_spacing(10)

        self.append_column(column_text)        
        
        # column
        renderer_toggle = Gtk.CellRendererToggle()
        renderer_toggle.connect("toggled", self.on_cell_visible_toggled)
        column_toggle = Gtk.TreeViewColumn("V", renderer_toggle, active=6, visible = 5)
        self.append_column(column_toggle)

        # column
        renderer_text = Gtk.CellRendererText()
        column_text = Gtk.TreeViewColumn("F", renderer_text, text=8, visible = 7)
        self.append_column(column_text)        
                
        self.connect('row-activated', self.on_select)
        self.connect('button-release-event', self.on_treeview_mouse_button_release_event )

    def on_cell_visible_toggled(self, widget, path):
        ##print(list(path))
        self.treestore[path][6] = not self.treestore[path][6]
        ##print(list(self.treestore[path]))
        
        self.selectedID  = int(self.treestore[path][1])
        vismol_object = self.main.vm_session.vm_objects_dic[self.selectedID]
        vismol_object.active = self.treestore[path][6]
        
        
        self.main.bottom_notebook.seqview.refresh_vobject_sequence_list(self.main.vm_session.vm_objects_dic)
       
        self.refresh_trajectory_scalebar()
        self.main.vm_session.vm_glcore.queue_draw()

    def set_active_system(self, e_id):
        """This function activates a system using the e_id as input."""
        for i, row in enumerate(self.treestore):
            # Only the row with the desired e_id will have the radio button set to True
            row[4] = (i == e_id)
        # Update the active system ID in the session
        self.main.p_session.active_id = e_id

    def on_cell_active_radio_toggled(self, widget, path):
        """Callback triggered when a radio button in the TreeView is toggled."""
        selected_path = Gtk.TreePath(path)

        for row in self.treestore:
            # Store the system ID before switching to the new one
            active_id_before = self.main.p_session.active_id

            # Update the radio button state: True only for the selected row
            row[4] = row.path == selected_path

            if row.path == selected_path:
                # Retrieve the system ID associated with the selected row
                system_e_id = row[0]

                # Retrieve the TreeIter of the selected system
                treeview_iter = self.main.system_treeview_iters[system_e_id]

                # Update the active system ID in the session
                self.main.p_session.active_id = system_e_id
                active_id_after = system_e_id

                # Update annotations in the text buffer with the new active system
                self.main.bottom_notebook.change_annotations_textbuffer(
                    active_id_before,
                    active_id_after
                )
            else:
                pass

        # Refresh the main status bar with updated information
        self.main.refresh_main_statusbar()

        # Update all interface windows and dialogs to reflect the change
        self.main.uptade_interface_windows_and_dialogs()
    
    def on_select(self, tree, path, selection):
        """Triggered when a row in the TreeView is selected."""
        # Get the currently selected row and the model
        model, iter = tree.get_selection().get_selected()        

        # Retrieve values from the selected row
        data2 = model.get_value(iter, 2)
        data1 = model.get_value(iter, 1)
        data0 = model.get_value(iter, 0)

        # Update the trajectory range in the player window
        self.main.trajectory_player_window.change_range(upper=self.treestore[path][8])

        # Store selected vobject index and system ID
        selection             = tree.get_selection()
        model                 = tree.get_model()
        self.vm_object_index  = int(model.get_value(iter, 1))
        self.system_e_id      = int(model.get_value(iter, 0))

        # Update the number of frames for each vobject in the TreeView
        for index, vobject in self.main.vm_session.vm_objects_dic.items():
            treeview_iter = vobject.e_treeview_iter
            vobj_id = self.treestore[treeview_iter][1]

            if vobj_id != -1:
                vobj = self.main.vm_session.vm_objects_dic[vobj_id]
                size = len(vobj.frames)
                self.treestore[treeview_iter][8] = size

    def on_treeview_mouse_button_release_event(self, tree, event):
        """Callback for mouse button release events inside the TreeView."""
        model, iter = tree.get_selection().get_selected()        
        selection   = tree.get_selection()
        model       = tree.get_model()

        try:
            self.vm_object_index = int(model.get_value(iter, 1))
            self.system_e_id     = int(model.get_value(iter, 0))
        except:
            return False
        
        # Right-click → open context menu
        if event.button == 3:
            self.treeview_menu.open_menu(self.system_e_id, self.vm_object_index)

        # Middle-click → center visualization on the selected vobject
        if event.button == 2:
            if self.vm_object_index != -1:  # Ensure it is a vismol object
                vismol_object = self.main.vm_session.vm_objects_dic[self.vm_object_index]
                self.main.vm_session.vm_glcore.center_on_coordinates(
                    vismol_object,
                    vismol_object.mass_center
                )

        # Left-click → currently does nothing
        if event.button == 1:
            pass

    def refresh_number_of_frames (self):
        """Refresh the number of frames in the TreeView for each vobject."""
         
        #.This function refreshes the number of frames  on the main treeview.  
        #.The self.tree_iters list contains all the "parents", or the treeview lines, in the TreeView
        #.vismol_object.e_treeview_iter_parent_key
        
        for index, vobject in self.main.vm_session.vm_objects_dic.items():
            if getattr(vobject, 'e_treeview_iter', False):
                treeview_iter = vobject.e_treeview_iter
                size = len(vobject.frames)

                # Update number of frames and active status
                self.treestore[treeview_iter][8] = size
                self.treestore[treeview_iter][6] = vobject.active

    def refresh_trajectory_scalebar(self):
        """Update the trajectory scalebar range based on the largest active trajectory."""
        higher = 1
        for index, vobject in self.main.vm_session.vm_objects_dic.items():
            #print(index, vobject.name, len(vobject.frames))
            
            if vobject.active:
                size = len(vobject.frames)
                if size > higher:
                    higher = size
        #print('higher', higher)
        # Update the trajectory player window with the new range
        self.main.trajectory_player_window.change_range(upper=higher)
        self.main.trajectory_player_window.upper = higher

    def refresh(self):
        """Rebuild the TreeView with the current systems and vobjects."""
        # Clear the TreeStore before repopulating
        self.treestore.clear()

        # Add all systems to the TreeView
        for e_id in self.main.p_session.psystem.keys():
            system = self.main.p_session.psystem[e_id]
            self.add_new_system_to_treeview(system)

        # Add all vobjects to the TreeView
        for v_obj_index in self.main.vm_session.vm_objects_dic.keys():
            vismol_object = self.main.vm_session.vm_objects_dic[v_obj_index]
            self.add_vismol_object_to_treeview(vismol_object)
