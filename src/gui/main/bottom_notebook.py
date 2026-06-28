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

class BottomNoteBook:
    
    def __init__ (self, main):
        """ Function doc """
        self.home    = main.home
        self.main    = main
        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(self.home,'src/gui/MainWindow_text_and_logs.glade'))
        
        self.widget = self.builder.get_object('notebook_text_and_logs')
        self._define_status_treeview()
    
        #-----------------------------------------------------------------------
        self.annotations_textviewer = self.builder.get_object('annontations_textviewer')
        self.annotations_textbuffer = self.annotations_textviewer.get_buffer()
        #-----------------------------------------------------------------------
        
        
        
        #-----------------------------------------------------------------------
        self.seqview =  GtkSequenceViewer(main = self.main)
        notebook = self.widget
        tab_label = Gtk.Label("-")
        tab_label1 = Gtk.Label("sequence")
        self.seqview.label = tab_label
        vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=3)
        #vbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=10)
        vbox.pack_start(tab_label, False, False, 0)
        vbox.pack_start(self.seqview, True, True, 0)
        
        #notebook.append_page(self.seqview, tab_label)
        notebook.append_page(vbox, tab_label1)
        #-----------------------------------------------------------------------
        self.seqview.add_events(Gdk.EventMask.POINTER_MOTION_MASK)
        self.seqview.add_events(Gdk.EventMask.BUTTON_PRESS_MASK)
        self.seqview.add_events(Gdk.EventMask.BUTTON_RELEASE_MASK)
        #
        self.seqview.connect("motion-notify-event" , self.seqview.on_motion)
        self.seqview.connect("button_press_event"  , self.seqview.button_press)
        self.seqview.connect("button-release-event", self.seqview.button_release)
        self.seqview.set_font_size(17)
        ##-----------------------------------------------------------------------

        self.main.vm_session.vm_widget.connect_after("button-release-event", self.meu_evento_personalizado)
    
    def meu_evento_personalizado (self, widget, event):
        """ Function doc """
        #print(widget, event)
        
        #self.main.bottom_notebook.seqview.text_drawing_area.queue_draw()
        #print(self.main.vm_session.selections[self.main.vm_session.current_selection].active)
    
    def _define_status_treeview (self):
        """   
        This function builds the status treeview
        """
        self.treeview = self.builder.get_object('treeview_status')
        self.status_liststore = Gtk.ListStore(str, # time 
                                              str, # message
                                              str, # logfile_path
                                              GdkPixbuf.Pixbuf,  
                                              )
        self.treeview.set_model(self.status_liststore)
        
        # column
        renderer_text = Gtk.CellRendererText()
        #renderer_text.connect("edited", self.text_edited, model, column)
        column_text = Gtk.TreeViewColumn("time", renderer_text, text=0)
        column_text.set_sizing(Gtk.TreeViewColumnSizing.AUTOSIZE)
        column_text.set_resizable(True)
        self.treeview.append_column(column_text)        
        

        
        # column
        renderer_pixbuf = Gtk.CellRendererPixbuf()
        renderer_text   = Gtk.CellRendererText()
        column_text     = Gtk.TreeViewColumn("message")#, renderer_text, text=2)
        
        column_text.pack_start(renderer_pixbuf, False)
        column_text.add_attribute(renderer_pixbuf, "pixbuf", 3)
        column_text.pack_start(renderer_text, True)
        column_text.add_attribute(renderer_text, "text", 1)
        self.treeview.append_column(column_text)
        
        
        
        #renderer_text = Gtk.CellRendererText()
        #column_text = Gtk.TreeViewColumn("message", renderer_text, text=1, visible = True)
        #self.treeview.append_column(column_text)        
                
        #self.treeview.connect('row-activated', self.on_select)
        #self.treeview.connect('button-release-event', self.on_treeview_mouse_button_release_event )
        self.treeview.set_headers_visible(False)

    def status_teeview_add_new_item (self, message = None , logfile = '', system = None):
        """ Function doc """
        sqr_color      = get_colorful_square_pixel_buffer(system)
        current_time   = time.time()
        formatted_time = time.strftime("%Y-%m-%d   %H:%M:%S", time.localtime(current_time))
        #print(sqr_color)
        
        # Add a new row to the ListStore model
        self.status_liststore.append([formatted_time, message, logfile, sqr_color])
        
        # Get the path of the newly added row
        path = Gtk.TreePath.new_from_indices([len(self.status_liststore)-1])

        # Scroll the TreeView to the newly added row
        self.treeview.scroll_to_cell(path, None, True, 0.5, 0.5)
        return path 
        
    def set_active_system_text_to_textbuffer (self):
        e_id = self.main.p_session.active_id
        new_text = self.main.p_session.psystem[e_id].e_annotations
        self.annotations_textbuffer.set_text(new_text, -1)
    
    def get_active_system_text_from_textbuffer (self):
        """ Function doc """
        start_iter = self.annotations_textbuffer.get_start_iter()
        end_iter = self.annotations_textbuffer.get_end_iter()
        extracted_text = self.annotations_textbuffer.get_text(start_iter, end_iter, include_hidden_chars=False)
        
        e_id = self.main.p_session.active_id
        self.main.p_session.psystem[e_id].e_annotations = extracted_text

    def change_annotations_textbuffer (self, before, after):
        """ Function doc """
        start_iter = self.annotations_textbuffer.get_start_iter()
        end_iter = self.annotations_textbuffer.get_end_iter()
        extracted_text = self.annotations_textbuffer.get_text(start_iter, end_iter, include_hidden_chars=False)
        
        self.main.p_session.psystem[before].e_annotations = extracted_text
        new_text = self.main.p_session.psystem[after].e_annotations
        self.annotations_textbuffer.set_text(new_text, -1)
