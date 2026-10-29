#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
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
#      to facilitate QC/MM partitioning and molecular simulations.
#
from util.debug import dprint
from logging import getLogger
logger = getLogger(__name__)
import os, sys, time
import gi 
import signal
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk, Pango
from gi.repository import GdkPixbuf
from gi.repository import GLib
import traceback


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
from gui.windows.setup.windows_and_dialogs import PrepareAmberSystemWindow
from gui.windows.setup.windows_and_dialogs import PrepareOPLSSystemWindow
from gui.windows.setup.windows_and_dialogs import ForceFieldParametersWindow
from gui.windows.setup.windows_and_dialogs import PrepareLigandAntechamberWindow
from gui.windows.setup.windows_and_dialogs import PrepareNamdRunWindow
from gui.windows.setup.windows_and_dialogs import PrepareSMDWindow
from gui.windows.setup.windows_and_dialogs import PreparePackmolWindow
from gui.windows.setup.windows_and_dialogs import PrepareVinaDockingWindow
from gui.windows.setup.windows_and_dialogs import PrepareAutoDockGPUWindow
from gui.windows.setup.windows_and_dialogs import PrepareAddHydrogensWindow
from gui.windows.setup.windows_and_dialogs import SimpleDialog
from gui.windows.setup.edit_frames_dialog import EditFrameDialog
from gui.windows.setup.edit_cell          import EditCellWindow

from gui.windows.setup.easyhybrid_terminal    import TerminalWindow
from gui.windows.setup.selection_list_window  import *
from gui.windows.setup.setup_interface        import EasyHybridPreferencesWindow
from gui.windows.setup.process_manager_window import ProcessManagerWindow

from gui.windows.simulation.single_point_window          import SinglePointWindow
from gui.windows.simulation.geometry_optimization_window import GeometryOptimization
from gui.windows.simulation.transition_state_search_window import TransitionStateSearchWindow
from gui.windows.simulation.reaction_path_window import ReactionPathWindow
from gui.windows.simulation.conjugate_peak_refinement_window import ConjugatePeakRefinementWindow
from gui.windows.simulation.PES_scan_window              import PotentialEnergyScanWindow
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
from gui.windows.analysis.RMSF_analysis_window                    import RMSFAnalysisWindow
from gui.windows.analysis.RDF_analysis_window                      import RDFAnalysisWindow
from gui.windows.analysis.density_analysis_window                 import DensityAnalysisWindow
from gui.windows.analysis.align_trajectory                        import AlignTrajectoryWindow #/home/fernando/programs/EasyHybrid3/src/gui/windows/analysis/RMSD_analysis_window.py
from gui.windows.analysis.reimaging_trajectory                    import ReimagingTrajectoryWindow #/home/fernando/programs/EasyHybrid3/src/gui/windows/analysis/RMSD_analysis_window.py

from util.geometric_analysis import get_simple_distance
from util.sequence_plot import GtkSequenceViewer
from gui.windows.analysis.ramachandran_analysis_window import RamachandranAnalysisWindow

from gui.windows.builder.builder_sidebar      import BuilderSidebarWindow
from gui.windows.builder.fragment_library_window import FragmentLibraryWindow
from gui.windows.builder.structure_library_window import StructureLibraryWindow
from gui.windows.builder.atom_types_window import AtomTypesWindow
from gui.windows.builder.dyff_parameters_window import DYFFParametersWindow
from gui.windows.builder.transform_selection_window import TransformSelectionWindow
from gui.windows.builder.solvate_window import SolvateWindow
from gui.windows.builder.dihedral_angle_window import DihedralAngleWindow


from pdynamo.pDynamo2EasyHybrid import pDynamoSession
import numpy as np



from pCore                     import Align                                        , \
                                      Clone                                        , \
                                      logFile                                      , \
                                      Selection                                    , \
                                      TestScript_InputDataPath                     , \
                                      TestScript_OutputDataPath                    , \
                                      XHTMLLogFileWriter

# --- imports entre modulos adicionados na refatoracao ---
from gui.main.main_treeview import EasyHybridMainTreeView
from gui.main.treeview_menu import TreeViewMenu
from gui.main.bottom_notebook import BottomNoteBook

class MainWindow:
    """
    The main application window for EasyHybrid, a graphical interface for pDynamo3.
    This class handles the initialization of all GUI components, including system
    and coordinate management, treeviews, notebooks, status bars, and auxiliary windows.
    """
    def __init__ (self, vm_session = None, home = None, version = None):
        """
        Initialize the MainWindow instance.

        Parameters:
        vm_session : object
            The virtual molecular session providing visualization and simulation tools.
        home : str
            Home directory path, used for loading resources such as UI files and icons.
        version : str
            The version string of the EasyHybrid application.
        """
        self.home               =  home
        self.EASYHYBRID_VERSION =  version


        '''Search home folder. Every time a file is loaded into memory 
        this attribute will be modified.'''
        self.current_search_folder = None

        


        '''#- - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -#'''
        self.system_liststore        = Gtk.ListStore(str, int, GdkPixbuf.Pixbuf)
        
        #.Connects a system to a treeview postion
        self.system_treeview_iters   = {} #. A dictionary: {e_id : treeview_iter, }
        self.system_liststore_iters  = {} #. A dictionary: {e_id : liststore_iter, }
        
        '''The "vobject_liststore_dict" is a dictionary where the access key is
        the e_id, which is the index of the system of interest generated in 
        "pDynamo2Easyhybrid/pDynamoSession/add_new_system_to_psession". 
        Each dictionary element contains a liststore that includes the respective 
        vobjects.'''
        self.vobject_liststore_dict  = {               
                                       0 : Gtk.ListStore(str,  int, int, GdkPixbuf.Pixbuf)  # name, object_index, e_id, pixel_buffer
                                       }                                 
        '''#- - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -#'''
        self.job_history_liststore = Gtk.ListStore(
            str,               # 0: system name
            str,               # 1: job type
            str,               # 2: potential
            str,               # 3: start time
            str,               # 4: end time
            str,               # 5: status
            GdkPixbuf.Pixbuf,  # 6: color/icon
            int,               # 7: e_id
            int,               # 8: step counter
            str                # 9: system tag (e_tag) -- appended at the
                                #    end deliberately, not inserted after
                                #    column 0, so every EXISTING column-index
                                #    reference elsewhere (process_manager_
                                #    window.py has many) stays correct; the
                                #    "Tag" column's VISUAL position (right
                                #    after "System Name") is controlled by
                                #    treeview column insertion order, which
                                #    is independent of liststore column order.
        )

        # -------------------- GTK BUILDER --------------------
        # Load Glade UI and connect signals
        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(self.home,'src/gui/MainWindow.glade'))
        self.builder.connect_signals(self)
        self.window = self.builder.get_object('window1')
        self.refresh_run_on_menu()     # . Simulate > Run Simulations On (2026-10-28)
        # Restore the dimensions saved from the last session, but only if
        # gl_parameters['save_window_size'] is on (the "Save window size"
        # checkbox in Preferences). If off, always use the fixed 1200x600
        # size -- window_resize() below also respects this flag, so nothing
        # is overwritten in gl_parameters in that case.
        # [BUG FIX] Uses the local PARAMETER 'vm_session', not 'self.vm_session'
        # -- self.vm_session is only assigned further down in this __init__
        # (line ~223), so accessing it here raised AttributeError.
        if vm_session.vm_config.gl_parameters.get('save_window_size', True):
            saved_w = vm_session.vm_config.gl_parameters.get('main_window_width', 1200)
            saved_h = vm_session.vm_config.gl_parameters.get('main_window_height', 600)
        else:
            saved_w, saved_h = 1200, 600
        # Stored to compute the initial position of paned_V further down in
        # this __init__ (see 'self.paned_V_position'), based on the same
        # window height restored here.
        self._restored_window_height = saved_h
        self.window.set_default_size(saved_w, saved_h)                          
        self.window.set_title('EasyHybrid {}'.format(self.EASYHYBRID_VERSION))                          
        
        
        # Drag and drop files
        #-------------------------------------------------------------------
        target = Gtk.TargetEntry.new(
                                         "text/uri-list",
                                         0,
                                         0
                                     )
        self.window.drag_dest_set(
                                     Gtk.DestDefaults.ALL,
                                     [target],
                                     Gdk.DragAction.COPY
                                 )
        self.window.connect("drag-data-received", self.on_drag_data_received)
        #-------------------------------------------------------------------
        
        
        # Set application icon
        self.pixbuf = GdkPixbuf.Pixbuf.new_from_file(os.path.join(self.home,"src/gui/icons/icon_2.png"))
        self.window.set_icon(self.pixbuf)
        
        # Status bar initialization
        self.statusbar_main = self.builder.get_object('statusbar1')
        self.statusbar_main.push(1,'Welcome to EasyHybrid version {}, a pDynamo3 graphical tool'.format(self.EASYHYBRID_VERSION))

        # Hover status bar: atom under the mouse cursor (see update_hover_statusbar)
        self.statusbar_hover = self.builder.get_object('statusbar_hover')
        self.statusbar_hover_context = self.statusbar_hover.get_context_id('hover')
        # right-aligned text: GtkStatusbar packs its label with expand=False
        # and halign=start, so it must be re-packed to fill the bar first
        message_area = self.statusbar_hover.get_message_area()
        for label in message_area.get_children():
            message_area.set_child_packing(label, True, True, 0, Gtk.PackType.START)
            label.set_halign(Gtk.Align.FILL)
            label.set_xalign(1.0)

        self.paned_V         = self.builder.get_object('paned_V')
        
        # -------------------- SESSION MANAGEMENT --------------------
        self.vm_session      = vm_session#( main = None)
        self.vm_session.main = self
        self.vm_session.hover_atom_callback = self.update_hover_statusbar
        
        self.vm_session.vm_object_counter = 0
        self.vm_session.insert_glmenu()
        
        # Connect key press/release events for 3D visualization widget
        self.window.connect("key-press-event",   self.vm_session.vm_widget.key_pressed)
        self.window.connect("key-release-event", self.vm_session.vm_widget.key_released)

        # [BUG FIX] ESC used to close the whole app instantly: vm_widget's
        # own _pressed_Escape (graphics_engine/vismol/gui/vismol_gtkwidget.py)
        # calls self.quit() -> Gtk.main_quit() directly, with no
        # confirmation and completely bypassing on_delete_event's "unsaved
        # changes" check below. Overridden here at the INSTANCE level
        # (vm_widget is a single object, not the class) instead of editing
        # the graphics engine itself, so ESC now goes through the exact
        # same close path as clicking the window's [x] button.
        # . 2026-10-24 (user request): ESC does nothing in the main window -- closing is
        #   done only from the window's [x] or File > Quit (ESC is the key people press to
        #   cancel things, and closed the program when there were no unsaved changes).
        self.vm_session.vm_widget._pressed_Escape = lambda: None


        # -------------------- SELECTION BOX --------------------
        self.menu_box = self.builder.get_object('toolbutton_selection_box')
        self.box2 = self.builder.get_object('box2')
    
        self.selection_box_frame = VismolSelectionTypeBox(vm_session = self.vm_session)
        self.vm_session.selection_box_frame = self.selection_box_frame 
        self.menu_box.add(self.selection_box_frame.box)

        '''#- - - - - - - - - - - -  pDynamo - - - - - - - - - - - - - - -#'''
        self.p_session = pDynamoSession(vm_session = vm_session)
        '''#- - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -#'''


        '''This gtk list is declared in the VismolGLWidget file 
           (it does not depend on the creation of Treeview)'''
        #self.Vismol_Objects_ListStore = self.vm_session.Vismol_Objects_ListStore
        #self.box2.pack_start(self.toolbar_builder, True, True, 1)
        
        #-------------------------------------------------------------------      
        #                         notebook_V1
        #-------------------------------------------------------------------
       
        #-------------------------------------------------------------------      
        #                         notebook_H1
        #-------------------------------------------------------------------
        self.notebook_H1 = Gtk.Notebook()
        self.ScrolledWindow_notebook_H1 = Gtk.ScrolledWindow()
        
        self.main_treeview                = EasyHybridMainTreeView()
        self.main_treeview.main           = self
        self.main_treeview.treeview_menu  = TreeViewMenu(self.main_treeview)
        
        
        self.ScrolledWindow_notebook_H1.add(self.main_treeview)
        
        #self.treeview = GtkEasyHybridMainTreeView(self, vm_session)
        #self.ScrolledWindow_notebook_H1.add(self.treeview)

        # . 2026-10-21: mini trajectory player at the bottom of the Objects treeview
        from gui.widgets.custom_widgets.mini_player import MiniPlayer
        self.mini_player = MiniPlayer(self)
        self.box_objects_page = Gtk.Box(orientation = Gtk.Orientation.VERTICAL, spacing = 0)
        self.box_objects_page.pack_start(self.ScrolledWindow_notebook_H1, True, True, 0)
        self.box_objects_page.pack_start(Gtk.Separator(orientation = Gtk.Orientation.HORIZONTAL), False, False, 0)
        self.box_objects_page.pack_start(self.mini_player, False, False, 0)
        self.notebook_H1.append_page(child = self.box_objects_page, tab_label= Gtk.Label(label = 'Objects'))
        


        # the label we use to show the selection
        self.label = Gtk.Label()
        self.label.set_text("")

        
        #-------------------------------------------------------------------
        #                         notebook_H2
        #-------------------------------------------------------------------
        self.notebook_H2 = Gtk.Notebook()
        #-------------------------------------------------------------------
        
        self.paned_H = Gtk.Paned(orientation = Gtk.Orientation.HORIZONTAL)
        self.button = Gtk.Button(label="Click Here")
        #-------------------------------------------------------------------
        self.vm_session = vm_session#( main = None)
        #-------------------------------------------------------------------
        
        
        self.container = Gtk.Box (orientation = Gtk.Orientation.VERTICAL)
        self.command_line_entry = Gtk.Entry()

        
        if self.vm_session is not None:
            # Add the main visualization widget to the container
            self.container.pack_start(self.vm_session.vm_widget, True, True, 0)
            self.notebook_H2.append_page(child = self.container, tab_label = Gtk.Label(label = 'view'))
            
            # Horizontal paned layout: notebook_H1 on left, notebook_H2 on right
            self.HBOX = Gtk.Box(orientation = Gtk.Orientation.VERTICAL, spacing = 0)
            self.HBOX.pack_start(self.notebook_H1, True, True, 0)

            self.paned_H.add(self.HBOX)
            self.paned_H.add(self.notebook_H2)
            self.paned_H.set_position(250)

            # . paned_V splits the GL/treeview area (top) from the bottom
            #   notebook (Status/Annotations/Sequence). Using plain .add()
            #   here would make the FIRST child (the top, containing the
            #   GL area) fixed-size and the SECOND child (the bottom
            #   notebook) absorb all resizing -- the opposite of what we
            #   want: the 3D view should resize with the window, and the
            #   bottom notebook should keep whatever size the user set for
            #   it. pack1/pack2 let us set "resize" explicitly instead of
            #   relying on add()'s implicit (and here backwards) default.
            self.bottom_notebook = BottomNoteBook(main = self)
            self.paned_V.pack1 (self.paned_H,               True,  True)
            self.paned_V.pack2 (self.bottom_notebook.widget, False, True)
            # [UPDATE] It used to be a fixed value (400px), which made no
            # sense at all if the restored window had a quite different
            # height (see self._restored_window_height, computed up top
            # from gl_parameters['main_window_height']).
            # Now it stores/uses the position as a PROPORTION of the window
            # height (gl_parameters['main_window_paned_v_ratio']) -- so the
            # split between the 3D/treeview area and the bottom notebook
            # follows the window size, instead of being "stuck" at a pixel
            # number that only made sense for the old default size
            # (1200x600). See window_resize() and on_paned_v_position_
            # changed() below, which keep this continuously updated.
            default_paned_v_ratio = 400.0 / 600.0
            paned_v_ratio = self.vm_session.vm_config.gl_parameters.get(
                'main_window_paned_v_ratio', default_paned_v_ratio)
            self.paned_V_position = int(round(paned_v_ratio * self._restored_window_height))
            self.paned_V.set_position(self.paned_V_position)
            self.paned_V.connect("notify::position", self.on_paned_v_position_changed)

            self.bottom_notebook.status_teeview_add_new_item(message = 'Welcome to EasyHybrid, have a happy simulation day!')

        ''' - - - - - - - - - - - - Terminal e Trajectory Play - - - - - - - - - - - - - - - - '''
        self.trajectory_player_button = self.builder.get_object('toolbutton_trajectory_tool1')
        self.cmd_terminal_button      = self.builder.get_object('toolbutton_terminal')
        ''' - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - '''
        
        
        '''Text buffer associated with the command terminal (cmd) output
           
           used in:
                self.terminal_window           = TerminalWindow  (main = self)
            
        '''
        self.terminal_text_buffer    = Gtk.TextBuffer()
        
        
        
        


        #             EASYHYBRID SESSION FILE
        self.session_filename = None
        
        # Timer periodico de autosave (criterio "a cada N minutos" -- o
        # criterio "a cada N eventos" fica em pDynamoSession.register_
        # change_and_maybe_autosave, chamado nas acoes do treeview). Os dois
        # disparam o mesmo p_session._do_autosave(), o que vier primeiro.
        self._autosave_timer_id = None
        self._restart_autosave_timer()
        

        
        '''#- - - - - - - - - - G T K  W I N D O W S - - - - - - - - - - -#'''
        self.window_list = []
        
        self.filechooser                  = FileChooser(main_window = self)
        self.NewSystemWindow              = ImportANewSystemWindow       ( main = self )
        
        self.setup_QCModel_window         = EasyHybridSetupQCModelWindow ( main = self )
        self.window_list.append(self.setup_QCModel_window)
        
        self.geometry_optimization_window = GeometryOptimization  ( main = self )
        self.window_list.append(self.geometry_optimization_window)

        self.transition_state_search_window = TransitionStateSearchWindow ( main = self )
        self.window_list.append(self.transition_state_search_window)

        self.reaction_path_window = ReactionPathWindow ( main = self )
        self.window_list.append(self.reaction_path_window)

        self.cpr_window = ConjugatePeakRefinementWindow ( main = self )
        self.window_list.append(self.cpr_window)
        
        # "Advanced Reaction Coordinate Scans" used to be a separate window
        # (AdvancedPotentialEnergyScanWindow) -- merged into PES_scan_window
        # itself via its checkbtn_advanced_mode toggle. See
        # menuitem_advanced_rc_scans below and PotentialEnergyScanWindow.
        # open_window(advanced=...).
        self.PES_scan_window              = PotentialEnergyScanWindow    ( main=  self)
        self.window_list.append(self.PES_scan_window)

        self.selection_list_window        = SelectionListWindow          ( main = self, system_liststore =  self.system_liststore)
        self.window_list.append(self.selection_list_window)

        
        self.process_manager_window = ProcessManagerWindow( main =  self)
        self.window_list.append(self.process_manager_window)
        
        
        self.go_to_atom_window            = EasyHybridGoToAtomWindow     ( main = self, system_liststore = self.system_liststore)
        self.window_list.append(self.go_to_atom_window)

        self.pDynamo_selection_window     = EasyHybridSelectionWindow    ( main = self)
        
        self.export_data_window           = ExportDataWindow             ( main = self)
        self.window_list.append(self.export_data_window)

        
        self.energy_refinement_window     = EnergyRefinementWindow       ( main = self)
        self.window_list.append(self.energy_refinement_window)
        
        self.single_point_window          = SinglePointWindow            ( main = self)
        self.window_list.append(self.single_point_window)
        #self.window_list.append(self.energy_refinement_window)

        self.import_trajectory_window     = ImportTrajectoryWindow       (main = self)
        self.window_list.append(self.import_trajectory_window)

        
        self.PES_analysis_window          = PotentialEnergyAnalysisWindow(main = self)#, coor_liststore = self.system_liststore)
        self.window_list.append(self.PES_analysis_window)


        self.trajectory_player_window  = TrajectoryPlayerWindow (main = self)
        self.terminal_window           = TerminalWindow  (main = self)
        self.builder_sidebar_window    = BuilderSidebarWindow (main = self)
        self.fragment_library_window   = FragmentLibraryWindow (main = self)
        self.structure_library_window  = StructureLibraryWindow (main = self)
        self.atom_types_window         = AtomTypesWindow (main = self)
        self.dyff_parameters_window    = DYFFParametersWindow (main = self)   # replaces Atom Types in the Builder sidebar
        self.transform_selection_window = TransformSelectionWindow (main = self)
        self.solvate_window             = SolvateWindow (main = self)
        self.dihedral_angle_window      = DihedralAngleWindow (main = self)
        
        self.molecular_dynamics_window  = MolecularDynamicsWindow(main = self)
        self.window_list.append(self.molecular_dynamics_window)
        
        self.umbrella_sampling_window   = UmbrellaSamplingWindow(main = self)
        self.window_list.append(self.umbrella_sampling_window)
        
        
        self.chain_of_states_opt_window = ChainOfStatesOptWindow(main = self)
        self.window_list.append(self.chain_of_states_opt_window)
        
        self.normal_modes_analysis_window =   NormalModesAnalysisWindow (main = self)
        self.surface_analysis_window =   SurfaceAnalysisWindow (main = self)

        self.distance_angle_dihedral_analysis_window = DistanceAngleDihedralAnalysisWindow (main = self)
        

        self.rmsd_analysis_window = RMSDAnalysisWindow (main = self, system_liststore = self.system_liststore)
        self.rmsf_analysis_window = RMSFAnalysisWindow (main = self)
        self.rdf_analysis_window  = RDFAnalysisWindow  (main = self)
        self.density_analysis_window = DensityAnalysisWindow (main = self)
        self.ramachandran_analysis_window = RamachandranAnalysisWindow (main = self)
        
        self.align_trajectory_window = AlignTrajectoryWindow (main = self, system_liststore = self.system_liststore)
        
        self.reimaging_trajectory_window = ReimagingTrajectoryWindow (main = self, system_liststore = self.system_liststore)
        
        
        self.normal_modes_window          =   NormalModesWindow (main = self)
        self.window_list.append(self.normal_modes_window)

        
        self.WHAM_window   =  WHAMWindow(main = self)
        self.window_list.append(self.WHAM_window)

        self.simple_dialog =  SimpleDialog(main = self)

        self.edit_frames_dialog = EditFrameDialog(main = self)
        self.merge_system_window = MergeSystemWindow(main = self)
        # . 2026-10-21: System > Solvate System... (replaces the old Edit > Solvate System)
        from gui.windows.setup.windows_and_dialogs.system_windows.solvate_ready_system import SolvateReadySystemWindow
        self.solvate_system_window = SolvateReadySystemWindow(main = self)
        self.prepare_amber_system_window = PrepareAmberSystemWindow(main = self)
        self.prepare_opls_system_window  = PrepareOPLSSystemWindow(main = self)
        # . 2026-10-17: Force Field Parameters (report + editing) replaces the OPLS Parameters window
        self.opls_parameters_window      = ForceFieldParametersWindow(main = self)
        self.prepare_ligand_antechamber_window = PrepareLigandAntechamberWindow(main = self)
        self.prepare_namd_run_window = PrepareNamdRunWindow(main = self)
        self.prepare_smd_window = PrepareSMDWindow(main = self)
        self.prepare_packmol_window = PreparePackmolWindow(main = self)
        self.prepare_vina_docking_window = PrepareVinaDockingWindow(main = self)
        self.prepare_autodock_gpu_window = PrepareAutoDockGPUWindow(main = self)
        self.prepare_add_hydrogens_window = PrepareAddHydrogensWindow(main = self)
        self.preferences_window = EasyHybridPreferencesWindow(main = self)

        self.make_solvent_box_window = MakeSolventBoxWindow(main = self)

        self.edit_cell_window = EditCellWindow(main = self)
        self.window_list.append(self.edit_cell_window)
        '''#- - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -#'''

        # -------------------- OPEN RECENT MENU --------------------
        # Populated dynamically (not in the .glade) from vm_config's
        # persisted 'recent_files' list -- see gui/config.py. Refreshed
        # on the submenu's own "show" signal so a file that was removed
        # from disk since it was recorded gets pruned right before display.
        self.recent_files_submenu = self.builder.get_object('menu_open_recent')
        self.recent_files_submenu.connect('show', self.on_recent_files_menu_show)
        self.refresh_recent_files_menu()

        # -------------------- WINDOW SIGNALS --------------------
        #self.window.connect("destroy", Gtk.main_quit)
        #self.window.connect("delete-event", Gtk.main_quit)
        #self.window.connect("destroy",      self.quit_easyhybrid)
        self.window.connect("delete-event", self.on_delete_event)
        self.window.connect("check-resize", self.window_resize)
        self.window.show_all()
        # [EN] macOS fix belt-and-suspenders: present() + an idle_add()
        # nudge for the very first frame. Harmless on every platform
        # (present() is a no-op if the window is already focused/shown;
        # the idle callback just asks the glArea/DrawingArea to redraw
        # once the main loop is idle) -- not gated to macOS.
        self.window.present()
        GLib.idle_add(self.vm_session.vm_widget.queue_draw)

    def on_drag_data_received(self, widget, drag_context,
                              x, y, data, info, time):

        for uri in data.get_uris():
            path, _ = GLib.filename_from_uri(uri)
            dprint(path)
            self.vm_session.load(path)

    def restart (self):
        """ Function doc """
        self.main_treeview.treestore.clear()
        self.system_liststore.clear()      
        self.vobject_liststore_dict.clear()
        #self.bottom_notebook.treeview.clear()
        self.session_filename = None
        
        #-------------------------------------------
        self.p_session.restart()
        
        #self.p_session.active_id = 0
        #self.p_session.psystem   =  {
        #                           0:None
        #                          }
        #self.p_session.psystem_name_list    = []
        #self.p_session.psystem_name_counter = 1
        #
        #self.p_session.counter      = 0
        #self.p_session.color_palette_counter = 0
        #-------------------------------------------

        
        self.vm_session.restart()
        self.job_history_liststore.clear()
        
    def window_resize (self, a, b =None, c=None):
        """ Function doc """
        # [BUG FIX] This function used to read w,h and do nothing with them --
        # no window dimension was persisted anywhere. It stores
        # em memoria aqui (barato, roda a cada evento 'check-resize'); a
        # the actual disk write happens in on_delete_event, not here, so as
        # not to write to the config file on every dragged pixel.
        # Respects gl_parameters['save_window_size']: if off, does not
        # atualiza nada (o tamanho fixo 1200x600 continua valendo na
        # proxima abertura).
        if not self.vm_session.vm_config.gl_parameters.get('save_window_size', True):
            return
        w, h = a.get_size()
        self.vm_session.vm_config.gl_parameters['main_window_width']  = int(w)
        self.vm_session.vm_config.gl_parameters['main_window_height'] = int(h)
        
        # Mantem o divisor paned_V (area 3D/treeview vs. notebook de baixo)
        # following the current proportion, instead of staying fixed at a
        # pixel number while the window grows/shrinks. set_position() triggers
        # "notify::position" -> on_paned_v_position_changed(), que so'
        # reconfirms the SAME proportion (does not diverge, just redundant).
        paned_v = getattr(self, 'paned_V', None)
        if paned_v is not None and h > 0:
            ratio = self.vm_session.vm_config.gl_parameters.get('main_window_paned_v_ratio', 400.0/600.0)
            new_pos = int(round(ratio * h))
            if paned_v.get_position() != new_pos:
                paned_v.set_position(new_pos)
    
    def on_paned_v_position_changed (self, paned, gparam):
        """ Atualiza gl_parameters['main_window_paned_v_ratio'] sempre que a
            posicao do paned_V muda -- seja por arrasto do usuario (o caso
            que realmente importa persistir) ou programaticamente (quando
            window_resize acima reaplica a mesma proporcao apos a janela
            mudar de tamanho -- nesse caso o valor recalculado e' o mesmo,
            so' um pouco redundante, sem problema).

            Respeita gl_parameters['save_window_size'] (mesma flag que
            controla o resto da persistencia de layout da janela -- ver
            window_resize/on_delete_event). """
        if not self.vm_session.vm_config.gl_parameters.get('save_window_size', True):
            return
        h = self.window.get_size()[1]
        if h > 0:
            self.paned_V_position = paned.get_position()
            self.vm_session.vm_config.gl_parameters['main_window_paned_v_ratio'] = self.paned_V_position / float(h)
    
    def _restart_autosave_timer (self):
        """ (Re)inicia o timer periodico de autosave, conforme gl_parameters
            ['autosave_interval_minutes']. Chamado uma vez no __init__ e de
            novo sempre que a janela de Preferences aplica mudancas (o
            intervalo pode ter sido alterado la'). """
        if getattr(self, '_autosave_timer_id', None) is not None:
            GLib.source_remove(self._autosave_timer_id)
            self._autosave_timer_id = None
        interval_min = self.vm_session.vm_config.gl_parameters.get('autosave_interval_minutes', 5)
        try:
            interval_min = float(interval_min)
        except (TypeError, ValueError):
            interval_min = 5.0
        # Piso de 10s: evita um valor absurdo (ex.: 0 ou negativo, digitado
        # por engano na Preferences) travar a interface com autosave a
        # cada tick do loop.
        interval_sec = max(int(interval_min * 60), 10)
        self._autosave_timer_id = GLib.timeout_add_seconds(interval_sec, self._on_autosave_timer_tick)

    def _on_autosave_timer_tick (self):
        """ Callback do timer periodico de autosave (criterio "a cada N
            minutos"). So' salva se autosave estiver ligado E houver
            mudancas nao salvas (self.p_session.changed). Retorna True
            sempre, pra o GLib manter o timer se repetindo -- retornar
            False cancelaria o timer. """
        gl_parameters = self.vm_session.vm_config.gl_parameters
        if gl_parameters.get('autosave', True) and self.p_session.changed:
            self.p_session._do_autosave()
        return True
    
    def add_vobject_to_vobject_liststore_dict (self, vismol_object):
        """ Adds a vobject to the vobject liststore. """
        e_id   = vismol_object.e_id
        system = self.p_session.psystem[e_id]
        
        #--------------------------------------------------
        #                  PIXEL BUFFER 
        #--------------------------------------------------
        sqr_color = get_colorful_square_pixel_buffer(system)
        
        
        if getattr(vismol_object,'is_surface', False):
            pass
        else:
            vismol_object.liststore_iter = self.vobject_liststore_dict[e_id].append([vismol_object.name, 
                                                                                    vismol_object.index, 
                                                                                    vismol_object.e_id,
                                                                                    sqr_color]
                                                                                    )

    def get_active_vobject (self, e_id = None):
        """ [EN] BUG FIX: replaces a made-up lookup
        (self.p_session.systems[e_id]['vobject']) that was used at three
        call sites in this file (two added by a fix earlier in this same
        project, one pre-existing at 'toolbutton_system_check') -- all
        three crashed with "AttributeError: 'pDynamoSession' object has
        no attribute 'systems'". self.p_session.systems never actually
        exists anywhere in pDynamoSession (confirmed by grepping the
        whole class); the earlier fix mistakenly copied this exact
        pattern from the pre-existing, apparently equally untested
        'toolbutton_system_check' handler, assuming it was a proven,
        working pattern.

        The only thing that IS real is self.p_session.psystem[e_id],
        which gives the pDynamo "system" wrapper object itself -- but
        that object has no back-reference to its own VismolObject either
        (confirmed: no "system.vobject" assignment exists anywhere in
        session.py). The correct way to get "the" VismolObject for a
        given system is the other direction: every VismolObject already
        stores its OWN system's e_id (vismol_object.e_id, set when the
        object was created/added -- see eSession._add_vismol_object() and
        main_treeview.add_vismol_object_to_treeview(), both of which
        already rely on exactly this attribute), so this searches
        vm_session.vm_objects_dic for the (normally unique) object
        belonging to that system, deliberately excluding is_surface
        objects (a system can have surface VismolObjects attached too,
        sharing the same e_id -- see the Surface Analysis Window work --
        which are not "the" molecule object callers here actually want).

        e_id defaults to the currently active system
        (self.p_session.active_id) if not given explicitly. Returns None
        if nothing matches (e.g. no system loaded yet). """
        if e_id is None:
            e_id = self.p_session.active_id
        for vobject in self.vm_session.vm_objects_dic.values ( ):
            if getattr ( vobject, "e_id", None ) == e_id and not getattr ( vobject, "is_surface", False ):
                return vobject
        return None

    def get_vobjects_for_system ( self, e_id, include_surfaces = False ):
        """ [EN] ALL VismolObjects belonging to system `e_id` -- unlike
        get_active_vobject() above (which returns a single, non-surface
        "the molecule" object), a system can genuinely have MORE than one:
        several docking poses loaded under the same e_id, or one or more
        molecular-surface child objects (is_surface=True, see util/
        molecular_surface.py). `include_surfaces=False` (default) matches
        get_active_vobject()'s own filtering, for callers that only care
        about real, editable molecules (e.g. the Builder's "Edit in
        Builder" own vobject picker, added 2026-09-24 -- see gui/windows/
        builder/builder_entry_dialog.py's choose_and_open_builder()/
        choose_vobject_and_edit()); pass True to also include surface
        objects. No reusable helper for this existed before -- every
        call site inlined its own ad hoc `[v for v in vm_objects_dic.
        values() if v.e_id == e_id ...]` variant. """
        return [ vobject for vobject in self.vm_session.vm_objects_dic.values ( )
                 if getattr ( vobject, "e_id", None ) == e_id
                 and ( include_surfaces or not getattr ( vobject, "is_surface", False ) ) ]

    def on_main_toolbar_clicked (self, button):
        """ Function doc """
        if button  == self.builder.get_object('toolbutton_new_system'):
            self.NewSystemWindow.open_window()
            '''
            from gui.windows.setup.xtb_fragment_charges_window import XtbFragmentChargesWindow
            system = self.p_session.psystem[self.p_session.active_id]
            
            win = XtbFragmentChargesWindow(main=self.window, system=system)
            win.open_window()
            #'''        
        if button  == self.builder.get_object('toolbutton_save'):
            
            if self.session_filename == None:
                self.gtk_save_as_file(button)
                
            else:
                self.p_session.save_easyhybrid_session( filename = self.session_filename)
        
        if button  == self.builder.get_object('toolbutton_save_as'):
            self.gtk_save_as_file (button)
        
        if button  == self.builder.get_object('_show_cell'):
            # [EN] Was "print('Under construction!')" -- the Builder GUI
            # didn't exist at all when this was written (only the
            # backend atom_ops.py/click_mode.py modules did, driven
            # entirely via the terminal's `new`/`placemode`/`add`/`tool`
            # commands -- see easyhybrid_terminal.py). Now that the
            # sidebar exists (gui/windows/builder/builder_sidebar.py),
            # this button opens it -- opening the sidebar IS entering
            # Builder editing mode, see that module's own design note.
            #
            # [EN] DESIGN CHANGE (user's own explicit request): this
            # generic toolbar trigger doesn't already know which molecule
            # (if any) the user wants -- unlike the treeview's own
            # per-row "Edit in Builder" menu item -- so it now asks New
            # vs Edit Existing first (builder_entry_dialog.py) instead of
            # always defaulting straight to a blank canvas.
            #
            # [EN] 2026-10-17, user's request: now an "Edit and Build"
            # TOGGLE. On -> New/Edit chooser; off -> close the session
            # (its own Yes/No/Cancel DYFF question may keep it open). The
            # button then mirrors the real session state, whatever way
            # it was opened/closed -- see BuilderSidebarWindow.
            # sync_main_toggle().
            if getattr ( self, "_syncing_builder_toggle", False ):
                return
            sidebar = self.builder_sidebar_window
            if button.get_active ( ) and not sidebar.visible:
                from gui.windows.builder.builder_entry_dialog import choose_and_open_builder
                choose_and_open_builder ( self )
            elif not button.get_active ( ) and sidebar.visible:
                sidebar.close_window ( allow_cancel = True )
            sidebar.sync_main_toggle ( )

        if button == self.builder.get_object('toolbutton_terminal'):
            if button.get_active ():
                self.terminal_window.open_window()
            else:
                self.terminal_window.close_window(None, None)
        
        if button == self.builder.get_object('toolbutton_trajectory_tool1'):
            if button.get_active ():
                self.trajectory_player_window.open_window()
                #self.traj_frame.hide()
            else:
                self.trajectory_player_window.close_window(button = None)
        
        if button == self.builder.get_object('button_go_to_atom'):
            self.go_to_atom_window.open_window()

        if button  == self.builder.get_object('selections'):
            self.selection_list_window.open_window()

        if button  == self.builder.get_object('toolbutton_energy'):
            #self.p_session.run_simulation (parameters = {'simulation_type' : 'Energy', 'system': self.p_session.psystem[self.p_session.active_id]})
            #self.energy_refinement_window.open_window()
            self.single_point_window.open_window()
            
        if button  == self.builder.get_object('toolbutton_setup_QCModel'):
            # [EN] BUG FIX: this used to call open_window() with no
            # argument at all -- qcmodel_window.py's open_window()
            # defaults vismol_object to None in that case, and
            # on_button_ok() (clicking "OK" in that window) later crashes
            # with "AttributeError: 'NoneType' object has no attribute
            # 'chains'" inside check_charge_fragmentation(), unrelated to
            # the Builder or anything else -- ANY use of this exact
            # toolbar button hit this. Fixed by passing the active
            # vismol_object, same pattern already used two lines below by
            # 'toolbutton_system_check' in this same handler.
            self.setup_QCModel_window.open_window(self.get_active_vobject())
        
        if button  == self.builder.get_object('toolbutton_system_check'): 
            self.get_active_vobject().get_backbone_indexes ()

        if button  == self.builder.get_object('toolbutton_geometry_optimization'):
            self.geometry_optimization_window.open_window()
        
        if button  == self.builder.get_object('toolbutton_pDynamo_selections'):
            self.pDynamo_selection_window.open_window()
            
            '''
            atom1 = self.vm_session.picking_selections.picking_selections_list[0]
            ##print (atom1.chain, atom1.resn, atom1.resi, atom1.name)
            #print ("%s:%s.%s:%s" %(atom1.chain, atom1.resn, atom1.resi, atom1.name))
            
            _centerAtom ="%s:%s.%s:%s" %(atom1.chain, atom1.resn, atom1.resi, atom1.name)
            _radius =  10.0
            
            self.p_session.selections (_centerAtom, _radius)
            '''
       
        if button  == self.builder.get_object('toolbutton_normal_modes'):
            self.normal_modes_window.open_window()
        
        if button  == self.builder.get_object('toolbutton_pes_scan'):
            self.PES_scan_window.open_window()
        
        #if button  == self.builder.get_object('toolbutton_pes_scan'):
        #    self.PES_scan_window.open_window()
        
        if button  == self.builder.get_object('toolbutton_molecular_dynamics'):
            self.molecular_dynamics_window.open_window()
            
        if button  == self.builder.get_object('toolbutton_umbrella_sampling'):
            ##print('toolbutton_umbrella_sampling')
            self.umbrella_sampling_window.open_window()
        
        if button  == self.builder.get_object('toolbutton_chain_of_states_opt'):
            ##print('toolbutton_umbrella_sampling')
            self.chain_of_states_opt_window.open_window()
        
        if button  == self.builder.get_object('toolbutton_monte_carlo'):
            ##print('toolbutton_umbrella_sampling')
            
            self.vm_session.get_dihedral()
            #lista = list(self.vm_session.vobject_names.keys())
            #print(lista)
            #obj = self.vm_session.vobject_names[lista[0]]
            #print(obj.topology)
        
        #if button  == self.builder.get_object('button_test'):
        if button  == self.builder.get_object('button_task_list'):
            ##print('toolbutton_umbrella_sampling')
            self.process_manager_window.open_window()

        if button  == self.builder.get_object('toolbutton_export_img'):
            
            ##print('toolbutton_umbrella_sampling')
            self.save_glArea_image(None)
            #self.vm_session.vm_widget.save_image("saida.png")
    
    def on_main_menu_activate (self, menuitem):
        """ Function doc """
        ##print(menuitem)
        
        if menuitem == self.builder.get_object('menuitem_new'):
            self.NewSystemWindow.open_window()
        
        elif menuitem == self.builder.get_object('menuitem_open'):
            self.open_gtk_load_files (menuitem)
            
        elif menuitem == self.builder.get_object('menuitem_save'):
            if self.session_filename == None:
                self.gtk_save_as_file(menuitem)
            else:
                self.p_session.save_easyhybrid_session( filename = self.session_filename)            

        elif menuitem == self.builder.get_object('menuitem_save_as'):
            self.gtk_save_as_file (menuitem)

        elif menuitem == self.builder.get_object('menuitem_export'):
            self.export_data_window.open_window(sys_selected = self.p_session.active_id)
        
        elif menuitem == self.builder.get_object('menuitem_import'):
            self.import_trajectory_window.open_window()
            
        elif menuitem == self.builder.get_object('menuitem_quit'):
            self.on_delete_event(self, None)
            #Gtk.main_quit()
            pass
        
        
        #----------------------------------------------------------------------
        #                           E D I T 
        #----------------------------------------------------------------------
        elif menuitem == self.builder.get_object('menuitem_make_solvent_box'):
            self.make_solvent_box_window.open_window()
        elif menuitem == self.builder.get_object('menuitem_go_to_atom'):
            self.go_to_atom_window.open_window()
        
        elif menuitem == self.builder.get_object('menuitem_preferences'):
            self.preferences_window.open_window()
            #print(menuitem, 'menuitem_preferences', self.vm_session.vm_glcore.bckgrnd_color)
            #self.vm_session.vm_glcore.bckgrnd_color = [1,1,1,1]
            #print(menuitem, 'menuitem_preferences', self.vm_session.vm_glcore.bckgrnd_color)
            #self.vm_session.vm_config.gl_parameters["line_width"] =20
            #self.vm_session.vm_glcore.light_position = [0, 10, 100.0]
        #----------------------------------------------------------------------
        
        
        #----------------------------------------------------------------------
        #                           S E L E C T I N G
        #----------------------------------------------------------------------
        elif menuitem == self.builder.get_object('menuitem_list'):
            self.selection_list_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_seltype_viewing'):
            if self.selection_box_frame:
                self.selection_box_frame.change_toggle_button_selecting_mode_status(False)
            else:
                self.picking_selection_mode = False
            self.vm_session.vm_glcore.queue_draw()
            
        elif menuitem == self.builder.get_object('menuitem_seltype_picking'):
            if self.selection_box_frame:
                self.selection_box_frame.change_toggle_button_selecting_mode_status(True)
            else:
                self.picking_selection_mode = True
            self.vm_session.vm_glcore.queue_draw()
        
        
        elif menuitem == self.builder.get_object('menuitem_send_to_sel_list'):
            sel_list, sel_resi_table = self.vm_session.build_index_list_from_atom_selection()
            if sel_list:
                
                self.p_session.add_a_new_item_to_selection_list (system_id = self.p_session.active_id, 
                                                                                   indexes = sel_list, 
                                                                                    )
            
            
        elif menuitem == self.builder.get_object('menuitem_extend_selection'):
            self.pDynamo_selection_window.open_window()
        
        elif menuitem == self.builder.get_object('menuitem_by_atom'):
            self.vm_session.viewing_selection_mode(sel_type = 'atom')
        
        elif menuitem == self.builder.get_object('menuitem_by_residue'):
            self.vm_session.viewing_selection_mode(sel_type = 'residue')
        
        elif menuitem == self.builder.get_object('menuitem_by_chain'):
            self.vm_session.viewing_selection_mode(sel_type = 'chain')
        
        elif menuitem == self.builder.get_object('menuitem_by_molecule'):
            self.vm_session.viewing_selection_mode(sel_type = 'molecule')
        
        elif menuitem == self.builder.get_object('menuitem_color'):
            selection               = self.vm_session.selections[self.vm_session.current_selection]
            self.colorchooserdialog = Gtk.ColorChooserDialog()
            
            if self.colorchooserdialog.run() == Gtk.ResponseType.OK:
                color = self.colorchooserdialog.get_rgba()
                #print(color.red,color.green, color.blue )
                new_color = [color.red, color.green, color.blue]

            self.colorchooserdialog.destroy()
            self.vm_session.set_color(color =new_color)
        
        elif menuitem == self.builder.get_object('menuitem_set_as_qc'):
            active_id = self.p_session.active_id
            qc_list, residue_dict, vismol_object = self.vm_session.build_index_list_from_atom_selection(return_vobject = True )
            if qc_list:
                self.p_session.psystem[active_id].e_qc_residue_table = residue_dict
                self.p_session.psystem[active_id].e_qc_table = qc_list
                self.run_dialog_set_QC_atoms(vismol_object = vismol_object)
            
        elif menuitem == self.builder.get_object('menuitem_set_as_fixed'):
            self.vm_session.set_as_fixed_atoms()
        
        elif menuitem == self.builder.get_object('menuitem_set_as_free'):
            self.vm_session.set_as_free_atoms()
        
        elif menuitem == self.builder.get_object('menuitem_prune'):
            #self.vm_session.prune_atoms()
            atomlist, resi_table = self.vm_session.build_index_list_from_atom_selection()

            if atomlist:
                atomlist = list(set(atomlist))
                
                num_of_atoms = len(atomlist)
                name = self.p_session.psystem[self.p_session.active_id].label
                tag  = self.p_session.psystem[self.p_session.active_id].e_tag
                
                dialog =  EasyHybridDialogPrune(self ,num_of_atoms, name, tag)

                if dialog.prune:
                    name         = dialog.name        
                    tag          = dialog.tag  
                    color        = dialog.color 
                    #print('color', color)
                    self.p_session.prune_system (selection = atomlist, name = name, summary = True, tag = tag, color = color)


        #----------------------------------------------------------------------
        #                            V I E W     M E N U
        #----------------------------------------------------------------------
        elif menuitem == self.builder.get_object('menuitem_terminal'):
            button = self.builder.get_object('toolbutton_terminal') 
            self.terminal_window.open_window()
            if self.terminal_window.visible:
                button.set_active (True) 
            else:
                self.terminal_window.close_window(button = None)



        elif menuitem == self.builder.get_object('menuitem_trajectory_tool'):
            button = self.builder.get_object('toolbutton_trajectory_tool1') 
            
            self.trajectory_player_window.open_window()
            if self.trajectory_player_window.Visible:
                button.set_active(True)
            else:
                self.trajectory_player_window.close_window(button = None)

        elif menuitem == self.builder.get_object('menuitem_file_tools'):
            if menuitem.get_active():
                self.builder.get_object('toolbar1').show()
            else:
                self.builder.get_object('toolbar1').hide()
        

        elif menuitem == self.builder.get_object('menuitem_seltype_box'):
            if menuitem.get_active():
                self.builder.get_object('toolbutton_selection_box').show()
            else:
                self.builder.get_object('toolbutton_selection_box').hide()

        
        
        elif menuitem == self.builder.get_object('menuitem_check_selection_toolbar'):
            if menuitem.get_active():
                self.builder.get_object('toolbar2_selections').show()
            else:
                self.builder.get_object('toolbar2_selections').hide()
        
        
        elif menuitem == self.builder.get_object('menuitem_check_pDynamo_tools_bar'):
            if menuitem.get_active():
                self.builder.get_object('toolbar4_pdynamo_tools').show()
            else:
                self.builder.get_object('toolbar4_pdynamo_tools').hide()
        
        elif menuitem == self.builder.get_object('menuitem_message_window'):
            if menuitem.get_active():
                self.bottom_notebook.widget.show()
            else:
                self.bottom_notebook.widget.hide()
        
        #----------------------------------------------------------------------
        elif menuitem == self.builder.get_object('menuitem_check_xyz_axes'):
            if menuitem.get_active():
                self.vm_session.show_axes()
            else:
                self.vm_session.hide_axes()
        
        
        elif menuitem == self.builder.get_object('menuitem_bg_color'):
            self.colorchooserdialog = Gtk.ColorChooserDialog() 
            
            if self.colorchooserdialog.run() == Gtk.ResponseType.OK:
                color = self.colorchooserdialog.get_rgba()
                new_color = [color.red, color.green, color.blue, 1]
            else:
                new_color = False
                self.colorchooserdialog.destroy()
        
            if new_color:
                self.vm_session.vm_glcore.bckgrnd_color = new_color
                self.colorchooserdialog.destroy()
                #system = self.p_session.psystem[self.p_session.active_id]
                #self.change_reference_color(system, new_color)
        #----------------------------------------------------------------------
        
        
        
        
        
        
        
        
        
        
        #----------------------------------------------------------------------
        #                            S Y S T E M 
        #----------------------------------------------------------------------
        elif menuitem == self.builder.get_object('menuitem_info'):
            system = self.p_session.psystem[self.p_session.active_id]
            window = InfoWindow(system)
        
        elif menuitem == self.builder.get_object('menuitem_rename'):
            e_id   = self.p_session.active_id
            system = self.p_session.psystem[e_id]
            # [EN] BUG FIX: this used to construct its own PreferencesWindow
            # directly (never calling .set_names(), so it opened showing the
            # glade's placeholder text instead of the system's actual current
            # name/tag) and stashed it in its own untracked self.preferences,
            # completely independent of treeview_menu's rename_window_visible
            # guard -- so a row-menu rename and a menu-bar rename could each
            # open their own window at the same time. Routing both through
            # the same open_rename_window() fixes both issues at once.
            self.main_treeview.treeview_menu.open_rename_window(
                e_id, -1, system.label, system.e_tag)

        elif menuitem == self.builder.get_object('menuitem_change_working_folder'):
            self.main_treeview.treeview_menu.change_working_folder_for_active_system()

        elif menuitem == self.builder.get_object('menuitem_qc_setup'):
            # [EN] Same bug/fix as 'toolbutton_setup_QCModel' above --
            # see that comment for the full traceback/reasoning.
            self.setup_QCModel_window.open_window(self.get_active_vobject())
            
        elif menuitem == self.builder.get_object('menuitem_show_cell'):
            system = self.p_session.psystem[self.p_session.active_id]
            for key, vobject in self.vm_session.vm_objects_dic.items():
                self.vm_session.show_cell (vobject)
            
            
        elif menuitem == self.builder.get_object('menuitem_hide_cell'):
            system = self.p_session.psystem[self.p_session.active_id]
            for key, vobject in self.vm_session.vm_objects_dic.items():
                self.vm_session.hide_cell (vobject)
            #if self.builder.get_object('toogle_show_cell').get_active():
            #    print('toogle_show_cell - show', self.builder.get_object('toogle_show_cell').get_active())
            #else:
            #    print('toogle_show_cell - hide', self.builder.get_object('toogle_show_cell').get_active())
        
        #elif menuitem == self.builder.get_object('menuitem_cell_and_symmetry'):
        #    print ('menuitem_cell_and_symmetry')

        elif menuitem == self.builder.get_object('menuitem_edit_cell'):
            self.edit_cell_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_center_in_cell'):
            from gui.windows.setup.crystal_tools import center_in_cell
            center_in_cell(self)

        elif menuitem == self.builder.get_object('menuitem_expand_p1'):
            from gui.windows.setup.crystal_tools import expand_to_p1
            expand_to_p1(self)
            
            
            
        elif menuitem == self.builder.get_object('menuitem_change_color'):
            #selection               = self.selections[self.current_selection]
            self.colorchooserdialog = Gtk.ColorChooserDialog()
            
            if self.colorchooserdialog.run() == Gtk.ResponseType.OK:
                color = self.colorchooserdialog.get_rgba()
                new_color = [color.red, color.green, color.blue]
            else:
                new_color = False
                self.colorchooserdialog.destroy()
        
            if new_color:
                self.colorchooserdialog.destroy()
                system = self.p_session.psystem[self.p_session.active_id]
                self.change_reference_color(system, new_color)
        
        elif menuitem == self.builder.get_object('menuitem_remove_restraints'): 
            #print('menuitem_remove_restraints')
            system = self.p_session.psystem[self.p_session.active_id]
            system.e_restraint_counter = 0
            system.e_restraints_dict   = {}
            # 2026-10-26: the dashed lines of the removed restraints stayed in
            # the 3D view and the Selections and Restraints window
            try:
                self.p_session.update_restaint_representation(e_id=system.e_id)
            except Exception:
                dprint(traceback.format_exc())
            try:
                if getattr(self.selection_list_window, 'Visible', False):
                    self.selection_list_window.update_window(selections=False, restraints=True)
            except Exception:
                dprint(traceback.format_exc())

        elif menuitem == self.builder.get_object('menuitem_custom_colors'): 
            system = self.p_session.psystem[self.p_session.active_id]
            name = system.label
            num_of_custom_colors = len(system.e_custom_colors)
            
            msg = 'You are about to remove {} customized atom color list(s) from the system:\n\n{}\n\nDo you want to continue?'.format(num_of_custom_colors,name )
            dialog = SimpleDialog(self )
            yes_or_no = dialog.question (msg)
            if yes_or_no:
                system.e_custom_colors = []
            else:
                pass



            
        elif menuitem == self.builder.get_object('menuitem_remove_fixed'): 
            system = self.p_session.psystem[self.p_session.active_id]
            
            for key, vobject in self.vm_session.vm_objects_dic.items():
                if vobject.e_id == self.p_session.active_id:
                    for index, atom in vobject.atoms.items():
                        atom.color = atom._init_color()
                    
                    #for index in system.e_fixed_table:
                    #    atom = vobject.atoms[index] 
                    #    atom.color = atom._init_color()
                    #    #print(atom, index, atom.index)
                    self.p_session._apply_fixed_representation_to_vobject(system_id = None, vismol_object = vobject)
                    self.p_session._apply_QC_representation_to_vobject   (system_id = None, vismol_object = vobject)                    
            
            self.vm_session.vm_glcore.queue_draw()
            system.e_fixed_table = []
            system.freeAtoms = None            
            self.refresh_widgets()



        elif menuitem == self.builder.get_object('menuitem_fragment_charges'):
            from gui.windows.setup.xtb_fragment_charges_window import XtbFragmentChargesWindow
            system = self.p_session.psystem[self.p_session.active_id]
            
            win = XtbFragmentChargesWindow(main=self, system=system)
            win.open_window()

        elif menuitem == self.builder.get_object('menuitem_esp_charges'):
            from gui.windows.setup.esp_charges_window import ESPChargesWindow
            if not hasattr(self, 'esp_charges_window'):
                self.esp_charges_window = ESPChargesWindow(main=self)
            self.esp_charges_window.open_window()

            
        elif menuitem == self.builder.get_object('menuitem_charge_inspection'): 
            # . 2026-10-24: a window with the diagnosis, the correction and a report (the direct
            #   correction is still p_session.check_for_fragmented_charges())
            from gui.windows.setup.charge_block_inspector_window import ChargeBlockInspectorWindow
            if getattr(self, 'charge_block_inspector_window', None) is None:
                self.charge_block_inspector_window = ChargeBlockInspectorWindow(main = self)
            self.charge_block_inspector_window.open_window()
        
        
        elif menuitem == self.builder.get_object('menuitem_nb_model'): 
            self.p_session.remove_NBModel()
            self.refresh_main_statusbar()
        
        elif menuitem == self.builder.get_object('menuitem_PES_data'): 
            msg = self.p_session.remove_PES_data()
            #if TrueFalse:
            self.simple_dialog.info(msg = msg )
            #else:
            #    self.simple_dialog.error(msg = msg )
        
        
        elif menuitem == self.builder.get_object('menuitem_new_nb_model'):
            self.p_session.define_NBModel()
            self.refresh_main_statusbar()

        elif menuitem == self.builder.get_object('menuitem_force_field_opls'):
            # . 2026-10-17: missing parameters -> offers the Extras > OPLS parametrization
            from gui.windows.setup.opls_assign import assign_opls
            assign_opls ( self, self.p_session.psystem.get ( self.p_session.active_id ) )

        elif menuitem == self.builder.get_object('menuitem_prepare_dyff_system'):
            # . 2026-10-21: System > Force Field > Define Force Field; asks the DYFF version first
            from gui.windows.setup.prepare_dyff_dialog import prepare_dyff_system
            prepare_dyff_system(self)

        elif menuitem == self.builder.get_object('menuitem_parametrization_report'):
            self.open_parametrization_window()

        elif menuitem == self.builder.get_object('menuitem_load_dyff_parametrization'):
            # . 2026-10-23: a saved DYFF modified parametrization (pdynamo/param_report.py) applied again
            self._load_dyff_parametrization_dialog()

        elif menuitem == self.builder.get_object('menuitem_remove_force_field'):
            system = self.p_session.psystem.get(self.p_session.active_id)
            if system is None:
                self.simple_dialog.error(msg="No active system.")
            elif getattr(system, "mmModel", None) is None:
                self.simple_dialog.info(msg="'{}' has no force field.".format(system.label), title="Remove Force Field")
            elif self.simple_dialog.question("Remove the force field (MM and nonbonded models) of '{}'?\n\n"
                                             "Atoms, bonds and coordinates are kept; energies cannot be computed "
                                             "until a force field is defined again.".format(system.label)):
                ok, msg = self.p_session.remove_MMModel(system)
                self.refresh_main_statusbar()
                if ok:
                    self.bottom_notebook.status_teeview_add_new_item(message=msg, system=system)
                else:
                    self.simple_dialog.error(msg=msg)
        elif menuitem == self.builder.get_object('menuitem_force_field_dyff_library'):
            self._dyff_residue_library_dialog()

        elif menuitem == self.builder.get_object('menuitem_force_field_dyff'):
            ok, msg = self.p_session.define_MMModel ( force_field = 'DYFF' )
            self.refresh_main_statusbar()
            if ok:
                self.bottom_notebook.status_teeview_add_new_item ( message = msg, system = self.p_session.psystem.get ( self.p_session.active_id ) )
            else:
                self.simple_dialog.error ( msg = msg )
            
            
            
        elif menuitem == self.builder.get_object('menuitem_remove_qc'): 
            system = self.p_session.psystem[self.p_session.active_id]
            '''
            # The QC region will be deleted. First we have to disable 
            # the QC representation and restore the MM representation 
            # to the atoms.
            '''
            self.p_session.clean_QC_representation_to_vobject()

            system.qcModel = None
            #for vismol_object in self.vm_session.vm_objects_dic.values():
            #    if vismol_object.e_id == system.e_id:
            #        self.p_session._apply_QC_representation_to_vobject   (system_id = None, vismol_object = vismol_object)                   
            self.vm_session.vm_glcore.queue_draw()
            self.refresh_widgets()

        
        elif menuitem == self.builder.get_object('menuitem_solvate'):
            self.solvate_system_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_prepare_amber_tleap'):
            self.prepare_amber_system_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_prepare_opls_system'):
            self.prepare_opls_system_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_opls_parameters'):
            self.opls_parameters_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_prepare_ligand_antechamber'):
            self.prepare_ligand_antechamber_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_run_namd'):
            self.prepare_namd_run_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_run_smd'):
            self.prepare_smd_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_run_packmol'):
            self.prepare_packmol_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_run_vina'):
            self.prepare_vina_docking_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_run_autodock_gpu'):
            self.prepare_autodock_gpu_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_add_missing_hydrogens'):
            self.prepare_add_hydrogens_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_merge'):
            """ Function doc """
            system = self.p_session.psystem[self.p_session.active_id]
            e_id = system.e_id
            self.merge_system_window.selected_system_id = e_id
            self.merge_system_window.open_window()
            
        
        
        
        elif menuitem == self.builder.get_object('menuitem_clone'): 
            system = self.p_session.psystem[self.p_session.active_id]
            
            backup = []
            backup.append(system.e_treeview_iter)
            backup.append(system.e_liststore_iter)
            system.e_treeview_iter   = None
            system.e_liststore_iter  = None
            
            new_system = Clone(system)
            system.e_treeview_iter   = backup[0]
            system.e_liststore_iter  = backup[1]
            #print('menuitem_clone')
            
            new_system = self.p_session.add_new_system_to_psession (system = new_system)
            self.main_treeview.add_new_system_to_treeview (new_system)
            ff  =  getattr(new_system.mmModel, 'forceField', "None")
        
            self.bottom_notebook.status_teeview_add_new_item(message = 'New System:  {} ({}) - Force Field:  {}'.format(new_system.label, new_system.e_tag, ff), system = new_system)
            self.p_session._add_vismol_object_to_easyhybrid_session (new_system, True) #, name = 'olha o  coco')
            #----------------------------------------------------------------------
        elif menuitem == self.builder.get_object('menuitem_prune_to_qc_region'): 
            system = self.p_session.psystem[self.p_session.active_id]
            
            new_system = system.PruneToQCRegion()
            new_system = self.p_session.add_new_system_to_psession (system = new_system)            
            self.main_treeview.add_new_system_to_treeview (new_system)
            # Determine the force field used
            # sys_type = {0: 'AMBER', 1: 'CHARMM'}        
            ff  =  getattr(new_system.mmModel, 'forceField', "None")
            
            
            # Add information about the new system to the status treeview
            self.bottom_notebook.status_teeview_add_new_item(message = 'New System:  {} ({}) - Force Field:  {}'.format(system.label, system.e_tag, ff), system =system)
            
            # Add the system as a vismol object to the easyhybrid session
            self.p_session._add_vismol_object_to_easyhybrid_session (new_system, True)  
        
        #---------------------------------------------------------------------------------------------------------
        
        
        elif menuitem == self.builder.get_object('menuitem_energy'): 
            self.single_point_window.open_window()
            
        elif menuitem == self.builder.get_object('menuitem_geometry_optimization'):
            self.geometry_optimization_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_transition_state_search'):
            self.transition_state_search_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_reaction_path'):
            self.reaction_path_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_cpr'):
            self.cpr_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_molecular_dynamics'):
            self.molecular_dynamics_window.open_window()
            
        elif menuitem == self.builder.get_object('menuitem_normal_modes'):
            self.normal_modes_window.open_window()
            
        elif menuitem == self.builder.get_object('menuitem_rection_coordinate_scans'):
            self.PES_scan_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_advanced_rc_scans'):
            self.PES_scan_window.open_window(advanced=True)

        elif menuitem == self.builder.get_object('menuitem_nudged_elastic_band'):
            self.chain_of_states_opt_window.open_window()
            
        elif menuitem == self.builder.get_object('menuitem_umbrella_sampling'):
            self.umbrella_sampling_window.open_window()
        


        #----------------------------------------------------------------------
        #                        A N A L Y S I S 
        #----------------------------------------------------------------------


        elif menuitem == self.builder.get_object('menuitem_energy_analysis'):
            self.PES_analysis_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_mormal_modes_analysis'):
            self.normal_modes_analysis_window.open_window()

        elif menuitem == self.builder.get_object('test_item'):
            # . the IR spectrum is plotted from the imported normal modes ('Plot IR Spectrum')
            self.normal_modes_analysis_window.open_window()
        
        elif menuitem == self.builder.get_object('menuitem_surface_analysis'):
            # . 2026-10-24: new Orbitals and Surfaces window (energy levels, LCAO); the previous
            #   SurfaceAnalysisWindow is kept in the code for comparison
            from gui.windows.analysis.orbitals_surfaces_window import OrbitalsSurfacesWindow
            if not hasattr(self, 'orbitals_surfaces_window'):
                self.orbitals_surfaces_window = OrbitalsSurfacesWindow(main=self)
            self.orbitals_surfaces_window.open_window()
            #self.surface_list_window.open_window()

        
        elif menuitem == self.builder.get_object('menuitem_energy_refinement'):
            self.energy_refinement_window.open_window()
        
        elif menuitem == self.builder.get_object('menuitem_align_trajectory'):
            #print('Hello!!')
            self.align_trajectory_window.open_window()
            #self.energy_refinement_window.open_window()
        
        elif menuitem == self.builder.get_object('menuitem_reimaging'):
            #print('reimaging_trajectory!!')
            self.reimaging_trajectory_window.open_window()
            #self.energy_refinement_window.open_window()
        
        elif menuitem == self.builder.get_object('menuitem_d_a_d_analysis'):
            self.distance_angle_dihedral_analysis_window.open_window()
        
        elif menuitem == self.builder.get_object('menuitem_rama'):
            self.ramachandran_analysis_window.open_window()
        
        elif menuitem == self.builder.get_object('menuitem_RMSD_tool'):
            #self.rmsd_tool_window.open_window()
            self.rmsd_analysis_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_rmsf_analysis'):
            self.rmsf_analysis_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_excited_states'):
            from gui.windows.simulation.excited_states_window import ExcitedStatesWindow
            if not hasattr(self, 'excited_states_window'):
                self.excited_states_window = ExcitedStatesWindow(main=self)
            self.excited_states_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_molecular_properties'):
            from gui.windows.analysis.molecular_properties_window import MolecularPropertiesWindow
            if not hasattr(self, 'molecular_properties_window'):
                self.molecular_properties_window = MolecularPropertiesWindow(main=self)
            self.molecular_properties_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_trajectory_analysis'):
            from gui.windows.analysis.trajectory_analysis_window import TrajectoryAnalysisWindow
            if not hasattr(self, 'trajectory_analysis_window'):
                self.trajectory_analysis_window = TrajectoryAnalysisWindow(main=self)
            self.trajectory_analysis_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_rdf_analysis'):
            self.rdf_analysis_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_density_analysis'):
            self.density_analysis_window.open_window()
        
        elif menuitem == self.builder.get_object('test_histograms'):
            from util.easyplot import ImagePlot, XYPlot
            import random
            self.plot = XYPlot()
            
            
            logs = [
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window0.ptRes_Values_simpler2.dat' ,
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window1.ptRes_Values_simpler2.dat' ,
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window2.ptRes_Values_simpler2.dat' ,
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window3.ptRes_Values_simpler2.dat' ,
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window4.ptRes_Values_simpler2.dat' ,
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window5.ptRes_Values_simpler2.dat' ,
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window6.ptRes_Values_simpler2.dat' ,
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window7.ptRes_Values_simpler2.dat' ,
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window8.ptRes_Values_simpler2.dat' ,
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window9.ptRes_Values_simpler2.dat' ,
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window10.ptRes_Values_simpler2.dat',
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window11.ptRes_Values_simpler2.dat',
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window12.ptRes_Values_simpler2.dat',
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window13.ptRes_Values_simpler2.dat',
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window14.ptRes_Values_simpler2.dat',
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window15.ptRes_Values_simpler2.dat',
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window16.ptRes_Values_simpler2.dat',
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window17.ptRes_Values_simpler2.dat',
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window18.ptRes_Values_simpler2.dat',
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window19.ptRes_Values_simpler2.dat',
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window20.ptRes_Values_simpler2.dat',
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window21.ptRes_Values_simpler2.dat',
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window22.ptRes_Values_simpler2.dat',
                #'/home/fernando/PROJECTS/OMP/HPC_19_09_2024/5-Hu_OMPD_QC196_AMBER_XTB QC MODEL_QC196_umb_sam/data_collection/window24.ptRes_Values_simpler2.dat',
       
                '/home/fernando/Documents/teste/output_window0.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window1.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window2.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window3.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window4.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window5.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window6.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window7.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window8.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window9.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window10.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window11.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window12.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window13.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window14.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window15.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window16.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window17.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window18.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window19.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window20.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window21.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window22.ptRes_histogram.dat',
                '/home/fernando/Documents/teste/output_window23.ptRes_histogram.dat',
       
             
            # '/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_0.ptRes_Values.dat'  ,
            # '/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_1.ptRes_Values.dat'  ,
            # '/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_2.ptRes_Values.dat'  ,
            # '/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_3.ptRes_Values.dat'  ,
            # '/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_4.ptRes_Values.dat'  ,
            # '/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_5.ptRes_Values.dat'  ,
            # '/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_6.ptRes_Values.dat'  ,
            # '/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_7.ptRes_Values.dat'  ,
            # '/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_8.ptRes_Values.dat'  ,
            # '/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_9.ptRes_Values.dat'  ,
            # '/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_10.ptRes_Values.dat' ,
            # '/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_11.ptRes_Values.dat' ,
            # '/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_12.ptRes_Values.dat' ,
            # '/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_13.ptRes_Values.dat' ,
            # '/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_14.ptRes_Values.dat' ,
            # '/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_15.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_16.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_17.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_18.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_19.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_20.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_21.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_22.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_23.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_24.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_25.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_26.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_27.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_28.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_29.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_30.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_31.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_32.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_33.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_34.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_35.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_0.ptRes_Values.dat'  ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_1.ptRes_Values.dat'  ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_2.ptRes_Values.dat'  ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_3.ptRes_Values.dat'  ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_4.ptRes_Values.dat'  ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_5.ptRes_Values.dat'  ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_6.ptRes_Values.dat'  ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_7.ptRes_Values.dat'  ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_8.ptRes_Values.dat'  ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_9.ptRes_Values.dat'  ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_10.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_11.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_12.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_13.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_14.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_15.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_16.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_17.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_18.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_19.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_20.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_21.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_22.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_23.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_24.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_25.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_26.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_27.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_28.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_29.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_30.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_31.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_32.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_33.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_34.ptRes_Values.dat' ,
            #'/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_psi_35.ptRes_Values.dat' ,           
            ]
            
            
            
            
            for i , log in enumerate(logs):
                
                X = [] 
                Y = [] 
                
                r = random.random()
                g = random.random()
                b = random.random()
                rgb = [r,g,b,]
                
                #print(rgb)
                data = open(log, 'r')
                for line in data:
                    line2 = line.split()
                    X.append(float(line2[0]) )
                    Y.append(float(line2[1]))
                #for i  in range(10):
                #    X.append(i ) 
                #    Y.append(i**2)
                self.plot.add ( X = X, Y = Y,
                                symbol = None, sym_color = [1,1,1], sym_fill = False, 
                                line = 'solid', line_color = rgb, energy_label = None)
            
            #X = []
            #Y = []
            #data = open('/home/fernando/programs/pDynamo3/scratch/examples-3.1.2/bAlaPhiPsiPMF/generatedFiles/bAla_phi_14.ptRes_Values.dat', 'r')
            #
            #for line in data:
            #    line2 = line.split()
            #    X.append(float(line2[0]) )
            #    Y.append(float(line2[1]))
            #self.plot.add ( X = X, Y = Y)
            #self.plot.Ymax = 100
            #print(self.plot.Ymax_list)
            
            self.plot.Ymax_list= [100]
            window =  Gtk.Window()
            window.add(self.plot)
            window.show_all()
        
        
        elif menuitem == self.builder.get_object('menuitem_WHAM'):
            self.WHAM_window.open_window()

        elif menuitem == self.builder.get_object('menuitem_umbrella_browser'):
            self.open_umbrella_browser()
        
        
        
        elif menuitem == self.builder.get_object('menuitem_about'):
            builder = Gtk.Builder()
            builder.add_from_file(os.path.join(self.home,'src/gui/windows/setup/about_dialog.glade'))

            dialog = builder.get_object('about_dialog')
            # . the version comes from _version.py (EASYHYBRID_VERSION), not from the .glade
            if self.EASYHYBRID_VERSION:
                dialog.set_version(str(self.EASYHYBRID_VERSION))
            dialog.set_transient_for(self.window)
            dialog.run()
            dialog.destroy()

    def gtk_save_as_file (self, button):
        """ Function doc """
        dialog = Gtk.FileChooserDialog(
            title="Save File", 
            parent=self.window, 
            action=Gtk.FileChooserAction.SAVE, 
            #buttons=(
            #    Gtk.STOCK_CANCEL, 
            #    Gtk.ResponseType.CANCEL, 
            #    Gtk.STOCK_SAVE, 
            #    Gtk.ResponseType.OK
            #)
            )
        
        dialog.add_buttons(
            Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
            Gtk.STOCK_SAVE, Gtk.ResponseType.OK
            )
        
        '''
        # Add file filters to limit file types
        filter_text = Gtk.FileFilter()
        filter_text.set_name("Text files")
        filter_text.add_mime_type("text/plain")
        dialog.add_filter(filter_text)

        filter_py = Gtk.FileFilter()
        filter_py.set_name("Python files")
        filter_py.add_mime_type("text/x-python")
        dialog.add_filter(filter_py)
        '''
        
        # Set default filename and directory
        dialog.set_current_name("untitled.easy")
        dialog.set_current_folder("~/Documents")
        response = dialog.run()

        if response == Gtk.ResponseType.OK:
            dprint("File selected: " + dialog.get_filename())
            self.p_session.save_easyhybrid_session( filename = dialog.get_filename())
            dialog.destroy()
            # Save file here...
        elif response == Gtk.ResponseType.CANCEL:
            dprint("Save operation canceled.")
            dialog.destroy()

    def open_gtk_load_files (self, button):
        '''Easyhybrid and pkl pdynamo file search '''
        filters = []

        file_filter = Gtk.FileFilter()
        file_filter.set_name("EasyHybrid / pDynamo files (*.easy, *.pkl, *.yaml)")
        file_filter.add_pattern("*.easy")
        file_filter.add_pattern("*.pkl")
        file_filter.add_pattern("*.yaml")

        filters.append(file_filter)

        
        # EasyHybrid session
        file_filter = Gtk.FileFilter()
        file_filter.set_name("EasyHybrid3 session files (*.easy)")
        file_filter.add_pattern("*.easy")
        filters.append(file_filter)

        # pDynamo PKL
        file_filter = Gtk.FileFilter()
        file_filter.set_name("pDynamo3 system files (*.pkl)")
        file_filter.add_pattern("*.pkl")
        filters.append(file_filter)

        # YAML
        file_filter = Gtk.FileFilter()
        file_filter.set_name("YAML files (*.yaml)")
        file_filter.add_pattern("*.yaml")
        filters.append(file_filter)

        # All files
        file_filter = Gtk.FileFilter()
        file_filter.set_name("All files")
        file_filter.add_pattern("*")
        filters.append(file_filter)

        filename = self.filechooser.open(filters=filters)

        if filename:
            self.vm_session.load(filename)

    def refresh_recent_files_menu (self):
        """ Rebuilds the "Open Recent" submenu from vm_config's persisted
            list. Called once at startup and every time the submenu is
            about to be shown (see on_recent_files_menu_show) -- cheap
            enough (at most max_entries items) to just rebuild from
            scratch rather than diffing. """
        for child in self.recent_files_submenu.get_children():
            self.recent_files_submenu.remove(child)

        recent = self.vm_session.vm_config.get_recent_files()

        if not recent:
            placeholder = Gtk.MenuItem(label="(no recent files)")
            placeholder.set_sensitive(False)
            placeholder.show()
            self.recent_files_submenu.append(placeholder)
            return

        for entry in recent:
            path  = entry.get('path', '')
            label = "{}  [{}]".format(os.path.basename(path), entry.get('type', '?'))
            item  = Gtk.MenuItem(label=label)
            item.set_tooltip_text(path)
            item.connect('activate', self.on_recent_file_item_activate, path)
            item.show()
            self.recent_files_submenu.append(item)

        separator = Gtk.SeparatorMenuItem()
        separator.show()
        self.recent_files_submenu.append(separator)

        clear_item = Gtk.MenuItem(label="Clear Recent")
        clear_item.connect('activate', self.on_recent_files_clear_activate)
        clear_item.show()
        self.recent_files_submenu.append(clear_item)

    def on_recent_files_menu_show (self, submenu):
        """ Function doc """
        self.refresh_recent_files_menu()

    def on_recent_file_item_activate (self, menuitem, path):
        """ Function doc """
        self.vm_session.load(path)

    def on_recent_files_clear_activate (self, menuitem):
        """ Function doc """
        self.vm_session.vm_config.clear_recent_files()
        self.refresh_recent_files_menu()

    def run_dialog_set_QC_atoms (self, _type = None, vismol_object = None):
        """ Function doc """
        dialog = EasyHybridDialogSetQCAtoms(self.window)
        response = dialog.run()

        if response == Gtk.ResponseType.YES:
            ##print("The OK button was clicked")
            self.setup_QCModel_window.open_window(vismol_object)
        elif response == Gtk.ResponseType.CANCEL:
            pass
            #print("The Cancel button was clicked")

        dialog.destroy()

    def uptade_interface_windows_and_dialogs (self, parameters = None):
        """ Function doc """
        for window in self.window_list:
            #try:
            window.update()
            #except:
            #    pass    
            #if window.Visible:
            #    window.update()
        
    def change_reference_color (self, system, new_color):
        """ Function doc """
        system.e_color_palette['C'] = new_color
        sqr_color                   = get_colorful_square_pixel_buffer (system)
        e_id =system.e_id
        #print(system.e_liststore_iter, sqr_color)
        self.system_liststore[self.system_liststore_iters[e_id]][2] = sqr_color 
        e_id = system.e_id
        self.main_treeview_model =  self.main_treeview.get_model() 
        #self.main_treeview_model[system.e_treeview_iter][9] = sqr_color
        self.main_treeview_model[self.system_treeview_iters[e_id]][9] = sqr_color
        
        for index, vobject in self.vm_session.vm_objects_dic.items():
            
           
            treeview_iter = vobject.e_treeview_iter
            if self.main_treeview_model[treeview_iter][0] == system.e_id:
                
                #If vobject is a surface, no pixelsqr on treeview
                is_surface = getattr(vobject, 'is_surface', False)
                if is_surface:                
                    self.main_treeview_model[treeview_iter][9] = None
                else:
                    self.main_treeview_model[treeview_iter][9] = sqr_color
                self.vobject_liststore_dict[system.e_id][vobject.liststore_iter][3]=sqr_color
            else:
                pass

    def rename_tag (self, e_id, tag):
        """ Function doc """
        #selection     = self.get_selection()
        #(model, iter) = selection.get_selected()
        #e_id     = model.get_value(iter, 0)
        self.p_session.psystem[e_id].e_tag = tag

    def rename (self, e_id = None, v_id = -1, name = None):
        """ Renames a system's own label (v_id == -1) or one of its
        vismol objects (v_id != -1). Returns True on success, False if
        the new name was rejected (name already taken by something else
        -- an error dialog is shown before returning in that case, so
        callers like PreferencesWindow.on_button_apply can just check the
        return value instead of assuming this always works). """
        if v_id == -1: #.change the header
            # [EN] BUG FIX: systems had NO uniqueness check at all
            # (unlike vobjects just below), so two systems could silently
            # end up with the exact same displayed name/label.
            for other_e_id, other_system in self.p_session.psystem.items():
                if other_e_id != e_id and other_system is not None and other_system.label == name:
                    self.simple_dialog.error(msg='A system named "{}" already exists.'.format(name))
                    return False

            _iter = self.system_treeview_iters[e_id]
            #_iter = self.p_session.psystem[e_id].e_treeview_iter
            self.main_treeview.treestore[_iter][2] = str(e_id)+' - '+ name
            self.p_session.psystem[e_id].label  = name

            liststore_iter = self.system_liststore_iters[e_id]
            #liststore_iter = self.p_session.psystem[e_id].e_liststore_iter
            self.system_liststore[liststore_iter][0] = str(e_id)+' - '+ name
            return True

        else:
            vm_object = self.vm_session.vm_objects_dic[v_id]

            # [EN] BUG FIX: this used to reject ANY name already present
            # in vobject_names -- including the object's OWN current
            # name, which is always already a key there. So clicking
            # Apply without actually changing the name (or re-confirming
            # the same name after editing something else) silently
            # "failed" (dprint only, no dialog) while the caller closed
            # the window anyway, looking like a successful no-op rename.
            if name == vm_object.name:
                return True

            if name in self.vm_session.vobject_names.keys():
                self.simple_dialog.error(msg='An object named "{}" already exists.'.format(name))
                return False

            _iter = vm_object.e_treeview_iter
            self.main_treeview.treestore[_iter][2] = name

            old_name = vm_object.name
            self.vm_session.vobject_names.pop(old_name)

            vm_object.name = name

            self.vm_session.vobject_names[name] = vm_object
            try:
                self.vobject_liststore_dict[e_id][vm_object.liststore_iter][0] = name
            except:
                #means that it is surface
                pass
            return True
    
    def delete_system (self, system_e_id = None ):
        """Remove a system and its associated vobjects from the session.

        This method deletes a system identified by its `system_e_id` from all
        relevant data structures, including:
            1. `vm_objects_dic` in `vm_session`.
            2. `vobject_liststore_dict` in the main object.
            3. `vobject_names` in `vm_session`.
            4. `system_liststore` and `treestore` in the main object.
            5. `psystem` in the `p_session` object.

        Args:
            system_e_id (int, optional): Identifier of the system (access key).
                This value can be obtained from `vobject.e_id`. If `None`, 
                nothing will be removed.

        Returns:
            None
        """      
        
        if system_e_id != None:
            
            removed_vobjects = [ v for v in self.vm_session.vm_objects_dic.values() if v.e_id == system_e_id ]
            '''organizing the list of vobjects that should be removed from vm_object_dic'''
            pop_list = []
            for index, vobject in self.vm_session.vm_objects_dic.items():
                if vobject.e_id == system_e_id:                   
                    self.main_treeview.treestore.remove(vobject.e_treeview_iter)
                    pop_list.append(index)
            
            '''removing from vm_object_dic'''
            for index in pop_list:
                dprint('removing sys = {}, obj = {}, name = {}'.format(system_e_id, 
                                                                      index,
                                                                      self.vm_session.vm_objects_dic[index].name))
                
                name  = self.vm_session.vm_objects_dic[index].name
                self.vm_session.vm_objects_dic.pop(index)
                self.vm_session.vobject_names.pop(name)
            
            '''removing vobject from vobject_liststore_dict'''
            self.vobject_liststore_dict.pop(system_e_id)
            
            
            #  - - - - - - - - removing system treeview items - - - - - - - - - - -
            system = self.p_session.get_system(index = system_e_id)
            e_id = system.e_id 
            self.system_liststore.remove(self.system_liststore_iters[e_id] )
            #self.main_treeview.treestore.remove(system.e_treeview_iter)
            self.main_treeview.treestore.remove(self.system_treeview_iters[e_id])
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
            # Remove the system from p_session and update the graphical view.
            a = self.p_session.delete_system(system_e_id)

            # [BUG FIX] active_id was never updated here -- if the system
            # just removed was the active one (the common case), active_id
            # kept pointing at a key delete_system() had just popped from
            # psystem. Nothing crashed immediately, but the very next
            # thing that touched "the active system" (e.g. File > Save As
            # -> get_active_system_text_from_textbuffer()) blew up with
            # KeyError instead. Point active_id at another remaining
            # system if there is one, otherwise reset to the same empty
            # placeholder state used by pDynamoSession.__init__/restart().
            if self.p_session.active_id not in self.p_session.psystem:
                remaining_ids = list(self.p_session.psystem.keys())
                if remaining_ids:
                    self.p_session.active_id = remaining_ids[-1]
                else:
                    self.p_session.active_id = 0
                    self.p_session.psystem[0] = None

            self._on_vobjects_removed(removed_vobjects)
            self.vm_session.vm_glcore.queue_draw()

    def delete_vm_object (self, vm_object_index = None):
        """ 
        
        vm_object_index = is the access key to the object. You can get it from vobject.index
        
        '''When an object is removed it has to be removed from the treeview and 
        vobject_liststore_dict, in addition to the vm_object_dic in the .vm_session.'''
        
        """
        if vm_object_index != None:
            #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            vobject = self.vm_session.vm_objects_dic[vm_object_index]
            #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
            name = vobject.name
            #  - - - - REMOVING plotting data FROM system.e_logfile_data - - - - - -
            system = self.p_session.psystem[vobject.e_id]
            #print(vobject.e_id, vm_object_index)
            if vm_object_index in system.e_logfile_data.keys():
                dprint ('deleting plotting data...')
                system.e_logfile_data[vm_object_index] = None
                system.e_logfile_data.pop(vm_object_index)
            #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
            #  - - - - - - REMOVING vobj FROM  vobject_liststore_dict - - - - - - -
            if getattr(vobject, 'is_surface', False):
                pass
            else:
                self.vobject_liststore_dict[vobject.e_id].remove(vobject.liststore_iter)
            #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
            #  - - - - - - REMOVING vobj FROM  treestore - - - - - - -
            self.main_treeview.treestore.remove(vobject.e_treeview_iter)
            #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
            #  - - - - - - - - REMOVING vobj FROM vm_object_dic - - - - - - - - - -
            self.vm_session.vm_objects_dic[vm_object_index] = None
            self.vm_session.vm_objects_dic.pop(vm_object_index)# = None
            self.vm_session.vobject_names.pop(name)# = None
            #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            self._on_vobjects_removed([vobject])
            self.vm_session.vm_glcore.queue_draw()

    def _on_vobjects_removed(self, vobjects):
        """ [EN] 2026-10-04 BUG FIX (user's report: a system deleted from the
        treeview while the Builder was editing it left "residual atoms and
        indices" in later Builder objects). Called by delete_system() and
        delete_vm_object() with the objects just removed:
          - their atoms leave vm_session.atom_dic_id (picking id -> Atom;
            never cleaned before, every deleted system leaked its atoms);
          - if the Builder was editing one of them, the session ends
            silently (builder_sidebar.close_window() sees the object is gone:
            no "incorporate" question, nothing rebuilt), which also closes
            the Builder's DYFF Parameters / torsion fit windows -- they kept
            showing and labelling the deleted object's atoms;
          - a DYFF Parameters window left open on one of them is closed. """
        vobjects = [v for v in vobjects if v is not None]
        if not vobjects:
            return
        atom_dic_id = self.vm_session.atom_dic_id
        for vobject in vobjects:
            for atom in vobject.atoms.values():
                if atom_dic_id.get(getattr(atom, 'unique_id', None)) is atom:
                    atom_dic_id.pop(atom.unique_id, None)
        removed = set(id(v) for v in vobjects)
        self._drop_removed_atoms_from_selections(removed)
        target = getattr(self.vm_session, 'builder_target_object', None)
        sidebar = getattr(self, 'builder_sidebar_window', None)
        if target is not None and id(target) in removed and sidebar is not None and getattr(sidebar, 'visible', False):
            sidebar.close_window()
            self.statusbar_main.push(1, "Builder closed: the object being edited ('{}') was deleted.".format(target.name))
        elif target is not None and id(target) in removed:
            self.vm_session.builder_target_object = None
        dyff_window = getattr(self, 'dyff_parameters_window', None)
        if dyff_window is not None and dyff_window.visible and id(getattr(dyff_window, 'target_object', None)) in removed:
            dyff_window.close_window()

    def show_sidebar_page(self, name, title=None):
        """ [EN] 2026-10-17: right side panel (sidebar_root_box, LibreOffice
        style) -- shows page `name` of sidebar_stack (e.g. "builder_add") with
        `title` in its header. Tools can add their own pages to
        self.builder.get_object('sidebar_stack') with add_titled(). """
        box   = self.builder.get_object('sidebar_root_box')
        stack = self.builder.get_object('sidebar_stack')
        page  = stack.get_child_by_name(name)
        if page is None:
            return
        if title is None:
            title = stack.child_get_property(page, 'title') or ''
        self.builder.get_object('sidebar_title_label').set_text(title)
        page.show_all()
        stack.set_visible_child(page)
        box.set_no_show_all(False)
        box.show_all()

    def hide_sidebar(self):
        """ [EN] 2026-10-17: hides the right side panel (see show_sidebar_page()). """
        box = self.builder.get_object('sidebar_root_box')
        box.hide()
        box.set_no_show_all(True)

    def _drop_removed_atoms_from_selections(self, removed):
        """ [EN] 2026-10-09 BUG FIX (user's report: Builder object created and
        deleted -> in the next object pk1/pk2, as projected by the DYFF
        Parameters window, pointed to the deleted object). The picking list
        (pk1..pk4) kept the deleted object's Atom objects and their marker
        spheres stayed on screen; a click on the new object filled only the
        EMPTY slots (pk3/pk4), so pk1/pk2 stayed on the deleted atoms (and
        dihedral/bond tools used that object's topology). Removed atoms now
        leave the picking list (the survivors are re-picked in order, which
        also redraws the markers) and every viewing selection.
        removed: set of id() of the removed vobjects. """
        session = self.vm_session
        picking = getattr(session, 'picking_selections', None)
        if picking is not None:
            pk_list = picking.picking_selections_list
            if any(a is not None and id(a.vm_object) in removed for a in pk_list):
                survivors = [a for a in pk_list if a is not None and id(a.vm_object) not in removed]
                picking.selection_function_picking(None)
                for atom in survivors:
                    picking.selection_function_picking(atom)
        for selection in getattr(session, 'selections', {}).values():
            atoms = getattr(selection, 'selected_atoms', None)
            if atoms:
                for atom in [a for a in atoms if id(a.vm_object) in removed]:
                    atoms.discard(atom)
            objects = getattr(selection, 'selected_objects', None)
            if objects:
                for vobject in [v for v in objects if id(v) in removed]:
                    objects.discard(vobject)


    def _dyff_residue_library_dialog(self, as_new=False):
        """ [EN] 2026-10-07: System > Add > Force Field > DYFF + residue library.
        Residues are matched to the precomputed xTB library (src/pdynamo/residue_library.py);
        the rest gets CM5 charges on a capped subsystem (worker thread); then DYFF is
        defined with the library edits and the charges on the system's object.
        as_new (2026-10-21): applied to a copy of the system (_dyff_copy_system). """
        import threading
        from pdynamo import residue_library as RL
        system = self.p_session.psystem.get(self.p_session.active_id)
        if system is None:
            self.simple_dialog.error(msg="No active system."); return
        if getattr(system, "qcModel", None) is not None:
            self.simple_dialog.error(msg="The system has a QC model: remove the QC region before redefining the MM model."); return
        # . the library matches residues by their bond graph and DYFF types atoms from bond
        #   orders: a PDB read without CONECT records has no bonds at all (every atom would
        #   become its own "rest" group, hundreds of xTB jobs)
        if self._dyff_needs_preparation(system):
            # . 2026-10-08: a raw structure goes through Prepare DYFF System (no OPLS typing:
            #   hydrogens, protonation, metal sites, truncated residues), which then applies
            #   the library itself
            from gui.windows.setup.prepare_dyff_dialog import PrepareDYFFDialog
            PrepareDYFFDialog(self).run(system)
            return
        vobject = next((v for v in self.vm_session.vm_objects_dic.values() if getattr(v, "e_id", None) == system.e_id), None)
        if vobject is None or len(vobject.atoms) != len(system.atoms):
            self.simple_dialog.error(msg="No object with the same atoms as the active system."); return
        try:
            p = RL.plan(system)
        except Exception as error:
            self.simple_dialog.error(msg="Residue library: {}".format(error)); return
        s = p["summary"]
        dialog = Gtk.Dialog(title="DYFF + residue library" + (" (new system)" if as_new else ""), transient_for=self.window, modal=True)
        dialog.add_buttons("Cancel", Gtk.ResponseType.CANCEL, "Apply", Gtk.ResponseType.OK)
        grid = Gtk.Grid(column_spacing=8, row_spacing=6, margin=10)
        text = ("{}: {} of {} residues matched in the library ({}).\n"
                "Not in the library: {}{}.\n{} atoms outside the library will get CM5 charges (capped subsystem, xTB)."
                ).format(system.label, s["matched"], s["residues"],
                         ", ".join("{} {}".format(n, t) for t, n in sorted(s["templates"].items())) or "none",
                         ", ".join(s["unmatched"][:12]) or "none", " ..." if len(s["unmatched"]) > 12 else "", s["rest_atoms"])
        if s.get("waters") or s.get("ions"):
            text += ("\n{} water molecules (TIP3P charges and oxygen LJ on DYFF bonded terms -- not a validated water model) "
                     "and {} monatomic ions (formal charge) get fixed parameters, no xTB.").format(s.get("waters", 0), s.get("ions", 0))
        label = Gtk.Label(label=text, xalign=0); label.set_line_wrap(True); label.set_max_width_chars(90)
        grid.attach(label, 0, 0, 2, 1)
        solvent = Gtk.CheckButton(label="CM5 of the rest in implicit water (ALPB)"); solvent.set_active(True)
        grid.attach(solvent, 0, 1, 2, 1)
        bonded = Gtk.CheckButton(label="Outside the library also: bonds/angles from the xTB Hessian + torsions from the DFT torsion library")
        bonded.set_tooltip_text("Same capped subsystem: modified Seminario on its xTB Hessian (calibrated to DFT by chemical environment), "
                                "then the per-bond torsion series of pdynamo/bond_library.py (found by chemical environment; DYFF rules "
                                "where none is found). Off: CM5 charges only, DYFF bonded rules.")
        bonded.set_active(bool(s["rest_atoms"])); bonded.set_sensitive(bool(s["rest_atoms"]))
        grid.attach(bonded, 0, 2, 2, 1)
        method_combo, scale_spin = self._dyff_method_and_scale_widgets(grid, 3)
        scan = self._dyff_scan_torsions_widget(bonded)
        grid.attach(scan, 0, 5, 2, 1)
        # . 2026-10-21, user's request: charge and multiplicity of the xTB calculations
        charge_spin, mult_spin, electrons = self._dyff_charge_multiplicity_widgets(grid, 6, system, vobject)
        grid.attach(Gtk.Label(label="Ligand charges:", xalign=0), 0, 8, 1, 1)
        ligand_combo = self._dyff_ligand_charges_widget()
        ligand_combo.set_sensitive(bool(s["rest_atoms"]))
        grid.attach(ligand_combo, 1, 8, 1, 1)
        charge_spin.set_sensitive(bool(s["rest_atoms"])); mult_spin.set_sensitive(bool(s["rest_atoms"]))
        dialog.get_content_area().add(grid); dialog.show_all()
        while True:
            response = dialog.run()
            if response != Gtk.ResponseType.OK or not s["rest_atoms"]: break
            q, m = int(charge_spin.get_value()), int(mult_spin.get_value())
            if (electrons - q) % 2 == (m - 1) % 2: break
            self.simple_dialog.error(msg="Charge {:+d} gives {} electrons: the multiplicity must be {} (1, 3, 5... for an even "
                                         "number of electrons; 2, 4, 6... for an odd one).".format(
                                         q, electrons - q, "odd" if (electrons - q) % 2 == 0 else "even"))
        use_water = solvent.get_active(); use_bonded = bonded.get_active()
        scan_torsions = scan.get_active() and use_bonded
        total_charge, multiplicity = int(charge_spin.get_value()), int(mult_spin.get_value())
        method, cm5_scale = method_combo.get_active_id(), scale_spin.get_value()
        ligand_charges = ligand_combo.get_active_id() or "cm5"
        dialog.destroy()
        if response != Gtk.ResponseType.OK: return
        if as_new:
            system, vobject = self._dyff_copy_system(system, vobject)
            if vobject is None:
                self.simple_dialog.error(msg="The copy of the system has no object."); return
            p = RL.plan(system)
        self._run_dyff_library(system, vobject, p, use_water=use_water, use_bonded=use_bonded, method=method, cm5_scale=cm5_scale,
                               scan_torsions=scan_torsions, total_charge=total_charge, multiplicity=multiplicity,
                               ligand_charges=ligand_charges, open_report=True)

    def _load_dyff_parametrization_dialog(self):
        """ [EN] 2026-10-23: System > Force Field > Define Force Field > Load DYFF Parametrization. """
        system = self.p_session.psystem.get(self.p_session.active_id)
        if system is None:
            self.simple_dialog.error(msg="No active system."); return
        if getattr(system, "qcModel", None) is not None:
            self.simple_dialog.error(msg="'{}' has a QC model: DYFF cannot be assigned.".format(system.label)); return
        chooser = Gtk.FileChooserDialog(title="Load DYFF Parametrization (parametrization.json)", parent=self.window,
                                        action=Gtk.FileChooserAction.OPEN)
        chooser.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
        flt = Gtk.FileFilter(); flt.set_name("DYFF parametrization (parametrization.json)"); flt.add_pattern("parametrization.json")
        chooser.add_filter(flt)
        folder = getattr(system, "e_working_folder", None)
        if folder and os.path.isdir(folder):
            chooser.set_current_folder(folder)
        response, path = chooser.run(), chooser.get_filename()
        chooser.destroy()
        if response != Gtk.ResponseType.OK or not path:
            return
        try:
            data = self.load_dyff_parametrization(path, system)
        except Exception as error:
            self.simple_dialog.error(msg="Load DYFF Parametrization: {}".format(error)); return
        report = os.path.join(data["_folder"], "report.html")
        if os.path.exists(report) and self.simple_dialog.question(
                "The parameters of '{}' are being applied.\n\nOpen the parametrization report?".format(system.label)):
            self._open_in_browser(report)

    def open_ai_assistant(self, *_):
        from gui.windows.setup.ai_assistant_window import AIAssistantWindow
        if getattr(self, "ai_assistant_window", None) is None:
            self.ai_assistant_window = AIAssistantWindow(main=self)
        self.ai_assistant_window.open_window()

    # ------------------------------------------------------------ remote machines
    def refresh_run_on_menu(self):
        """ Rebuilds the radio list of Simulate > Run Simulations On (menu_run_on,
            in the glade file) from the profiles of Remote Machines. """
        from remote import profiles as remote_profiles
        submenu = self.builder.get_object('menu_run_on')
        for child in submenu.get_children():
            submenu.remove(child)
        current = getattr(getattr(self, "p_session", None), "remote_target", None)
        group = None
        names = [None] + list(remote_profiles.load_all())
        if current not in names:
            current = None
        for name in names:
            label = "This Computer" if name is None else "{} (remote)".format(name)
            item = Gtk.RadioMenuItem.new_with_label_from_widget(group, label)
            group = item
            item.set_active(name == current)
            item.connect("toggled", self._on_run_on_toggled, name)
            submenu.append(item)
        submenu.show_all()

    def _on_run_on_toggled(self, item, name):
        if not item.get_active():
            return
        self.p_session.remote_target = name
        if name is None:
            self.statusbar_main.push(1, "Simulations run on this computer.")
        else:
            self.statusbar_main.push(1, "Simulations run on '{}': results come back to their output folder "
                                        "(see the Process Manager).".format(name))

    def open_remote_machines_window(self, *_):
        from gui.windows.setup.remote_machines_window import RemoteMachinesWindow
        if getattr(self, "remote_machines_window", None) is None:
            self.remote_machines_window = RemoteMachinesWindow(main=self)
        self.remote_machines_window.open_window()

    def on_menuitem_control_server_toggled(self, menuitem):
        """ [EN] 2026-10-28: Edit > AI Control Server -- starts (after the options
            dialog) or stops the control server of src/automation (docs/AI_CONTROL.md). """
        if getattr(self, "_syncing_control_menu", False):
            return
        server = getattr(self, "control_server", None)
        if menuitem.get_active():
            if server is not None and server.running:
                return
            from gui.windows.setup.ai_control_dialog import AIControlDialog, save_and_apply
            settings = AIControlDialog(parent=self.window).run()
            if settings is None:
                self.sync_control_menu(False)
                return
            status = save_and_apply(settings)
            try:
                from automation.server import ControlServer
                self.control_server = ControlServer(self, settings=settings)
                path = self.control_server.start()
            except Exception as error:
                self.sync_control_menu(False)
                self.simple_dialog.error(msg="Could not start the AI control server:\n{}".format(error))
                return
            self.statusbar_main.push(1, "AI control server on: {}. {}".format(path, status))
        else:
            if server is not None and server.running:
                server.stop()
            self.statusbar_main.push(1, "AI control server off.")

    def on_menuitem_control_settings_activate(self, menuitem):
        """ Edit > AI Control Settings...: the options, also while the server runs. """
        from gui.windows.setup.ai_control_dialog import AIControlDialog, save_and_apply
        server = getattr(self, "control_server", None)
        running = server is not None and server.running
        settings = AIControlDialog(parent=self.window, settings=server.settings if running else None,
                                   starting=False).run()
        if settings is None:
            return
        status = save_and_apply(settings)
        if running:
            server.apply_settings(settings)
        self.statusbar_main.push(1, status)

    def sync_control_menu(self, active):
        """ Shows the state of the control server in Edit > AI Control Server
            without starting or stopping anything (--control, terminal 'control'). """
        item = self.builder.get_object("menuitem_control_server")
        if item is None or item.get_active() == bool(active):
            return
        self._syncing_control_menu = True
        try:
            item.set_active(bool(active))
        finally:
            self._syncing_control_menu = False

    def open_umbrella_browser(self, folder=None):
        """ [EN] 2026-10-23: Analysis > Umbrella Sampling Browser. """
        from gui.windows.analysis.umbrella_browser_window import UmbrellaBrowserWindow
        if getattr(self, "umbrella_browser_window", None) is None:
            self.umbrella_browser_window = UmbrellaBrowserWindow(main=self)
        self.umbrella_browser_window.open_window(folder=folder)

    def open_parametrization_window(self, folder=None, system=None):
        """ [EN] 2026-10-23: System > Force Field > Parametrization Report (DYFF modified ligands). """
        from gui.windows.analysis.parametrization_window import ParametrizationWindow
        if getattr(self, "parametrization_window", None) is None:
            self.parametrization_window = ParametrizationWindow(main=self)
        self.parametrization_window.open_window(folder=folder, system=system)

    def _open_in_browser(self, path):
        """ [EN] 2026-10-23: a local HTML file in the default browser, without blocking. """
        import threading
        import webbrowser
        threading.Thread(target=lambda: webbrowser.open("file://" + os.path.abspath(path)), daemon=True).start()

    def load_dyff_parametrization(self, folder, system=None, vobject=None):
        """ [EN] 2026-10-23: applies a saved DYFF parametrization (pdynamo/param_report.py
        output folder) to the active system -- the same atoms, elements and bonds -- with the
        residue library and the stored parameters of the atoms outside it, no xTB/ORCA. """
        from pdynamo import residue_library as RL
        from pdynamo import param_report as PR
        data = PR.load(folder)
        if system is None:
            system = self.p_session.psystem[self.p_session.active_id]
        if vobject is None:
            vobject = self.get_active_vobject(getattr(system, "e_id", None))
        if vobject is None:
            raise ValueError("no object of the system '{}'".format(system.label))
        reason = PR.check_system(data, system)
        if reason:
            raise ValueError("this parametrization does not fit '{}': {}".format(system.label, reason))
        p = RL.plan(system)
        if sorted(p["rest_atoms"]) != sorted(data["rest_atoms"]):
            raise ValueError("the residue library now covers different atoms than when the parametrization was made")
        o = data.get("settings") or {}
        self._run_dyff_library(system, vobject, p, use_water=o.get("water", True), use_bonded=o.get("bonded", True),
                               method=o.get("method", "gfn2"), cm5_scale=o.get("cm5_scale", 1.0),
                               scan_torsions=o.get("scan_torsions", False), ligand_charges=o.get("ligand_charges", "cm5"),
                               total_charge=o.get("total_charge"), multiplicity=o.get("multiplicity_override"),
                               precomputed=data)
        return data

    def reapply_dyff_library(self, vobject, extra_message=""):
        """ [EN] 2026-10-08: DYFF + residue library again on the system of `vobject`, with
        the options it was parametrized with (vobject.dyff_library_options) -- after the
        structure was edited in the Builder. Residues still matching the library get the
        same parameters; the edited region (e.g. a metal site with a bonded hydroxide) is
        outside the library and is refined with xTB; metal sites come from the structure. """
        from pdynamo import residue_library as RL
        system = self.p_session.psystem.get(getattr(vobject, "e_id", None))
        if system is None:
            self.simple_dialog.error(msg="No pDynamo system for '{}'.".format(vobject.name)); return False
        if getattr(system, "qcModel", None) is not None:
            self.simple_dialog.error(msg="'{}' has a QC model: DYFF was not assigned.".format(system.label)); return False
        if len(vobject.atoms) != len(system.atoms):
            self.simple_dialog.error(msg="The object and its system have different atoms."); return False
        options = getattr(vobject, "dyff_library_options", None) or {}
        try:
            p = RL.plan(system)
        except Exception as error:
            self.simple_dialog.error(msg="Residue library: {}".format(error)); return False
        self._run_dyff_library(system, vobject, p, use_water=options.get("water", True), use_bonded=options.get("bonded", True),
                               scan_torsions=options.get("scan_torsions", False),
                               ligand_charges=options.get("ligand_charges", "cm5"),
                               total_charge=options.get("total_charge"), multiplicity=options.get("multiplicity_override"),
                               method=options.get("method", "gfn2"), cm5_scale=options.get("cm5_scale", 1.0),
                               extra_message=extra_message)
        return True

    def _dyff_copy_system(self, system, vobject=None):
        """ [EN] 2026-10-21, user's request: DYFF preparation makes a NEW system by
        default. A copy of `system` (coordinates, fixed atoms, cell; tag DYFF, name
        <label>_DYFF) is added to the session; what DYFF reads from the object --
        Builder partial charges and formal charges per atom, manual and metal type
        overrides, the previous library options -- goes to the copy's object.
        Returns ( new system, its object ). """
        import copy
        new_system = copy.deepcopy(system)
        new_system.e_step_counter = 0
        # . a system read from a file is labelled with its path: <file name>_DYFF
        label = os.path.splitext(os.path.basename(system.label or "system"))[0].replace(" ", "_") or "system"
        new_system = self.p_session.add_prepared_pdynamo_system(
            new_system, name="{}_DYFF".format(label), tag="DYFF",
            working_folder=getattr(system, "e_working_folder", None),
            input_files=dict(getattr(system, "e_input_files", None) or {}))
        new_vobject = next((v for v in self.vm_session.vm_objects_dic.values()
                            if getattr(v, "e_id", None) == new_system.e_id), None)
        if vobject is not None and new_vobject is not None and len(vobject.atoms) == len(new_vobject.atoms):
            for attribute in ("manual_atom_type_overrides", "dyff_type_overrides", "dyff_library_options"):
                value = getattr(vobject, attribute, None)
                if value is not None:
                    setattr(new_vobject, attribute, copy.deepcopy(value))
            for i, atom in vobject.atoms.items():
                target = new_vobject.atoms.get(i)
                if target is None: continue
                for attribute in ("mm_charge", "formal_charge"):
                    if getattr(atom, attribute, None) is not None:
                        setattr(target, attribute, getattr(atom, attribute))
        return new_system, new_vobject

    def _ask_dyff_charge_multiplicity(self, system, vobject, intro):
        """ [EN] 2026-10-21, user's request: charge and multiplicity of the xTB
        calculations when the dialog that started the parametrization could not ask
        (Prepare DYFF System: the hydrogens, hence the electrons, only exist after the
        build). Returns ( charge, multiplicity or None = keep the metal sites' ) or
        None when cancelled. """
        dialog = Gtk.Dialog(title="DYFF: charge and multiplicity", transient_for=self.window, modal=True)
        dialog.add_buttons("Cancel", Gtk.ResponseType.CANCEL, "Continue", Gtk.ResponseType.OK)
        dialog.set_default_response(Gtk.ResponseType.OK)
        grid = Gtk.Grid(column_spacing=8, row_spacing=8, margin=12)
        label = Gtk.Label(label=intro, xalign=0); label.set_line_wrap(True); label.set_max_width_chars(80)
        grid.attach(label, 0, 0, 2, 1)
        charge_spin, mult_spin, electrons = self._dyff_charge_multiplicity_widgets(grid, 1, system, vobject)
        default_mult = int(mult_spin.get_value())
        # . unchanged with metal sites: keep each site's own multiplicity (several metals)
        from pdynamo import dyff_prep as DP
        try:
            has_sites = bool(DP.sites_from_system(system, (getattr(vobject, "dyff_library_options", None) or {}).get("multiplicities"))[1])
        except Exception:
            has_sites = False
        dialog.get_content_area().add(grid); dialog.show_all()
        result = None
        while True:
            response = dialog.run()
            if response != Gtk.ResponseType.OK: break
            q, m = int(charge_spin.get_value()), int(mult_spin.get_value())
            if (electrons - q) % 2 == (m - 1) % 2:
                result = (q, None if (has_sites and m == default_mult) else m); break
            self.simple_dialog.error(msg="Charge {:+d} gives {} electrons: the multiplicity must be {} (1, 3, 5... for an even "
                                         "number of electrons; 2, 4, 6... for an odd one).".format(
                                         q, electrons - q, "odd" if (electrons - q) % 2 == 0 else "even"))
        dialog.destroy()
        return result

    def _dyff_charge_multiplicity_widgets(self, grid, row, system, vobject=None):
        """ [EN] 2026-10-21, user's request: "Charge" and "Multiplicity" of the xTB
        calculations (CM5 charges, Hessians, torsion scans) of the atoms outside the
        library. Defaults: the sum of the formal charges; the metal site's multiplicity
        (as kept on the object or the DYFF default for the element) or, without metals,
        1/2 from the electron count -- always made consistent with the electron count.
        Returns ( charge spin, multiplicity spin, number of electrons at charge 0 ). """
        from pdynamo import dyff_prep as DP
        electrons = sum(int(a.atomicNumber) for a in system.atoms)
        charge = sum(int(getattr(a, "formalCharge", 0) or 0) for a in system.atoms)
        previous = (getattr(vobject, "dyff_library_options", None) or {}) if vobject is not None else {}
        mult = previous.get("multiplicity_override")
        if mult is None:
            try:
                _, states = DP.sites_from_system(system, previous.get("multiplicities"))
            except Exception:
                states = {}
            mult = int(states[sorted(states)[0]]["multiplicity"]) if states else 1
        if previous.get("total_charge") is not None:
            charge = int(previous["total_charge"])
        if (electrons - charge) % 2 != (mult - 1) % 2:
            mult = max(1, mult - 1) if mult > 1 else 2
        label = Gtk.Label(label="Charge / multiplicity (xTB, atoms outside the library):", xalign=0)
        tip = ("Total charge and spin multiplicity used by the xTB calculations (CM5 charges, Hessians, torsion scans) of "
               "the atoms outside the residue library. Charge: default = sum of the formal charges. Multiplicity: 2S+1; "
               "default = the metal site's (the DYFF default for the element) or 1/2 from the electron count. With one "
               "metal it is the metal's spin state (e.g. Fe(II) 1 low spin, 5 high spin); with several metals it goes on "
               "the first one. Must agree with the number of electrons (odd for an even count).")
        label.set_tooltip_text(tip)
        grid.attach(label, 0, row, 1, 1)
        box = Gtk.Box(spacing=6)
        charge_spin = Gtk.SpinButton.new_with_range(-20, 20, 1); charge_spin.set_value(charge); charge_spin.set_tooltip_text(tip)
        mult_spin = Gtk.SpinButton.new_with_range(1, 9, 1); mult_spin.set_value(mult); mult_spin.set_tooltip_text(tip)
        box.pack_start(charge_spin, False, False, 0); box.pack_start(Gtk.Label(label="/"), False, False, 0)
        box.pack_start(mult_spin, False, False, 0)
        grid.attach(box, 1, row, 1, 1)
        return charge_spin, mult_spin, electrons

    @staticmethod
    def _dyff_needs_preparation(system):
        """ [EN] 2026-10-21 BUG FIX (user's report: a complex drawn from the structure
        library failed with "Metal-containing cofactors are not supported yet"): a system
        was sent to Prepare DYFF System (the raw-PDB route) whenever it had no double,
        triple or aromatic bond -- so every saturated molecule built in the Builder (43 of
        the 173 metal complexes of the library, alkanes, sugars) was treated as a raw PDB.
        Raw now means: no bonds, or a main-group atom whose bonds (and formal charge) do
        not complete its valence (missing hydrogens or bond orders). Atoms bonded to a
        metal are not checked (dative bonds). """
        connectivity = getattr(system, "connectivity", None)
        bonds = list(connectivity.bonds) if connectivity is not None else []
        if not bonds:
            return True
        order_of = {"Single": 1, "Double": 2, "Triple": 3}
        valence = {6: lambda q: 4 - abs(q), 7: lambda q: 3 + q, 8: lambda q: 2 + q, 9: lambda q: 1, 17: lambda q: 1,
                   35: lambda q: 1, 53: lambda q: 1, 1: lambda q: 1}
        main_group = set(valence) | {5, 14, 15, 16, 34}
        used = [0.0] * len(system.atoms)
        next_to_metal = set()
        for b in bonds:
            i, j = b.node1.index, b.node2.index
            order = 1.5 if getattr(b, "isAromatic", False) else order_of.get(getattr(getattr(b, "type", None), "name", None), 1)
            used[i] += order; used[j] += order
            zi, zj = system.atoms[i].atomicNumber, system.atoms[j].atomicNumber
            if zi not in main_group: next_to_metal.add(j)
            if zj not in main_group: next_to_metal.add(i)
        for k, atom in enumerate(system.atoms):
            rule = valence.get(atom.atomicNumber)
            if rule is None or k in next_to_metal:
                continue
            if used[k] + 0.25 < rule(int(getattr(atom, "formalCharge", 0) or 0)):
                return True
        return False

    def _dyff_scan_torsions_widget(self, bonded=None):
        """ [EN] 2026-10-21, user's request: optional xTB scans + fits of the flexible
        rotatable bonds that the DFT torsion library does not cover (residue_library.
        flexible_scan_targets / scan_flexible_torsions). Off by default: one relaxed
        xTB scan (24 points + single points) per bond. """
        check = Gtk.CheckButton(label="Scan with xTB and fit the flexible bonds not in the torsion library (slow)")
        check.set_tooltip_text("For the atoms outside the residue library: every rotatable single bond whose torsion is not "
                               "in the DFT torsion library gets a relaxed GFN2-xTB scan (24 points, on the same capped "
                               "fragment as the Hessian) and a fitted torsion. Rigid bonds are skipped: rings, double/aromatic "
                               "bonds, amide-like C(=O)-N bonds and symmetric CH3/NH2/NH3+ rotors. Each scan takes from "
                               "seconds to a few minutes; scans run in parallel. Off: those bonds keep the DYFF rule.")
        check.set_active(False)
        if bonded is not None:
            check.set_sensitive(bonded.get_active() and bonded.get_sensitive())
            bonded.connect("toggled", lambda b: (check.set_sensitive(b.get_active()), b.get_active() or check.set_active(False)))
        return check

    def _dyff_ligand_charges_widget(self):
        """ [EN] 2026-10-23, user's decision: the charge model of the ligands (the whole
        molecules among the atoms outside the residue library): CM5 (default) or RESP-2
        on the GFN2-xTB or HF/6-31G* potential (residue_library.resp_group_charges).
        Groups bonded to the rest of the system and metal sites keep CM5. """
        from pdynamo.residue_library import LIGAND_CHARGE_MODELS
        combo = Gtk.ComboBoxText()
        for key, text in LIGAND_CHARGE_MODELS:
            combo.append(key, text)
        combo.set_active_id("cm5")
        combo.set_tooltip_text("Charges of the ligands (whole molecules outside the residue library). RESP-2 is fitted to the "
                               "gas-phase electrostatic potential (System > Charges > ESP / RESP Charges, identical to AmberTools "
                               "resp). On the NCI Atlas RESP-HF was slightly better than CM5 and RESP-xTB slightly worse; on "
                               "protein-ligand interaction energies RESP-xTB was better. Termini, modified residues, covalent "
                               "ligands and metal sites keep CM5.")
        return combo

    def _dyff_method_and_scale_widgets(self, grid, row):
        """ [EN] 2026-10-08: the Hessian method (GFN2-xTB / g-xTB / DFT, the latter not yet
        available) and the CM5 charge scale, shared by the DYFF + residue library dialog
        and Prepare DYFF System. Returns ( combo, spin ). """
        grid.attach(Gtk.Label(label="Hessian method (atoms outside the library, metal sites):", xalign=0), 0, row, 1, 1)
        combo = Gtk.ComboBoxText()
        combo.append("gfn2", "GFN2-xTB (implicit water allowed)")
        combo.append("gxtb", "g-xTB (more accurate; gas phase only)")
        combo.append("dft", "DFT (ORCA) -- planned, not available yet")
        combo.set_active_id("gfn2")
        def no_dft(widget):
            if widget.get_active_id() == "dft":
                widget.set_active_id("gfn2")
                self.simple_dialog.info(msg="DFT (ORCA) refinement is planned but not implemented yet; GFN2-xTB is used.",
                                        title="DFT refinement")
        combo.connect("changed", no_dft)
        combo.set_tooltip_text("g-xTB runs through tblite, which has no implicit solvent: with charged sites prefer GFN2-xTB + ALPB.")
        grid.attach(combo, 1, row, 1, 1)
        grid.attach(Gtk.Label(label="CM5 charge scale:", xalign=0), 0, row + 1, 1, 1)
        spin = Gtk.SpinButton.new_with_range(0.5, 1.5, 0.05); spin.set_digits(2); spin.set_value(1.0)
        spin.set_tooltip_text("Charges scaled within every residue / group, totals kept: q' = q_mean + s (q - q_mean). On gas-phase DFT "
                              "dipeptide energies and forces (SPICE) 0.7-0.8 is best; condensed-phase work has not been benchmarked: "
                              "keep 1.0 unless you test it.")
        grid.attach(spin, 1, row + 1, 1, 1)
        return combo, spin

    def _run_dyff_library(self, system, vobject, p, use_water=True, use_bonded=True, method="gfn2", site_states=None,
                          cm5_scale=1.0, extra_message="", on_success=None, scan_torsions=False,
                          total_charge=None, multiplicity=None, ligand_charges="cm5", open_report=False, precomputed=None):
        """ [EN] 2026-10-08: applies DYFF + residue library to `system` (plan p from
        residue_library.plan): the atoms outside the library refined on capped subsystems
        in a worker thread (method gfn2/gxtb; site_states = spin of metal sites), the
        charges optionally scaled, then the MM model defined.
        2026-10-23: every run with atoms outside the library writes an output folder
        (pdynamo/param_report.py: parametrization.json, report.html, the original xTB/ORCA
        logs) in the system's working folder; open_report: show report.html when done.
        precomputed: a loaded parametrization.json (param_report.load) whose parameters are
        applied instead of computing them again. """
        import threading
        import time
        from pdynamo import param_report as PR
        from pdynamo import residue_library as RL
        from pdynamo import dyff_prep as DP
        s = p["summary"]
        if method == "dft": method = "gfn2"
        previous = getattr(vobject, "dyff_library_options", None) or {}
        if site_states is None:
            # . 2026-10-08: metal sites read from the structure itself (a system edited in the
            #   Builder after Prepare DYFF: hydroxide bonded to Zn, donors changed, ...); the
            #   spin states chosen at preparation are kept by atom path
            overrides, site_states = DP.sites_from_system(system, previous.get("multiplicities"))
            if overrides:
                vobject.dyff_type_overrides = overrides
        if multiplicity is not None and site_states:
            # . 2026-10-21: the multiplicity chosen in the dialog goes on the metal site
            #   (several metals: the first one, the others singlet)
            for k, m in enumerate(sorted(site_states)):
                site_states[m]["multiplicity"] = int(multiplicity) if k == 0 else 1
        # . kept on the object: Edit in Builder carries them, and closing the Builder can
        #   re-apply the library with the same choices (reapply_dyff_library)
        vobject.dyff_library_options = {"water": use_water, "bonded": use_bonded, "method": method, "cm5_scale": cm5_scale,
                                        "scan_torsions": bool(scan_torsions and use_bonded),
                                        "ligand_charges": ligand_charges,
                                        "total_charge": total_charge, "multiplicity_override": multiplicity,
                                        "multiplicities": {system.atoms[i].path: st["multiplicity"] for i, st in (site_states or {}).items()}}
        self.statusbar_main.push(1, "DYFF + residue library: {} atoms outside the library...".format(s["rest_atoms"]))
        def say(text):
            GLib.idle_add(lambda: self.statusbar_main.push(1, "DYFF + residue library: " + text) and False)
        folder = None
        if p["rest_atoms"] and precomputed is None:
            try:
                folder = PR.new_folder(getattr(system, "e_working_folder", None) or os.getcwd(), system.label)
            except Exception as error:
                print("DYFF parametrization report: no output folder ({})".format(error))
        started = time.time()
        def work():
            if precomputed is not None:
                GLib.idle_add(done, PR.rest_result(precomputed) if p["rest_atoms"] else None, None); return
            try:
                rest = RL.rest_parameters(system, p["rest_atoms"], solvent="water" if use_water else None, log_dir=os.path.join(folder, "logs") if folder else None,
                                          hessian=use_bonded, torsions=use_bonded, say=say, method=method,
                                          site_states=site_states, scan_torsions=bool(scan_torsions and use_bonded),
                                          known_charges=p["charges"], known_edits=p["edits"],
                                          total_charge=total_charge, multiplicity=multiplicity,
                                          ligand_charges=ligand_charges) if p["rest_atoms"] else None
                GLib.idle_add(done, rest, None)
            except Exception as error:
                traceback.print_exc()
                GLib.idle_add(done, None, str(error))
        def done(rest, error):
            if error:
                self.simple_dialog.error(msg="Parameters for the atoms outside the library failed:\n{}".format(error[-600:])); return False
            charges = dict(p["charges"]); edits = {k: dict(v) for k, v in p["edits"].items()}
            if rest:
                charges.update(rest["charges"])
                for key in ("bond_terms", "angle_terms", "dihedral_bonds"):
                    edits.setdefault(key, {}).update(rest["edits"][key])
            if abs(cm5_scale - 1.0) > 1e-9:
                symbols, neighbours = RL.system_graph(system)
                charges = RL.scale_charges(charges, RL.charge_groups(p["matches"], p["rest_atoms"], neighbours), cm5_scale)
            total = RL.store_on_object(vobject, system, edits, charges)
            ok, msg = self.p_session.define_MMModel(force_field='DYFF', system=system)
            if ok and on_success is not None:
                # . 2026-10-17: e.g. the Builder's Force Field Settings (nonbonded model)
                msg += on_success() or ""
            self.refresh_main_statusbar()
            if ok:
                msg += " Residue library: {}/{} residues; total charge {:+.3f}.".format(s["matched"], s["residues"], total)
                if rest and use_bonded:
                    r = rest["report"]
                    msg += " Outside the library ({}): {} bonds, {} angles from the Hessian; {} of {} rotatable bonds from the DFT torsion library.".format(
                        {"gxtb": "g-xTB"}.get(method, "GFN2-xTB"), r["bonds"], r["angles"], r["torsions"], r["rotatable"])
                    if "scan_targets" in r:
                        # . 2026-10-21: automatic xTB scans of the flexible bonds the library does not have
                        msg += " xTB torsion scans: {} of {} flexible bond(s) fitted".format(r.get("scanned", 0), r["scan_targets"])
                        details = [d for d in (r.get("scan_details") or []) if d.get("accepted", True)
                                   and d.get("fit_rmse") is not None and d.get("before_rmse") is not None]
                        if details:
                            msg += " (RMSE vs xTB, mean: {:.2f} -> {:.2f} kcal/mol)".format(
                                sum(d["before_rmse"] for d in details) / len(details), sum(d["fit_rmse"] for d in details) / len(details))
                            poor = [d for d in details if d.get("poor")]
                            if poor:
                                msg += " WARNING: {} fitted torsion(s) still above 2 kcal/mol RMSE (bonds {}): the scan profile is not a pure torsion (internal hydrogen bonds?); check them in Force Field Parameters.".format(
                                    len(poor), ", ".join(d["bond"] for d in poor[:5]))
                        failed = r.get("scan_failed") or []
                        msg += (" ({} failed: {}).".format(len(failed), "; ".join(failed[:3])) if failed else ".")
                if rest and rest["report"].get("ligand_charges"):
                    # . 2026-10-23: RESP charges of the ligands
                    groups = rest["report"]["ligand_charges"]
                    resp = [g for g in groups if g["status"].startswith("RESP")]
                    kept = [g for g in groups if g["status"].startswith("CM5")]
                    msg += " Ligand charges: {} molecule(s) with {}".format(len(resp), resp[0]["status"] if resp else "RESP")
                    msg += ("; {} group(s) kept CM5 ({}).".format(len(kept), "; ".join(sorted({g["status"][10:] for g in kept})))
                            if kept else ".")
                if abs(cm5_scale - 1.0) > 1e-9:
                    msg += " CM5 charges scaled by {:.2f}.".format(cm5_scale)
                msg += extra_message
                if precomputed is not None:
                    vobject.dyff_param_folder = precomputed.get("_folder")
                    msg += " Parameters read from {} (no new calculation).".format(precomputed.get("_folder", "parametrization.json"))
                elif folder:
                    try:
                        data = PR.collect(system, p, rest, dict(vobject.dyff_library_options), charges, edits,
                                          messages=[msg], total_seconds=time.time() - started)
                        path = PR.write(folder, data)
                        vobject.dyff_param_folder = folder
                        msg += " Report and logs: {}".format(path)
                        if open_report:
                            # . 2026-10-23: the DYFF Parametrization window (report.html is written as well)
                            GLib.idle_add(lambda: self.open_parametrization_window(folder, system) and False)
                    except Exception as error:
                        traceback.print_exc()
                        msg += " (the parametrization report could not be written: {})".format(error)
                self.bottom_notebook.status_teeview_add_new_item(message=msg, system=system)
                self.statusbar_main.push(1, msg)
            else:
                self.simple_dialog.error(msg=msg)
            return False
        threading.Thread(target=work, daemon=True).start()

    def refresh_main_statusbar(self, message = None, psystem = None):
        """ Function doc """
        psystem = self.p_session.psystem[self.p_session.active_id]
            
        if message:
            string = message
            #print(string)
            #self.statusbar_main.push(1,string)
        else:

            #string = 'System: {}  Size: {}  '.format()
            name    = psystem.label
            size    = len(psystem.atoms)
            string = 'system: {}    atoms: {}    '.format(name, size)

            if psystem.qcModel:
                hamiltonian   = getattr(psystem.qcModel, 'hamiltonian', False)
                if hamiltonian:
                    pass
                else:
                    try:
                        itens = psystem.qcModel.SummaryItems()
                        #print(itens)
                        hamiltonian = itens[0][0]
                    except:
                        hamiltonian = 'external'
                    
                    
                    
                
                n_QC_atoms    = len(list(psystem.qcState.pureQCAtoms))
                
                
                summary_items = psystem.electronicState.SummaryItems()
                
                string += 'hamiltonian: {}    QC atoms: {}    QC charge: {}    spin multiplicity {}    '.format(  hamiltonian, 
                                                                                                                  n_QC_atoms,
                                                                                                                  summary_items[1][1],
                                                                                                                  summary_items[2][1],
                                                                                                                 )
                
            n_fixed_atoms = len(psystem.e_fixed_table)
            string += 'fixed atoms: {}    '.format(n_fixed_atoms)
            
            if psystem.mmModel:
                forceField = psystem.mmModel.forceField
                string += 'force field: {}    '.format(forceField)
            
                if psystem.nbModel:
                    nbmodel = psystem.mmModel.forceField
                    string += 'nbModel: True    '
                    
                    summary_items = psystem.nbModel.SummaryItems()
                    
                
                else:
                    string += 'nbModel: False    '
                
                
                if psystem.symmetry:
                    #nbmodel = psystem.mmModel.forceField
                    string += 'PBC: True    symmetry: {}    '.format( psystem.symmetry.crystalSystem.label)
                    #print(psystem.symmetry)
                    #print(psystem.symmetryParameters)
                    #summary_items = psystem.nbModel.SummaryItems()
                    
                
                else:
                    string += 'PBC: False    '
            
            
            
            
            
            '''
            color_palette = self.p_session.get_color_palette()
            color     = color_palette['C']
            res_color = [int(color[0]*255),int(color[1]*255),int(color[2]*255)] 
            #swatch    = self.getColouredPixmap( res_color[0], res_color[1], res_color[2] )
            
            r, g, b, a = int(color[0]*255), int(color[1]*255), int(color[2]*255), 255 
            
            
            CHANNEL_BITS=8
            WIDTH=10
            HEIGHT=10
            swatch = GdkPixbuf.Pixbuf.new( GdkPixbuf.Colorspace.RGB, True, CHANNEL_BITS, WIDTH, HEIGHT ) 
            swatch.fill( (r<<24) | (g<<16) | (b<<8) | a ) # RGBA
            
            self.statusbar_main.push(0, string, swatch)
            '''
        self.statusbar_main.push(1,string)

    def update_hover_statusbar(self, atom = None):
        """ Shows the atom under the mouse cursor in the hover status bar,
        as  Resname(Resnumber)  Symbol/AtomName/AtomType (Id)  -- AtomType is
        the MM force-field type, left out when the system has no MM model.
        Called by VismolGLCore whenever the hovered atom changes; None
        clears it. """
        self.statusbar_hover.remove_all(self.statusbar_hover_context)
        if atom is None:
            return
        
        vobject = atom.vm_object
        index   = atom.atom_id
        fields  = [atom.symbol, atom.name]
        
        psystem = self.p_session.psystem.get(vobject.e_id) if vobject.e_id is not None else None
        if psystem is not None and len(psystem.atoms) == len(vobject.atoms):
            mmState = getattr(psystem, 'mmState', None)
            if mmState is not None and getattr(mmState, 'atomTypes', None) is not None:
                try:
                    fields.append(str(mmState.atomTypes[mmState.atomTypeIndices[index]]))
                except (IndexError, TypeError):
                    pass
        
        string = '{} ({})'.format('/'.join(fields), atom.index if atom.index is not None else index+1)
        if atom.residue is not None:
            string = '{}({})  {}'.format(atom.residue.name, atom.residue.index, string)
        self.statusbar_hover.push(self.statusbar_hover_context, string)

    def refresh_widgets (self, statusbar = True):
        """ Function doc """
        #if self.main_treeview:
        #    self.main_treeview.refresh_trajectory_scalebar()
        if statusbar:
            self.refresh_main_statusbar()

    def print_btn (self, widget):
        """ Function doc """
        #print(self.box_reac)
        parm = self.box_reac.get_rc_data()

    def run_menu_item_test (self, widget):
        """ Function doc """
        pklist = self.vm_session.picking_selections.picking_selections_list
        vobj =  pklist[0].vm_object
        #print(vobj.residues)
        #print(vobj.chains)
        
        #print(self.vm_session.picking_selections.picking_selections_list)
        pklist = self.vm_session.picking_selections.picking_selections_list
        if pklist[0]  and pklist[1]:
            #print(pklist[0], pklist[1])
            vobject = pklist[0].vm_object
            size = len(vobject.frames)
            
            dist1 = None
            dist2 = None
            dist3 = None
            
            for i in range(size):

                a1_coord = pklist[0].coords(frame=i)
                a2_coord = pklist[1].coords(frame=i)
                
                dist1 = get_simple_distance(a1_coord, a2_coord)
                
                
                if pklist[1] and pklist[2]:
                    a1_coord = pklist[1].coords(frame=i)
                    a2_coord = pklist[2].coords(frame=i)
                    dist2 = get_simple_distance(a1_coord, a2_coord)
                
                if pklist[2] and pklist[3]:
                    a1_coord = pklist[2].coords(frame=i)
                    a2_coord = pklist[3].coords(frame=i)
                    dist3 = get_simple_distance(a1_coord, a2_coord)
                dprint(dist1, dist2, dist3)
                
    def run_test (self, widget):
        '''Test'''
        pass
        
    def save_glArea_image (self, widget):
        """ Opens the Preview/Export window: a live (screen-resolution)
            preview of the current OpenGL view where the user can toggle
            and tweak the cartoon filter, save/load named style presets,
            and export the final PNG at a chosen resolution multiplier.
            The export re-captures the scene at full quality through
            VismolGLCore.render_to_image() rather than upscaling the
            preview image - see VismolGTKWidget.open_export_preview().
        """
        glwidget = self.vm_session.vm_glcore.parent_widget
        glwidget.open_export_preview()

    def opengl_to_povray_camera(self,  view_matrix, projection_matrix, camera_position, use_clipping=True):
        """
        Convert OpenGL matrices to a POV-Ray camera block.

        Args:
            view_matrix (np.ndarray): 4x4 OpenGL view matrix.
            projection_matrix (np.ndarray): 4x4 OpenGL projection matrix.
            camera_position (list/np.ndarray): Camera position [x, y, z].
            use_clipping (bool): Se True, adiciona o clipping (near/far).

        Returns:
            str: POV-Ray camera block.
        """
        import math
        # Vertical field of view (in degrees)
        fov_y_rad = 2 * math.atan(1 / projection_matrix[1,1])
        fov_y_deg = math.degrees(fov_y_rad)
        
        # near/far planes from the projection matrix
        if use_clipping:
            p22 = projection_matrix[2,2]
            p32 = projection_matrix[3,2]
            # resolver f e n: 
            # p22 = -(f+n)/(f-n), p32 = -2*f*n/(f-n)
            denom = p22 - 1e-12  # avoid division by zero
            far = (-p32)/(p22 + 1)
            near = far * (p22 - 1)/(p22 + 1)
        else:
            near, far = None, None

        # POV-Ray default look_at (origin)
        look_at = np.array([0,0,0], dtype=float)
        
        # Assemble the camera block
        pov_lines = []
        pov_lines.append("camera {")
        pov_lines.append(f"    location <{camera_position[0]}, {camera_position[1]}, {camera_position[2]}>")
        pov_lines.append(f"    look_at <{look_at[0]}, {look_at[1]}, {look_at[2]}>")
        pov_lines.append(f"    angle {fov_y_deg:.4f}")
        
        if use_clipping and near is not None and far is not None:
            pov_lines.append(f"    clipping {{ {near:.4f}, {far:.4f} }}")
        
        pov_lines.append("}")
        
        return "\n".join(pov_lines)

    def _persist_window_and_config (self):
        """ Salva gl_parameters (incluindo as dimensoes atuais da janela,
            atualizadas a cada resize em window_resize()) em disco antes de
            fechar o EasyHybrid, para restaurar na proxima abertura (ver
            MainWindow.__init__). Silencioso em caso de falha -- nao deve
            impedir o encerramento do programa. """
        try:
            self.vm_session.vm_config.save_easyhybrid_config()
        except Exception as e:
            logger.warning('Could not save preferences on exit: %s', e)

    def on_delete_event(self, widget, event):
        if self.p_session.changed:
            if self.session_filename == None:
                filename = 'untitled'
                #dialog.set_text('''Save changes to the project "{}" before closing?\n Your changes will be lost if you don't save them'''.format(filename))
            else:
                #dialog.set_text('''Save changes to the project "{}" before closing?\n Your changes will be lost if you don't save them'''.format(self.session_filename))
                filename = os.path.split(self.session_filename)
                filename = filename[-1]
            
            
            dialog = Gtk.MessageDialog(
                transient_for=self.window,
                modal=True,
                message_type=Gtk.MessageType.QUESTION,
                text='''Save changes to the EasyHybrid session "{}" before closing?\n Your changes will be lost if you don't save them'''.format(filename)
            )
            
            # .Save changes to the project "Untitled" before closing?
            # .Your changes will be lost if you don' tsave them
            
            dialog.add_button("Don't Save", Gtk.ResponseType.NO)
            dialog.add_button("Cancel", Gtk.ResponseType.CANCEL)
            dialog.add_button("Save", Gtk.ResponseType.YES)
            response = dialog.run()
            dialog.destroy()
            #print (response) 
            

            
            if response == -9:   # Don't Save and close
                self._persist_window_and_config()
                Gtk.main_quit()
            
            elif response == -6: # Cancel - don't close
                return True 
            
            elif response == -8: # Save and close
                if self.session_filename == None:
                    self.gtk_save_as_file(widget)
                else:
                    self.p_session.save_easyhybrid_session( filename = self.session_filename)
                #print('Saving' )
                self._persist_window_and_config()
                Gtk.main_quit() 
            
            else:
                return True 
        else:
            self._persist_window_and_config()
            Gtk.main_quit()
