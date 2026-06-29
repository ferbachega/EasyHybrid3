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
from gi.repository import GLib


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

# --- imports entre modulos adicionados na refatoracao ---
from gui.main.main_treeview import EasyHybridMainTreeView
from gui.main.treeview_menu import TreeViewMenu
from gui.main.bottom_notebook import BottomNoteBook
from gui.main.preferences_window import PreferencesWindow

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
            int                # 8: step counter
        )

        # -------------------- GTK BUILDER --------------------
        # Load Glade UI and connect signals
        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(self.home,'src/gui/MainWindow.glade'))
        self.builder.connect_signals(self)
        self.window = self.builder.get_object('window1')
        self.window.set_default_size(1200, 600)                          
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
        
        self.paned_V         = self.builder.get_object('paned_V')
        
        # -------------------- SESSION MANAGEMENT --------------------
        self.vm_session      = vm_session#( main = None)
        self.vm_session.main = self
        
        self.vm_session.vm_object_counter = 0
        self.vm_session.insert_glmenu()
        
        # Connect key press/release events for 3D visualization widget
        self.window.connect("key-press-event",   self.vm_session.vm_widget.key_pressed)
        self.window.connect("key-release-event", self.vm_session.vm_widget.key_released)
        
                
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

        self.notebook_H1.append_page(child = self.ScrolledWindow_notebook_H1, tab_label= Gtk.Label(label = 'Objects'))
        


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
            self.paned_V.add(self.paned_H)
            
            # Bottom notebook for logs, status, and command terminal
            self.bottom_notebook = BottomNoteBook(main = self)
            self.paned_V.add(self.bottom_notebook.widget)
            self.paned_V_position = 400
            self.paned_V.set_position(self.paned_V_position)

            self.bottom_notebook.status_teeview_add_new_item(message = 'Welcome to EasyHybrid 3.0, have a happy simulation day!')

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
        

        
        '''#- - - - - - - - - - G T K  W I N D O W S - - - - - - - - - - -#'''
        self.window_list = []
        
        self.filechooser                  = FileChooser(main_window = self)
        self.NewSystemWindow              = ImportANewSystemWindow       ( main = self )
        
        self.setup_QCModel_window         = EasyHybridSetupQCModelWindow ( main = self )
        self.window_list.append(self.setup_QCModel_window)
        
        self.geometry_optimization_window = GeometryOptimization  ( main = self )
        self.window_list.append(self.geometry_optimization_window)
        
        self.PES_scan_window              = PotentialEnergyScanWindow    ( main=  self)
        self.window_list.append(self.PES_scan_window)
        
        self.PES_advanced_scan_window              = AdvancedPotentialEnergyScanWindow ( main=  self)
        self.window_list.append(self.PES_advanced_scan_window)

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
        
        self.align_trajectory_window = AlignTrajectoryWindow (main = self, system_liststore = self.system_liststore)
        
        self.reimaging_trajectory_window = ReimagingTrajectoryWindow (main = self, system_liststore = self.system_liststore)
        
        
        self.normal_modes_window          =   NormalModesWindow (main = self)
        self.window_list.append(self.normal_modes_window)

        
        self.WHAM_window   =  WHAMWindow(main = self)
        self.window_list.append(self.WHAM_window)

        self.simple_dialog =  SimpleDialog(main = self)

        self.edit_frames_dialog = EditFrameDialog(main = self)
        self.merge_system_window = MergeSystemWindow(main = self)
        self.solvate_system_window = SolvateSystemWindow(main = self)
        self.preferences_window = EasyHybridPreferencesWindow(main = self)

        self.make_solvent_box_window = MakeSolventBoxWindow(main = self)
        '''#- - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -#'''

        # -------------------- WINDOW SIGNALS --------------------
        #self.window.connect("destroy", Gtk.main_quit)
        #self.window.connect("delete-event", Gtk.main_quit)
        #self.window.connect("destroy",      self.quit_easyhybrid)
        self.window.connect("delete-event", self.on_delete_event)
        self.window.connect("check-resize", self.window_resize)       
        self.window.show_all()


    def on_drag_data_received(self, widget, drag_context,
                              x, y, data, info, time):

        for uri in data.get_uris():
            path, _ = GLib.filename_from_uri(uri)
            print(path)
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
        w, h = a.get_size()
    
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

    def on_main_toolbar_clicked (self, button):
        """ Function doc """
        if button  == self.builder.get_object('toolbutton_new_system'):
            self.NewSystemWindow.open_window()
        
        if button  == self.builder.get_object('toolbutton_save'):
            
            if self.session_filename == None:
                self.gtk_save_as_file(button)
                
            else:
                self.p_session.save_easyhybrid_session( filename = self.session_filename)
        
        if button  == self.builder.get_object('toolbutton_save_as'):
            self.gtk_save_as_file (button)
        
        if button  == self.builder.get_object('_show_cell'):
            print('Under construction!')
            #self.run_test(None)
            
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
            self.setup_QCModel_window.open_window()
        
        if button  == self.builder.get_object('toolbutton_system_check'): 
            self.p_session.systems[self.p_session.active_id]['vobject'].get_backbone_indexes ()

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
        
        if button  == self.builder.get_object('button_test'):
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
            e_id = self.p_session.active_id
            v_id = -1
            self.preferences = PreferencesWindow(main = self , 
                                         e_id = e_id     ,
                                         v_id = v_id     )
        
        
        elif menuitem == self.builder.get_object('menuitem_qc_setup'):
            self.setup_QCModel_window.open_window()
            
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



        elif menuitem == self.builder.get_object('menuitem_charge_inspection'): 
            TrueFalse, msg = self.p_session.check_for_fragmented_charges()
            #if true_or_false:
            if TrueFalse:
                self.simple_dialog.info(msg = msg )
            else:
                self.simple_dialog.error(msg = msg )
        
        
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
            
        elif menuitem == self.builder.get_object('menuitem_molecular_dynamics'):
            self.molecular_dynamics_window.open_window()
            
        elif menuitem == self.builder.get_object('menuitem_normal_modes'):
            self.normal_modes_window.open_window()
            
        elif menuitem == self.builder.get_object('menuitem_rection_coordinate_scans'):
            self.PES_scan_window.open_window()
        
        elif menuitem == self.builder.get_object('menuitem_advanced_rc_scans'):
            self.PES_advanced_scan_window.open_window()
            
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
        
        elif menuitem == self.builder.get_object('menuitem_surface_analysis'):
            self.surface_analysis_window.open_window()
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
            rama = RamachandranWindow()
        
        elif menuitem == self.builder.get_object('menuitem_RMSD_tool'):
            #self.rmsd_tool_window.open_window()
            self.rmsd_analysis_window.open_window()
        
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
        
        
        
        elif menuitem == self.builder.get_object('menuitem_about'):
            builder = Gtk.Builder()
            builder.add_from_file(os.path.join(self.home,'src/gui/windows/setup/about_dialog.glade'))

            dialog = builder.get_object('about_dialog')
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
            print("File selected: " + dialog.get_filename())
            self.p_session.save_easyhybrid_session( filename = dialog.get_filename())
            dialog.destroy()
            # Save file here...
        elif response == Gtk.ResponseType.CANCEL:
            print("Save operation canceled.")
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
        #print(name, v_id,  e_id)
        #print(name)
        if v_id == -1: #.change the header
            _iter = self.system_treeview_iters[e_id] 
            #_iter = self.p_session.psystem[e_id].e_treeview_iter
            self.main_treeview.treestore[_iter][2] = str(e_id)+' - '+ name
            self.p_session.psystem[e_id].label  = name
            
            liststore_iter = self.system_liststore_iters[e_id]
            #liststore_iter = self.p_session.psystem[e_id].e_liststore_iter
            self.system_liststore[liststore_iter][0] = str(e_id)+' - '+ name
  
        else:
            #print(self.vm_session.vobject_names.values())
            if name in self.vm_session.vobject_names.keys():
                print('Invalid name.')
                return False

            else:
                _iter = self.vm_session.vm_objects_dic[v_id].e_treeview_iter
                self.main_treeview.treestore[_iter][2] = name
                
                old_name = self.vm_session.vm_objects_dic[v_id].name          
                self.vm_session.vobject_names.pop(old_name)
                

                self.vm_session.vm_objects_dic[v_id].name = name
                
                self.vm_session.vobject_names[name] = self.vm_session.vm_objects_dic[v_id]
                #print('aqui')
                try:
                    self.vobject_liststore_dict[e_id][self.vm_session.vm_objects_dic[v_id].liststore_iter][0] = name
                except:
                    #means that it is surface 
                    pass
    
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
            
            '''organizing the list of vobjects that should be removed from vm_object_dic'''
            pop_list = []
            for index, vobject in self.vm_session.vm_objects_dic.items():
                if vobject.e_id == system_e_id:                   
                    self.main_treeview.treestore.remove(vobject.e_treeview_iter)
                    pop_list.append(index)
            
            '''removing from vm_object_dic'''
            for index in pop_list:
                print('removing sys = {}, obj = {}, name = {}'.format(system_e_id, 
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
                print ('deleting plotting data...')
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

            self.vm_session.vm_glcore.queue_draw()
        
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
                print(dist1, dist2, dist3)
                
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
        Converte matrizes OpenGL para um bloco de câmera POV-Ray.

        Args:
            view_matrix (np.ndarray): 4x4 matriz de visualização OpenGL.
            projection_matrix (np.ndarray): 4x4 matriz de projeção OpenGL.
            camera_position (list/np.ndarray): Posição da câmera [x, y, z].
            use_clipping (bool): Se True, adiciona o clipping (near/far).

        Returns:
            str: Bloco de câmera POV-Ray.
        """
        import math
        # Campo de visão vertical (em graus)
        fov_y_rad = 2 * math.atan(1 / projection_matrix[1,1])
        fov_y_deg = math.degrees(fov_y_rad)
        
        # Planos near/far da matriz de projeção
        if use_clipping:
            p22 = projection_matrix[2,2]
            p32 = projection_matrix[3,2]
            # resolver f e n: 
            # p22 = -(f+n)/(f-n), p32 = -2*f*n/(f-n)
            denom = p22 - 1e-12  # evitar divisão por zero
            far = (-p32)/(p22 + 1)
            near = far * (p22 - 1)/(p22 + 1)
        else:
            near, far = None, None

        # POV-Ray look_at padrão (origem)
        look_at = np.array([0,0,0], dtype=float)
        
        # Monta o bloco da câmera
        pov_lines = []
        pov_lines.append("camera {")
        pov_lines.append(f"    location <{camera_position[0]}, {camera_position[1]}, {camera_position[2]}>")
        pov_lines.append(f"    look_at <{look_at[0]}, {look_at[1]}, {look_at[2]}>")
        pov_lines.append(f"    angle {fov_y_deg:.4f}")
        
        if use_clipping and near is not None and far is not None:
            pov_lines.append(f"    clipping {{ {near:.4f}, {far:.4f} }}")
        
        pov_lines.append("}")
        
        return "\n".join(pov_lines)

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
                Gtk.main_quit()
            
            elif response == -6: # Cancel - don't close
                return True 
            
            elif response == -8: # Save and close
                if self.session_filename == None:
                    self.gtk_save_as_file(widget)
                else:
                    self.p_session.save_easyhybrid_session( filename = self.session_filename)
                #print('Saving' )
                Gtk.main_quit() 
            
            else:
                return True 
        else:
            Gtk.main_quit()
