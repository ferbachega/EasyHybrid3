#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  main.py
#  
#  Copyright 2023 Fernando <fernando@fernando-Victus-by-HP-Gaming-Laptop-15-fb1xxx>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
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
#  
import os, sys, time
import gi 
import signal
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
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

from gui.windows.setup.easyhybrid_terminal    import TerminalWindow
from gui.windows.setup.selection_list_window  import *
from gui.windows.setup.setup_interface        import EasyHybridPreferencesWindow

from gui.windows.simulation.single_point_window          import SinglePointWindow
from gui.windows.simulation.geometry_optimization_window import GeometryOptimization
from gui.windows.simulation.PES_scan_window              import PotentialEnergyScanWindow 
from gui.windows.simulation.molecular_dynamics_window    import MolecularDynamicsWindow 
from gui.windows.simulation.umbrella_sampling_window     import UmbrellaSamplingWindow 
from gui.windows.simulation.chain_of_states_opt_window   import ChainOfStatesOptWindow 
from gui.windows.simulation.normal_modes_window          import NormalModesWindow 

from gui.windows.analysis.WHAM_analysis_window           import WHAMWindow 
from gui.windows.analysis.normal_modes_analysis_window   import NormalModesAnalysisWindow 
from gui.windows.analysis.surface_analysis_window        import SurfaceAnalysisWindow 
from gui.windows.analysis.energy_refinement_window       import EnergyRefinementWindow
from gui.windows.analysis.PES_analysis_window            import PotentialEnergyAnalysisWindow

from util.geometric_analysis import get_simple_distance
from util.sequence_plot import GtkSequenceViewer


from pdynamo.pDynamo2EasyHybrid import pDynamoSession
import numpy as np



from pCore                     import Align                                        , \
                                      Clone                                        , \
                                      logFile                                      , \
                                      Selection                                    , \
                                      TestScript_InputDataPath                     , \
                                      TestScript_OutputDataPath                    , \
                                      XHTMLLogFileWriter







class MainWindow:
    """ Class doc """
    def __init__ (self, vm_session = None, home = None, version = None):
        """ Class initialiser """
        self.home               =  home
        self.EASYHYBRID_VERSION =  version
        #signal.signal(signal.SIGINT, self.on_sigint_received)
        #print (self.home, self.EASYHYBRID_VERSION)

        '''Search home folder. Every time a file is loaded into memory 
        this attribute will be modified.'''
        self.current_search_folder = None

        


        '''#- - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -#'''
        self.system_liststore        = Gtk.ListStore(str, int, GdkPixbuf.Pixbuf)

        
        '''The "vobject_liststore_dict" is a dictionary where the access key is
        the e_id, which is the index of the system of interest generated in 
        "pDynamo2Easyhybrid/pDynamoSession/append_system_to_pdynamo_session". 
        Each dictionary element contains a liststore that includes the respective 
        vobjects.'''
        self.vobject_liststore_dict  = {               
                                       0 : Gtk.ListStore(str,  int, int, GdkPixbuf.Pixbuf)  # name, object_index, e_id, pixel_buffer
                                       }                                 
        '''#- - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -#'''



        
        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(self.home,'src/gui/MainWindow.glade'))
        self.builder.connect_signals(self)
        self.window = self.builder.get_object('window1')
        self.window.set_default_size(1200, 600)                          
        self.window.set_title('EasyHybrid {}'.format(self.EASYHYBRID_VERSION))                          


        self.statusbar_main = self.builder.get_object('statusbar1')
        self.statusbar_main.push(1,'Welcome to EasyHybrid version {}, a pDynamo3 graphical tool'.format(self.EASYHYBRID_VERSION))
        
        self.paned_V         = self.builder.get_object('paned_V')
        self.vm_session      = vm_session#( main = None)
        self.vm_session.main = self
        
        self.vm_session.vm_object_counter = 0
        self.vm_session.insert_glmenu()
        
        self.window.connect("key-press-event",   self.vm_session.vm_widget.key_pressed)
        self.window.connect("key-release-event", self.vm_session.vm_widget.key_released)
        
                
        
        self.menu_box = self.builder.get_object('toolbutton_selection_box')
        self.box2 = self.builder.get_object('box2')
        
        #self.selection_box = self.vm_session.selection_box
        self.selection_box_frame = VismolSelectionTypeBox(vm_session = self.vm_session)
        self.vm_session.selection_box_frame = self.selection_box_frame 
        #print('\n\n\n',self.selection_box_frame,self.vm_session.selection_box_frame ,'\n\n')
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
            #player

            self.container.pack_start(self.vm_session.vm_widget, True, True, 0)
            
            #self.traj_frame = self.vm_session.trajectory_frame
            #self.container.pack_start(self.traj_frame, False, False, 1)
            #self.container.pack_start(self.command_line_entry, False, False, 0)

            self.notebook_H2.append_page(child = self.container, tab_label = Gtk.Label(label = 'view'))
            #self.notebook_H2.append_page(Gtk.TextView(), Gtk.Label('logs'))
            
            
            #self.HBOX = Gtk.Box(orientation = Gtk.Orientation.VERTICAL, spacing = 6)
            self.HBOX = Gtk.Box(orientation = Gtk.Orientation.VERTICAL, spacing = 0)
            self.HBOX.pack_start(self.notebook_H1, True, True, 0)
            #self.label_frame = Gtk.Label('Frame Number: 0')
            #self.HBOX.pack_start(self.label_frame, False, False, 0)
            #self.HBOX.pack_start(self.traj_frame, False, False, 1)

            #self.paned_H.add(self.notebook_H1)
            self.paned_H.add(self.HBOX)
            self.paned_H.add(self.notebook_H2)
            self.paned_H.set_position(250)

            self.paned_V.add(self.paned_H)
            
            self.bottom_notebook = BottomNoteBook(main = self)
            #self.paned_V.add(self.builder.get_object('notebook_text_and_logs'))
            self.paned_V.add(self.bottom_notebook.widget)
            self.paned_V_position = 400
            self.paned_V.set_position(self.paned_V_position)

            self.bottom_notebook.status_teeview_add_new_item(message = 'Welcome to EasyHybrid 3.0, have a happy simulation day!')
            
            #for i in range (10):
            #    self.bottom_notebook.status_teeview_add_new_items([str(i),str(i),str(i)])
            
            #self.paned_V.add(Gtk.TextView())
            
            #self.paned_V.add(self.traj_frame)
        ''' - - - - - - - - - - - - Terminal e Trajectory Play - - - - - - - - - - - - - - - - '''
        self.trajectory_player_button = self.builder.get_object('toolbutton_trajectory_tool1')
        self.cmd_terminal_button      = self.builder.get_object('toolbutton_terminal')
        #self.player_frame = self.vm_session.player_frame
        #self.player_frame.show_all()
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

        self.selection_list_window        = SelectionListWindow          ( main = self, system_liststore =  self.system_liststore)
        self.window_list.append(self.selection_list_window)

        
        self.go_to_atom_window            = EasyHybridGoToAtomWindow     ( main = self, system_liststore = self.system_liststore)
        self.window_list.append(self.go_to_atom_window)

        
        #self.pDynamo_selection_window     = PDynamoSelectionWindow       ( main = self)
        self.pDynamo_selection_window     = EasyHybridSelectionWindow       ( main = self)
        
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
        
        
        self.normal_modes_window          =   NormalModesWindow (main = self)
        self.window_list.append(self.normal_modes_window)

        
        self.WHAM_window   =  WHAMWindow(main = self)
        self.window_list.append(self.WHAM_window)

        self.simple_dialog =  SimpleDialog(main = self)


        self.merge_system_window = MergeSystemWindow(main = self)
        self.solvate_system_window = SolvateSystemWindow(main = self)
        self.preferences_window = EasyHybridPreferencesWindow(main = self)

        self.make_solvent_box_window = MakeSolventBoxWindow(main = self)
        '''#- - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -#'''

        #self.window.connect("destroy", Gtk.main_quit)
        #self.window.connect("delete-event", Gtk.main_quit)
        #self.window.connect("destroy",      self.quit_easyhybrid)
        self.window.connect("delete-event", self.on_delete_event)
        self.window.connect("check-resize", self.window_resize)
        
        
        #self.builder.get_object('button_test').connect("clicked",    self.run_test)
        
        self.window.show_all()

    def restart (self):
        """ Function doc """
        self.main_treeview.treestore.clear()
        self.system_liststore.clear()      
        self.vobject_liststore_dict.clear()
        #self.bottom_notebook.treeview.clear()
        self.session_filename = None
        self.vm_session.restart()

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
        
        
        
        vismol_object.liststore_iter = self.vobject_liststore_dict[e_id].append([vismol_object.name, 
                                                                                 vismol_object.index, 
                                                                                 vismol_object.e_id,
                                                                                 sqr_color]
                                                                                 )

    def on_main_toolbar_clicked (self, button):
        """ Function doc """
        if button  == self.builder.get_object('toolbutton_new_system'):
            self.NewSystemWindow.OpenWindow()
        
        if button  == self.builder.get_object('toolbutton_save'):
            
            if self.session_filename == None:
                self.gtk_save_as_file(button)
                
            else:
                self.p_session.save_easyhybrid_session( filename = self.session_filename)
        
        if button  == self.builder.get_object('toolbutton_save_as'):
            self.gtk_save_as_file (button)
        
        if button  == self.builder.get_object('_show_cell'):
            print('aqui')
            self.run_test(None)
            
        if button == self.builder.get_object('toolbutton_terminal'):
            if button.get_active ():
                self.terminal_window.OpenWindow()
            else:
                self.terminal_window.CloseWindow(None, None)
        
        if button == self.builder.get_object('toolbutton_trajectory_tool1'):
            if button.get_active ():
                self.trajectory_player_window.OpenWindow()
                #self.traj_frame.hide()
            else:
                self.trajectory_player_window.CloseWindow(button = None)
        
        if button == self.builder.get_object('button_go_to_atom'):
            self.go_to_atom_window.OpenWindow()

        if button  == self.builder.get_object('selections'):
            self.selection_list_window.OpenWindow()

        if button  == self.builder.get_object('toolbutton_energy'):
            #self.p_session.run_simulation (parameters = {'simulation_type' : 'Energy', 'system': self.p_session.psystem[self.p_session.active_id]})
            #self.energy_refinement_window.OpenWindow()
            self.single_point_window.OpenWindow()
            
        if button  == self.builder.get_object('toolbutton_setup_QCModel'):
            self.setup_QCModel_window.OpenWindow()
        
        if button  == self.builder.get_object('toolbutton_system_check'): 
            self.p_session.systems[self.p_session.active_id]['vobject'].get_backbone_indexes ()

        if button  == self.builder.get_object('toolbutton_geometry_optimization'):
            self.geometry_optimization_window.OpenWindow()
        
        if button  == self.builder.get_object('toolbutton_pDynamo_selections'):
            self.pDynamo_selection_window.OpenWindow()
            
            '''
            atom1 = self.vm_session.picking_selections.picking_selections_list[0]
            ##print (atom1.chain, atom1.resn, atom1.resi, atom1.name)
            #print ("%s:%s.%s:%s" %(atom1.chain, atom1.resn, atom1.resi, atom1.name))
            
            _centerAtom ="%s:%s.%s:%s" %(atom1.chain, atom1.resn, atom1.resi, atom1.name)
            _radius =  10.0
            
            self.p_session.selections (_centerAtom, _radius)
            '''
       
        if button  == self.builder.get_object('toolbutton_normal_modes'):
            self.normal_modes_window.OpenWindow()
        
        if button  == self.builder.get_object('toolbutton_pes_scan'):
            self.PES_scan_window.OpenWindow()
        
        #if button  == self.builder.get_object('toolbutton_pes_scan'):
        #    self.PES_scan_window.OpenWindow()
        
        if button  == self.builder.get_object('toolbutton_molecular_dynamics'):
            self.molecular_dynamics_window.OpenWindow()
            
        if button  == self.builder.get_object('toolbutton_umbrella_sampling'):
            ##print('toolbutton_umbrella_sampling')
            self.umbrella_sampling_window.OpenWindow()
        
        if button  == self.builder.get_object('toolbutton_chain_of_states_opt'):
            ##print('toolbutton_umbrella_sampling')
            self.chain_of_states_opt_window.OpenWindow()
        
        if button  == self.builder.get_object('toolbutton_monte_carlo'):
            ##print('toolbutton_umbrella_sampling')
            self.normal_modes_analysis_window.OpenWindow()

    def on_main_menu_activate (self, menuitem):
        """ Function doc """
        ##print(menuitem)
        
        if menuitem == self.builder.get_object('menuitem_new'):
            self.NewSystemWindow.OpenWindow()
        
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
            self.export_data_window.OpenWindow()
        
        elif menuitem == self.builder.get_object('menuitem_import'):
            self.import_trajectory_window.OpenWindow()
            
        elif menuitem == self.builder.get_object('menuitem_quit'):
            self.on_delete_event(self, None)
            #Gtk.main_quit()
            pass
        
        
        #----------------------------------------------------------------------
        #                           E D I T 
        #----------------------------------------------------------------------
        elif menuitem == self.builder.get_object('menuitem_make_solvent_box'):
            self.make_solvent_box_window.OpenWindow()
        elif menuitem == self.builder.get_object('menuitem_go_to_atom'):
            self.go_to_atom_window.OpenWindow()
        
        elif menuitem == self.builder.get_object('menuitem_preferences'):
            self.preferences_window.OpenWindow()
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
            self.selection_list_window.OpenWindow()

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
            self.pDynamo_selection_window.OpenWindow()
        
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
                
                dialog =  EasyHybridDialogPrune(self.home ,num_of_atoms, name, tag)

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
            self.terminal_window.OpenWindow()
            if self.terminal_window.visible:
                button.set_active (True) 
            else:
                self.terminal_window.CloseWindow(button = None)



        elif menuitem == self.builder.get_object('menuitem_trajectory_tool'):
            button = self.builder.get_object('toolbutton_trajectory_tool1') 
            
            self.trajectory_player_window.OpenWindow()
            if self.trajectory_player_window.Visible:
                button.set_active(True)
            else:
                self.trajectory_player_window.CloseWindow(button = None)

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
            self.setup_QCModel_window.OpenWindow()
            
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
            self.solvate_system_window.OpenWindow()
        
        elif menuitem == self.builder.get_object('menuitem_merge'): 
            """ Function doc """
            system = self.p_session.psystem[self.p_session.active_id]
            e_id = system.e_id
            self.merge_system_window.selected_system_id = e_id
            self.merge_system_window.OpenWindow()
            
        
        
        
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
            
            new_system = self.p_session.append_system_to_pdynamo_session (system = new_system)
            self.main_treeview.add_new_system_to_treeview (new_system)
            ff  =  getattr(new_system.mmModel, 'forceField', "None")
        
            self.bottom_notebook.status_teeview_add_new_item(message = 'New System:  {} ({}) - Force Field:  {}'.format(new_system.label, new_system.e_tag, ff), system = new_system)
            self.p_session._add_vismol_object_to_easyhybrid_session (new_system, True) #, name = 'olha o  coco')
            #----------------------------------------------------------------------
        elif menuitem == self.builder.get_object('menuitem_prune_to_qc_region'): 
            system = self.p_session.psystem[self.p_session.active_id]
            
            new_system = system.PruneToQCRegion()
            new_system = self.p_session.append_system_to_pdynamo_session (system = new_system)            
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
            self.single_point_window.OpenWindow()
            
        elif menuitem == self.builder.get_object('menuitem_geometry_optimization'):
            self.geometry_optimization_window.OpenWindow()
            
        elif menuitem == self.builder.get_object('menuitem_molecular_dynamics'):
            self.molecular_dynamics_window.OpenWindow()
            
        elif menuitem == self.builder.get_object('menuitem_normal_modes'):
            self.normal_modes_window.OpenWindow()
            
        elif menuitem == self.builder.get_object('menuitem_rection_coordinate_scans'):
            self.PES_scan_window.OpenWindow()
            
        elif menuitem == self.builder.get_object('menuitem_nudged_elastic_band'):
            self.chain_of_states_opt_window.OpenWindow()
            
        elif menuitem == self.builder.get_object('menuitem_umbrella_sampling'):
            self.umbrella_sampling_window.OpenWindow()
        


        elif menuitem == self.builder.get_object('menuitem_energy_analysis'):
            self.PES_analysis_window.OpenWindow()

        elif menuitem == self.builder.get_object('menuitem_mormal_modes_analysis'):
            self.normal_modes_analysis_window.OpenWindow()
        
        elif menuitem == self.builder.get_object('menuitem_surface_analysis'):
            self.surface_analysis_window.OpenWindow()


        
        elif menuitem == self.builder.get_object('menuitem_energy_refinement'):
            self.energy_refinement_window.OpenWindow()
        
        
        elif menuitem == self.builder.get_object('menuitem_WHAM'):
            self.WHAM_window.OpenWindow()
        
        
        
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
        ''' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '''
        
        #'''
        filter = Gtk.FileFilter()  
        filter.set_name("EasyHybrid3 session files - *.easy")
        filter.add_mime_type("Easy files files")
        filter.add_pattern("*.easy")
        filters.append(filter)
        #'''
        

        filter = Gtk.FileFilter()  
        filter.set_name("pDynamo3 system files - *.pkl")
        filter.add_mime_type("PKL files")
        filter.add_pattern("*.pkl")
        filters.append(filter)
        
        filter = Gtk.FileFilter()  
        filter.set_name("pDynamo3 system files - *.yaml")
        filter.add_mime_type("YAML files")
        filter.add_pattern("*.yaml")
        filters.append(filter)
        
        filter = Gtk.FileFilter()
        filter.set_name("All files")
        filter.add_pattern("*")
        filters.append(filter)
        ''' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '''
        
        
        ''' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '''
        filename = self.filechooser.open(filters = filters)
        ''' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '''
    
        if filename:
            if filename[-4:] == 'easy':
                if os.path.exists(filename+'~'):
                    msg = 'There is a newer temporary file for the project you are loading.\n Would you like to load the most current file?'
                    dialog = SimpleDialog(self)
                    yes_or_no = dialog.question (msg)
                    if yes_or_no:
                        self.p_session.load_easyhybrid_serialization_file(filename+'~', tmp = True)
                    else:
                        self.p_session.load_easyhybrid_serialization_file(filename)
                else:
                    self.p_session.load_easyhybrid_serialization_file(filename)

            elif filename[-5:] == 'easy~':
                #print('ehf file')            
                #self.save_vismol_file = filename
                self.p_session.load_easyhybrid_serialization_file(filename)            
            else:
                files = {'coordinates': filename}
                systemtype = 3
                self.p_session.load_a_new_pDynamo_system_from_dict(files, systemtype)
        else:
            pass
        
    def run_dialog_set_QC_atoms (self, _type = None, vismol_object = None):
        """ Function doc """
        dialog = EasyHybridDialogSetQCAtoms(self.window)
        response = dialog.run()

        if response == Gtk.ResponseType.YES:
            ##print("The OK button was clicked")
            self.setup_QCModel_window.OpenWindow(vismol_object)
        elif response == Gtk.ResponseType.CANCEL:
            pass
            #print("The Cancel button was clicked")

        dialog.destroy()

    def uptade_interface_windows_and_dialogs (self, parameters = None):
        """ Function doc """
        for window in self.window_list:
            window.update()
        
    def change_reference_color (self, system, new_color):
        """ Function doc """
        system.e_color_palette['C'] = new_color
        sqr_color                   = get_colorful_square_pixel_buffer (system)
        
        
        self.system_liststore[system.e_liststore_iter][2] = sqr_color 

        self.main_treeview_model =  self.main_treeview.get_model() 
        self.main_treeview_model[system.e_treeview_iter][9] = sqr_color
        
        for index, vobject in self.vm_session.vm_objects_dic.items():
            treeview_iter = vobject.e_treeview_iter
            if self.main_treeview_model[treeview_iter][0] == system.e_id:
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
        if v_id == -1:
            _iter = self.p_session.psystem[e_id].e_treeview_iter
            self.main_treeview.treestore[_iter][2] = str(e_id)+' - '+ name
            self.p_session.psystem[e_id].label  = name
            liststore_iter = self.p_session.psystem[e_id].e_liststore_iter
            self.system_liststore[liststore_iter][0] = str(e_id)+' - '+ name
        
        else:
            _iter = self.vm_session.vm_objects_dic[v_id].e_treeview_iter
            self.main_treeview.treestore[_iter][2] = name
            self.vm_session.vm_objects_dic[v_id].name = name
            self.vobject_liststore_dict[e_id][self.vm_session.vm_objects_dic[v_id].liststore_iter][0] = name
            
    def delete_system (self, system_e_id = None ):
        """ 
        system_e_id = is the access key to the object. You can get it from vobject.e_id

        1) remove vobjects from vm_object_dic (self.vm_object_dic in vm_session object)
        2) remove vobjects from vobject_liststore_dict (self.vobject_liststore_dict in main object)
        
        3) remove system from system_liststore (self.system_liststore in main object)
        4) remove system from treestore (system.e_treeview_iter)
        
        5) remove system from p_session (p_session.psystem[sys_e_id] in p_session object)
        """
        #print(system_e_id)
        #parent_key = self.treeview.main.p_session.psystem[system_e_id].e_treeview_iter_parent_key       
        
        if system_e_id != None:
            
            '''organizing the list of vobjects that should be removed from vm_object_dic'''
            pop_list = []
            for index, vobject in self.vm_session.vm_objects_dic.items():
                if vobject.e_id == system_e_id:                   
                    self.main_treeview.treestore.remove(vobject.e_treeview_iter)
                    pop_list.append(index)
            '''removing from vm_object_dic'''
            for index in pop_list:
                self.vm_session.vm_objects_dic.pop(index)
            
            '''removing vobject from vobject_liststore_dict'''
            self.vobject_liststore_dict.pop(system_e_id)
            
            
            #  - - - - - - - - removing system treeview items - - - - - - - - - - -
            system = self.p_session.get_system(index = system_e_id)
            self.system_liststore.remove(system.e_liststore_iter)
            self.main_treeview.treestore.remove(system.e_treeview_iter)
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
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
            
            
            #  - - - - REMOVING plotting data FROM system.e_logfile_data - - - - - -
            system = self.p_session.psystem[vobject.e_id]
            #print(vobject.e_id, vm_object_index)
            if vm_object_index in system.e_logfile_data.keys():
                print ('deleting plotting data...')
                system.e_logfile_data[vm_object_index] = None
                system.e_logfile_data.pop(vm_object_index)
            #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
            #  - - - - - - REMOVING vobj FROM  vobject_liststore_dict - - - - - - -
            self.vobject_liststore_dict[vobject.e_id].remove(vobject.liststore_iter)
            #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
            #  - - - - - - REMOVING vobj FROM  treestore - - - - - - -
            self.main_treeview.treestore.remove(vobject.e_treeview_iter)
            #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
            #  - - - - - - - - REMOVING vobj FROM vm_object_dic - - - - - - - - - -
            self.vm_session.vm_objects_dic[vm_object_index] = None
            self.vm_session.vm_objects_dic.pop(vm_object_index)# = None
            #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            self.vm_session.vm_glcore.queue_draw()
        
    def refresh_main_statusbar(self, message = None, psystem = None):
        """ Function doc """
        psystem = self.p_session.psystem[self.p_session.active_id]
            
        if message:
            string = message
            print(string)
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
        if statusbar:
            self.refresh_main_statusbar()

    def print_btn (self, widget):
        """ Function doc """
        #print(self.box_reac)
        parm = self.box_reac.get_rc_data()

    def run_menu_item_test (self, widget):
        """ Function doc """
        #print(self.vm_session.picking_selections.picking_selections_list)
        pklist = self.vm_session.picking_selections.picking_selections_list
        if pklist[0]  and pklist[1]:
            print(pklist[0], pklist[1])
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

        #'''
        #print('aloowww')
        #print(self.vm_session.vm_glcore.glcamera)
        print(self.vm_session.vm_glcore.glcamera.z_near          ) #= self.z_far
        print(self.vm_session.vm_glcore.glcamera.z_far        ) #= self.fog_end - self.min_zfar
        print('\nview matrixes: \n', self.vm_session.vm_glcore.glcamera.view_matrix      ) #= self._get_view_matrix(pos)
        print('\nprojection_matrix: \n',self.vm_session.vm_glcore.glcamera.projection_matrix) #= self._get_projection_matrix()
        print('\ncamera position: \n', self.vm_session.vm_glcore.glcamera.get_position()) #= self._get_projection_matrix()
        
        camera = self.vm_session.vm_glcore.glcamera.get_position()
        for key, vobj in self.vm_session.vm_objects_dic.items():
            print(vobj.model_mat)
            print(vobj.trans_mat)
        
        
        #inv_mat = mop.get_inverse_matrix(self.vm_session.vm_glcore.glcamera.view_matrix )
        #print ('inverse of nview', inv_mat)
        #
        #print ('inverse of nview', inv_mat[3])
        
        
        view_matrix = self.vm_session.vm_glcore.glcamera.view_matrix
        forward = view_matrix[2]
        
        
        print(vobj.colors[1][0])
        pov_file  = open('temp.pov', 'w')

        text = '// Generated by EasyHybrid \n'

        text+= '#include "colors.inc"\n'
        text+= '#include "textures.inc"\n'
        text+= '#include "shapes.inc"\n'
        text+= '#include "stones1.inc"\n'

        text+= '// Rotation matrix\n'
        text+= '#declare myTransforms = transform {\n'
        text+= '    matrix <{}, {} ,{}, {}, {}, {}, {}, {}, {}, 0.000000, 0.000000 ,0.000000>\n'.format(vobj.model_mat[0][0],  vobj.model_mat[0][1], vobj.model_mat[0][2],
                                                                                                        vobj.model_mat[1][0],  vobj.model_mat[1][1], vobj.model_mat[1][2],
                                                                                                        vobj.model_mat[2][0],  vobj.model_mat[2][1], vobj.model_mat[2][2],
                                                                                                        )
        text+= '}\n'



        text+= '// CAMERA                                                          \n'
        text+= 'camera                                                             \n'
        text+= '{	//orthographic                                                 \n'
        #text+= '    right     1.26887871853547 *x                                  \n'
        #text+= '    up        y                                                    \n'
        text+= '    direction -z                                                   \n'
        text+= '    location  < {}, {},  {} >                                      \n'.format(camera[0], camera[1], camera[2])
        text+= '    angle   {}                                      \n'.format(self.vm_session.vm_glcore.glcamera.field_of_view*3)
        text+= '    look_at   < {}, {},  {} >                                      \n'.format(forward[0],forward[1],forward[2]) #(vobj.model_mat[3][0],  vobj.model_mat[3][1], vobj.model_mat[3][2])
        #text+= '     scale     20.91588996061168                                   \n'
        text+= '    translate < {} , {} , {} >                                     \n'.format(-vobj.model_mat[3][0], -vobj.model_mat[3][1], -vobj.model_mat[3][2] ) 
        text+= '}                                                                  \n'
        text+= '                                                                   \n'
        text+= '  // LIGHT 1                                                       \n'
        text+= 'light_source                                                       \n'
        text+= '{                                                                  \n'
        text+= '    <  0.000000,  0.000000,4000>                             \n'
        text+= '    color 0.8*White                                                \n'
        text+= '}                                                                  \n'
        text+= '// BACKGROUND                                                      \n'
        text+= 'background                                                         \n'
        text+= '{                                                                  \n'
        text+= '    color rgb < 1.000000, 1.000000, 1.000000 >                     \n'
        text+= '}                                                                  \n'
        text+= '                                                                   \n'
        text+= '# declare molecule = union {                                       \n'
        text+= '// ATOMS                                                           \n'
        
                                                                
        
        for key, vobj in self.vm_session.vm_objects_dic.items():
            for i, atom in vobj.atoms.items():
                #print(i, atom)
                xyz = atom.coords(frame = 0)

                text+= 'sphere                                                             \n'
                text+= '{                                                                  \n'
                text+= '    <      {},    {},     {}>       {}  \n'.format(xyz[0], xyz[1], xyz[2], atom.vdw_rad*0.2)
                text+= '    texture { finish { Dull } }                                    \n'
                text+= '    pigment { '+'rgb<     {},      {},      {}'.format(vobj.colors[i][0], vobj.colors[i][1],vobj.colors[i][2])+'> }  \n'
                text+= '}                                                                  \n'

            for bond in vobj.bonds:
                xyz_i = bond.atom_i.coords(frame = 0)
                xyz_j = bond.atom_j.coords(frame = 0)
                
                text+= 'cylinder                                                          \n'
                text+= '{                                                                 \n'
                text+= '    <      {},      {},     {}>,               \n'.format(xyz_i[0], xyz_i[1], xyz_i[2])
                text+= '    <      {},      {},     {}>                \n'.format(xyz_j[0], xyz_j[1], xyz_j[2])
                text+= '          0.15                                                \n'
                text+= '    texture { finish { Dull } }                                   \n'
                text+= '    pigment { rgb<      1,      1,      1> } \n'
                text+= '}                                                                 \n'


        text+='transform { myTransforms }                                              \n'
        text+='}                                                                       \n'
        text+='                                                                        \n'
        #text+='object {molecule}                                                       \n'
 
        text+='intersection{                                              \n'                          
        text+='    object {molecule}                                      \n'                    
        text+='    plane{z, %s pigment{rgbf<1,1,1,0>}}                    \n'% str(-self.vm_session.vm_glcore.glcamera.z_near)
        text+='}                                                          \n'
        #text+= 'plane                                                              \n'
        #text+= '{                                                                  \n'
        #text+= '	z, {}                                                          \n'.format( -self.vm_session.vm_glcore.glcamera.z_far )
        #text+= '    pigment { rgb <   1.000000, 1.000000, 1.000000 > }                  \n'#.format(0.0, 0.0 , self.vm_session.vm_glcore.glcamera.fog_end )
        #text+= '} \n' 

        pov_file.write(text)
        #print(text)
        #os.system('povray +A0.3 -UV +W1000 +H1000 +Itemp.pov +Otemp2.png')
        print('povray +A0.3 -UV +W1000 +H1000 +Itemp.pov +Otemp2.png')
        
        
        print(self.vm_session.vm_glcore.glcamera.z_near, self.vm_session.vm_glcore.glcamera.z_far)
        #print(self.vm_session.vm_glcore.glcamera._get_view_matrix([0,0,10]))
        print(self.vm_session.vm_glcore.glcamera.fog_end  )
        print(self.vm_session.vm_glcore.glcamera.fog_start)
        print(self.vm_session.vm_glcore.zero_reference_point)
        #'''
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        '''
        print('model matrix:',      self.vm_session.vm_glcore.model_mat                 , type(self.vm_session.vm_glcore.model_mat                 )  )
        print('projection matrix:', self.vm_session.vm_glcore.glcamera.projection_matrix, type(self.vm_session.vm_glcore.glcamera.projection_matrix)  )
        print('view matrix:',       self.vm_session.vm_glcore.glcamera.view_matrix      , type(self.vm_session.vm_glcore.glcamera.view_matrix      )  )
        
        proj_mat = self.vm_session.vm_glcore.glcamera.projection_matrix
        view_mat = self.vm_session.vm_glcore.glcamera.view_matrix
        #model_mat= self.vm_session.vm_glcore.model_mat
        # Multiplica as matrizes na ordem correta

        

        pov_file  = open('temp.pov', 'w')
        
        text = '// Generated by EasyHybrid \n'
        
        text+= '#include "colors.inc"\n'
        text+= '#include "textures.inc"\n'
        text+= '#include "shapes.inc"\n'
        text+= '#include "stones1.inc"\n'
        
        text+= '\n'
        text+= 'camera {\n'
        text+= '    location <0,0,-100>\n'
        #text+= '    location <{},{},{}>\n'.format(camera_position_povray[0],camera_position_povray[1],camera_position_povray[2])
        text+= '    look_at  <0,0,0>\n'
        #text+= '    look_at  <{},{},{}>\n'.format(look_at_point_povray[0],look_at_point_povray[1],look_at_point_povray[2])
        text+= '    angle 45\n'
        text+= '}\n'

        text+= '                                                                   \n'
        text+= '  // LIGHT 1                                                       \n'
        text+= 'light_source                                                       \n'
        text+= '{                                                                  \n'
        text+= '    <  0.000000,  0.000000,-109.953541>                             \n'
        text+= '    color 0.8*White                                                \n'
        text+= '}                                                                  \n'
        text+= '// BACKGROUND                                                      \n'
        text+= 'background                                                         \n'
        text+= '{                                                                  \n'
        text+= '    color rgb < 1.000000, 1.000000, 1.000000 >                     \n'
        text+= '}                                                                  \n'
        text+= '                                                                   \n'
        text+= '# declare molecule = union {                                       \n'
        text+= '// ATOMS                                                           \n'
        

        for key, vobj in self.vm_session.vm_objects_dic.items():
            if vobj.active:
                model_mat = vobj.model_mat
                result = np.matmul(proj_mat, view_mat)
                result = np.matmul(result, model_mat)
                
                for i, atom in vobj.atoms.items():
                    #print(i, atom)
                    xyz = atom.coords(frame = 0)
                    #result = np.matmul(result, vert_coord)
                    xyz = np.append(xyz, 1.0)
                    xyz = np.matmul(result, xyz)
                    text+= 'sphere                                                             \n'
                    text+= '{                                                                  \n'
                    text+= '    <      {},    {},     {}>       {}  \n'.format(xyz[0], xyz[1], xyz[2], atom.vdw_rad*0.3)
                    text+= '    texture { finish { Dull } }                                    \n'
                    text+= '    pigment { '+'rgb<     {},      {},      {}'.format(vobj.colors[i][0], vobj.colors[i][1],vobj.colors[i][2])+'> }  \n'
                    text+= '}                                                                  \n'
        
                for bond in vobj.bonds:
                    xyz_i = bond.atom_i.coords(frame = 0)
                    xyz_j = bond.atom_j.coords(frame = 0)
                    
                    text+= 'cylinder                                                          \n'
                    text+= '{                                                                 \n'
                    text+= '    <      {},      {},     {}>,               \n'.format(xyz_i[0], xyz_i[1], xyz_i[2])
                    text+= '    <      {},      {},     {}>                \n'.format(xyz_j[0], xyz_j[1], xyz_j[2])
                    text+= '          0.20                                                \n'
                    text+= '    texture { finish { Dull } }                                   \n'
                    text+= '    pigment { rgb<      0.193545,      0.193545,      0.193545> } \n'
                    text+= '}                                                                 \n'
        text+= '}                                                                 \n'
        text+='object {molecule}'
        pov_file.write(text)
        
        width  = self.vm_session.vm_glcore.width 
        height = self.vm_session.vm_glcore.height
        
        os.system('povray +A0.1 -UV +W{} +H{} +Itemp.pov +Otemp2.png'.format(width, height))
        print('povray +A0.1 -UV +W{} +H{} +Itemp.pov +Otemp2.png'.format(width, height))
        #'''
 

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
    


class EasyHybridMainTreeView(Gtk.TreeView):
    
    def __init__(self):
        super().__init__( )

        self.treestore = Gtk.TreeStore(int , #0 system_e_id           
                                       int , #1 vobject 
                                       str , #2 name 
                                       
                                       bool, #3 is the radiobutton visible?
                                       bool, #4 is the radiobutton active?
                                       
                                       bool, #5 is the Toggle box visible?
                                       bool, #6 is the toggle boz active?
                                       
                                       bool, #7 is the Frames visible?
                                       int , #8 Number of frames
                                       
                                       GdkPixbuf.Pixbuf #9 Color
                                       )
        
        
        self.set_model(self.treestore)
        self.tree_iters         = []                               
        
        '''self.tree_iters_dict stores the path to 
        access the vobject in the treeviewm. 
        Every vobject has a key to access the "iter"'''
        self.tree_iters_dict    = {}
        self.tree_iters_counter = 0
        self._create_treeview ()
        


    
    def add_new_system_to_treeview (self, system):
        """ Function doc """
        sqr_color   = get_colorful_square_pixel_buffer(system)
        #color      =  system.e_color_palette['C']
        #res_color  = [int(color[0]*255),int(color[1]*255),int(color[2]*255)] 
        #sqr_color  =  getColouredPixmap( res_color[0], res_color[1], res_color[2] )
        
        for row in self.treestore:
            row[4] = False
        parent = self.treestore.append(None, [system.e_id, -1,str(system.e_id)+' - '+ system.label, True, True, False, False, False, 0, sqr_color])
        system.e_treeview_iter = parent
        
        
        '''To improve organization and accessibility, we will add a GtkListStore 
        to a dictionary that will be accessed by all windows. Each GtkListStore 
        in the dictionary will contain the vobjects for a particular system or project'''


       
        system.e_liststore_iter = self.main.system_liststore.append([str(system.e_id)+' - '+ system.label, system.e_id, sqr_color])
        self.main.vobject_liststore_dict[system.e_id] = Gtk.ListStore(str, int, int,sqr_color)
        self.main.bottom_notebook.set_active_system_text_to_textbuffer()
    
    def add_vismol_object_to_treeview(self, vismol_object):
        """ Function doc """
        e_id   = vismol_object.e_id
        system = self.main.p_session.get_system(e_id)
        sqr_color   = get_colorful_square_pixel_buffer (system)
        #color      =  system.e_color_palette['C']
        #res_color  = [int(color[0]*255),int(color[1]*255),int(color[2]*255)] 
        #sqr_color  =  getColouredPixmap( res_color[0], res_color[1], res_color[2] )
        
        parent = system.e_treeview_iter #self.tree_iters_dict[system.e_treeview_iter_parent_key]
        
        size   = len(vismol_object.frames)
        _iter  = self.treestore.append(parent, [e_id,  vismol_object.index , vismol_object.name, False, False , True, vismol_object.active, True, size, sqr_color])
        
        self.tree_iters.append(_iter)
        self.tree_iters_dict[self.tree_iters_counter] = parent
        vismol_object.e_treeview_iter_parent_key = self.tree_iters_counter
        
        vismol_object.e_treeview_iter = _iter
        
        self.tree_iters_counter += 1
        
        for row in self.treestore:
            pass

        self.expand_row(row.path, True)
        self.refresh_trajectory_scalebar()
    
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
        
        self.refresh_trajectory_scalebar()
        '''
        higher = 1
        for index, vobject in self.main.vm_session.vm_objects_dic.items():
            treeview_iter = vobject.e_treeview_iter
            if vobject.active:
                size = len(vobject.frames)
                if size > higher:
                    higher = size
                else:
                    pass
            else:
                pass
            
            #self.treestore[treeview_iter][8] = size
            #print(index, self.treestore[treeview_iter][2], 'frames', len(vobject.frames))
        self.main.trajectory_player_window.change_range(upper = higher)
        self.main.trajectory_player_window.upper = higher
        #'''
        self.main.vm_session.vm_glcore.queue_draw()

    
    def on_cell_active_radio_toggled(self, widget, path):
        ##print (path)
        selected_path = Gtk.TreePath(path)
        
        ##print(self.treestore[path][1], path, self.treestore[path][0])
        #
        for row in self.treestore:
            ##print(row.path, selected_path)
            
            active_id_before = self.main.p_session.active_id
            
            row[4] = row.path == selected_path
            if row.path == selected_path:
                system_e_id = row[0]
                self.main.p_session.active_id = system_e_id
                active_id_after = system_e_id
                
                self.main.bottom_notebook.change_annotations_textbuffer(
                                                                  active_id_before,
                                                                  active_id_after
                                                                  )
                #print (system_e_id, self.main.p_session.psystem[system_e_id].label)
            else:
                pass
        #self.main.refresh_active_system_liststore ()
        self.main.refresh_main_statusbar ()
        self.main.uptade_interface_windows_and_dialogs()
    
    def on_select(self, tree, path, selection):
        '''---------------------- Row information ---------------------'''
        # Get the current selected row and the model.
        model, iter = tree.get_selection().get_selected()        
        
        # Look up the current value on the selected row and get
        # a new value to change it to.
        data2 = model.get_value(iter, 2)
        data1 = model.get_value(iter, 1)
        data0 = model.get_value(iter, 0)
        #print(data0, data1, data2, self.treestore[path][8])
        self.main.trajectory_player_window.change_range(upper = self.treestore[path][8])

        selection             = tree.get_selection()
        model                 = tree.get_model()
        self.vm_object_index  = int(model.get_value(iter, 1))
        self.system_e_id      = int(model.get_value(iter, 0))
        '''------------------------------------------------------------'''
        
        
        # - - - - - - printing the information form all rows - - - - - - - 
        #'''
        for index , vobject  in self.main.vm_session.vm_objects_dic.items():
            treeview_iter = vobject.e_treeview_iter
            
            vobj_id = self.treestore[treeview_iter][1]
            if vobj_id == -1:
                pass
                #print ('vobj_id', vobj_id)
            else:
                vobj = self.main.vm_session.vm_objects_dic[vobj_id]
                size = len(vobj.frames)
                self.treestore[treeview_iter][8] = size
                #print(vobj_id, self.treestore[treeview_iter][2], 'frames', len(vobj.frames))
        #'''

    def on_treeview_mouse_button_release_event (self, tree, event):
        """ Function doc """
        
        model, iter = tree.get_selection().get_selected()        
        selection   = tree.get_selection()
        model       = tree.get_model()
        try:
            self.vm_object_index  = int(model.get_value(iter, 1))
            self.system_e_id      = int(model.get_value(iter, 0))
        except:
            return False
        
        
        if event.button == 3:
            self.treeview_menu.open_menu(self.system_e_id, self.vm_object_index)

        if event.button == 2:
            #print('event.button 2',self.vm_object_index)
            if self.vm_object_index == -1:
                #means that is not a vismol object
                pass
            else:
                vismol_object = self.main.vm_session.vm_objects_dic[self.vm_object_index]
                self.main.vm_session.vm_glcore.center_on_coordinates(vismol_object, vismol_object.mass_center)

        if event.button == 1:
            #print('event.button 1')
            pass
    
    
    def refresh_number_of_frames (self):
        """ 
        This function refreshes the number of frames  on the main treeview.  
        The self.tree_iters list contains all the "parents", or the treeview lines, in the TreeView
        vismol_object.e_treeview_iter_parent_key
        """
        for index, vobject in self.main.vm_session.vm_objects_dic.items():
            
            treeview_iter = vobject.e_treeview_iter
            size = len(vobject.frames)
            self.treestore[treeview_iter][8] = size
            #print(index, self.treestore[treeview_iter][2], 'frames', len(vobject.frames))


    def refresh_trajectory_scalebar (self):
        """ Function doc """
        higher = 1
        for index, vobject in self.main.vm_session.vm_objects_dic.items():
            treeview_iter = vobject.e_treeview_iter
            if vobject.active:
                size = len(vobject.frames)
                if size > higher:
                    higher = size
                else:
                    pass
            else:
                pass
        self.main.trajectory_player_window.change_range(upper = higher)
        self.main.trajectory_player_window.upper = higher


    def refresh (self):
        """ Function doc """
        self.treestore.clear()
        
        for e_id in self.main.p_session.psystem.keys():
            system = self.main.p_session.psystem[e_id]
            self.add_new_system_to_treeview (system)
        
        
        
        for v_obj_index in self.main.vm_session.vm_objects_dic.keys():
            vismol_object = self.main.vm_session.vm_objects_dic[v_obj_index]
            self.add_vismol_object_to_treeview(vismol_object)


class TreeViewMenu:
    """ Class doc """
    
    def __init__ (self, treeview):
        """ Class initialiser """
        self.treeview = treeview
        self.main     = treeview.main 
        self.filechooser   = FileChooser()
        self.rename_window_visible =  False
        
        vobject_menu_items = {
                                'header'                : None    ,

                                '_separator'            : ''      ,
                                #'Info'                  : self.f2 ,
                                'Rename'                : self._menu_rename ,
                                '_separator'            : ''      ,
                                'Show / Hide Cell'      : self.show_or_hide_cell,
                                #'Load Data Into System' : self._menu_load_data_to_system,
                                #'Change Color'  : self.change_system_color_palette ,
                                'Go To Atom'            : self._menu_go_to_atom ,
                                '_separator'            : ''      ,

                                'Export As...'          : self._menu_export_data_window ,
                                #'Merge System With...'  : self.f3 ,
                                '_separator'            : ''      ,

                                'Delete'                : self._menu_delete_vm_object ,
                                #'test'  : self.f1 ,
                                #'f1'    : self.f1 ,
                                #'f2'    : self.f2 ,
                                #'gordão': self.f3 ,
                                #'delete': self.f3 ,
                                }

        system_menu_items = {
                                
                                'header'                  : None                            ,
                                                          
                                '_separator'              : ''                              ,
                                                          
                                'Info'                    : self._show_info                  ,
                                                          
                                '_separator'              : ''                              ,
                                                          
                                'Rename'                  : self._menu_rename               ,
                                'Import Data...'          : self._menu_load_data_to_system  ,
                                'Reference Color'           : self._menu_change_color_palette ,
                                #'Edit Parameters'         : self.f2                         ,
                                'Export As...'            : self._menu_export_data_window    ,
                                
                                '_separator'              : ''                              ,
                                'Merge With...'           : self._menu_merge_system         ,
                                'Clone'                   : self._menu_clone_system         ,
                                
                                '_separator'              : ''                              ,
                                
                                'Delete'                  : self._menu_delete_system        ,
                                #'test'  : self.f1 ,
                                #'f1'    : self.f1 ,
                                #'f2'    : self.f2 ,
                                #'gordão': self.f3 ,
                                #'delete': self.f3 ,
                                }
                    
                    
                    
        self.tree_view_vobj_menu  , self.tree_header_vobj_menu    = self.build_tree_view_menu(vobject_menu_items)
        self.tree_view_sys_menu   , self.tree_header_sys_menu     = self.build_tree_view_menu(system_menu_items)

    def show_or_hide_cell (self, widget):
        """ Function doc """
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()

        old_name = model.get_value(iter, 2)
        v_id     = model.get_value(iter, 1)
        e_id     = model.get_value(iter, 0)
        #----------------------------------------------------------------------        
        #system = self.main.p_session.psystem[e_id]
        
        vobject = self.main.vm_session.vm_objects_dic[v_id]
        if "cell_lines" in vobject.representations.keys():
            if vobject.representations["cell_lines"].active:
                self.main.vm_session.hide_cell (vobject)
            
            else:
                self.main.vm_session.show_cell (vobject)
        else:
            self.main.vm_session.show_cell (vobject)
            



    def _show_info (self, widget):
        """ Function doc """
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        e_id          = int(model.get_value(iter, 0))  # @+
        #----------------------------------------------------------------------        
        system = self.main.p_session.psystem[e_id]
        window = InfoWindow(system)
    

    def _menu_export_data_window (self,vobject = None ):
        """ Function doc """
        self.treeview.main.export_data_window.OpenWindow()
    
    
    def _menu_load_data_to_system (self, vobject = None ):
        """ Function doc """
        selection        = self.treeview.get_selection()
        model, iter      = selection.get_selected()
        #print (list(model))
        self.main.import_trajectory_window.OpenWindow(sys_selected = model.get_value(iter, 0))


    def _menu_change_color_palette (self, widget):
        """ Function doc """
        #selection               = self.selections[self.current_selection]
        self.colorchooserdialog = Gtk.ColorChooserDialog()
        
        if self.colorchooserdialog.run() == Gtk.ResponseType.OK:
            color = self.colorchooserdialog.get_rgba()
            #print(color.red,color.green, color.blue )
            new_color = [color.red, color.green, color.blue]

        else:
            new_color = False
        
        
    
        if new_color:
            self.colorchooserdialog.destroy()

            #----------------------------------------------------------------------
            selection     = self.treeview.get_selection()
            (model, iter) = selection.get_selected()
            self.selectedID  = int(model.get_value(iter, 0))  # @+
            #----------------------------------------------------------------------
            
            system = self.main.p_session.psystem[self.selectedID]
            
            self.main.change_reference_color(system, new_color)

    
    def _menu_merge_system (self, widget):
        """ Function doc """
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        e_id          = int(model.get_value(iter, 0)) 
        self.main.merge_system_window.selected_system_id = e_id
        self.main.merge_system_window.OpenWindow()
    
    
    def _menu_clone_system (self, widget):
        """ Function doc """
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        e_id          = int(model.get_value(iter, 0))
        self.main.p_session.clone_system(e_id)
        self._save_backup_file()
    
    def _menu_go_to_atom (self, vobject = None):
        """ Function doc """
        ##print('f2')
        #self._show_lines(vobject = self.vobjects[0], indices = [0,1,2,3,4] )
        self.treeview.main.go_to_atom_window.OpenWindow()
        #self.treeview.vm_session.go_to_atom_window.OpenWindow()

    def f3 (self, vobject = None):
        """ Function doc """
        
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()


        self.selectedID  = int(model.get_value(iter, 1))  # @+
        
        
        
        del self.treeview.vm_session.vobjects_dic[self.selectedID]
        '''
        vobject = self.treeview.vm_session.vobjects_dic.pop(self.selectedID)
        del vobject
        '''
        self.treeview.store.clear()
        for vobj_index ,vis_object in self.treeview.vm_session.vobjects_dic.items():
            data = [vis_object.active          , 
                    str(vobj_index),
                    vis_object.name            , 
                    str(len(vis_object.atoms)) , 
                    str(len(vis_object.frames)),
                   ]
            model.append(data)
        self.treeview.vm_session.glwidget.queue_draw()

    
    def _menu_rename (self, menu_item = None ):
        """  
        menu_item = Gtk.MenuItem object at 0x7fbdcc035700 (GtkMenuItem at 0x37cf6c0)
        
        """
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()

        old_name = model.get_value(iter, 2)
        v_id     = model.get_value(iter, 1)
        e_id     = model.get_value(iter, 0)
        tag      = self.main.p_session.psystem[e_id].e_tag 
        
        old_name = old_name.split("- ")
        old_name = old_name[-1]
        
        if self.rename_window_visible:
            self.preferences.set_names (old_name, tag)
            pass
        else:
            
            self.preferences = PreferencesWindow(main = self.main, 
                                                 e_id = e_id     ,
                                                 v_id = v_id     )
            self.preferences.set_names (old_name, tag)
        self._save_backup_file()
        
    def destroy (self, widget):
        """ Function doc """
        self.rename_window_visible = False
    
    def _menu_delete_vm_object (self, widget):
        """ Function doc """
        self.main.delete_vm_object ( vm_object_index = self.vobject_index)
        self._save_backup_file()
    
    def _menu_delete_system (self, widget):
        """ Function doc """
        self.main.delete_system (system_e_id = self.system_e_id )
        self._save_backup_file()
        #self.save_easyhybrid_session( filename = self.main.session_filename, tmp = True)
        
    def build_tree_view_menu (self, menu_items = None):
        """ Function doc """
        tree_view_menu = Gtk.Menu()
        menu_header    = None
        
        for label in menu_items:
            
            
            if menu_items[label] == None:
                # just a label
                
                mitem = Gtk.MenuItem(label = label)
                if label == 'header':
                    menu_header    = mitem
                
                
            elif  label == '_separator':
                mitem = Gtk.SeparatorMenuItem()
            
            else:
                mitem = Gtk.MenuItem(label = label)
                mitem.connect('activate', menu_items[label])
            
            tree_view_menu.append(mitem)
            #mitem = Gtk.SeparatorMenuItem()
            #self.tree_view_menu.append(mitem)

        tree_view_menu.show_all()
        return tree_view_menu, menu_header


    def open_menu (self, system_e_id = None , vobject_index = None):
        """ Function doc """
        self.system_e_id     = system_e_id    
        self.vobject_index = vobject_index
        ##print(system_e_id, vobject_index)
        
        system = self.treeview.main.p_session.psystem[self.system_e_id] 
        self.tree_header_sys_menu.set_label(system.label)
        
        if vobject_index == -1:
            
            self.tree_view_sys_menu.popup(None, None, None, None, 0, 0)

        if vobject_index != None and vobject_index != -1:
            
            vismol_object = self.treeview.main.vm_session.vm_objects_dic[vobject_index]
            self.tree_header_vobj_menu.set_label(vismol_object.name)
            
            self.tree_view_vobj_menu.popup(None, None, None, None, 0, 0)
                

    def _save_backup_file (self):
        """ Function doc """
        self.main.p_session.save_easyhybrid_session( filename = self.main.session_filename, tmp = True)
        
        
class BottomNoteBook:
    
    def __init__ (self, main):
        """ Function doc """
        self.home    = main.home
        self.main    = main
        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(self.home,'src/gui/MainWindow_text_and_logs.glade'))
        
        self.widget = self.builder.get_object('notebook_text_and_logs')

        self._define_status_treeview()
    

        self.annotations_textviewer = self.builder.get_object('annontations_textviewer')
        self.annotations_textbuffer = self.annotations_textviewer.get_buffer()

        #-----------------------------------------------------------------------
        seqview =  SequenceViewerBox(main = self.main)
        notebook = self.widget
        tab_label = Gtk.Label("sequence")
        notebook.append_page(seqview, tab_label)
        #-----------------------------------------------------------------------
        
        
        
        #seqview.add_events(Gdk.EventMask.POINTER_MOTION_MASK)
        #seqview.add_events(Gdk.EventMask.BUTTON_PRESS_MASK)
        #seqview.add_events(Gdk.EventMask.BUTTON_RELEASE_MASK)
        #
        #seqview.connect("motion-notify-event" , seqview.on_motion)
        #seqview.connect("button_press_event"  , seqview.button_press)
        #seqview.connect("button-release-event", seqview.button_release)
        ##-----------------------------------------------------------------------




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
        #label.set_text("Extracted Text: " + extracted_text)


class PreferencesWindow:
    """ Class doc """
    
    def __init__ (self, main = None, e_id = None , v_id = None):
        """ Class initialiser """
        
        self.xml='''
<?xml version="1.0" encoding="UTF-8"?>
<!-- Generated with glade 3.22.2 -->
<interface>
  <requires lib="gtk+" version="3.20"/>
  <object class="GtkWindow" id="window">
    <property name="can_focus">False</property>
    <property name="default_width">300</property>
    <child type="titlebar">
      <placeholder/>
    </child>
    <child>
      <object class="GtkAlignment">
        <property name="visible">True</property>
        <property name="can_focus">False</property>
        <property name="top_padding">5</property>
        <property name="bottom_padding">5</property>
        <property name="left_padding">5</property>
        <property name="right_padding">5</property>
        <child>
          <object class="GtkBox">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="orientation">vertical</property>
            <property name="spacing">5</property>
            <child>
              <object class="GtkGrid">
                <property name="visible">True</property>
                <property name="can_focus">False</property>
                <property name="row_spacing">5</property>
                <property name="column_spacing">5</property>
                <child>
                  <object class="GtkEntry" id="entry_name">
                    <property name="visible">True</property>
                    <property name="can_focus">True</property>
                    <property name="hexpand">True</property>
                    <property name="text" translatable="yes">new pDynamo system</property>
                  </object>
                  <packing>
                    <property name="left_attach">1</property>
                    <property name="top_attach">0</property>
                  </packing>
                </child>
                <child>
                  <object class="GtkLabel" id="label_name">
                    <property name="visible">True</property>
                    <property name="can_focus">False</property>
                    <property name="halign">end</property>
                    <property name="label" translatable="yes">Name:</property>
                  </object>
                  <packing>
                    <property name="left_attach">0</property>
                    <property name="top_attach">0</property>
                  </packing>
                </child>
                <child>
                  <object class="GtkLabel" id="label_tag">
                    <property name="visible">True</property>
                    <property name="can_focus">False</property>
                    <property name="halign">end</property>
                    <property name="label" translatable="yes">Tag: </property>
                  </object>
                  <packing>
                    <property name="left_attach">0</property>
                    <property name="top_attach">1</property>
                  </packing>
                </child>
                <child>
                  <object class="GtkEntry" id="entry_tag">
                    <property name="visible">True</property>
                    <property name="can_focus">True</property>
                    <property name="max_length">15</property>
                    <property name="width_chars">10</property>
                    <property name="text" translatable="yes">molsys</property>
                  </object>
                  <packing>
                    <property name="left_attach">1</property>
                    <property name="top_attach">1</property>
                  </packing>
                </child>
                <child>
                  <object class="GtkLabel" id="label_color">
                    <property name="visible">True</property>
                    <property name="can_focus">False</property>
                    <property name="label" translatable="yes">Color:</property>
                  </object>
                  <packing>
                    <property name="left_attach">0</property>
                    <property name="top_attach">2</property>
                  </packing>
                </child>
                <child>
                  <object class="GtkColorButton" id="button_color">
                    <property name="visible">True</property>
                    <property name="can_focus">True</property>
                    <property name="receives_default">True</property>
                    <property name="rgba">rgb(138,226,52)</property>
                  </object>
                  <packing>
                    <property name="left_attach">1</property>
                    <property name="top_attach">2</property>
                  </packing>
                </child>
              </object>
              <packing>
                <property name="expand">False</property>
                <property name="fill">True</property>
                <property name="position">0</property>
              </packing>
            </child>
            <child>
              <object class="GtkButtonBox" id="dialog-action_area2">
                <property name="visible">True</property>
                <property name="can_focus">False</property>
                <property name="layout_style">end</property>
                <child>
                  <object class="GtkButton" id="button_cancel">
                    <property name="label" translatable="yes">Cancel</property>
                    <property name="visible">True</property>
                    <property name="can_focus">True</property>
                    <property name="receives_default">True</property>
                  </object>
                  <packing>
                    <property name="expand">False</property>
                    <property name="fill">False</property>
                    <property name="position">0</property>
                  </packing>
                </child>
                <child>
                  <object class="GtkButton" id="button_apply">
                    <property name="label" translatable="yes">Apply</property>
                    <property name="visible">True</property>
                    <property name="can_focus">True</property>
                    <property name="receives_default">True</property>
                  </object>
                  <packing>
                    <property name="expand">False</property>
                    <property name="fill">False</property>
                    <property name="position">1</property>
                  </packing>
                </child>
              </object>
              <packing>
                <property name="expand">False</property>
                <property name="fill">False</property>
                <property name="position">2</property>
              </packing>
            </child>
          </object>
        </child>
      </object>
    </child>
  </object>
</interface>
'''        
        
        self.main       = main
        self.vm_session = main.vm_session
        self.p_session  = main.p_session
        
        self.e_id  = e_id 
        self.v_id  = v_id 
         
        
        self.builder = Gtk.Builder()
        self.builder.add_from_string(self.xml)

        self.window = self.builder.get_object('window')

        self.entry_name   = self.builder.get_object('entry_name')
        self.entry_tag    = self.builder.get_object('entry_tag')
        self.button_color = self.builder.get_object('button_color')
        
        self.button_cancel = self.builder.get_object('button_cancel')
        self.button_apply  = self.builder.get_object('button_apply')
        
        self.button_apply.connect('clicked', self.on_button_apply )
        self.button_cancel.connect('clicked', self.on_button_cancel )
        self.main.main_treeview.treeview_menu.rename_window_visible = True
        
        #self.window.connect("destroy", self.on_button_cancel)
        
        self.window.set_resizable(False)
        self.window.show_all()
        self.button_color.hide()
        self.builder.get_object('label_color').hide()
    
    def set_names (self, name, tag):
        """ Function doc """
        self.entry_name.set_text(name)
        self.entry_tag .set_text(tag)
        
    
    def on_button_apply (self, button):
        """ Function doc """
        name  = self.entry_name.get_text()
        tag   = self.entry_tag .get_text()
    
        self.rename_window_visible = False
        self.main.rename_tag(e_id = self.e_id, tag = tag)
        self.main.rename(e_id  = self.e_id, 
                         v_id  = self.v_id, 
                         name  = name)
        self.window.destroy()
        self.main.main_treeview.treeview_menu.rename_window_visible = False
    
    def on_button_cancel (self, button):
        """ Function doc """
        self.window.destroy()
        self.main.main_treeview.treeview_menu.rename_window_visible = False
