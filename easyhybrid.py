#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  simple.py
#
#  Copyright 2022 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
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

import logging
import gi 
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
from gi.repository import GdkPixbuf

import os, sys, time


from vismol.core.vismol_session import VismolSession
import vismol.utils.matrix_operations  as mop

from gEngine.eSession import EasyHybridSession
from gui.gtk_widgets import VismolSelectionTypeBox
from gui.gtk_widgets import FileChooser

from gui.windows.windows_and_dialogs import ImportANewSystemWindow
from gui.windows.windows_and_dialogs import EasyHybridDialogSetQCAtoms
from gui.windows.windows_and_dialogs import EasyHybridSetupQCModelWindow
from gui.windows.windows_and_dialogs import EasyHybridGoToAtomWindow
from gui.windows.windows_and_dialogs import PDynamoSelectionWindow
from gui.windows.windows_and_dialogs import ExportDataWindow
#from gui.windows.windows_and_dialogs import EnergyRefinementWindow
from gui.windows.windows_and_dialogs import ImportTrajectoryWindow
from gui.windows.windows_and_dialogs import TrajectoryPlayerWindow
from gui.windows.windows_and_dialogs import PotentialEnergyAnalysisWindow

from gui.windows.easyhybrid_terminal import TerminalWindow
from gui.windows.geometry_optimization_window import *
from gui.windows.selection_list_window      import *
from gui.windows.PES_scan_window            import PotentialEnergyScanWindow 

from gui.windows.molecular_dynamics          import MolecularDynamicsWindow 
from gui.windows.umbrella_sampling_window    import UmbrellaSamplingWindow 
from gui.windows.chain_of_states_opt_window  import ChainOfStatesOptWindow 
from gui.windows.normal_modes_analysis       import NormalModesAnalysisWindow 
from gui.windows.normal_modes_window         import NormalModesWindow 


from pdynamo.pDynamo2EasyHybrid import pDynamoSession

logger = logging.getLogger(__name__)
EASYHYBRID_VERSION = '3.0'
EASYHYBRID_HOME    = os.environ.get('EASYHYBRID_HOME')#'/home/fernando/programs/EasyHybrid3'#os.environ.get('EASYHYBRID_HOME')





class MainWindow:
    """ Class doc """
    def __init__ (self, vm_session = None):
        """ Class initialiser """
        self.home               =  EASYHYBRID_HOME
        self.EASYHYBRID_VERSION =  EASYHYBRID_VERSION
        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(self.home,'gui/MainWindow.glade'))
        self.builder.add_from_file(os.path.join(self.home,'gui/MainWindow_text_and_logs.glade'))

        self.builder.connect_signals(self)
        self.window = self.builder.get_object('window1')
        self.window.set_default_size(1200, 600)                          
        self.window.set_title('EasyHybrid {}'.format(EASYHYBRID_VERSION))                          


        self.statusbar_main = self.builder.get_object('statusbar1')
        self.statusbar_main.push(1,'Welcome to EasyHybrid version {}, a pDynamo3 graphical tool'.format(EASYHYBRID_VERSION))
        
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
        self.selection_box_frane = VismolSelectionTypeBox(vm_session = self.vm_session)
        self.vm_session.selection_box_frane = self.selection_box_frane 
        self.menu_box.add(self.selection_box_frane.box)

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
            
            
            #self.HBOX.pack_start(self.traj_frame, False, False, 1)

            #self.paned_H.add(self.notebook_H1)
            self.paned_H.add(self.HBOX)
            self.paned_H.add(self.notebook_H2)
            self.paned_H.set_position(250)

            self.paned_V.add(self.paned_H)
            
            self.bottom_notebook = BottonNoteBook(main = self)
            #self.paned_V.add(self.builder.get_object('notebook_text_and_logs'))
            self.paned_V.add(self.bottom_notebook.widget)
            self.paned_V_position = 400
            self.paned_V.set_position(self.paned_V_position)

            self.bottom_notebook.status_teeview_add_new_item(message = 'This is EasyHybrid 3.0, have a happy simulation day!')
            
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

        
        
        
        '''#- - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -#'''
        self.system_liststore        = Gtk.ListStore(str, int)
        #self.active_system_liststore = Gtk.ListStore(str, int)
        
        self.vobject_liststore_dict  = { 
                                       0 : Gtk.ListStore(str, int, int)  
                                       # sys_id : GtkListiStore 
                                       }
        '''#- - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -#'''




        '''#- - - - - - - - - - - -  pDynamo - - - - - - - - - - - - - - -#'''
        self.p_session = pDynamoSession(vm_session = vm_session)
        '''#- - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -#'''

        
        '''#- - - - - - - - - - G T K  W I N D O W S - - - - - - - - - - -#'''
        self.window_list = []
        
        self.filechooser                  = FileChooser()
        self.NewSystemWindow              = ImportANewSystemWindow       ( main = self )
        
        self.setup_QCModel_window         = EasyHybridSetupQCModelWindow ( main = self )
        self.window_list.append(self.setup_QCModel_window)
        
        self.geometry_optimization_window = GeometryOptimizatrionWindow  ( main = self )
        self.window_list.append(self.geometry_optimization_window)
        
        self.PES_scan_window              = PotentialEnergyScanWindow    ( main=  self)
        self.window_list.append(self.PES_scan_window)

        self.selection_list_window        = SelectionListWindow          ( main = self, system_liststore =  self.system_liststore)
        self.window_list.append(self.selection_list_window)

        
        self.go_to_atom_window            = EasyHybridGoToAtomWindow     ( main = self, system_liststore = self.system_liststore)
        self.window_list.append(self.go_to_atom_window)

        
        self.pDynamo_selection_window     = PDynamoSelectionWindow       ( main = self)
        
        self.export_data_window           = ExportDataWindow             ( main = self)
        self.window_list.append(self.export_data_window)

        
        #self.energy_refinement_window     = EnergyRefinementWindow       ( main = self)
        #self.window_list.append(self.energy_refinement_window)

        self.import_trajectory_window     = ImportTrajectoryWindow       (main = self)
        self.window_list.append(self.import_trajectory_window)

        
        self.PES_analysis_window          = PotentialEnergyAnalysisWindow(main = self)#, coor_liststore = self.system_liststore)
        self.window_list.append(self.PES_analysis_window)


        self.trajectory_player_window  = TrajectoryPlayerWindow (main = self)
        self.terminal_window           = TerminalWindow  (main = self)
        
        self.molecular_dynamics_window  = MolecularDynamicsWindow(main = self)
        
        self.umbrella_sampling_window   = UmbrellaSamplingWindow(main = self)
        
        self.chain_of_states_opt_window = ChainOfStatesOptWindow(main = self)
        
        self.normal_modes_analysis_window =   NormalModesAnalysisWindow (main = self)
        self.normal_modes_window          =   NormalModesWindow (main = self)
        '''#- - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -#'''

        self.window.connect("destroy", Gtk.main_quit)
        self.window.connect("check-resize", self.window_resize)
        self.window.connect("delete-event",    Gtk.main_quit)
        
        self.builder.get_object('button_test').connect("clicked",    self.run_test)
        
        self.window.show_all()

    def window_resize (self, a, b =None, c=None):
        """ Function doc """
        w, h = a.get_size()
        #self.paned_V_position += h/100
        #self.paned_V.set_position(self.paned_V_position)
        
        #print ( a.get_size(), b, c)
    
    def add_vobject_to_vobject_liststore_dict (self, vismol_object):
        """ Function doc """
        e_id = vismol_object.e_id
        vismol_object.liststore_iter = self.vobject_liststore_dict[e_id].append([vismol_object.name, 
                                                                              vismol_object.index, 
                                                                              vismol_object.e_id])

        print('\n\n\n\n')
        print(self.vobject_liststore_dict[e_id])
        print(list(self.vobject_liststore_dict[e_id]))
        print('\n\n\n\n')


    def clear_vobject_liststore_dict (self, e_id = 'all'):
        """ Function doc """
        if e_id == 'all':
            for e_id, liststore in self.vobject_liststore_dict.items(): 
                liststore.clear()
            else:
                self.vobject_liststore_dict[e_id].clear()
    
    
    def refresh_vobject_liststore_dict  (self, e_id = 'all'):
        """ Function doc """
        if e_id == 'all':
            self.clear_vobject_liststore_dict()
            for index, vobject in self.vm_session.vm_objects_dic.items():
                self.add_vobject_to_vobject_liststore_dict(vobject)
        
        else:
            self.clear_vobject_liststore_dict(e_id = e_id)
            for index, vobject in self.vm_session.vm_objects_dic.items:
                if e_id == index:
                    self.add_vobject_to_vobject_liststore_dict(vobject)
                else:
                    pass


    def on_main_toolbar_clicked (self, button):
        """ Function doc """
        if button  == self.builder.get_object('toolbutton_new_system'):
            self.NewSystemWindow.OpenWindow()
        
        if button  == self.builder.get_object('toolbutton_save'):
            self.p_session.save_easyhybrid_session( filename = 'session.easy')
        
        if button  == self.builder.get_object('toolbutton_save_as'):
            self.gtk_save_as_file (button)
            
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
            self.energy_refinement_window.OpenWindow()
            
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
            self.gtk_load_files (menuitem)
            
        elif menuitem == self.builder.get_object('menuitem_save'):
            self.gtk_save_file (menuitem)
            
        elif menuitem == self.builder.get_object('menuitem_save_as'):
            self.gtk_save_as_file (menuitem)

        elif menuitem == self.builder.get_object('menuitem_export'):
            self.export_data_window.OpenWindow()
            
        elif menuitem == self.builder.get_object('menuitem_quit'):
            ##print(menuitem, 'menu_item_merge_system')
            pass
        
        
        
        
        elif menuitem == self.builder.get_object('menuitem_energy'):
            self.gtk_get_energy(button)
            
        elif menuitem == self.builder.get_object('menuitem_geometry_optimization'):
            self.geometry_optimization_window.OpenWindow()
            
        elif menuitem == self.builder.get_object('menuitem_molecular_dynamics'):
            self.molecular_dynamics_window.OpenWindow()
            
        elif menuitem == self.builder.get_object('menuitem_normal_modes'):
            pass
            
        elif menuitem == self.builder.get_object('menuitem_rection_coordinate_scans'):
            self.PES_scan_window.OpenWindow()
            
        elif menuitem == self.builder.get_object('menuitem_nudged_elastic_band'):
            pass
            
        elif menuitem == self.builder.get_object('menuitem_umbrella_sampling'):
            self.umbrella_sampling_window.OpenWindow()
        
        elif menuitem == self.builder.get_object('menuitem_check_pDynamo_tools_bar'):
            if menuitem.get_active():
                self.builder.get_object('toolbar4_pdynamo_tools').show()
            else:
                self.builder.get_object('toolbar4_pdynamo_tools').hide()
        
        
        elif menuitem == self.builder.get_object('menuitem_check_selection_toolbar'):
            if menuitem.get_active():
                self.builder.get_object('toolbar2_selections').show()
            else:
                self.builder.get_object('toolbar2_selections').hide()

        elif menuitem == self.builder.get_object('menuitem_energy_analysis'):
            self.PES_analysis_window.OpenWindow()

        elif menuitem == self.builder.get_object('menuitem_about'):
            about_dialog = Gtk.AboutDialog()
            text = '''EasyHybrid is a free and open source graphical environment for the pDynamo3 package. It is developed using a combination of Python3, Cython3, GTK, and modern OpenGL. It utilizes the VisMol graphical engine to render 3D structures. EasyHybrid is a graphical extension of pDynamo that allows users to perform most of the basic pDynamo routines within its interface, as well as inspect, edit, and export pDynamo systems for further simulations using Python scripting in text mode. '''
            # Set the properties of the about dialog
            about_dialog.set_program_name("EasyHybrid")
            about_dialog.set_version("3.0")
            #about_dialog.set_copyright("Copyright 2022 My Company")
            #about_dialog.set_comments("A simple program for demonstrating about dialogs")
            about_dialog.set_comments(text)
            about_dialog.set_website("https://github.com/ferbachega/EasyHybrid3")
            about_dialog.set_website_label("Visit the EasyHybrid website")
            about_dialog.set_authors(["Fernando Bachega", "Igor Barden", 'Luis Timmers', "Martin Field"])
            about_dialog.set_documenters(["Fernando Bachega", "Igor Barden"])

            # Connect the "response" signal to the hide() method of the about dialog,
            # so that the dialog will be closed when the user clicks one of the buttons
            about_dialog.connect("response", lambda d, r: d.hide())

            # Show the about dialog
            about_dialog.show()


    def gtk_load_files (self, button):
        '''Easyhybrid and pkl pdynamo file search '''
        filters = []        
        ''' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '''
        
        '''
        filter = Gtk.FileFilter()  
        filter.set_name("EasyHybrid3 files - *.easy")
        filter.add_mime_type("Easy files files")
        filter.add_pattern("*.easy")
        filters.append(filter)
        #'''
        
        filter = Gtk.FileFilter()  
        filter.set_name("pDynamo3 files - *.pkl")
        filter.add_mime_type("PKL files")
        filter.add_pattern("*.pkl")
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
                #print('ehf file')            
                self.save_vismol_file = filename
                self.vm_session.load_easyhybrid_serialization_file(filename)            
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
                hamiltonian   = psystem.qcModel.hamiltonian
                n_QC_atoms    = len(list(psystem.qcState.pureQCAtoms))
                
                
                summary_items = psystem.electronicState.SummaryItems()
                
                string += 'hamiltonian: {}    QC atoms: {}    QC charge: {}    spin multiplicity {}    '.format(  hamiltonian, 
                                                                                                                  n_QC_atoms,
                                                                                                                  summary_items[1][1],
                                                                                                                  summary_items[2][1],
                                                                                                                 )
                
            n_fixed_atoms = len(psystem.e_fixed_table)
            string += 'fixed_atoms: {}    '.format(n_fixed_atoms)
            
            if psystem.mmModel:
                forceField = psystem.mmModel.forceField
                string += 'forceField: {}    '.format(forceField)
            
                if psystem.nbModel:
                    nbmodel = psystem.mmModel.forceField
                    string += 'nbModel: True    '
                    
                    summary_items = psystem.nbModel.SummaryItems()
                    
                
                else:
                    string += 'nbModel: False    '
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
    
    
    def run_test (self, widget):
        """ Function doc """
        #print('aloowww')
        #print(self.vm_session.vm_glcore.glcamera)
        print(self.vm_session.vm_glcore.glcamera.fog_end          ) #= self.z_far
        print(self.vm_session.vm_glcore.glcamera.fog_start        ) #= self.fog_end - self.min_zfar
        print('\nview matrixes: \n', self.vm_session.vm_glcore.glcamera.view_matrix      ) #= self._get_view_matrix(pos)
        print('\nprojection_matrix: \n',self.vm_session.vm_glcore.glcamera.projection_matrix) #= self._get_projection_matrix()
        print('\ncamera position: \n', self.vm_session.vm_glcore.glcamera.get_position()) #= self._get_projection_matrix()
        
        camera = self.vm_session.vm_glcore.glcamera.get_position()
        for key, vobj in self.vm_session.vm_objects_dic.items():
            print(vobj.model_mat)
            print(vobj.trans_mat)
        
        
        inv_mat = mop.get_inverse_matrix(self.vm_session.vm_glcore.glcamera.view_matrix )
        print ('inverse of nview', inv_mat)
        
        print ('inverse of nview', inv_mat[3])
        
        
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
        text+= '    <  0.000000,  0.000000,109.953541>                             \n'
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


        text+='transform { myTransforms }                                              \n'
        text+='}                                                                       \n'
        text+='                                                                        \n'
        text+='object {molecule}                                                       \n'

        pov_file.write(text)
        #print(text)
        os.system('povray +A0.3 -UV +W2218 +H1748 +Itemp.pov +Otemp2.png')

        
        

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
        
        for row in self.treestore:
            row[4] = False
        parent = self.treestore.append(None, [system.e_id, -1,str(system.e_id)+' - '+ system.label, True, True, False, False, False, 0])
        system.e_treeview_iter = parent
        
        
        '''To improve organization and accessibility, we will add a GtkListStore 
        to a dictionary that will be accessed by all windows. Each GtkListStore 
        in the dictionary will contain the vobjects for a particular system or project'''
        system.e_liststore_iter = self.main.system_liststore.append([str(system.e_id)+' - '+ system.label, system.e_id])
        self.main.vobject_liststore_dict[system.e_id] = Gtk.ListStore(str, int, int)
        
    
    def add_vismol_object_to_treeview(self, vismol_object):
        """ Function doc """
        e_id   = vismol_object.e_id
        system = self.main.p_session.get_system(e_id)
        parent = system.e_treeview_iter #self.tree_iters_dict[system.e_treeview_iter_parent_key]
        
        size   = len(vismol_object.frames)
        _iter  = self.treestore.append(parent, [e_id,  vismol_object.index , vismol_object.name, False, False , True, vismol_object.active, True, size])
        
        self.tree_iters.append(_iter)
        self.tree_iters_dict[self.tree_iters_counter] = parent
        vismol_object.e_treeview_iter_parent_key = self.tree_iters_counter
        
        vismol_object.e_treeview_iter = _iter
        
        self.tree_iters_counter += 1
        
        for row in self.treestore:
            pass

        self.expand_row(row.path, True)
        
    
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
        renderer_text = Gtk.CellRendererText()
        #renderer_text.connect("edited", self.text_edited, model, column)
        column_text = Gtk.TreeViewColumn("System", renderer_text, text=2)
        column_text.set_sizing(Gtk.TreeViewColumnSizing.AUTOSIZE)
        column_text.set_resizable(True)
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
        self.main.vm_session.vm_glcore.queue_draw()

    
    def on_cell_active_radio_toggled(self, widget, path):
        ##print (path)
        selected_path = Gtk.TreePath(path)
        
        ##print(self.treestore[path][1], path, self.treestore[path][0])
        #
        for row in self.treestore:
            ##print(row.path, selected_path)
            row[4] = row.path == selected_path
            if row.path == selected_path:
                system_e_id = row[0]
                self.main.p_session.active_id = system_e_id
                #print (system_e_id, self.main.p_session.psystem[system_e_id].label)
            else:
                pass
        #self.main.refresh_active_system_liststore ()
        self.main.refresh_main_statusbar ()
        self.main.uptade_interface_windows_and_dialogs()




    #def text_edited( self, w, row, new_text, model, column):
    #  model[row][column] = new_text

    def on_select(self, tree, path, selection):
        '''---------------------- Row information ---------------------'''
        # Get the current selected row and the model.
        model, iter = tree.get_selection().get_selected()        
        
        # Look up the current value on the selected row and get
        # a new value to change it to.
        data2 = model.get_value(iter, 2)
        data1 = model.get_value(iter, 1)
        data0 = model.get_value(iter, 0)
        print(data0, data1, data2, self.treestore[path][8])
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
                print(vobj_id, self.treestore[treeview_iter][2], 'frames', len(vobj.frames))
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
            print(index, self.treestore[treeview_iter][2], 'frames', len(vobject.frames))


    def refresh (self):
        """ Function doc """
        self.treestore.clear()
        
        for e_id in self.main.p_session.psystem.keys():
            system = self.main.p_session.psystem[e_id]
            self.add_new_system_to_treeview (system)
        
        
        
        for v_obj_index in self.main.vm_session.vm_objects_dic.keys():
            vismol_object = self.main.vm_session.vm_objects_dic[v_obj_index]
            self.add_vismol_object_to_treeview(vismol_object)

    def rename (self, name):
        selection     = self.get_selection()
        (model, iter) = selection.get_selected()
        self.key      = model.get_value(iter, 0)
        sys           = model.get_value(iter, 1)
        
        

        old_name = model.get_value(iter, 2)
        v_id     = model.get_value(iter, 1)
        e_id     = model.get_value(iter, 0)
        if v_id == -1:
            #rename  system
            self.treestore[iter][2] = str(e_id)+' - '+ name
            self.main.p_session.psystem[e_id].label  = name
            liststore_iter = self.main.p_session.psystem[e_id].e_liststore_iter
            self.main.system_liststore[liststore_iter][0] = str(e_id)+' - '+ name
        
        else:
            self.treestore[iter][2] = name
            self.main.vm_session.vm_objects_dic[v_id].name = name
            self.main.vobject_liststore_dict[e_id][self.main.vm_session.vm_objects_dic[v_id].liststore_iter][0] = name
            
        #self.p_session.psystem[self.e_id].e_selections[new_name] = self.p_session.psystem[self.e_id].e_selections[self.key]
        #self.p_session.psystem[self.e_id].e_selections.pop(self.key)
        #self.sele_window.update_window()


class TreeViewMenu:
    """ Class doc """
    
    def __init__ (self, treeview):
        """ Class initialiser """
        self.treeview = treeview
        self.main     = treeview.main 
        #self.main     = self.treeview.main
        self.filechooser   = FileChooser()
        self.rename_window_visible =  False
        
        vobject_menu_items = {
                                'header'                : None    ,

                                '_separator'            : ''      ,
                                'Info'                  : self.f2 ,
                                'Rename Object'                : self.menu_rename ,
                                '_separator'            : ''      ,

                                'Load Data Into System' : self.load_data_to_a_system ,
                                'Change Color'  : self.change_system_color_palette ,
                                'Edit Parameters'       : self.f2 ,
                                '_separator'            : ''      ,

                                'Export As...'          : self.menu_export_data_window ,
                                'Merge System With...'  : self.f3 ,
                                '_separator'            : ''      ,

                                'Delete'                : self._menu_delete_vm_object ,
                                #'test'  : self.f1 ,
                                #'f1'    : self.f1 ,
                                #'f2'    : self.f2 ,
                                #'gordão': self.f3 ,
                                #'delete': self.f3 ,
                                }

        system_menu_items = {
                                
                                'header'                : None    ,
                                '_separator'            : ''      ,
                                'Rename System'        : self.menu_rename ,
                                '_separator'            : ''      ,
                                'Info'                  : self.f2 ,
                                'Load Data Into System' : self.load_data_to_a_system ,
                                'Change Color '  : self.change_system_color_palette ,
                                'Edit Parameters'       : self.f2 ,
                                'Export As...'          : self.menu_export_data_window ,
                                'Merge System With...'  : self.f3 ,
                                'Delete'                : self._menu_delete_system ,
                                #'test'  : self.f1 ,
                                #'f1'    : self.f1 ,
                                #'f2'    : self.f2 ,
                                #'gordão': self.f3 ,
                                #'delete': self.f3 ,
                                }
                    
                    
                    
        self.tree_view_vobj_menu  , self.tree_header_vobj_menu    = self.build_tree_view_menu(vobject_menu_items)
        self.tree_view_sys_menu   , self.tree_header_sys_menu     = self.build_tree_view_menu(system_menu_items)


    def menu_export_data_window (self,vobject = None ):
        """ Function doc """
        self.treeview.main.export_data_window.OpenWindow()
    
    
    def load_data_to_a_system (self, vobject = None ):
        """ Function doc """
        selection        = self.treeview.get_selection()
        model, iter      = selection.get_selected()
        #print (list(model))
        self.treeview.main.import_trajectory_window.OpenWindow(sys_selected = model.get_value(iter, 0))

    def change_system_color_palette (self, widget):
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
            #print('self.selectedID',self.selectedID)
            print('self.selectedID',type(self.main.p_session.psystem[self.selectedID].e_color_palette ))
            
            self.main.p_session.psystem[self.selectedID].e_color_palette['C'] = new_color

            #self.set_color(color =new_color)
    
    def f2 (self, vobject = None):
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

    
    
    def menu_rename (self, menu_item = None ):
        """  
        menu_item = Gtk.MenuItem object at 0x7fbdcc035700 (GtkMenuItem at 0x37cf6c0)
        
        """
        if self.rename_window_visible:
            pass
        else:
            #
            #self.e_id     = self.sele_window.system_names_combo.get_active()
            #self.e_id     = system_e_id
            

            
            #print(self.key, self.e_id)
            self.window = Gtk.Window()
            self.window.connect('destroy', self.destroy)
            self.window.set_keep_above(True)
            self.entry  = Gtk.Entry()
            
            self.entry.connect('activate', self.get_new_name)
            self.window.add(self.entry)
            self.rename_window_visible = True
            self.window.show_all()
            print(menu_item)

    def get_new_name (self, menu_item):
        """ Function doc """
        print(self.entry.get_text())
        
        new_name = self.entry.get_text()
        self.treeview.rename(new_name)
    
        self.window.destroy()
        self.rename_window_visible = False
    
    def destroy (self, widget):
        """ Function doc """
        self.rename_window_visible = False
    
    def _menu_delete_vm_object (self, widget):
        """ Function doc """
        self.delete_vm_object ( vm_object_index = self.vobject_index)
    
    
    def _menu_delete_system (self, widget):
        """ Function doc """
        self.delete_system (system_e_id = self.system_e_id )
    
    
    def delete_vm_object (self, vm_object_index = None):
        """ 
        
        vm_object_index = is the access key to the object. You can get it from vobject.index
        
        '''When an object is removed it has to be removed from the treeview and 
        vobject_liststore_dict, in addition to the vm_object_dic in the .vm_session.'''
        
        """
        if vm_object_index != None:
            #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            vobject = self.treeview.main.vm_session.vm_objects_dic[vm_object_index]
            #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
            
            #  - - - - - - REMOVING vobj FROM  vobject_liststore_dict - - - - - - -
            self.treeview.main.vobject_liststore_dict[vobject.e_id].remove(vobject.liststore_iter)
            #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
            #  - - - - - - REMOVING vobj FROM  treestore - - - - - - -
            self.treeview.treestore.remove(vobject.e_treeview_iter)
            #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
            #  - - - - - - - - REMOVING vobj FROM vm_object_dic - - - - - - - - - -
            self.treeview.main.vm_session.vm_objects_dic[vm_object_index] = None
            self.treeview.main.vm_session.vm_objects_dic.pop(vm_object_index)# = None
            #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            self.treeview.main.vm_session.vm_glcore.queue_draw()
            
        #self.treeview.main.refresh_active_system_liststore()
        
        #if self.treeview.main.selection_list_window.visible:
        #    self.treeview.main.selection_list_window.update_window(system_names = True, coordinates = True,  selections = True)
    
    
    def delete_system (self, system_e_id = None ):
        """ 
        system_e_id = is the access key to the object. You can get it from vobject.index

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
            for index, vobject in self.treeview.main.vm_session.vm_objects_dic.items():
                if vobject.e_id == system_e_id:                   
                    self.treeview.treestore.remove(vobject.e_treeview_iter)
                    pop_list.append(index)
            '''removing from vm_object_dic'''
            for index in pop_list:
                self.treeview.main.vm_session.vm_objects_dic.pop(index)
            
            '''removing vobject from vobject_liststore_dict'''
            self.treeview.main.vobject_liststore_dict.pop(system_e_id)
            
            
            #  - - - - - - - - removing system treeview items - - - - - - - - - - -
            system = self.treeview.main.p_session.get_system(index = system_e_id)
            self.treeview.main.system_liststore.remove(system.e_liststore_iter)
            self.treeview.treestore.remove(system.e_treeview_iter)
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
            a = self.treeview.main.p_session.delete_system(system_e_id)
            self.treeview.main.vm_session.vm_glcore.queue_draw()
            
        
        #self.treeview.refresh()
        #self.treeview.main.refresh_vobject_liststore_dict()
        
        #if self.treeview.main.selection_list_window.visible:
        #    self.treeview.main.selection_list_window.update_window(system_names = True, coordinates = True,  selections = True)


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
                

class BottonNoteBook:
    
    def __init__ (self, main):
        """ Function doc """
        self.home    = main.home
        self.main    = main
        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(self.home,'gui/MainWindow_text_and_logs.glade'))
        
        self.widget = self.builder.get_object('notebook_text_and_logs')

        self._define_status_treeview()

    
    def _define_status_treeview (self):
        """ Function doc """
        self.treeview = self.builder.get_object('treeview_status')
        self.status_liststore = Gtk.ListStore(str, # time 
                                              str, # message
                                              str, # logfile_path
                                              
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
        renderer_text = Gtk.CellRendererText()
        column_text = Gtk.TreeViewColumn("message", renderer_text, text=1, visible = True)
        self.treeview.append_column(column_text)        
                
        #self.treeview.connect('row-activated', self.on_select)
        #self.treeview.connect('button-release-event', self.on_treeview_mouse_button_release_event )
        self.treeview.set_headers_visible(False)

    def status_teeview_add_new_item (self, message = None , logfile = ''):
        """ Function doc """
        current_time   = time.time()
        formatted_time = time.strftime("%Y-%m-%d   %H:%M:%S", time.localtime(current_time))
        
        self.status_liststore.append([formatted_time, message, logfile])



def main():
    logging.basicConfig(format="%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s",
                        datefmt="%Y-%m-%d:%H:%M:%S", level=logging.DEBUG)
    
    vm_session = EasyHybridSession( a = True)
    vm_session.vm_widget.insert_glmenu()
    
    main_window = MainWindow(vm_session = vm_session)
    #window = Gtk.Window(title="Vismol window")
    #container = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
    #container.pack_start(vm_session.vm_widget, True, True, 0)
    #window.connect("key-press-event", vm_session.vm_widget.key_pressed)
    #window.connect("key-release-event", vm_session.vm_widget.key_released)
    #window.add(container)
    #window.connect("delete-event", Gtk.main_quit)
    #window.show_all()
    try:
        filein = sys.argv[-1]
        vm_session.load_molecule(filein)
    except:
        pass
    Gtk.main()
    return 0

if __name__ == "__main__":
    main()

