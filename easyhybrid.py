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
import gi, sys, os
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
from vismol.core.vismol_session import VismolSession

from gEngine.eSession import EasyHybridSession
from gui.gtk_widgets import VismolSelectionTypeBox
from gui.gtk_widgets import FileChooser

from gui.windows.windows_and_dialogs import ImportANewSystemWindow

from pdynamo.pDynamo2EasyHybrid import pDynamoSession

logger = logging.getLogger(__name__)
EASYHYBRID_VERSION = '3.0'
VISMOL_HOME = '/home/fernando/programs/EasyHybrid3'#os.environ.get('VISMOL_HOME')

class MainWindow:
    """ Class doc """
    def __init__ (self, vm_session = None):
        """ Class initialiser """
        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(VISMOL_HOME,'gui/MainWindow.glade'))
        #self.builder.add_from_file(os.path.join(VISMOL_HOME,'GTKGUI/toolbar_builder.glade'))
        self.builder.connect_signals(self)
        self.window = self.builder.get_object('window1')
        self.window.set_default_size(1200, 600)                          
        self.window.set_title('EasyHybrid {}'.format(EASYHYBRID_VERSION))                          


        self.statusbar_main = self.builder.get_object('statusbar1')
        #self.statusbar_main.push(0,'wellcome to EasyHydrid')
        self.statusbar_main.push(1,'Welcome to EasyHybrid version {}, a pDynamo3 graphical tool'.format(EASYHYBRID_VERSION))
        
        self.paned_V = self.builder.get_object('paned_V')
        #self.nootbook  =  self.builder.get_object('notebook2')
        #self.window = Gtk.Window(title="VisMol window")
        #self.main_box = Gtk.Box()
        self.vm_session = vm_session#( main_session = None)
        self.vm_session.main_session = self
        
        self.window.connect("key-press-event",   self.vm_session.vm_widget.key_pressed)
        self.window.connect("key-release-event", self.vm_session.vm_widget.key_released)
        
                
        
        self.menu_box = self.builder.get_object('toolbutton_selection_box')
        self.box2 = self.builder.get_object('box2')
        
        #self.selection_box = self.vm_session.selection_box
        self.selection_box = VismolSelectionTypeBox(vm_session = self.vm_session)
        self.menu_box.add(self.selection_box.box)

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
        

        #self.treeview = GtkEasyHybridMainTreeView(self, vm_session)
        #self.ScrolledWindow_notebook_H1.add(self.treeview)

        self.notebook_H1.append_page(self.ScrolledWindow_notebook_H1, Gtk.Label('Objects'))
        


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
        self.vm_session = vm_session#( main_session = None)
        #self.filechooser   = FileChooser()
        #-------------------------------------------------------------------
        
        
        self.container = Gtk.Box (orientation = Gtk.Orientation.VERTICAL)
        self.command_line_entry = Gtk.Entry()

        
        if self.vm_session is not None:
            #player

            self.container.pack_start(self.vm_session.vm_widget, True, True, 0)
            
            #self.traj_frame = self.vm_session.trajectory_frame
            #self.container.pack_start(self.traj_frame, False, False, 1)
            #self.container.pack_start(self.command_line_entry, False, False, 0)

            self.notebook_H2.append_page(self.container, Gtk.Label('view'))
            self.notebook_H2.append_page(Gtk.TextView(), Gtk.Label('logs'))
            
            
            #self.HBOX = Gtk.Box(orientation = Gtk.Orientation.VERTICAL, spacing = 6)
            self.HBOX = Gtk.Box(orientation = Gtk.Orientation.VERTICAL, spacing = 0)
            self.HBOX.pack_start(self.notebook_H1, True, True, 0)
            #self.HBOX.pack_start(self.traj_frame, False, False, 1)

            #self.paned_H.add(self.notebook_H1)
            self.paned_H.add(self.HBOX)
            self.paned_H.add(self.notebook_H2)

            self.paned_V.add(self.paned_H)
            #self.paned_V.add(Gtk.TextView())
            
            #self.paned_V.add(self.traj_frame)
        

        #self.player_frame = self.vm_session.player_frame
        #self.player_frame.show_all()
        '''#- - - - - - - - - - G T K  W I N D O W S - - - - - - - - - - -#'''

        self.filechooser                  = FileChooser()
        self.NewSystemWindow              = ImportANewSystemWindow(main = self, home  = VISMOL_HOME)
        '''#- - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -#'''



        '''#- - - - - - - - - - - -  pDynamo - - - - - - - - - - - - - - -#'''
        self.p_session = pDynamoSession(vm_session = vm_session)
        '''#- - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - -#'''


        self.window.connect("destroy", Gtk.main_quit)
        self.window.connect("delete-event",    Gtk.main_quit)
        self.window.show_all()


    def on_main_toolbar_clicked (self, button):
        """ Function doc """
        if button  == self.builder.get_object('toolbutton_new_system'):
            self.NewSystemWindow.OpenWindow()
        
        if button  == self.builder.get_object('toolbutton_save'):
            self.gtk_save_file (button)
             
        if button  == self.builder.get_object('toolbutton_save_as'):
            self.gtk_save_as_file (button)
            
        if button == self.builder.get_object('toolbutton_terminal'):
            if button.get_active ():
                self.terminal_window.OpenWindow()
            else:
                self.terminal_window.CloseWindow(None, None)
        
        if button == self.builder.get_object('toolbutton_trajectory_tool1'):
            if button.get_active ():
                self.traj_frame.hide()
                #window = Gtk.Window()
                #window.add(self.traj_frame)
                #window.show_all()
                #self.traj_frame
                #print('ativo')
            else:
                print('desativo')
                self.traj_frame.show()
        if button == self.builder.get_object('button_go_to_atom'):
            self.treeview.main_session.go_to_atom_window.OpenWindow()

        if button  == self.builder.get_object('selections'):
            #print('OpenWindow')
            self.selection_list_window.OpenWindow()

        if button  == self.builder.get_object('toolbutton_energy'):
            self.energy_refinement_window.OpenWindow()
            #self.gtk_get_energy(button)
            
        if button  == self.builder.get_object('toolbutton_setup_QCModel'):
            #self.dialog_import_a_new_systen = EasyHybridImportANewSystemDialog(self.p_session, self)
            #self.dialog_import_a_new_systen.run()
            #self.dialog_import_a_new_systen.hide()
            #self.NewSystemWindow.OpenWindow()
            self.setup_QCModel_window.OpenWindow()
        
        if button  == self.builder.get_object('toolbutton_system_check'): 
            self.p_session.systems[self.p_session.active_id]['vobject'].get_backbone_indexes ()
            print(self.p_session.systems[self.p_session.active_id]['vobject'].c_alpha_bonds)          
            print(self.p_session.systems[self.p_session.active_id]['vobject'].c_alpha_atoms)
        
        
        if button  == self.builder.get_object('toolbutton_geometry_optimization'):
            self.geometry_optimization_window.OpenWindow()
        
        if button  == self.builder.get_object('toolbutton_pDynamo_selections'):
            #print('toolbutton_pDynamo_selections')
            #print('self.p_session.picking_selections_list', self.vm_session.picking_selections.picking_selections_list)
            self.pDynamo_selection_window.OpenWindow()
            
            '''
            atom1 = self.vm_session.picking_selections.picking_selections_list[0]
            #print (atom1.chain, atom1.resn, atom1.resi, atom1.name)
            print ("%s:%s.%s:%s" %(atom1.chain, atom1.resn, atom1.resi, atom1.name))
            
            _centerAtom ="%s:%s.%s:%s" %(atom1.chain, atom1.resn, atom1.resi, atom1.name)
            _radius =  10.0
            
            self.p_session.selections (_centerAtom, _radius)
            '''
       
        if button  == self.builder.get_object('toolbutton_pes_scan'):
            self.PES_scan_window.OpenWindow()
        
        if button  == self.builder.get_object('toolbutton_molecular_dynamics'):
            self.molecular_dynamics_window.OpenWindow()

        if button  == self.builder.get_object('toolbutton_umbrella_sampling'):
            #print('toolbutton_umbrella_sampling')
            self.umbrella_sampling_window.OpenWindow()


def main():
    logging.basicConfig(format="%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s",
                        datefmt="%Y-%m-%d:%H:%M:%S", level=logging.DEBUG)
    # logging.basicConfig(level=logging.DEBUG)
    
    #vm_session = VismolSession(toolkit="Gtk_3.0")
    
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

