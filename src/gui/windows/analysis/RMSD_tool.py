#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  easyhybrid_pDynamo_selection.py
#  
#  Copyright 2022 Fernando <fernando@winter>
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
import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
import gc
import os

#import math
from gi.repository import Gtk, Pango
from gi.repository import Gdk
from gi.repository import GdkPixbuf
#from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.widgets.custom_widgets import FolderChooserButton

from gui.windows.setup.windows_and_dialogs import TextWindow

#from util.geometric_analysis import get_simple_distance
#from util.geometric_analysis import get_simple_angle
#from util.geometric_analysis import get_simple_dihedral

from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox


from pprint import pprint
import numpy as np

class RMSDToolWindow:

    def __init__(self, main = None , system_liststore = None):
        """ Class initialiser """
        self.main                = main
        self.home                = main.home
        self.p_session           = self.main.p_session
        self.vm_session          = self.main.vm_session
        self.Visible             = False  
        
        self.system_liststore      = system_liststore
        self.coordinates_liststore = Gtk.ListStore(str, int)
        
        self.parameters = []
        

    def OpenWindow (self, vobject = None):
        """ Function doc """
        if self.Visible  ==  False:

            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.main.home,'src/gui/windows/analysis/RMSD_tool.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('window')
            self.window.set_title('RMSD Analysis')
            #self.window.set_keep_above(True)            
            #self.window.set_default_size(100, 500)
            
            self.box_system      = self.builder.get_object('box_system')
            self.box_coordinates = self.builder.get_object('box_coordinates')
            
            self.combobox_systems     = SystemComboBox     (self.main)
            self.combobox_coordinates = CoordinatesComboBox(self.main.vobject_liststore_dict[self.main.p_session.active_id])
            
            self.combobox_systems.connect("changed", self.on_combobox_systems_changed)
            
            
            self.box_system.pack_start       (self.combobox_systems, False, False, 0)
            self.box_coordinates.pack_start  (self.combobox_coordinates, False, False, 0)
            
            self.build_treeview()
            '''
            #------------------------------------------------------------------------
            # define_measurement
            #------------------------------------------------------------------------
            self.btn_distance = self.builder.get_object('btn_distance')
            self.btn_distance.connect("button_press_event"  ,self.define_measurement) 
            self.btn_distance = self.builder.get_object('btn_angle')
            self.btn_distance.connect("button_press_event"  ,self.define_measurement) 
            self.btn_distance = self.builder.get_object('btn_dihedral')
            self.btn_distance.connect("button_press_event"  ,self.define_measurement) 
            #------------------------------------------------------------------------
            self.btn_cancel = self.builder.get_object('btn_cancel')
            self.btn_cancel.connect("button_press_event"  ,self.CloseWindow)
            #------------------------------------------------------------------------
            #------------------------------------------------------------------------
            self.btn_clear = self.builder.get_object('btn_clear')
            self.btn_clear.connect("button_press_event"  ,self.on_btn_clear)
            #------------------------------------------------------------------------
            
            
            self.folder_chooser_button = FolderChooserButton(main =  self.main, sel_type = 'folder', home =  self.home)
            self.builder.get_object('folder_chooser_box').pack_start(self.folder_chooser_button.btn, True, True, 0)
            system_id      = self.p_session.active_id

            
            self.build_treeview()
            '''
            
            self.window.show_all()
            #self.column_type.set_spacing(40)
            #self.column_1   .set_spacing(40)
            #self.column_2   .set_spacing(40)
            #self.column_3   .set_spacing(40)
            #self.column_4   .set_spacing(40)
            
            self.Visible  = True
            
      
    def on_combobox_systems_changed (self, widget):
        """ Function doc """
        cb_id = widget.get_system_id()
        
        if cb_id is not None:
            
            self.update_window (coordinates = True)
            
            key =  self.combobox_coordinates.get_vobject_id()
            #_, key = self.coordinates_liststore[cb_id]
            
            self.VObj = self.vm_session.vm_objects_dic[key]

            
    def update_window (self, system_names = False, coordinates = True,  selections = False ):
        """ Function doc """

        if self.Visible:
            
            _id = self.combobox_systems.get_active()
            if _id == -1:
                '''_id = -1 means no item inside the combobox'''
                return None
            else:    
                _, system_id, pixbuf = self.system_liststore[_id]
            
            if system_names:
                self.refresh_system_liststore ()
                self.combobox_systems.set_active(_id)
            
            if coordinates:
                self.refresh_coordinates_liststore ()
    

    def refresh_coordinates_liststore(self, system_id = None):
        """ Function doc """
        system_id = self.combobox_systems.get_system_id()
        #print(2313, system_id,self.main.vobject_liststore_dict )
        self.combobox_coordinates.set_model(self.main.vobject_liststore_dict[system_id])
        self.combobox_coordinates.set_active_vobject(-1)
            

    def build_treeview (self):
        """ Function doc """
       
        self.liststore = Gtk.ListStore(str  , #0 system_e_id           
                                       str  , #1 vobject 
                                       str  , #2 vobject 
                                       str  , #3 name 
                                       str  , #4 vobject 
                                       str  , #5 name 
                                       )
        
 
    
        
        self.treeview = self.builder.get_object('rmsd_treeview')
        #self.treeview.set_spacing(40)
        #connect("button_press_event"  , seqview.button_press)
        
        ## column
        #renderer_radio = Gtk.CellRendererToggle()
        #renderer_radio.set_radio(True)
        #renderer_radio.connect("toggled", self.on_cell_active_radio_toggled)
        #column_radio = Gtk.TreeViewColumn("A", renderer_radio, active=0)
        #self.treeview.append_column(column_radio)

        # column
        renderer_text = Gtk.CellRendererText()
        column_text = Gtk.TreeViewColumn("Run", renderer_text, text=0)
        #column_text.set_sizing(Gtk.TreeViewColumnSizing.AUTOSIZE)
        column_text.set_resizable(True)
        column_text.set_spacing(40)
        self.treeview.append_column(column_text)        
        self.column_type = column_text
        
        # column
        renderer_text = Gtk.CellRendererText()
        column_text = Gtk.TreeViewColumn("vobject", renderer_text, text=1)
        column_text.set_sizing(Gtk.TreeViewColumnSizing.AUTOSIZE)
        column_text.set_resizable(True)
        #column_text.set_spacing(40)
        self.treeview.append_column(column_text)        
        self.column_1 = column_text

        # column
        renderer_text = Gtk.CellRendererText()
        column_text = Gtk.TreeViewColumn("min", renderer_text, text=2)
        column_text.set_sizing(Gtk.TreeViewColumnSizing.AUTOSIZE)
        column_text.set_resizable(True)
        #column_text.set_spacing(40)
        self.treeview.append_column(column_text)  
        self.column_2 = column_text
        
        # column
        renderer_text = Gtk.CellRendererText()
        column_text = Gtk.TreeViewColumn("max", renderer_text, text=3)
        column_text.set_sizing(Gtk.TreeViewColumnSizing.AUTOSIZE)
        column_text.set_resizable(True)
        self.treeview.append_column(column_text)  
        #column_text.set_spacing(40)
        self.column_3 = column_text
        
        # column
        renderer_text = Gtk.CellRendererText()
        column_text = Gtk.TreeViewColumn("average", renderer_text, text=4)
        column_text.set_sizing(Gtk.TreeViewColumnSizing.AUTOSIZE)
        column_text.set_resizable(True)
        column_text.set_spacing(40)
        self.treeview.append_column(column_text)  
        self.column_4 = column_text



        ## column
        #renderer_text = Gtk.CellRendererText()
        #column_text = Gtk.TreeViewColumn("Atoms", renderer_text, text=2)
        #column_text.set_sizing(Gtk.TreeViewColumnSizing.AUTOSIZE)
        #column_text.set_resizable(True)
        #self.treeview.append_column(column_text)  

        self.treeview.set_model(self.liststore)
        self.treeview.connect('button-release-event', self.on_treeview_mouse_button_release_event )
        self.treeview.connect("key-press-event", self.on_key_press)


    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False
        self.on_btn_clear(None, None)


    def define_measurement(self, widget, event):
        """ Function doc """
        
        atom1 = self.vm_session.picking_selections.picking_selections_list[0]
        atom2 = self.vm_session.picking_selections.picking_selections_list[1]
        atom3 = self.vm_session.picking_selections.picking_selections_list[2]
        atom4 = self.vm_session.picking_selections.picking_selections_list[3]

        #print(atom1, atom2, atom3, atom4)
        
        atom_list = [atom1, atom2, atom3, atom4]
        
        
        if widget ==  self.builder.get_object('btn_distance'):
            
            liststore_item = ['distance']
            
            if atom1 and atom2:
                
                idx1 = str(atom1.index)
                idx2 = str(atom2.index)
                
                sym1 = atom1.symbol
                sym2 = atom2.symbol
                
                self.parameters.append(['distance',atom1, atom2])
                
                item1 = '  {} / {}  '.format(idx1, sym1)
                liststore_item.append(item1)
                item2 = '  {} / {}  '.format(idx2, sym2)
                liststore_item.append(item2)
                
                liststore_item.append('-')
                liststore_item.append('-')
                liststore_item.append('-')
                self.liststore.append(liststore_item)
                
        if widget ==  self.builder.get_object('btn_angle'):
            
            liststore_item = ['angle']
            if atom1 and atom2 and atom3:
                
                idx1 = str(atom1.index)
                idx2 = str(atom2.index)
                idx3 = str(atom3.index)
                
                sym1 = atom1.symbol
                sym2 = atom2.symbol
                sym3 = atom3.symbol
                
                self.parameters.append(['angle',atom1, atom2, atom3])
                
                item1 = '  {} / {}  '.format(idx1, sym1)
                liststore_item.append(item1)
                
                item2 = '  {} / {}  '.format(idx2, sym2)
                liststore_item.append(item2)
                
                item3 = '  {} / {}  '.format(idx3, sym3)
                liststore_item.append(item3)
                
                #liststore_item.append('-')
                liststore_item.append('-')
                liststore_item.append('-')
                self.liststore.append(liststore_item)
                
        if widget ==  self.builder.get_object('btn_dihedral'):
            
            liststore_item = ['dihedral']
            if atom1 and atom2 and atom3 and atom4:
                
                idx1 = str(atom1.index)
                idx2 = str(atom2.index)
                idx3 = str(atom3.index)
                idx4 = str(atom4.index)
                
                sym1 = atom1.symbol
                sym2 = atom2.symbol
                sym3 = atom3.symbol
                sym4 = atom4.symbol
                
                self.parameters.append(['dihedral',atom1, atom2, atom3, atom4])
                
                item1 = '  {} / {}  '.format(idx1, sym1)
                liststore_item.append(item1)
                
                item2 = '  {} / {}  '.format(idx2, sym2)
                liststore_item.append(item2)
                
                item3 = '  {} / {}  '.format(idx3, sym3)
                liststore_item.append(item3)
                
                item4 = '  {} / {}  '.format(idx4, sym4)
                liststore_item.append(item4)
                
                liststore_item.append('-')
                #liststore_item.append('-')
                self.liststore.append(liststore_item)
    

    def on_select(self, tree, path, selection):
        '''---------------------- Row information ---------------------'''
        # Get the current selected row and the model.
        model, iter = tree.get_selection().get_selected()        
        
        # Look up the current value on the selected row and get
        # a new value to change it to.
        data2 = model.get_value(iter, 2)
        data1 = model.get_value(iter, 1)
        data0 = model.get_value(iter, 0)
        print (data0,data1,data2)


    def on_treeview_mouse_button_release_event (self, tree, event):
        """ Function doc """
        
        model, iter = tree.get_selection().get_selected()        
        selection   = tree.get_selection()
        model       = tree.get_model()
        #try:
        #    self.vm_object_index  = int(model.get_value(iter, 1))
        #    self.system_e_id      = int(model.get_value(iter, 0))
        #except:
        #    return False
        
        
        if event.button == 3:
            #self.treeview_menu.open_menu(self.system_e_id, self.vm_object_index)
            print('button == 3')
        if event.button == 2:
            #print('event.button 2',self.vm_object_index)
            #if self.vm_object_index == -1:
            #    #means that is not a vismol object
            #    pass
            #else:
            #    vismol_object = self.main.vm_session.vm_objects_dic[self.vm_object_index]
            #    self.main.vm_session.vm_glcore.center_on_coordinates(vismol_object, vismol_object.mass_center)
            print('button == 2')
        if event.button == 1:
            print('button == 1')


    def on_btn_clear (self, widget, event ):
        """ Function doc """
        self.liststore.clear()
        self.parameters = []
    
    
    def on_key_press(self, widget, event):
        """ Deletes the key pressed and deletes an item if it is 'Delete' """
        if event.keyval == Gdk.KEY_Delete:  # Verifica se a tecla pressionada é Delete
            selection = self.treeview.get_selection()
            model, treeiter = selection.get_selected()
            
            if treeiter is not None:
                path = model.get_path(treeiter)
                
                deleted_line = path.get_indices()[0] # this is an int
                
                self.parameters.pop(deleted_line)
                
                print(f"Selected line to be deleted: {path.get_indices()[0]}")
                model.remove(treeiter)  # Remove the selected line


    def on_btn_run_analysis (self, widget, event = None):
        """ Function doc """

        measures = []
        jobsize   = len(self.parameters)
        print('\n\n')
        
        for job in self.parameters:
            
            vobject      = job[1].vm_object
            size         = len(vobject.frames)
            job_distance = []
            
            for i in range(size):
                if job[0] == 'distance':    
                    a1_coord = job[1].coords(frame=i)
                    a2_coord = job[2].coords(frame=i)
                    
                    dist1 = get_simple_distance(a1_coord, a2_coord)            
                    job_distance.append(dist1)
                
                elif job[0] == 'angle':
                    a1_coord = job[1].coords(frame=i)
                    a2_coord = job[2].coords(frame=i)
                    a3_coord = job[3].coords(frame=i)
                    
                    angle = get_simple_angle( a1_coord, a2_coord, a3_coord)            
                    job_distance.append(angle)
                
                
                elif  job[0] == 'dihedral':
                    a1_coord = job[1].coords(frame=i)
                    a2_coord = job[2].coords(frame=i)
                    a3_coord = job[3].coords(frame=i)
                    a4_coord = job[4].coords(frame=i)
                    
                    dangle = get_simple_dihedral( a1_coord, a2_coord, a3_coord, a4_coord)
                    job_distance.append(dangle)
                
                else:
                    pass
                    
            measures.append(job_distance)
        
        #print(measures)
        '''
        # printing the results 
        '''
        text = ''
        for i in range(size):
            text += ' {:<4} '.format(i)
            for j in range(jobsize):
                text += ' {:^.7f} '.format(measures[j][i])
            
            text += '\n' 
        
        print(text)
        textwindow = TextWindow(text)

    
    
    

