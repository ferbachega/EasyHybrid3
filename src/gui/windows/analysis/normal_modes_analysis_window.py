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
#from GTKGUI.gtkWidgets.filechooser import FileChooser
#from easyhybrid.pDynamoMethods.pDynamo2Vismol import *
import gc
import os
import threading
import time

from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox


VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')


class NormalModesAnalysisWindow(Gtk.Window):
    """ Class doc """
    def __init__(self, main = None, system_liststore = None ):
        """ Class initialiser """
        self.main          = main#self.main.system_liststore
        self.home                  = main.home
        self.visible               = False        
        self.p_session             = main.p_session
        self.vm_session            = main.vm_session
        
        self.modes = None
        self.running = False
        self.stop_thread = True
        
        self.lower = 0
        self.upper = 0
        
        self.selection_liststore       = Gtk.ListStore(str, int)
        self.selection_liststore_dict  =   {
                                           # system_e_id : Gtk.ListStore(str, int)
                                           }
    
    def OpenWindow (self):
        """  """
        if self.visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'src/gui/windows/analysis/normal_modes_analysis_window.glade')) #/home/fernando/programs/EasyHybrid3/gui/windows/normal_mode_analysis.glade
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('window')
            self.window.set_default_size(200, 600)  
            self.window.set_title('Normal Modes')  
            self.window.set_keep_above(True)
            
            
            
            
            # - - - - - - - systems combobox - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            self.box1 = self.builder.get_object('box_system')
            self.system_names_combo = SystemComboBox (self.main)
            #self.system_names_combo =self.builder.get_object('systems_combobox')
            #self.system_names_combo.set_model(self.main.system_liststore)

            self.system_names_combo.connect("changed", self.on_system_names_combobox_changed)
            #renderer_text = Gtk.CellRendererText()
            #self.system_names_combo.pack_start(renderer_text, True)
            #self.system_names_combo.add_attribute(renderer_text, "text", 0)
            '''--------------------------------------------------------------------------------------------'''
            self.box1.pack_start(self.system_names_combo, False, False, 0)


            self.coordinates_liststore = Gtk.ListStore(str, int, int) 
            self.coordinates_liststore.append(['New Visual Object', -1 , -1])
            # - - - - - - - coordinates combobox - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            self.coordinates_combobox =self.builder.get_object('coordinates_combobox')
            renderer_text2 = Gtk.CellRendererText()
            self.coordinates_combobox.pack_start(renderer_text2, True)
            self.coordinates_combobox.add_attribute(renderer_text2, "text", 0)
            self.coordinates_combobox.connect("changed", self.on_coordinates_combobox_changed)
            '''--------------------------------------------------------------------------------------------'''
            
            
            system  = self.main.p_session.get_system()
            for row in self.main.vobject_liststore_dict[self.main.p_session.active_id]:
                
                a = list (row)
                self.coordinates_liststore.append([a[0], a[1], a[2]])

            
            self.coordinates_combobox.set_model(self.coordinates_liststore)
            




            self.liststore = Gtk.ListStore(bool , #0 system_e_id           
                                           str  , #1 vobject 
                                           str  , #2 name 
                    
                                           )
            
            #for i in range(0,10):
            #    self.liststore.append([False, str(i), str(i**3)])
        
            
            self.treeview = self.builder.get_object('selection_treeview')



            # column
            renderer_radio = Gtk.CellRendererToggle()
            renderer_radio.set_radio(True)
            renderer_radio.connect("toggled", self.on_cell_active_radio_toggled)
            column_radio = Gtk.TreeViewColumn("A", renderer_radio, active=0)
            self.treeview.append_column(column_radio)

            # column
            renderer_text = Gtk.CellRendererText()
            column_text = Gtk.TreeViewColumn("Mode", renderer_text, text=1)
            column_text.set_sizing(Gtk.TreeViewColumnSizing.AUTOSIZE)
            column_text.set_resizable(True)
            self.treeview.append_column(column_text)        



            # column
            renderer_text = Gtk.CellRendererText()
            column_text = Gtk.TreeViewColumn("Frequency", renderer_text, text=2)
            column_text.set_sizing(Gtk.TreeViewColumnSizing.AUTOSIZE)
            column_text.set_resizable(True)
            self.treeview.append_column(column_text)  

            self.treeview.set_model(self.liststore)
                    
            #self.connect('row-activated', self.on_select)
            #self.connect('button-release-event', self.on_treeview_mouse_button_release_event )






            self.spin_button =self.builder.get_object('spin_button')






            # - - - - - - - Selection Treeview - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            '''
            self.treeview = self.builder.get_object('selection_treeview')
            for i, column_title in enumerate(['Selection',"Number of Atoms"]):
                renderer = Gtk.CellRendererText()
                column = Gtk.TreeViewColumn(column_title, renderer, text=i)
                self.treeview.append_column(column)
            self.treeview.set_model(self.selection_liststore)

            self.treeview.connect('row_activated', self.on_treeview_Objects_row_activated )
            self.treeview.connect('button-release-event', self.on_treeview_Objects_button_release_event )
            '''
            '''--------------------------------------------------------------------------------------------'''


            #self.refresh_system_liststore()
            self.treeview_menu         = TreeViewMenu(self)
            self.window.show_all()                                               
            #self.system_names_combo.set_active(0)
            self.visible    =  True
            '''--------------------------------------------------------------------------------------------'''

    def CloseWindow (self, button, data  = None):
        """ Function doc """

        self.stop(None)
        #self.running = False  


        self.window.destroy()
        self.visible    =  False
        #print('self.visible',self.visible)
    
    def on_cell_active_radio_toggled (self, widget, path):
        """ Function doc """
        selected_path = Gtk.TreePath(path)
        for row in self.liststore:
            ##print(row.path, selected_path)
            row[0] = row.path == selected_path
            
            if row.path == selected_path:
                #print (row[0],  row[1],  row[2], self.modes[int(row[1])][1])
                
                self.lower = int(self.modes[int(row[1])][1][0])
                self.upper = int(self.modes[int(row[1])][1][1])
                self.value = int(self.lower) 
            else:
                pass    
    
    def on_button_import (self, button):
        """ Function doc """
        dialog = Gtk.FileChooserDialog("Select Folders", None,
                                       Gtk.FileChooserAction.SELECT_FOLDER,
                                       (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                                        Gtk.STOCK_OPEN, Gtk.ResponseType.OK))
        dialog.set_select_multiple(True)
        
        parameters = {}
        
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            # Get the list of selected folder paths
            folder_paths = dialog.get_filenames()
            #print(folder_paths)
            
            #-----------------------------------------------------------------
            tree_iter = self.system_names_combo.get_active_iter()
            if tree_iter != None:
                model = self.system_names_combo.get_model()
                parameters['system_id'] = model[tree_iter][1]
            
            else:
                parameters['system_id'] = 0
            #-----------------------------------------------------------------
            
            
            
            #-----------------------------------------------------------------
            tree_iter = self.coordinates_combobox.get_active_iter()
            if tree_iter != None:
                model = self.coordinates_combobox.get_model()
                parameters['vobject_id'] = model[tree_iter][1]
            
            else:
                parameters['vobject_id'] = -1
            #-----------------------------------------------------------------
            
            
            
            parameters['new_vobj_name'] = 'normal modes'
            parameters['data_type']     = 'normal_modes'

            parameters['trajectories'] = folder_paths
            parameters['logfile'] = False
            self.modes = self.p_session.import_data (parameters)
            
            
            self.liststore.clear()# = Gtk.ListStore(bool , #0 system_e_id           
                                  #         str  , #1 vobject 
                                  #         str  , #2 name 
                                  #
                                  #         )
            
            for i, data in self.modes.items():
                self.liststore.append([False, str(i), '{:.3f}'.format(float(data[0]))])
            
            
            
            
        dialog.destroy()
    
    
    def _coordinates_model_update (self, e_id):
        """ Function doc """
        #------------------------------------------------------------------------------------
        '''The combobox accesses, according to the id of the active system, 
        listostore of the dictionary object_list state_dict'''
        if self.Visible:
            #e_id = self.main.p_session.active_id 
            self.combobox_starting_coordinates.set_model(e_id)
            #------------------------------------------------------------------------------------
            size = len(self.main.vobject_liststore_dict[e_id])
            self.combobox_starting_coordinates.set_active(size-1)
            #------------------------------------------------------------------------------------
        else:
            pass
    
    
    def update_window (self, system_names = True, coordinates = False,  selections = True):
        """ Function doc """
        if self.visible:
            
            _id = self.system_names_combo.get_active()
            if _id == -1:
                '''_id = -1 means no item inside the combobox'''
                #self.selection_liststore.clear()
                #self.coordinates_liststore.clear()
                return None
            else:    
                _, system_id = self.main.system_liststore[_id]
            
            
            #if system_names:
            #    self.refresh_system_liststore ()
            #    self.system_names_combo.set_active(_id)
            #
            #if coordinates:
            #    self.refresh_coordinates_liststore ()
            #
            #
            if selections:
                _, system_id = self.main.system_liststore[_id]
                self.refresh_selection_liststore(system_id)
        else:
            pass
    
    
    def update (self, system_names = True, coordinates = False,  selections = True ):
        """ Function doc """
        pass
    
    def refresh_system_liststore (self):
        """ Function doc """
        self.main.refresh_system_liststore()
    
    def refresh_selection_liststore (self, system_id = None ):
        """ Function doc """
        self.selection_liststore.clear()
        
        
        ''' QC atoms'''
        if self.p_session.psystem[system_id].qcModel:
            self.p_session.psystem[system_id].e_selections["QC atoms"] = list(self.p_session.psystem[system_id].qcState.pureQCAtoms)


        '''Fixed atoms'''
        if self.p_session.psystem[system_id].freeAtoms is None:
            pass
        
        else:
            self.p_session.psystem[system_id].e_selections["fixed atoms"] = self.p_session.get_fixed_atoms_from_system(self.p_session.psystem[system_id])
  
        
        for selection , indexes in self.p_session.psystem[system_id].e_selections.items():
            self.selection_liststore.append([selection, len(indexes)])
    
    
    def refresh_coordinates_liststore(self, system_id = None):
        """ Function doc """
        cb_id = self.coordinates_combobox.get_active()
        
        if system_id:
            pass
        else:
            _id = self.system_names_combo.get_active()
            if _id == -1:
                return False
            else:
                #print('_id', _id)
                _, system_id = self.main.system_liststore[_id]
        
        self.coordinates_liststore.clear()
        for key , vobject in self.vm_session.vm_objects_dic.items():
            if vobject.e_id == system_id:
                self.coordinates_liststore.append([vobject.name, key])
        
        self.coordinates_combobox.set_active(len(self.coordinates_liststore)-1)
        
    
    def on_system_names_combobox_changed (self, widget):
        """ Function doc """
        index = self.system_names_combo.get_system_id()
       
        #if index == -1:
        #    '''_id = -1 means no item inside the combobox'''
        #    return None
        #
        #else:    
        #    _, system_id = self.main.system_liststore[index]
        
    
    def on_coordinates_combobox_changed(self, widget):
        """ Function doc """
        index = self.coordinates_combobox.get_active()
       
        if index == -1:
            '''_id = -1 means no item inside the combobox'''
            return None
        
        else:    
            name, vobject_id, system_id = self.coordinates_liststore[index]
            
            #print(name, vobject_id, system_id)
            
            try:
                vobject = self.vm_session.vm_objects_dic[vobject_id]
                self.liststore.clear()
                for i, data in vobject.normal_modes_dict.items():
                    self.liststore.append([False, str(i), data[0]])
            except:
                print('vobject has no Normal Modes data')
                pass
            
            #self.coordinates_liststore = Gtk.ListStore(str, int, int) 
            #self.coordinates_liststore.append(['New Visual Object', -1 , -1])
            #
            ##system  = self.main.p_session.get_system()
            #for row in self.main.vobject_liststore_dict[system_id]:
            #    
            #    a = list (row)
            #    self.coordinates_liststore.append([a[0], a[1], a[2]])
            #
            #self.coordinates_combobox.set_model(self.coordinates_liststore)
            #print (_, system_id)
            #
            ##self.refresh_selection_liststore (system_id)
            #
            ##size  =  len(list(self.main.vobject_liststore_dict[system_id]))
            #self.coordinates_combobox.set_active(0)

        
    def on_treeview_Objects_row_activated(self, tree, event, data):
        
        _id = self.system_names_combo.get_active()
        _, system_id = self.main.system_liststore[_id]
        
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        
        '''
        "key" is the acesses key to the dictionary containg the selection lists
        there is no two selection lists with the same name.
        indexes =  A list of atoms for selection
        '''
        key     = model.get_value(iter, 0)
        indexes = self.p_session.psystem[system_id].e_selections[key] 
    
        
        '''
        first we have to get the indexes from the comboboxes. These indexes are 
        not the same access keys as systema or vobject, they are just the 
        positions in the respective liststores
        '''
        sys_names_cb_index = self.system_names_combo.get_active()
        coords_cb_index    = self.coordinates_combobox.get_active()
        
        
        _, e_id          = self.main.system_liststore[sys_names_cb_index]
        _, vobj_id, e_id = self.main.vobject_liststore_dict[e_id][coords_cb_index]
        

        '''----------------------------- Applying the selection ------------------------------------'''
        vobject = self.vm_session.vm_objects_dic[vobj_id]
        self.vm_session.selections[self.vm_session.current_selection].selecting_by_indexes (vismol_object = vobject, indexes = indexes, clear = True)
        self.vm_session.selections[self.vm_session.current_selection].active = True
        '''-----------------------------------------------------------------------------------------'''
                
        self.vm_session.vm_glcore.queue_draw()
        
    
    def on_treeview_Objects_button_release_event(self, tree, event):
        '''
         str  ,   #                                   # 0
         bool ,   # toggle active=1                   # 1
         bool ,   # toggle visible = 3                # 2 
                                                      
         bool ,   # radio  active  = 2                # 3 
         bool ,   # radio  visible = 4                # 4 
                                                      
         bool  ,  # traj radio  active = 5            # 5 
         bool  ,  # is trajectory radio visible?      # 6 
                                                      
         int,     #                                   # 7
         int,     # pdynamo system index              # 8
         int,)    # frames  # 9
        '''
        
        
        _id = self.system_names_combo.get_active()
        if _id == -1:
            '''_id = -1 means no item inside the combobox'''
            return None
        else:    
            _, system_id = self.main.system_liststore[_id]
            
            
            
        if event.button == 3:
            selection     = self.treeview.get_selection()
            (model, iter) = selection.get_selected()
            for item in model:
                pass
                #print (item[0], model[iter][0])
            if iter != None:
                self.treeview_menu.open_menu(iter, system_id)

        if event.button == 1:
            print ('event.button == 1')


    def play (self, button):
        # Create and start a new thread when the button is clicked
        if self.running:
            pass
        
        else:
            thread = threading.Thread(target=self.long_task)
            thread.start()
            
            self.running = True
        
        
    def long_task(self):
        # This flag is used to stop the thread
        self.stop_thread = False
        i = 0
        
        self.lower
        self.upper
        
        self.value = self.lower
        while self.stop_thread == False:
            if self.stop_thread:
                return

            self.vm_session.set_frame(int(self.value))
            self.value += 1
            #print(self.value)

            if self.value >= self.upper:
                self.value = self.lower
                self.vm_session.set_frame(self.value)
            #time.sleep(0.01)
            time.sleep(1/self.spin_button.get_value())
            
    def stop (self, button):
        """ Function doc """
        self.stop_thread = True
        self.running = False    

    def forward (self, button):
        """ Function doc """
        value =  int(self.scale.get_value())
        value = value+1
        self.scale.set_value(int(value))
        self.vm_session.set_frame(int(value))
        #print(value)
        return value
        
        


class TreeViewMenu:
    """ Class doc """
    
    def __init__ (self, sele_window):
        """ Class initialiser """
        pass
        self.treeview = sele_window.treeview
        self.p_session = sele_window.p_session
        self.sele_window = sele_window 
        functions = {
                    'Rename'                : self.print_test ,
                    'Delete'                : self.delete_system ,
                    }
        self.build_tree_view_menu(functions)
        self.rename_window_visible = False

    def print_test (self, menu_item = None ):
        """  
        menu_item = Gtk.MenuItem object at 0x7fbdcc035700 (GtkMenuItem at 0x37cf6c0)
        
        """
        if self.rename_window_visible:
            pass
        else:
            #
            self.e_id     = self.sele_window.system_names_combo.get_active()
            selection     = self.treeview.get_selection()
            (model, iter) = selection.get_selected()
            self.key      = model.get_value(iter, 0)
            sys           = model.get_value(iter, 1)
            
            print('key: ', self.key, 'e_id: ',self.e_id)
            self.window = Gtk.Window()
            self.window.connect('destroy', self.destroy)
            self.window.set_keep_above(True)
            self.entry  = Gtk.Entry()
            
            self.entry.connect('activate', self.rename)
            self.window.add(self.entry)
            self.rename_window_visible = True
            self.window.show_all()
            #print(menu_item)

    def rename (self, menu_item):
        """ Function doc """
        print('New name: ', self.entry.get_text())
        new_name = self.entry.get_text()
        pass
        
        
        self.p_session.psystem[self.e_id].e_selections[new_name] = self.p_session.psystem[self.e_id].e_selections[self.key]
        self.p_session.psystem[self.e_id].e_selections.pop(self.key)
        self.sele_window.update_window()
        self.window.destroy()
        self.rename_window_visible = False
        
    def destroy (self, widget):
        """ Function doc """
        self.rename_window_visible = False

    def delete_system (self, menu_item = None ):
        """ Function doc """
        selection = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        #print(model[iter][0])


        sele = self.p_session.psystem[self.system_id].e_selections.pop(model[iter][0])
        #print ('deleting',sele)
        #print ('selections', self.p_session.systems[self.system_id]['selections'])
        self.sele_window.update_window (system_names = False, coordinates = False,  selections = True )


    def build_tree_view_menu (self, menu_items = None):
        """ Function doc """
        self.tree_view_menu = Gtk.Menu()
        for label in menu_items:
            mitem = Gtk.MenuItem(label)
            mitem.connect('activate', menu_items[label])
            self.tree_view_menu.append(mitem)
            #mitem = Gtk.SeparatorMenuItem()
            #self.tree_view_menu.append(mitem)

        self.tree_view_menu.show_all()

    
    def open_menu (self, vobject = None, system_id = None):
        """ Function doc """
        self.system_id = system_id
        #print (vobject)
        self.tree_view_menu.popup(None, None, None, None, 0, 0)


