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

from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')


class SelectionListWindow(Gtk.Window):
    """ Class doc """
    def __init__(self, main = None, system_liststore = None ):
        """ Class initialiser """
        self.main_session          = main#self.main_session.system_liststore
        self.home                  = main.home
        self.visible               = False        
        self.p_session             = main.p_session
        self.vm_session            = main.vm_session
        

        
        self.selection_liststore       = Gtk.ListStore(str, int)
        self.selection_liststore_dict  =   {
                                           # system_e_id : Gtk.ListStore(str, int)
                                           }
    
        #-------------------------------------------
        # - - - - - - - - Restraints - - - - - - - - 
        #-------------------------------------------
        self.restraint_types = {0 : 'distance',
                                #1 : 'angle',
                                #2 : 'dihedral',
                                3 : 'all'
                                }
        self.rest_type_liststore = Gtk.ListStore(str)
        
        self.restraint_liststore       = Gtk.ListStore(bool, # True  - / False toggle
                                                        str, # name
                                                        str, # type  -  distance / angle / dihedral
                                                        str, # atoms -
                                                        str, # values
                                                        str, # force constant
                                                        int,) # e_id) 
        #--------------------------------------------
        
        
    def OpenWindow (self):
        """ Function doc /home/fernando/programs/VisMol/easyhybrid/gui/selection_list.glade"""
        if self.visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'src/gui/windows/setup/selection_and_restraint_list.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('window')
            self.window.set_default_size(200, 600)  
            self.window.set_title('Selections')  
            self.window.set_keep_above(True)
            
            # - - - - - - - systems combobox - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            self.box = self.builder.get_object('box_system')
            self.combobox_systems = SystemComboBox(self.main_session)
            self.combobox_systems.connect("changed", self.on_combobox_systemsbox_changed)
            '''--------------------------------------------------------------------------------------------'''
            self.box.pack_start(self.combobox_systems, False, False, 0)



            # - - - - - - - coordinates combobox - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            self.box_coordinates = self.builder.get_object('box_coordinates')
            self.coordinates_combobox = CoordinatesComboBox() #self.builder.get_object('coordinates_combobox')
            self.box_coordinates.pack_start(self.coordinates_combobox, False, False, 0)
            '''--------------------------------------------------------------------------------------------'''


            system  = self.main_session.p_session.get_system()
            self.combobox_systems.set_active_iter(system.e_liststore_iter)
            


            # - - - - - - - Selection Treeview - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            self.treeview = self.builder.get_object('selection_treeview')
            for i, column_title in enumerate(['Selection',"Number of Atoms"]):
                renderer = Gtk.CellRendererText()
                column = Gtk.TreeViewColumn(column_title, renderer, text=i)
                self.treeview.append_column(column)
            self.treeview.set_model(self.selection_liststore)

            self.treeview.connect('row_activated', self.on_treeview_Objects_row_activated )
            self.treeview.connect('button-release-event', self.on_treeview_Objects_button_release_event )
            '''--------------------------------------------------------------------------------------------'''

            
            
            #---------------------------------------------------------------------------------
            #                - - - - - - - Restraint Treeview - - - - - - -
            #---------------------------------------------------------------------------------
            self.restraint_treeview = self.builder.get_object('restraint_treeview')
            
            renderer_toggle = Gtk.CellRendererToggle()
            renderer_toggle.connect("toggled", self._on_cell_visible_toggled)
            column_toggle = Gtk.TreeViewColumn("A", renderer_toggle, active=0)#, visible = True)
            self.restraint_treeview.append_column(column_toggle)
            
            for i, column_title in enumerate(["Index", 'Type', 'Atoms', "Dist/Angle", "Force Constant"]):
                renderer = Gtk.CellRendererText()
                column = Gtk.TreeViewColumn(column_title, renderer, text=i+1)
                self.restraint_treeview.append_column(column)
            self.restraint_treeview.set_model(self.restraint_liststore)
            #---------------------------------------------------------------------------------
            
            
            
            
            self.rest_type_liststore.clear()
            #---------------------------------------------------------------------------------
            #                           Restraints COMBOBOX
            #---------------------------------------------------------------------------------
            self.comobobox_restraints = self.builder.get_object('combobox_restraints')
            for key, rest_type in self.restraint_types.items():
                self.rest_type_liststore.append([rest_type])
            
            self.comobobox_restraints.set_model(self.rest_type_liststore)
            self.comobobox_restraints.connect("changed", self._on_rest_types_changed)
            renderer_text = Gtk.CellRendererText()
            self.comobobox_restraints.pack_start(renderer_text, True)
            self.comobobox_restraints.add_attribute(renderer_text, "text", 0)
            self.comobobox_restraints.set_active(0)
            #self.restraint_treeview.connect('row_activated', self.on_treeview_Objects_row_activated )
            #self.restraint_treeview.connect('button-release-event', self.on_treeview_Objects_button_release_event )
            '''--------------------------------------------------------------------------------------------'''



            #self.refresh_system_liststore()
            self.treeview_menu         = TreeViewMenu(self)
            self.window.show_all()                                               
            #self.combobox_systems.set_active(0)
            self.visible    =  True
            '''--------------------------------------------------------------------------------------------'''
            self.update_window()
        else:
            self.window.present()
            
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        #self.BackUpWindowData()
        self.window.destroy()
        self.visible    =  False
        #print('self.visible',self.visible)
    
    def _on_rest_types_changed (self, widget):
        """ Function doc """
        print(self.comobobox_restraints.get_active())
    
    
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
    
     
    #def update_restraint_representation (_):
    #    """ Function doc """




     
        
    def _on_cell_visible_toggled (self, widget, path):
        """ Function doc """
        #print(widget)
        ##print(list(path))
        self.restraint_liststore[path][0] = not self.restraint_liststore[path][0]
        
        print(list(self.restraint_liststore[path]))
        e_id = self.restraint_liststore[path][6]
        name = self.restraint_liststore[path][1]
        system = self.main_session.p_session.psystem[e_id]
        system.e_restraints_dict[name][0] =  self.restraint_liststore[path][0]
        
        self.p_session.update_restaint_representation(e_id)
        

            

        
    def update_window (self, system_names = True, coordinates = False,  selections = True, restraints = True):
        """ Function doc """
        if self.visible:
            
            system_id = self.combobox_systems.get_system_id()
            if system_id is not None:
                if selections:
                    self.refresh_selection_liststore(system_id)
                if restraints:
                    self.refresh_restraint_liststore(system_id)
        else:
            self.refresh_restraint_liststore(self.main_session.p_session.active_id)
    
    def update (self, system_names = True, coordinates = False,  selections = True ):
        """ Function doc """
        pass
        #if self.visible:
        #    
        #    _id = self.combobox_systems.get_active()
        #    if _id == -1:
        #        '''_id = -1 means no item inside the combobox'''
        #        #self.selection_liststore.clear()
        #        #self.coordinates_liststore.clear()
        #        return None
        #    else:    
        #        _, system_id = self.main_session.system_liststore[_id]
        #    
        #    if system_names:
        #        self.refresh_system_liststore ()
        #        self.combobox_systems.set_active(_id)
        #    
        #    if coordinates:
        #        self.refresh_coordinates_liststore ()
        #    
        #    
        #    if selections:
        #        _, system_id = self.main_session.system_liststore[_id]
        #        self.refresh_selection_liststore(system_id)
        #else:
        #    pass

    
    def refresh_system_liststore (self):
        """ Function doc """
        self.main_session.refresh_system_liststore()
    
    
    def refresh_restraint_liststore (self, system_id = None):
        """ Function doc """
        self.restraint_liststore.clear()
        

        for name, restraint in self.p_session.psystem[system_id].e_restraints_dict.items():
            _bool = restraint[0]
            name  = restraint[1]
            _type = restraint[2]
            if _type == 'distance':
                atons = '{} / {}'.format(restraint[3][0],restraint[3][1]) 
            dist_or_angle = '{:.4f}'.format(restraint[4])
            force_const   = str(restraint[5])
            e_id          =  restraint[6] 
                                           #(bool,  str,   str,   str,          str,       str    , int    )
            self.restraint_liststore.append([_bool, name, _type, atons, dist_or_angle, force_const, e_id   ])
            
        self.p_session.update_restaint_representation(system_id)
        self.main_session.vm_session.vm_glcore.queue_draw()

                
    
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
            _id = self.combobox_systems.get_active()
            if _id == -1:
                return False
            else:
                #print('_id', _id)
                _, system_id = self.main_session.system_liststore[_id]
        
        self.coordinates_liststore.clear()
        for key , vobject in self.vm_session.vm_objects_dic.items():
            if vobject.e_id == system_id:
                self.coordinates_liststore.append([vobject.name, key])
        
        self.coordinates_combobox.set_active(len(self.coordinates_liststore)-1)
        
        
    def on_combobox_systemsbox_changed(self, widget):
        """ Function doc """
        system_id = self.combobox_systems.get_system_id()
       
        if system_id is not None:
            self.coordinates_combobox.set_model(self.main_session.vobject_liststore_dict[system_id])
            self.refresh_selection_liststore (system_id)            
            size  =  len(list(self.main_session.vobject_liststore_dict[system_id]))
            self.coordinates_combobox.set_active(size-1)
            
            self.update_window ( selections = False, restraints = True)
        
    def on_treeview_Objects_row_activated(self, tree, event, data):
        
        system_id = self.combobox_systems.get_system_id()
        
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        
        '''
        "key" is the acesses key to the dictionary containg the selection lists
        there is no two selection lists with the same name.
        indexes =  A list of atoms for selection
        '''
        key     = model.get_value(iter, 0)
        indexes = self.p_session.psystem[system_id].e_selections[key] 
        #print(self.p_session.psystem[system_id].e_qc_table)
        #print(self.p_session.psystem[system_id].e_selections)
        
        pymol_indexes = []
        for i in indexes:
            pymol_indexes.append(i+1)
        print(pymol_indexes)
        '''
        first we have to get the indexes from the comboboxes. These indexes are 
        not the same access keys as system or vobject, they are just the 
        positions in the respective liststores
        '''
        sys_names_cb_index = self.combobox_systems.get_active()
        coords_cb_index    = self.coordinates_combobox.get_active()
        
        
        #_, e_id          = self.main_session.system_liststore[sys_names_cb_index]
        _, vobj_id, e_id, pixbuf = self.main_session.vobject_liststore_dict[system_id][coords_cb_index]
        

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
        
        
        system_id = self.combobox_systems.get_system_id()
            
            
        if event.button == 3:
            selection     = self.treeview.get_selection()
            (model, iter) = selection.get_selected()
            for item in model:
                pass
                #print (item[0], model[iter][0])
            if iter != None:
                self.treeview_menu.open_menu(iter, system_id)

        if event.button == 1:
            print ('event.button == 1:')



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
            self.e_id     = self.sele_window.combobox_systems.get_active()
            selection     = self.treeview.get_selection()
            (model, iter) = selection.get_selected()
            self.key      = model.get_value(iter, 0)
            sys           = model.get_value(iter, 1)
            
            print(self.key, self.e_id)
            self.window = Gtk.Window()
            self.window.connect('destroy', self.destroy)
            self.window.set_keep_above(True)
            self.entry  = Gtk.Entry()
            
            self.entry.connect('activate', self.rename)
            self.window.add(self.entry)
            self.rename_window_visible = True
            self.window.show_all()
            print(menu_item)

    def rename (self, menu_item):
        """ Function doc """
        print(self.entry.get_text())
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


