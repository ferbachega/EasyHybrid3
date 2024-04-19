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
import numpy as np

from copy import deepcopy
import multiprocessing



from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from vismol.libgl.representations import SurfaceRepresentation

#----------------------------------------------------------------------
from pSimulation import QCGridPropertyGenerator
from pCore       import *
#----------------------------------------------------------------------

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')




class SurfaceAnalysisWindow(Gtk.Window):
    """ Class doc """
    def __init__(self, main = None, system_liststore = None ):
        """ Class initialiser """
        self.main       = main#self.main.system_liststore
        self.home       = main.home
        self.visible    = False        
        self.p_session  = main.p_session
        self.vm_session = main.vm_session
        
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
            self.builder.add_from_file(os.path.join(self.home,'src/gui/windows/analysis/surface_analysis_window.glade')) #/home/fernando/programs/EasyHybrid3/gui/windows/normal_mode_analysis.glade
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('window')
            self.window.set_default_size(200, 600)  
            self.window.set_title('Surfaces')  
            self.window.set_keep_above(True)
            
            
            
            
            # - - - - - - - systems combobox - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            self.box1 = self.builder.get_object('box_system')
            self.system_names_combo = SystemComboBox (self.main)
            self.system_names_combo.connect("changed", self.on_system_names_combobox_changed)
            '''--------------------------------------------------------------------------------------------'''
            self.box1.pack_start(self.system_names_combo, False, False, 0)

            self.coordinates_liststore = Gtk.ListStore(str, int, int)
            self.box2 = self.builder.get_object('box_coordinates')
            self.coordinates_combobox = CoordinatesComboBox(self.main.vobject_liststore_dict[self.p_session.active_id]) 
            
            self.box2.pack_start(self.coordinates_combobox, False, False, 0)
            
            
      
            system  = self.main.p_session.get_system()
            for row in self.main.vobject_liststore_dict[self.main.p_session.active_id]:
                
                a = list (row)
                self.coordinates_liststore.append([a[0], a[1], a[2]])

            
            self.coordinates_combobox.set_model(self.coordinates_liststore)
            



            '''
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

            #'''




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


            columns = [' ', 'Orbital', 'Occ.', 'Energy', 'visible']
            
            self.liststore = Gtk.ListStore(int, str, int, float, bool)
            #for item in data:
            #    self.liststore.append(item)

            self.treeview = self.builder.get_object('selection_treeview')#Gtk.TreeView(model=self.liststore)
            self.treeview.set_model(self.liststore)
             
            for i, column_title in enumerate(columns):
                renderer = Gtk.CellRendererText()
                if column_title == 'visible':
                    renderer_toggle = Gtk.CellRendererToggle()
                    renderer_toggle.connect("toggled", self.on_cell_toggled)
                    column = Gtk.TreeViewColumn(column_title, renderer_toggle, active=4)
                else:    
                    column = Gtk.TreeViewColumn(column_title, renderer, text=i)
                self.treeview.append_column(column)


            #self.refresh_system_liststore()
            #self.treeview_menu         = TreeViewMenu(self)
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

    def on_cell_toggled(self, widget, path):
        self.liststore[path][4] = not self.liststore[path][4]

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
        system_id = self.system_names_combo.get_system_id()
        system  = self.main.p_session.get_system(system_id)
        
        if system_id is not None:
            self.coordinates_combobox.set_model(self.main.vobject_liststore_dict[system_id])
            #self.refresh_selection_liststore (system_id)            
            size  =  len(list(self.main.vobject_liststore_dict[system_id]))
            self.coordinates_combobox.set_active(size-1)
        
        
        system.Energy  ( )
        orbitalsP = system.scratch.orbitalsP
        energies  = orbitalsP.energies
        
        LUMO      = orbitalsP.occupancyHandler.numberOccupied
        HOMO      = LUMO - 1
        
        
        
        orbitals = []
        for i,energy in enumerate(energies):
            if i >= LUMO:
                if i-LUMO == 0:
                    label = 'LUMO ' 
                else:
                    label = 'LUMO +'+str(i-LUMO)
            else:
                
                if i-HOMO == 0:
                    label = 'HOMO' 
                else:
                    label = 'HOMO '+str(i-HOMO)
            orbitals.append([i, label, orbitalsP.occupancies[i], energy, False ])
        
        for i in range(len(orbitals)):
            reverse_index = -i-1 #- len(orbitals)
            print(reverse_index, orbitals[reverse_index ])
            self.liststore.append(orbitals[reverse_index ])
        
        print(index)
    
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

    def on_treeview_Objects_row_activated(self, tree, event, data):

        vobject_id    = self.coordinates_combobox.get_vobject_id()
        vismol_object = self.main.vm_session.vm_objects_dic[vobject_id]

        #_id = self.system_names_combo.get_active()
        #_, system_id = self.main.system_liststore[_id]
        _isovalue    = float(self.builder.get_object('entry_isovalue').get_text())
        _GridSpacing = float(self.builder.get_object('entry_spacing') .get_text())
        
        system  = self.main.p_session.get_system()
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        #'''
        backup = []
        try:
            backup.append(system.e_treeview_iter)
            backup.append(system.e_liststore_iter)
            system.e_treeview_iter   = None
            system.e_liststore_iter  = None
        except:
            pass
        '''
        "key" is the acesses key to the dictionary containg the selection lists
        there is no two selection lists with the same name.
        indexes =  A list of atoms for selection
        '''
        key     = model.get_value(iter, 0)
        #indexes = self.p_session.psystem[system_id].e_selections[key] 
        #vismol_object = self.vm_session.vm_objects_dic[0]
        print(key, vismol_object.frames.shape[0])
        
        
        #_GridSpacing = 0.6
        _OrbitalTag    = "Grid Orbitals"
        _IsosurfaceTag = "Isosurface"
        
        trajectory = [None]*vismol_object.frames.shape[0]
        joblist = []
        for frame in range(vismol_object.frames.shape[0]):
            #self.p_session.get_coordinates_from_vobject_to_pDynamo_system(vobject)
            #'''
            #print('AQUIIIII', frame)
            self.p_session.get_coordinates_from_vobject_to_pDynamo_system( vobject = vismol_object, 
                                                                          system_id = None, 
                                                                          frame = frame)
            
            parameters = {
            '_GridSpacing'   : _GridSpacing,
            '_OrbitalTag'    : _OrbitalTag,
            '_isovalue'      : _isovalue,
            '_IsosurfaceTag' : _IsosurfaceTag,
            'orbital_key'    : key,
            }

            coords = self.p_session.get_coordinates_from_vobject (vobject = vismol_object, frame = frame)
            
            joblist.append([frame, system, coords, parameters])
            

        

        
        
        p = multiprocessing.Pool(processes = multiprocessing.cpu_count())
        results = p.map(generate_grid_parallel, joblist)
        
                #if interface:
        try:
            system.e_treeview_iter   = backup[0]
            system.e_liststore_iter  = backup[1]
        except:
            pass
        vismol_object.surface_trajectory = results # trajectory
        #-----------------------------------------------------------------------
        #generator.ExportProperty ( "/home/fernando/programs/EasyHybrid3/examples/scripts/tmp", _IsosurfaceTag )
        vismol_object.representations["surface1"] =  SurfaceRepresentation(vismol_object = vismol_object             ,
                                                                           vismol_glcore = self.vm_session.vm_glcore ,  
                                                                           name          = 'surface'                 ,
                                                                           active        = True                      ,
                                                                           indexes       = []                        ,
                                                                           is_dynamic    = False                     ,
                                                                           iso_color     = [1,0,0]                   ,
                                                                           surface_name  = 'obital_plus'                )
                                                     

        #-----------------------------------------------------------------------
        #generator.ExportProperty ( "/home/fernando/programs/EasyHybrid3/examples/scripts/tmp", _IsosurfaceTag )
        vismol_object.representations["surface2"] =  SurfaceRepresentation(vismol_object = vismol_object             ,
                                                                           vismol_glcore = self.vm_session.vm_glcore ,  
                                                                           name          = 'surface'                 ,
                                                                           active        = True                      ,
                                                                           indexes       = []                        ,
                                                                           is_dynamic    = False                     ,
                                                                           iso_color     = [0,0,1]                   ,
                                                                           surface_name  = 'obital_minus'           )
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



def surfece_parser ( surface, iso_color):
    """ Function doc """
    normals   = surface.polygonNormals
    polygons  = surface.polygons
    vertices  = surface.vertices
    colors    = []
    vertices2 = []
    normals2  = []
    
    for p in range ( polygons.rows ):
        #print ( "facet normal " )
        for c in range ( 3 ): 
            normals2.append (normals[p,c])
        #print ( "\n    outer loop" )
        for v in polygons[p,:]:
            #text = "\n        vertex "
            for c in range ( 3 ): 
                vertices2.append((vertices[v,c])/1.889725989 ) # convert from Bohr to angstrom
                #text += " {} ".format(vertices[v,c])
            for rgb in iso_color: 
                #vertices2.append(rgb)
                colors.append(rgb)
    
    vertices = np.array(vertices2, dtype=np.float32)
    colors   = np.array(colors, dtype=np.float32)
    indexes  = np.array(range(len(vertices)*3), dtype=np.uint32)
    return vertices, colors, indexes
    #active   = True
def apply_coords_to_system (system, coords):
    """ Function doc """
    
    for i, xyz in enumerate(coords):
        system.coordinates3[i][0] = xyz[0]
        system.coordinates3[i][1] = xyz[1]
        system.coordinates3[i][2] = xyz[2]

def generate_grid_parallel (job):
    """ Function doc 
    
    [frame, system, coords, parameters]
    
    """
    i          = job[0]
    system     = job[1]
    coords     = job[2]
    parameters = job[3]
    
    '''
    parameters = {
            '_GridSpacing'   : _GridSpacing
            '_OrbitalTag'    : _OrbitalTag
            '_IsosurfaceTag' : _IsosurfaceTag
            }
    '''
    
    _GridSpacing   = parameters['_GridSpacing']
    _OrbitalTag    = parameters['_OrbitalTag']
    _isovalue      = parameters['_isovalue']  
    _IsosurfaceTag = parameters['_IsosurfaceTag']
    key            = parameters['orbital_key']
    
    apply_coords_to_system(system, coords)
    system.Energy()
    
    
    #-----------------------------------------------------------------------
    # . Calculate the system grid properties.
    #-----------------------------------------------------------------------
    generator = QCGridPropertyGenerator.FromSystem (system )
    generator.DefineGrid    ( gridSpacing = _GridSpacing ) # . Some value in atomic units - e.g. 0.2
    

    
    orbital_iso = {}

    generator.GridOrbitals  ( [ key ]    ,       tag = _OrbitalTag    ) # . List of orbital indices (can be one only)    
    
    generator.Isosurface    ( _OrbitalTag, _isovalue, tag = _IsosurfaceTag )
    surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
    isosurface_p = surfaceProperty.isosurface # . This is the surface you can display.
    
    
    
    vertices, colors, indexes = surfece_parser ( surface = isosurface_p , iso_color = [1,0,0] )
    
    orbital_iso['obital_plus'] = [vertices, colors, indexes]
    
    generator.Isosurface    ( _OrbitalTag, _isovalue*-1, tag = _IsosurfaceTag )
    surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
    isosurface_n = surfaceProperty.isosurface # . This is the surface you can display.
    
    vertices, colors, indexes = surfece_parser ( surface = isosurface_n , iso_color = [0,0,1] )
    orbital_iso['obital_minus'] = [vertices, colors, indexes]
    
    return orbital_iso
    #trajectory[frame] = orbital_iso










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


