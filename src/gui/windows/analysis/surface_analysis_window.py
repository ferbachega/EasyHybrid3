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
from vismol.core.vismol_object import VismolObject

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
        self.frame      = self.vm_session.frame
        
        self.orbital_liststore_dict = {
                                       #vm_object_id:[liststore1, liststore2...]
                                       }
        self.wave_function_dict = {
                                       #vm_object_id:[orbitals, gridgenerator]
                                      }
        self.counter = 1000
        
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
            self.box1.pack_start(self.system_names_combo, False, False, 0)
            '''--------------------------------------------------------------------------------------------'''

            
            
            '''--------------------------------------------------------------------------------------------'''
            self.coordinates_liststore = Gtk.ListStore(str, int, int)
            self.box2 = self.builder.get_object('box_coordinates')
            self.coordinates_combobox = CoordinatesComboBox(self.main.vobject_liststore_dict[self.p_session.active_id]) 
            self.box2.pack_start(self.coordinates_combobox, False, False, 0)
            '''--------------------------------------------------------------------------------------------'''


            self.btn_import_wfunction =  self.builder.get_object('btn_import_wavefunction')
            self.btn_import_wfunction.connect('clicked', self.on_button_import_wavefunction)


            #                       SURFACE TYPE COMBOBOX
            #'''--------------------------------------------------------------------------------------------'''
            self.box_surface_type = self.builder.get_object('box_surface_type')
            
            self.cbx_surface_type = Gtk.ComboBoxText()
            self.cbx_surface_type.connect("changed", self.surface_combobox_change)

            self.surface_options = ["Orbitals"               , 
                                    "Electrostatic Potential",
                                    "Density"              , 
                                    'External'               ]
            for i, item in  enumerate(self.surface_options):
                self.cbx_surface_type.insert(i, str(i), item )
                
            self.box_surface_type.pack_start(self.cbx_surface_type, False, False, 0)
            self.cbx_surface_type.set_active(0)
            #'''--------------------------------------------------------------------------------------------'''


            self.label_frame = self.builder.get_object('label_frame')
                 
            system  = self.main.p_session.get_system()
            for row in self.main.vobject_liststore_dict[self.main.p_session.active_id]:
                
                a = list (row)
                self.coordinates_liststore.append([a[0], a[1], a[2]])

            self.coordinates_combobox.set_model(self.coordinates_liststore)
            

            columns = [' ', 'Orbital', 'Occ.', 'Energy']#, 'visible']
            
            self.liststore = Gtk.ListStore(int, str, int, float, bool)

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
            
            #-----------------------------------------------------------
            self.btn_render = self.builder.get_object('btn_render')
            self.btn_render.connect('clicked', self.on_render_button)
            #-----------------------------------------------------------
            
            #self.refresh_system_liststore()
            #self.treeview_menu         = TreeViewMenu(self)
            self.window.show_all()                                               
            
            
            self.builder.get_object('btn_external_file').hide()
            self.builder.get_object('label_external_file').hide()
            #self.system_names_combo.set_active(0)
            self.visible    =  True
            '''--------------------------------------------------------------------------------------------'''

    def CloseWindow (self, button, data  = None):
        """ Function doc """

        #self.stop(None)
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
            if vobject.is_surface:
                print(vobject.is_surface)
                pass
            else:
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
        
        
    def surface_combobox_change (self, widget):
        """ Function doc """
        index = self.cbx_surface_type.get_active()
        print(index)
        if index in [1,2]:
            #self.builder.get_object('btn_import_wavefunction').set_sensitive(False)
            #self.builder.get_object('selection_treeview')     .set_sensitive(False)
            self.builder.get_object('btn_import_wavefunction').hide()
            self.builder.get_object('selection_treeview')     .hide()
            self.builder.get_object('label_external_file').hide()
            self.builder.get_object('btn_external_file')  .hide()
        
        elif index == 3:
            self.builder.get_object('label_external_file').show()
            self.builder.get_object('btn_external_file').show()
            self.builder.get_object('btn_import_wavefunction').set_sensitive(False)
            self.builder.get_object('selection_treeview'     ).set_sensitive(False)
            self.builder.get_object('btn_import_wavefunction').hide()
            self.builder.get_object('selection_treeview'     ).hide()
        
        else:
            self.builder.get_object('btn_import_wavefunction').set_sensitive(True)
            self.builder.get_object('selection_treeview')     .set_sensitive(True)
            self.builder.get_object('btn_import_wavefunction').show()
            self.builder.get_object('selection_treeview')     .show()
            self.builder.get_object('label_external_file').hide()
            self.builder.get_object('btn_external_file')  .hide()

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

        system_id = self.system_names_combo.get_system_id()
        system    = self.main.p_session.psystem[system_id]

        vobject_id    = self.coordinates_combobox.get_vobject_id()
        vismol_object = self.main.vm_session.vm_objects_dic[vobject_id]

        _isovalue    = float(self.builder.get_object('entry_isovalue').get_text())
        _GridSpacing = float(self.builder.get_object('entry_spacing') .get_text())
        
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        

        #system  = self.main.p_session.get_system()
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
        print(key, vismol_object.frames.shape[0], model, model[iter][1])
        name = str(key) +' '+model[iter][1]#+' '+ str(model[iter][3])
        #_GridSpacing = 0.6
        _OrbitalTag    = "Grid Orbitals"
        _IsosurfaceTag = "Isosurface"
        
        
        trajectory = [None]*vismol_object.frames.shape[0]
        joblist = []
        
        for frame in range(vismol_object.frames.shape[0]):
            #self.p_session.get_coordinates_from_vobject_to_pDynamo_system(vobject)
            #'''
            print(vismol_object, frame)
            self.p_session.get_coordinates_from_vobject_to_pDynamo_system( vobject = vismol_object, 
                                                                           system_id = None, 
                                                                           frame = frame)
            
            parameters = {
            'type'           : 'orbital',
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
        
        #vobject_tmp = VismolObject(name="UNK", index=-1,
        vobject_tmp = VismolObject(name= name , index=-1,
                                   vismol_session        = self.vm_session,
                                   trajectory            = [],
                                   bonds_pair_of_indexes = [0,1])
        
        
        vobject_tmp.model_mat = vismol_object.model_mat 
        vobject_tmp.trans_mat = vismol_object.trans_mat 
        
        #vismol_object.surface_trajectory = results # trajectory
        vobject_tmp.surface_trajectory = results # trajectory
        #-----------------------------------------------------------------------
        #generator.ExportProperty ( "/home/fernando/programs/EasyHybrid3/examples/scripts/tmp", _IsosurfaceTag )
        #vismol_object.representations["surface1"] =  SurfaceRepresentation(vismol_object = vismol_object             ,
        vobject_tmp.representations["surface1"] =  SurfaceRepresentation(vismol_object = vobject_tmp             ,
                                                                           vismol_glcore = self.vm_session.vm_glcore ,  
                                                                           name          = 'surface'                 ,
                                                                           active        = True                      ,
                                                                           indexes       = []                        ,
                                                                           is_dynamic    = False                     ,
                                                                           iso_color     = [1,0,0]                   ,
                                                                           surface_name  = 'obital_plus'                )
                                                     

        #-----------------------------------------------------------------------
        #generator.ExportProperty ( "/home/fernando/programs/EasyHybrid3/examples/scripts/tmp", _IsosurfaceTag )
        #vismol_object.representations["surface2"] =  SurfaceRepresentation(vismol_object = vismol_object             ,
        vobject_tmp.representations["surface2"] =  SurfaceRepresentation(vismol_object = vobject_tmp             ,
                                                                           vismol_glcore = self.vm_session.vm_glcore ,  
                                                                           name          = 'surface'                 ,
                                                                           active        = True                      ,
                                                                           indexes       = []                        ,
                                                                           is_dynamic    = False                     ,
                                                                           iso_color     = [0,0,1]                   ,
                                                                           surface_name  = 'obital_minus'           )
        
        vobject_tmp.parameters = parameters
        
        vobject_tmp.frames = vismol_object.frames
        vobject_tmp.active = True
        vobject_tmp.is_surface = True
        vobject_tmp.e_id = system.e_id
        self.vm_session._add_vismol_object(vobject_tmp, show_molecule=False, autocenter=False)
        
        
        self.main.main_treeview.add_vismol_object_to_treeview(vobject_tmp,vismol_object.e_treeview_iter )
        self.main.add_vobject_to_vobject_liststore_dict(vobject_tmp)
        self.main.refresh_widgets()
        self.vm_session.vm_glcore.queue_draw()
        self.counter +=1
        
        
        
        '''
        for frame , data in enumerate(self.wave_function_dict[vobject_id]):
            
            
            generator = data[2]
            #system.Energy()
            parameters = {
            '_GridSpacing'   : _GridSpacing,
            '_OrbitalTag'    : _OrbitalTag,
            '_isovalue'      : _isovalue,
            '_IsosurfaceTag' : _IsosurfaceTag,
            'orbital_key'    : key,
            }
            
            joblist.append([frame, generator, parameters])
        
        p = multiprocessing.Pool(processes = multiprocessing.cpu_count())
        #'''
        '''
        results = p.map(generate_grid_parallel, joblist)
        
        print (results)
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
        #'''
    
    
    
    def _update_liststore (self):
        """ Function doc """
        
        vobject_id = self.coordinates_combobox.get_vobject_id()
        print(vobject_id)
        
        model = self.treeview.get_model()
        if model is not None:
            # Remove todos os itens do modelo
            model.clear()
       
        if self.frame > len(self.orbital_liststore_dict[vobject_id]):
            self.treeview.set_model(self.orbital_liststore_dict[vobject_id][-1])
        else:
            
            #for frame , data in enumerate(self.wave_function_dict[vobject_id]):
            orbitals = self.wave_function_dict[vobject_id][self.frame][0]
            for i in range(len(orbitals)):
                reverse_index = -i-1 #- len(orbitals)
                model.append(orbitals[reverse_index ])
            
        
    def set_frame (self ):
        """ Function doc """
        self.frame =  self.vm_session.frame
        self.label_frame.set_text('Frame = {}'.format(self.frame))
        self._update_liststore()

    
    def on_render_button (self, widget):
        """ Function doc """
        index = self.cbx_surface_type.get_active()
        print(index)

        system_id = self.system_names_combo.get_system_id()
        system    = self.main.p_session.psystem[system_id]

        vobject_id    = self.coordinates_combobox.get_vobject_id()
        vismol_object = self.main.vm_session.vm_objects_dic[vobject_id]

        _isovalue    = float(self.builder.get_object('entry_isovalue').get_text())
        _GridSpacing = float(self.builder.get_object('entry_spacing') .get_text())
        
        
        backup = []
        try:
            backup.append(system.e_treeview_iter)
            backup.append(system.e_liststore_iter)
            system.e_treeview_iter   = None
            system.e_liststore_iter  = None
        except:
            pass
        
        
        if index == 2:
            joblist = []
            for frame in range(vismol_object.frames.shape[0]):
                #'''
                self.p_session.get_coordinates_from_vobject_to_pDynamo_system( vobject = vismol_object, 
                                                                               system_id = None, 
                                                                               frame = frame)
                parameters = {
                'type'           : 'density',
                '_GridSpacing'   : _GridSpacing,
                '_OrbitalTag'    : 'density',
                '_isovalue'      : _isovalue,
                '_IsosurfaceTag' : 'density',
                'orbital_key'    : 0,
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
            name ='Density'
            #vobject_tmp = VismolObject(name="UNK", index=-1,
            vobject_tmp = VismolObject(name= name , index=-1,
                                        vismol_session        = self.vm_session,
                                        trajectory            = [],
                                        bonds_pair_of_indexes = [0,1])
            
            
            vobject_tmp.model_mat = vismol_object.model_mat 
            vobject_tmp.trans_mat = vismol_object.trans_mat 
            
            #vismol_object.surface_trajectory = results # trajectory
            vobject_tmp.surface_trajectory = results # trajectory
            #-----------------------------------------------------------------------
            #generator.ExportProperty ( "/home/fernando/programs/EasyHybrid3/examples/scripts/tmp", _IsosurfaceTag )
            #vismol_object.representations["surface1"] =  SurfaceRepresentation(vismol_object = vismol_object             ,
            vobject_tmp.representations["surface1"] =  SurfaceRepresentation(vismol_object = vobject_tmp             ,
                                                                                vismol_glcore = self.vm_session.vm_glcore ,  
                                                                                name          = 'surface'                 ,
                                                                                active        = True                      ,
                                                                                indexes       = []                        ,
                                                                                is_dynamic    = False                     ,
                                                                                iso_color     = [0,0,1]                   ,
                                                                                surface_name  = 'obital_plus'                )
                                                            
        
            #-----------------------------------------------------------------------
            #generator.ExportProperty ( "/home/fernando/programs/EasyHybrid3/examples/scripts/tmp", _IsosurfaceTag )
            #vismol_object.representations["surface2"] =  SurfaceRepresentation(vismol_object = vismol_object             ,
            
            #vobject_tmp.representations["surface2"] =  SurfaceRepresentation(vismol_object = vobject_tmp             ,
            #                                                                    vismol_glcore = self.vm_session.vm_glcore ,  
            #                                                                    name          = 'surface'                 ,
            #                                                                    active        = True                      ,
            #                                                                    indexes       = []                        ,
            #                                                                    is_dynamic    = False                     ,
            #                                                                    iso_color     = [0,0,1]                   ,
            #                                                                    surface_name  = 'obital_minus'           )
            #
            
            vobject_tmp.parameters = parameters
            
            vobject_tmp.frames = vismol_object.frames
            vobject_tmp.active = True
            vobject_tmp.is_surface = True
            vobject_tmp.e_id = system.e_id
            self.vm_session._add_vismol_object(vobject_tmp, show_molecule=False, autocenter=False)
            
            
            self.main.main_treeview.add_vismol_object_to_treeview(vobject_tmp,vismol_object.e_treeview_iter )
            # Add the VisMol object to the vobject liststore dictionary
            self.main.add_vobject_to_vobject_liststore_dict(vobject_tmp)
            # Refresh the widgets in the main window
            self.main.refresh_widgets()
            
            
            
            
            #print(vobject_tmp, vobject_tmp.surface_trajectory)
            #self.vm_session.vm_objects_dic[self.counter] = vobject_tmp
            self.vm_session.vm_glcore.queue_draw()
            self.counter +=1

        elif index ==1:
            return False
            joblist = []
            for frame in range(vismol_object.frames.shape[0]):
                #'''
                self.p_session.get_coordinates_from_vobject_to_pDynamo_system( vobject = vismol_object, 
                                                                               system_id = None, 
                                                                               frame = frame)
                parameters = {
                'type'           : 'potential',
                '_GridSpacing'   : _GridSpacing,
                '_OrbitalTag'    : 'potential',
                '_isovalue'      : _isovalue,
                '_IsosurfaceTag' : 'potential',
                'orbital_key'    : 0,
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
            name ='Potential'
            vobject_tmp = VismolObject(name= name , index=-1,
                                        vismol_session        = self.vm_session,
                                        trajectory            = [],
                                        bonds_pair_of_indexes = [0,1])
            
            
            vobject_tmp.model_mat = vismol_object.model_mat 
            vobject_tmp.trans_mat = vismol_object.trans_mat 
            vobject_tmp.surface_trajectory = results # trajectory
            #-----------------------------------------------------------------------
            vobject_tmp.representations["surface1"] =  SurfaceRepresentation(vismol_object = vobject_tmp             ,
                                                                                vismol_glcore = self.vm_session.vm_glcore ,  
                                                                                name          = 'surface'                 ,
                                                                                active        = True                      ,
                                                                                indexes       = []                        ,
                                                                                is_dynamic    = False                     ,
                                                                                iso_color     = [1,0,0]                   ,
                                                                                surface_name  = 'obital_plus'             )
                                                            
        
            #-----------------------------------------------------------------------
            vobject_tmp.representations["surface2"] =  SurfaceRepresentation(vismol_object = vobject_tmp             ,
                                                                                vismol_glcore = self.vm_session.vm_glcore ,  
                                                                                name          = 'surface'                 ,
                                                                                active        = True                      ,
                                                                                indexes       = []                        ,
                                                                                is_dynamic    = False                     ,
                                                                                iso_color     = [0,0,1]                   ,
                                                                                surface_name  = 'obital_minus'           )
            
            
            vobject_tmp.parameters = parameters
            
            vobject_tmp.frames = vismol_object.frames
            vobject_tmp.active = True
            vobject_tmp.is_surface = True
            vobject_tmp.e_id = system.e_id
            self.vm_session._add_vismol_object(vobject_tmp, show_molecule=False, autocenter=False)
            
            
            self.main.main_treeview.add_vismol_object_to_treeview(vobject_tmp,vismol_object.e_treeview_iter )
            # Add the VisMol object to the vobject liststore dictionary
            self.main.add_vobject_to_vobject_liststore_dict(vobject_tmp)
            # Refresh the widgets in the main window
            self.main.refresh_widgets()
            self.vm_session.vm_glcore.queue_draw()
            self.counter +=1

        elif index == 0:
            selection     = self.treeview.get_selection()
            (model, iter) = selection.get_selected()
        
            
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
            print(key, vismol_object.frames.shape[0], model, model[iter][1])
            name = str(key) +' '+model[iter][1]#+' '+ str(model[iter][3])
            #_GridSpacing = 0.6
            _OrbitalTag    = "Grid Orbitals"
            _IsosurfaceTag = "Isosurface"
            
            
            trajectory = [None]*vismol_object.frames.shape[0]
            joblist = []
            
            for frame in range(vismol_object.frames.shape[0]):
                #self.p_session.get_coordinates_from_vobject_to_pDynamo_system(vobject)
                #'''
                self.p_session.get_coordinates_from_vobject_to_pDynamo_system( vobject = vismol_object, 
                                                                               system_id = None, 
                                                                               frame = frame)
                
                parameters = {
                'type'           : 'orbital',
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
            
            vobject_tmp = VismolObject(name= name , index=-1,
                                       vismol_session        = self.vm_session,
                                       trajectory            = [],
                                       bonds_pair_of_indexes = [0,1])
            
            vobject_tmp.model_mat = vismol_object.model_mat 
            vobject_tmp.trans_mat = vismol_object.trans_mat 
            vobject_tmp.surface_trajectory = results # trajectory
            #-----------------------------------------------------------------------
            vobject_tmp.representations["surface1"] =  SurfaceRepresentation(vismol_object = vobject_tmp             ,
                                                                               vismol_glcore = self.vm_session.vm_glcore ,  
                                                                               name          = 'surface'                 ,
                                                                               active        = True                      ,
                                                                               indexes       = []                        ,
                                                                               is_dynamic    = False                     ,
                                                                               iso_color     = [1,0,0]                   ,
                                                                               surface_name  = 'obital_plus'                )
                                                         

            #-----------------------------------------------------------------------
            vobject_tmp.representations["surface2"] =  SurfaceRepresentation(vismol_object = vobject_tmp             ,
                                                                               vismol_glcore = self.vm_session.vm_glcore ,  
                                                                               name          = 'surface'                 ,
                                                                               active        = True                      ,
                                                                               indexes       = []                        ,
                                                                               is_dynamic    = False                     ,
                                                                               iso_color     = [0,0,1]                   ,
                                                                               surface_name  = 'obital_minus'           )
            
            vobject_tmp.parameters = parameters
            
            vobject_tmp.frames = vismol_object.frames
            vobject_tmp.active = True
            vobject_tmp.is_surface = True
            vobject_tmp.e_id = system.e_id
            self.vm_session._add_vismol_object(vobject_tmp, show_molecule=False, autocenter=False)
            
            
            self.main.main_treeview.add_vismol_object_to_treeview(vobject_tmp,vismol_object.e_treeview_iter )
            # Add the VisMol object to the vobject liststore dictionary
            self.main.add_vobject_to_vobject_liststore_dict(vobject_tmp)
            # Refresh the widgets in the main window
            self.main.refresh_widgets()
            self.vm_session.vm_glcore.queue_draw()
            self.counter +=1













    def on_button_import_wavefunction (self, widget):
        """ Function doc """
        print('on_button_import_wavefunction')
        
        system_id = self.system_names_combo.get_system_id()
        system    = self.main.p_session.psystem[system_id]
        
        vobject_id    = self.coordinates_combobox.get_vobject_id()
        vismol_object = self.main.vm_session.vm_objects_dic[vobject_id]
        print(system_id, vobject_id, system, vismol_object)
        backup = []
        try:
            backup.append(system.e_treeview_iter)
            backup.append(system.e_liststore_iter)
            system.e_treeview_iter   = None
            system.e_liststore_iter  = None
        except:
            pass
        
        
        #frame = self.vm_session.frame
        #coords = self.p_session.get_coordinates_from_vobject (vobject = vismol_object, frame = frame)
        #joblist = [[frame, system, coords]] 
        
        #'''
        trajectory = [None]*vismol_object.frames.shape[0]
        joblist = []
        for frame in range(vismol_object.frames.shape[0]):
            self.p_session.get_coordinates_from_vobject_to_pDynamo_system( vobject = vismol_object, 
                                                                          system_id = system_id, 
                                                                          frame = frame)
            parameters = None
            coords = self.p_session.get_coordinates_from_vobject (vobject = vismol_object, frame = frame)
            joblist.append([frame, system, coords])
        #'''    
        p = multiprocessing.Pool(processes = multiprocessing.cpu_count())
        results = p.map(generate_wavefunction_parallel, joblist)
        
        self.wave_function_dict[vobject_id] = results
        

        try:
            system.e_treeview_iter   = backup[0]
            system.e_liststore_iter  = backup[1]
        except:
            pass

        #print(self.wave_function_dict)


        
        #'''
        
        self.orbital_liststore_dict[vobject_id]= []
        
        for frame , data in enumerate(self.wave_function_dict[vobject_id]):
            orbitals = data[0]
            
            self.liststore = Gtk.ListStore(int, str, int, float, bool)
            for i in range(len(orbitals)):
                reverse_index = -i-1 #- len(orbitals)
                #print(reverse_index, orbitals[reverse_index ])
        
                self.liststore.append(orbitals[reverse_index ])
                self.orbital_liststore_dict[vobject_id].append(self.liststore)
        
        print()
        self.treeview.set_model(self.orbital_liststore_dict[vobject_id][self.frame])


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



def surface_parser ( surface, iso_color):
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
    #_type      = job[4]
    '''
    parameters = {
            '_GridSpacing'   : _GridSpacing
            '_OrbitalTag'    : _OrbitalTag
            '_IsosurfaceTag' : _IsosurfaceTag
            }
    '''
    _type          = parameters['type']
    _GridSpacing   = parameters['_GridSpacing']
    _OrbitalTag    = parameters['_OrbitalTag']
    _isovalue      = parameters['_isovalue']  
    _IsosurfaceTag = parameters['_IsosurfaceTag']
    key            = parameters['orbital_key']
    #print(parameters, type(system))
    #-----------------------------------------------------------------------
    # . Calculate the system grid properties.
    #-----------------------------------------------------------------------
    #system    = system.Energy()
    #energies  = orbitalsP.energies
    
    apply_coords_to_system(system, coords)
    system.Energy()
    
    #-----------------------------------------------------------------------
    # . Calculate the system grid properties.
    #-----------------------------------------------------------------------
    generator = QCGridPropertyGenerator.FromSystem (system )
    generator.DefineGrid    ( gridSpacing = _GridSpacing ) # . Some value in atomic units - e.g. 0.2
    

    
    orbital_iso = {}
    if _type == 'orbital':
        generator.GridOrbitals  ( [ key ]    ,       tag = _OrbitalTag    ) # . List of orbital indices (can be one only)    
        generator.Isosurface    ( _OrbitalTag, _isovalue, tag = _IsosurfaceTag )
    
        surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
        isosurface_p = surfaceProperty.isosurface # . This is the surface you can display.
        
        vertices, colors, indexes = surface_parser ( surface = isosurface_p , iso_color = [1,0,0] )
        
        orbital_iso['obital_plus'] = [vertices, colors, indexes]
        generator.Isosurface    ( _OrbitalTag, _isovalue*-1, tag = _IsosurfaceTag )
        surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
        isosurface_n = surfaceProperty.isosurface # . This is the surface you can display.
        
        vertices, colors, indexes = surface_parser ( surface = isosurface_n , iso_color = [0,0,1] )
        orbital_iso['obital_minus'] = [vertices, colors, indexes]
    
    elif _type == 'potential':
        generator.GridPotential  (                   tag = _IsosurfaceTag   ) # . List of orbital indices (can be one only)    
        generator.Isosurface    ( 'potential', _isovalue, tag = _IsosurfaceTag)
        surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
        isosurface_p = surfaceProperty.isosurface # . This is the surface you can display.
        
        vertices, colors, indexes = surface_parser ( surface = isosurface_p , iso_color = [1,0,0] )
        
        orbital_iso['obital_plus'] = [vertices, colors, indexes]
        generator.Isosurface    ( _OrbitalTag, _isovalue*-1, tag = _IsosurfaceTag )
        surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
        isosurface_n = surfaceProperty.isosurface # . This is the surface you can display.
        
        vertices, colors, indexes = surface_parser ( surface = isosurface_n , iso_color = [0,0,1] )
        orbital_iso['obital_minus'] = [vertices, colors, indexes]
    
    
    
    
    
    elif _type == 'density':
        generator.GridDensity  (                   tag = 'density'   ) # . List of orbital indices (can be one only)    
        generator.Isosurface    ( 'density', _isovalue, tag = _IsosurfaceTag)
        surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
        isosurface_p = surfaceProperty.isosurface # . This is the surface you can display.
        vertices, colors, indexes = surface_parser ( surface = isosurface_p , iso_color = [0,0,1] )
        orbital_iso['obital_plus'] = [vertices, colors, indexes]
        
    else:
        pass
        
   
    #vertices, colors, indexes = surface_parser ( surface = isosurface_p , iso_color = [1,0,0] )
    #
    #orbital_iso['obital_plus'] = [vertices, colors, indexes]
    #
    #generator.Isosurface    ( _OrbitalTag, _isovalue*-1, tag = _IsosurfaceTag )
    #surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
    #isosurface_n = surfaceProperty.isosurface # . This is the surface you can display.
    #
    #vertices, colors, indexes = surface_parser ( surface = isosurface_n , iso_color = [0,0,1] )
    #orbital_iso['obital_minus'] = [vertices, colors, indexes]
    
    return orbital_iso
    
    #generator.DefineGrid    ( gridSpacing = _GridSpacing ) # . Some value in atomic units - e.g. 0.2
    #
    #orbital_iso = {}
    #print ('key, _OrbitalTag:', key, _OrbitalTag)
    #generator.GridOrbitals  ( [ key ], tag = _OrbitalTag) # . List of orbital indices (can be one only)    
    #
    #generator.Isosurface    ( _OrbitalTag, _isovalue, tag = _IsosurfaceTag )
    #surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
    #isosurface_p = surfaceProperty.isosurface # . This is the surface you can display.
    #
    #
    #
    #vertices, colors, indexes = surface_parser ( surface = isosurface_p , iso_color = [1,0,0] )
    #
    #orbital_iso['obital_plus'] = [vertices, colors, indexes]
    #
    #generator.Isosurface    ( _OrbitalTag, _isovalue*-1, tag = _IsosurfaceTag )
    #surfaceProperty = generator.GetProperty ( _IsosurfaceTag )
    #isosurface_n = surfaceProperty.isosurface # . This is the surface you can display.
    #
    #vertices, colors, indexes = surface_parser ( surface = isosurface_n , iso_color = [0,0,1] )
    #orbital_iso['obital_minus'] = [vertices, colors, indexes]
    #
    #return orbital_iso



def generate_wavefunction_parallel(job):
    """ Function doc """
    generate_grid_parallel

    i          = job[0]
    system     = job[1]
    coords     = job[2]
    #parameters = job[3]
    
    apply_coords_to_system(system, coords)
    system.Energy()

    orbitalsP = system.scratch.orbitalsP
    energies  = orbitalsP.energies
    
    LUMO      = orbitalsP.occupancyHandler.numberOccupied
    HOMO      = LUMO - 1
    
    #generator = QCGridPropertyGenerator.FromSystem (system )
    #generator.orbitalsP = orbitalsP
    orbitals  = []
    generator = QCGridPropertyGenerator.FromSystem ( system )
    
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
        
    return orbitals, system, generator
    
    
    
    
    
    '''
    for i in range(len(orbitals)):
        reverse_index = -i-1 #- len(orbitals)
        print(reverse_index, orbitals[reverse_index ])
        self.liststore.append(orbitals[reverse_index ])
    '''











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


