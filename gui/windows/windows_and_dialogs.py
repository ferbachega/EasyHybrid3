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
#  self.home,'gui/windows/export_system_window.glade')
#  
import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
from gi.repository import Gdk
from gi.repository import GdkPixbuf
#from GTKGUI.gtkWidgets.filechooser import FileChooser
#from easyhybrid.pDynamoMethods.pDynamo2Vismol import *
from gui.gtk_widgets import FolderChooserButton
from gui.gtk_widgets import VismolTrajectoryFrame

from util.file_parser import get_file_type  
from util.file_parser import read_MOL2  

import gc
import os

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')


#from GTKGUI.gtkWidgets.main_treeview import GtkMainTreeView


class EnergyRefinementWindow():

    def __init__(self, main = None ):
        """ Class initialiser """
        self.main     = main
        self.home     = main.home
        #self.p_session           = self.easyhybrid_main.p_session
        #self.vm_session          = main.vm_session
        self.Visible             =  False        
        
        self.vobject_liststore   = Gtk.ListStore(str, int)
        self.data_liststore      = Gtk.ListStore(str, int)

    def OpenWindow (self, vobject = None):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'gui/windows/energy_refinement_window.glade'))
            self.builder.connect_signals(self)
            #self.vobject = vobject
            self.window = self.builder.get_object('window')
            self.window.set_title('Energy Refinement Window')
            self.window.connect('destroy', self.CloseWindow)
            
            self.builder.get_object('button_cancel').connect('clicked', self.CloseWindow)
            self.builder.get_object('button_run').connect('clicked', self.on_button_run_clicked)
            
            '''--------------------------------------------------------------------------------------------''' 
            self.combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
            #self.starting_coords_liststore = self.main.active_system_liststore
            self.combobox_starting_coordinates.set_model(self.main.vobject_liststore_dict[self.main.p_session.active_id])
            renderer_text = Gtk.CellRendererText()
            self.combobox_starting_coordinates.pack_start(renderer_text, True)
            self.combobox_starting_coordinates.add_attribute(renderer_text, "text", 0)
            size = len(self.main.vobject_liststore_dict[self.main.p_session.active_id])
            self.combobox_starting_coordinates.set_active(size-1)
            '''--------------------------------------------------------------------------------------------''' 
            
            #'''--------------------------------------------------------------------------------------------------
            self.folder_chooser_button2 = FolderChooserButton(main = self.main, sel_type = 'file')
            self.builder.get_object('folder_chooser_box2').pack_start(self.folder_chooser_button2.btn, True, True, 0)
            #'''--------------------------------------------------------------------------------------------------
            
            #'''--------------------------------------------------------------------------------------------------
            self.folder_chooser_button1 = FolderChooserButton(main =  self.main)
            self.builder.get_object('folder_chooser_box1').pack_start(self.folder_chooser_button1.btn, True, True, 0)
            #'''--------------------------------------------------------------------------------------------------
            
            
            
            self.builder.get_object('combobox_coordinate_type').set_active(0)
            #'''--------------------------------------------------------------------------------------------------
            self.builder.get_object('box_reaction_coordinate').set_sensitive(False)
            self.builder.get_object('label_CPUs').set_sensitive(False)
            self.builder.get_object('n_CPUs_spinbutton').set_sensitive(False)
            self.builder.get_object('label_input_logfile').set_sensitive(False)
            self.builder.get_object('file_chooser_btn_logfile').set_sensitive(False)
            self.builder.get_object('folder_chooser_box2').set_sensitive(False)
            self.builder.get_object('combobox_coordinate_type').set_sensitive(False)
            self.builder.get_object('label_coordinate_type').set_sensitive(False)
            self.builder.get_object('label_file_or_folder').set_sensitive(False)
            self.builder.get_object('label_input_logfile').set_sensitive(False)
            self.builder.get_object('frame_output').set_sensitive(False)
            #'''--------------------------------------------------------------------------------------------------
            
            self.on_combobox_coordinate_type_changed(None)
 
            self.window.show_all()
            self.Visible  = True
    
        else:
            pass

    def on_combobox_coordinate_type_changed (self, widget):
        """ Function doc """
        if self.builder.get_object('combobox_coordinate_type').get_active() == 0:
            
            #self.builder.get_object('folder_chooser_box2').set_sensitive(False)
            #self.builder.get_object('label_file_or_folder').set_sensitive(False)
            #
            #self.builder.get_object('combobox_starting_coordinates').set_sensitive(True)
            #self.builder.get_object('label_coordinates').set_sensitive(True)
            self.builder.get_object('folder_chooser_box2').hide()
            self.builder.get_object('label_file_or_folder').hide()

            self.builder.get_object('combobox_starting_coordinates').show()
            self.builder.get_object('label_coordinates').show()
        
        else:
            #self.builder.get_object('combobox_starting_coordinates').set_sensitive(False)
            #self.builder.get_object('label_coordinates').set_sensitive(False)
            #
            #self.builder.get_object('folder_chooser_box2').set_sensitive(True)
            #self.builder.get_object('label_file_or_folder').set_sensitive(True)
            
            self.builder.get_object('combobox_starting_coordinates').hide()
            self.builder.get_object('label_coordinates').hide()
            
            self.builder.get_object('folder_chooser_box2').show()
            self.builder.get_object('label_file_or_folder').show()

    def on_coordinates_combobox_change (self, widget):
        """ Function doc """
        _id = self.coordinates_combobox.get_active()
        print(_id)
        vobject_index = None
        #-----------------------------------------------------------------------------
        _iter = self.coordinates_combobox.get_active_iter()
        if _iter is not None:
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.coordinates_combobox.get_model()
            _name, vobject_index = model[_iter][:2]
            print ('\n\n\_name, vobject_index:', _name, vobject_index, '\n\n')
        #-----------------------------------------------------------------------------
        self.vobject = self.main.vm_session.vm_objects_dic[vobject_index]

        self.data_liststore.clear()
        for index , data in enumerate(self.main.p_session.systems[self.vobject.easyhybrid_system_id]['logfile_data'][vobject_index]):
            #print(data)
            self.data_liststore.append([data['name'], index])
        
        
        #self.data_liststore.append(['all', 2])
        
        self.data_combobox.set_active(0)

    def on_data_combobox_change (self, widget):
        """ Function doc """

        _iter = self.data_combobox.get_active_iter()
        if _iter is not None:
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.data_combobox.get_model()
            _name, index = model[_iter][:2]
            #print ('\n\n\_name, index:', _name,  index, '\n\n')
        
        #self.vobject = self.main.vm_session.vm_objects_dic[vobject_index]
        self.data = self.main.p_session.systems[self.vobject.easyhybrid_system_id]['logfile_data'][self.vobject.index][index] 
        print(self.data)
        self._draw_data(cla = True)

    def radiobutton_single_point_toggled_cb (self, widget):
        """ Function doc """
        if self.builder.get_object('radiobutton_single_point').get_active():
            self.builder.get_object('box_reaction_coordinate').set_sensitive(False)
            self.builder.get_object('label_CPUs').set_sensitive(False)
            self.builder.get_object('n_CPUs_spinbutton').set_sensitive(False)
            self.builder.get_object('label_input_logfile').set_sensitive(False)
            self.builder.get_object('file_chooser_btn_logfile').set_sensitive(False)
            self.builder.get_object('folder_chooser_box2').set_sensitive(False)
            self.builder.get_object('combobox_coordinate_type').set_sensitive(False)
            self.builder.get_object('label_coordinate_type').set_sensitive(False)
            self.builder.get_object('label_file_or_folder').set_sensitive(False)
            self.builder.get_object('label_input_logfile').set_sensitive(False)
            self.builder.get_object('frame_output').set_sensitive(False)
            
            self.builder.get_object('combobox_coordinate_type').set_active(0)
            
        else:
            self.builder.get_object('box_reaction_coordinate').set_sensitive(True)
            self.builder.get_object('label_CPUs').set_sensitive(True)
            self.builder.get_object('n_CPUs_spinbutton').set_sensitive(True)
            self.builder.get_object('label_input_logfile').set_sensitive(True)
            self.builder.get_object('file_chooser_btn_logfile').set_sensitive(True)
            self.builder.get_object('folder_chooser_box2').set_sensitive(True)
            self.builder.get_object('combobox_coordinate_type').set_sensitive(True)
            self.builder.get_object('label_coordinate_type').set_sensitive(True)
            self.builder.get_object('label_file_or_folder').set_sensitive(True)
            self.builder.get_object('label_input_logfile').set_sensitive(True)
            self.builder.get_object('frame_output').set_sensitive(True)
            
            self.on_combobox_coordinate_type_changed (None)

    def change_check_button_reaction_coordinate (self, widget):
        """ Function doc """
        #radiobutton_bidimensional = self.builder.get_object('radiobutton_bidimensional')
        if self.builder.get_object('label_check_button_reaction_coordinate2').get_active():
            self.box_reaction_coordinate2.set_sensitive(True)
            self.builder.get_object('n_CPUs_spinbutton').set_sensitive(True)
            self.builder.get_object('n_CPUs_label').set_sensitive(True)
        else:
            self.box_reaction_coordinate2.set_sensitive(False)
            self.builder.get_object('n_CPUs_spinbutton').set_sensitive(False)
            self.builder.get_object('n_CPUs_label')     .set_sensitive(False)

    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False

    def on_button_run_clicked (self, button):
        """ Function doc """
        
        #----------------------------------------------------------------------------------------------
        #                            S I N G L E    P O I N T 
        if self.builder.get_object('radiobutton_single_point').get_active():
            
            combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
            tree_iter = combobox_starting_coordinates.get_active_iter()
            if tree_iter is not None:
                
                '''selecting the vismol object from the content that is in the combobox '''
                model = combobox_starting_coordinates.get_model()
                name, vobject_id = model[tree_iter][:2]
                vobject = self.main.vm_session.vm_objects_dic[vobject_id]
                
                '''This function imports the coordinates of a vobject into the dynamo system in memory.''' 
                print('vobject:', vobject.name, len(vobject.frames) )
                self.main.p_session.get_coordinates_from_vobject_to_pDynamo_system(vobject)

            energy = self.main.p_session.get_energy()
            
            self.CloseWindow ( button =None, data  = None)
            
            dialog = EasyHybridDialogEnergy(parent = self.window, energy = energy, name = vobject.name)
            response = dialog.run()
            dialog.destroy()
        #----------------------------------------------------------------------------------------------

        else:
            
            
            if self.builder.get_object('combobox_coordinate_type').get_active() == 0:
                
                #data = {
                #   'name': "energy refinement",
                #   'type': "plot2D",
                #   'RC1' : RC1,
                #   'RC2' : RC2,
                #   'Z'   : Z
                #   }
                
                
                combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
                tree_iter = combobox_starting_coordinates.get_active_iter()
                if tree_iter is not None:
                    
                    '''selecting the vismol object from the content that is in the combobox '''
                    model = combobox_starting_coordinates.get_model()
                    name, vobject_id = model[tree_iter][:2]
                    vobject = self.main.vm_session.vm_objects_dic[vobject_id]
                    
                
                for frame_number in range(0, len(vobject.frames)):
                    
                    '''This function imports the coordinates of a vobject into the dynamo system in memory.''' 
                    #print('vobject:', vobject.name, len(vobject.frames) )
                
                    self.main.p_session.get_coordinates_from_vobject_to_pDynamo_system(vobject = vobject, system_id =  None, frame = frame_number)

                    energy = self.main.p_session.get_energy(log = True)
                    
                    #print(vobject.idx_2D_xy, energy)
                    #except:
                    print(frame_number, energy)
                
                self.CloseWindow ( button =None, data  = None)

    def update (self):
        """ Function doc """
        #print('VismolGoToAtomWindow2 update')
        pass


class ExportDataWindow:
    """ Class doc """
    def __init__(self, main = None):
        """ Class initialiser """
        self.main    = main
        self.home    = main.home
        self.Visible =  False        
        self.starting_coords_liststore = Gtk.ListStore(str, int)
        self.is_single_frame  = True
        
    def OpenWindow (self, sys_selected = None):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'gui/windows/export_system_window.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('export_data_window')
            self.window.set_title('Export Data Window')
            self.window.set_keep_above(True)
            '''--------------------------------------------------------------------------------------------'''
            
            
            '''--------------------------------------------------------------------------------------------'''
            self.format_store = Gtk.ListStore(str)
            self.formats = {
                       0 : 'pkl - pDynamo system'         ,
                       1 : 'pkl - pDynamo coordinates'    ,
                       2 : 'pdb - Protein Data Bank'      ,
                       3 : 'xyz - cartesian coordinates'  ,
                       4 : 'mol'                          ,
                       5 : 'mol2'                         ,
                       }
                       
            for key, _format in self.formats.items():
                self.format_store.append([_format])
                #print (format)
            self.combobox_fileformat = self.builder.get_object('combobox_fileformat')
            self.combobox_fileformat.set_model(self.format_store)
            #self.formats_combo.connect("changed", self.on_name_combo_changed)
            self.combobox_fileformat.set_model(self.format_store)
            
            renderer_text = Gtk.CellRendererText()
            self.combobox_fileformat.pack_start(renderer_text, True)
            self.combobox_fileformat.add_attribute(renderer_text, "text", 0)
            '''--------------------------------------------------------------------------------------------'''
            
            
            
            '''--------------------------------------------------------------------------------------------'''
            #self.psystem_liststore = Gtk.ListStore(str,int)
            #names = [ ]
            #for key , system in self.main.p_session.psystem.items():
            #    try:
            #        name = system.label
            #        print ([name, int(key)])
            #        self.psystem_liststore.append([name, int(key)])
            #    except:
            #        print(system)
            self.psystem_combo = self.builder.get_object('combobox_pdynamo_system')
            self.psystem_combo.set_model(self.main.system_liststore)
            #self.psystem_combo.connect("changed", self.on_name_combo_changed)
            #self.psystem_combo.set_model(self.psystem_liststore)
            
            renderer_text = Gtk.CellRendererText()
            self.psystem_combo.pack_start(renderer_text, True)
            self.psystem_combo.add_attribute(renderer_text, "text", 0)
            if sys_selected:
                self.psystem_combo.set_active(sys_selected)
            else:
                #sys_selected = self.main
                self.psystem_combo.set_active(0)
                _, sys_selected = self.main.system_liststore[0]
                #sys_selected = self.main.system_liststore
            '''--------------------------------------------------------------------------------------------'''
            
            
            #'''--------------------------------------------------------------------------------------------------
            self.folder_chooser_button = FolderChooserButton(main =  self.main)
            self.builder.get_object('folder_chooser_box').pack_start(self.folder_chooser_button.btn, True, True, 0)
            #'''--------------------------------------------------------------------------------------------------

            #------------------------------------------------------------------------------------
            
            self.combobox_starting_coordinates = self.builder.get_object('vobjects_combobox')
            self.combobox_starting_coordinates.set_model(self.main.vobject_liststore_dict[sys_selected])
            
            renderer_text = Gtk.CellRendererText()
            self.combobox_starting_coordinates.pack_start(renderer_text, True)
            self.combobox_starting_coordinates.add_attribute(renderer_text, "text", 0)
            
            size = len(self.main.vobject_liststore_dict[sys_selected])
            self.combobox_starting_coordinates.set_active(size-1)
            #------------------------------------------------------------------------------------
        
            #------------------------------------------------------------------------------------
            self.folder_chooser_button.set_folder(HOME)
            #------------------------------------------------------------------------------------

            self.combobox_fileformat.set_active(0)
            self.window.show_all()
            self.Visible  = True
    
    def combobox_pdynamo_system (self, widget):
        """ Function doc """
        #print('combobox_pdynamo_system aqui')
        if self.Visible:
            tree_iter = self.psystem_combo.get_active_iter()
            if tree_iter is not None:
                
                '''selecting the vismol object from the content that is in the combobox '''
                model = self.psystem_combo.get_model()
                name, sys_id = model[tree_iter][:2]
                print (name, sys_id)

            self.combobox_starting_coordinates.set_model(self.main.vobject_liststore_dict[sys_id])
            self.combobox_starting_coordinates.set_active(0)
    
    def on_vobject_combo_changed (self, widget):
        '''this combobox has the reference to the starting coordinates of a simulation'''
        #combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]
            print (name, model[tree_iter][:2])
            #name, vobject_id = model[tree_iter][:2]
        
        
        if len(self.main.vm_session.vm_objects_dic[vobject_id].frames) > 1:
            print(self.main.vm_session.vm_objects_dic[vobject_id].name,
                  len(self.main.vm_session.vm_objects_dic[vobject_id].frames),'True')
            self.is_single_frame = True
            
            if self.combobox_fileformat.get_active( ) == 0:
                self.builder.get_object('entry_first').set_sensitive(False)
                self.builder.get_object('label_first').set_sensitive(False)
                self.builder.get_object('entry_stride').set_sensitive(False)
                self.builder.get_object('label_stride').set_sensitive(False)
            else:
                self.builder.get_object('entry_first').set_sensitive(True)
                self.builder.get_object('label_first').set_sensitive(True)            
                self.builder.get_object('entry_stride').set_sensitive(True)
                self.builder.get_object('label_stride').set_sensitive(True)
            
            
        else:
            print(self.main.vm_session.vm_objects_dic[vobject_id].name,
                  len(self.main.vm_session.vm_objects_dic[vobject_id].frames),'False')
            
            self.is_single_frame = False
            self.builder.get_object('entry_first').set_sensitive(False)
            self.builder.get_object('label_first').set_sensitive(False)
            self.builder.get_object('entry_stride').set_sensitive(False)
            self.builder.get_object('label_stride').set_sensitive(False)
        
    def on_combobox_fileformat (self, widget):
        """ Function doc """
        print('on_combobox_fileformat')
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]
            print (name, model[tree_iter][:2])
            
        if len(self.main.vm_session.vm_objects_dic[vobject_id].frames) > 1:
            print(self.main.vm_session.vm_objects_dic[vobject_id].name,
                  len(self.main.vm_session.vm_objects_dic[vobject_id].frames),'True')
            self.is_single_frame = True
            if self.combobox_fileformat.get_active( ) == 0:
                self.builder.get_object('entry_first').set_sensitive(False)
                self.builder.get_object('label_first').set_sensitive(False)
                self.builder.get_object('entry_stride').set_sensitive(False)
                self.builder.get_object('label_stride').set_sensitive(False)
            else:
                self.builder.get_object('entry_first').set_sensitive(True)
                self.builder.get_object('label_first').set_sensitive(True)            
                self.builder.get_object('entry_stride').set_sensitive(True)
                self.builder.get_object('label_stride').set_sensitive(True)
        else:
            print(self.main.vm_session.vm_objects_dic[vobject_id].name,
                  len(self.main.vm_session.vm_objects_dic[vobject_id].frames),'False')
            
            self.is_single_frame = False
            self.builder.get_object('entry_first').set_sensitive(False)
            self.builder.get_object('label_first').set_sensitive(False)
            self.builder.get_object('entry_stride').set_sensitive(False)
            self.builder.get_object('label_stride').set_sensitive(False)
   
    def on_name_combo_changed (self, widget):
        """ Function doc """
        if  self.combox.get_active() == 0:
            self.folder_chooser_button.sel_type = 'folder'
        else:
            self.folder_chooser_button.sel_type = 'file'

    def on_radio_button_change (self, widget):
        """ Function doc """
        if self.builder.get_object('radiobutton_singlefile').get_active():
            self.is_single_file = True
        else:
            self.is_single_file = False

    def export_data (self, button):
        """ Function doc """
        
        parameters = {
                      'system_id'       : None,
                      'vobject_id'      : None, 
                      'format'          : None,
                      'is_single_file'  : True,
                      
                      'fist'            : 0   ,
                      'last'            : -1  ,
                      'stride'          : 1   ,
                      
                      'filename'        :'exported_system',
                      'folder'          :None, 
                      
                      }
        
        
        '''--------------------------------------------------------------------------'''
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]
        parameters['vobject_id'] = vobject_id
        '''--------------------------------------------------------------------------'''


        
        '''---------------------  Getting the correct format  ---------------------------'''
        _format  = self.combobox_fileformat.get_active()
        parameters['format'] = _format
        '''------------------------------------------------------------------------------'''
        
        
        
        '''------------------------   Fomrat  ------------------------------------------'''
        folder   = self.folder_chooser_button.get_folder()
        filename = self.builder.get_object('entry_filename').get_text()
        parameters['folder'] = folder
        parameters['filename'] = filename
        '''------------------------------------------------------------------------------'''
        
        
        
        '''------------------------  pDynam System   ------------------------------------'''
        tree_iter = self.psystem_combo.get_active_iter()
        if tree_iter is not None:
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.psystem_combo.get_model()
            name, system_id = model[tree_iter][:2]
        parameters['system_id'] = system_id
        '''------------------------------------------------------------------------------'''
        
        
        '''--------------------------- Is single file ? ---------------------------------'''
        if self.builder.get_object('radiobutton_singlefile').get_active():
            parameters['is_single_file'] = True          
        else:
            parameters['is_single_file'] = False
        '''------------------------------------------------------------------------------'''
        
        
        
        '''----------------------   FIRST  LAST and STRIDE   -----------------------------'''
        parameters['first']  = int(self.builder.get_object('entry_first').get_text() )
        parameters['last']   = int(self.builder.get_object('entry_last').get_text()  )
        parameters['stride'] = int(self.builder.get_object('entry_stride').get_text())
        '''-------------------------------------------------------------------------------'''

        print(parameters)
        '''------------------------------------------------------------------------------'''
        self.main.p_session.export_system (parameters)
        '''------------------------------------------------------------------------------'''
    
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False

    def update (self):
        """ Function doc """
        #print('VismolGoToAtomWindow2 update')
        pass


class PDynamoSelectionWindow:
    """ Class doc """
    def __init__(self, main = None):
        """ Class initialiser """
        self.main = main
        self.vm_session      = main.vm_session
        self.Visible         = False        
        self.home            = main.home
        self.p_session       = main.p_session
        
        
        self.chain = ''
        self.resn  = ''
        self.resi  = ''
        self.atom  = ''

 
    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home, 'gui/windows/pDynamo_selection.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('pdynamo_selection_window')
            self.window.set_title('pDynamo Selection Window')
            self.window.set_keep_above(True)
            
            #self.chain_entry  = self.builder.get_object('chain_entry')
            #self.resn_entry   = self.builder.get_object('resn_entry' )
            #self.resi_entry   = self.builder.get_object('resi_entry' )
            #self.atom_entry   = self.builder.get_object('atom_entry' )
            self.builder.get_object('chain_entry').set_text(str(self.chain) )
            self.builder.get_object('resn_entry' ).set_text(str(self.resn ) )
            self.builder.get_object('resi_entry' ).set_text(str(self.resi ) )
            self.builder.get_object('atom_entry' ).set_text(str(self.atom ) )
            
            
            self.box_combo_methods = self.builder.get_object('box_combo_methods')
            
            #self.self.method_combobox   = self.builder.get_object('self.method_combobox'   )
            #self.radius_spinbutton = self.builder.get_object('radius_spinbutton' )
            #self.action_combobox   = self.builder.get_object('action_combobox'   )

            method_store = Gtk.ListStore(str)
            method_store.append(["ByComponent"])
            method_store.append(["Complement"])
            method_store.append(["ByAtom"])


            #------------------------------------------------------------------
            self.method_combo = Gtk.ComboBox.new_with_model(method_store)
            renderer_text = Gtk.CellRendererText()
            self.method_combo.pack_start(renderer_text, True)
            self.method_combo.add_attribute(renderer_text, "text", 0)            
            
            self.method_combo.set_entry_text_column(0)
            self.method_combo.set_active(0)
            self.box_combo_methods.pack_start(self.method_combo, False, False, 0)
            #------------------------------------------------------------------

            
            
            
            
            self.radius_spinbutton  = self.builder.get_object('radius_spinbutton' )
            #------------------------------------------------------------------
            self.fps_adjustment = Gtk.Adjustment(value          = 24 , 
                                                 upper          = 100, 
                                                 step_increment = 1  , 
                                                 page_increment = 10 )

            self.radius_spinbutton.set_adjustment ( self.fps_adjustment)
            #------------------------------------------------------------------




            '''
            self.box_combo_action = self.builder.get_object('box_combo_action' )
            
            action_store = Gtk.ListStore(str)
            action_store.append(["ByComponent"])
            action_store.append(["Complement"])
            action_store.append(["ByLinearPolymer"])
            #------------------------------------------------------------------
            self.action_combo = Gtk.ComboBox.new_with_model(action_store)
            renderer_text = Gtk.CellRendererText()
            self.action_combo.pack_start(renderer_text, True)
            self.action_combo.add_attribute(renderer_text, "text", 0)            
            
            self.action_combo.set_entry_text_column(0)
            self.action_combo.set_active(0)
            
            self.box_combo_action.pack_start(self.action_combo, False, False, 0)
            #------------------------------------------------------------------
            '''
                
                
                

            self.window.show_all()
            self.Visible  = True
    
        else:
            self.window.present()
    
    def import_data (self, button):
        """ Function doc """
        entry_name    = None
        idnum     = self.combobox_pdynamo_system.get_active()
        text      = self.combobox_pdynamo_system.get_active_text()
        
        #print(idnum, text )
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False
    
    
    def run_selection (self, button):
        """ Function doc """
        #print('run_selection', self.method_combo.get_active(), self.radius_spinbutton.get_value ())
        self.chain = self.builder.get_object('chain_entry').get_text()
        self.resn  = self.builder.get_object('resn_entry' ).get_text()
        self.resi  = self.builder.get_object('resi_entry' ).get_text()
        self.atom  = self.builder.get_object('atom_entry' ).get_text()
        
        #print(chain,resn,resi, atom)
        #'''
        #atom1 = self.vm_session.picking_selections.picking_selections_list[0]
        #print (atom1.chain, atom1.resn, atom1.resi, atom1.name)
        #print ("%s:%s.%s:%s" %(chain,resn,resi,atom))
        
        _centerAtom = "%s:%s.%s:%s" %(self.chain, 
                                      self.resn,
                                      self.resi,
                                      self.atom)
        _radius     =  self.radius_spinbutton.get_value ()
        _method     =  self.method_combo.get_active()
        
        #self.main.p_session.selections (_centerAtom, _radius, _method )
        indexes = self.main.p_session.p_selections (system_id = None, 
                                                  _centerAtom = _centerAtom, 
                                                      _radius = _radius, 
                                                      _method = _method)
        
        #print (indexes)
        self.main.vm_session.selections[self.vm_session.current_selection].selecting_by_indexes(self.atom1.vm_object, 
                                                           indexes, 
                                                           clear=True)
        if self.main.vm_session.selection_box_frane:
            self.main.vm_session.selection_box_frane.change_toggle_button_selecting_mode_status(False)
        else:
            self.main.vm_session._picking_selection_mode = False
        
        self.main.vm_session.selections[self.main.vm_session.current_selection].active = True
        self.main.vm_session.vm_glcore.queue_draw()
        
        
    def get_info_from_selection (self, button):
        """ Function doc """
        #chain = self.builder.get_object('chain_entry').get_text()
        #resn  = self.builder.get_object('resn_entry' ).get_text()
        #resi  = self.builder.get_object('resi_entry' ).get_text()
        #atom  = self.builder.get_object('atom_entry' ).get_text()
        
        atom1 = self.vm_session.picking_selections.picking_selections_list[0]
        self.atom1 = atom1
        if atom1:
            #atom1.chain, atom1.resn, atom1.resi, atom1.name
            
            if atom1.chain.name =='':
                chain = '*'
            else:
                chain = atom1.chain.name
            #residue = atom1.vm_object.residues[atom1.residue]
            self.builder.get_object('chain_entry').set_text(str(chain )      )
            self.builder.get_object('resn_entry' ).set_text(str(atom1.residue.name) )
            self.builder.get_object('resi_entry' ).set_text(str(atom1.residue.index) )
            self.builder.get_object('atom_entry' ).set_text(str(atom1.name) )
        else:
            print('use picking selection to chose the central atom')
            
    def update (self):
        """ Function doc """
        #print('VismolGoToAtomWindow2 update')
        pass
    
            
class EasyHybridSetupQCModelWindow:
    """ Class doc """
    
    def __init__(self, main = None):
        """ Class initialiser """
        self.main_session     = main
        self.home             = main.home
        self.Visible          = False        
        self.vismol_object    = None 
        
        self.methods_liststore = Gtk.ListStore(str, str, str)
        
        self.method_id    = 0    # 0 for am1, 1 for pm3  and...
        self.charge       = 0
        self.multiplicity = 1
        self.restricted   = True
        
        self.adjustment_charge = Gtk.Adjustment(value=0,
                                           lower=-100,
                                           upper=100,
                                           step_increment=1,
                                           page_increment=1,
                                           page_size=0)
        
        self.adjustment_multiplicity = Gtk.Adjustment(value=1,
                                           lower=1,
                                           upper=100,
                                           step_increment=1,
                                           page_increment=1,
                                           page_size=0)
        
        self.methods_id_dictionary = {
                                      0 : 'am1'     ,
                                      1 : 'am1dphot',
                                      2 : 'pm3'     ,
                                      3 : 'pm6'     ,
                                      4 : 'mndo'    ,
                                
                                      }
        
      
    def refresh (self, option = 'all'):
        """ Function doc """
        self.update_number_of_qc_atoms()
       
        
    def OpenWindow (self, vismol_object = None):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home, 'gui/windows/setup_qc_model_window2.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('SetupQCModelWindow')
            self.window.set_keep_above(True)

            self.vismol_object = vismol_object
            
            '''--------------------------------------------------------------------------------------------'''
            self.methods_type_store = Gtk.ListStore(str)
            methods_types = [
                "am1",
                "am1dphot",
                "pm3",
                "pm6",
                "mndo",
                "ab initio - ORCA",
                "DFT / DFTB",
                ]
            for method_type in methods_types:
                self.methods_type_store.append([method_type])
                #print (method_type)
            
            self.methods_combo = self.builder.get_object('QCModel_methods_combobox')
            self.methods_combo.connect("changed", self.on_name_combo_changed)
            self.methods_combo.set_model(self.methods_type_store)
            renderer_text = Gtk.CellRendererText()
            self.methods_combo.pack_start(renderer_text, True)
            self.methods_combo.add_attribute(renderer_text, "text", 0)
            self.methods_combo.set_active(self.method_id)
            '''--------------------------------------------------------------------------------------------'''


            self.spinbutton_charge       = self.builder.get_object('spinbutton_charge'      )
            self.spinbutton_multiplicity = self.builder.get_object('spinbutton_multiplicity')
            self.spinbutton_charge      .set_adjustment(self.adjustment_charge)
            self.spinbutton_multiplicity.set_adjustment(self.adjustment_multiplicity)
            
            self.window.show_all()                                               
            self.builder.connect_signals(self)                                   
            
            
            ''' Updating the number of atoms '''
            self.update_number_of_qc_atoms ()
            
            self.Visible  =  True

        else:
            ''' Updating the number of atoms '''
            self.update_number_of_qc_atoms ()
            self.window.present()
            
    def update_number_of_qc_atoms (self):
        """ Function doc """
        self.entry_number_of_qc_atoms = self.builder.get_object('entry_number_of_qc_atoms')
        
        '''   Estiamting the QC charge  '''
        '''----------------------------------------------------------------------------------------------'''
        system  = self.main_session.p_session.psystem[self.main_session.p_session.active_id]
        estimated_charge = 0.0
        
        if len(system.e_qc_table) > 0:
            
            for index in system.e_qc_table:
                estimated_charge += system.mmState.charges[index]
            
            estimated_charge = int(round(estimated_charge))
            self.spinbutton_charge.set_value (estimated_charge)
            
            number_of_qc_atoms = len(system.e_qc_table)
            self.entry_number_of_qc_atoms.set_text(str(number_of_qc_atoms))
        
        else:
            
            if system.mmModel:
                for charge in system.mmState.charges:
                    estimated_charge += charge
                
                estimated_charge = int(round(estimated_charge))
                self.spinbutton_charge.set_value (estimated_charge)
            else:
                print('mmModel not found, pure QC system.')
            number_of_qc_atoms = len(system.atoms)
            self.entry_number_of_qc_atoms.set_text(str(number_of_qc_atoms)+ ' (all)')
        '''----------------------------------------------------------------------------------------------'''

    
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False
    
    #----------------------------------------------------------------
    def on_spian_button_change (self, widget):
        """ Function doc """
        self.charge       = self.spinbutton_charge.get_value_as_int()
        self.multiplicity = self.spinbutton_multiplicity.get_value_as_int()
        
    def on_name_combo_changed (self, combobox):
        """ Function doc """
        self.method_id = self.builder.get_object('QCModel_methods_combobox').get_active()
    
    def on_button_ok (self, button):
        """ Function doc """
        
        if self.builder.get_object('radio_button_restricted').get_active():
            #print("%s is active" % (self.builder.get_object('radio_button_restricted').get_label()))
            self.restricted = True
        else:
            #print("%s is not active" % (self.builder.get_object('radio_button_restricted').get_label()))
            self.restricted = False
        

        
        parameters = {
                    'qcmodel'          : 'mndo'           ,
                    'charge'           : self.charge      ,
                    'multiplicity'     : self.multiplicity,
                    'method'           : self.methods_id_dictionary[self.method_id]   ,
                    'isSpinRestricted' : self.restricted  ,
                     }


        parameters['energyTolerance'  ] = float(self.builder.get_object('entry_energyTolerance').get_text())
        parameters['densityTolerance' ] = float(self.builder.get_object('entry_densityTolerance').get_text())
        parameters['maximumIterations'] = int(self.builder.get_object('entry_maximumIterations').get_text())


        self.main_session.p_session.define_a_new_QCModel(system = None,  parameters = parameters, vismol_object =  self.vismol_object)
        #self.main_session.update_gui_widgets ()
        self.window.destroy()
        self.Visible    =  False

    def update (self):
        """ Function doc """
        #print('VismolGoToAtomWindow2 update')
        pass


class EasyHybridGoToAtomWindow(Gtk.Window):
    def __init__(self, main = None, system_liststore = None):
        """ Class initialiser """
        
        self.main  = main
        self.vm_session    = self.main.vm_session
        self.p_session     = self.main.p_session
        self.system_liststore      = system_liststore
        self.coordinates_liststore = Gtk.ListStore(str, int)
        
        
        self.residue_liststore = Gtk.ListStore(GdkPixbuf.Pixbuf, int, str, str, int)
        self.atom_liststore    = Gtk.ListStore(GdkPixbuf.Pixbuf, int, str, str, float, int, )
        self.residue_filter    = False
        self.visible           = False



        self.residues_dictionary = {
                               'WAT': [165,42,42], 
                               'SOL': [165,42,42], 
                               'HOH': [165,42,42], 
                               
                               'CYS': [207,40,168], 
                               'ASP': [128,0,128], 
                               'SER': [0, 128, 0], 
                               'GLN': [0, 128, 0], 
                               'LYS': [255, 0, 0],
                               'ILE': [30, 144, 255], 
                               'PRO': [255, 255,0], 
                               'THR': [0, 128, 0], 
                               'PHE': [6,176,176], 
                               'ASN': [0, 128, 0], 
                               'GLY': [255, 165,0], 
                               'HIS': [6,176,176], 
                               
                               # amber
                               "HID": [6,176,176],
                               "HIE": [6,176,176],
                               "HIP": [6,176,176],
                               "ASH": [0, 128, 0],
                               "GLH": [0, 128, 0],
                               "CYX": [207,40,168],
                               
                               # charmm
                               "HSD": [6,176,176], 
                               "HSE": [6,176,176], 
                               "HSP": [6,176,176], 
                               
                               
                               'LEU': [30, 144, 255], 
                               'ARG': [255, 0, 0], 
                               'TRP': [6,176,176], 
                               'ALA': [30, 144, 255], 
                               'VAL': [30, 144, 255], 
                               'GLU': [128,0,128], 
                               'TYR': [6,176,176], 
                               'MET': [30, 144, 255]}



    def OpenWindow (self):
        """ Function doc """
        if self.visible  ==  False:
            
            #self.vm_session.Vismol_Objects_ListStore
            
            #------------------------------------------------------------------#
            #                  SYSTEM combobox and Label
            #------------------------------------------------------------------#
            self.box_vertical    = Gtk.Box(orientation = Gtk.Orientation.VERTICAL,   spacing = 10)
            self.box_horizontal1 = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 10)
             

            
            self.label1  = Gtk.Label()
            self.label1.set_text('System:')
            self.box_horizontal1.pack_start(self.label1, False, False, 0)

            
            self.combobox_systems = Gtk.ComboBox.new_with_model(self.system_liststore)
            self.combobox_systems.connect("changed", self.on_combobox_systems_changed)
            renderer_text = Gtk.CellRendererText()
            self.combobox_systems.pack_start(renderer_text, True)
            self.combobox_systems.add_attribute(renderer_text, "text", 0)
            
            self.box_horizontal1.pack_start(self.combobox_systems, False, False, 0)
            #------------------------------------------------------------------#
            
            
            
            
            
            #------------------------------------------------------------------#
            #                  COORDINATES combobox and Label
            #------------------------------------------------------------------#
            self.label1  = Gtk.Label()
            self.label1.set_text('Coordinates:')
            self.box_horizontal1.pack_start(self.label1, False, False, 0)
            self.coordinates_combobox = Gtk.ComboBox.new_with_model(self.main.vobject_liststore_dict[self.main.p_session.active_id])
            #self.coordinates_combobox.connect("changed", self.on_self.coordinates_combobox_changed)
            renderer_text = Gtk.CellRendererText()
            self.coordinates_combobox.pack_start(renderer_text, True)
            self.coordinates_combobox.add_attribute(renderer_text, "text", 0)
            self.box_horizontal1.pack_start(self.coordinates_combobox, False, False, 0)
            #------------------------------------------------------------------#
            
            


            
            #------------------------------------------------------------------#
            #                  CHAIN combobox and Label
            #------------------------------------------------------------------#
            self.box_horizontal2 = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
            
            
            self.label2  = Gtk.Label()
            self.label2.set_text('Chain:')
            self.box_horizontal2.pack_start(self.label2, False, False, 0)

            self.liststore_chains = Gtk.ListStore(str)
            
            self.combobox_chains = Gtk.ComboBox.new_with_model(self.liststore_chains)
            self.combobox_chains.connect("changed", self.on_combobox_chains_changed)
            renderer_text = Gtk.CellRendererText()
            self.combobox_chains.pack_start(renderer_text, True)
            self.combobox_chains.add_attribute(renderer_text, "text", 0)
            #vbox.pack_start(self.combobox_chains, False, False, True)
            self.box_horizontal2.pack_start(self.combobox_chains, False, False, 0)
            
            
            
            
            #------------------------------------------------------------------#
            #                  RESIDUES combobox and Label
            #------------------------------------------------------------------#
            
            self.label3  = Gtk.Label()
            self.label3.set_text('Residue type:')
            self.box_horizontal2.pack_start(self.label3, False, False, 0)

            self.liststore_residues = Gtk.ListStore(str)
            
            self.combobox_residues = Gtk.ComboBox.new_with_model(self.liststore_residues)
            self.combobox_residues.connect("changed", self.on_combobox_residues_changed)
            renderer_text = Gtk.CellRendererText()
            self.combobox_residues.pack_start(renderer_text, True)
            self.combobox_residues.add_attribute(renderer_text, "text", 0)
            #vbox.pack_start(self.combobox_chains, False, False, True)
            self.box_horizontal2.pack_start(self.combobox_residues, False, False, 0)
            
            
            #------------------------------------------------------------------#
            
            
            
            
            
            self.treeviewbox_horizontal = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
            
            #------------------------------------------------------------------------------------------
            #self.treeview = Gtk.TreeView(model =self.residue_liststore)
            
            
            #-----------------------------------------------------------------------------------------
            #                                      Chain filter
            #-----------------------------------------------------------------------------------------
            self.current_filter_chain = None
            # Creating the filter, feeding it with the liststore model
            self.chain_filter = self.residue_liststore.filter_new()
            # setting the filter function, note that we're not using the
            self.chain_filter.set_visible_func(self.chain_filter_func)
            #-----------------------------------------------------------------------------------------
            
            
            self.treeview = Gtk.TreeView(model = self.chain_filter)
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.MULTIPLE)
            
            self.treeview.connect("button-release-event", self.on_treeview_Objects_button_release_event)
            self.treeview.connect("row-activated", self.on_treeview_row_activated_event)
            
            for i, column_title in enumerate(
                ['', "index", "Residue",  "Chain", 'size']
            ):
                if i == 0:
                    #cell = Gtk.CellRendererToggle()
                    #cell.set_property('activatable', True)
                    #cell.connect('toggled', self.on_chk_renderer_toggled, self.residue_liststore)
                    
                    column = Gtk.TreeViewColumn(column_title, Gtk.CellRendererPixbuf(), pixbuf=0)
                    
                    #column = Gtk.TreeViewColumn(column_title, cell )
                    #column.add_attribute(cell, 'active', 0)
                    self.treeview.append_column(column)
                    #print ('aqui')
                else:
                    renderer = Gtk.CellRendererText()
                    #renderer.connect('toggled', self.on_chk_renderer_toggled, self.residue_liststore)
                    column = Gtk.TreeViewColumn(column_title, renderer, text=i)
                    self.treeview.append_column(column)

            
            self.current_filter_chain = None
            

            self.scrollable_treelist = Gtk.ScrolledWindow()
            self.scrollable_treelist.set_vexpand(True)
            self.scrollable_treelist.add(self.treeview)
            #------------------------------------------------------------------------------------------
            


            
            
            #------------------------------------------------------------------------------------------
            self.treeview_atom = Gtk.TreeView(model =self.atom_liststore)
            self.treeview_atom.connect("button-release-event", self.on_treeview_atom_button_release_event)
            self.treeview_atom.connect("row-activated", self.on_treeview_atom_row_activated_event)

            for i, column_title in enumerate(
                ['', "index", "name", "MM atom", 'MM charge']
            ):
                if i == 0:
                    #cell = Gtk.CellRendererToggle()
                    #cell.set_property('activatable', True)
                    #cell.connect('toggled', self.on_chk_renderer_toggled, self.atom_liststore)
                    #column = Gtk.TreeViewColumn(column_title, cell )
                    #column.add_attribute(cell, 'active', 0)
                    column = Gtk.TreeViewColumn(column_title, Gtk.CellRendererPixbuf(), pixbuf=0)
                    self.treeview_atom.append_column(column)
                    #print ('aqui')
                else:
                    renderer = Gtk.CellRendererText()
                    column = Gtk.TreeViewColumn(column_title, renderer, text=i)
                    self.treeview_atom.append_column(column)

            self.scrollable_treelist2 = Gtk.ScrolledWindow()
            self.scrollable_treelist2.set_vexpand(True)
            self.scrollable_treelist2.add(self.treeview_atom)
            #------------------------------------------------------------------------------------------
            
            
            
            
            
            
            
            
            
            
            
            self.box_vertical.pack_start(self.box_horizontal1, False, True, 0)
            self.box_vertical.pack_start(self.box_horizontal2, False, True, 0)
            self.treeviewbox_horizontal.pack_start(self.scrollable_treelist, True, True, 0)
            self.treeviewbox_horizontal.pack_start(self.scrollable_treelist2, True, True, 0)
            
            self.box_vertical.pack_start(self.treeviewbox_horizontal, False, True, 0)
            
            #self.refresh_system_liststore()
            #self.update_window (system_names = True, coordinates = True)
            
            self.combobox_systems.set_active(0)
            self.window =  Gtk.Window()
            self.window.set_border_width(10)
            self.window.set_default_size(600, 600)  
            self.window.add(self.box_vertical)
            self.window.connect("destroy", self.CloseWindow)
            self.window.set_title('Go to Atom Window') 
            self.window.show_all() 
                                                          
            #                                                                
            #self.builder.connect_signals(self)                                   
            
            self.visible  =  True
            #self.PutBackUpWindowData()
            #gtk.main()
            #----------------------------------------------------------------
        else:
            self.window.present()

    def CloseWindow (self, button):
        """ Function doc """
        #self.BackUpWindowData()
        self.window.destroy()
        self.visible    =  False
        #print('self.visible',self.visible)


    def getColouredPixmap( self, r, g, b, a=255 ):
        """ Given components, return a colour swatch pixmap """
        CHANNEL_BITS=8
        WIDTH=10
        HEIGHT=10
        swatch = GdkPixbuf.Pixbuf.new( GdkPixbuf.Colorspace.RGB, True, CHANNEL_BITS, WIDTH, HEIGHT ) 
        swatch.fill( (r<<24) | (g<<16) | (b<<8) | a ) # RGBA
        return swatch


    def refresh (self, option = 'all'):
        """ Function doc """
        self.update_window()
    
    def update_window (self, system_names = False, coordinates = True,  selections = False ):
        """ Function doc """

        if self.visible:
            
            _id = self.combobox_systems.get_active()
            if _id == -1:
                '''_id = -1 means no item inside the combobox'''
                return None
            else:    
                _, system_id = self.system_liststore[_id]
            
            if system_names:
                self.refresh_system_liststore ()
                self.combobox_systems.set_active(_id)
            
            if coordinates:
                self.refresh_coordinates_liststore ()
                
            
            #if selections:
            #    _, system_id = self.system_liststore[_id]
            #    self.refresh_selection_liststore(system_id)
        else:
            if system_names:
                self.refresh_system_liststore ()
            if coordinates:
                self.refresh_coordinates_liststore ()
            
    def refresh_coordinates_liststore(self, system_id = None):
        """ Function doc """
        cb_id =  self.coordinates_combobox.get_active()
        if system_id:
            pass
        else:
            _id = self.combobox_systems.get_active()
            if _id == -1:
                return False
            else:
                #print('_id', _id)
                _, system_id = self.system_liststore[_id]
        
        self.coordinates_liststore.clear()
        n = 0
        for key , vobject in self.vm_session.vm_objects_dic.items():
            if vobject.e_id == system_id:
                self.coordinates_liststore.append([vobject.name, key])
                n += 1
        
        self.coordinates_combobox.set_active(n-1)
        
    def refresh_system_liststore (self):
        """ Function doc """
        self.main.refresh_system_liststore()
        #self.system_liststore     .clear()
        ##self.selection_liststore  .clear()
        #'''--------------------------------------------------------------------------------------------'''
        ##self.combobox_systems =self.builder.get_object('systems_combobox')
        #for key, system  in self.p_session.systems.items():
        #    try:
        #        self.system_liststore.append([system['name'], key])
        #    except:
        #        print(system)

    def on_combobox_residues_changed (self, widget):
        """ Function doc """
        tree_iter = widget.get_active_iter()
        if tree_iter is not None:
            model = widget.get_model()
            residue = model[tree_iter][0]
            #print("Selected: country=%s" % country)
        
            self.current_filter_residue = residue
            
            #print("%s Chain selected!" % self.current_filter_residue)
            
            # we update the filter, which updates in turn the view
            if self.residue_filter:
                self.residue_filter.refilter()
        
        
    def on_combobox_chains_changed (self, widget):
        """ Function doc """
        ##---------------------------------------------------------------
        #self.current_filter_chain = None
        ## Creating the filter, feeding it with the liststore model
        #self.chain_filter = self.residue_liststore.filter_new()
        ## setting the filter function, note that we're not using the
        #self.chain_filter.set_visible_func(self.chain_filter_func)
        ##---------------------------------------------------------------
        
        tree_iter = self.combobox_chains.get_active_iter()
        if tree_iter is not None:
            model = self.combobox_chains.get_model()
            chain = model[tree_iter][0]
            #print("Selected: country=%s" % country)
        
        self.current_filter_chain = chain
        #print("%s Chain selected!" % self.current_filter_chain)
        # we update the filter, which updates in turn the view
        self.chain_filter.refilter()
    
    
    def on_combobox_systems_changed (self, widget):
        """ Function doc """
        #print(widget)
        #print(widget.get_active())
        
        #print(widget.get_active_id())
        #print(widget.get_active_iter())
        
        #self.vm_session.vm_objects_dic.items()
        
        #self.vm_session.vm_objects_dic.items()
        cb_id = widget.get_active()
        if cb_id == -1:
            return None
        else:
            
            self.update_window (coordinates = True)
            
            cb_id =  self.coordinates_combobox.get_active()
            _, key = self.coordinates_liststore[cb_id]
            
            self.VObj = self.vm_session.vm_objects_dic[key]
            
            
            self.liststore_chains = Gtk.ListStore(str)
            self.liststore_chains.append(['all'])
            chains = self.VObj.chains.keys()
            
            #self.chain_liststore = Gtk.ListStore(str)
            
            for chain in chains:
                self.liststore_chains.append([chain])
            self.combobox_chains.set_model(self.liststore_chains)
            self.combobox_chains.set_active(0)
            
            
            self.residue_liststore.clear() #= Gtk.ListStore(GdkPixbuf.Pixbuf, int, str, str, int)
            for chain in self.VObj.chains:
                for resi in self.VObj.chains[chain].residues:
                    #print(res.resi, res.resn, chain,  len(res.atoms) ) 
                    
                    #res_color = self.residues_dictionary[self.VObj.residues[resi].name]
                    color  =  self.VObj.color_palette['C']
                    res_color  = [int(color[0]*255),int(color[1]*255),int(color[2]*255)] 
                    swatch = self.getColouredPixmap( res_color[0], res_color[1], res_color[2] )
                    #self.colour_list.append( [ 'Sarcoline (literally mea…', 184, 134,  11, swatch ] )
                    #self.colour_list.append( [ 'Sarcoline (literally meaning: "flesh coloured", so varies...)', 184, 134,  11, swatch ] )

                    print(swatch, self.VObj.residues[resi].index, self.VObj.residues[resi].name , chain,  len(self.VObj.residues[resi].atoms)) 
                    self.residue_liststore.append(list([swatch, self.VObj.residues[resi].index, self.VObj.residues[resi].name , chain,  len(self.VObj.residues[resi].atoms)]))
            #-----------------------------------------------------------------------------------------
            #                                      Chain filter
            #-----------------------------------------------------------------------------------------
            self.current_filter_chain = None
            # Creating the filter, feeding it with the liststore model
            self.chain_filter = self.residue_liststore.filter_new()
            # setting the filter function, note that we're not using the
            self.chain_filter.set_visible_func(self.chain_filter_func)
            #-----------------------------------------------------------------------------------------
            
            
            
            
            
            
            #-----------------------------------------------------------------------------------------
            #                                      Residue combobox
            #-----------------------------------------------------------------------------------------
            self.liststore_residues = Gtk.ListStore(str)
            self.liststore_residues.append(['all'])
            
            resn_labels = {}
            
            for resi, residue in self.VObj.residues.items():
                resn_labels[residue.name] = True
            
            for resn in resn_labels.keys():
                #print (resn)
                self.liststore_residues.append([resn])
            
            self.combobox_residues.set_model(self.liststore_residues)
            self.combobox_residues.set_active(0)
            
            #-----------------------------------------------------------------------------------------
            #                                      Residue filter
            #-----------------------------------------------------------------------------------------
            self.current_filter_residue = None
            # Creating the filter, feeding it with the liststore model
            self.residue_filter = self.chain_filter.filter_new()
            # setting the filter function, note that we're not using the
            self.residue_filter.set_visible_func(self.residue_filter_func)
            #-----------------------------------------------------------------------------------------        
            
            
            #self.treeview.set_model(self.residue_liststore)
            #self.treeview.set_model(self.chain_filter)
            self.treeview.set_model(self.residue_filter)
            
        
    def on_treeview_atom_button_release_event(self, tree, event):
        if event.button == 2:
            selection     = tree.get_selection()
            model         = tree.get_model()
            (model, iter) = selection.get_selected()
            
            
            if iter != None:
                self.selectedID  = int(model.get_value(iter, 1))-1  # @+
                atom = self.VObj.atoms[self.selectedID]
                self.vm_session.vm_glcore.center_on_atom(atom)
       
    def on_treeview_atom_row_activated_event (self, tree, rowline , column):
        """ Function doc """
        selection   = tree.get_selection()
        model       = tree.get_model()
        data        = list(model[rowline])
        pickedID    = data[-1]
                
                
        cb_id =  self.coordinates_combobox.get_active()
        _, key = self.coordinates_liststore[cb_id]
        self.VObj = self.vm_session.vm_objects_dic[key]
                
        atom_picked = self.VObj.atoms[pickedID]
                    
        #atom_picked = self.vm_session.atom_dic_id[pickedID]

        #self.vm_session.selections[self.vm_session.current_selection].selection_function_viewing_set( selected= [atom_picked], _type= "atom")
        self.vm_session.selections[self.vm_session.current_selection].selecting_by_atom([atom_picked],  True)
        self.vm_session.selections[self.vm_session.current_selection]._build_selected_atoms_coords_and_selected_objects_from_selected_atoms()
        self.vm_session.vm_glcore.queue_draw()

    
    def on_treeview_row_activated_event(self, tree, rowline , column ):
        #print (A,B,C)
        selection     = tree.get_selection()
        model         = tree.get_model()
        
        #print(model)
        #print(rowline, list(model[rowline]))
        
        data  = list(model[rowline])
        self.selectedID  = int(data[1])  # @+
        self.selectedObj = str(data[2])
        self.selectedChn = str(data[3])
        
        cb_id =  self.coordinates_combobox.get_active()
        _, key = self.coordinates_liststore[cb_id]
        self.VObj = self.vm_session.vm_objects_dic[key]
        
        res = self.VObj.chains[self.selectedChn].residues[self.selectedID]
        
        
        '''centering and selecting'''
        frame = self.vm_session.get_frame ()
        res.get_center_of_mass(frame = frame)
        
        #self.vm_session.vm_glcore.center_on_coordinates(res.vm_object, res.mass_center)
        
        atom_keys = list(res.atoms.values())
        self.vm_session._selection_function_set({atom_keys[0]})
        self.vm_session.vm_glcore.updated_coords = True
        #self.vm_session.selections[self.vm_session.current_selection].selection_function_viewing_set( selected= {atom_keys[0]}, _type= "residue")
        
        #self.vm_session.selections[self.vm_session.current_selection].selecting_by_residue( selected_atoms = {atom_keys[0]} )
        #self.vm_session.selections[self.vm_session.current_selection]._build_selected_atoms_coords_and_selected_objects_from_selected_atoms()
        
        self.vm_session.vm_glcore.queue_draw()
        
        #self.atom_liststore.clear()
        
        #for atom in res.atoms.values():
        #    #res_color = self.residues_dictionary[self.VObj.residues[resi].name]
        #    color  =  self.VObj.color_palette[atom.symbol]
        #    swatch = self.getColouredPixmap( int(color[0]*255), int(color[1]*255), int(color[2]*255) )
        #    self.atom_liststore.append(list([swatch, int(atom.index), atom.name, atom.symbol ,  charges[atom.index-1] , int(atom.atom_id) ]))

        #self.treeview_atom.set_model(self.atom_liststore)
    def on_treeview_Objects_button_release_event(self, tree, event):
        #print ( tree, event)
        
        if event.button == 3:
            print ("button 3")
            #selection     = tree.get_selection()
            #model         = tree.get_model()
            #(model, iter) = selection.get_selected()
            #
            #
            #
            #if iter != None:
            #    self.selectedID  = int(model.get_value(iter, 1))  # @+
            #    self.selectedObj = str(model.get_value(iter, 2))
            #    print(self.selectedID, self.selectedObj, self.VObj.residues[self.selectedID])
            #    #self.builder.get_object('TreeViewObjLabel').set_label('- ' +self.selectedObj+' -' )
            #
            #    #widget = self.builder.get_object('treeview_menu')
            #    #widget.popup(None, None, None, None, event.button, event.time)
            #    #print ('button == 3')


        if event.button == 2:
            #print ('button == 2')
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.SINGLE)

            
            selection     = tree.get_selection()
            model         = tree.get_model()
            (model, iter) = selection.get_selected()

            if iter != None:
                self.selectedID  = int(model.get_value(iter, 1))  # @+
                self.selectedObj = str(model.get_value(iter, 2))
                self.selectedChn = str(model.get_value(iter, 3))
                res = self.VObj.chains[self.selectedChn].residues[self.selectedID]
                frame = self.vm_session.get_frame ()
                res.get_center_of_mass(frame = frame)
                
                
                #atomTypes       =      self.p_session.psystem[self.VObj.e_id].mmState.atomTypes
                charges         = list(self.p_session.psystem[self.VObj.e_id].mmState.charges)
                atomTypes       =      self.p_session.psystem[self.VObj.e_id].mmState.atomTypes
                atomTypeIndices = list(self.p_session.psystem[self.VObj.e_id].mmState.atomTypeIndices)
                
                self.vm_session.vm_glcore.center_on_coordinates(res.vm_object, res.mass_center)
        
                self.atom_liststore.clear()
                for atom in res.atoms.values():
                    #self.atom_liststore.append(list([True, int(atom.index), atom.name, atom.symbol, atom.charge, atom.atom_id ]))
                    #self.atom_liststore.append(list([True, int(atom.index), atom.name, atom.symbol ,  charges[atom.index-1] , int(atom.atom_id)]))
                    color  =  self.VObj.color_palette[atom.symbol]
                    swatch = self.getColouredPixmap( int(color[0]*255), int(color[1]*255), int(color[2]*255) )
                    self.atom_liststore.append(list([swatch, int(atom.index), atom.name, atom.symbol ,  charges[atom.index-1] , int(atom.atom_id) ]))
            
            
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.MULTIPLE)

        
        if event.button == 1:
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.SINGLE)

            selection     = tree.get_selection()
            model         = tree.get_model()
            (model, iter) = selection.get_selected()
            
            
            if iter != None:
                self.selectedID  = int(model.get_value(iter, 1))  # @+
                self.selectedObj = str(model.get_value(iter, 2))
                self.selectedChn = str(model.get_value(iter, 3))
                res = self.VObj.chains[self.selectedChn].residues[self.selectedID]
                
                
                self.atom_liststore.clear()
                
                '''when we don't have an MM (molecular mechanics) system, the "system" object doesn't 
                have the mmState attribute, so we need to get the information of atom names and indexes
                 in another way.'''
                #--------------------------------------------------------------------------------------
                try:
                    charges         = list(self.p_session.psystem[self.VObj.e_id].mmState.charges)
                    atomTypes       =      self.p_session.psystem[self.VObj.e_id].mmState.atomTypes
                    atomTypeIndices = list(self.p_session.psystem[self.VObj.e_id].mmState.atomTypeIndices)
                except:
                    charges         = [0.0]*len(self.p_session.psystem[self.VObj.e_id].atoms)
                    atomTypes       = []
                    atomTypeIndices = []
                    for atom in self.p_session.psystem[self.VObj.e_id].atoms.items:
                        atomTypes.append('-')
                        atomTypeIndices.append(atom.index)
                #--------------------------------------------------------------------------------------
                for atom in res.atoms.values():
                     
                    color  =  self.VObj.color_palette[atom.symbol]
                    #print (color)
                    swatch = self.getColouredPixmap( int(color[0]*255), int(color[1]*255), int(color[2]*255) )
                    self.atom_liststore.append(list([swatch, int(atom.index), str(atom.name), str(atomTypes[atomTypeIndices[atom.index-1]]) ,  charges[atom.index-1] , int(atom.atom_id)]))
                    #self.atom_liststore.append(list([True, int(atom.index), atom.name, atomTypes[atomTypeIndices[atom.index-1]] ,  charges[atom.index-1] , int(atom.atom_id)]))

                
                
                #self.treeview_atom.set_model(self.atom_liststore)
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.MULTIPLE)

    
    def on_chk_renderer_toggled(self, cell, path, model):
        print(model[path][0])


    def residue_filter_func(self, model, iter, data):
        """Tests if the language in the row is the one in the filter"""
        if (
            self.current_filter_residue is None
            or self.current_filter_residue == "all"
        ):
            return True
        else:
            return model[iter][2] == self.current_filter_residue
   
            
    def chain_filter_func(self, model, iter, data):
        """Tests if the language in the row is the one in the filter"""
        if (
            self.current_filter_chain is None
            or self.current_filter_chain == "all"
        ):
            return True
        else:
            return model[iter][3] == self.current_filter_chain

    #def on_selection_button_clicked(self, widget):
    #    """Called on any of the button clicks"""
    #    # we set the current language filter to the button's label
    #    self.current_filter_language = widget.get_label()
    #    print("%s language selected!" % self.current_filter_language)
    #    # we update the filter, which updates in turn the view
    #    self.language_filter.refilter()
    #
    #def on_button1_clicked(self, widget):
    #    print("Hello")
    #
    #def on_button2_clicked(self, widget):
    #    print("Goodbye")
    def update (self):
        """ Function doc """
        #print('VismolGoToAtomWindow2 update')
        pass
        #self.self.combobox_systems.set_active(-1)

class EasyHybridDialogSetQCAtoms(Gtk.Dialog):
    def __init__(self, parent):
        super().__init__(title="New QC list", transient_for=parent, flags=0)
        self.add_buttons(
            Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_YES, Gtk.ResponseType.YES
        )

        self.set_default_size(150, 100)

        label = Gtk.Label(label="A new quantum region has been defined. Would you like to set up your QC parameters now?")

        box = self.get_content_area()
        box.add(label)
        self.show_all()

    def update (self):
        """ Function doc """
        #print('VismolGoToAtomWindow2 update')
        pass


class EasyHybridDialogPrune:
    """ Class doc """

    def __init__ (self, home = None, num_of_atoms =  0,  name = 'Unknow', tag = 'UNK'):
        """ Class initialiser """
        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(home,'gui/windows/windows_and_dialogs.glade'))
        #self.builder.connect_signals(self)
        
        self.dialog       = self.builder.get_object('dialog_prune')

        self.builder.get_object('entry_number_of_atoms').set_text(str(num_of_atoms))
        self.builder.get_object('entry_name').set_text(name + '_pruned')
        self.builder.get_object('entry_tag').set_text(tag)
        
        self.builder.get_object('prune_dialog_button_prune').connect("clicked", self.on_click_button_prune)
        self.builder.get_object('prune_dialog_button_cancel').connect("clicked", self.on_click_button_cancel)
        
        self.prune = False
        self.name  = None
        self.tag   = None
        answer = self.dialog.run()
        #print ('answer', answer)
    
    def on_click_button_prune (self, widget):
        """ Function doc """
        num_of_atoms = self.builder.get_object('entry_number_of_atoms').get_text( )
        self.name    = self.builder.get_object('entry_name').get_text( )
        self.tag     = self.builder.get_object('entry_tag').get_text( )
        self.prune   = True

        self.dialog.destroy()

    def on_click_button_cancel (self, widget):
        """ Function doc """
        self.dialog.destroy()
        self.prune   = False

    def update (self):
        """ Function doc """
        #print('VismolGoToAtomWindow2 update')
        pass


class ImportANewSystemWindow(Gtk.Window):
    """ Class doc """
    def __init__(self, main = None, home = None):
        """ Class initialiser """
        
        self.easyhybrid_main     = main
        self.home                = main.home
        self.Visible             = False        
        
        self.residue_liststore = Gtk.ListStore(str, str, str)
        #self.atom_liststore    = Gtk.ListStore(bool, int, str, str, int, int)
        #self.residue_filter    = False
        
        self.charmm_txt = '''When using the traditional CHARMM/PSF format for system preparation, it is necessary to have three types of files: parameter files (in formats such as prm or par), topology files in the form of psf, and coordinate files (which can be in various formats including chm, crd, pdb, xyz, etc. It's worth noting that multiple parameter files may be required depending on the system being simulated.'''
        self.amber_txt = '''In order to properly run simulations utilizing the AMBER force field, it is necessary to have two types of files: topologies (in the form of either top or prmtop files) and coordinates (which can be in the form of crd, pdb, xyz, or other similar file types).'''
        self.OPLS_txt = '''When preparing systems natively in pDynamo utilizing the OPLS force field, it is necessary to have two components: a folder containing the OPLS parameters, and a topology and coordinate file (in formats such as pdb, mol2, or mol).'''
        self.DYFF_txt = '''The DYFF force field is a generic force field specifically designed for use within the pDynamo program. When utilizing DYFF to prepare systems natively in pDynamo, it is essential to have two components: a folder containing the necessary force field parameters, and a coordinate file in formats such as pdb, mol2, or mol.'''
        self.gmx_txt = '''Systems prepared natively in GROMACS using the CHARMM(or AMBER) force field require: A parameter/topology (top) file and coordinate file (pdb, mol2, mol, ...) .  '''
        
    
    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'gui/windows/import_system_window.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('ImportNewSystemWindow')
            self.window.set_border_width(10)
            self.window.set_default_size(500, 370)  

            
            '''--------------------------------------------------------------------------------------------'''
            self.system_type_store = Gtk.ListStore(str, int)
            system_types = [
                ["AMBER"                        , 0],
                ["CHARMM"                       , 1],
                ["OPLS"                         , 2],
                ['DYFF'                         , 5], 
                ["pdynamo files (*.pkl, *.yaml)", 3],
                ["other (*.pdb, *.xyz, *.mol2)" , 4],
                ]
            for system_type in system_types:
                self.system_type_store.append(system_type)
                #print (system_type)
            
            self.system_types_combo = Gtk.ComboBox.new_with_model(self.system_type_store)
            #self.system_types_combo = self.builder.get_object('system_type_combox_from_import_a_new_system')
            self.box_combo = self.builder.get_object('box')
            self.box_combo.pack_start(self.system_types_combo, True, True, 0)
            
            self.system_types_combo.connect("changed", self.on_name_combo_changed)
            self.system_types_combo.set_model(self.system_type_store)
            
            renderer_text = Gtk.CellRendererText()
            self.system_types_combo.pack_start(renderer_text, True)
            self.system_types_combo.add_attribute(renderer_text, "text", 0)
            '''--------------------------------------------------------------------------------------------'''
           
            
            self.treeview = self.builder.get_object('gtktreeview_import_system')
            for i, column_title in enumerate(['file', "type", "number of atoms"]):
                renderer = Gtk.CellRendererText()
                column = Gtk.TreeViewColumn(column_title, renderer, text=i)
                self.treeview.append_column(column)

            
            self.window.show_all()                                               
            self.builder.connect_signals(self)                                   
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()

            self.Visible  =  True
            
            self.files    = {
                            'amber_prmtop': None,
                            'charmm_par'  : [],
                            'charmm_psf'  : None,
                            'charmm_extra': None, 
                            'prm_folder' : [],
                            'coordinates' : None,
                            }
            self.system_types_combo.set_active(0)

            self.window.connect('destroy', self.CloseWindow)
            self.builder.get_object('button_load_files')        .connect('clicked', self.on_button_load_files_clicked)
            self.builder.get_object('button_remove_files')      .connect('clicked', self.on_button_delete_files_clicked)
            #self.builder.get_object('gtkbox_OPLS_folderchooser').connect('clicked',)
            self.builder.get_object('button_cancel')            .connect('clicked', self.CloseWindow)
            self.builder.get_object('import_import_system')     .connect('clicked', self.on_button_import_system_clicked)

        else:
            self.window.present()
            #----------------------------------------------------------------
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        #self.BackUpWindowData()
        self.window.destroy()
        self.Visible    =  False
        #print('self.Visible',self.Visible)
    
    def on_name_combo_changed(self, widget):
        """ Function doc """
        #fftype = self.system_types_combo.get_active()
        
        active_iter = self.system_types_combo.get_active_iter()
        if active_iter:
            model = self.system_types_combo.get_model()
            fftype = model[active_iter][1]
            #print(model[active_iter][0])
        
        #print (fftype)
        self.files    = {
            'amber_prmtop': None,
            'charmm_par'  : [],
            'charmm_psf'  : None,
            'charmm_extra': None, 
            'prm_folder' : [],
            'coordinates' : None,
            'charges'     : None,
            }
        self.residue_liststore = Gtk.ListStore(str, str, str)
        self.treeview.set_model(self.residue_liststore)
            
        if fftype == 0: #AMBER
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()
            self.builder.get_object('gtk_label_fftype').set_text(self.amber_txt)
 
            
        if fftype == 1: # "CHARMM":
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()
            self.builder.get_object('gtk_label_fftype').set_text(self.charmm_txt)

            
        if fftype == 10: #"GROMACS":
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()
            self.builder.get_object('gtk_label_fftype').set_text(self.gmx_txt)

            
        if fftype == 2:#"OPLS":
            self.builder.get_object('gtkbox_OPLS_folderchooser').show()
            self.builder.get_object('gtk_label_fftype').set_text(self.OPLS_txt)
            
            '''Eventually, the user may have another set of parameters, but 
            by default, the pDynamo3 parameter directories are first searched.'''
            path = os.environ.get('PDYNAMO3_PARAMETERS')
            path = os.path.join(path,'forceFields/opls/protein')
            self.builder.get_object('OPLS_folderchooserbutton').set_filename(path)
            
            
        if fftype == 3: #"pDynamo files(*.pkl,*.yaml)":
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()
            
        if fftype == 4: #"Other(*.pdb,*.xyz,*.mol2...)":
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()

        if fftype == 5:#"DYFF":
            self.builder.get_object('gtkbox_OPLS_folderchooser').show()
            self.builder.get_object('gtk_label_fftype').set_text(self.DYFF_txt)
            #try:
            '''Eventually, the user may have another set of parameters, but 
            by default, the pDynamo3 parameter directories are first searched.'''
            path = os.environ.get('PDYNAMO3_PARAMETERS')            
            path = os.path.join(path,'forceFields/dyff/dyff-1.0')
            self.builder.get_object('OPLS_folderchooserbutton').set_filename(path)
            #except:
            #    pass
                
    def filetype_parser(self, filein, systemtype):
        filetype = get_file_type(filein)
        #print (filetype, systemtype)
        if filetype in ['top', 'prmtop', 'TOP', 'PRMTOP', ]:
            if systemtype == 0:
                self.files['amber_prmtop'] = filein
                return 'amber parameters/topologies'
                
   
        
        elif filetype in ['par', 'prmtop', 'prm', 'PAR', 'PRM', 'str']:
            if systemtype == 1:
                self.files['charmm_par'].append(filein)
                return 'charmm parameters'
        
        
        
        
        elif filetype in ['psf', 'PSF','psfx', 'PSFX']:
            if systemtype == 1:
                self.files['charmm_psf'] = filein
                return 'charmm topologies'
        
        
        
        elif filetype in ['pdb', 'PDB','mol','MOL','mol2','MOL2', 'xyz', 'XYZ', 'crd', 'inpcrd', 'chm']:
            #if systemtype == 1:
            self.files['coordinates'] = filein
            if filetype == 'mol2':
                atoms, bonds = read_MOL2(filein)
                self.files['charges'] = atoms['charges']
                pass
            return 'coordinates'
        
        
        
        elif filetype in ['pkl', 'PKL']:
            #if systemtype == 1:
            self.files['coordinates'] = filein
            return 'pDynamo coordinates'
        else:
            return 'unknow'
        
    #def GetFileType(self, filename):
    #    file_type = filename.split('.')
    #    return file_type[-1]

    def on_button_delete_files_clicked (self, button):
        """ Function doc """
        self.residue_liststore = Gtk.ListStore(str, str, str)
        self.treeview.set_model(self.residue_liststore)
    
    def on_button_load_files_clicked (self, button):
        """ Function doc """
        files = self.easyhybrid_main.filechooser.open(select_multiple = True)
        #print(files)

        for _file in files:
            #for res in self.VObj.chains[chain].residues:
                ##print(res.resi, res.resn, chain,  len(res.atoms) ) 
            
            active_iter = self.system_types_combo.get_active_iter()
            if active_iter:
                model = self.system_types_combo.get_model()
                #print(model[active_iter][1]) # access the value of the active item
                systemtype = model[active_iter][1]
            
            filetype = self.filetype_parser( _file, systemtype)
            self.residue_liststore.append(list([_file, filetype, 'unk' ]))
        self.treeview.set_model(self.residue_liststore)
        self.files['prm_folder'] =  self.builder.get_object('OPLS_folderchooserbutton').get_filename()
        #print(self.files)

    def on_button_import_a_new_system_clicked (self, button):
        """ Function doc """
        
        if button == self.builder.get_object('ok_button_import_a_new_system'):
            print('ok_button_import_a_new_system')
            #self.on_button1_clicked_create_new_project(button)
            #self.dialog.hide()
        if button == self.builder.get_object('cancel_button_import_a_new_system'):
            print('cancel_button_import_a_new_system')
            self.dialog.hide()
            
    def on_button_import_system_clicked (self, button):
        #print('ok_button_import_a_new_system')
        #systemtype = self.system_types_combo.get_active()
        active_iter = self.system_types_combo.get_active_iter()
        if active_iter:
            model = self.system_types_combo.get_model()
            #print(model[active_iter][1]) # access the value of the active item
            systemtype = model[active_iter][1]
        else:
            return None
        self.files['prm_folder'] =  self.builder.get_object('OPLS_folderchooserbutton').get_filename()
        name =  self.builder.get_object('entry_system_name').get_text()
        tag  =  self.builder.get_object('entry_system_tag').get_text()
        print(self.files, systemtype)
        
        #'''
        self.easyhybrid_main.p_session.load_a_new_pDynamo_system_from_dict(input_files = self.files, 
                                                                           system_type = systemtype, 
                                                                           name        = name      ,
                                                                           tag         = tag
                                                                           )
        #'''
        #if systemtype == 2:
        #    self.files['prm_folder'] =  self.builder.get_object('OPLS_folderchooserbutton').get_filename()
        
        #print ('systemtype',systemtype, self.files )
        ##self.easyhybrid_main.p_session.get_bonds_from_pDynamo_system()
        #
        #vobject = self.easyhybrid_main.p_session.build_vobject_from_pDynamo_system (name = name)
        self.CloseWindow(button, data  = None)

    def update (self):
        """ Function doc """
        #print('VismolGoToAtomWindow2 update')
        pass






class ImportTrajectoryWindow:
    """ Class doc """
    def __init__(self, main = None):
        """ Class initialiser """
        self.main                = main
        self.p_session           = main.p_session
        self.home                = main.home
        
        self.Visible             =  False        
        self.starting_coords_liststore = Gtk.ListStore(str, int)
        
        self.data_type_dict = {
                        0:'pklfile'    , #  single file
                        1:'pklfolder'  , #  trajectory
                        2:'pklfolder2D', #  2d trajectory  
                        3:'pdbfile'    ,
                        4:'pdbfolder'  ,
                        5:'dcd',
                        6:'crd',
                        7:'xyz',
                        8:'mol2',
                        9:'netcdf',
                       10:'log_file'  }
        
        self.folder_type_list = ['pklfolder', 'pklfolder2D', 'pdbfolder']

    def OpenWindow (self, sys_selected = None):
        """ Function doc """
        if self.Visible  ==  False:
            '''--------------------------------------------------------------------------------------------'''
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'gui/windows/import_trajectory_window.glade'))
            self.builder.connect_signals(self)
            self.window = self.builder.get_object('import_trajectory_window')
            self.window.set_title('Import Trajectory Window')
            self.window.set_keep_above(True)
            '''--------------------------------------------------------------------------------------------'''

            #------------------------------------------------------------------------------------
            self.combobox_pdynamo_system = self.builder.get_object('combobox_pdynamo_system')
            renderer_text = Gtk.CellRendererText()
            self.combobox_pdynamo_system.add_attribute(renderer_text, "text", 0)
            

            
            #------------------------------------------------------------------------------------
            """This block of instructions is used to organize the system combobox. 
            The menu function provides the system ID, but the position in the combobox 
            list must be calculated dynamically, as the index of items will not match 
            the system ID when a system has been deleted."""
            #self._define_system_liststore()
            #self.system_liststore = Gtk.ListStore(str, int)
            #names      = [ ]
            
            counter    = 0
            set_active = 0
            for key , system in self.p_session.psystem.items():
                if system:
                    #name = system.label
                    #self.system_liststore.append([name, int(key)])
                    if sys_selected == int(key):
                        set_active = counter
                        counter += 1
                    else:
                        counter += 1
                else:
                    pass
            #------------------------------------------------------------------------------------
            
            
            #self.combobox_pdynamo_system.set_model(self.system_liststore)
            self.combobox_pdynamo_system.set_model(self.main.system_liststore)

            if sys_selected:
                self.combobox_pdynamo_system.set_active(set_active)
            else:
                self.combobox_pdynamo_system.set_active(set_active)
            #------------------------------------------------------------------------------------
            
            
            
            
            
            #'''--------------------------------------------------------------------------------------------------
            self.folder_chooser_button = FolderChooserButton(main =  self.main)
            self.builder.get_object('folder_chooser_box').pack_start(self.folder_chooser_button.btn, True, True, 0)
            #'''--------------------------------------------------------------------------------------------------
            
            
            
            #------------------------------------------------------------------------------------
            self.combobox_starting_coordinates = self.builder.get_object('vobjects_combobox')
            #self.combobox_starting_coordinates.set_model(self.starting_coords_liststore)
            
            renderer_text = Gtk.CellRendererText()
            self.combobox_starting_coordinates.pack_start(renderer_text, True)
            self.combobox_starting_coordinates.add_attribute(renderer_text, "text", 0)
            
            size = len(self.starting_coords_liststore)
            self.combobox_starting_coordinates.set_active(size-1)
            #------------------------------------------------------------------------------------
            
            
            #------------------------------------------------------------------------------------
            self.builder.get_object('vobjects_combobox').set_sensitive(False)
            #------------------------------------------------------------------------------------

            #'''--------------------------------------------------------------------------------------------'''
            self.combox = self.builder.get_object('combobox_coordinate_type')
            self.combox.connect("changed", self.on_combobox_coordinate_type)
    
            #self.folder_chooser_button.btn.connect('clicked', self.print_test)
            self.on_combobox_pdynamo_system(None)
            self.combox.set_active(0)
            self.window.show_all()
            self.Visible  = True
        
        else:
            self.window.present()
    
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False
    
    def print_test (self, widget):
        """ Function doc """
        #print(self.folder_chooser_button.folder)
        trajfolder  = self.folder_chooser_button.folder
        basename  = os.path.basename(trajfolder)
        logfile   = basename[:-5]+'log'
        
        logfile   = os.path.join(trajfolder, logfile)
        
        if exists(logfile):
            self.builder.get_object('file_chooser_btn_logfile').set_filename(logfile)
        else:
            pass
        
    def on_combobox_pdynamo_system (self, widget):
        """ Function doc """
        tree_iter = self.combobox_pdynamo_system.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.combobox_pdynamo_system.get_model()
            name, sys_id = model[tree_iter][:2]
            print (name, sys_id)
        else:
            return False

        self.combobox_starting_coordinates = self.builder.get_object('vobjects_combobox')
        self.combobox_starting_coordinates.set_model(self.main.vobject_liststore_dict[sys_id])

    def on_combobox_coordinate_type (self, widget):
        """ Function doc """
        data_type = self.builder.get_object('combobox_coordinate_type').get_active() 

        
        
        if data_type == 10:
            self.builder.get_object('frame_stride_box').set_sensitive(False) 
            self.builder.get_object('folder_chooser_box').set_sensitive(False) 
            self.builder.get_object('entry_create_a_new_vobj').set_sensitive(False)
            self.builder.get_object('vobjects_combobox').set_sensitive(True)

        
        elif data_type == 0:
            self.builder.get_object('frame_stride_box').set_sensitive(False) 
            self.builder.get_object('folder_chooser_box').set_sensitive(True) 
            self.builder.get_object('entry_create_a_new_vobj').set_sensitive(True)
            self.builder.get_object('vobjects_combobox').set_sensitive(True)

        elif data_type in [1,3,5,6,7,8,9]:
            self.builder.get_object('frame_stride_box').set_sensitive(True) 
            self.builder.get_object('folder_chooser_box').set_sensitive(True) 
            self.builder.get_object('entry_create_a_new_vobj').set_sensitive(True)
            self.builder.get_object('vobjects_combobox').set_sensitive(True)

        elif data_type in [2,4]:
            self.builder.get_object('frame_stride_box').set_sensitive(False) 
            self.builder.get_object('folder_chooser_box').set_sensitive(True) 
            self.builder.get_object('entry_create_a_new_vobj').set_sensitive(True)
            self.builder.get_object('vobjects_combobox').set_sensitive(False)
        else:
            self.builder.get_object('frame_stride_box').set_sensitive(True) 
            self.builder.get_object('folder_chooser_box').set_sensitive(True) 
            self.builder.get_object('entry_create_a_new_vobj').set_sensitive(True)
            self.builder.get_object('vobjects_combobox').set_sensitive(True)        
        
        print (data_type, self.data_type_dict[data_type])
        data_type = self.data_type_dict[data_type]
        
        if  data_type in self.folder_type_list:
            self.folder_chooser_button.sel_type = 'folder'
        else:
            self.folder_chooser_button.sel_type = 'file'        

    def on_vobject_combo_changed (self, widget):
        '''this combobox has the reference to the starting coordinates of a simulation'''
        #combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]
            print (name, model[tree_iter][:2])
            #name, vobject_id = model[tree_iter][:2]
    def on_name_combo_changed (self, widget):
        """ Function doc """
        traj_type = self.builder.get_object('combobox_coordinate_type').get_active() 
        print (traj_type, self.data_type_dict[traj_type])
        traj_type = self.data_type_dict[traj_type]
        
        if  traj_type in self.folder_type_list:
            self.folder_chooser_button.sel_type = 'folder'
        else:
            self.folder_chooser_button.sel_type = 'file'

    def on_radio_button_change (self, widget):
        """ Function doc """
        if self.builder.get_object('radiobutton_import_as_new_object').get_active():
            self.builder.get_object('vobjects_combobox').set_sensitive(False)
            self.builder.get_object('entry_create_a_new_vobj').set_sensitive(True)

        if self.builder.get_object('radiobutton_append_to_a_vobject').get_active():
            self.builder.get_object('entry_create_a_new_vobj').set_sensitive(False)
            self.builder.get_object('vobjects_combobox').set_sensitive(True)

    def import_data (self, button):
        """ Function doc """
        
        parameters = {
                      'system_id'    : None,
                                     
                      'data_path'    : None,
                      'data_type'    : 0   ,
                      
                      'new_vobj_name': None, 
                      'vobject_id'   : None,
                      'vobject'      : None,

                      'logfile'      : None,
                      
                      'first'        : None,
                      'last'         : None,
                      'stride'       : None,
                     }
        
        #entry_name     = None
        # - - - system - - - 
        tree_iter = self.combobox_pdynamo_system.get_active()
        model = self.combobox_pdynamo_system.get_model()
        name, sys_id = model[tree_iter][:2]
        parameters['system_id'] = sys_id
        
        # data, a file or a folder - many files
        parameters['data_path'] = self.folder_chooser_button.get_folder()
        
        # A log file - optional
        parameters['logfile']   = self.builder.get_object('file_chooser_btn_logfile').get_filename()
        
        
        data_type = self.builder.get_object('combobox_coordinate_type').get_active()        
        parameters['data_type'] = self.data_type_dict[data_type]
        #-----------------------------------------------------------------------------
        #tree_iter = self.combobox_pdynamo_system.get_active_iter()
        #if tree_iter is not None:
        #    '''selecting the vismol object from the content that is in the combobox '''
        #    model = self.combobox_pdynamo_system.get_model()
        #    _name, system_id = model[tree_iter][:2]
        #
        #else:
        #    print ('Error: Please select a pDynamo system before continuing. ')
        #-----------------------------------------------------------------------------
        
        
        
        #----------------------------------------------------------------------------------------------
        '''Here, we determine whether it is necessary to create a new object 
        or append the coordinates to an existing one.'''
        if self.builder.get_object('radiobutton_import_as_new_object').get_active():
            parameters['new_vobj_name'] = self.builder.get_object('entry_create_a_new_vobj').get_text()
        else:
            #-----------------------------------------------------------------------------
            tree_iter = self.combobox_starting_coordinates.get_active_iter()
            if tree_iter is not None:
                
                '''selecting the vismol object from the content that is in the combobox '''
                model = self.combobox_starting_coordinates.get_model()
                name, vobject_id = model[tree_iter][:2]
                vobject = self.main.vm_session.vm_objects_dic[vobject_id]
                
                parameters['vobject_id'] = vobject_id
                parameters['vobject']    = vobject
            #-------------------------------------------------------------------------------------------

        print(parameters)
        self.main.p_session.import_data ( parameters ) 

    def update (self):
        """ Function doc """
        #print('VismolGoToAtomWindow2 update')
        pass



class TrajectoryPlayerWindow:
    """ Class doc """
    
    def __init__ (self, main = None):
        """ Class initialiser """
        self.main = main
        self.vm_session = main.vm_session
        self.p_session  = main.p_session
            
        self.Visible = False
    
    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.window = Gtk.Window()
            self.window.connect('destroy-event', self.CloseWindow)
            self.window.connect('delete-event', self.CloseWindow)
            self.vm_traj_obj = VismolTrajectoryFrame(self.vm_session)
            
            self.window.add(self.vm_traj_obj.get_box())
            
            self.window.set_default_size(300, 50)  
            self.window.set_title('EasyHybrid Player')  
            self.window.set_keep_above(True)
            
            self.window.show_all()
            self.Visible  = True
    
        else:
            self.window.present()
        

    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.vm_traj_obj.stop(None)

        self.window.destroy()
        self.Visible    =  False
        self.main.trajectory_player_button.set_active(False)


    def change_range (self, upper = 100):
        """ Function doc """
        if self.Visible:
            print('upper =', upper)
            self.vm_traj_obj.change_range (upper = upper)
        








from matplotlib.backends.backend_gtk3agg import FigureCanvas  # or gtk3cairo.
from matplotlib.figure import Figure
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import copy as cp
import gc
import os


class PotentialEnergyAnalysisWindow():

    def __init__(self, main = None ):
        """ Class initialiser """
        self.main                = main
        self.p_session           = self.main.p_session
        self.vm_session          = self.main.vm_session
        self.Visible             =  False        
        
        self.vobject_liststore   = Gtk.ListStore(str, int)
        self.data_liststore      = Gtk.ListStore(str, int)
        
        #self.interpolate = False#True
        self.interpolate =  True
        self.opt_methods = { 0 : 'ConjugatedGradient',
                             1 : 'SteepestDescent'   ,
                             2 : 'LFBGS'             ,
                             3 : 'QuasiNewton'       ,
                             4 : 'FIRE'              }
        
        self.xdata = []
        self.ydata = []
        self.zdata = []
        self.xy_traj = []
        self.pcm = None
        self.color_bar = None
        self.vobject = None
        self.traj_export_index = 1

    def OpenWindow (self, vobject = None):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()

            self.builder.add_from_file(os.path.join(self.main.home,'gui/windows/PES_analysis_window.glade'))
            self.builder.connect_signals(self)
            self.vobject = vobject
            self.window = self.builder.get_object('window')
            self.window.set_title('Analysis Window')
            self.window.set_keep_above(True)            
            self.window.set_default_size(700, 450)
            
            self.hbox           = self.builder.get_object('hbox_matplot_figures')
            self.scale_traj     = self.builder.get_object('scale_trajectory_from_PES')
            self.adjustment     = Gtk.Adjustment(value         = 0,
                                                 lower         = 0,
                                                 upper         = 100,
                                                 step_increment= 1,
                                                 page_increment= 1,
                                                 page_size     = 1)

            self.scale_traj.set_adjustment ( self.adjustment)
            self.scale_traj.set_digits(0)
            self.reaction_coord_label = self.builder.get_object('RC1_RC2_label')
            
            
            
            
            
            '''-------------------------------------------------------------'''
            self.fig = Figure(figsize=(6, 4))#,constrained_layout=True)
            self.canvas = FigureCanvas(self.fig)  # a Gtk.DrawingArea
            self.canvas.mpl_connect('button_press_event', self.on_pick)
            self.canvas.mpl_connect('motion_notify_event', self.on_motion_notify_event)

            self.hbox.pack_start(self.canvas, True, True, 0)
            self.ax = self.fig.add_subplot(1,1,1)
            '''-------------------------------------------------------------'''




            '''-------------------------------------------------------------'''
            self.fig2 = Figure(figsize=(6, 4), dpi=100)
            self.ax2 = self.fig2.add_subplot(1,1,1)
            self.ax3 = self.fig2.add_subplot(1,1,1)
            
            #self.line2, = self.ax2.plot([], [], '-o')
            #self.secax = self.ax2.secondary_xaxis('top', functions=(self.seconday_xaxis_function))

            self.canvas2 = FigureCanvas(self.fig2)  # a Gtk.DrawingArea
            self.hbox.pack_start(self.canvas2, True, True, 0)        
            '''-------------------------------------------------------------'''





            
            self.grid = self.builder.get_object('grid_setup')
            
            sys_selected = 0
            
            
            
            #'''------------------------------------------------------------------------------------
            self.vobject_liststore = Gtk.ListStore(str, int, int) # name, vobj_id, sys_id
            names = [ ]
            
            for system_id , system in self.p_session.psystem.items():
                for vobject_id in system.e_logfile_data.keys():
                    _vobject = self.vm_session.vm_objects_dic[vobject_id]
                    self.vobject_liststore.append([_vobject.name, vobject_id, system_id])
            
            
            #for key , system in self.p_session.psystem.items():
            #    #print ()
            #    try:
            #        print("system data", system.e_logfile_data)
            #        if system.e_logfile_data:
            #            for vobject_id in  system.e_logfile_data.keys():
            #                print(['vobject_id:', vobject_id ,systeme_logfile_data.keys()])
            #                
            #                try:
            #                    _vobject = self.main.vm_session.vm_objects_dic[vobject_id]
            #                    print(['_vobject:', _vobject.name,_vobject.index])
            #                    self.vobject_liststore.append([_vobject.name, _vobject.index])
            #                except:
            #                    print('self.vobject_liststore.append([_vobject.name, _vobject.index])', 'failed!')
            #    except:
            #        pass
            #self.vobject_liststore.append(['all', _vobject.index])

            self.coordinates_combobox = Gtk.ComboBox()
            renderer_text = Gtk.CellRendererText()
            self.coordinates_combobox.pack_start(renderer_text, True)
            self.coordinates_combobox.add_attribute(renderer_text, "text", 0)
            
            self.coordinates_combobox.set_model(self.vobject_liststore)
            
            self.coordinates_combobox.connect('changed', self.on_coordinates_combobox_change)
            self.grid.attach (self.coordinates_combobox, 1, 0, 1, 1)
            #------------------------------------------------------------------------------------
            
            
            
            
            
            #------------------------------------------------------------------------------------
            self.data_combobox = Gtk.ComboBox()
            renderer_text = Gtk.CellRendererText()
            self.data_combobox.pack_start(renderer_text, True)
            self.data_combobox.add_attribute(renderer_text, "text", 0)
            self.data_combobox.set_model(self.data_liststore)
            self.data_combobox.connect('changed', self.on_data_combobox_change)
            self.grid.attach (self.data_combobox, 3, 0, 1, 1)
            #------------------------------------------------------------------------------------
            
            
            
            #------------------------------------------------------------------------------------
            if sys_selected:
                self.coordinates_combobox.set_active(sys_selected)
            else:
                self.coordinates_combobox.set_active(0)
            #------------------------------------------------------------------------------------'''
 
            self.window.show_all()
            self.Visible  = True
    
        else:
            
            _id = self.coordinates_combobox.get_active()
            
            self.vobject_liststore = Gtk.ListStore(str, int, int) # name, vobj_id, sys_id
            
            for system_id , system in self.p_session.psystem.items():
                for vobject_id in system.e_logfile_data.keys():
                    _vobject = self.vm_session.vm_objects_dic[vobject_id]
                    self.vobject_liststore.append([_vobject.name, vobject_id, system_id])
            
            
            #self.vobject_liststore = Gtk.ListStore(str, int)
            #names = [ ]
            #for key , system in self.p_session.psystem.items():
            #    try:
            #        for vobject_id in  system.e_logfile_data.keys():
            #            try:
            #                _vobject = self.main.vm_session.vm_objects_dic[vobject_id]
            #                print(['_vobject:', _vobject.name,_vobject.index])
            #                self.vobject_liststore.append([_vobject.name, _vobject.index])
            #            except:
            #                print('self.vobject_liststore.append([_vobject.name, _vobject.index])', 'failed!')
            #    except:
            #        pass
            #self.vobject_liststore.append(['all', _vobject.index])
            self.coordinates_combobox.set_model(self.vobject_liststore)
            self.coordinates_combobox.set_active(_id)

    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False

    def seconday_xaxis_function (self):
        """ Function doc """
        pass

    def _draw_data (self, cla = True, refresh = True):
        """ Function doc """
        print('self.pcm.cla()')
        if cla:
            if self.ax:
                self.ax.cla()
                #if self.pcm:
                self.fig.gca().clear()
                self.fig.canvas.draw()
            else:
                pass

        self.pcm = self.ax.pcolormesh(range(len(self.data['Z'][0])), 
                                      range(len(self.data['Z'])), 
                                      self.data['Z'], 
                                      cmap='jet', 
                                      #vmin=0, #shading='auto'
                                      )
                                      
        
        am = self.ax.contour(range(len(self.data['Z'][0])), 
                             range(len(self.data['Z']))   , 
                             self.data['Z']               ,
                             colors='k')
        if self.color_bar:
            self.color_bar.remove() 
        self.color_bar = self.fig.colorbar(self.pcm,  ax=self.ax)#, extend='both')
        self.fig.canvas.draw()
    
    def on_coordinates_combobox_change (self, widget):
        """ Function doc """
        _id = self.coordinates_combobox.get_active()
        print(_id)
        vobject_index = None
        #-----------------------------------------------------------------------------
        _iter = self.coordinates_combobox.get_active_iter()
        if _iter is not None:
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.coordinates_combobox.get_model()
            _name, vobject_index = model[_iter][:2]
            print ('\n\n\_name, vobject_index:', _name, vobject_index, '\n\n')
        #-----------------------------------------------------------------------------
        self.vobject = self.main.vm_session.vm_objects_dic[vobject_index]

        self.data_liststore.clear()
        for index , data in enumerate(self.p_session.psystem[self.vobject.e_id].e_logfile_data[vobject_index]):
            #print(data)
            self.data_liststore.append([data['name'], index])
        
        
        #self.data_liststore.append(['all', 2])
        
        self.data_combobox.set_active(0)

    def on_data_combobox_change (self, widget):
        """ Function doc """

        _iter = self.data_combobox.get_active_iter()
        if _iter is not None:
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.data_combobox.get_model()
            _name, index = model[_iter][:2]
            #print ('\n\n\_name, index:', _name,  index, '\n\n')
        
        #self.vobject = self.main.vm_session.vm_objects_dic[vobject_index]
        self.data = self.p_session.psystem[self.vobject.e_id].e_logfile_data[self.vobject.index][index] 
        print(self.data)
        self._draw_data(cla = True)
        
    def on_motion_notify_event (self, event):
        """ Function doc """
        if event.xdata is None or event.ydata is None:
            pass
        else:
            x, y = int(event.xdata), int(event.ydata)
            try:
                rc1 = self.data['RC1'][y][x]
                rc2 = self.data['RC2'][y][x]
                energy = self.data['Z'][y][x]
                
                label_text =  '({:3d}) RC1: {:.4f}     ({:3d}) RC2: {:.4f}  '.format(x, rc1, y, rc2) 
                
                self.reaction_coord_label.set_text(label_text)
                print('{:10.5f} {:10.5f} {:4} {:4} {:13.5f}'.format(rc1, rc2, x,y, energy))
            except:
                pass
            
    def on_pick(self, event):
        ''' Function doc '''
        
        #self.ax.pcolormesh(self.X, self.Y, self.Z, cmap='jet', vmin=0, shading='gouraud')
        #self.ax.pcolormesh(self.X, self.Y, self.Z, cmap='jet', vmin=0, shading='gouraud')
        
        #self.pcm = self.ax.pcolormesh(self.X, self.Y, self.Z, cmap='jet', vmin=0)#, shading='gouraud')
        #am = self.ax.contour(self.X, self.Y, self.Z,lspacing = 10, colors='k')
        button_number = int(event.button)
        
        if button_number == 1:
            self.ax.cla()
            self.pcm = self.ax.pcolormesh(range(len(self.data['Z'][0])), 
                                          range(len(self.data['Z'])), 
                                          self.data['Z'], 
                                          cmap='jet',# shading='auto', 
                                          )#vmin=0)
                                          
            
            am = self.ax.contour(range(len(self.data['Z'][0])), 
                                 range(len(self.data['Z']))   , 
                                 self.data['Z']               ,
                                 colors='k')
            
            
            
            if event.xdata is None or event.ydata is None:
                self.xdata  = []
                self.ydata  = []
                self.indexes = []
                self.zdata   = []
                self.ax2.cla()
                
                self.xy_traj = []
                self.scale_traj_new_definitions()

            else:
                
                #print('you pressed', event.button)
                #print('you pressed', event.xdata, event.ydata)
                x, y = int(event.xdata), int(event.ydata)
                #print (x, y)
                print('you pressed', x, y, self.data['Z'][y][x])
                #print('you pressed', event)
                
                if self.interpolate:
                    if len(self.xy_traj) > 0 :
                        print([self.xy_traj[-1],[x,y] ])
                        xy_list = build_chain_of_states( [self.xy_traj[-1],[x,y] ])
                    else:
                        xy_list = [[x,y]]
                else: xy_list = [[x,y]]
                
                for xy in xy_list:
                    x = xy[0]
                    y = xy[1]
                    
                    #self.xdata.append(event.xdata)
                    #self.ydata.append(event.ydata)
                    self.xdata.append(x+0.5)
                    self.ydata.append(y+0.5)
                    self.xy_traj.append([x,y])
                    
                    self.scale_traj_new_definitions()
                    
                    self.zdata.append(self.data['Z'][y][x])
                #self.line, = self.ax.plot(event.xdata, event.ydata, '-o')  # plot the first row
            
            self.line, = self.ax.plot(self.xdata, self.ydata, '-ok')  # plot the first row
            self.ax2.cla()
            self.line2, = self.ax2.plot(range(0, len(self.zdata)), self.zdata, '-ob')  # plot the first row
            
            self.canvas.draw()
            self.canvas2.draw()

    def scale_traj_new_definitions(self):
        #self.scale_traj
        self.scale_traj.set_range(0, len(self.xy_traj))
        #self.scale_traj.set_value(self.vm_session.get_frame())

    def on_scaler_frame_value_changed (self, hscale, text= None,  data=None):
        """ Function doc """
        value = self.scale_traj.get_value()
        pos   = self.scale_traj.get_value_pos ()
        
        
        #print(self.xy_traj[int(value)])
        xy = self.xy_traj[int(value)]
        print(xy, self.zdata[int(value)])
        
        self.ax2.cla()
        self.ax3.cla()
        self.ax2.plot(range(0, len(self.zdata)), self.zdata, '-ob')
        self.ax3.plot( [int(value)], [self.zdata[int(value)]], '-or')
        
        if self.vobject:
            frame = self.vobject.idx_2D_xy[(xy[0], xy[1])]
            self.main.vm_session.frame = int(frame)
            #self.main.vm_session.set_frame(int(frame)) 
            self.main.vm_session.vm_glcore.updated_coords = True
            self.main.vm_session.vm_widget.queue_draw()
        self.canvas2.draw()
        #self.scale_traj.set_value(int(value))

    def on_button_export_trajectory (self, widget):
        """ Function doc """
        print (
        self.xdata   ,
        self.ydata   ,
        self.zdata   ,
        self.xy_traj ,
        self.vobject ,
        )
        
        #print (self.main.p_session.systems.keys())
        #print (self.vobject.easyhybrid_system_id)
        active_id = self.main.p_session.active_id
        
        new_vismol_object = self.main.p_session.generate_new_empty_vismol_object(system_id = self.vobject.e_id , 
                                                                                 name      = 'new_coordinates' )
        
        
        #self.main.p_session.active_id = self.vobject.e_id
        #system_id = self.vobject.e_id

        #frames     = []
        #frame      = []

        for xy in self.xy_traj:
            frame_number = self.vobject.idx_2D_xy[(xy[0], xy[1])]
            #frames.append(self.vobject.frames[frame_number])
            print(self.vobject.frames[frame_number])
            new_vismol_object.frames = np.vstack((new_vismol_object.frames, 
                                                  [self.vobject.frames[frame_number]]))
        
        self.traj_export_index += 1
        self.vm_session.main.main_treeview.refresh_number_of_frames()


    def change_check_button_reaction_coordinate (self, widget):
        """ Function doc """
        #radiobutton_bidimensional = self.builder.get_object('radiobutton_bidimensional')
        if self.builder.get_object('label_check_button_reaction_coordinate2').get_active():
            self.box_reaction_coordinate2.set_sensitive(True)
            self.builder.get_object('n_CPUs_spinbutton').set_sensitive(True)
            self.builder.get_object('n_CPUs_label').set_sensitive(True)
        else:
            self.box_reaction_coordinate2.set_sensitive(False)
            self.builder.get_object('n_CPUs_spinbutton').set_sensitive(False)
            self.builder.get_object('n_CPUs_label')     .set_sensitive(False)

    #======================================================================================

    def  update (self):
        """ Function doc """
        pass


#def parse_1D_scan_logfile (logfile):
#    """ Function doc """
#    
#
#def parse_2D_scan_logfile (logfile):
#    """ Function doc """
#    
#    data  =  open(logfile, 'r')
#    lines =  data.readlines()
#    xlist = []
#    ylist = []
#    zlist = []
#    
#    lastline = lines[-1].split()
#    x_size = int(lastline[0])
#    y_size = int(lastline[1])
#
#    
#    rows = y_size+1
#    cols = x_size+1
#     
#    Z       = [[0]*cols for _ in range(rows)]
#    RC1     = [[0]*cols for _ in range(rows)]
#    RC2     = [[0]*cols for _ in range(rows)]
#    
#    #print (zlist)
#    
#    for line in lines[1:]:
#        line2 = line.split()
#        x = int(line2[0])
#        y = int(line2[1])
#        #print(x,y, line2[-1])
#        
#        
#    
#        Z[y][x]       = float(line2[-1]) 
#        RC1[y][x]     = float(line2[-3]) 
#        RC2[y][x]     = float(line2[-2]) 
#    data.close()
#
#    data = {
#           'RC1': RC1,
#           'RC2': RC2,
#           'Z': Z
#           }
#    return data
#
def find_the_midpoint (coord1 , coord2):
	""" Function doc """
	#print (coord1, coord2)
	
	x = float(coord2[0] - coord1[0])
	x = (x/2)
	#print ('x', x)
	x = coord1[0] + x


	y = float(coord2[1] - coord1[1])
	y = (y/2)
	#print ('y', y)

	y = coord1[1] + y

	#print (x, y)
	#return [int(x), int(round(y))]
	return [int(x), int( y )]

def build_chain_of_states( input_coord):
 
    #print (input_coord)
    inset_point = True


    while inset_point == True:
        a = 0
        counter = 0

        while a == 0:

            try:
                coord1 =  input_coord[ counter  ]
                coord2 =  input_coord[ counter+1]

                inset_point = check_distance (coord1 , coord2)
                
                if inset_point == False:
                    counter += 1
                    #print ('inset_point == False')
                else:
                    midpoint = find_the_midpoint (coord1 , coord2)
                    #print counter, counter+1, midpoint, input_coord
                    input_coord.insert(counter+1, 0 )
                    input_coord[counter+1] = midpoint
            except:
                a = True
                #print input_coord
    return input_coord[1:]


def check_distance (coord1 , coord2):
	""" Function doc """
	x = float(coord2[0] - coord1[0])
	y = float(coord2[1] - coord1[1])
	d =  (x**2 + y**2)**0.5
	if d < 1.42:
		return False
	else:
		#print 'not too close'
		return True





def main(args):
    
    win = PotentialEnergyAnalysisWindow()
    win.OpenWindow()
    Gtk.main()
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))





























