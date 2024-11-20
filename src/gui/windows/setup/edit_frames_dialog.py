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
from gi.repository import Gtk, GdkPixbuf, Gdk
#from GTKGUI.gtkWidgets.filechooser import FileChooser
#from easyhybrid.pDynamoMethods.pDynamo2Vismol import *
import gc
import os
import numpy as np

from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')


class EditFrameDialog():
    """ Class doc """
    def __init__(self, main = None):
        """ Class initialiser """
        self.main_session          = main#self.main_session.system_liststore
        self.home                  = main.home
        self.visible               = False        
        self.p_session             = main.p_session
        self.vm_session            = main.vm_session
        
        self.gl_parameters = self.vm_session.vm_config.gl_parameters

    def OpenWindow (self, vobj_id):
        """ Function doc /home/fernando/programs/VisMol/easyhybrid/gui/selection_list.glade"""
        self.vobj_id = vobj_id
        
        if self.visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'src/gui/windows/setup/edit_frames_dialog.glade'))
            self.window = self.builder.get_object('dialog')
            self.window.connect('destroy', self.CloseWindow) 
            
            self.entry_delete_frames = self.builder.get_object('entry_delete_frames')
            self.button_cancel = self.builder.get_object('button_cancel')
            self.button_apply  = self.builder.get_object('button_apply')
            
            self.button_cancel.connect('clicked', self.CloseWindow)
            self.button_apply.connect('clicked', self.apply)
            
            self.check_button_reverse = self.builder.get_object('check_button_reverse')
            self.check_create_new_vobject = self.builder.get_object('check_create_new_vobject')
            
            self.window.show_all()                                               
            self.visible    =  True
            '''--------------------------------------------------------------------------------------------'''

        else:
            self.window.present()
            
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        #self.BackUpWindowData()
        self.window.destroy()
        self.visible    =  False
        #print('self.visible',self.visible)
    
    def apply (self, widget):
        """ Function doc """
        text = self.entry_delete_frames.get_text()
        print(text,  self.vobj_id )
        vobject = self.main_session.vm_session.vm_objects_dic[self.vobj_id]
        frames  = vobject.frames
        
        
        indices_to_remove = []
        
        #processing the text data (frames to be deleted)
        rawdata = text.split(',')
        
        for data in rawdata:
           
            # este parte  processa um range de indexes que deve ser deletado
            # exemplot  2-5 (intervalo entre 2 e 5).
            if ':' in data:
                index_range = data.split(':')
                
                if len(index_range) == 2:
                    first   = int(index_range[0])
                    last    = int(index_range[1])
                    indexes = range(first, last)
                    
                    for index in indexes:
                        if index in indices_to_remove:
                            pass
                        else:
                            indices_to_remove.append(index)
                else:
                    pass
            
            # esta parte processa os indexes que deve ser deletados.
            else:
                data = data.strip()
                
                if data == '':
                    pass
                else:
                    index = int(data)
                    if index in indices_to_remove:
                        pass
                    else:
                        indices_to_remove.append(index)
                
        
        if self.check_create_new_vobject.get_active():
            system = self.main_session.p_session.psystem[vobject.e_id]
            vobject = self.main_session.p_session._add_vismol_object_to_easyhybrid_session (system, show_molecule=True, name = 'edited_coords')


        coords = np.delete(frames, indices_to_remove, axis=0)
        #print( frames.shape )
        vobject.frames = coords

        if self.check_button_reverse.get_active():
            #reverser (invert order)
            coords_inverted = coords[::-1]
            vobject.frames = coords_inverted
            print('do reverse')
        
        
        # Apply fixed representation to the VisMol object
        self.main_session.p_session._apply_fixed_representation_to_vobject(vismol_object =vobject)
        
        # Apply QC representation to the VisMol object
        self.main_session.p_session._apply_QC_representation_to_vobject(vismol_object =vobject)
        
        # Refresh the widgets in the main window
        self.main_session.main_treeview.refresh_number_of_frames()
        #self.main_session.p_session.main.refresh_widgets()        
        
        
        
        
        self.CloseWindow(None, None)
