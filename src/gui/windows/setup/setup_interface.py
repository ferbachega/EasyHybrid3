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



#'''
def get_colorful_square_pixel_buffer (system = None,  atomic_symbol = 'C'):
    """ Function doc """
    if system is not None:
        color        =  system.e_color_palette[atomic_symbol]
        res_color    = [int(color[0]*255),int(color[1]*255),int(color[2]*255)] 
        pixelbuffer  =  getColouredPixmap( res_color[0], res_color[1], res_color[2] )
    
    else:
        res_color    = [255,255,255] 
        pixelbuffer  =  getColouredPixmap( res_color[0], res_color[1], res_color[2] )
    return pixelbuffer


def getColouredPixmap( r, g, b, a=255 ):
    """ Given components, return a colour swatch pixmap """
    CHANNEL_BITS=8
    WIDTH=10
    HEIGHT=10
    swatch = GdkPixbuf.Pixbuf.new( GdkPixbuf.Colorspace.RGB, True, CHANNEL_BITS, WIDTH, HEIGHT ) 
    swatch.fill( (r<<24) | (g<<16) | (b<<8) | a ) # RGBA
    return swatch
#'''

class EasyHybridPreferencesWindow():
    """ Class doc """
    def __init__(self, main = None):
        """ Class initialiser """
        self.main_session          = main#self.main_session.system_liststore
        self.home                  = main.home
        self.visible               = False        
        self.p_session             = main.p_session
        self.vm_session            = main.vm_session

    def OpenWindow (self):
        """ Function doc /home/fernando/programs/VisMol/easyhybrid/gui/selection_list.glade"""
        if self.visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'src/gui/windows/setup/setup_interface.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('window')
            self.window.set_default_size(600, 600)  
            self.window.set_title('Preferences')  
            self.window.set_keep_above(True)
            
            
            self.elements_gtk_scrolled = self.builder.get_object('elements_gtk_scrolled')
            self.built_element_treeview()
            
            
            #       background
            self.colorbutton = self.builder.get_object('gtk_color_button_background')
            bgcolor = self.vm_session.vm_config.gl_parameters["background_color"]
            rgba = Gdk.RGBA(bgcolor[0], bgcolor[1], bgcolor[2])
            self.colorbutton.set_rgba(rgba)
            
            
            #-------------------------------------------------------------------------------------
            #     viewing colors
            rgba = Gdk.RGBA(0, 1, 1)
            self.colorbutton_viewing_selections = self.builder.get_object('color_btn_view_sel_dots')
            self.colorbutton_viewing_selections.set_rgba(rgba)
            self.colorbutton_viewing_selections.connect('color-set', self.color_set_viewing_selections)
        
            self.entry_viewing_dot_size = self.builder.get_object('entry_view_sel_dots')
            #-------------------------------------------------------------------------------------


            #-------------------------------------------------------------------------------------
            #     picking selections
            #bgcolor = self.vm_session.vm_config.gl_parameters["background_color"]
            #rgba = Gdk.RGBA(bgcolor[0], bgcolor[1], bgcolor[2])
            self.color_btn_picking_labels  = self.builder.get_object('color_btn_pk_label')
            self.color_btn_picking_labels.connect('color-set', self.color_set_viewing_selections)

            self.color_btn_pk_dist_label     = self.builder.get_object('color_btn_pk_dist_label')
            self.color_btn_pk_dist_label.connect('color-set', self.color_set_viewing_selections)

            self.color_btn_pk_dist_lines     = self.builder.get_object('color_btn_pk_dist_lines')
            self.color_btn_pk_dist_lines.connect('color-set', self.color_set_viewing_selections)
            #-------------------------------------------------------------------------------------
            
            
            #-------------------------------------------------------------------------------------
            self.color_button_pk1_sphr = self.builder.get_object('color_btn_pk_sphr_1')
            self.color_button_pk1_sphr.connect('color-set', self.color_set_viewing_selections)

            self.color_button_pk2_sphr = self.builder.get_object('color_btn_pk_sphr_2')
            self.color_button_pk2_sphr.connect('color-set', self.color_set_viewing_selections)

            self.color_button_pk3_sphr = self.builder.get_object('color_btn_pk_sphr_3')
            self.color_button_pk3_sphr.connect('color-set', self.color_set_viewing_selections)

            self.color_button_pk4_sphr = self.builder.get_object('color_btn_pk_sphr_4')
            self.color_button_pk4_sphr.connect('color-set', self.color_set_viewing_selections)
            #-------------------------------------------------------------------------------------
            
            
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
    def get_color_pixbuf(self, rgb_values):
        rgb = rgb_values
        
        rgb
        #pixbuf = getColouredPixmap(r =155, g=155, b=155, a=255)
        pixbuf = getColouredPixmap(r = round(rgb[0]*255), g=round(rgb[1]*255), b=round(rgb[2]*255), a=255)
        
        #pixel_data = bytes(rgb_values * 4)  # RGBA format
        #pixbuf = GdkPixbuf.Pixbuf.new_from_data(pixel_data, GdkPixbuf.Colorspace.RGB, True, 8, 1, 1, 4)
        return pixbuf

    def render_color_square(self, column, cell, model, iter, data):
        color_pixbuf = model.get_value(iter, 2)
        cell.set_property("pixbuf", color_pixbuf)

    def built_element_treeview (self):
        """ Function doc """
         # Criando uma TreeView
        self.treeview = Gtk.TreeView()

        # Criando as colunas
        self.numero_atomico_column = Gtk.TreeViewColumn("Número Atômico")
        self.simbolo_column = Gtk.TreeViewColumn("Símbolo")
        self.cor_column = Gtk.TreeViewColumn("Cor")

        # Adicionando as colunas à TreeView
        self.treeview.append_column(self.numero_atomico_column)
        self.treeview.append_column(self.simbolo_column)
        self.treeview.append_column(self.cor_column)

        # Criando os CellRenderers para exibir os dados nas colunas
        self.numero_atomico_cell = Gtk.CellRendererText()
        self.simbolo_cell = Gtk.CellRendererText()
        self.cor_cell = Gtk.CellRendererPixbuf()

        # Adicionando os CellRenderers às colunas
        self.numero_atomico_column.pack_start(self.numero_atomico_cell, True)
        self.simbolo_column.pack_start(self.simbolo_cell, True)
        self.cor_column.pack_start(self.cor_cell, True)

        # Definindo os atributos dos CellRenderers para exibir os dados corretos
        self.numero_atomico_column.add_attribute(self.numero_atomico_cell, "text", 0)
        self.simbolo_column.add_attribute(self.simbolo_cell, "text", 1)
        self.cor_column.set_cell_data_func(self.cor_cell, self.render_color_square, None)

        # Criando um modelo para os dados
        self.liststore = Gtk.ListStore(int, str, GdkPixbuf.Pixbuf)
        
        # Preenchendo o modelo com os dados dos átomos
        for symbol, data in self.vm_session.periodic_table.elements_by_symbol.items():
            numero_atomico = data[0]
            cor = self.get_color_pixbuf(data[2])
            self.liststore.append([numero_atomico, symbol, cor])

        # Conectando o modelo à TreeView
        self.treeview.set_model(self.liststore)
        self.elements_gtk_scrolled.add(self.treeview)



    def color_set_viewing_selections (self, widget):
        """ Function doc """
        color = widget.get_rgba()
        print("Selected color: ", list(color))
        color = list(color)
        
        #-----------------------------------------------------------------------
        #                            viewing  dots
        #-----------------------------------------------------------------------
        if widget == self.colorbutton_viewing_selections:
            #color = list(color)
            color = color[:-1]
            self.vm_session.vm_config.gl_parameters["picking_dots_color"] = color
            for vm_object in self.vm_session.vm_objects_dic.values():
                vm_object.core_representations["picking_dots"] = None
        #-----------------------------------------------------------------------
        
        
        
        
        #-----------------------------------------------------------------------
        #                     pk labels  / distance labels
        #-----------------------------------------------------------------------
        elif widget == self.color_btn_picking_labels:
            self.vm_session.vm_config.gl_parameters["pk_label_color"] = color
            self.vm_session.vm_glcore.vm_font.vao = None
        
        elif widget == self.color_btn_pk_dist_label:
            self.vm_session.vm_config.gl_parameters["pk_dist_label_color"] = color
            self.vm_session.vm_glcore.vm_font_dist.vao = None
            pass
        
        elif widget == self.color_btn_pk_dist_lines:
            print(color)
            self.vm_session.vm_config.gl_parameters["dashed_dist_lines_color"] = color
            #self.vm_session.vm_glcore.vm_font_dist.vao = None
            pass
        #-----------------------------------------------------------------------
        
        
        
        
        #-------------------------------------------------------------------------
        #                         pk Spheres
        #-------------------------------------------------------------------------
        elif widget == self.color_button_pk1_sphr:
            self.vm_session.set_pk_sphr_selection_color(color, 1)
            pass
        elif widget == self.color_button_pk2_sphr:
            self.vm_session.set_pk_sphr_selection_color(color, 2)
            pass
        elif widget == self.color_button_pk3_sphr:
            self.vm_session.set_pk_sphr_selection_color(color, 3)
            pass
        elif widget == self.color_button_pk4_sphr:
            self.vm_session.set_pk_sphr_selection_color(color, 4)
            pass
        #-------------------------------------------------------------------------


        return color
