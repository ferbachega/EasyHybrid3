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
    WIDTH= 20
    HEIGHT=20
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
        
        self.gl_parameters = self.vm_session.vm_config.gl_parameters

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
            self.window.connect('destroy-event', self.CloseWindow)

            
            self.elements_gtk_scrolled = self.builder.get_object('elements_gtk_scrolled')
            self.built_element_treeview()
            
            self.tmp_autosave             = self.builder.get_object('checkbox_tmp_autosave') 
            self.ask_autosave_and_unsaved = self.builder.get_object('checkbox_ask_autosave_and_unsaved')            
            #----------------------------------------------------------------------
            #. Paths
            #----------------------------------------------------------------------
            self.entry_startup_path = self.builder.get_object('entry_startup_path')
            self.entry_scratch_path = self.builder.get_object('entry_scratch_path')
            self.entry_tmp_path     = self.builder.get_object('entry_tmp_path')
            #----------------------------------------------------------------------

            #-------------------------------------------------------------------------------------
            self.set_interface_startup_shutdown_paramters()
            #-------------------------------------------------------------------------------------
            
            #-------------------------------------------------------------------------------------
            self.set_general_parameters()
            
            self.set_selection_parameters()
                        
            self.set_light_parameters()
            
            self.set_line_parameters()
            
            self.set_bond_parameters()
            
            self.set_stick_parameters()
            
            self.set_sphere_parameters()
            
            self.set_paths_and_folders()
            #-------------------------------------------------------------------------------------
            
            

            
            self.btn_apply_changes = self.builder.get_object('btn_apply_changes')
            self.btn_apply_changes.connect('clicked', self.on_btn_apply_all_changes)
            
            self.btn_apply_and_save_changes = self.builder.get_object('btn_apply_and_save_changes')
            self.btn_apply_and_save_changes.connect('clicked', self.on_btn_apply_and_save_changes)
            
            self.btn_reset_parms = self.builder.get_object('btn_reset_parameters')
            self.btn_reset_parms.connect('clicked', self.on_btn_reset_parms)
            
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
        print('self.visible',self.visible)
    
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
        self.numero_atomico_column = Gtk.TreeViewColumn("Number")
        self.simbolo_column = Gtk.TreeViewColumn("symbol")
        self.cor_column = Gtk.TreeViewColumn("Color")
        self.mass_column = Gtk.TreeViewColumn("Mass")
        self.rcov_column = Gtk.TreeViewColumn("r (cov)")
        self.rvdw_column = Gtk.TreeViewColumn("r (vdw)")

        # Adicionando as colunas à TreeView
        self.treeview.append_column(self.numero_atomico_column)
        self.treeview.append_column(self.simbolo_column)
        self.treeview.append_column(self.cor_column)
        self.treeview.append_column(self.mass_column)
        self.treeview.append_column(self.rcov_column)
        self.treeview.append_column(self.rvdw_column)


        # Criando os CellRenderers para exibir os dados nas colunas
        self.numero_atomico_cell = Gtk.CellRendererText()
        self.simbolo_cell = Gtk.CellRendererText()
        self.cor_cell = Gtk.CellRendererPixbuf()

        self.mass_cell = Gtk.CellRendererText()
        self.rcov_cell = Gtk.CellRendererText()
        self.rvdw_cell = Gtk.CellRendererText()

        # Adicionando os CellRenderers às colunas
        self.numero_atomico_column.pack_start(self.numero_atomico_cell, True)
        self.simbolo_column.pack_start(self.simbolo_cell, True)
        self.cor_column.pack_start(self.cor_cell, True)

        self.mass_column.pack_start(self.mass_cell, True)
        self.rcov_column.pack_start(self.rcov_cell, True)
        self.rvdw_column.pack_start(self.rvdw_cell, True)

        # Definindo os atributos dos CellRenderers para exibir os dados corretos
        self.numero_atomico_column.add_attribute(self.numero_atomico_cell, "text", 0)
        self.simbolo_column.add_attribute(self.simbolo_cell, "text", 1)
        self.cor_column.set_cell_data_func(self.cor_cell, self.render_color_square, None)


        self.mass_column.add_attribute(self.mass_cell, "text", 3)
        self.rcov_column.add_attribute(self.rcov_cell, "text", 4)
        self.rvdw_column.add_attribute(self.rvdw_cell, "text", 5)


        # Criando um modelo para os dados
        self.liststore = Gtk.ListStore(int, str, GdkPixbuf.Pixbuf, str, str, str, str)
        
        # Preenchendo o modelo com os dados dos átomos
        for symbol, data in self.vm_session.periodic_table.elements_by_symbol.items():
            numero_atomico = data[0]
            cor = self.get_color_pixbuf(data[2])
            mass = str(data[3])
            rdis = str(data[4])
            rcov = str(data[5])
            rvdw = str(data[6])
            name = data[1]
            self.liststore.append([numero_atomico, symbol, cor, mass, rcov, rvdw, name])
        
        
        # Conectando o modelo à TreeView
        self.treeview.set_model(self.liststore)
        #self.treeview.set_tooltip_cell(Gtk.Tooltip(), None, self.simbolo_cell, None, 1)
        self.elements_gtk_scrolled.add(self.treeview)

    def __apply_all_changes (self):
        """ Function doc """
        self.__apply_light_parameters()
        self.__apply_sphere_parameters()
        self.__apply_stick_parameters()
        self.__apply_lines_parameters() 
        self.__apply_bond_parameters() 
        self.__apply_viewer_selections_parameters() 
        self.__apply_viewer_general_parameters() 
        self.__apply_interface_general_parameters()
        self.vm_session.vm_glcore.queue_draw()
    
    def on_btn_reset_parms (self, widget):
        """ Function doc """
        #print('on_btn_reset_parms')
        isOK = self.vm_session.vm_config.reset_parameters()
     
        if isOK:
            #-------------------------------------------------------------------------------------
            self.set_interface_startup_shutdown_paramters()
            #-------------------------------------------------------------------------------------

            self.set_general_parameters()
            self.set_selection_parameters()
            self.set_light_parameters()
            self.set_line_parameters()
            self.set_bond_parameters()
            self.set_stick_parameters()
            self.set_sphere_parameters()
            #-------------------------------------------------------------------------------------
            self.vm_session.vm_glcore.queue_draw()
    
    def on_btn_apply_all_changes (self, widget):
        """ Function doc """
        self.__apply_all_changes()
    
    def on_btn_apply_and_save_changes (self, widget):
        """ Function doc """
        self.__apply_all_changes()
        #self.__apply_light_parameters()
        #self.__apply_sphere_parameters()
        #self.__apply_stick_parameters()
        #self.__apply_lines_parameters() 
        #self.__apply_bond_parameters() 
        #self.__apply_viewer_selections_parameters() 
        #self.__apply_viewer_general_parameters() 
        #self.__apply_interface_general_parameters()
        self.vm_session.vm_config.save_easyhybrid_config()
    
    def color_set_viewing_selections (self, widget):
        """ 
        not used anymore       
        """
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

    
    
    
    def set_interface_startup_shutdown_paramters (self):
        """ Function doc """
        a = self.vm_session.vm_config.gl_parameters['autosave']      
        b = self.vm_session.vm_config.gl_parameters['askSaveUnsave'] 
        
        self.tmp_autosave            .set_active(a)
        self.ask_autosave_and_unsaved.set_active(b)
        pass
    
    def set_general_parameters (self):
        """ Function doc """
        #       background
        self.colorbutton = self.builder.get_object('btn_background_color')
        bgcolor = self.vm_session.vm_config.gl_parameters["background_color"]
        rgba = Gdk.RGBA(bgcolor[0], bgcolor[1], bgcolor[2])
        self.colorbutton.set_rgba(rgba)

        
        #   scroll step
        self.entry_scroll_step = self.builder.get_object('entry_scroll_step')
        scroll_step = self.vm_session.vm_config.gl_parameters["scroll_step"]
        self.entry_scroll_step.set_text(str(scroll_step))
        
        self.entry_sleep_time_coc = self.builder.get_object('entry_sleep_time_coc')
        sleep_time_coc = self.vm_session.vm_config.gl_parameters["center_on_coord_sleep_time"]
        self.entry_sleep_time_coc.set_text(str(sleep_time_coc))
        
        
        
        self.entry_field_of_view = self.builder.get_object('entry_field_of_view')
        field_of_view            = self.vm_session.vm_config.gl_parameters["field_of_view"]
        self.entry_field_of_view.set_text(str(field_of_view))

    def set_selection_parameters (self):
        """ Function doc """
        #-------------------------------------------------------------------------------------
        #                          viewing colors
        #-------------------------------------------------------------------------------------
        color = self.vm_session.vm_config.gl_parameters["picking_dots_color"]
        rgba = Gdk.RGBA(color[0], color[1], color[2])
        self.colorbutton_viewing_selections = self.builder.get_object('color_btn_view_sel_dots')
        self.colorbutton_viewing_selections.set_rgba(rgba)   
        
        self.entry_viewing_dot_size = self.builder.get_object('entry_view_sel_dots')
        #-------------------------------------------------------------------------------------
        
        #-------------------------------------------------------------------------------------
        self.entry_view_sel_dot_size = self.builder.get_object('entry_view_sel_dot_size')
        dot_sel_size = self.vm_session.vm_config.gl_parameters['dot_sel_size']
        self.entry_view_sel_dot_size.set_text(str(dot_sel_size))
        #-------------------------------------------------------------------------------------

        
        #-------------------------------------------------------------------------------------
        #                          picking selections
        #-------------------------------------------------------------------------------------
        color = self.vm_session.vm_config.gl_parameters["pk_label_color"]
        rgba = Gdk.RGBA(color[0], color[1], color[2])
        self.color_btn_picking_labels  = self.builder.get_object('color_btn_pk_label')
        self.color_btn_picking_labels.set_rgba(rgba)
        #-------------------------------------------------------------------------------------
        
        
        #-------------------------------------------------------------------------------------
        #                    picking selections dashed list
        #-------------------------------------------------------------------------------------
        color = self.vm_session.vm_config.gl_parameters["pk_dist_label_color"]
        rgba = Gdk.RGBA(color[0], color[1], color[2])
        self.color_btn_pk_dist_label     = self.builder.get_object('color_btn_pk_dist_label')
        self.color_btn_pk_dist_label.set_rgba(rgba)
        #-------------------------------------------------------------------------------------
        color = self.vm_session.vm_config.gl_parameters["dashed_dist_lines_color"]
        rgba = Gdk.RGBA(color[0], color[1], color[2])
        self.color_btn_pk_dist_lines     = self.builder.get_object('color_btn_pk_dist_lines')
        self.color_btn_pk_dist_lines.set_rgba(rgba)
        #-------------------------------------------------------------------------------------



        #-------------------------------------------------------------------------------------
        #                              picking spheres
        #-------------------------------------------------------------------------------------
        color =  self.vm_session.picking_selections.pk_scolor['pk1']
        rgba = Gdk.RGBA(color[0], color[1], color[2])
        self.color_button_pk1_sphr = self.builder.get_object('color_btn_pk_sphr_1')
        self.color_button_pk1_sphr.set_rgba(rgba)
        #-------------------------------------------------------------------------------------
        color =  self.vm_session.picking_selections.pk_scolor['pk2']
        rgba = Gdk.RGBA(color[0], color[1], color[2])
        self.color_button_pk2_sphr = self.builder.get_object('color_btn_pk_sphr_2')
        self.color_button_pk2_sphr.set_rgba(rgba)
        #-------------------------------------------------------------------------------------
        color =  self.vm_session.picking_selections.pk_scolor['pk3']
        rgba = Gdk.RGBA(color[0], color[1], color[2])
        self.color_button_pk3_sphr = self.builder.get_object('color_btn_pk_sphr_3')
        self.color_button_pk3_sphr.set_rgba(rgba)
        #-------------------------------------------------------------------------------------
        color =  self.vm_session.picking_selections.pk_scolor['pk4']
        rgba = Gdk.RGBA(color[0], color[1], color[2])
        self.color_button_pk4_sphr = self.builder.get_object('color_btn_pk_sphr_4')
        self.color_button_pk4_sphr.set_rgba(rgba)
        #-------------------------------------------------------------------------------------
        


        self.entry_pk_label_size      =  self.builder.get_object('entry_pk_label_size')
        self.entry_pk_dist_label_size =  self.builder.get_object('entry_pk_dist_label_size')
        self.entry_pk_dist_line_size  =  self.builder.get_object('entry_pk_dist_line_size')
        
    def set_line_parameters (self):
        """ Function doc """
        gridsize             = self.vm_session.vm_config.gl_parameters['gridsize']              
        line_width_selection = self.vm_session.vm_config.gl_parameters['line_width_selection']
        line_type            = self.vm_session.vm_config.gl_parameters['line_type']
        
        self.entry_lines_with            =  self.builder.get_object('entry_lines_with')
        self.entry_lines_with_selections =  self.builder.get_object('entry_lines_with_selections')
        self.entry_lines_type            =  self.builder.get_object('entry_lines_type')

        self.entry_lines_with           .set_text(str(gridsize))   
        self.entry_lines_with_selections.set_text(str(line_width_selection))   
        self.entry_lines_type           .set_text(str(line_type))   
    
    def set_stick_parameters (self):
        """ Function doc """
        sticks_radius = self.vm_session.vm_config.gl_parameters['sticks_radius']
        sticks_type   = self.vm_session.vm_config.gl_parameters['sticks_type']
        
        self.entry_stick_radius = self.builder.get_object('entry_stick_radius')         
        self.entry_stick_type   = self.builder.get_object('entry_stick_type')          
        
        self.entry_stick_radius.set_text(str(sticks_radius))
        self.entry_stick_type  .set_text(str(sticks_type))
        
        
        self.btn_stick_unique_color       = self.builder.get_object('btn_stick_unique_color')     
        self.checkbox_stick_unique_color  = self.builder.get_object('checkbox_stick_unique_color')
        
        self.entry_dbond_radius           = self.builder.get_object('entry_dbond_radius')
        self.entry_dbond_type             = self.builder.get_object('entry_dbond_type')
        self.btn_dbond_unique_color       = self.builder.get_object('btn_dbond_unique_color')
        self.checkbox_dbond_unique_color  = self.builder.get_object('checkbox_dbond_unique_color')
        
        self.entry_dbond_radius.set_text(str(sticks_radius))
        self.entry_dbond_type  .set_text(str(sticks_type))

    def set_light_parameters (self):
        """ Function doc """
        
 
        self.entry_light_position         = self.builder.get_object('entry_light_position')
        light_position         = self.vm_session.vm_config.gl_parameters['light_position']
        self.entry_light_position.set_text(str(light_position))
        
        
        color            = self.vm_session.vm_config.gl_parameters['light_color']
        self.btn_light_color              = self.builder.get_object('btn_light_color')
        rgba = Gdk.RGBA(color[0], color[1], color[2])
        self.btn_light_color.set_rgba(rgba)

        
        light_ambient_coef     = self.vm_session.vm_config.gl_parameters['light_ambient_coef']
        self.entry_light_amb_coef           = self.builder.get_object('entry_light_amb_coef')
        self.entry_light_amb_coef.set_text(str(light_ambient_coef))
        
        
        self.entry_light_shiness  = self.builder.get_object('entry_light_shiness')
        light_shininess           = self.vm_session.vm_config.gl_parameters['light_shininess']
        self.entry_light_shiness.set_text(str(light_shininess))
        
        
        self.btn_light_spec_color         = self.builder.get_object('btn_light_spec_color')
        color = self.vm_session.vm_config.gl_parameters['light_specular_color']
        rgba = Gdk.RGBA(color[0], color[1], color[2])
        self.btn_light_spec_color.set_rgba(rgba)
        
        
        self.entry_light_intensity = self.builder.get_object('entry_light_intensity')
        light_intensity        = self.vm_session.vm_config.gl_parameters['light_intensity']
        self.entry_light_intensity.set_text(str(light_intensity))
        
    def set_bond_parameters (self):
        """ Function doc """
        gridsize        = self.vm_session.vm_config.gl_parameters['gridsize']
        maxbond         = self.vm_session.vm_config.gl_parameters['maxbond']
        bond_tolerance  = self.vm_session.vm_config.gl_parameters['bond_tolerance']

        self.entry_grid_size     = self.builder.get_object('entry_grid_size')
        self.entry_max_bond_size = self.builder.get_object('entry_max_bond_size')
        self.entry_bond_tol      = self.builder.get_object('entry_bond_tol')
        
        self.entry_grid_size    .set_text(str(gridsize))
        self.entry_max_bond_size.set_text(str(maxbond))
        self.entry_bond_tol     .set_text(str(bond_tolerance))
        
    def set_sphere_parameters (self):
        """ Function doc """

        sphere_type      = self.vm_session.vm_config.gl_parameters['sphere_type']
        sphere_quality   = self.vm_session.vm_config.gl_parameters['sphere_quality']
        sphere_scale     = self.vm_session.vm_config.gl_parameters['sphere_scale']

        self.entry_sphere_scale   = self.builder.get_object('entry_sphere_scale')
        self.entry_sphere_quality = self.builder.get_object('entry_sphere_quality')
        self.entry_sphere_type    = self.builder.get_object('entry_sphere_type')

        self.entry_sphere_scale  .set_text(str(sphere_scale))
        self.entry_sphere_quality.set_text(str(sphere_quality))
        self.entry_sphere_type   .set_text(str(sphere_type))

    def set_paths_and_folders (self):
        """ Function doc """
        startup_path = self.vm_session.vm_config.gl_parameters['startup_path']
        self.entry_startup_path.set_text(startup_path)
    
    
    
    def __apply_light_parameters (self):
        #---------------------------------------------------------------
        light_pos = self.entry_light_position.get_text()
        light_pos = light_pos.replace('[', '')
        light_pos = light_pos.replace(']', '')
        light_pos = light_pos.split(',')
        
        new_pos = []
        for pos in light_pos:
            new_pos.append(float(pos))
        self.gl_parameters['light_position'] = new_pos
        #---------------------------------------------------------------

        
        #---------------------------------------------------------------
        color = self.btn_light_color.get_rgba()
        color = list(color)
        color = color[:-1]
        self.gl_parameters["light_color"] = color
        #---------------------------------------------------------------

        
        #---------------------------------------------------------------
        amb_coef = float(self.entry_light_amb_coef.get_text())
        self.gl_parameters["light_ambient_coef"] = amb_coef
        #---------------------------------------------------------------

        
        #---------------------------------------------------------------
        light_shiness = float(self.entry_light_shiness.get_text())
        self.gl_parameters["light_shininess"] = light_shiness
        #---------------------------------------------------------------
        
        #---------------------------------------------------------------
        light_intensity =  self.entry_light_intensity.get_text() 
        light_intensity = light_intensity.replace('[', '')
        light_intensity = light_intensity.replace(']', '')
        light_intensity = light_intensity.split(',')
        
        new_int = []
        for pos in light_intensity:
            new_int.append(float(pos))
        self.gl_parameters['light_intensity'] = new_int
        #---------------------------------------------------------------
        
        #---------------------------------------------------------------
        color = self.btn_light_spec_color.get_rgba()
        color = list(color)
        color = color[:-1]
        self.gl_parameters["light_specular_color"] = color
        #---------------------------------------------------------------
    def __apply_sphere_parameters (self):
        """ Function doc """

        sphere_type    = int(self.entry_sphere_type   .get_text())
        sphere_quality = int(self.entry_sphere_quality.get_text())
        sphere_scale   = float(self.entry_sphere_scale .get_text())
        
        self.gl_parameters['sphere_type']    = sphere_type    
        self.gl_parameters['sphere_quality'] = sphere_quality 
        self.gl_parameters['sphere_scale']   = sphere_scale   

    def __apply_stick_parameters (self):
        """ Function doc """
        stick_radius = float(self.entry_stick_radius.get_text())
        stick_type   = float(self.entry_stick_type  .get_text())
        
        dbond_radius = float(self.entry_dbond_radius.get_text())
        dbond_type   = float(self.entry_dbond_type  .get_text())
        
        self.gl_parameters['sticks_radius'] = stick_radius
        self.gl_parameters['sticks_type']   = stick_type  
        #self.gl_parameters[]
        #self.gl_parameters[]

    def __apply_lines_parameters (self):
        lines_with            = float(self.entry_lines_with           .get_text())
        lines_with_selections = float(self.entry_lines_with_selections.get_text())
        lines_type            = float(self.entry_lines_type           .get_text())
        
        self.gl_parameters['line_width']           = lines_with           
        self.gl_parameters['line_width_selection'] = lines_with_selections
        self.gl_parameters['line_type']            = lines_type           

    def __apply_bond_parameters (self):
        """ Function doc """
        #try:
        grid_size     = float(self.entry_grid_size    .get_text())
        max_bond_size = float(self.entry_max_bond_size.get_text())
        bond_tol      = float(self.entry_bond_tol     .get_text())
        
        self.gl_parameters['gridsize']       = grid_size    
        self.gl_parameters['maxbond']        = max_bond_size
        self.gl_parameters['bond_tolerance'] = bond_tol     
        
    def __apply_viewer_selections_parameters (self):
        """ Function doc """
        #-----------------------------------------------------------------------
        #                            viewing  dots
        #-----------------------------------------------------------------------        
        widget = self.builder.get_object('color_btn_view_sel_dots')
        color = widget.get_rgba()
        color = list(color)
        color = color[:-1]
        self.vm_session.vm_config.gl_parameters["picking_dots_color"] = color
        for vm_object in self.vm_session.vm_objects_dic.values():
            vm_object.core_representations["picking_dots"] = None
        #-----------------------------------------------------------------------
        
        
        
        
        #-----------------------------------------------------------------------
        #                     pk labels  / distance labels
        #-----------------------------------------------------------------------
        btn_pk_label = self.builder.get_object('color_btn_pk_label')
        color = btn_pk_label.get_rgba()
        color = list(color)
        #color = color[:-1]
        self.vm_session.vm_config.gl_parameters["pk_label_color"] = color
        self.vm_session.vm_glcore.vm_font.vao = None
        #-----------------------------------------------------------------------

        #-----------------------------------------------------------------------
        widget = self.builder.get_object('color_btn_pk_dist_label')
        color = widget.get_rgba()
        color = list(color)
        #color = color[:-1]
        self.vm_session.vm_config.gl_parameters["pk_dist_label_color"] = color
        self.vm_session.vm_glcore.vm_font_dist.vao = None
        #-----------------------------------------------------------------------

        #-----------------------------------------------------------------------
        widget = self.builder.get_object('color_btn_pk_dist_lines')
        color = widget.get_rgba()
        color = list(color)
        color = color[:-1]
        self.vm_session.vm_config.gl_parameters["dashed_dist_lines_color"] = color
        #self.vm_session.vm_glcore.vm_font_dist.vao = None
        
        
        #-----------------------------------------------------------------------
        

        #-------------------------------------------------------------------------
        #                         pk Spheres
        #-------------------------------------------------------------------------
        widget = self.builder.get_object('color_btn_pk_sphr_1')
        color = widget.get_rgba()
        color = list(color)
        color = color[:-1]
        self.vm_session.set_pk_sphr_selection_color(color, 1)

        widget = self.builder.get_object('color_btn_pk_sphr_2')
        color = widget.get_rgba()
        color = list(color)
        color = color[:-1]
        self.vm_session.set_pk_sphr_selection_color(color, 2)
        
        widget = self.builder.get_object('color_btn_pk_sphr_3')
        color = widget.get_rgba()
        color = list(color)
        color = color[:-1]
        self.vm_session.set_pk_sphr_selection_color(color, 3)
        
        widget = self.builder.get_object('color_btn_pk_sphr_4')
        color = widget.get_rgba()
        color = list(color)
        color = color[:-1]
        self.vm_session.set_pk_sphr_selection_color(color, 4)
        #-------------------------------------------------------------------------
    def __apply_viewer_general_parameters (self):
        #---------------------------------------------------------------
        self.btn_background_color = self.builder.get_object('btn_background_color')
        color = self.btn_background_color.get_rgba()
        color = list(color)
        #color = color[:-1]
        self.gl_parameters["background_color"] = color
        self.vm_session.vm_glcore.bckgrnd_color = color
        #---------------------------------------------------------------

        
        
        #---------------------------------------------------------------
        scroll_step    = float(self.builder.get_object('entry_scroll_step')   .get_text() )
        sleep_time_coc = float(self.builder.get_object('entry_sleep_time_coc').get_text() )
        field_of_view  = float(self.builder.get_object('entry_field_of_view') .get_text() )
        
        self.gl_parameters["scroll_step"]                = scroll_step    
        self.gl_parameters["center_on_coord_sleep_time"] = sleep_time_coc 
        self.gl_parameters["field_of_view"]              = field_of_view  
        #---------------------------------------------------------------
    def __apply_interface_general_parameters (self):
        """ Function doc """
        path = self.entry_startup_path.get_text()
        
        if os.path.isdir(path):
            self.vm_session.vm_config.gl_parameters['startup_path'] = path
            print('Defining New Startup Path:', path)
        else:
            dialog = Gtk.MessageDialog(
                                flags=0,
                                message_type=Gtk.MessageType.ERROR,
                                buttons=Gtk.ButtonsType.OK,
                                text="ERROR: Folder not found.",
                            )

            dialog.format_secondary_text(
                                    "If the desired path is correct, create the folder using your file manager."
                                        )
            dialog.run()
            #print("ERROR dialog closed")
            dialog.destroy()
        
        
        a = self.builder.get_object('checkbox_tmp_autosave').get_active()
        b = self.builder.get_object('checkbox_ask_autosave_and_unsaved').get_active()
        c = self.builder.get_object('checkbox_save_window_size').get_active()
        
        
        self.vm_session.vm_config.gl_parameters['autosave']      = a
        self.vm_session.vm_config.gl_parameters['askSaveUnsave'] = b
        print(a,b)
        #self.vm_session.vm_config.gl_parameters['startup_path'] =
        
        
        
        #
        #b = self.builder.get_object('entry_startup_path')
        #b = self.builder.get_object('entry_scratch_path')
        #b = self.builder.get_object('entry_tempfiles_path')






















