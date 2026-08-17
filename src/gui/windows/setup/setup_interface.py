#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Selection utilities for pDynamo systems
#
#  Copyright 2022-2025 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
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
#  Maintainer:
#      Fernando Bachega <ferbachega@gmail.com> or <easyhybrid3@gmail.com>
#
#  Description:
#      Provides functions for selecting atoms and residues in pDynamo systems
#      to facilitate QM/MM partitioning and molecular simulations.
#
from util.debug import dprint
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
from vismol.libgl.vismol_font import list_available_fonts
from vismol.libgl.vismol_font import resolve_font_path
from vismol.libgl.vismol_font import DEFAULT_FONT_FILE, DEFAULT_FONT_SIZE
from vismol.libgl.representations import compute_atom_label_text
from pprint import pprint
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
        self.set_paths_and_folders_in_parameters()
        
        
        
    def open_window (self):
        """ Function doc /home/fernando/programs/VisMol/easyhybrid/gui/selection_list.glade"""
        if self.visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'src/gui/windows/setup/setup_interface.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('window')
            self.window.set_default_size(600, 600)  
            self.window.set_title('Preferences')  
            self.window.set_keep_above(True)
            self.window.connect('destroy-event', self.close_window)

            
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

            # [NOVO] estado inicial dos controles de cor unica das dynamic bonds
            self.set_dynamic_bonds_parameters()
            
            self.set_sphere_parameters()
            
            self.set_paths_and_folders()
            #-------------------------------------------------------------------------------------
            
            

            
            self.btn_apply_changes = self.builder.get_object('btn_apply_changes')
            self.btn_apply_changes.connect('clicked', self.on_btn_apply_all_changes)
            
            self.btn_apply_and_save_changes = self.builder.get_object('btn_apply_and_save_changes')
            self.btn_apply_and_save_changes.connect('clicked', self.on_btn_apply_and_save_changes)
            
            self.btn_reset_parms = self.builder.get_object('btn_reset_parameters')
            self.btn_reset_parms.connect('clicked', self.on_btn_reset_parms)
            
            
            self.btn_cancel = self.builder.get_object('btn_cancel')
            self.btn_cancel.connect('clicked', self.close_window)

            self.window.show_all()                                               
            self.visible    =  True
            '''--------------------------------------------------------------------------------------------'''

        else:
            self.window.present()
            
    def close_window (self, button, data  = None):
        """ Function doc """
        #self.BackUpWindowData()
        self.window.destroy()
        self.visible    =  False
        dprint('self.visible',self.visible)
    
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
        dprint("Selected color: ", list(color))
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
            dprint(color)
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

    # ----------------------------------------------------------------------- #
    #  [NOVO] Dynamic bonds: cor unica opcional                                 #
    # ----------------------------------------------------------------------- #
    # ----------------------------------------------------------------------- #
    #  [NOVO] Dynamic bonds: cor unica opcional                                 #
    # ----------------------------------------------------------------------- #
    def set_dynamic_bonds_parameters(self):
        """Configura os controles de cor unica das dynamic bonds na abertura.

        Le as preferencias atuais e reflete nos widgets do Glade:
          - checkbox_dynamic_bonds_single_color (GtkCheckButton)
          - colorbutton_dynamic_bonds          (GtkColorButton)
        Se os widgets nao existirem no glade (ainda nao adicionados), sai sem erro.
        """
        try:
            from gi.repository import Gdk
            gp = self.vm_session.vm_config.gl_parameters

            self.checkbox_dynamic_bonds_single_color = self.builder.get_object(
                'checkbox_dbond_unique_color')
            self.colorbutton_dynamic_bonds = self.builder.get_object(
                'btn_dbond_unique_color')

            if self.checkbox_dynamic_bonds_single_color is not None:
                self.checkbox_dynamic_bonds_single_color.set_active(
                    gp.get("dynamic_bonds_single_color", False))
                self.checkbox_dynamic_bonds_single_color.connect(
                    'toggled', self.on_dynamic_bonds_single_color_toggled)

            if self.colorbutton_dynamic_bonds is not None:
                c = gp.get("dynamic_bonds_color", [1.0, 1.0, 1.0, 1.0])
                rgba = Gdk.RGBA()
                rgba.red, rgba.green, rgba.blue = c[0], c[1], c[2]
                rgba.alpha = c[3] if len(c) > 3 else 1.0
                self.colorbutton_dynamic_bonds.set_rgba(rgba)
                self.colorbutton_dynamic_bonds.connect(
                    'color-set', self.on_dynamic_bonds_color_set)
        except Exception as e:
            dprint("set_dynamic_bonds_parameters skipped:", e)

    def on_dynamic_bonds_single_color_toggled(self, widget):
        """Liga/desliga o uso de cor unica para as ligacoes dinamicas (QC).

        Conectar no Glade ao 'toggled' de um GtkCheckButton chamado
        'checkbox_dynamic_bonds_single_color'. Ao mudar, atualiza a preferencia
        e reconstroi a representacao 'dynamic' dos objetos para refletir na tela.
        """
        active = widget.get_active()
        self.vm_session.vm_config.gl_parameters["dynamic_bonds_single_color"] = active
        self._refresh_dynamic_bonds_representation()

    def on_dynamic_bonds_color_set(self, widget):
        """Define a cor unica das ligacoes dinamicas.

        Conectar no Glade ao 'color-set' de um GtkColorButton chamado
        'colorbutton_dynamic_bonds'. Guarda a cor (RGBA em [0,1]) e, se a opcao
        de cor unica estiver ligada, atualiza a representacao imediatamente.
        """
        color = list(widget.get_rgba())  # [r, g, b, a] em [0,1]
        self.vm_session.vm_config.gl_parameters["dynamic_bonds_color"] = color
        if self.vm_session.vm_config.gl_parameters.get("dynamic_bonds_single_color", False):
            self._refresh_dynamic_bonds_representation()

    def _refresh_dynamic_bonds_representation(self):
        """Reconstroi a representacao 'dynamic' de todos os objetos visuais.

        Zera a representacao 'dynamic' existente para que ela seja recriada
        (create_representation lera as preferencias novas e aplicara/limpara a
        cor unica), e redesenha a cena.
        """
        try:
            for vm_object in self.vm_session.vm_objects_dic.values():
                if "dynamic" in vm_object.representations and \
                   vm_object.representations["dynamic"] is not None:
                    # Aplica a mudanca na representacao ja existente, sem recriar:
                    rep = vm_object.representations["dynamic"]
                    gp = self.vm_session.vm_config.gl_parameters
                    if gp.get("dynamic_bonds_single_color", False):
                        rep.set_uniform_color(gp.get("dynamic_bonds_color", [1.0, 1.0, 1.0, 1.0]))
                    else:
                        rep.clear_uniform_color()
            self.vm_session.vm_glcore.queue_draw()
        except Exception as e:
            dprint("dynamic bonds color refresh failed:", e)

    def set_interface_startup_shutdown_paramters (self):
        """ Function doc """
        a = self.vm_session.vm_config.gl_parameters['autosave']      
        b = self.vm_session.vm_config.gl_parameters['askSaveUnsave'] 
        
        self.tmp_autosave            .set_active(a)
        self.ask_autosave_and_unsaved.set_active(b)
        
        # Criterio de autosave (timer em minutos + contador de eventos --
        # dispara pelo que vier primeiro. Ver MainWindow._restart_autosave_
        # timer e pDynamoSession.register_change_and_maybe_autosave).
        self.entry_autosave_interval    = self.builder.get_object('entry_autosave_interval')
        self.entry_autosave_event_count = self.builder.get_object('entry_autosave_event_count')
        interval = self.vm_session.vm_config.gl_parameters.get('autosave_interval_minutes', 5)
        count    = self.vm_session.vm_config.gl_parameters.get('autosave_event_count', 20)
        self.entry_autosave_interval   .set_text(str(interval))
        self.entry_autosave_event_count.set_text(str(count))
        
        # "Save window size" -- ja existia no glade, era lido em __apply_
        # interface_general_parameters mas o valor nunca era usado (nem
        # setado aqui na abertura). Ver gl_parameters['save_window_size'].
        self.checkbox_save_window_size = self.builder.get_object('checkbox_save_window_size')
        save_window_size = self.vm_session.vm_config.gl_parameters.get('save_window_size', True)
        self.checkbox_save_window_size.set_active(bool(save_window_size))
        pass
    
    def set_general_parameters (self):
        """ Function doc """
        #       background
        self.colorbutton = self.builder.get_object('btn_background_color')
        bgcolor = self.vm_session.vm_config.gl_parameters["background_color"]
        rgba = Gdk.RGBA(bgcolor[0], bgcolor[1], bgcolor[2])
        self.colorbutton.set_rgba(rgba)

        
        # mouse_rotation_sensibility
        self.entry_rot_sensibililty = self.builder.get_object('entry_rot_sensibililty')
        entry_rot_sensibililty = self.vm_session.vm_config.gl_parameters["mouse_rotation_sensibility"]
        self.entry_rot_sensibililty.set_text(str(entry_rot_sensibililty))
        
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

        #-------------------------------------------------------------------------------------
        #    Labels: single "scale with zoom" option for ALL glArea labels
        #-------------------------------------------------------------------------------------
        # Replaces what used to be three separate checkboxes (picking,
        # distance and atom labels each had their own). zoom_sensitivity
        # is a 0.0..1.0 float per VismolFont (see vismol_font.py); this
        # single checkbox only offers the two extremes and is applied to
        # every label font at once in __apply_viewer_selections_parameters().
        self.chk_labels_scale_with_zoom = self.builder.get_object('chk_labels_scale_with_zoom')
        if self.chk_labels_scale_with_zoom is not None:
            zoom_sensitivity = self.vm_session.vm_config.gl_parameters.get('labels_zoom_sensitivity', 1.0)
            self.chk_labels_scale_with_zoom.set_active(zoom_sensitivity >= 0.5)
        #-------------------------------------------------------------------------------------

    def _populate_font_combo(self, combo, spin, current_font, current_size):
        """ Fills a font-family GtkComboBoxText with the bundled .ttf
            fonts (see vismol/libgl/fonts/) and selects `current_font`,
            and sets `spin` (a GtkSpinButton) to `current_size`. Shared
            helper used for both the Picking/Distance labels font and the
            Atom Labels (index/charge/residue.../glArea) font, which are
            configured independently of each other.
        """
        if combo is not None:
            combo.remove_all()
            available_fonts = list_available_fonts()
            if current_font not in available_fonts:
                available_fonts = [current_font] + available_fonts
            active_index = 0
            for i, font_name in enumerate(available_fonts):
                # Nice display name: drop the .ttf extension
                display_name = os.path.splitext(font_name)[0]
                combo.append(font_name, display_name)
                if font_name == current_font:
                    active_index = i
            combo.set_active(active_index)

        if spin is not None:
            spin.set_value(current_size)

    def _iter_all_vm_objects_atoms(self):
        """ Yields every (vm_object, atom) pair currently loaded, across
            every molecule in the session. Used by the "Show for All
            Atoms" / "Hide All" atom-label buttons.
        """
        for vm_object in self.vm_session.vm_objects_dic.values():
            for atom in vm_object.atoms.values():
                yield vm_object, atom

    def _current_label_content(self):
        """ Reads the 'Label Content' combo (falls back to gl_parameters,
            then to 'name' if the widget isn't available for some reason).
        """
        combo = self.builder.get_object('combo_label_content')
        if combo is not None:
            content = combo.get_active_id()
            if content:
                return content
        return self.vm_session.vm_config.gl_parameters.get('label_content', 'name')

    def on_btn_show_all_atom_labels (self, widget):
        """ Turns the "Labels" representation (atom index/charge/residue
            name/.../chain, per the Label Content selector) ON immediately
            for every atom of every loaded molecule, using the currently
            selected content option. This talks directly to the
            representation's *own* font object
            (vm_object.representations['labels'].vm_font) -- the one
            actually used for drawing -- not vm_object.vm_font, which the
            "labels" representation does not read.
        """
        content = self._current_label_content()
        self.vm_session.vm_config.gl_parameters['label_content'] = content

        for vm_object in self.vm_session.vm_objects_dic.values():
            if len(vm_object.atoms) == 0:
                continue
            if vm_object.representations.get('labels') is None:
                vm_object.create_representation(rep_type='labels')
            rep = vm_object.representations['labels']

            indexes = []
            for atom in vm_object.atoms.values():
                atom.label_text = compute_atom_label_text(atom, content)
                atom.labels = True
                indexes.append(atom.atom_id)

            rep.define_new_indexes_to_vbo(indexes)
            rep.active = True

        self.vm_session.vm_glcore.queue_draw()

    def on_btn_hide_all_atom_labels (self, widget):
        """ Turns the "Labels" representation OFF for every atom of every
            loaded molecule.
        """
        for vm_object in self.vm_session.vm_objects_dic.values():
            for atom in vm_object.atoms.values():
                atom.labels = False
            rep = vm_object.representations.get('labels')
            if rep is not None:
                rep.active = False

        self.vm_session.vm_glcore.queue_draw()

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
        #    Picking / Distance Labels: SHARED font family, INDEPENDENT sizes
        #-------------------------------------------------------------------------------------
        # Picking labels (#1 #2 #3 #4) and distance labels use the same
        # font family (one combo box, 'combo_pk_dist_font_family') but
        # keep their own, independently adjustable font sizes -- picking
        # via 'label_font_size', distance via 'pk_dist_label_font_size'.
        # "Scale with zoom" used to be a per-label checkbox here; it's
        # now a single option for ALL labels in the glArea, on the
        # Viewer > General tab (see set_general_parameters()).
        gp = self.vm_session.vm_config.gl_parameters

        self.combo_pk_dist_font_family = self.builder.get_object('combo_pk_dist_font_family')
        self.spin_pk_dist_font_size    = self.builder.get_object('spin_pk_dist_font_size')
        self._populate_font_combo(self.combo_pk_dist_font_family,
                                   self.spin_pk_dist_font_size,
                                   gp.get('label_font_file', DEFAULT_FONT_FILE),
                                   gp.get('label_font_size', DEFAULT_FONT_SIZE))

        self.spin_dist_label_font_size = self.builder.get_object('spin_dist_label_font_size')
        self._populate_font_combo(None,
                                   self.spin_dist_label_font_size,
                                   None,
                                   gp.get('pk_dist_label_font_size', gp.get('label_font_size', DEFAULT_FONT_SIZE)))
        #-------------------------------------------------------------------------------------

        #-------------------------------------------------------------------------------------
        #    Atom Labels (glArea): content (index/charge/resname/...) + family + size
        #-------------------------------------------------------------------------------------
        self.combo_atom_label_font_family = self.builder.get_object('combo_atom_label_font_family')
        self.spin_atom_label_font_size    = self.builder.get_object('spin_atom_label_font_size')
        self._populate_font_combo(self.combo_atom_label_font_family,
                                   self.spin_atom_label_font_size,
                                   gp.get('atom_label_font_file', DEFAULT_FONT_FILE),
                                   gp.get('atom_label_font_size', DEFAULT_FONT_SIZE))

        self.combo_label_content = self.builder.get_object('combo_label_content')
        if self.combo_label_content is not None:
            self.combo_label_content.remove_all()
            # (internal value, display label) -- internal value is what
            # gets stored in gl_parameters['label_content'] and passed to
            # representations.compute_atom_label_text().
            content_options = [
                ('name',          'Atom Name'),
                ('symbol',        'Atom Symbol'),
                ('index',         'Atom Index'),
                ('mm_charge',     'MM Charge'),
                ('residue_name',  'Residue Name'),
                ('residue_index', 'Residue Index'),
                ('chain',         'Chain'),
            ]
            current_content = gp.get('label_content', 'name')
            active_index = 0
            for i, (value, display) in enumerate(content_options):
                self.combo_label_content.append(value, display)
                if value == current_content:
                    active_index = i
            self.combo_label_content.set_active(active_index)

        self.btn_show_all_atom_labels = self.builder.get_object('btn_show_all_atom_labels')
        if self.btn_show_all_atom_labels is not None:
            self.btn_show_all_atom_labels.connect('clicked', self.on_btn_show_all_atom_labels)

        self.btn_hide_all_atom_labels = self.builder.get_object('btn_hide_all_atom_labels')
        if self.btn_hide_all_atom_labels is not None:
            self.btn_hide_all_atom_labels.connect('clicked', self.on_btn_hide_all_atom_labels)
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
        lines_with           = self.vm_session.vm_config.gl_parameters['line_width']              
        line_width_selection = self.vm_session.vm_config.gl_parameters['line_width_selection']
        line_type            = self.vm_session.vm_config.gl_parameters['line_type']
        
        self.entry_lines_with            =  self.builder.get_object('entry_lines_with')
        self.entry_lines_with_selections =  self.builder.get_object('entry_lines_with_selections')
        self.entry_lines_type            =  self.builder.get_object('entry_lines_type')

        self.entry_lines_with           .set_text(str(lines_with))   
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
        
        # Liga/desliga a representacao de ligacoes duplas/triplas nos sticks
        # (gl_parameters['multiple_bonds'] -- ver SticksRepresentation._get_
        # bond_order_per_bond em representations.py). Quando desligado, todas
        # as ligacoes sao desenhadas como simples, independente da ordem
        # percebida.
        self.checkbox_multiple_bonds = self.builder.get_object('checkbox_multiple_bonds')
        multiple_bonds = self.vm_session.vm_config.gl_parameters.get('multiple_bonds', True)
        self.checkbox_multiple_bonds.set_active(bool(multiple_bonds))

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
        self.set_paths_and_folders_in_parameters()
        parameters = self.vm_session.vm_config.gl_parameters      
        self.builder.get_object('entry_tmp_path').set_text(parameters['tmp_path'])
        self.builder.get_object('entry_workspace_path').set_text(parameters['workspace_path'])
        self.builder.get_object('entry_startup_path').set_text(parameters['startup_path'])

    def set_paths_and_folders_in_parameters (self):
        """ Function doc """
        parameters = self.vm_session.vm_config.gl_parameters      
        if 'tmp_path' in parameters.keys():         
            if os.path.isdir(parameters['tmp_path']):
                pass
            else:
                dprint('Folder not found:', parameters['tmp_path'])
                PDYNAMO3_SCRATCH = os.environ.get('PDYNAMO3_SCRATCH')
                parameters['tmp_path'] = PDYNAMO3_SCRATCH
        else:            
            PDYNAMO3_SCRATCH = os.environ.get('PDYNAMO3_SCRATCH')
            parameters['tmp_path'] = PDYNAMO3_SCRATCH
        
        
        if 'workspace_path'  in parameters.keys():
            if os.path.isdir(parameters['workspace_path']):
                pass
            else:
                dprint('Folder not found:', parameters['workspace_path'])
                workspace_path = os.path.join(self.home, 'workspace')
                parameters['workspace_path'] = workspace_path
        else:
            workspace_path = os.path.join(self.home, 'workspace')
            parameters['workspace_path'] = workspace_path
       
        
        if 'startup_path' in parameters.keys():
            if os.path.isdir(parameters['startup_path']):
                pass
            else:
                dprint('Folder not found:', parameters['startup_path'])
                parameters['startup_path'] = self.home
        else:
            parameters['startup_path'] = self.home

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
        
        # Liga/desliga ligacoes duplas/triplas nos sticks. Nao precisa
        # reconstruir nenhuma representacao: SticksRepresentation._get_bond_
        # order_per_bond le' este valor de gl_parameters a cada frame (ver
        # _refresh_bond_order_tbo), entao o efeito aparece no proximo redraw
        # (ja disparado por __apply_all_changes logo apos todos os
        # __apply_*_parameters, via vm_glcore.queue_draw()).
        self.gl_parameters['multiple_bonds'] = self.checkbox_multiple_bonds.get_active()

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

        #-----------------------------------------------------------------------
        #    Labels: single "scale with zoom" option for ALL glArea labels
        #-----------------------------------------------------------------------
        # Read once here (Viewer > General tab) and applied below to every
        # label font: picking, distance and atom labels. zoom_sensitivity
        # is just a uniform value read every frame by the shader, so it
        # can be set directly on each VismolFont with no VAO rebuild.
        chk_scale_zoom = self.builder.get_object('chk_labels_scale_with_zoom')
        labels_zoom_sensitivity = 1.0 if (chk_scale_zoom is not None and chk_scale_zoom.get_active()) else 0.0
        self.vm_session.vm_config.gl_parameters['labels_zoom_sensitivity'] = labels_zoom_sensitivity
        #-----------------------------------------------------------------------

        #-----------------------------------------------------------------------
        #    Picking / Distance Labels: SHARED font family, INDEPENDENT sizes
        #-----------------------------------------------------------------------
        # Picking labels (#1 #2 #3 #4, on vm_glcore.vm_font/vm_font_static)
        # and distance labels (vm_glcore.vm_font_dist) share ONE font
        # family (a single combo box) but keep their own font sizes.
        # apply_settings() marks each font's VAO for regeneration
        # (vao=None) so the new font/size take effect on the next draw,
        # without needing to reload the molecules.
        combo_font = self.builder.get_object('combo_pk_dist_font_family')
        spin_pk_size   = self.builder.get_object('spin_pk_dist_font_size')
        spin_dist_size = self.builder.get_object('spin_dist_label_font_size')
        if combo_font is not None and spin_pk_size is not None:
            font_file = combo_font.get_active_id() or DEFAULT_FONT_FILE
            pk_size   = spin_pk_size.get_value()
            dist_size = spin_dist_size.get_value() if spin_dist_size is not None else pk_size

            self.vm_session.vm_config.gl_parameters['label_font_file']         = font_file
            self.vm_session.vm_config.gl_parameters['label_font_size']         = pk_size
            self.vm_session.vm_config.gl_parameters['pk_dist_label_font_file'] = font_file
            self.vm_session.vm_config.gl_parameters['pk_dist_label_font_size'] = dist_size

            vm_glcore = self.vm_session.vm_glcore
            for font_obj in (vm_glcore.vm_font, vm_glcore.vm_font_static):
                font_obj.apply_settings(font_file=font_file, size=pk_size)
                font_obj.zoom_sensitivity = labels_zoom_sensitivity
            vm_glcore.vm_font_dist.apply_settings(font_file=font_file, size=dist_size)
            vm_glcore.vm_font_dist.zoom_sensitivity = labels_zoom_sensitivity
        #-----------------------------------------------------------------------

        #-----------------------------------------------------------------------
        #    Atom Labels (glArea): content (index/charge/resname/...) + family + size
        #-----------------------------------------------------------------------
        # This is a SEPARATE font from the one above: the "labels"
        # representation (atom index/MM charge/residue name/residue
        # index/chain -- also settable per-selection via the glArea
        # right-click "Show" menu) owns its OWN VismolFont instance, one
        # per vm_object, at vm_object.representations['labels'].vm_font --
        # NOT vm_object.vm_font, which that representation never reads.
        combo_font = self.builder.get_object('combo_atom_label_font_family')
        spin_size  = self.builder.get_object('spin_atom_label_font_size')
        combo_content = self.builder.get_object('combo_label_content')
        if combo_font is not None and spin_size is not None:
            font_file = combo_font.get_active_id() or DEFAULT_FONT_FILE
            font_size = spin_size.get_value()
            content = (combo_content.get_active_id() if combo_content is not None else None) or 'name'

            self.vm_session.vm_config.gl_parameters['atom_label_font_file'] = font_file
            self.vm_session.vm_config.gl_parameters['atom_label_font_size'] = font_size
            self.vm_session.vm_config.gl_parameters['label_content']       = content

            for vm_object in self.vm_session.vm_objects_dic.values():
                rep = vm_object.representations.get('labels')
                if rep is None:
                    # Representation not created yet for this object (no
                    # atom has ever been labeled): nothing to refresh.
                    continue
                rep.vm_font.apply_settings(font_file=font_file, size=font_size)
                rep.vm_font.zoom_sensitivity = labels_zoom_sensitivity
                # Re-derive the text for every atom currently labeled
                # (atom.labels == True) so a Label Content change (e.g.
                # "Atom Name" -> "MM Charge") is reflected immediately on
                # Apply, without needing to re-pick the selection.
                relabeled_indexes = []
                for atom in vm_object.atoms.values():
                    if getattr(atom, 'labels', False):
                        atom.label_text = compute_atom_label_text(atom, content)
                        relabeled_indexes.append(atom.atom_id)
                if relabeled_indexes:
                    rep.define_new_indexes_to_vbo(relabeled_indexes)
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
        self.btn_seqviewer_bg_color = self.builder.get_object('btn_seqviewer_bg_color')
        color = self.btn_seqviewer_bg_color.get_rgba()
        color = list(color)
        self.vm_session.main.bottom_notebook.seqview.text_drawing_area.bg_color = color

        self.btn_seqviewer_label_color = self.builder.get_object('btn_seqviewer_label_color')
        color = self.btn_seqviewer_label_color.get_rgba()
        color = list(color)
        self.vm_session.main.bottom_notebook.seqview.text_drawing_area.marker_color = color
        #---------------------------------------------------------------
        
        
        
        #---------------------------------------------------------------
        scroll_step      = float(self.builder.get_object('entry_scroll_step')   .get_text() )
        sleep_time_coc   = float(self.builder.get_object('entry_sleep_time_coc').get_text() )
        field_of_view    = float(self.builder.get_object('entry_field_of_view') .get_text() )
        rot_sensibililty = float(self.builder.get_object('entry_rot_sensibililty') .get_text() )
        
        self.gl_parameters["scroll_step"]                = scroll_step    
        self.gl_parameters["center_on_coord_sleep_time"] = sleep_time_coc 
        self.gl_parameters["field_of_view"]              = field_of_view  
        self.gl_parameters["mouse_rotation_sensibility"] = rot_sensibililty  
        #---------------------------------------------------------------
    def __apply_interface_general_parameters (self):
        """ Function doc """
        
        #---------------------------------------------------------------
        #                       Startup Path
        #---------------------------------------------------------------
        path = self.entry_startup_path.get_text()
        if os.path.isdir(path):
            self.vm_session.vm_config.gl_parameters['startup_path'] = path
            dprint('Defining New Startup Path:', path)

        else:
            dialog = Gtk.MessageDialog(
                                flags=0,
                                message_type=Gtk.MessageType.ERROR,
                                buttons=Gtk.ButtonsType.OK,
                                text="ERROR: Startup folder not found.",
                            )

            dialog.format_secondary_text(
                                    "If the desired path is correct, create the folder using your file manager."
                                        )
            dialog.run()
            dialog.destroy()
        #---------------------------------------------------------------
        
        
        #---------------------------------------------------------------
        #                       Workspace 
        #---------------------------------------------------------------
        workspace_path = self.builder.get_object('entry_workspace_path').get_text()
        if os.path.isdir(workspace_path):
            self.vm_session.vm_config.gl_parameters['workspace_path'] = workspace_path
            dprint('Defining workspace path:', workspace_path)
        else:
            dialog = Gtk.MessageDialog(
                                flags=0,
                                message_type=Gtk.MessageType.ERROR,
                                buttons=Gtk.ButtonsType.OK,
                                text="ERROR: Workspace folder not found.",
                            )
            dialog.format_secondary_text(
                                    "If the desired path is correct, create the folder using your file manager."
                                        )
            dialog.run()
            dialog.destroy()
        #---------------------------------------------------------------

        #---------------------------------------------------------------
        #                       Temp 
        #---------------------------------------------------------------
        tmp_path       = self.builder.get_object('entry_tmp_path').get_text()
        if os.path.isdir(tmp_path):
            self.vm_session.vm_config.gl_parameters['tmp_path'] = tmp_path
            dprint('Defining temporary path:', tmp_path)
        else:
            dialog = Gtk.MessageDialog(
                                flags=0,
                                message_type=Gtk.MessageType.ERROR,
                                buttons=Gtk.ButtonsType.OK,
                                text="ERROR: Temporary folder not found.",
                            )
            dialog.format_secondary_text(
                                    "If the desired path is correct, create the folder using your file manager."
                                        )
            dialog.run()
            dialog.destroy()
        #---------------------------------------------------------------
        
        

        a = self.builder.get_object('checkbox_tmp_autosave').get_active()
        b = self.builder.get_object('checkbox_ask_autosave_and_unsaved').get_active()
        c = self.builder.get_object('checkbox_save_window_size').get_active()

        self.vm_session.vm_config.gl_parameters['autosave']      = a
        self.vm_session.vm_config.gl_parameters['askSaveUnsave'] = b
        # [BUG FIX] 'c' era lido mas nunca usado -- o checkbox "Save window
        # size" nao tinha efeito nenhum. Agora liga/desliga a persistencia
        # de dimensoes da janela (ver MainWindow.__init__/window_resize/
        # on_delete_event).
        self.vm_session.vm_config.gl_parameters['save_window_size'] = c
        dprint(a,b,c)
        
        # Criterio de autosave (timer em minutos + contador de eventos).
        # Aceita virgula OU ponto decimal; cai pro valor atual em caso de
        # entrada invalida, em vez de quebrar o Apply inteiro.
        try:
            interval = float(self.entry_autosave_interval.get_text().replace(',', '.'))
            self.vm_session.vm_config.gl_parameters['autosave_interval_minutes'] = interval
        except ValueError:
            dprint('Invalid autosave interval, keeping previous value.')
        try:
            count = int(float(self.entry_autosave_event_count.get_text().replace(',', '.')))
            self.vm_session.vm_config.gl_parameters['autosave_event_count'] = count
        except ValueError:
            dprint('Invalid autosave event count, keeping previous value.')
        
        # O intervalo pode ter mudado: reinicia o timer periodico de
        # autosave da janela principal com o novo valor. self.main_session
        # e' a propria instancia de MainWindow (ver EasyHybridPreferencesWindow.
        # __init__: self.main_session = main).
        restart_timer = getattr(self.main_session, '_restart_autosave_timer', None)
        if restart_timer is not None:
            restart_timer()
        
        
        '''
        a = self.builder.get_object('checkbox_output_tag').get_active()
        b = self.builder.get_object('checkbox_output_sys_name').get_active()
        c = self.builder.get_object('checkbox_output_ff_model').get_active()        
        d = self.builder.get_object('checkbox_output_qc_model').get_active()
        e = self.builder.get_object('checkbox_output_qc_charge').get_active()
        f = self.builder.get_object('checkbox_output_qc_multiplicity').get_active()
        g = self.builder.get_object('checkbox_output_qc_size').get_active()
        h = self.builder.get_object('checkbox_output_simtype').get_active()
        
        self.vm_session.vm_config.gl_parameters['fname_output_tag']             = a
        self.vm_session.vm_config.gl_parameters['fname_output_sys_name']        = b
        self.vm_session.vm_config.gl_parameters['fname_output_ff_model']        = c
        self.vm_session.vm_config.gl_parameters['fname_output_qc_model']        = d
        self.vm_session.vm_config.gl_parameters['fname_output_qc_charge']       = e
        self.vm_session.vm_config.gl_parameters['fname_output_qc_multiplicity'] = f
        self.vm_session.vm_config.gl_parameters['fname_output_qc_size']         = g
        self.vm_session.vm_config.gl_parameters['fname_output_simtype']         = h
        '''
        
        
        
        
        #self.vm_session.vm_config.gl_parameters['startup_path'] =
        #pprint(self.vm_session.vm_config.gl_parameters)
        #startup_path   = self.builder.get_object('entry_startup_path').get_text()
        #workspace_path = self.builder.get_object('entry_workspace_path').get_text()
        #tmp_path       = self.builder.get_object('entry_tmp_path').get_text()
        #
        #self.vm_session.vm_config.gl_parameters['startup_path'  ] = startup_path  
        #self.vm_session.vm_config.gl_parameters['workspace_path'] = workspace_path
        #self.vm_session.vm_config.gl_parameters['tmp_path'      ] = tmp_path      






















