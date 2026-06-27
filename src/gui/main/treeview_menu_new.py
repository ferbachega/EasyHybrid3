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
import os, sys, time
import gi 
import signal
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk, Pango
from gi.repository import GdkPixbuf


from gui.widgets.custom_widgets  import VismolSelectionTypeBox
from gui.widgets.custom_widgets  import FileChooser
from gui.widgets.custom_widgets  import get_colorful_square_pixel_buffer
from gui.widgets.custom_widgets  import ReactionCoordinateBox
from gui.widgets.custom_widgets  import SequenceViewerBox

from gui.windows.setup.windows_and_dialogs import ImportANewSystemWindow
from gui.windows.setup.windows_and_dialogs import EasyHybridDialogSetQCAtoms
from gui.windows.setup.windows_and_dialogs import EasyHybridSetupQCModelWindow
from gui.windows.setup.windows_and_dialogs import EasyHybridGoToAtomWindow
#from gui.windows.setup.windows_and_dialogs import PDynamoSelectionWindow
from gui.windows.setup.windows_and_dialogs import EasyHybridSelectionWindow
from gui.windows.setup.windows_and_dialogs import ExportDataWindow
from gui.windows.setup.windows_and_dialogs import EasyHybridDialogPrune
from gui.windows.setup.windows_and_dialogs import MakeSolventBoxWindow


from gui.windows.setup.windows_and_dialogs import ImportTrajectoryWindow
from gui.windows.setup.windows_and_dialogs import TrajectoryPlayerWindow
from gui.windows.setup.windows_and_dialogs import InfoWindow
from gui.windows.setup.windows_and_dialogs import MergeSystemWindow
from gui.windows.setup.windows_and_dialogs import SolvateSystemWindow
from gui.windows.setup.windows_and_dialogs import SimpleDialog
from gui.windows.setup.edit_frames_dialog import EditFrameDialog

from gui.windows.setup.easyhybrid_terminal    import TerminalWindow
from gui.windows.setup.selection_list_window  import *
from gui.windows.setup.setup_interface        import EasyHybridPreferencesWindow
from gui.windows.setup.process_manager_window import ProcessManagerWindow

from gui.windows.simulation.single_point_window          import SinglePointWindow
from gui.windows.simulation.geometry_optimization_window import GeometryOptimization
from gui.windows.simulation.PES_scan_window              import PotentialEnergyScanWindow 
from gui.windows.simulation.PES_advanced_scan_window     import AdvancedPotentialEnergyScanWindow 
from gui.windows.simulation.molecular_dynamics_window    import MolecularDynamicsWindow 
from gui.windows.simulation.umbrella_sampling_window     import UmbrellaSamplingWindow 
from gui.windows.simulation.chain_of_states_opt_window   import ChainOfStatesOptWindow 
from gui.windows.simulation.normal_modes_window          import NormalModesWindow 

from gui.windows.analysis.WHAM_analysis_window                    import WHAMWindow 
from gui.windows.analysis.normal_modes_analysis_window            import NormalModesAnalysisWindow 
from gui.windows.analysis.surface_analysis_window                 import SurfaceAnalysisWindow 
#from gui.windows.analysis.surface_list_window                     import SurfaceListWindow 
from gui.windows.analysis.energy_refinement_window                import EnergyRefinementWindow
from gui.windows.analysis.PES_analysis_window                     import PotentialEnergyAnalysisWindow
from gui.windows.analysis.distance_angle_dihedral_analysis_window import DistanceAngleDihedralAnalysisWindow
from gui.windows.analysis.RMSD_tool                               import RMSDToolWindow
from gui.windows.analysis.RMSD_analysis_window                    import RMSDAnalysisWindow #/home/fernando/programs/EasyHybrid3/src/gui/windows/analysis/RMSD_analysis_window.py
from gui.windows.analysis.align_trajectory                        import AlignTrajectoryWindow #/home/fernando/programs/EasyHybrid3/src/gui/windows/analysis/RMSD_analysis_window.py
from gui.windows.analysis.reimaging_trajectory                    import ReimagingTrajectoryWindow #/home/fernando/programs/EasyHybrid3/src/gui/windows/analysis/RMSD_analysis_window.py

from util.geometric_analysis import get_simple_distance
from util.sequence_plot import GtkSequenceViewer
from util.rama_plot import RamachandranWindow


from pdynamo.pDynamo2EasyHybrid import pDynamoSession
import numpy as np



from pCore                     import Align                                        , \
                                      Clone                                        , \
                                      logFile                                      , \
                                      Selection                                    , \
                                      TestScript_InputDataPath                     , \
                                      TestScript_OutputDataPath                    , \
                                      XHTMLLogFileWriter

# --- imports entre modulos adicionados na refatoracao ---
from gui.main.preferences_window import PreferencesWindow

class TreeViewMenu_new:
    """ Class doc """
    
    def __init__ (self, treeview):
        """ Class initialiser """
        self.treeview = treeview
        self.main     = treeview.main 
        self.filechooser   = FileChooser()
        self.rename_window_visible =  False
        
        
        new_vobject_menu_items = {
                                'header'                : ['MenuItem', None]   ,
                                '_separator'            : ['separator', None],
                                'Rename'                : ['MenuItem',self._menu_rename] ,
                                '_separator'            : ['separator', None],
                                'Show / Hide Cell'      : ['MenuItem',self.show_or_hide_cell],

                                'Edit Frames'           : ['MenuItem',self.call_editframe_window],
                                'Go To Atom'            : ['MenuItem',self._menu_go_to_atom] ,
                                '_separator'            : ['separator', None],

                                'Export As...'          : ['MenuItem',self._menu_export_data_window] ,
                                '_separator'            : ['separator', None],
                                'Delete'                : ['MenuItem',self._menu_delete_vm_object] ,
                                }

        new_system_menu_items = {
                                
                                'header'                  : None                            ,
                                                          
                                '_separator'              : ''                              ,
                                                          
                                'Info'                    : self._show_info                  ,
                                                          
                                '_separator'              : ''                              ,
                                                          
                                'Rename'                  : self._menu_rename               ,
                                'Import Data...'          : self._menu_load_data_to_system  ,
                                'Reference Color'         : self._menu_change_color_palette ,
                                #'Edit Parameters'        : self.f2                         ,
                                'Export As...'            : self._menu_export_data_window    ,
                                
                                '_separator'              : ''                              ,
                                'Merge With...'           : self._menu_merge_system         ,
                                'Clone System'            : self._menu_clone_system         ,
                                
                                '_separator'              : ''                              ,
                                
                                'Delete'                  : self._menu_delete_system        ,
                                #'test'  : self.f1 ,
                                #'f1'    : self.f1 ,
                                #'f2'    : self.f2 ,
                                #'gordão': self.f3 ,
                                #'delete': self.f3 ,
                                }
                    

        menu_items = {
            "Abrir": self._menu_rename,
            "Salvar": self._menu_rename,
            "_separator": None,

            "Exportar": {
                "Como PNG": self._menu_rename,
                "Como JPG": self._menu_rename,
                "_separator": None,
                "Avançado": {
                    "Alta Qualidade": self._menu_rename,
                    "Baixa Qualidade": self._menu_rename,
                }
            },

            "header": None
        }


        vobject_menu_items = {
                                'header'                : None    ,

                                '_separator'            : ''      ,
                                #'Info'                  : self.f2 ,
                                'Rename'                : self._menu_rename ,
                                '_separator'            : ''      ,
                                'Show / Hide Cell'      : self.show_or_hide_cell,
                                
                                
                                'Edit Frames'           : self.call_editframe_window,
                                #'Load Data Into System' : self._menu_load_data_to_system,
                                #'Change Color'  : self.change_system_color_palette ,
                                'Go To Atom'            : self._menu_go_to_atom ,
                                '_separator'            : ''      ,

                                'Export As...'          : self._menu_export_data_window ,
                                #'Merge System With...'  : self.f3 ,
                                '_separator'            : ''      ,

                                'Delete'                : self._menu_delete_vm_object ,
                                #'test'  : self.f1 ,
                                #'f1'    : self.f1 ,
                                #'f2'    : self.f2 ,
                                #'gordão': self.f3 ,
                                #'delete': self.f3 ,
                                }

        system_menu_items = {
                                
                                'header'                  : None                            ,
                                                          
                                '_separator'              : ''                              ,
                                                          
                                'Info'                    : self._show_info                  ,
                                                          
                                '_separator'              : ''                              ,
                                                          
                                'Rename'                  : self._menu_rename               ,
                                'Import Data...'          : self._menu_load_data_to_system  ,
                                'Reference Color'         : self._menu_change_color_palette ,
                                #'Edit Parameters'        : self.f2                         ,
                                'Export As...'            : self._menu_export_data_window    ,
                                
                                '_separator'              : ''                              ,
                                'Merge With...'           : self._menu_merge_system         ,
                                'Clone System'            : self._menu_clone_system         ,
                                
                                '_separator'              : ''                              ,
                                
                                'Delete'                  : self._menu_delete_system        ,
                                #'test'  : self.f1 ,
                                #'f1'    : self.f1 ,
                                #'f2'    : self.f2 ,
                                #'gordão': self.f3 ,
                                #'delete': self.f3 ,
                                }
                    
                    
                    
        #self.tree_view_vobj_menu  , self.tree_header_vobj_menu    = self._build_treemenu_from_dicts(new_vobject_menu_items)
        self.tree_view_vobj_menu  , self.tree_header_vobj_menu    = self.build_tree_view_menu(vobject_menu_items)
        self.tree_view_sys_menu   , self.tree_header_sys_menu     = self.build_tree_view_menu(system_menu_items)

    def show_or_hide_cell (self, widget):
        """ Function doc """
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()

        old_name = model.get_value(iter, 2)
        v_id     = model.get_value(iter, 1)
        e_id     = model.get_value(iter, 0)
        #----------------------------------------------------------------------        
        #system = self.main.p_session.psystem[e_id]
        
        vobject = self.main.vm_session.vm_objects_dic[v_id]
        if "cell_lines" in vobject.representations.keys():
            if vobject.representations["cell_lines"].active:
                self.main.vm_session.hide_cell (vobject)
            
            else:
                self.main.vm_session.show_cell (vobject)
        else:
            self.main.vm_session.show_cell (vobject)
            
    def call_editframe_window (self, widget):
        """ Function doc """
        selection        = self.treeview.get_selection()
        (model, iter)    = selection.get_selected()
        e_id             = int(model.get_value(iter, 0))
        vm_object_index  = int(model.get_value(iter, 1))
        self.main.edit_frames_dialog.open_window (vm_object_index)
        
    def _show_info (self, widget):
        """ Function doc """
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        e_id          = int(model.get_value(iter, 0))  # @+
        #----------------------------------------------------------------------        
        system = self.main.p_session.psystem[e_id]
        window = InfoWindow(system)
        
    def _menu_export_data_window (self,vobject = None ):
        """ Function doc """
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        e_id          = int(model.get_value(iter, 0)) 
        #print(e_id)
        self.treeview.main.export_data_window.open_window(sys_selected = e_id)
        
    def _menu_load_data_to_system (self, vobject = None ):
        """ Function doc """
        selection        = self.treeview.get_selection()
        model, iter      = selection.get_selected()
        #print (list(model))
        self.main.import_trajectory_window.open_window(sys_selected = model.get_value(iter, 0))

    def _menu_change_color_palette (self, widget):
        """ Function doc """
        #selection               = self.selections[self.current_selection]
        self.colorchooserdialog = Gtk.ColorChooserDialog()
        
        if self.colorchooserdialog.run() == Gtk.ResponseType.OK:
            color = self.colorchooserdialog.get_rgba()
            #print(color.red,color.green, color.blue )
            new_color = [color.red, color.green, color.blue]

        else:
            new_color = False
        
        
    
        if new_color:
            self.colorchooserdialog.destroy()

            #----------------------------------------------------------------------
            selection     = self.treeview.get_selection()
            (model, iter) = selection.get_selected()
            self.selectedID  = int(model.get_value(iter, 0))  # @+
            #----------------------------------------------------------------------
            
            system = self.main.p_session.psystem[self.selectedID]
            
            self.main.change_reference_color(system, new_color)

    def _menu_merge_system (self, widget):
        """ Function doc """
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        e_id          = int(model.get_value(iter, 0)) 
        self.main.merge_system_window.selected_system_id = e_id
        self.main.merge_system_window.open_window()
    
    def _menu_clone_system (self, widget):
        """ Function doc """
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        e_id          = int(model.get_value(iter, 0))
        name          = model.get_value(iter, 2) 
        #tag           = self.main.p_session.psystem[e_id].tag
        #--------------------------------------------------------------                                                            
        dialog = EasyHybridDialogPrune(main         = self.main,    
                                       num_of_atoms = 'all'    ,    
                                       name         = name     ,    
                                       tag          = 'UNK'    ,    
                                       e_id         = e_id     ,    
                                       _type        = 1        )    
        name         = dialog.name        
        tag          = dialog.tag  
        color        = dialog.color 
        vobject_id   = dialog.vobject_id
        #--------------------------------------------------------------                                                            

        
        vobject      = self.main.vm_session.vm_objects_dic[vobject_id]
        self.main.p_session.set_psystem_coordinates_from_vobject(vobject   = vobject,
                                                                           system_id =  e_id  )
        #print(e_id)
        self.main.p_session.clone_system( e_id    = e_id, 
                                          vobject = vobject, 
                                          name    = name, 
                                          tag     = tag, 
                                          color   = color)
        self._save_backup_file()
    
    def _menu_go_to_atom (self, vobject = None):
        """ Function doc """
        ##print('f2')
        #self._show_lines(vobject = self.vobjects[0], indices = [0,1,2,3,4] )
        self.treeview.main.go_to_atom_window.open_window()
        #self.treeview.vm_session.go_to_atom_window.open_window()
    def f3 (self, vobject = None):
        """ Function doc """
        
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()


        self.selectedID  = int(model.get_value(iter, 1))  # @+
        
        
        
        del self.treeview.vm_session.vobjects_dic[self.selectedID]
        '''
        vobject = self.treeview.vm_session.vobjects_dic.pop(self.selectedID)
        del vobject
        '''
        self.treeview.store.clear()
        for vobj_index ,vis_object in self.treeview.vm_session.vobjects_dic.items():
            data = [vis_object.active          , 
                    str(vobj_index),
                    vis_object.name            , 
                    str(len(vis_object.atoms)) , 
                    str(len(vis_object.frames)),
                   ]
            model.append(data)
        self.treeview.vm_session.glwidget.queue_draw()

    def _menu_rename (self, menu_item = None ):
        """  
        menu_item = Gtk.MenuItem object at 0x7fbdcc035700 (GtkMenuItem at 0x37cf6c0)
        
        """
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()

        old_name = model.get_value(iter, 2)
        v_id     = model.get_value(iter, 1)
        e_id     = model.get_value(iter, 0)
        tag      = self.main.p_session.psystem[e_id].e_tag 
        
        old_name = old_name.split("- ")
        old_name = old_name[-1]
        
        if self.rename_window_visible:
            self.preferences.set_names (old_name, tag)
            pass
        else:
            
            self.preferences = PreferencesWindow(main = self.main, 
                                                 e_id = e_id     ,
                                                 v_id = v_id     )
            self.preferences.set_names (old_name, tag)
        self._save_backup_file()
        
    def destroy (self, widget):
        """ Function doc """
        self.rename_window_visible = False
    
    def _menu_delete_vm_object (self, widget):
        """ Function doc """
        self.main.delete_vm_object ( vm_object_index = self.vobject_index)
        self._save_backup_file()
    
    def _menu_delete_system (self, widget):
        """ Function doc """
        self.main.delete_system (system_e_id = self.system_e_id )
        self._save_backup_file()
 
    def build_tree_view_menu_old(self, menu_items = None):
        """ Function doc """
        tree_view_menu = Gtk.Menu()
        menu_header    = None
        
        for label in menu_items:
            
            
            if menu_items[label] == None:
                # just a label
                
                mitem = Gtk.MenuItem(label = label)
                if label == 'header':
                    menu_header    = mitem
                
                
            elif  label == '_separator':
                mitem = Gtk.SeparatorMenuItem()
            
            else:
                mitem = Gtk.MenuItem(label = label)
                mitem.connect('activate', menu_items[label])
            
            tree_view_menu.append(mitem)
            #mitem = Gtk.SeparatorMenuItem()
            #self.tree_view_menu.append(mitem)

        tree_view_menu.show_all()
        return tree_view_menu, menu_header

    def build_tree_view_menu (self, menu_items):
        """Cria menus e submenus a partir de um dicionário."""
        menu = Gtk.Menu()
        menu_header = None

        for label, value in menu_items.items():

            # --- Separador ---
            if label == "_separator":
                mitem = Gtk.SeparatorMenuItem()

            # --- Header (item destacado) ---
            elif label == "header":
                mitem = Gtk.MenuItem(label=label)
                mitem.set_sensitive(False)      # desabilita
                menu_header = mitem

            # --- Submenu (value é um dicionário) ---
            elif isinstance(value, dict):
                mitem = Gtk.MenuItem(label=label)
                # cria o submenu recursivamente
                submenu, _ = self.build_tree_view_menu(value)
                mitem.set_submenu(submenu)

            # --- Item simples sem callback ---
            elif value is None:
                mitem = Gtk.MenuItem(label=label)

            # --- Item com callback ---
            else:
                mitem = Gtk.MenuItem(label=label)
                mitem.connect("activate", value)

            menu.append(mitem)

        menu.show_all()
        return menu, menu_header

        def _build_submenus_from_dicts(self, menu_dict):
            """ Function doc """
            menu = Gtk.Menu()
            for key in menu_dict:
                mitem = Gtk.MenuItem(key)
                
                
                if menu_dict[key][0] == "submenu":
                    menu2 = self._build_submenus_from_dicts(menu_dict[key][1])
                    mitem.set_submenu(menu2)
                
                elif menu_dict[key][0] == "separator":
                    mitem = Gtk.SeparatorMenuItem()
                
                else:
                    if menu_dict[key][1] != None:
                        mitem.connect("activate", menu_dict[key][1])
                    else:
                        pass
                menu.append(mitem)
            return menu
        
        def _build_treemenu_from_dicts(self, menu_dict):
            """ Function doc """
            tree_view_menu = Gtk.Menu()
            menu_header    = None
            
            for key in menu_dict:
                mitem = Gtk.MenuItem(label=key)
                
                if menu_dict[key][0] == "submenu":
                    menu2 = self._build_submenus_from_dicts(menu_dict[key][1])
                    mitem.set_submenu(menu2)
                
                elif key == 'header':
                        menu_header    = mitem
                
                elif menu_dict[key][0] == "separator":
                    mitem = Gtk.SeparatorMenuItem()
                else:
                    if menu_dict[key][1] != None:
                        mitem.connect("activate", menu_dict[key][1])
                    else:
                        pass
                tree_view_menu.append(mitem)
            return tree_view_menu, menu_header


        def open_menu (self, system_e_id = None , vobject_index = None):
            """ Function doc """
            self.system_e_id     = system_e_id    
            self.vobject_index = vobject_index
            ##print(system_e_id, vobject_index)
            
            system = self.treeview.main.p_session.psystem[self.system_e_id] 
            self.tree_header_sys_menu.set_label(system.label)
            
            if vobject_index == -1:
                
                self.tree_view_sys_menu.popup(None, None, None, None, 0, 0)

            if vobject_index != None and vobject_index != -1:
                
                vismol_object = self.treeview.main.vm_session.vm_objects_dic[vobject_index]
                self.tree_header_vobj_menu.set_label(vismol_object.name)
                
                self.tree_view_vobj_menu.popup(None, None, None, None, 0, 0)
                    
        def _save_backup_file (self):
            """ Function doc """
            self.main.p_session.save_easyhybrid_session( filename = self.main.session_filename, tmp = True)
