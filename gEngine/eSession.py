#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  simple.py
#
#  Copyright 2022 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
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

import logging
import gi, sys
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
from vismol.core.vismol_session import VismolSession
from gui.windows.windows_and_dialogs import EasyHybridDialogPrune
from vismol.core.vismol_selections import VismolViewingSelection as VMSele
import numpy as np
logger = logging.getLogger(__name__)


class GLMenu:
    """ Class doc """
    def insert_glmenu (self, bg_menu  = None, 
                            sele_menu = None, 
                             obj_menu = None, 
                            pick_menu = None):
        """ Function doc """
        



        def _viewing_selection_mode_atom (_):
            """ Function doc """
            self.viewing_selection_mode(sel_type = 'atom')
        
        def _viewing_selection_mode_residue (_):
            """ Function doc """
            self.viewing_selection_mode(sel_type = 'residue')
        
        def _viewing_selection_mode_chain (_):
            """ Function doc """
            self.viewing_selection_mode(sel_type = 'chain')

        def _selection_type_picking(self, widget):
            if self.selection_box_frame:
                self.selection_box_frame.change_toggle_button_selecting_mode_status(True)
            else:
                self.vm_session.picking_selection_mode = True
            self.vm_glcore.queue_draw()
        
        def _selection_type_viewing(self, widget):
            if self.selection_box_frame:
                self.selection_box_frame.change_toggle_button_selecting_mode_status(False)
            else:
                self.vm_session.picking_selection_mode = False
            self.vm_glcore.queue_draw()

        if sele_menu is None:
            ''' Standard Sele Menu '''
            
            def menu_show_dynamic_bonds (_):
                """ Function doc """
               ##print('dynamic_test')
                self.show_or_hide( rep_type = 'dynamic_bonds', show = True)
            
            def menu_hide_dynamic_bonds (_):
                """ Function doc """
               ##print('dynamic_test')
                self.show_or_hide( rep_type = 'dynamic_bonds', show = False)
            
            def menu_show_ribbons (_):
                """ Function doc """
               ##print('dynamic_test')
                self.show_or_hide( rep_type = 'ribbons', show = True)
            
            def menu_hide_ribbons (_):
                """ Function doc """
               ##print('dynamic_test')
                self.show_or_hide( rep_type = 'ribbons', show = False)
            
            def select_test (_):
                """ Function doc """
                self.select(indexes = 'all')
            
            def menu_show_lines (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'lines', show = True)
            
            def menu_show_dotted_lines (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'dotted_lines', show = True)

            def menu_hide_lines (_):
                """ Function doc """
                ##print('hide')
                self.show_or_hide( rep_type = 'lines', show = False)
            
            def menu_hide_dotted_lines (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'dotted_lines', show = False)
            
            def menu_show_sticks (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'sticks', show = True)
            
            def menu_show_nonbonded (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'nonbonded', show = True)
            
            def menu_hide_nonbonded (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'nonbonded', show = False)

            def menu_hide_sticks (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'sticks', show = False)

            def menu_show_spheres (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'spheres', show = True)

            def menu_hide_spheres (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'spheres', show = False)
            
            def menu_show_dots (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'dots', show = True)

            def menu_hide_dots (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'dots', show = False)
            
            
            
            def menu_set_color_grey (_):
                """ Function doc """
                self.set_color(color = [0.3     , 0.3     , 0.5 ] )
            
            def menu_set_color_green (_):
                """ Function doc """
                self.set_color(color = [0.0     , 1.0     , 0.0 ] )
            
            def menu_set_color_yellow (_):
                """ Function doc """
                self.set_color(color = [1.0     , 1.0     , 0.0 ] )
            
            def menu_set_color_light_blue (_):
                """ Function doc """
                self.set_color(color = [0.5     , 0.5     , 1.0 ] )
            
            def menu_set_color_light_red (_):
                """ Function doc """
                self.set_color(color = [1.0     , 0.5     , 0.5 ] )
            
            def menu_set_color_purple (_):
                """ Function doc """
                self.set_color(color = [1.0     , 0.0     , 1.0 ] )
            
            def menu_set_color_orange (_):
                """ Function doc """
                self.set_color(color = [1.0     , 0.5     , 0.0 ] )

            def menu_set_color_magenta (_):
                """ Function doc """

            def menu_color_change (_):
                """ Function doc """
                selection               = self.selections[self.current_selection]
                self.colorchooserdialog = Gtk.ColorChooserDialog()
                
                if self.colorchooserdialog.run() == Gtk.ResponseType.OK:
                    color = self.colorchooserdialog.get_rgba()
                    #print(color.red,color.green, color.blue )
                    new_color = [color.red, color.green, color.blue]

                self.colorchooserdialog.destroy()
                self.set_color(new_color)

            def set_as_qc_atoms (_):
                """ Function doc """
                #selection = self.selections[self.current_selection]
                active_id = self.main.p_session.active_id
                qc_list, residue_dict, vismol_object = self.build_index_list_from_atom_selection(return_vobject = True )
                 
                if qc_list:
                    
                    self.main.p_session.psystem[active_id].e_qc_residue_table = residue_dict
                    self.main.p_session.psystem[active_id].e_qc_table = qc_list
                    self.main.run_dialog_set_QC_atoms(vismol_object = vismol_object)

            def set_as_free_atoms (_):
                """ Function doc """
                selection = self.selections[self.current_selection]
                
                freelist = []                
                for atom in selection.selected_atoms:
                    ##print(atom.index, atom.name, atom.color) 
                    
                    '''checks if the selected atoms belong to the active project'''
                    true_or_false = self.check_selected_atom(atom, dialog = True)
                    if true_or_false:
                        freelist.append(atom.index -1)
                    else:
                        return False
                    
                '''here we are returning the original color of the selected atoms'''
                for key, vobject in self.vm_objects_dic.items():
                    if vobject.e_id == self.main.p_session.active_id:
                        #'''
                        for index in freelist:
                           ##print(index,vobject. atoms[index])
                            atom = vobject.atoms[index]
                            atom.color = atom._init_color()
                
                ##print(atom.index, atom.name, atom.color) 
                #----------------------------------------------
                pdmsys_active =   self.main.p_session.active_id
                
                a = set(self.main.p_session.psystem[pdmsys_active].e_fixed_table)
                b = set(freelist)
                
                c = a - b

                fixedlist =  set(self.main.p_session.psystem[pdmsys_active].e_fixed_table) -set(freelist)
                #guarantee that the atom index appears only once in the list
                fixedlist = list(c) 
                
 
                refresh = self.main.p_session.define_free_or_fixed_atoms_from_iterable (fixedlist)
                if refresh:
                    for key, vobject in self.vm_objects_dic.items():
                        if vobject.e_id == self.main.p_session.active_id:
                            self.main.p_session._apply_fixed_representation_to_vobject(system_id = None, vismol_object = vobject)
                            self.main.p_session._apply_QC_representation_to_vobject(system_id = None, vismol_object = vobject)                    
                            self.vm_glcore.queue_draw()
            
            
            def prune_atoms (_):
                """ Function doc """
                
                
                #builder = Gtk.Builder()
                #builder.add_from_file(os.path.join(self.home,'gui/windows/windows_and_dialogs.glade'))
                


                atomlist, resi_table = self.build_index_list_from_atom_selection()
                
                if atomlist:
                    atomlist = list(set(atomlist))
                    
                    num_of_atoms = len(atomlist)
                    name = self.main.p_session.psystem[self.main.p_session.active_id].label
                    tag  = self.main.p_session.psystem[self.main.p_session.active_id].e_tag
                    
                    dialog =  EasyHybridDialogPrune(self.main.home ,num_of_atoms, name, tag)
                    
       

                    if dialog.prune:
                        #print ("Prune")
                        name         = dialog.name        
                        tag          = dialog.tag  
                        
                        #for row in self.treestore:
                        #    #row[2] = row.path == selected_path
                        #    row[3] =  False
                        #
                        self.main.p_session.prune_system (selection = atomlist, name = name, summary = True, tag = tag)
            
            def set_as_fixed_atoms (_):
                """ Function doc """
                
                fixedlist, sel_resi_table = self.build_index_list_from_atom_selection()
                
                if fixedlist:
                    pdmsys_active = self.main.p_session.active_id
                    fixedlist = list(fixedlist) + list(self.main.p_session.psystem[pdmsys_active].e_fixed_table)
                    #guarantee that the atom index appears only once in the list
                    fixedlist = list(set(fixedlist)) 
                    #print ('fixedlist',fixedlist)
                    #sending to pDynamo
                    
                    refresh = self.main.p_session.define_free_or_fixed_atoms_from_iterable (fixedlist)
                    if refresh:
                        
                        for v_id, vobject in self.vm_objects_dic.items():
                            if vobject.e_id == self.main.p_session.active_id:
                                self.main.p_session._apply_fixed_representation_to_vobject(system_id = None, vismol_object = vobject)
                                self.main.p_session._apply_QC_representation_to_vobject(system_id = None, vismol_object = vobject)
                                self.vm_glcore.queue_draw()
                    #self.main.p_session.vismol_selection_qc = selection.copy()
            
            
            def add_selection_to_sel_list (_):
                """ Function doc """
                ##print('self.selections[self.current_selection].invert_selection()')
                sel_list, sel_resi_table = self.build_index_list_from_atom_selection()
                if sel_list:
                    
                    self.main.p_session.add_a_new_item_to_selection_list (system_id = self.main.p_session.active_id, 
                                                                                       indexes = sel_list, 
                                                                                        )

                
                #self.selections[self.current_selection].invert_selection()
            
            def invert_selection (_):
                """ Function doc """
                ##print('self.selections[self.current_selection].invert_selection()')
                self.selections[self.current_selection].invert_selection()
            
            
            sele_menu = { 
                    #'header' : ['MenuItem', None],
                    
                    'separator0':['separator', None],
                    
                    'Send to Selection List':['MenuItem', add_selection_to_sel_list],
                    
                    'separator1':['separator', None],
                    
                    
                    'show'   : [
                                'submenu' ,{
                                            
                                            'lines'         : ['MenuItem', menu_show_lines],
                                            #'dotted_lines'  : ['MenuItem', menu_show_dotted_lines],
                                            'sticks'        : ['MenuItem', menu_show_sticks],
                                            'spheres'       : ['MenuItem', menu_show_spheres],
                                            #'dots'          : ['MenuItem', menu_show_dots],
                                            'dynamic bonds' : ['MenuItem', menu_show_dynamic_bonds],
                                            'ribbons'       : ['MenuItem', menu_show_ribbons],
                                            'separator2'    : ['separator', None],
                                            'nonbonded'     : ['MenuItem', menu_show_nonbonded],
                    
                                           }
                               ],
                    
                    
                    'hide'   : [
                                'submenu',  {
                                            'lines'         : ['MenuItem', menu_hide_lines],
                                            #'dotted_lines'  : ['MenuItem', menu_hide_dotted_lines],
                                            'sticks'        : ['MenuItem', menu_hide_sticks],
                                            'spheres'       : ['MenuItem', menu_hide_spheres],
                                            #'dots'          : ['MenuItem', menu_hide_dots],
                                            'dynamic bonds' : ['MenuItem', menu_hide_dynamic_bonds],
                                            'ribbons' : ['MenuItem', menu_hide_ribbons],

                                            'separator2'    : ['separator', None],
                                            'nonbonded'     : ['MenuItem', menu_hide_nonbonded],
                                            }
                                ],
                    
                    'color'   : [
                                'submenu',  {
                                            'grey'          : ['MenuItem', menu_set_color_grey],
                                            'yellow'        : ['MenuItem', menu_set_color_yellow],
                                            'green'         : ['MenuItem', menu_set_color_green],
                                            'light_blue'    : ['MenuItem', menu_set_color_light_blue],
                                            'light_red'     : ['MenuItem', menu_set_color_light_red],
                                            'purple'        : ['MenuItem', menu_set_color_purple],
                                            'orange'        : ['MenuItem', menu_set_color_orange],
                                            'custon'        : ['MenuItem', menu_color_change],
                                            #'dotted_lines'  : ['MenuItem', menu_hide_dotted_lines],
                                            }
                                ],
                    
                    'Invert Selection':['MenuItem', invert_selection],
                    
                    

                    
                    'separator2':['separator', None],

                    
                    
                    'Selection type'   : [
                                'submenu' ,{
                                            
                                            'viewing'   :  ['MenuItem', _selection_type_viewing],
                                            'picking'   :  ['MenuItem', _selection_type_picking],
                                            #'separator2':['separator', None],
                                            #'nonbonded' : ['MenuItem', None],
                    
                                           }
                                        ],
                    
                    'Selection Mode'   : [
                                'submenu' ,{
                                            
                                            'Atoms'     :  ['MenuItem', _viewing_selection_mode_atom],
                                            'Residue'   :  ['MenuItem', _viewing_selection_mode_residue],
                                            'Chain'     :  ['MenuItem', _viewing_selection_mode_chain],
                                            #'separator2':['separator', None],
                                            #'nonbonded' : ['MenuItem', None],
                    
                                           }
                               ],
                    
                    'separator3':['separator', None],
                    
                    'Set as QC atoms'      :  ['MenuItem', set_as_qc_atoms],
                    
                    'separator4':['separator', None],

                    'Set as fixed atoms'   :  ['MenuItem', set_as_fixed_atoms],
                    'Set as free atoms'   :  ['MenuItem', set_as_free_atoms],
                    
                    'separator5':['separator', None],
                    'prune to selection'  :  ['MenuItem', prune_atoms],

                    'separator6':['separator', None],

                    
                    #'Label Mode':  ['submenu' , {
                    #                        'Atom'         : [
                    #                                           'submenu', {
                    #                                                       'lines'    : ['MenuItem', None],
                    #                                                       'sticks'   : ['MenuItem', None],
                    #                                                       'spheres'  : ['MenuItem', None],
                    #                                                       'nonbonded': ['MenuItem', None],
                    #                                                       }
                    #                                          ],
                    #                        
                    #                        'Atom index'   : ['MenuItem', None],
                    #                        'residue name' : ['MenuItem', None],
                    #                        'residue_index': ['MenuItem', None],
                    #                       },
                    #          ]
                    }
      
        
        if bg_menu is None:
            ''' Standard Bg Menu'''
            
            def open_structure_data (_):
                """ Function doc """
                ##print('ebaaaa')
                #self.filechooser   = FileChooser()
                #filename = self.filechooser.open()
                #self.load (filename, widget = None, autocenter = True)
                self.main.selection_list_window.OpenWindow()
            
            def import_system_menu (_):
                """ Function doc """
                self.main.NewSystemWindow.OpenWindow()
            
            def active_selection (_):
                """ Function doc """
                self.selections[self.current_selection].active = True
                self.vm_glcore.queue_draw()
                
            bg_menu = { 

                    'separator3'       : ['separator', None],
                    'Active Selection' : ['MenuItem', active_selection],

                    
                    'separator0'          : ['separator', None],
                    'Show Selection list' : ['MenuItem', open_structure_data],
                    
                    
                    #'funcao teste' : ['MenuItem', self.teste],                  
                    #'funcao teste2': ['MenuItem', self.teste2], 

                    'separator1':['separator', None],
                    'Selection type'   : [
                                'submenu' ,{
                                            
                                            'viewing'   :  ['MenuItem', _selection_type_viewing],
                                            'picking'   :  ['MenuItem', _selection_type_picking],
                                            #'separator2':['separator', None],
                                            #'nonbonded' : ['MenuItem', None],
                    
                                           }
                                        ],
                    
                    'Selection Mode'   : [
                                'submenu' ,{
                                            
                                            'atoms'     :  ['MenuItem', _viewing_selection_mode_atom],
                                            'residue'   :  ['MenuItem', _viewing_selection_mode_residue],
                                            'chain'     :  ['MenuItem', _viewing_selection_mode_chain],
                                            #'separator2':['separator', None],
                                            #'nonbonded' : ['MenuItem', None],
                    
                                           }
                               ],
                    
                    
                    #'hide'   : [
                    #            'submenu',  {
                    #                        'lines'    : ['MenuItem', menu_hide_lines],
                    #                        'sticks'   : ['MenuItem', menu_hide_sticks],
                    #                        'spheres'  : ['MenuItem', menu_hide_spheres],
                    #                        'nonbonded': ['MenuItem', None],
                    #                        }
                    #            ],
                    
                    
                    'separator2':['separator', None],
                    'Import System' : ['MenuItem', import_system_menu],
                    'separator3':['separator', None],

                    
                    
                    #'label':  ['submenu' , {
                    #                        'Atom'         : [
                    #                                           'submenu', {
                    #                                                       'lines'    : ['MenuItem', None],
                    #                                                       'sticks'   : ['MenuItem', None],
                    #                                                       'spheres'  : ['MenuItem', None],
                    #                                                       'nonbonded': ['MenuItem', None],
                    #                                                       }
                    #                                          ],
                    #                        
                    #                        'Atom index'   : ['MenuItem', None],
                    #                        'residue name' : ['MenuItem', None],
                    #                        'residue_index': ['MenuItem', None],
                    #                       },
                    #           ]
                    }

        
        if obj_menu is None:
            ''' Standard Obj Menu'''
            obj_menu = { 
                    'OBJ menu' : ['MenuItem', None],
                    
                    
                    'separator1':['separator', None],
                    
                    
                    'show'   : [
                                'submenu' ,{
                                            
                                            'lines'    : ['MenuItem', menu_show_lines],
                                            'sticks'   : ['MenuItem', menu_show_sticks],
                                            'spheres'  : ['MenuItem', menu_show_spheres],
                                            'separator2':['separator', None],
                                            'nonbonded': ['MenuItem', None],
                    
                                           }
                               ],
                    
                    
                    'hide'   : [
                                'submenu',  {
                                            'lines'    : ['MenuItem', menu_hide_lines],
                                            'sticks'   : ['MenuItem', menu_hide_sticks],
                                            'spheres'  : ['MenuItem', menu_hide_spheres],
                                            'nonbonded': ['MenuItem', None],
                                            }
                                ],
                    
                    
                    'separator2':['separator', None],

                    
                    
                    'label':  ['submenu' , {
                                            'Atom'         : [
                                                               'submenu', {
                                                                           'lines'    : ['MenuItem', None],
                                                                           'sticks'   : ['MenuItem', None],
                                                                           'spheres'  : ['MenuItem', None],
                                                                           'nonbonded': ['MenuItem', None],
                                                                           }
                                                              ],
                                            
                                            'atomic index' : ['MenuItem', None],
                                            'residue name' : ['MenuItem', None],
                                            'residue_index': ['MenuItem', None],
                                           },
                               ]
                    }


        if pick_menu is None:
            ''' Standard Sele Menu '''
            pick_menu = { 
                    'header' : ['MenuItem', None],
                    
                    
                    
                    'separator1':['separator', None],
                    
                    
                    'show'   : [
                                'submenu' ,{
                                            
                                            'lines'         : ['MenuItem', menu_show_lines],
                                            'sticks'        : ['MenuItem', menu_show_sticks],
                                            'spheres'       : ['MenuItem', menu_show_spheres],
                                            'dynamic bonds' : ['MenuItem', menu_show_dynamic_bonds],
                                            'separator2'    : ['separator', None],
                                            'nonbonded'     : ['MenuItem', None],
                    
                                           }
                               ],
                    
                    
                    'hide'   : [
                                'submenu',  {
                                            'lines'    : ['MenuItem', menu_hide_lines],
                                            'sticks'   : ['MenuItem', menu_hide_sticks],
                                            'spheres'  : ['MenuItem', menu_hide_spheres],
                                            'dynamic bonds' : ['MenuItem', menu_hide_dynamic_bonds],

                                            'nonbonded': ['MenuItem', None],
                                            }
                                ],
                    
                    
                    'separator2':['separator', None],

                    }


        self.vm_widget.insert_glmenu( bg_menu   = bg_menu, 
                                     sele_menu = sele_menu, 
                                     obj_menu  = obj_menu, 
                                     pick_menu = pick_menu)




class EasyHybridSession(VismolSession, GLMenu):
    """ Class doc """
    
    def __init__ (self, a = None):
        """ Class initialiser """
        super().__init__(toolkit="Gtk_3.0")
        ##print('\n\n\n',a, '\n\n\n')
    
    def _selection_function (self, selected, _type = None, disable = True):
        #"""     P I C K I N G     S E L E C T I O N S     """
        ##print('_selection_function')
        if self._picking_selection_mode:
            self.picking_selections.selection_function_picking(selected)
        
        #"""     V I E W I N G     S E L E C T I O N S     """
        else:
            self.selections[self.current_selection].selection_function_viewing(selected, _type, disable)

    
    def viewing_selection_mode(self, sel_type="atom"):
        """ Function doc """
        ##print(self.selection_box_frane)
        
        if self.selection_box_frane:
            #print('Selection mode:', sel_type)
            self.selection_box_frane.change_sel_type_in_combobox(sel_type)
        self.selections[self.current_selection].selection_mode = sel_type

    
    def _add_vismol_object(self, vismol_object, show_molecule=True, autocenter=True):
        """ Function doc """
        if vismol_object.index in self.vm_objects_dic.keys():
            logger.warning("The VismolObject with id {} already exists. \
                The data will be overwritten.".format(vismol_object.index))
        
        self.vm_objects_dic[self.vm_object_counter] = vismol_object
        vismol_object.index = self.vm_object_counter
        self.vm_object_counter +=1
        #self.atom_id_counter += len(vismol_object.atoms)
        
        for atom in vismol_object.atoms.values():
            self.atom_dic_id[atom.unique_id] = atom
        
        if show_molecule:
            vismol_object.create_representation(rep_type="lines")
            vismol_object.create_representation(rep_type="nonbonded")
            #vismol_object.create_representation(reprep_type="sticks")
            if autocenter:
                self.vm_glcore.center_on_coordinates(vismol_object, vismol_object.mass_center)
            else:
                self.vm_glcore.queue_draw()
        return vismol_object


    def build_index_list_from_atom_selection (self, return_vobject = False ):
        """  
        returns the index_list and residue_dict
        
        Residue_dict is a dictionary where the access key is the index 
        of the residue. The key gives access to a list of atoms 
        """
        selection    = self.selections[self.current_selection]
        
        index_list   = []  
        residue_dict = {} 
        vismol_object = None
        
        for atom in selection.selected_atoms:
            ##print(atom.vobject.easyhybrid_system_id , pdmsys_active)
            true_or_false = self.check_selected_atom (atom)
            if true_or_false:
                
                index_list.append(atom.index -1)
                
                vismol_object = atom.vm_object
                
                ##print(atom.residue, atom.index -1)
                
                if atom.residue.index in residue_dict:
                    residue_dict[atom.residue.index].append(atom.index -1)
                else:
                    residue_dict[atom.residue.index] = [atom.index -1]

                
            else:
                return False
        
        if return_vobject:
            return index_list, residue_dict, vismol_object
        else:
            return index_list, residue_dict


    def check_selected_atom(self, atom, dialog = True):
        '''checks if selected atoms belong to the dynamo system in memory'''
        if atom.vm_object.e_id != self.main.p_session.active_id:
            ##print(atom.index-1, atom.name, atom.resn)
            
            name = self.main.p_session.psystem[self.main.p_session.active_id].label
            
            dialog = Gtk.MessageDialog(
                        transient_for = self.main.window,
                        flags=0,
                        message_type=Gtk.MessageType.INFO,
                        buttons=Gtk.ButtonsType.OK,
                        text="Invalid Atom Selection",
                        )
            dialog.format_secondary_text(
"""Your atom selection does not belong to the active pDynamo system:
 
{} ({}) 

You can choose the active pDynamo system by changing the radio 
button position in the main treeview (active column).""".format(name,self.main.p_session.active_id)
            )
            dialog.run()
            #print("INFO dialog closed")
            dialog.destroy()
            return False
        else:
            return True
        

    
    #def set_color (self, color = [0.5 , 0.5 , 0.5]):
    #    """ Function doc """
    #    selection         = self.selections[self.current_selection]
    #
    #    vm_object_list = []
    #
    #    for atom in selection.selected_atoms:
    #        if atom.symbol == 'C':
    #            atom.color = color
    #            if atom.vm_object in vm_object_list:
    #                pass
    #            else:
    #                vm_object_list.append(atom.vm_object)
    #    
    #    self.vm_glcore.queue_draw()
    #    for vm_object in vm_object_list:
    #        vm_object._generate_color_vectors(True)
    #        
    #        for rep  in vm_object.representations.keys():
    #            if vm_object.representations[rep]:
    #                try:
    #                    vm_object.representations[rep]._load_color_vbo(vm_object.colors)
    #                #
    #                except:
    #                    pass
    #                    print("VisMol/vModel/Representations.py, line 123, in _set_colors_to_buffer GL.glBindBuffer(GL.GL_ARRAY_BUFFER, ctypes.ArgumentError: argument 2: <class 'TypeError'>: wrong type'")
    #                    
    #
    #    self.vm_glcore.queue_draw()
    #    return True


    def set_color_by_index (self, vobject = None, indexes = [ ], color = [0.9, 0.9, 0.9] ):
        """ Function doc """
        for atom_index in indexes:
            vobject.atoms[atom_index].color = color    

        vobject._generate_color_vectors ( 
                                          True

                                               )
        self.vm_glcore.queue_draw()
        for rep  in vobject.representations.keys():
            if vobject.representations[rep]:
                #try:
                try:
                    vobject.representations[rep]._load_color_vbo(vobject.colors)
                #
                except:
                    pass
                    print("VisMol/vModel/Representations.py, line 123, in _set_colors_to_buffer GL.glBindBuffer(GL.GL_ARRAY_BUFFER, ctypes.ArgumentError: argument 2: <class 'TypeError'>: wrong type'")
                    
        return True


    def create_new_selection (self):
        """ Function doc """
        return VMSele(self)




    def forward_frame(self):
        """ Function doc """
        frame = self.frame + 1
        #print('frame',frame)
        for i, vm_object in enumerate(self.vm_objects_dic.values()):
            if frame < vm_object.frames.shape[0]:
                self.frame += 1
                self.vm_glcore.updated_coords = True
                break
            else:
                pass
        else:
            self.vm_glcore.updated_coords = False
    
    def reverse_frame(self):
        """ Function doc """
        if self.frame - 1 >= 0:
            self.frame -= 1
            self.vm_glcore.updated_coords = True
        else:
            self.vm_glcore.updated_coords = False
    
    def set_frame(self, frame=0):
        """ Function doc """
        assert frame >= 0
        self.frame = np.uint32(frame)
        self.vm_glcore.updated_coords = True
        self.vm_widget.queue_draw()
