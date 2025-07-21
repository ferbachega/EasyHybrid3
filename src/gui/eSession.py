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
import vismol.utils.selectors as selectors

#from gui.widgets.custom_widgets import get_distance

from gui.windows.setup.windows_and_dialogs import EasyHybridDialogPrune
from gui.windows.setup.windows_and_dialogs import AddHarmonicRestraintDialog
      
from vismol.core.vismol_selections import VismolViewingSelection as VMSele
from vismol.core.vismol_selections import VismolPickingSelection as VMPick
import numpy as np
logger = logging.getLogger(__name__)

from util.geometric_analysis            import get_distance 
from util.geometric_analysis            import get_dihedral 
from util.geometric_analysis            import get_angle 

class CommandLine:
    """ Class doc """
    
    def __init__ (self, vm_session):
        """ Class initialiser """
        pass
        self.vm_session  = vm_session
        #self.main        = self.vm_session.main
        #self.p_session   = self.vm_session.main.p_session
        '''''' 
        

    def command_line_parser (self, cmd_text):
        """             
            command:
            
            func arg1 = value1, arg2 = value2+value3... , arg3 = ... 
        
        """
        #obtaing the function name
        func = cmd_text.split()
        func = func[0]
        
        #obtating the arguments 
        args_raw = cmd_text.replace(func, '')
        args_raw = args_raw.replace(' ', '')
        #args_raw = args_raw.strip()
        args_raw = args_raw.split(',')
        
        print (args_raw)
        args = []
        for argr in args_raw:
            argr = argr.split('=')
            for arg in argr:
                #print(arg.strip())
                args.append(arg.strip())
        
        
        #func =arg_raw data[0]
        #args = data[1:]     
        #print ('data', func, args)
        
        return func, args
        
    
    def run(self, cmd ):
        """ Function doc"""
        
        func, args = self.command_line_parser(cmd)
        #print(func, args)
        _func = getattr(self, func)
        #print(_func)
        try:
            log = _func(args)
        except AttributeError as ae:
            logger.error("Command '{}' not implemented".format(func))
            logger.error(ae)
            log = "Command '{}' not implemented".format(func)
        return log
    
            
    def list (self, args = None):
        """ Function doc """
        text = 'here is list'
        if args == [] or args == ['']:
            for e_id, system in self.vm_session.main.p_session.psystem.items():
                print('system: {} - {}'.format(e_id, system.label))               
                
                text += '\nsystem: {} - {}'.format(e_id, system.label)
                for index,  vobject in self.vm_session.vm_objects_dic.items():
                    if vobject.e_id == e_id:
                        print('          {} - {}'.format(index, vobject.name))
                        text += '\n         {} - {}'.format(index, vobject.name)
        else:
            pass
            
        
        return text
            
        
        '''
        #if system_id:
        #    
        #    print('system {} {}:'.format(system_id, self.vm_session.main.p_session.psystem[system_id].label))
        #    for index,  vobject in self.vm_session.vm_objects_dic.items():
        #        if vobject.e_id == system_id:
        #            print('     {} {}:'.format(index, vobject.name))
        #        
        else:
            for e_id, system in self.vm_session.main.p_session.psystem.items():
                
                print('system {} {}:'.format(e_id, system.label))
                
                for index,  vobject in self.vm_session.vm_objects_dic.items():
                    if vobject.e_id == e_id:
                        print('     {} {}:'.format(index, vobject.name))
        '''


    def select (sys = None, obj = None, name = [], symbol = [], resi = [], resn = [], chain = []):
        """ Function doc 
        
        select 
        
        """

    
    def show (self, args):
        """ Function doc 
        
        example:
        
        show type = lines, sele = sel01 
        
        """
        print('print show', args)
        

    def clear (self, args = None):
        """ Function doc """
        self.vm_session.main.terminal_text_buffer.set_text("")
    
    

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
        
        def _viewing_selection_mode_molecule (_):
            """ Function doc """
            self.viewing_selection_mode(sel_type = 'molecule')
        


        #def _selection_type_picking(self, widget):
        def _selection_type_picking(_):
            
            if self.selection_box_frame:
                self.selection_box_frame.change_toggle_button_selecting_mode_status(True)
            else:
                self.picking_selection_mode = True
            self.vm_glcore.queue_draw()
        
        #def _selection_type_viewing(self, widget):
        def _selection_type_viewing(_):
            if self.selection_box_frame:
                self.selection_box_frame.change_toggle_button_selecting_mode_status(False)
            else:
                self.picking_selection_mode = False
            self.vm_glcore.queue_draw()

        if sele_menu is None:
            ''' Standard Sele Menu '''
            
            
            def menu_show_atom_name (_):
                #print('menu_show_atom_name') 
                selection = self.show_or_hide( rep_type = 'labels', show = True)
                
                for atom in selection.selected_atoms:
                    atom.label_text = atom.name
            
            
            def menu_show_atom_symbol (_):   
                selection = self.show_or_hide( rep_type = 'labels', show = True)
                
                for atom in selection.selected_atoms:
                    atom.label_text = atom.symbol
            
            
            def menu_show_atom_index (_):    
                selection = self.show_or_hide( rep_type = 'labels', show = True)
                for atom in selection.selected_atoms:
                    atom.label_text = str(atom.index)
                    
            def menu_show_atom_MM_charge (_):               
                selection = self.show_or_hide( rep_type = 'labels', show = True)
                for atom in selection.selected_atoms:
                    VObj = atom.vm_object
                    atom.label_text  = '%4.3f'%(float(self.main.p_session.psystem[VObj.e_id].mmState.charges[atom.index-1]))            

            def menu_show_residue_name  (_): 
                selection = self.show_or_hide( rep_type = 'labels', show = True)
                for atom in selection.selected_atoms:
                    #VObj = atom.vm_object
                    atom.label_text = atom.residue.name              
            
            def menu_show_residue_index  (_):
                selection = self.show_or_hide( rep_type = 'labels', show = True)
                for atom in selection.selected_atoms:
                    #VObj = atom.vm_object
                    atom.label_text = str(atom.residue.index)             
            
            def menu_show_chain (_):
                selection = self.show_or_hide( rep_type = 'labels', show = True)
                for atom in selection.selected_atoms:
                    #VObj = atom.vm_object
                    atom.label_text = str(atom.chain.name)  
            
            def menu_hide_label (_):
                """ Function doc """
                selection = self.show_or_hide( rep_type = 'labels', show = False)
            
            
            
            def menu_show_dynamic_bonds (_):
                """ Function doc """
               ##print('dynamic_test')
                sele = self.show_or_hide( rep_type = 'dynamic', show = True)
                #print(sele.get_selection_info())
            def menu_hide_dynamic_bonds (_):
                """ Function doc """
               ##print('dynamic_test')
                self.show_or_hide( rep_type = 'dynamic', show = False)
            
            def menu_show_ribbons (_):
                """ Function doc """
               ##print('dynamic_test')
                self.show_or_hide( rep_type = 'ribbons', show = True)
                self.show_or_hide( rep_type = 'ribbon_sphere', show = True)
                print('ribbon_sphere')
            def menu_hide_ribbons (_):
                """ Function doc """
               ##print('dynamic_test')
                self.show_or_hide( rep_type = 'ribbons', show = False)
                self.show_or_hide( rep_type = 'ribbon_sphere', show = False)

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
                self.show_or_hide( rep_type = 'stick_spheres', show = True)
                #self.show_or_hide( rep_type = 'dash', show = True)
            
            def menu_show_nonbonded (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'nonbonded', show = True)
            
            def menu_hide_nonbonded (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'nonbonded', show = False)

            def menu_hide_sticks (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'sticks', show = False)
                self.show_or_hide( rep_type = 'stick_spheres', show = False)

            def menu_show_spheres (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'spheres', show = True)

            def menu_hide_spheres (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'spheres', show = False)
            
            def menu_show_vdw_spheres (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'vdw_spheres', show = True)

            def menu_hide_vdw_spheres (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'vdw_spheres', show = False)
            
            def menu_show_dots (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'dots', show = True)

            def menu_hide_dots (_):
                """ Function doc """
                self.show_or_hide( rep_type = 'dots', show = False)
            
            def menu_hide_hydrogens (_):
                """ Function doc """
                #if selection is None:
                selection = self.selections[self.current_selection]
                #for a
                
                vobjects = {}
                for atom in selection.selected_atoms:
                    # is atom a hydrogen
                    if atom.symbol == 'H':
                        if atom.vm_object.index in vobjects.keys():
                            vobjects[atom.vm_object.index].append(atom.index-1)
                        else:
                            vobjects[atom.vm_object.index] = []
                            vobjects[atom.vm_object.index].append(atom.index-1)
                
                #print(vobjects)
                
                self.selections[self.current_selection]= VMSele(self)
                
                for vobj_id, indexes in vobjects.items():
                    '''----------------------------- Applying the selection ------------------------------------'''
                    vobject = self.vm_objects_dic[vobj_id]
                    #print ('1 - ',self.selections[self.current_selection].selected_atoms)
                    self.selections[self.current_selection].selecting_by_indexes (vismol_object = vobject, indexes = indexes, clear = True)
                    self.selections[self.current_selection].active = True
                    '''-----------------------------------------------------------------------------------------'''
                    #print ('2 - ',self.selections[self.current_selection].selected_atoms)

                    self.show_or_hide( rep_type = 'lines', show = False)
                    self.show_or_hide( rep_type = 'sticks', show = False)
                    self.show_or_hide( rep_type = 'spheres', show = False)

            
            #def menu_set_color_grey (_):
            #    """ Function doc """
            #    self.set_color(color = [0.3     , 0.3     , 0.5 ] )
            #
            #def menu_set_color_green (_):
            #    """ Function doc """
            #    self.set_color(color = [0.0     , 1.0     , 0.0 ] )
            #
            #def menu_set_color_yellow (_):
            #    """ Function doc """
            #    self.set_color(color = [1.0     , 1.0     , 0.0 ] )
            #
            #def menu_set_color_light_blue (_):
            #    """ Function doc """
            #    self.set_color(color = [0.5     , 0.5     , 1.0 ] )
            #
            #def menu_set_color_light_red (_):
            #    """ Function doc """
            #    self.set_color(color = [1.0     , 0.5     , 0.5 ] )
            #
            #def menu_set_color_purple (_):
            #    """ Function doc """
            #    self.set_color(color = [1.0     , 0.0     , 1.0 ] )
            #
            #def menu_set_color_orange (_):
            #    """ Function doc """
            #    self.set_color(color = [1.0     , 0.5     , 0.0 ] )
            #
            #def menu_set_color_magenta (_):
            #    """ Function doc """

            def menu_color_change (_):
                """ Function doc """
                selection               = self.selections[self.current_selection]
                self.colorchooserdialog = Gtk.ColorChooserDialog()
                
                if self.colorchooserdialog.run() == Gtk.ResponseType.OK:
                    color = self.colorchooserdialog.get_rgba()
                    #print(color.red,color.green, color.blue )
                    new_color = [color.red, color.green, color.blue]
                
                
                
                
                self.colorchooserdialog.destroy()
                self.set_color(color = new_color)

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
        


                atomlist, resi_table = self.build_index_list_from_atom_selection()
                
                if atomlist:
                    atomlist = list(set(atomlist))
                    
                    num_of_atoms = len(atomlist)
                    name = self.main.p_session.psystem[self.main.p_session.active_id].label
                    tag  = self.main.p_session.psystem[self.main.p_session.active_id].e_tag
                    
                    dialog =  EasyHybridDialogPrune(self.main  ,num_of_atoms, name, tag)
                    
       

                    if dialog.prune:
                        #print ("Prune")
                        name         = dialog.name        
                        tag          = dialog.tag  
                        color        = dialog.color 
                        vobject_id   = dialog.vobject_id
                        
                        vobject = self.main.vm_session.vm_objects_dic[vobject_id]
                        self.main.p_session.get_coordinates_from_vobject_to_pDynamo_system(vobject)
                        
                        
                        #for row in self.treestore:
                        #    #row[2] = row.path == selected_path
                        #    row[3] =  False
                        #print('color', color)

                        self.main.p_session.prune_system (selection = atomlist, name = name, summary = True, tag = tag, color = color)
            
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
            
            def call_selection_modify_window (_):
                """ Function doc """
                self.main.pDynamo_selection_window.OpenWindow()
            
            
            sele_menu = { 
                    #'header' : ['MenuItem', None],
                    
                    '- Selection Menu -':['MenuItem', None],
                    
                    'separator0':['separator', None],
                    
                    'Send to Selection List':['MenuItem', add_selection_to_sel_list],
                    
                    'separator1':['separator', None],
                    
                    
                    
                    'Extend Selection' :['MenuItem', call_selection_modify_window],
                    
                    #'Modify' :['submenu',  
                    #
                    #                        {
                    #                        'Around'          : ['submenu',
                    #                                                        {
                    #                                                        'Atom within 4 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Atom within 5 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Atom within 6 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Atom within 8 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Atom within 12 A'         : ['MenuItem', menu_show_lines],
                    #                                                        'Atom within 15 A'         : ['MenuItem', menu_show_lines],
                    #                                                        'separator9'    : ['separator', None],
                    #                                                        'Residues within 4 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Residues within 5 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Residues within 6 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Residues within 8 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Residues within 12 A'         : ['MenuItem', menu_show_lines],
                    #                                                        'Residues within 15 A'         : ['MenuItem', menu_show_lines],
                    #                                                        'separator2'    : ['separator', None],
                    #                                                        'Molecules within 4 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Molecules within 5 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Molecules within 6 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Molecules within 8 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Molecules within 12 A'         : ['MenuItem', menu_show_lines],
                    #                                                        'Molecules within 15 A'         : ['MenuItem', menu_show_lines],
                    #                                                        }
                    #
                    #                                            ],
                    #                        
                    #                        'Expand'          : ['submenu',
                    #                                                        {
                    #                                                        'Atom within 4 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Atom within 5 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Atom within 6 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Atom within 8 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Atom within 12 A'         : ['MenuItem', menu_show_lines],
                    #                                                        'Atom within 15 A'         : ['MenuItem', menu_show_lines],
                    #                                                        'separator9'    : ['separator', None],
                    #                                                        'Residues within 4 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Residues within 5 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Residues within 6 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Residues within 8 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Residues within 12 A'         : ['MenuItem', menu_show_lines],
                    #                                                        'Residues within 15 A'         : ['MenuItem', menu_show_lines],
                    #                                                        'separator2'    : ['separator', None],
                    #                                                        'Molecules within 4 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Molecules within 5 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Molecules within 6 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Molecules within 8 A'         : ['MenuItem', menu_show_lines] ,
                    #                                                        'Molecules within 12 A'         : ['MenuItem', menu_show_lines],
                    #                                                        'Molecules within 15 A'         : ['MenuItem', menu_show_lines],
                    #                                                        }
                    #
                    #                                            ],
                    #                        
                    #                        
                    #                        'Invert'        : ['MenuItem', menu_set_color_yellow],
                    #
                    #                       }
                    #          ],
                    
                    'Show'   : [
                                'submenu' ,{
                                            
                                            'labels'        : [
                                                               'submenu' ,{
                                                                           'Name'          : ['MenuItem', menu_show_atom_name     ],      
                                                                           'Symbol'        : ['MenuItem', menu_show_atom_symbol   ],      
                                                                           'Index'         : ['MenuItem', menu_show_atom_index    ],      
                                                                           'Charge(MM)'    : ['MenuItem', menu_show_atom_MM_charge],      
                                                                           'Residue Name'  : ['MenuItem', menu_show_residue_name  ],      
                                                                           'Residue Index' : ['MenuItem', menu_show_residue_index ],      
                                                                           'Chain'         : ['MenuItem', menu_show_chain         ],      
                                                                           
                                                                           }
                                                              ],
                                            
                                            'lines'         : ['MenuItem', menu_show_lines],
                                            #'dotted_lines'  : ['MenuItem', menu_show_dotted_lines],
                                            'sticks'        : ['MenuItem', menu_show_sticks],
                                            'dynamic bonds' : ['MenuItem', menu_show_dynamic_bonds],
                                            'separator1'    : ['separator', None],
                                            'spheres'       : ['MenuItem', menu_show_spheres],
                                            'vdw_spheres'   : ['MenuItem', menu_show_vdw_spheres],
                                            #'dots'          : ['MenuItem', menu_show_dots],
                                            
                                            'separator2'    : ['separator', None],

                                            'ribbons'       : ['MenuItem', menu_show_ribbons],
                                            'separator3'    : ['separator', None],
                                            'nonbonded'     : ['MenuItem', menu_show_nonbonded],
                    
                                           }
                               ],
                    
                    
                    'Hide'   : [
                                'submenu',  {
                                            'labels'        : ['MenuItem', menu_hide_label],
                                                              #[
                                                              # 'submenu' ,{
                                                              #             'Name'          : ['MenuItem', menu_hide_label],      
                                                              #             'Symbol'        : ['MenuItem', menu_hide_label],      
                                                              #             'Index'         : ['MenuItem', menu_hide_label],      
                                                              #             'Charge(MM)'    : ['MenuItem', menu_hide_label],      
                                                              #             'Residue Name'  : ['MenuItem', menu_hide_label],      
                                                              #             'Residue Index' : ['MenuItem', menu_hide_label],      
                                                              #             'Chain'         : ['MenuItem', menu_hide_label],      
                                                              #             
                                                              #             }
                                                              #],
                                            
                                            
                                            'lines'         : ['MenuItem', menu_hide_lines],
                                            #'dotted_lines'  : ['MenuItem', menu_hide_dotted_lines],
                                            'sticks'        : ['MenuItem', menu_hide_sticks],
                                            'dynamic bonds' : ['MenuItem', menu_hide_dynamic_bonds],
                                            'separator1'    : ['separator', None],

                                            'spheres'       : ['MenuItem', menu_hide_spheres],
                                            'vdw_spheres'   : ['MenuItem', menu_hide_vdw_spheres],
                                            #'dots'          : ['MenuItem', menu_hide_dots],
                                            'separator2'    : ['separator', None],

                                            'ribbons' : ['MenuItem', menu_hide_ribbons],

                                            'separator3'    : ['separator', None],
                                            'nonbonded'     : ['MenuItem', menu_hide_nonbonded],
                                            'separator4'    : ['separator', None],
                                            'hydrogens'     : ['MenuItem', menu_hide_hydrogens],

                                            }
                                ],
                    
                    'Change Color'   : ['MenuItem', menu_color_change],
                                #['submenu',  {
                                #            'grey'          : ['MenuItem', menu_set_color_grey],
                                #            'yellow'        : ['MenuItem', menu_set_color_yellow],
                                #            'green'         : ['MenuItem', menu_set_color_green],
                                #            'light_blue'    : ['MenuItem', menu_set_color_light_blue],
                                #            'light_red'     : ['MenuItem', menu_set_color_light_red],
                                #            'purple'        : ['MenuItem', menu_set_color_purple],
                                #            'orange'        : ['MenuItem', menu_set_color_orange],
                                #            'separator1'    : ['separator', None],
                                #
                                #            'custon'        : ['MenuItem', menu_color_change],
                                #            'separator2'    : ['separator', None],
                                #
                                #            #'dotted_lines'  : ['MenuItem', menu_hide_dotted_lines],
                                #            }
                                #],
                    

                    

                    
                    'separator2':['separator', None],

                    
                    
                    'Selection Type'   : [
                                'submenu' ,{
                                            
                                            'Viewing'   :  ['MenuItem', _selection_type_viewing],
                                            'Picking'   :  ['MenuItem', _selection_type_picking],
                                            #'separator2':['separator', None],
                                            #'nonbonded' : ['MenuItem', None],
                    
                                           }
                                        ],
                    
                    'Selecting By:'   : [
                                'submenu' ,{
                                            
                                            'Atoms'     :  ['MenuItem', _viewing_selection_mode_atom],
                                            'Residue'   :  ['MenuItem', _viewing_selection_mode_residue],
                                            'Chain'     :  ['MenuItem', _viewing_selection_mode_chain],
                                            'Molecule'  :  ['MenuItem', _viewing_selection_mode_molecule],
                                            #'separator2':['separator', None],
                                            #'nonbonded' : ['MenuItem', None],
                    
                                           }
                               ],
                    
                    'separator3':['separator', None],
                    
                    'Set as QC Atoms'      :  ['MenuItem', set_as_qc_atoms],
                    
                    'separator4':['separator', None],

                    'Set as Fixed Atoms'   :  ['MenuItem', set_as_fixed_atoms],
                    'Set as Free Atoms'   :  ['MenuItem', set_as_free_atoms],
                    
                    'separator5':['separator', None],
                    'Prune to Selection'  :  ['MenuItem', prune_atoms],

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
                    'Selection Type'   : [
                                'submenu' ,{
                                            
                                            'viewing'   :  ['MenuItem', _selection_type_viewing],
                                            'picking'   :  ['MenuItem', _selection_type_picking],
                                            #'separator2':['separator', None],
                                            #'nonbonded' : ['MenuItem', None],
                    
                                           }
                                        ],
                    
                    'Selecting By:'   : [
                                'submenu' ,{
                                            
                                            'atoms'     :  ['MenuItem', _viewing_selection_mode_atom],
                                            'residue'   :  ['MenuItem', _viewing_selection_mode_residue],
                                            'chain'     :  ['MenuItem', _viewing_selection_mode_chain],
                                            'molecule'  :  ['MenuItem',_viewing_selection_mode_molecule],

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
            
            def center_on_atom (_):
                """ Function doc """
                #print('center')
                #print(self.vm_glcore.info_atom)
                self.vm_glcore.center_on_atom(self.vm_glcore.info_atom)
            obj_menu = { 
                    'OBJ menu' : ['MenuItem', None],
                    
                    
                    'separator1':['separator', None],
                    
                    'Center' : ['MenuItem', center_on_atom],
                    #'Zoom' : ['MenuItem', active_selection],
                    #'Orient' : ['MenuItem', active_selection],
                    
                    'separator2':['separator', None],
                    
                    
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
                    
                    
                    'separator3':['separator', None],

                    
                    
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
            
            def add_harmonic_restraint(_):
                """ Function doc """
                print('add_harmonic_restraint')
                atom1 = self.picking_selections.picking_selections_list[0]
                atom2 = self.picking_selections.picking_selections_list[1]
                atom3 = self.picking_selections.picking_selections_list[2]
                atom4 = self.picking_selections.picking_selections_list[3]
                if atom1:
                    self.vobject = atom1.vm_object
                else:
                    return None
                
                #--------------------------------------------------------------
                if atom1:
                    print(str(atom1.index-1),str(atom1.name) )
                else: print('use picking selection to chose the central atom')            
                #--------------------------------------------------------------
                if atom2:
                    print(str(atom2.index-1),str(atom2.name) )
                else: print('use picking selection to chose the central atom')            
                #--------------------------------------------------------------
                if atom2 and atom1:
                    parameters = {}
                    dist  = get_distance(atom1.vm_object, atom1.index-1, atom2.index-1)
                    parameters['atom1']  = atom1.index-1
                    parameters['atom2']  = atom2.index-1
                    parameters['system'] = self.main.p_session.psystem[self.vobject.e_id]
                    add_harmonic_restraint_dialog =  AddHarmonicRestraintDialog(self.main, atom1,  atom2, dist )
                    
                    #print(add_harmonic_restraint_dialog.ok)
                    
                    if add_harmonic_restraint_dialog.ok:
                        #print(add_harmonic_restraint_dialog.dist)
                        parameters['distance']       = float(add_harmonic_restraint_dialog.dist )
                        parameters['force_constant'] = float(add_harmonic_restraint_dialog.force)
                    
                    
                    
                        #if atom2 and atom1:
                        self.main.p_session.add_new_harmonic_restraint(parameters)
                        self.main.selection_list_window.update_window (selections = False, restraints = True)
            
            
            pick_menu = { 
                    'header' : ['MenuItem', None],
                    
                    
                    
                    'separator1'              :['separator', None],
                    'Add Harmonic Restraint'  :['MenuItem', add_harmonic_restraint],
                    
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
    
    def __init__ (self, vm_config = None):
        
        """ Class initialiser """
        super().__init__(toolkit="Gtk_3.0", vm_config = vm_config)
        
        ##print('\n\n\n',a, '\n\n\n')
        self.selection_box_frame = None
        self.cmd = CommandLine(self)
        #self.vm_widget.connect_after("button-press-event", self.meu_evento_personalizado)
        
        
        
    
    
    def restart (self):
        """ Function doc """
        self.picking_selection_mode = False # True/False  - interchange between viewing  and picking mode
        self.selections = {"sel_00": VMSele(self)}
        self.current_selection = "sel_00"
        self.picking_selections = VMPick(self)
        
        for key, vobject in self.vm_objects_dic.items():
            vobject.active = False
        
        self.vm_objects_dic = {}
        self.vm_glcore.queue_draw()  
        
          
    def _selection_function (self, selected, _type = None, disable = True):
        #"""     P I C K I N G     S E L E C T I O N S     """
        ##print('_selection_function')
        if self._picking_selection_mode:
            self.picking_selections.selection_function_picking(selected)
        
        #"""     V I E W I N G     S E L E C T I O N S     """
        else:
            self.selections[self.current_selection].selection_function_viewing(selected, _type, disable)


    def _selection_function_set(self, selected, _type=None, disable=True):
        """ Function doc """
        
        #print('selected', selected)
        if self.picking_selection_mode: # True for picking mode
            if selected:
                #assert len(selected) == 1 # bachega 06 / 18 /2025
                selected = list(selected)[0]
                self.picking_selections.selection_function_picking(selected)
            else:
                self.picking_selections.selection_function_picking(None)
        else: # False for viewing mode
            self.selections[self.current_selection].selection_function_viewing_set(selected, _type, disable)
        #this will refresh the sequence canvas
        self.main.bottom_notebook.seqview.text_drawing_area.queue_draw()
    
    def viewing_selection_mode(self, sel_type="atom"):
        """ Function doc """
        ##print(self.selection_box_frame)
        
        if self.selection_box_frame:
            #print('Selection mode:', sel_type)
            self.selection_box_frame.change_sel_type_in_combobox(sel_type)
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
            #vismol_object.create_representation(rep_type="surface")
            
            #vismol_object.create_representation(rep_type="dash")
            #vismol_object.create_representation(rep_type="sticks")

            #'''
            #from vismol.libgl.representations import LabelRepresentation
            ##print (vismol_object.cell_parameters)
            ##vismol_object.representations["labels"] =  LabelRepresentation(vismol_object, self.vm_glcore, name  = 'lines', active=True, indexes = vismol_object.cell_bonds)
            #vismol_object.representations["labels"] =  LabelRepresentation(vismol_object = vismol_object  ,  
            #                                                               vismol_glcore = self.vm_glcore , 
            #                                                               indexes       = [0,1,2] , 
            #                                                               labels        = None     , 
            #                                                               color         = [1, 1, 0, 1])
            #'''
            
            
            #from vismol.libgl.representations import DashedLinesRepresentation
            #vismol_object.representations["restraints"] = DashedLinesRepresentation(vismol_object, self.vm_glcore,
            #                                                                  active=True, indexes= [5, 17])
            #
            #vismol_object.representations["restraints"].define_new_indexes_to_vbo([5, 17 ])
            #
            
            self.main.p_session.update_restaint_representation(vismol_object.e_id)
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
                    #print("VisMol/vModel/Representations.py, line 123, in _set_colors_to_buffer GL.glBindBuffer(GL.GL_ARRAY_BUFFER, ctypes.ArgumentError: argument 2: <class 'TypeError'>: wrong type'")
                    
        return True


    def set_color (self, symbol = 'C', color = [0.9, 0.9, 0.0] ):
        """ Function doc """
        
        selection = self.selections[self.current_selection]
                
        atomlist = []
        
        vobjects = {}
                        
        for atom in selection.selected_atoms:
            if atom.symbol == symbol:
                #atomlist.append(atom.index-1)
                
                if atom.vm_object.index in vobjects.keys():
                    vobjects[atom.vm_object.index].append(atom.index-1)
                else:
                    vobjects[atom.vm_object.index] = [atom.index-1]
                    
                #vobject = atom.vm_object
        
        for vob_id, atomlist in vobjects.items():
            vobject = self.vm_objects_dic[vob_id]
                
            self.set_color_by_index (vobject = vobject, indexes = atomlist, color = color ) 
            self.main.p_session._apply_fixed_representation_to_vobject(system_id = None, vismol_object = vobject)
            self.main.p_session._apply_QC_representation_to_vobject(system_id = None, vismol_object = vobject)
    
    
    def create_new_selection (self):
        """ Function doc """
        return VMSele(self)


    def forward_frame(self):
        """ Function doc """
        #if self.main_session:
        #    print('here')
        #    self.main_session.trajectory_player_window.self.vm_traj_obj.forward()
        #else:
        frame = self.frame + 1
        self.set_frame(frame=frame)

        '''
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
        #'''
  

    def reverse_frame(self):
        """ Function doc """
        
        #if self.main_session.trajectory_player_window.Visible:
        #    self.main_session.trajectory_player_window.self.vm_traj_obj.reverse()
        #else:
       
        frame = self.frame
        if self.frame - 1 >= 0:
            frame -= 1
        else:
            frame  = 0
        self.set_frame(frame=frame)
        
        '''
        if self.frame - 1 >= 0:
            self.frame -= 1
            self.vm_glcore.updated_coords = True
        else:
            self.vm_glcore.updated_coords = False
        #'''

    
    def set_frame(self, frame=0):
        """ Function doc """
        assert frame >= 0
        self.frame = np.uint32(frame)
        self.vm_glcore.updated_coords = True
        
        if self.picking_selection_mode:
            self.picking_selections.update_pki_pkj_rep_coordinates()
            self.picking_selections.print_pk_distances()    
        
        #self.main_session.label_frame.set_text('Frame Number: ' + str(frame)) 
        self.vm_widget.queue_draw()
        #print('aqui')
        #self.main.surface_analysis_window.set_frame()
        '''
        for vobj_id, vobject in self.vm_objects_dic.items():
            for index, atom in vobject.atoms.items():
                for bond in atom.bonds:
                    print(index-1, atom.name, bond.get_indexes() )
        '''
  
        
    def define_vismol_object_molecules (self, vobject):
        """ Function doc """
        
        self.verified = set()
        
        self.molecules = {0 : []}
        

    def _define_inner_box (self, selection, grid_size):
        """  
        Establishes the inner box (smallest box in Cartesian space) that encompasses the selected atoms
        
        inputs:
                selection =  vismol selection object 
                grid_size =  float
                
        return :
                vobjects,             = list of involved object
                grid_min,             = list (xyz) 
                grid_max,             = list (xyz)
                selected_indexes_dict = dict, the keys ate the vobject.index -> list os selected indexes 
        
                                        max(xyz) = grid_max 
                  |-------|-------|-------|
                  |       |       |      x|
                  |    x  |       | x x x |
                  |     xx|       |x      |
                  |-------|-------|-------|
                  |       |xx  xxx|       |
                  |       |   x   |       |
                  |       | x   x |       |
                  |-------|-------|-------|
                  |      x|x      |x      |
                  |   x x |       | x     |
                  |x      |       |       |
                  |-------|-------|-------|
                   ^
        min(xyz) = grid_min

        
        """
        print ('Defining the selection inner box')
        i = []
        j = []
        k = []
        selected_indexes = set()
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        vobject  = None
        vobjects = []
        
        selected_indexes_dict = {}
        selected_indexes_dict = {}
        
        for atom in selection.selected_atoms:
            '''checks if the selected atoms belong to the active project'''
            true_or_false = self.check_selected_atom(atom, dialog = True)
            
            if true_or_false:
                xyz = atom.coords()
                
                vobject = atom.vm_object
                
                if vobject in vobjects:
                    pass
                else:
                    vobjects.append(vobject)
                
                #grid_pos = (int(xyz[0]/grid_size), int(xyz[1]/grid_size), int(xyz[2]/grid_size))
                #i.append(grid_pos[0])
                #j.append(grid_pos[1])
                #k.append(grid_pos[2])
                i.append(xyz[0])
                j.append(xyz[1])
                k.append(xyz[2])
                
                #selected_indexes.add(atom.index-1)
                
                if vobject.index in selected_indexes_dict.keys():
                    selected_indexes_dict[vobject.index].append(atom.index-1)
                else:
                    selected_indexes_dict[vobject.index]=[]
                    selected_indexes_dict[vobject.index].append(atom.index-1)


                #if grid_pos in selected_grip.keys():
                #    selected_grip[grid_pos].append(atom)
                #else:
                #    selected_grip[grid_pos] = []
                #    selected_grip[grid_pos].append(atom)
            else:
                print('invalid selection list')
                return False
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        #print ('\nkey (grid pos) / num of atoms',
        #'\n max i: ', max(i) , 'min i: ', min(i),
        #'\n max j: ', max(j) , 'min j: ', min(j),
        #'\n max k: ', max(k) , 'min k: ', min(k))
        
        grid_min = [min(i), min(j), min(k)]
        grid_max = [max(i), max(j), max(k)]
        return vobjects, grid_min, grid_max, selected_indexes_dict
    
    
    def _define_selectable_indexes (self, vobject, grid_min, grid_max, radius, grid_size = 10):
        """  
        Establishes the list of selectable atoms (big box)
        Sets up a box that contains selectable atoms (which will be checked later). 
        The box size is established by the chosen search radius.
        
        inputs:
                vobject   =  vismol  object 
                grid_min  =  list
                grid_max  =  list
                radius    =  float
                grid_size = float
                
        return :
                selectable_indexes = set()
        
        
                                        max(xyz) = grid_max 
           |-------|-------|-------|-------|-------|-------|
           |                                               |
           |      ||-------------------------------||      |
           |      ||       |       |       |       ||      |
           |      ||-------|-box of selectable-----||      |
           |      ||       |     atoms     |       ||      |
           |      ||       |       |       |       ||      |
           |      ||       |-------------- |       ||      |
           |      ||------||-------|-------||------||      |
           |      ||      ||       |       ||      ||      |
           |      ||      || inner box     ||      ||      |
           |      ||      ||  (selected    ||      ||      |
           |      ||------||----atoms)-----||------||      |
           |      ||      ||       |       ||      ||      |                  ^
           |      ||      ||       |       ||      ||      |
           |      ||      ||       |       ||      ||      |
           |      ||------||-------|-------||------||      |
           |      ||       |---------------|       ||      |
           |      ||       |       |       |       ||      |
           |      ||       |       |       |       ||      |
           |      ||-------|-------|-------|-------||      |
           |      ||       |       |       |       ||      |                  ^
           |      ||-------------------------------||      |
           |                                               |
           |-------|-------|-------|-------|-------|-------|


        min(xyz) = grid_min

        
        
        """
        
        print ('Defining the selectable indexes')
        
        selectable_indexes = set()
        grid = {}
        radius = radius 

        for index, atom in vobject.atoms.items():
            xyz = atom.coords()
            
            grid_pos = (int(xyz[0]/grid_size), 
                        int(xyz[1]/grid_size), 
                        int(xyz[2]/grid_size))
            
            
            if grid_pos[0] >= grid_min[0]-radius and grid_pos[0] <= grid_max[0]+radius:
                if grid_pos[1] >= grid_min[1]-radius and grid_pos[1] <= grid_max[1]+radius:                
                    if grid_pos[2] >= grid_min[2]-radius and grid_pos[2] <= grid_max[2]+radius:
                        #----------------------------
                        selectable_indexes.add(index)
                        #----------------------------
            
            
            if grid_pos in grid.keys():
                grid[grid_pos].append(atom)
            else:
                grid[grid_pos] = []
                grid[grid_pos].append(atom)
        
        return selectable_indexes

        
    def advanced_selection (self, selection = None, _type = 'Around' ,selecting_by = 'Residue',   radius = 10, grid_size = 10):
        """ Function doc """
        
        grid_size = radius / grid_size
        
        if selection is None:
            selection = self.selections[self.current_selection]
        else:
            pass        
        
        try:
            vobjects, grid_min, grid_max, selected_indexes_dict = self._define_inner_box(selection, grid_size)       
        except:
            return False, 'Selection error!\n\nThe initial selection must contain at least one selected atom.'
        
        self._selection_function_set(None)

        for vobject in vobjects:
            selected_indexes = set(selected_indexes_dict[vobject.index])
            selectable_indexes = self._define_selectable_indexes(vobject, grid_min, grid_max, radius, grid_size )
            selectable_indexes = selectable_indexes - selected_indexes
            #print('vobject.c_alpha_bonds', vobject.c_alpha_bonds )
            #print('vobject.c_alpha_atoms', vobject.c_alpha_atoms )
    
            
            try:
                coordinates = vobject.frames[self.frame]
            except:
                coordinates = vobject.frames[-1]
            
            new_selected_indexes, selectable_indexes = selectors.selection_spherical_expansion( 
                                                                                            selected_indexes, 
                                                                                            selectable_indexes, 
                                                                                            coordinates, 
                                                                                            radius )
            
            #-------------------------------------------------------------------------------------
            if _type == 'Around':
                new_selected_indexes = new_selected_indexes - selected_indexes
            
            
            elif _type == 'Expand':
                for index in selected_indexes:
                    new_selected_indexes.add(index)
            
            elif _type == 'ByComplement':
                keys = vobject.atoms.keys()
                keys = set(keys)
                
                new_selected_indexes = keys - new_selected_indexes
                #for index in selected_indexes:
                #    new_selected_indexes.add(index)
            
            
            else:
                pass
            #print ('_type', _type, 'selecting_by', selecting_by)
            self.selections[self.current_selection].selecting_by_indexes( vobject, 
                                                                new_selected_indexes, 
                                                                clear=True)
            #-------------------------------------------------------------------------------------
            
            
            
            
            #-------------------------------------------------------------------------------------
            if selecting_by == 'Residue':
                new_selected_indexes = self._complement_by_residue( new_selected_indexes,vobject )
            
            elif selecting_by == 'Molecule':
                new_selected_indexes = self._complement_by_molecule( new_selected_indexes, vobject )
            
            else:
                pass
            #-------------------------------------------------------------------------------------
            
            
    
                
            self.selections[self.current_selection].active = True
            self.vm_glcore.queue_draw()

        return True, 'Selection done!'
 
    
    def _complement_by_residue (self, new_selected_indexes, vobject ):
        """ Function doc """
        # - - - - - - - - - - selecting all the residues! - - - - - - - - - - - - - - - - - - - -  
        for atom in self.selections[self.current_selection].selected_atoms:
            keys = atom.residue.atoms.keys()
            for key in keys:
                new_selected_indexes.add(key)
        
        self.selections[self.current_selection].selecting_by_indexes( vobject, 
                                                           new_selected_indexes, 
                                                           clear=True)
        return new_selected_indexes
 
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
    def _complement_by_molecule (self, new_selected_indexes, vobject ):
        """ Function doc """
        # - - - - - - - - - - selecting all the residues! - - - - - - - - - - - - - - - - - - - -  
        for atom in self.selections[self.current_selection].selected_atoms:
            keys = atom.molecule.atoms.keys()
            for key in keys:
                new_selected_indexes.add(key)
        
        self.selections[self.current_selection].selecting_by_indexes( vobject, 
                                                           new_selected_indexes, 
                                                           clear=True)
        return new_selected_indexes
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         

    def selection_around(self, selection = None, radius = 23, grid_size = 10): 
        """ Function doc """
        selection = {}
        for atom in self.selections[self.current_selection].selected_atoms:
            selection[atom.index-1] = atom
        
        bonds = self.vm_objects_dic[0].find_bonded_and_nonbonded_atoms( selection=selection, frame=self.frame)
                                        

    def show_cell (self, vismol_object):
        """ Function doc """
        rep_labels = vismol_object.representations.keys()
        
        if  "cell_lines" in rep_labels:
            vismol_object.representations["cell_lines"].active = True
        
        else:
            if vismol_object.cell_parameters:
                from vismol.libgl.representations import CellLineRepresentation
                print (vismol_object.cell_parameters)
                vismol_object.representations["cell_lines"] =  CellLineRepresentation(vismol_object, self.vm_glcore,name  = 'lines', active=True, indexes = vismol_object.cell_bonds)
        self.vm_glcore.queue_draw()
   
    
    def hide_cell (self, vismol_object):
        if  "cell_lines" in  vismol_object.representations.keys():
            vismol_object.representations["cell_lines"].active = False
        self.vm_glcore.queue_draw()
