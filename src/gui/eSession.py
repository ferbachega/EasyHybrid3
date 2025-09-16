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
import gi, sys, os
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
from vismol.core.vismol_session import VismolSession
import vismol.utils.selectors as selectors

#from gui.widgets.custom_widgets import get_distance

from gui.windows.setup.windows_and_dialogs import EasyHybridDialogPrune
from gui.windows.setup.windows_and_dialogs import AddHarmonicRestraintDialog
from gui.windows.setup.windows_and_dialogs import SimpleDialog

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
        
        #print (args_raw)
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
    """     
    Class responsible for handling graphical layer (GL) menus.

    This includes selection mode changes (atom, residue, chain, molecule),
    picking/viewing modes, and standard label display options. 
    """
    def insert_glmenu (self, bg_menu  = None, 
                            sele_menu = None, 
                             obj_menu = None, 
                            pick_menu = None):
        """
        Insert and configure GL menus for visualization and selection.

        Parameters
        ----------
        bg_menu : Gtk.Menu, optional
            Background menu (context menu for background actions).
        sele_menu : Gtk.Menu, optional
            Selection menu. If None, a standard default selection menu is created.
        obj_menu : Gtk.Menu, optional
            Object menu for object-specific actions.
        pick_menu : Gtk.Menu, optional
            Picking menu for atom/molecule picking.
        """
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
            '''Standard atom labeling menu.'''
            
            
            def menu_show_atom_name (_):
                """Show atom names as labels for the selected atoms.""" 
                selection = self.show_or_hide( rep_type = 'labels', show = True)
                
                for atom in selection.selected_atoms:
                    atom.label_text = atom.name
            
            
            def menu_show_atom_symbol (_):   
                """Show atom symbols as labels for the selected atoms."""
                selection = self.show_or_hide( rep_type = 'labels', show = True)
                
                for atom in selection.selected_atoms:
                    atom.label_text = atom.symbol
            
            
            def menu_show_atom_index (_):    
                """Show atom indices as labels for the selected atoms."""
                selection = self.show_or_hide( rep_type = 'labels', show = True)
                for atom in selection.selected_atoms:
                    atom.label_text = str(atom.index)
                    
            def menu_show_atom_MM_charge (_):               
                """Show the MM charge of each selected atom as its label."""
                selection = self.show_or_hide( rep_type = 'labels', show = True)
                for atom in selection.selected_atoms:
                    VObj = atom.vm_object
                    atom.label_text  = '%4.3f'%(float(self.main.p_session.psystem[VObj.e_id].mmState.charges[atom.index-1]))            

            def menu_show_residue_name  (_): 
                """Show the residue name for each selected atom."""
                selection = self.show_or_hide( rep_type = 'labels', show = True)
                for atom in selection.selected_atoms:
                    #VObj = atom.vm_object
                    atom.label_text = atom.residue.name              
            
            def menu_show_residue_index  (_):
                """Show the residue index for each selected atom."""
                selection = self.show_or_hide( rep_type = 'labels', show = True)
                for atom in selection.selected_atoms:
                    #VObj = atom.vm_object
                    atom.label_text = str(atom.residue.index)             
            
            def menu_show_chain (_):
                """Show the chain name for each selected atom."""
                selection = self.show_or_hide( rep_type = 'labels', show = True)
                for atom in selection.selected_atoms:
                    #VObj = atom.vm_object
                    atom.label_text = str(atom.chain.name)  
            
            def menu_hide_label (_):
                """Hide all labels for the selected atoms."""
                selection = self.show_or_hide( rep_type = 'labels', show = False)
            
            
            
            def menu_show_dynamic_bonds (_):
                """Show dynamic bond representations."""
                #print('dynamic_test')
                sele = self.show_or_hide( rep_type = 'dynamic', show = True)
                
            def menu_hide_dynamic_bonds (_):
                """Hide dynamic bond representations."""
                #print('dynamic_test')
                self.show_or_hide( rep_type = 'dynamic', show = False)
            
            def menu_show_ribbons (_):
                """Show ribbon and ribbon-sphere representations."""
                self.show_or_hide( rep_type = 'ribbons', show = True)
                self.show_or_hide( rep_type = 'ribbon_sphere', show = True)
                #print('ribbon_sphere')
            
            def menu_hide_ribbons (_):
                """Hide ribbon and ribbon-sphere representations."""
                self.show_or_hide( rep_type = 'ribbons', show = False)
                self.show_or_hide( rep_type = 'ribbon_sphere', show = False)

            def select_test (_):
                """
                Select all atoms.
                Not truly used
                """
                self.select(indexes = 'all')
            
            def menu_show_lines (_):
                """Show line representations."""
                self.show_or_hide( rep_type = 'lines', show = True)

            def menu_hide_lines (_):
                """Hide line representations."""
                self.show_or_hide( rep_type = 'lines', show = False)
            
            def menu_show_dotted_lines (_):
                """Show dotted line representations."""
                self.show_or_hide( rep_type = 'dotted_lines', show = True)

            def menu_hide_dotted_lines (_):
                """Hide dotted line representations."""
                self.show_or_hide( rep_type = 'dotted_lines', show = False)
            
            def menu_show_sticks (_):
                """Show stick and stick-sphere representations."""
                self.show_or_hide( rep_type = 'sticks', show = True)
                self.show_or_hide( rep_type = 'stick_spheres', show = True)
            
            def menu_show_nonbonded (_):
                """Show nonbonded atom representations."""
                self.show_or_hide( rep_type = 'nonbonded', show = True)
            
            def menu_hide_nonbonded (_):
                """Hide nonbonded atom representations."""
                self.show_or_hide( rep_type = 'nonbonded', show = False)

            def menu_hide_sticks (_):
                """Hide stick and stick-sphere representations."""
                self.show_or_hide( rep_type = 'sticks', show = False)
                self.show_or_hide( rep_type = 'stick_spheres', show = False)

            def menu_show_spheres (_):
                """Show sphere representations."""
                self.show_or_hide( rep_type = 'spheres', show = True)

            def menu_hide_spheres (_):
                """Hide sphere representations."""
                self.show_or_hide( rep_type = 'spheres', show = False)
            
            def menu_show_vdw_spheres (_):
                """Show van der Waals (vdW) sphere representations."""
                self.show_or_hide( rep_type = 'vdw_spheres', show = True)

            def menu_hide_vdw_spheres (_):
                """Hide van der Waals (vdW) sphere representations."""
                self.show_or_hide( rep_type = 'vdw_spheres', show = False)
            
            def menu_show_dots (_):
                """Show dot surface representations."""
                self.show_or_hide( rep_type = 'dots', show = True)

            def menu_hide_dots(_):
                """Hide dot surface representations."""
                self.show_or_hide(rep_type='dots', show=False)
            
            def menu_hide_hydrogens (_):
                """
                Hide all hydrogen atoms from the current selection in 
                different representations (lines, sticks, spheres)."""
                # Get the current selection
                selection = self.selections[self.current_selection]
                
                # Dictionary to group hydrogen atom indices by their corresponding vismol object (vm_object)
                vobjects = {}
                for atom in selection.selected_atoms:
                    # Check if the atom is a hydrogen
                    if atom.symbol == 'H':
                        vm_index = atom.vm_object.index

                        # Initialize the entry if not present
                        if vm_index not in vobjects:
                            vobjects[vm_index] = []

                        # Store atom index (subtract 1 to match internal indexing convention)
                        vobjects[vm_index].append(atom.index - 1)
                

                # Reset the current selection before applying hydrogen hiding
                self.selections[self.current_selection]= VMSele(self)
                
                for vobj_id, indexes in vobjects.items():
                    '''----------------------------- Applying the selection ------------------------------------'''
                    vobject = self.vm_objects_dic[vobj_id]
                    
                    # Apply selection by indexes
                    self.selections[self.current_selection].selecting_by_indexes (vismol_object = vobject, indexes = indexes, clear = True)
                    self.selections[self.current_selection].active = True
                    '''-----------------------------------------------------------------------------------------'''
                    # Hide hydrogens in all visual representation types
                    self.show_or_hide( rep_type = 'lines', show = False)
                    self.show_or_hide( rep_type = 'sticks', show = False)
                    self.show_or_hide( rep_type = 'stick_spheres', show = False)
                    
                    self.show_or_hide( rep_type = 'spheres', show = False)
                    self.show_or_hide( rep_type = 'vdw_spheres', show = False)

            def menu_color_change (_):
                """  
                Open a color chooser dialog and apply a new color to the selected atoms.

                The user can decide whether to make this color customization permanent
                for each system involved in the selection.
                """
                # Get the current selection
                selection               = self.selections[self.current_selection]
                
                # Create the GTK color chooser dialog
                self.colorchooserdialog = Gtk.ColorChooserDialog()
                
                
                # If the user confirms a color choice
                if self.colorchooserdialog.run() == Gtk.ResponseType.OK:
                    # Retrieve the RGBA color values
                    color = self.colorchooserdialog.get_rgba()
                    new_color = [color.red, color.green, color.blue]
                    
                    # Apply the temporary color to the visualization
                    self.set_color(color = new_color)
                    

                    '''
                    his dictionary will temporarily store the list 
                    of indices and the new color, with the access key 
                    being the e_id of each system (remember: a viewing 
                    selection can select more than one system at the 
                    same time).
                    '''
                    e_system = {}
                    
                    '''We will iterate through all atoms in the selection 
                    (it may include more than one vobject and more than 
                    one system).'''
                    # Iterate through all selected atoms across all vobjects/systems
                    for atom in selection.selected_atoms:
                        
                        vobject = atom.vm_object
                        e_id = getattr(vobject, 'e_id', None)

                        if e_id not in e_system:
                            # Initialize entry for this system
                            e_system[e_id] = {
                                'indexes': [],
                                'color': new_color
                            }

                        # Store atom index (subtract 1 to match internal indexing convention)
                        e_system[e_id]['indexes'].append(atom.index - 1)
                    
                    
                    # Ask the user if the customization should be permanent
                    msg = 'You have just selected a new list of atoms with customized colors. \nWould you like to make this color customization permanent?'
                    dialog = SimpleDialog(self.main_session)
                    yes_or_no = dialog.question (msg)
                    
                    
                    if yes_or_no:
                        '''Adding to each system the new customized colors 
                        and the lists of atoms (indices).'''
                        for e_id in e_system.keys():
                            system = self.main.p_session.psystem[e_id]
                            custom_colors = e_system[e_id]
                            system.e_custom_colors.append(custom_colors)
                    else:
                        pass
                
                # Close the color chooser dialog regardless of the result
                self.colorchooserdialog.destroy()
                

            def set_as_qc_atoms (_):
                """ 
                Mark the currently selected atoms as quantum chemistry (QC) atoms.
    
                This updates the QC atom and residue tables in the active system
                and opens a dialog for further QC atom configuration.
                """

                active_id = self.main.p_session.active_id
                
                # Build list of selected atom indices, residue mapping, and vobject
                qc_list, residue_dict, vismol_object = self.build_index_list_from_atom_selection(return_vobject = True )
                 
                if qc_list:
                    # Store QC information in the active system
                    system = self.main.p_session.psystem[active_id]
                    system.e_qc_residue_table = residue_dict
                    system.e_qc_table = qc_list
                    
                    # Open a dialog to confirm/adjust QC atom selection
                    self.main.run_dialog_set_QC_atoms(vismol_object = vismol_object)

            def set_as_free_atoms (_):
                """ 
                Mark the currently selected atoms as free (non-fixed).
    
                This updates the fixed/free atom tables in the active system and
                restores the original colors of the freed atoms.
                """
                selection = self.selections[self.current_selection]
                
                freelist = []                
                
                # Validate and collect indices of selected atoms
                for atom in selection.selected_atoms:
                    
                    # Ensure atom belongs to the active project/system 
                    true_or_false = self.check_selected_atom(atom, dialog = True)
                    if true_or_false:
                        freelist.append(atom.index -1)
                    else:
                        # Abort if invalid atom is selected
                        return False
                    
                '''here we are returning the original color of the selected atoms'''
                # Restore original colors for freed atoms
                for key, vobject in self.vm_objects_dic.items():
                    if vobject.e_id == self.main.p_session.active_id:
                        for index in freelist:
                            atom = vobject.atoms[index]
                            atom.color = atom._init_color()
                
                
                # Update fixed atom table in the active system
                pdmsys_active =   self.main.p_session.active_id
                system = self.main.p_session.psystem[pdmsys_active]
                
                '''
                a = set(self.main.p_session.psystem[pdmsys_active].e_fixed_table)
                b = set(freelist)
                
                c = a - b

                fixedlist =  set(self.main.p_session.psystem[pdmsys_active].e_fixed_table) -set(freelist)
                #guarantee that the atom index appears only once in the list
                fixedlist = list(c) 
                '''
                # Remove freed atoms from fixed table
                fixedlist = set(system.e_fixed_table) - set(freelist)
                fixedlist = list(fixedlist)  # Ensure unique indices
                
                
                # Apply updated free/fixed configuration
                refresh = self.main.p_session.define_free_or_fixed_atoms_from_iterable (fixedlist)
                
                if refresh:
                    # Apply visual updates to the corresponding vobject
                    for key, vobject in self.vm_objects_dic.items():
                        if vobject.e_id == self.main.p_session.active_id:
                            self.main.p_session._apply_fixed_representation_to_vobject(system_id = None, vismol_object = vobject)
                            self.main.p_session._apply_QC_representation_to_vobject(system_id = None, vismol_object = vobject)                    
                            
                            # Redraw the visualization
                            self.vm_glcore.queue_draw()
            
            
            def prune_atoms (_):
                """ 
                Prune (remove) the currently selected atoms from the active system.
    
                This opens a prune dialog where the user can confirm the operation,
                assign a name, tag, and color, and specify the target vobject.
                """
        

                # Build list of selected atom indices and residue mapping
                atomlist, resi_table = self.build_index_list_from_atom_selection()
                
                if atomlist:
                    # Ensure uniqueness of atom indices
                    atomlist = list(set(atomlist))
                    num_of_atoms = len(atomlist)
                    
                    # Get active system details
                    active_id = self.main.p_session.active_id
                    system = self.main.p_session.psystem[active_id]
                    name, tag = system.label, system.e_tag
                    
                    # Open prune confirmation dialog
                    dialog = EasyHybridDialogPrune(self.main, num_of_atoms, name, tag)

                    if dialog.prune:
                        # Collect user-defined prune parameters
                        name         = dialog.name        
                        tag          = dialog.tag  
                        color        = dialog.color 
                        vobject_id   = dialog.vobject_id
                        
                        # Sync vobject coordinates with pDynamo system
                        vobject = self.main.vm_session.vm_objects_dic[vobject_id]
                        self.main.p_session.get_coordinates_from_vobject_to_pDynamo_system(vobject)
                        
                        # Perform the prune operation
                        self.main.p_session.prune_system (selection = atomlist, name = name, summary = True, tag = tag, color = color)
            
            def set_as_fixed_atoms (_):
                """
                Mark the currently selected atoms as fixed in the active system.
    
                Updates the fixed atom table and applies the corresponding
                fixed/QC visual representations to the active vobject.
                """
                
                # Collect indices of selected atoms
                fixedlist, sel_resi_table = self.build_index_list_from_atom_selection()
                
                if fixedlist:

                    pdmsys_active = self.main.p_session.active_id
                    system = self.main.p_session.psystem[pdmsys_active]
                    
                    # Merge new selection with existing fixed atom table
                    fixedlist = list(set(fixedlist).union(system.e_fixed_table))
                    
                    # Update the fixed/free status in the pDynamo system
                    refresh = self.main.p_session.define_free_or_fixed_atoms_from_iterable (fixedlist)
                    if refresh:
                        # Apply visual updates to the active vobject
                        for v_id, vobject in self.vm_objects_dic.items():
                            if vobject.e_id == self.main.p_session.active_id:
                                self.main.p_session._apply_fixed_representation_to_vobject(system_id = None, vismol_object = vobject)
                                self.main.p_session._apply_QC_representation_to_vobject(system_id = None, vismol_object = vobject)
                                # Redraw the visualization
                                self.vm_glcore.queue_draw()
           
            
            def add_selection_to_sel_list (_):
                """Add the currently selected atoms to the selection list of the active system."""
                sel_list, sel_resi_table = self.build_index_list_from_atom_selection()

                if sel_list:
                    
                    self.main.p_session.add_a_new_item_to_selection_list (system_id = self.main.p_session.active_id, 
                                                                                       indexes = sel_list, 
                                                                                        )

                            
            def invert_selection (_):
                """Invert the current selection of atoms (selected ↔ unselected)."""
                self.selections[self.current_selection].invert_selection()
            
            def call_selection_modify_window (_):
                """ Function doc """
                self.main.pDynamo_selection_window.OpenWindow()
            
            
            sele_menu = {
                        # ---------------- Standard Selection Menu ----------------
                        '- Selection Menu -': ['MenuItem', None],

                        # Separator for visual grouping
                        'separator0': ['separator', None],

                        # Add current selection to the selection list
                        'Send to Selection List': ['MenuItem', add_selection_to_sel_list],

                        'separator1': ['separator', None],

                        # Extend or modify the current selection
                        'Extend Selection': ['MenuItem', call_selection_modify_window],

                        # ---------------- Show Submenu ----------------
                        'Show': [
                            'submenu', {
                                # Labels submenu for displaying atom/residue information
                                'labels': [
                                    'submenu', {
                                        'Name': ['MenuItem', menu_show_atom_name],
                                        'Symbol': ['MenuItem', menu_show_atom_symbol],
                                        'Index': ['MenuItem', menu_show_atom_index],
                                        'Charge(MM)': ['MenuItem', menu_show_atom_MM_charge],
                                        'Residue Name': ['MenuItem', menu_show_residue_name],
                                        'Residue Index': ['MenuItem', menu_show_residue_index],
                                        'Chain': ['MenuItem', menu_show_chain],
                                    }
                                ],

                                # Graphical representations
                                'lines': ['MenuItem', menu_show_lines],
                                'sticks': ['MenuItem', menu_show_sticks],
                                'dynamic bonds': ['MenuItem', menu_show_dynamic_bonds],

                                'separator1': ['separator', None],

                                'spheres': ['MenuItem', menu_show_spheres],
                                'vdw_spheres': ['MenuItem', menu_show_vdw_spheres],

                                'separator2': ['separator', None],

                                'ribbons': ['MenuItem', menu_show_ribbons],

                                'separator3': ['separator', None],

                                'nonbonded': ['MenuItem', menu_show_nonbonded],
                            }
                        ],

                        # ---------------- Hide Submenu ----------------
                        'Hide': [
                            'submenu', {
                                'labels': ['MenuItem', menu_hide_label],
                                'lines': ['MenuItem', menu_hide_lines],
                                'sticks': ['MenuItem', menu_hide_sticks],
                                'dynamic bonds': ['MenuItem', menu_hide_dynamic_bonds],

                                'separator1': ['separator', None],

                                'spheres': ['MenuItem', menu_hide_spheres],
                                'vdw_spheres': ['MenuItem', menu_hide_vdw_spheres],

                                'separator2': ['separator', None],

                                'ribbons': ['MenuItem', menu_hide_ribbons],

                                'separator3': ['separator', None],

                                'nonbonded': ['MenuItem', menu_hide_nonbonded],

                                'separator4': ['separator', None],

                                'hydrogens': ['MenuItem', menu_hide_hydrogens],
                            }
                        ],

                        # Change color for the selected atoms
                        'Change Color': ['MenuItem', menu_color_change],

                        'separator2': ['separator', None],

                        # ---------------- Selection Type Submenu ----------------
                        'Selection Type': [
                            'submenu', {
                                'Viewing': ['MenuItem', _selection_type_viewing],
                                'Picking': ['MenuItem', _selection_type_picking],
                            }
                        ],

                        # ---------------- Selecting Mode Submenu ----------------
                        'Selecting By:': [
                            'submenu', {
                                'Atoms': ['MenuItem', _viewing_selection_mode_atom],
                                'Residue': ['MenuItem', _viewing_selection_mode_residue],
                                'Chain': ['MenuItem', _viewing_selection_mode_chain],
                                'Molecule': ['MenuItem', _viewing_selection_mode_molecule],
                            }
                        ],

                        'separator3': ['separator', None],

                        # ---------------- QC and Fixed/Free Atom Operations ----------------
                        'Set as QC Atoms': ['MenuItem', set_as_qc_atoms],

                        'separator4': ['separator', None],

                        'Set as Fixed Atoms': ['MenuItem', set_as_fixed_atoms],
                        'Set as Free Atoms': ['MenuItem', set_as_free_atoms],

                        'separator5': ['separator', None],

                        'Prune to Selection': ['MenuItem', prune_atoms],

                        'separator6': ['separator', None],
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
                # ---------------- Separator and Active Selection ----------------
                'separator3': ['separator', None],
                'Active Selection': ['MenuItem', active_selection],

                'separator0': ['separator', None],

                # Show the list of all selections
                'Show Selection list': ['MenuItem', open_structure_data],

                # ---------------- Optional test functions (commented out) ----------------
                'separator1': ['separator', None],

                # ---------------- Selection Type Submenu ----------------
                'Selection Type': [
                    'submenu', {
                        'viewing': ['MenuItem', _selection_type_viewing],
                        'picking': ['MenuItem', _selection_type_picking],
                        # 'separator2': ['separator', None],
                        # 'nonbonded': ['MenuItem', None],
                    }
                ],

                # ---------------- Selecting Mode Submenu ----------------
                'Selecting By:': [
                    'submenu', {
                        'atoms': ['MenuItem', _viewing_selection_mode_atom],
                        'residue': ['MenuItem', _viewing_selection_mode_residue],
                        'chain': ['MenuItem', _viewing_selection_mode_chain],
                        'molecule': ['MenuItem', _viewing_selection_mode_molecule],
                    }
                ],

                # ---------------- Separator and Import ----------------
                'separator2': ['separator', None],
                'Import System': ['MenuItem', import_system_menu],
                'separator3': ['separator', None],

                # ---------------- Placeholder for label submenu (commented out) ----------------
                # 'label': ['submenu', {
                #     'Atom': [
                #         'submenu', {
                #             'lines': ['MenuItem', None],
                #             'sticks': ['MenuItem', None],
                #             'spheres': ['MenuItem', None],
                #             'nonbonded': ['MenuItem', None],
                #         }
                #     ],
                #     'Atom index': ['MenuItem', None],
                #     'residue name': ['MenuItem', None],
                #     'residue_index': ['MenuItem', None],
                # }]
            }
        
        if obj_menu is None:
            ''' Standard Obj Menu'''
            
            def center_on_atom (_):
                """ Function doc """
                #print('center')
                #print(self.vm_glcore.info_atom)
                self.vm_glcore.center_on_atom(self.vm_glcore.info_atom)
            
            obj_menu = {
                        # Header for the object menu
                        'OBJ menu': ['MenuItem', None],

                        # Separator for visual grouping
                        'separator1': ['separator', None],

                        # Center view on the selected atom
                        'Center': ['MenuItem', center_on_atom],

                        # Zoom and Orient could be added here
                        # 'Zoom': ['MenuItem', active_selection],
                        # 'Orient': ['MenuItem', active_selection],

                        'separator2': ['separator', None],

                        # ---------------- Show submenu ----------------
                        'show': [
                            'submenu', {
                                'lines': ['MenuItem', menu_show_lines],
                                'sticks': ['MenuItem', menu_show_sticks],
                                'spheres': ['MenuItem', menu_show_spheres],
                                'separator2': ['separator', None],
                                'nonbonded': ['MenuItem', None],  # Placeholder
                            }
                        ],

                        # ---------------- Hide submenu ----------------
                        'hide': [
                            'submenu', {
                                'lines': ['MenuItem', menu_hide_lines],
                                'sticks': ['MenuItem', menu_hide_sticks],
                                'spheres': ['MenuItem', menu_hide_spheres],
                                'nonbonded': ['MenuItem', None],  # Placeholder
                            }
                        ],

                        'separator3': ['separator', None],

                        # ---------------- Labels submenu ----------------
                        'label': [
                            'submenu', {
                                'Atom': [
                                    'submenu', {
                                        'lines': ['MenuItem', None],       # Placeholder
                                        'sticks': ['MenuItem', None],      # Placeholder
                                        'spheres': ['MenuItem', None],     # Placeholder
                                        'nonbonded': ['MenuItem', None],   # Placeholder
                                    }
                                ],
                                'atomic index': ['MenuItem', None],       # Placeholder
                                'residue name': ['MenuItem', None],       # Placeholder
                                'residue_index': ['MenuItem', None],      # Placeholder
                            }
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
    """
    Main session class for EasyHybrid.
    Inherits from VismolSession and GLMenu.
    Manages loading of project/session files and system initialization.
    """
    
    def __init__(self, vm_config=None):
        """
        Initialize an EasyHybridSession.

        Parameters
        ----------
        vm_config : dict, optional
            Configuration dictionary for the Vismol session.
        """
        super().__init__(toolkit="Gtk_3.0", vm_config=vm_config)
        
        # Frame for selection box (not initialized yet)
        self.selection_box_frame = None
        
        # Dictionary to store unique names for VObjects
        # Important to allow selections by name (e.g. from the command line)
        self.vobject_names = {} 

    def load(self, filename=None):
        """
        Load a file into the EasyHybrid session.
        The behavior depends on the file extension:

        - *.easy   : Load EasyHybrid project file. If a temporary file (*.easy~) exists,
                     prompt the user whether to load the most recent one.
        - *.easy~  : Load directly from the temporary project file.
        - Other    : Assume coordinate/system file and load into a new pDynamo system.

        Parameters
        ----------
        filename : str, optional
            Path to the file to load. If None, the method does nothing.
        """
        if not filename:
            return  # Nothing to load

        # Handle EasyHybrid project files (*.easy)
        if filename.endswith(".easy"):
            temp_file = filename + "~"
            if os.path.exists(temp_file):
                msg = ("There is a newer temporary file for the project you are loading.\n"
                       "Would you like to load the most current file?")
                dialog = SimpleDialog(self.main_session)
                yes_or_no = dialog.question(msg)

                target_file = temp_file if yes_or_no else filename
                self.main_session.p_session.load_easyhybrid_serialization_file(
                    target_file, tmp=yes_or_no
                )
            else:
                self.main_session.p_session.load_easyhybrid_serialization_file(filename)

        # Handle temporary EasyHybrid project files (*.easy~)
        elif filename.endswith(".easy~"):
            self.main_session.p_session.load_easyhybrid_serialization_file(filename)

        # Handle all other file types (assumed coordinate/system files)
        else:
            files = {"coordinates": filename}
            systemtype = 3  # Hardcoded system type (could be parameterized later)
            self.main_session.p_session.load_a_new_pDynamo_system_from_dict(files, systemtype)


    def show (self, obj = None, rep = 'lines', sele = None):
        """
        Display the specified atoms of a VisMol object.

        Parameters:
        - obj : str or None
            Name of the VisMol object to show, or 'all' to show all objects.
        - rep : str
            Representation type to show ('lines', 'sticks', 'spheres', etc.).
        - sele : list[int] or None
            List of atom indices to show. If None, all atoms are selected.
        """
        vismol_objects = []
        
        # Determine which VisMol objects to operate on
        if obj in self.vobject_names.keys():
            vismol_objects.append(self.vobject_names[obj])
        
        elif obj == 'all':
            vismol_objects = self.vobject_names.values()
        else:
            # If object not found, return False
            return False
        
        # Iterate over selected VisMol objects
        for vismol_object in vismol_objects:
            # If no specific atom selection provided, select all atoms
            
            #if sele:
            #    pass
            #else:
            #   size =  len(vismol_object.atoms)
            #   sele =  list(range(size))
            
            if not sele:
                sele = list(range(len(vismol_object.atoms)))
            
            # Create a new selection object
            selection = self.create_new_selection()
            
            # Select atoms by index in the VisMol object
            selection.selecting_by_indexes(vismol_object, sele, clear=True)
            
            # Show the selected atoms using the specified representation
            self.show_or_hide(rep_type  = rep,
                              selection = selection, 
                              show      = True )
        
    def hide (self, obj = None, rep = 'lines', sele = None):
        """
        Hide the specified atoms of a VisMol object.

        Parameters:
        - obj : str or None
            Name of the VisMol object to hide, or 'all' to hide all objects.
        - rep : str
            Representation type to hide ('lines', 'sticks', 'spheres', etc.).
        - sele : list[int] or None
            List of atom indices to hide. If None, all atoms are selected.
        """
        
        vismol_objects = []
                
        if obj in self.vobject_names.keys():
            vismol_objects.append(self.vobject_names[obj])
        
        elif obj == 'all':
            vismol_objects = self.vobject_names.values()
        else:
            return False
        
        for vismol_object in vismol_objects:
            if not sele:
                sele = list(range(len(vismol_object.atoms)))
            '''
            if sele:
                pass
            else:
               size =  len(vismol_object.atoms)
               sele =  list(range(size))
            '''
            
            selection = self.create_new_selection()
            
            selection.selecting_by_indexes(vismol_object, sele, clear=True)
            
            self.show_or_hide(rep_type  = rep,
                              selection = selection, 
                              show      = False )
   
    def active (self, obj = None):
        """
        Set the specified VisMol objects as active.

        Parameters:
        - obj : str or None
            Name of the VisMol object to activate, or 'all' to activate all objects.
        """
        vismol_objects = []
        
        if obj in self.vobject_names.keys():
            vismol_objects.append(self.vobject_names[obj])
        
        elif obj == 'all':
            vismol_objects = self.vobject_names.values()
        else:
            return False
        
        for vismol_object in vismol_objects:
            vismol_object.active = True
        #.change the treeview status (active) self.treestore[treeview_iter][6] = vobject.active
        self.main_session.main_treeview.refresh_number_of_frames()
        self.vm_glcore.queue_draw()
    
    def deactivate (self, obj = None):
        """
        Set the specified VisMol objects as inactive.

        Parameters:
        - obj : str or None
            Name of the VisMol object to deactivate, or 'all' to deactivate all objects.
        """
        vismol_objects = []
        
        if obj in self.vobject_names.keys():
            vismol_objects.append(self.vobject_names[obj])
        
        elif obj == 'all':
            vismol_objects = self.vobject_names.values()
        else:
            return False
        
        for vismol_object in vismol_objects:
            vismol_object.active = False
        #.change the treeview status (active) self.treestore[treeview_iter][6] = vobject.active
        self.main_session.main_treeview.refresh_number_of_frames()
        self.vm_glcore.queue_draw()
    
    def color ( obj = None, sele = None):
        """ Function doc """
    
    def set_system_color ( obj = None, sele = None):
        """ Function doc """

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
        
    
    #-------------------------------------------------------------------
    #                        restricted methods
    #-------------------------------------------------------------------
          
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
        """
        Add a VismolObject to the current session.

        This method registers the object, ensures unique naming, 
        updates internal dictionaries (atoms and objects), 
        and creates default visual representations if requested.

        Parameters
        ----------
        vismol_object : VismolObject
            The object to be added to the session.
        show_molecule : bool, optional, default=True
            Whether to create visual representations for the object.
        autocenter : bool, optional, default=True
            Whether to automatically center the view on the new object.

        Returns
        -------
        VismolObject
            The VismolObject added to the session (with updated index and name).
        """

        # Check if object with same ID already exists
        if vismol_object.index in self.vm_objects_dic:
            logger.warning(
                f"The VismolObject with id {vismol_object.index} already exists. "
                "The data will be overwritten."
            )

        # Assign new index and register object
        self.vm_objects_dic[self.vm_object_counter] = vismol_object
        vismol_object.index = self.vm_object_counter
        self.vm_object_counter += 1

        # Ensure object name has no spaces
        vismol_object.name = vismol_object.name.replace(" ", "_")

        # ------------------------------------------------------------------
        # Ensure unique naming for the VismolObject
        # ------------------------------------------------------------------
        name_parts = vismol_object.name.split("_")
        prefix_tag = name_parts[0]

        # If the name doesn't start with the object's e_id, prepend it
        if prefix_tag != str(vismol_object.e_id):
            vismol_object.name = f"{vismol_object.e_id}_{vismol_object.name}"

        # If the name is already in use, append suffixes until unique
        while vismol_object.name in self.vobject_names:
            vismol_object.name = f"{vismol_object.name}_X"

        # Register the object in the names dictionary
        self.vobject_names[vismol_object.name] = vismol_object

        # ------------------------------------------------------------------
        # Register all atoms by unique ID
        # ------------------------------------------------------------------
        for atom in vismol_object.atoms.values():
            self.atom_dic_id[atom.unique_id] = atom

        # ------------------------------------------------------------------
        # Create default molecular representations if requested
        # ------------------------------------------------------------------
        if show_molecule:
            vismol_object.create_representation(rep_type="lines")
            vismol_object.create_representation(rep_type="nonbonded")

            # Update restraints and custom colors
            self.main.p_session.update_restaint_representation(vismol_object.e_id)
            self.main.p_session._apply_custom_colors_to_vobject(vismol_object)

            # Adjust viewport
            if autocenter:
                self.vm_glcore.center_on_coordinates(vismol_object, vismol_object.mass_center)
            else:
                self.vm_glcore.queue_draw()

        return vismol_object
        
    def build_index_list_from_atom_selection(self, return_vobject=False):
        """
        Build an index list and residue dictionary from the current atom selection.

        This method extracts indices of selected atoms and organizes them by residue. 
        Optionally, it also returns the associated VismolObject.

        Parameters
        ----------
        return_vobject : bool, optional, default=False
            If True, the returned tuple also includes the VismolObject of the selection.

        Returns
        -------
        tuple
            If return_vobject is False:
                (index_list, residue_dict)
            If return_vobject is True:
                (index_list, residue_dict, vismol_object)

            - index_list : list of int
                List of atom indices (0-based).
            - residue_dict : dict[int, list[int]]
                Dictionary mapping residue indices to lists of atom indices.
            - vismol_object : VismolObject, optional
                The VismolObject associated with the selected atoms.
        
        Notes
        -----
        - Atom indices are converted to 0-based (`atom.index - 1`).
        - If any atom fails `check_selected_atom`, the method returns False immediately.
        """
        selection    = self.selections[self.current_selection]
        
        index_list   = []  
        residue_dict = {} 
        vismol_object = None
        
        for atom in selection.selected_atoms:
            # Validate atom with a custom check
            if not self.check_selected_atom(atom):
                return False  # Early exit if validation fails

            # Convert to 0-based index
            atom_index = atom.index - 1
            index_list.append(atom_index)

            # Keep track of the VismolObject (all atoms should belong to the same one)
            vismol_object = atom.vm_object

            # Group atom indices by residue
            residue_index = atom.residue.index
            residue_dict.setdefault(residue_index, []).append(atom_index)

        if return_vobject:
            return index_list, residue_dict, vismol_object
        return index_list, residue_dict


    def check_selected_atom(self, atom, dialog = True):
        """Checks if the selected atom belongs to the active pDynamo system."""
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
        """Sets the color of specific atoms in a Vismol object and updates the OpenGL representation."""
        for atom_index in indexes:
            vobject.atoms[atom_index].color = color    

        vobject._generate_color_vectors ( 
                                          True

                                            )
        try:
            self.vm_glcore.queue_draw()
        except:
            pass
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


    def set_color (self, symbol = 'C', color = [0.9, 0.9, 0.0], selection = None ):
        """
        Sets the color of all atoms of a given element symbol in a giver selection
        
            selection is a vismol selection -->
        
            if selection = None, selection = the current selection.
        
        """
        if selection:
            pass
        else:
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
        """Creates a new selection object."""
        return VMSele(self)


    def forward_frame(self):
        """Moves the visualization forward by one frame."""

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
        """Moves the visualization backward by one frame."""
        frame = self.frame
        if self.frame - 1 >= 0:
            frame -= 1
        else:
            frame  = 0
        self.set_frame(frame=frame)


    def set_frame(self, frame=0):
        """Sets the current frame for all Vismol objects and updates representations."""
        assert frame >= 0
        self.frame = np.uint32(frame)
        self.vm_glcore.updated_coords = True
        
        if self.picking_selection_mode:
            self.picking_selections.update_pki_pkj_rep_coordinates()
            self.picking_selections.print_pk_distances()    
        
        #self.main_session.label_frame.set_text('Frame Number: ' + str(frame)) 
        self.vm_widget.queue_draw()
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
        
        Determines the smallest Cartesian box that contains the selected atoms.

        This function establishes the "inner box" (smallest box in Cartesian space) 
        encompassing all selected atoms in the provided selection. It returns the 
        involved VisMol objects, the minimum and maximum coordinates of the box, 
        and a dictionary mapping each VisMol object's index to its selected atom indices.


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



        Parameters:
        - selection : VisMolSelection
            The selection object containing atoms to consider.
        - grid_size : float
            Grid size (currently not used in computation but could be used for further processing).

        Returns:
        - vobjects : list
            List of involved VisMol objects.
        - grid_min : list[float]
            Minimum coordinates [x_min, y_min, z_min] of the inner box.
        - grid_max : list[float]
            Maximum coordinates [x_max, y_max, z_max] of the inner box.
        - selected_indexes_dict : dict
            Dictionary mapping each vobject.index to the list of selected atom indices.
            Example: {vobject_index: [atom_index_1, atom_index_2, ...]}

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
        
        # Iterate over each selected atom
        for atom in selection.selected_atoms:
            '''checks if the selected atoms belong to the active project'''
            true_or_false = self.check_selected_atom(atom, dialog = True)
            
            if true_or_false:
                xyz = atom.coords()
                vobject = atom.vm_object

                # Add vobject to list if not already included
                if vobject in vobjects:
                    pass
                else:
                    vobjects.append(vobject)
                # Collect coordinates
                i.append(xyz[0])
                j.append(xyz[1])
                k.append(xyz[2])
                
                # Update selected indexes dictionary
                if vobject.index in selected_indexes_dict.keys():
                    selected_indexes_dict[vobject.index].append(atom.index-1)
                else:
                    selected_indexes_dict[vobject.index]=[]
                    selected_indexes_dict[vobject.index].append(atom.index-1)

            else:
                print('invalid selection list')
                return False
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Compute the min and max coordinates for the inner box
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
