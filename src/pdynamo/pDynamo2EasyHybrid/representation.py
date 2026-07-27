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
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, GLib
import multiprocessing

import glob, math, os, os.path, sys, shutil
import pickle
import threading
from util.file_parser import read_MOL2  
from util.file_parser import read_SIMPLE_txt  
from util.file_parser import read_MOPAC_aux  

from datetime import date
import time

import numpy as np
import copy

import random
import string

from pprint import pprint

#VISMOL_HOME = os.environ.get('VISMOL_HOME')

#path fo the core python files on your machine
#sys.path.append(os.path.join(VISMOL_HOME,"easyhybrid/pDynamoMethods") )
#sys.path.append(os.path.join(VISMOL_HOME,"easyhybrid/gui"))

#from LogFile import LogFileReader

#from gEngine.vismol_object import EVismolObject

#from vismol.model.atom import Atom
from vismol.model.residue import Residue
from vismol.model.chain import Chain
from vismol.core.vismol_object import VismolObject
from vismol.model.atom import Atom
#print ('\n\n\n\n\n\nATOM',Atom,'\n\n\n\n\n\nATOM')
from logging import getLogger
logger = getLogger(__name__)

#---------------------------------------
from pBabel                    import*                                     
from pCore                     import*  
#---------------------------------------
from pMolecule                 import*                              
from pMolecule.MMModel         import*
from pMolecule.NBModel         import*                                     
from pMolecule.QCModel         import*
#---------------------------------------
from pScientific               import*                                     
from pScientific.Arrays        import*                                     
from pScientific.Geometry3     import*                                     
from pScientific.RandomNumbers import*                                     
from pScientific.Statistics    import*
from pScientific.Symmetry      import*
#---------------------------------------                              
from pSimulation               import*
#---------------------------------------


import numpy as np
#from vismol.model.molecular_properties import ATOM_TYPES
from vismol.libgl.representations import DashedLinesRepresentation

from util.colorpalette import CUSTOM_COLOR_PALETTE

from pdynamo.p_methods import GeometryOptimization
from pdynamo.p_methods import RelaxedSurfaceScan
from pdynamo.p_methods import AdvancedRelaxedSurfaceScan
from pdynamo.p_methods import MolecularDynamics
from pdynamo.p_methods import ChainOfStatesOptimizePath
from pdynamo.p_methods import NormalModes
from pdynamo.p_methods import EnergyCalculation
from pdynamo.p_methods import EnergyRefinement
from pdynamo.p_methods import UmbrellaSampling

from pdynamo.p_methods import WHAMAnalysis
from pdynamo.LogFileWriter import LogFileReader

from gui.windows.setup.windows_and_dialogs import call_message_dialog

class ModifyRepInVismol:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass

    def clean_QC_representation_to_vobject (self, system_id = None):
        """ Function doc """
        if system_id:
            system = self.psystem[system_id]
        else:
            system = self.psystem[self.active_id]
        
        if system.qcModel:
            qc_table = list(system.qcState.pureQCAtoms)
            for vismol_object in self.vm_session.vm_objects_dic.values():
                if vismol_object.e_id == system.e_id:
                    if getattr(vismol_object, 'is_surface', False):
                        pass
                    else:
                        
                        selection = self.vm_session.create_new_selection()
                        selection.selecting_by_indexes(vismol_object, qc_table, clear=True)
                        self.vm_session.show_or_hide(rep_type = 'lines',selection= selection, show = True )
                        self.vm_session.show_or_hide(rep_type = 'dynamic',selection= selection , show = False )
                        
                        #self.vm_session.show_or_hide(rep_type = 'stick',selection= selection , show = False )
                        for atom in vismol_object.atoms.values():
                            atom.spheres = False
                            atom.sticks = False
                        
                        vismol_object.representations['spheres'] = None
                        vismol_object.representations['sticks'] = None
                        #self._apply_QC_representation_to_vobject   (system_id = None, vismol_object = vismol_object)                   
            

    def _apply_fixed_representation_to_vobject (self, system_id = None, vismol_object = None):
        """ Function doc """
        
        if system_id:
            pass
        else:
            system_id = vismol_object.e_id
        
        self.get_fixed_table_from_pdynamo_system(system_id)
        
        
        #if self.psystem[system_id].freeAtoms is None:
        #    pass
        #else:
        #    if self.psystem[system_id].e_fixed_table == []:
        #        freeAtoms = self.psystem[system_id].freeAtoms
        #        freeAtoms                              = Selection.FromIterable (freeAtoms)
        #        selection_fixed                        = freeAtoms.Complement( upperBound = len (self.psystem[system_id].atoms ) )
        #        self.psystem[system_id].e_fixed_table  = list(selection_fixed)
        


        
        indexes = np.array(self.psystem[system_id].e_fixed_table, dtype=np.int32)    
        color   = np.array(self.fixed_color, dtype=np.float32)    
        self.vm_session.set_color_by_index(vobject       = vismol_object , 
                                           indexes       = indexes, 
                                           color         = color)

    def _apply_QC_representation_to_vobject (self, system_id = None, vismol_object = None, static = True):
        """ Function doc """
        if vismol_object:
            pass
        else:
            return False
            
        
        if system_id:
            pass
        else:
            system_id = vismol_object.e_id
        

        if self.psystem[system_id].qcModel:
            
            #if self.psystem[system_id].mmModel:
            self.psystem[system_id].e_qc_table = list(self.psystem[system_id].qcState.pureQCAtoms)
          
            '''
                This part of the code identifies which atoms in the QC 
            region are connected to atoms in the MM region (MM_QC_atoms). 
            Then build a list containing only atoms from the QC region 
            with no connection to the MM part. the line representation 
            will be erased for the atoms in this new list.
            '''
            #------------------------------------------------------------------
            MM_QC_atoms = []
            
            #print(self.psystem[system_id].e_qc_table)
            #psystem = self.p_session.psystem[self.p_session.active_id]
            if len(self.psystem[system_id].atoms) == len(self.psystem[system_id].e_qc_table ):
                pass
            else:
            
                for index in self.psystem[system_id].e_qc_table:
                    for bond in vismol_object.atoms[index].bonds:
                        a1 = bond.atom_i.index -1
                        a2 = bond.atom_j.index-1
                        #print (index, bond.atom_i.index-1, bond.atom_j.index-1 )
                        
                        if a1 in self.psystem[system_id].e_qc_table :
                            pass
                        else:
                            MM_QC_atoms.append(index)
                        
                        
                        if a2 in self.psystem[system_id].e_qc_table :
                            pass
                        else:
                            MM_QC_atoms.append(index)
                
            
            # list containing only atoms from the QC region with no connection to the MM part
            delete_lines = [] 
            
            for index in self.psystem[system_id].e_qc_table:
                if index in MM_QC_atoms:
                    #print('if index',index)
                    pass
                else:
                    delete_lines.append(index)
                    #print('else index',index)
            #------------------------------------------------------------------
            
            
            for atom in vismol_object.atoms.values():
                atom.spheres = False
                atom.sticks = False
                    
                
            selection = self.vm_session.create_new_selection()
            
            #print(self.psystem[system_id].e_qc_table)
            #e_qc_table = self.psystem[system_id].e_qc_table[]
            #selection.selecting_by_indexes(vismol_object,  e_qc_table, clear=True)
            
            selection.selecting_by_indexes(vismol_object, self.psystem[system_id].e_qc_table, clear=True)
            
            '''
            for _id , vismol_object in self.vm_session.vm_objects_dic.items():
                selection.selecting_by_indexes(vismol_object, self.psystem[system_id].e_qc_table, clear=False)
            
            #'''

            self.vm_session.show_or_hide (rep_type = 'spheres',selection= selection ,  show = True )
            #self.vm_session.show_or_hide (rep_type = 'spheres', show = True )

            if static:

                selection2 = self.vm_session.create_new_selection()
                selection2.selecting_by_indexes(vismol_object, delete_lines, clear=True)

                #self.vm_session.show_or_hide(rep_type = 'sticks',selection= selection , show = True )
                self.vm_session.show_or_hide(rep_type = 'dynamic',selection= selection , show = True )
                self.vm_session.show_or_hide(rep_type = 'lines',selection= selection2 , show = False )
                
                #self.vm_session.show_or_hide(rep_type = 'sticks' , show = True )

            else:
                pass
                #self.vm_session.show_or_hide(_type = 'dynamic_bonds' , vobject = vobject, selection_table = self.systems[system_id]['qc_table'] , show = True )

            for atom in vismol_object.atoms.values():
                atom.selected  =  False
                atom.vm_object.selected_atom_ids.discard(atom.atom_id)
            
    def  _apply_custom_colors_to_vobject (self,vismol_object = None):
        """ Function doc """
        #return False
        #e_qc_table = list(range(100))
        
        system = self.psystem[vismol_object.e_id]
        for custom_color in system.e_custom_colors:
            indexes = custom_color['indexes']
            color   = custom_color['color']
            
            selection  = self.vm_session.create_new_selection()
            selection.selecting_by_indexes(vismol_object, 
                                           indexes, 
                                           clear=True)

            self.vm_session.set_color (symbol = 'C', 
                                       color = color, 
                                       selection = selection )
