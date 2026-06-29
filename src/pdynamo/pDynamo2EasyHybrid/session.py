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

# --- imports entre modulos adicionados na refatoracao ---
from pdynamo.pDynamo2EasyHybrid.simulations_mixin import pSimulations
from pdynamo.pDynamo2EasyHybrid.analysis_mixin import pAnalysis
from pdynamo.pDynamo2EasyHybrid.representation import ModifyRepInVismol
from pdynamo.pDynamo2EasyHybrid.io_data import LoadAndSaveData
from pdynamo.pDynamo2EasyHybrid.import_trajectory import EasyHybridImportTrajectory
from pdynamo.pDynamo2EasyHybrid.restraints import Restraints
from pdynamo.pDynamo2EasyHybrid.helpers import Atom, generate_random_code, export_special_PDB

class pDynamoSession (pSimulations, pAnalysis, ModifyRepInVismol, LoadAndSaveData, EasyHybridImportTrajectory, Restraints):
    """ Class doc """
    
    def __init__ (self, vm_session = None):
        """ Class initialiser """
        self.vm_session  = vm_session
        self.main        = self.vm_session.main
        self.name        = 'p_session'
        
        self.nbModel_default         = NBModelCutOff.WithDefaults ( )
        self.fixed_color             = [0.5, 0.5, 0.5]
        self.pdynamo_distance_safety = 0.5
        
        '''self.active_id is the attribute that tells which 
        system is active for calculations in pdynamo 
        (always an integer value)'''
        self.active_id       = 0
        
        
        '''Now we can open more than one pdynamo system. 
        Each new system loaded into memory is stored in 
        the form of a dictionary, which has an int as 
        an access key.'''
        self.psystem =  {
                         0:None
                        }
        self.psystem_name_list    = []
        self.psystem_name_counter = 1
        
        self.counter      = 0
        self.color_palette_counter = 0
        
        '''This attribute is checked by the session before 
        closing it. If it is set to False, the session is 
        closed without question. If it is set to True, the 
        GUI asks the user if they want to save the session.
        "'''
        self.random_code = generate_random_code(10)
        self.changed = False
        
        '''
        There is now a self.process_pool attribute (a dictionary),
        which should store all generated processes, identified by the e_id
        of each system. This allows multiple systems to run processes
        simultaneously and be canceled if necessary.
        
        self.process_pool = {e_id_1 :  process_object_1, 
                             e_id_2 :  process_object_2,}
        
        #'''
        self.process_pool = {}
        self.GLib_monitor =None
        
    
    def set_active(self, system_e_id = None):
        """ Function doc """
        if system_e_id != None:
            self.active_id = system_e_id

        
    def get_system (self, index = None):
        """ Function doc """
        if index != None:
            return self.psystem[index]
        else:
            return self.psystem[self.active_id]
        
    
    def _check_name (self, name):
        """ Function doc """
        
        stop = False
        original_name = name
        counter = 0
        
        if name:
            
            while stop == False:
                if name in self.psystem_name_list:
                    name = original_name+' '+str(counter)
                    #print('_check_name', name)
                else:
                    stop = True
                
        else:
            name = 'new system'
        
        return name
        
    
    def restart (self):
        """ Function doc """
        #-------------------------------------------
        self.active_id = 0
        self.psystem   =  {
                        0:None
                        }
        self.psystem_name_list    = []
        self.psystem_name_counter = 1
        
        self.counter      = 0
        self.color_palette_counter = 0
        #-------------------------------------------
    
    
    def get_hamiltonian_type (self, e_id = None):
        """ Function doc """
        psystem = self.psystem[e_id]
        
        qc = False
        mm = False
        
        if psystem.qcModel:
            qc = 'QC'
            #hamiltonian   = getattr(psystem.qcModel, 'hamiltonian', False)
            #if hamiltonian:
            #    pass
            #else:
            #    try:
            #        itens = psystem.qcModel.SummaryItems()
            #        #print(itens)
            #        hamiltonian = itens[0][0]
            #    except:
            #        hamiltonian = 'external'    
            
        if psystem.mmModel:
            mm = 'MM'
            #forceField = psystem.mmModel.forceField
            #string += 'force field: {}    '.format(forceField)
            #
            #if psystem.nbModel:
            #    nbmodel = psystem.mmModel.forceField
            #    string += 'nbModel: True    '
            #    
            #    summary_items = psystem.nbModel.SummaryItems()
            #    
            #
            #else:
            #    string += 'nbModel: False    '
        if qc and mm:
            return 'QC/MM'
        elif qc:
            return 'QC'
        elif mm:
            return 'MM'
        else:
            pass
            
    
    def load_a_new_pDynamo_system_from_dict (self, input_files = {}, system_type = 0, name = None, tag = None, color = None, working_folder = None):
        """ Function doc """
        
        system = None 
        if system_type == 0:
            system              = ImportSystem       ( input_files['amber_prmtop'] )
            self.main.bottom_notebook.status_teeview_add_new_item(message = 'loading file:  {} '.format( input_files['amber_prmtop']), system = None)
            system.coordinates3 = ImportCoordinates3 ( input_files['coordinates'] )
            self.main.bottom_notebook.status_teeview_add_new_item(message = 'loading file:  {} '.format( input_files['coordinates']), system = None)

            self.define_NBModel(_type = 1, system = system)                      
        
        elif system_type == 1:
            parameters          = CHARMMParameterFileReader.PathsToParameters (input_files['charmm_par'])
            self.main.bottom_notebook.status_teeview_add_new_item(message = 'loading file:  {} '.format( input_files['charmm_par']), system = None)
            
            system              = ImportSystem       ( input_files['charmm_psf'] , isXPLOR = True, parameters = parameters )
            self.main.bottom_notebook.status_teeview_add_new_item(message = 'loading file:  {} '.format( input_files['charmm_psf']), system = None)
            
            system.coordinates3 = ImportCoordinates3 ( input_files['coordinates'] )
            self.main.bottom_notebook.status_teeview_add_new_item(message = 'loading file:  {} '.format( input_files['coordinates']), system = None)
            
            self.define_NBModel(_type = 1, system = system)        
        
        elif system_type == 2:
            mmModel        = MMModelOPLS.WithParameterSet ( input_files['prm_folder'] )       
            self.main.bottom_notebook.status_teeview_add_new_item(message = 'loading MMModel:  {} '.format( input_files['prm_folder']), system = None)
            
            system         = ImportSystem ( input_files['coordinates'])
            self.main.bottom_notebook.status_teeview_add_new_item(message = 'loading file:  {} '.format(input_files['coordinates']), system = None)
            
            system.DefineMMModel ( mmModel )
            self.define_NBModel(_type = 1, system = system)          
        
        elif system_type == 5:
            mmModel        = MMModelDYFF.WithParameterSet ( input_files['prm_folder'] )       
            self.main.bottom_notebook.status_teeview_add_new_item(message = 'loading MMModel:  {} '.format( input_files['prm_folder']), system = None)
            
            system         = ImportSystem ( input_files['coordinates'])
            self.main.bottom_notebook.status_teeview_add_new_item(message = 'loading file:  {} '.format(input_files['coordinates']), system = None)
            #system.BondsFromCoordinates3 ( )
            system.DefineMMModel ( mmModel )
            
            self.define_NBModel(_type = 1, system = system)          
            if input_files['charges']:
                print('\nGetting atomic charges from mol2 file!\n')
                self.main.bottom_notebook.status_teeview_add_new_item(message = 'Getting atomic charges from mol2 file!', system = None)
                for index, chg in enumerate(input_files['charges']):
                    system.mmState.charges[index] = chg
                
            #for i, atom in enumerate(system.atoms):
            #    print(i, atom.atomicNumber, atom.label)
            
        elif system_type == 3 or system_type == 4 :
            system = ImportSystem (input_files['coordinates'])
            self.main.bottom_notebook.status_teeview_add_new_item(message = 'loading file:  {} '.format(input_files['coordinates']), system = None)
        else:
            pass
        system.Summary()
        
        
        
        
        
        if working_folder:
            system.e_working_folder = working_folder
        else:
            pass
            
        ''''
        If the name is too long we have problems with the "system 
        comboboxes", they get huge and the use of the interface is impaired.
        '''
        if name == None:
            name = getattr (system, 'label', None)
            tag  = getattr (system,'e_tag', None)
        if len(name) > 70:
            name = name[-70:]
        
        
        self.remove_PES_data (system = system )
        
        '''Storing the system source files'''
        system.e_input_files = input_files
        system = self.add_new_system_to_psession (system = system,  name = name, tag = tag, color = color, changed = True )

        '''-----------------------------------------------------------------'''

        # Add the system to the main treeview
        self.main.main_treeview.add_new_system_to_treeview (system)

        # Determine the force field used
        # sys_type = {0: 'AMBER', 1: 'CHARMM'}        
        ff  =  getattr(system.mmModel, 'forceField', "None")
        
        
        # Add information about the new system to the status treeview
        self.main.bottom_notebook.status_teeview_add_new_item(message = 'New System:  {} ({}) - Force Field:  {}'.format(system.label, system.e_tag, ff), system =system)
        
        # Add the system as a vismol object to the easyhybrid session
        self._add_vismol_object_to_easyhybrid_session (system, True)  
        
        
        #self.main.refresh_active_system_liststore()
        #self.main.refresh_system_liststore ()
        ''' '''
   
    
    def _add_vismol_object_to_easyhybrid_session (self, system, show_molecule=True, name = 'new_coords', key6 = None):
        """ Function doc """
        # Create a VisMol object from the given pDynamo system
        vm_object = self._build_vobject_from_pdynamo_system ( system = system, name = name ) 
        
        if key6:
            vm_object.key6 = key6
        
        # Add the VisMol object to the VisMol session
        self.vm_session._add_vismol_object(vm_object, show_molecule=True)
        
        # Add the VisMol object to the main treeview
        self.main.main_treeview.add_vismol_object_to_treeview(vm_object)
        
        # Add the VisMol object to the vobject liststore dictionary
        self.main.add_vobject_to_vobject_liststore_dict(vm_object)
        
        # Apply fixed representation to the VisMol object
        self._apply_fixed_representation_to_vobject(vismol_object =vm_object)
        
        # Apply QC representation to the VisMol object
        self._apply_QC_representation_to_vobject(vismol_object =vm_object)
        
        # Refresh the widgets in the main window
        self.main.refresh_widgets()
        
        # Return the VisMol object
        return vm_object

    
    def add_new_system_to_psession (self, 
                                          system         = None            ,
                                          name           = 'pDynamo system',
                                          tag            = None            ,
                                          working_folder = None            ,
                                          color          = None            ,
                                          changed        = False ):
        """ Function doc """  
    
        
        # - - - - - - - Name and Tag - - - - - - - -
        if name:
            system.label = name
        else:
            pass
        
        if tag:
            system.e_tag = tag
        else:
            system.e_tag ='MolSys'
        # - - - - - - - - - - - - - - - - - - - - - -
        
        
        

        # - - - - - - - - - - - - - Working Folder - - - - - - - - - - - - - 
        # . When no name is provided.
        if working_folder is None:
            # . Checks whether the pkl file (hence system object) contains an associated workbook
            is_wf_already_set = getattr ( system, "e_working_folder", False )            
            
            if is_wf_already_set is False:
                # . If not, the pdynamo scratch folder will be used
                try:  
                    if os.path.isdir(self.vm_session.gl_parameters["tmp_path"]):
                        system.e_working_folder = self.vm_session.gl_parameters["tmp_path"]
                    else:
                        system.e_working_folder = os.environ.get('PDYNAMO3_SCRATCH')
                except:
                    system.e_working_folder =  os.environ.get('HOME')
            
            else:
                # . If there is a path to the working folder, verify that the folder exists.
                isExist = os.path.exists(is_wf_already_set)
                if isExist:
                    pass
                else:
                    system.e_working_folder = os.environ.get('PDYNAMO3_SCRATCH')
        else:
            # . When a working folder is provided by easyhybrid.
            system.e_working_folder = working_folder
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        #                            CHARGES
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        try:
            system.e_charges_backup           = list(system.AtomicCharges()).copy()
        except:
            system.e_charges_backup           = []
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        

        
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        system.e_active                   = False
        system.e_date                     = date.today()               # Time     
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        

        
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        #                                  COLORS
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if color is not None:
            system.e_color_palette            = self.vm_session.periodic_table.get_color_palette(custom = {"C":color}) #COLOR_PALETTE[0].copy() 
            #system.e_color_palette['C']       = color
        #'''
        else:
            e_color_palette = getattr(system, 'e_color_palette', None)
            if type(e_color_palette) == dict:
                pass
            else:
                #system.e_color_palette            = COLOR_PALETTE[self.color_palette_counter].copy() 
                custom = CUSTOM_COLOR_PALETTE[self.color_palette_counter].copy() 
                system.e_color_palette            = self.vm_session.periodic_table.get_color_palette(custom = custom) 
                #When the number of available colors runs out, we have to reset the counter 
                self.color_palette_counter        += 1
                if self.color_palette_counter > len(CUSTOM_COLOR_PALETTE.keys())-1:
                    self.color_palette_counter = 0
                else:
                    pass
        #'''
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        
        
        
        
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        system.e_active                   = False          
        system.e_bonds                    = None           
        system.e_sequence                 = None           

        #system.e_QCMM                     = False
        system.e_qc_table                 = []           
        system.e_qc_residue_table         = {}             #yes, it's a dict           
        system.e_fixed_table              = []             
        
        system.e_color_segments           = {}   # {1: [[idx1, idx2, idx3, ... ], [red, green, blue]]}
        
        
        if getattr ( system, "e_selections", False ):
            pass
            #system.e_selections               = {}             
        else:
            system.e_selections               = {}             
        
        
        if getattr ( system, "e_custom_colors", False ):
            pass
            #system.e_selections               = []             
        else:
            system.e_custom_colors             = [] 
        
        
        system.e_selections_counter       = 0             
        system.e_vm_objects               = {}                  
        
        
        if getattr ( system, "e_job_history", False ):
            pass
        else:
            system.e_job_history              = {}
        
        
        if getattr ( system, "e_logfile_data", False ):
            pass
        else:
            system.e_logfile_data             = {}             # <--- vobject_id : [data0, data1, data2], ...]  obs: each "data" is a dictionary 
        
        
        
        
        if getattr ( system, "e_step_counter", False ):
            pass
        else:
            system.e_step_counter             = 0              
        
        
        
        if getattr ( system, "e_restraints_dict", False ):
            pass
        else:
            system.e_restraint_counter        = 0
            system.e_restraints_dict          = {
                                                  #[True, 
                                                  # rest_name ,
                                                  # 'distance', 
                                                  # [parameters['atom1'],parameters['atom2']], 
                                                  # parameters['distance'],  
                                                  # parameters['force_constant']] 
                                                
                                                }
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        if getattr ( system, "e_annotations", False ):
            pass
        else:
            system.e_annotations = ''

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        
        
        '''
        Now each pdynamo system object and vismol object has a 
        "treeview_iter" attribute, which is a reference to access 
        the treeview elements
        
        the method "add_new_system_to_treeview" 
        
        '''
        system.e_treeview_iter    = None  #the method "add_new_system_to_treeview" modify this attribute 
        system.e_liststore_iter   = None  #the method "add_new_system_to_treeview" modify this attribute
        
        system.e_id               =  self.counter
        self.psystem[system.e_id] =  system 
        self.changed              = changed  
        self.active_id            = self.counter  
        self.counter    += 1
        return system


    def _get_sequence_from_pdynamo_system (self, system = None):
        """ Function doc """
        
        if system:
            sequence = getattr ( system, "sequence", None )
            
            if sequence is None: 
                '''If there is no sequence, easyhybrid will look for 
                the sequence_from_mol2 attribute that comes from the 
                modified MOL2FileReader.'''
                is_from_mol2 = getattr(system, 'sequence_from_mol2', False)
                if is_from_mol2:
                    sequence = system.sequence_from_mol2
                else:
                    sequence = Sequence.FromAtoms (system.atoms, componentLabel = "UNK.1" )
            else:
                #print('if sequence is None: - else')
                pass
        else:
            sequence = None
        
        return sequence


    def _get_atom_info_from_pdynamo_atom_obj (self, sequence = None, atom = None, is_from_mol2 = False):
        """
        To extract information from atom objects, 
        belonging to pdynamo, and to  organize it as a list
        that will be delivered later to build the 
        vismolObj
        
        Now chemical symbols are defined based on atomic numbers 
        (pdynamo's task of providing atomic numbers).
        
        """


        
        if is_from_mol2:
            #':DFX.1:C1'
            resdata = sequence[atom.index].split(':')
            resdata = resdata[1].split('.')
            
            resName = resdata[0]
            resSeq  = resdata[1]
            chainID = 'X'
        else:
            
            entityLabel = atom.parent.parent.label
            #entityLabel = atom.label
            useSegmentEntityLabels = False
            if useSegmentEntityLabels:
                chainID = ""
                segID   = entityLabel[0:4]
            else:
                chainID = entityLabel[0:1]
                segID   = ""            


            if sequence:
                resName, resSeq, iCode = sequence.ParseLabel ( atom.parent.label, fields = 3 )
            else:
                sequence = None
        
        #print(atom.index, atom.label,resName, resSeq)#, iCode,chainID, segID,  atom.atomicNumber, atom.connections)#, xyz[0], xyz[1], xyz[2] )
        
        index        = atom.index+1
        at_name      = atom.label
        at_resi      = int(resSeq)
        at_resn      = resName
        at_ch        = chainID
        
        
        '''
        Now chemical symbols are defined based on atomic numbers 
        (pdynamo's task of providing atomic numbers).
        '''
        try:
            at_symbol    = self.vm_session.periodic_table.check_symbol (symbol =None, number = atom.atomicNumber)
        except:
            at_symbol = "X"
        
        at_occup     = 0.0
        at_bfactor   = 0.0
        at_charge    = 0.0
        atom         = {
              'index'      : index      , 
              'name'       : at_name    , 
              'resi'       : at_resi    , 
              'resn'       : at_resn    , 
              'chain'      : at_ch      , 
              'symbol'     : at_symbol  , 
              'occupancy'  : at_occup   , 
              'bfactor'    : at_bfactor , 
              'charge'     : at_charge   
              }
        #print(at_symbol)
        return atom


    def _get_qc_table_from_pDynamo_system(self, system):
        system.e_qc_table = list(system.qcState.pureQCAtoms)

    
    def interpolate_frames_of_a_vobject (self,system, vobject):
        """ Function doc """
        #for 
        #p_coords = ImportCoordinates3 ( parameters['data_path'] )
        #v_coords = self._convert_pdynamo_coords_to_vismol(p_coords)
        #print (parameters['vobject'].frames)
        #coords = np.vstack((coords, f))
        
        
        atom_qtty = len(system.atoms.items)
        size = len(vobject.frames)
        coords    = np.empty([(size*2)-1, atom_qtty, 3], dtype=np.float32)
        i = 0
        for index in range(0, len(vobject.frames)-1,2):
            frame1 = self.get_coordinates_from_vobject ( vobject, index)
            frame2 = self.get_coordinates_from_vobject ( vobject, index+1)
            newframe = self._interpolate_frames (frame1, frame2)
        
            j = 0
            for xyz in frame1:
                x = np.float32(xyz[0])
                y = np.float32(xyz[1])
                z = np.float32(xyz[2])
                coords[i,j,:] = x, y, z
                j += 1
            
            i += 1
            
            j = 0
            for xyz in newframe:
                x = np.float32(xyz[0])
                y = np.float32(xyz[1])
                z = np.float32(xyz[2])
                coords[i,j,:] = x, y, z
                j += 1
            
            i += 1
            
            j = 0
            for xyz in frame2:
                x = np.float32(xyz[0])
                y = np.float32(xyz[1])
                z = np.float32(xyz[2])
                coords[i,j,:] = x, y, z
                j += 1
            i += 1
        
        vobject.frames = coords
        
        return vobject
        
        #vobject.frames = np.vstack((parameters['vobject'].frames, v_coords))
    
    def _interpolate_frames (self, frame1, frame2):
        """ Function doc """
        newframe = []
        
        for index in range(len(frame1)):
            xyz1 = frame1[index]
            xyz2 = frame2[index]
            
            x = xyz2[0] - xyz1[0]
            y = xyz2[1] - xyz1[1]
            z = xyz2[2] - xyz1[2]
            newframe.append([x,y,z])
        return newframe
            

    def get_coordinates_from_vobject (self, vobject = None, frame = -1):
        """ Function doc """
        coords = []
        for i, atom in vobject.atoms.items():
            #print(i, atom)
            xyz = atom.coords(frame = frame)
            coords.append(xyz)
        return coords


    def set_psystem_coordinates_from_vobject (self, vobject = None, system_id =  None, frame = -1):
        """ Function doc """
        if system_id != None:
            pass
        else:
            system_id = self.active_id
        
        ##print('Loading coordinates from', vobject.name, 'frame', frame)
        ##print(vobject.atoms)
        for i, atom in vobject.atoms.items():
            #print(i, atom)
            xyz = atom.coords(frame = frame)
            self.psystem[system_id].coordinates3[i][0] = xyz[0]
            self.psystem[system_id].coordinates3[i][1] = xyz[1]
            self.psystem[system_id].coordinates3[i][2] = xyz[2]
        
        #self.psystem[system_id].Summary()
        #print(vobject.cell_coordinates)
        #print(vobject.cell_coordinates)
        if self.psystem[system_id].symmetry:
            a = vobject.cell_coordinates[frame][7][0]
            b = vobject.cell_coordinates[frame][7][1]
            c = vobject.cell_coordinates[frame][7][2]
            
            alpha = vobject.cell_parameters['alpha']
            beta  = vobject.cell_parameters['beta']
            gamma = vobject.cell_parameters['gamma']
            
            print('setting cell parameters vobj to psys:')
            print('a = ', a)
            print('b = ', b)
            print('c = ', c)
            print('alpha = ', alpha)
            print('beta  = ', beta )
            print('gamma = ', gamma)
            
            self.psystem[system_id].symmetryParameters.SetCrystalParameters(a,b,c, alpha ,beta, gamma)
            
        
    def get_fixed_atoms_from_system (self, system):
        """
        Get a list of fixed atoms from the given system object.

        Args:
            system (YourSystemClass): The system object containing atoms and freeAtoms attribute.

        Returns:
            list: A list of fixed atoms (atoms not present in freeAtoms).
        """
        if system.freeAtoms is None:
            pass
        
        else:
            freeAtoms = system.freeAtoms
            freeAtoms = Selection.FromIterable (freeAtoms)
            selection_fixed = freeAtoms.Complement( upperBound = len (system.atoms ) )
            return list(selection_fixed)
    
    
    def get_fixed_table_from_pdynamo_system (self, system_id = None):
        """ Function doc """
        if system_id != None:
            pass
        else:
            system_id = self.active_id



        if self.psystem[system_id].freeAtoms is None:
            pass
        else:
            if self.psystem[system_id].e_fixed_table == []:
                freeAtoms = self.psystem[system_id].freeAtoms
                freeAtoms                              = Selection.FromIterable (freeAtoms)
                selection_fixed                        = freeAtoms.Complement( upperBound = len (self.psystem[system_id].atoms ) )
                self.psystem[system_id].e_fixed_table  = list(selection_fixed)
        
        return self.psystem[system_id].e_fixed_table
    
    
    def get_output_filename_from_system (self, stype = None):
        """ Function doc """
        psystem = self.psystem[self.active_id]
            
        name    = psystem.label
        size    = len(psystem.atoms)
        string = 'system: {}    atoms: {}    '.format(name, size)
        
        
        qc_string = ''
        
        if psystem.qcModel:
            hamiltonian   = getattr(psystem.qcModel, 'hamiltonian', False)
            if hamiltonian:
                pass
            else:
                try:
                    itens = psystem.qcModel.SummaryItems()
                    #print(itens)
                    hamiltonian = itens[0][0]
                except:
                    hamiltonian = 'external'
                
                
                
            
            n_QC_atoms    = len(list(psystem.qcState.pureQCAtoms))
           
            summary_items = psystem.electronicState.SummaryItems()
            
            qc_string += '{}_QC{}_'.format(  hamiltonian, n_QC_atoms)
        else:
            pass
            
            #string += 'hamiltonian: {}    QC atoms: {}    QC charge: {}    spin multiplicity {}    '.format(  hamiltonian, 
            #                                                                                                  n_QC_atoms,
            #                                                                                                  summary_items[1][1],
            #                                                                                                  summary_items[2][1],
            #                                                                                                 )
        
        
        mm_string = ''    
        n_fixed_atoms = len(psystem.e_fixed_table)
        string += 'fixed atoms: {}    '.format(n_fixed_atoms)
        
        if psystem.mmModel:
            forceField = psystem.mmModel.forceField
            string += 'force field: {}    '.format(forceField)
        
            if psystem.nbModel:
                nbmodel = psystem.mmModel.forceField
                string += 'nbModel: True    '
                
                summary_items = psystem.nbModel.SummaryItems()
                
            
            else:
                string += 'nbModel: False    '
            
            
            if psystem.symmetry:
                #nbmodel = psystem.mmModel.forceField
                string += 'PBC: True    symmetry: {}    '.format( psystem.symmetry.crystalSystem.label)
                #print(psystem.symmetry)
                #print(psystem.symmetryParameters)
                #summary_items = psystem.nbModel.SummaryItems()
                
            
            else:
                string += 'PBC: False    '
        
            #mm_string += '{}_F{}'.format(forceField, n_fixed_atoms)
            mm_string += '{}_'.format(forceField)
        
        else:
            pass
            
        
        
        filename = str(psystem.e_step_counter)+'-'+psystem.e_tag +'_'+mm_string + qc_string.upper() + stype
        #filename = '{}-{}_{}_{}_{}'str(psystem.e_step_counter)+'-'+system.e_tag +'_'+mm_string+ '_'+qc_string+' '+ stype
        return filename
    
    
    def get_color_palette (self, system_id = None):
        """ Function doc """
        if system_id:
            system = self.psystem[system_id]
            return system.e_color_palette
        else:
            system = self.psystem[self.active_id]
            return  system.e_color_palette 
    
    
    def _build_vobject_from_pdynamo_system(
            self,
            system=None,
            name='a_new_vismol_obj',
            e_id=None,
            is_vobject_active=True,
            autocenter=True,
            # refresh_qc_and_fixed=True,
            # add_vobject_to_vm_session=True,
            frames=None):
        """
        Build a VismolObject from a pDynamo System object.

        This function translates a molecular system represented internally
        in pDynamo into a VismolObject, which can be visualized, manipulated,
        and integrated with the EasyHybrid / Vismol session.

        Parameters
        ----------
        system : pDynamo System
            The pDynamo System instance to be converted into a visualizable
            VismolObject.

        name : str
            The name assigned to the new VismolObject.

        e_id : int, optional
            External identifier for linking this object with higher-level
            EasyHybrid structures.

        is_vobject_active : bool
            Whether the new VismolObject should be set as active.

        autocenter : bool
            If True, the object will be centered when first created.

        frames : np.ndarray, optional
            Pre-existing coordinates (frames) to assign. If None, the
            coordinates are extracted from the pDynamo system.

        Returns
        -------
        VismolObject
            The newly built VismolObject containing atoms, residues, chains,
            bonds, coordinates, and other metadata derived from the
            pDynamo System.
        """

        # -------------------------------------------------------------------------
        # 1. Extract sequence information from the pDynamo system
        # -------------------------------------------------------------------------
        sequence = self._get_sequence_from_pdynamo_system(system)

        # -------------------------------------------------------------------------
        # 2. Prepare atom list and coordinates array
        # -------------------------------------------------------------------------
        atoms = []
        atom_qtty = len(system.atoms.items)

        # Allocate a single frame of coordinates [1 x N_atoms x 3]
        coords = np.empty([1, atom_qtty, 3], dtype=np.float32)

        # Fill coordinates and atom metadata
        j = 0
        for atom in system.atoms.items:
            # Extract XYZ coordinates from the pDynamo system
            xyz = system.coordinates3[atom.index]
            x, y, z = np.float32(xyz[0]), np.float32(xyz[1]), np.float32(xyz[2])
            coords[0, j, :] = x, y, z

            # Check whether system was read from MOL2 file (affects labels)
            is_from_mol2 = getattr(system, 'sequence_from_mol2', False)

            # Build atom metadata dictionary from the pDynamo atom object
            atoms.append(
                self._get_atom_info_from_pdynamo_atom_obj(
                    sequence=sequence,
                    atom=atom,
                    is_from_mol2=is_from_mol2
                )
            )
            j += 1

        # -------------------------------------------------------------------------
        # 3. Create the VismolObject to represent the system visually
        # -------------------------------------------------------------------------
        vm_object = VismolObject(
            self.vm_session,
            len(self.vm_session.vm_objects_dic),  # assign unique index
            name=name,
            color_palette=system.e_color_palette
        )

        # Set model transformation matrix (used for rendering in OpenGL)
        vm_object.set_model_matrix(self.vm_session.vm_glcore.model_mat)

        # -------------------------------------------------------------------------
        # 4. Build hierarchy: Chains → Residues → Atoms
        # -------------------------------------------------------------------------
        initial = time.time()
        atom_id = 0

        for _atom in atoms:
            # Ensure chain exists in the VismolObject
            if _atom["chain"] not in vm_object.chains.keys():
                vm_object.chains[_atom["chain"]] = Chain(vm_object, name=_atom["chain"])
            _chain = vm_object.chains[_atom["chain"]]

            # Ensure residue exists inside the chain
            if _atom["resi"] not in _chain.residues.keys():
                _r = Residue(
                    vm_object,
                    name=_atom["resn"],
                    index=_atom["resi"],
                    chain=_chain
                )
                _chain.residues[_atom["resi"]] = _r

            _residue = _chain.residues[_atom["resi"]]

            # Create the Atom object inside the VismolObject
            atom = Atom(
                vismol_object=vm_object,
                name=_atom["name"],
                index=_atom["index"],
                symbol=_atom["symbol"],
                residue=_residue,
                chain=_chain,
                atom_id=atom_id,
                occupancy=_atom["occupancy"],
                bfactor=_atom["bfactor"],
                charge=_atom["charge"]
            )

            # Assign unique identifiers and update dictionaries
            atom.unique_id = self.vm_session.atom_id_counter
            atom._generate_atom_unique_color_id()

            _residue.atoms[atom_id] = atom
            vm_object.atoms[atom_id] = atom

            atom_id += 1
            self.vm_session.atom_id_counter += 1

        logger.debug(
            "Time used to build the tree: {:>8.5f} secs".format(time.time() - initial)
        )

        # -------------------------------------------------------------------------
        # 5. Assign coordinates and mass center
        # -------------------------------------------------------------------------
        vm_object.frames = coords
        vm_object.mass_center = np.mean(vm_object.frames[0], axis=0)

        # -------------------------------------------------------------------------
        # 6. Extract chemical bonds
        # -------------------------------------------------------------------------
        is_mmState = getattr(system, 'mmState', None)
        index_bonds = None
        bond_orders = None
        if is_mmState:
            # If system contains MM terms, extract harmonic bond terms
            #'''
            for term in system.mmState.mmTerms:
                if term.label == 'Harmonic Bond':
                    print('Bonds defined from pDynamo system topology.')
                    index_bonds = term.Get12Indices()
            #'''
            '''
            for term in system.mmState.mmTerms:
                if term.label == 'Harmonic Bond':
                    index_bonds = []
                    bond_orders = []
                    for k, i_j_bond_indexes in enumerate(system.connectivity.bondIndices):
                        #print(k, i_j_bond_indexes, i_j_bond_indexes[0], i_j_bond_indexes[1])
                        
                        #index_bonds.append(i_j_bond_indexes[0])
                        #index_bonds.append(i_j_bond_indexes[1])
                        pbonds= system.connectivity.bonds
                        #print(pbonds)
                        pbond = pbonds[k]
                        
                        i = pbond.node1.index
                        j = pbond.node2.index
                        index_bonds.append(i)
                        index_bonds.append(j)
                        #print(pbond.type.value[0])
                        bond_orders.append(pbond.type.value[0])
                    print(len(bond_orders), len(index_bonds))
            #'''
                
            if index_bonds:
                vm_object.define_bonds_from_external(index_bonds=index_bonds, bond_orders = bond_orders)
            else:
                vm_object.find_bonded_and_nonbonded_atoms()
                print('Bonds defined from distance.')
        else:
            # No topology: fallback to distance-based bond assignment
            vm_object.find_bonded_and_nonbonded_atoms()
            print('Bonds defined from distance.')

        # -------------------------------------------------------------------------
        # 7. Final metadata and flags
        # -------------------------------------------------------------------------
        vm_object.e_id = system.e_id
        vm_object._generate_color_vectors(self.vm_session.atom_id_counter)
        vm_object.active = True

        # TreeView / ListStore integration (used in GTK interface)
        vm_object.e_treeview_iter = None
        vm_object.liststore_iter = None

        # -------------------------------------------------------------------------
        # 8. Unit Cell (if periodic system)
        # -------------------------------------------------------------------------
        if system.symmetry:
            a = system.symmetryParameters.a
            b = system.symmetryParameters.b
            c = system.symmetryParameters.c
            alpha = system.symmetryParameters.alpha
            beta = system.symmetryParameters.beta
            gamma = system.symmetryParameters.gamma

            vm_object.set_cell(a, b, c, alpha, beta, gamma, color=[0.7, 0.7, 0.2])

        # -------------------------------------------------------------------------
        # Return the fully constructed VismolObject
        # -------------------------------------------------------------------------
        return vm_object


    def delete_system (self, system_e_id = None):
        """ Function doc """
        if system_e_id != None:
            
            if self.psystem[system_e_id].label in self.psystem_name_list:
                index = self.psystem_name_list.index(self.psystem[system_e_id].label)
                self.psystem_name_list.pop(index)
            self.psystem[system_e_id] = None
            self.psystem.pop(system_e_id)
            return True
        else:
            return False


    def define_free_or_fixed_atoms_from_iterable (self, fixedlist = None):
        """ Function doc """
        if fixedlist == []:
            self.psystem[self.active_id].e_fixed_table = []
            self.psystem[self.active_id].freeAtoms = None
            selection_fixed = []
            #self.refresh_qc_and_fixed_representations()

        else:
            selection_fixed                             = Selection.FromIterable (fixedlist)
            ##print( list(selection_fixed))
            self.psystem[self.active_id].e_fixed_table  = list(selection_fixed)
            selection_free                              = selection_fixed.Complement( upperBound = len (self.psystem[self.active_id].atoms ) )
        
            self.psystem[self.active_id].freeAtoms = selection_free
            #self.refresh_qc_and_fixed_representations()
            #self.psystem[self.active_id]['selections']["fixed atoms"] = list(selection_fixed)
            
            
        self.psystem[self.active_id].Summary()
        self.add_a_new_item_to_selection_list (system_id = self.active_id, 
                                                indexes = list(selection_fixed), 
                                                name    = 'fixed atoms')
        
        self.main.refresh_widgets()
        #self.refresh_qc_and_fixed_representations(QC_atoms = False)
        return True


    def get_symmetry_parameters (self, system_e_id = None):
        """ Function doc """
        if system_e_id != None:
            
            if self.psystem[system_e_id].symmetry:

                a = self.psystem[system_e_id].symmetryParameters.a
                b = self.psystem[system_e_id].symmetryParameters.b
                c = self.psystem[system_e_id].symmetryParameters.c
                
                alpha = self.psystem[system_e_id].symmetryParameters.alpha 
                beta  = self.psystem[system_e_id].symmetryParameters.beta  
                gamma = self.psystem[system_e_id].symmetryParameters.gamma
                
                cell  = [a,b,c, alpha, beta, gamma]
                return cell 
                
            else:
                
                return False
            
            
        else:
            return False
    
    
    def set_symmetry_parameters (self, system_e_id = None, 
                                       a     = 0.0,
                                       b     = 0.0,
                                       c     = 0.0,
                                       alpha = 0.0,
                                       beta  = 0.0,
                                       gamma = 0.0
                                       ):
        """ Function doc """
        if system_e_id != None:
            if self.psystem[system_e_id].symmetry:
                self.psystem[system_e_id].symmetryParameters.SetCrystalParameters(a,b,c, alpha,beta,gamma)
            else:
                return False
        else:
            return False


    def make_solvent_box (self, parameters):
        """ Function doc """
        # . Parameters.
        _Density = parameters['_Density']# 1000.0 # . Density of water (kg m^-3).
        _Refine  = parameters['_Refine']# True
        _Steps   = parameters['_Steps']#10000

        # . Box sizes.
        _XBox = parameters['_XBox'] #28.0
        _YBox = parameters['_YBox'] #28.0
        _ZBox = parameters['_ZBox'] #28.0
        
        molecule = parameters['molecule']
        
        # . Define the solvent MM and NB models.
        #mmModel = MMModelOPLS.WithParameterSet ( "bookSmallExamples" )
        #nbModel = NBModelCutOff.WithDefaults ( )
        #
        ## . Define the solvent molecule.
        #molecule = ImportSystem ( os.path.join ( dataPath, "water.mol" ) )
        #molecule.Summary ( )

        # . Create a symmetry parameters instance with the correct dimensions.
        symmetryParameters = SymmetryParameters ( )
        symmetryParameters.SetCrystalParameters ( _XBox, _YBox, _ZBox, 90.0, 90.0, 90.0 )

        # . Create the basic solvent box.
        solvent = BuildSolventBox ( CrystalSystemCubic ( ), symmetryParameters, molecule, _Density )
        solvent.label = "Solvent Box"
        solvent.DefineMMModel ( molecule.mmModel )
        solvent.DefineNBModel ( molecule.nbModel )
        if _Refine:
            ConjugateGradientMinimize_SystemGeometry(solvent                    ,                
                                                     logFrequency           = 10,
                                                     maximumIterations      = parameters['_Steps'],
                                                     rmsGradientTolerance   = 0.01)


        new_system = self.add_new_system_to_psession (system = solvent)
        self.main.main_treeview.add_new_system_to_treeview (new_system)
        ff  =  getattr(new_system.mmModel, 'forceField', "None")

        self.main.bottom_notebook.status_teeview_add_new_item(message = 'New System:  {} ({}) - Force Field:  {}'.format(new_system.label, new_system.e_tag, ff), system = new_system)
        self._add_vismol_object_to_easyhybrid_session (new_system, True) #, name = 'olha o  coco')
        return solvent


    def solvate_system (self, e_id = None, parameters = None):
        """ Function doc """
        # . Retrieve the system.
        _XBox      =  parameters['XBox']
        _YBox      =  parameters['YBox']
        _ZBox      =  parameters['ZBox']
        
        
        #----------------------------------------------------------------------
        _NPositive =  parameters['NPositive']
        if _NPositive > 0:
            cation     =  Unpickle ( parameters['cation'])
        else:
            cation     = None
        #----------------------------------------------------------------------
        
        
        #----------------------------------------------------------------------
        _NNegative =  parameters['NNegative']
        if _NNegative > 0:
            anion      =  Unpickle ( parameters['anion'])
        else:
            anion      =  None
        #----------------------------------------------------------------------
        
        
        solvent    =  parameters['solvent']
        system     =  self.psystem[e_id]
        system.Summary ( )
        
        if parameters['reorient']: masses = Array.FromIterable ( [ atom.mass for atom in system.atoms ] )
        
        # . Reorient the system if necessary (see the results of GetSolvationInformation.py).
        if parameters['reorient']: system.coordinates3.ToPrincipalAxes ( weights = masses )
    
        # . Add the counterions.
        
        if anion is None and cation is None:
            solute = system
        else:
            #print('\n\n\n AddCounterIons \n\n\n\n')
            #print( system, 
            #                          _NNegative, anion, 
            #                          _NPositive, cation, 
            #                          ( _XBox, _YBox, _ZBox )
            #                          )
            
            
            
            solute = AddCounterIons ( system, 
                                      _NNegative, anion, 
                                      _NPositive, cation, 
                                      ( _XBox, _YBox, _ZBox )
                                      )
        solute.Summary ( )
        
        #--------------------------------------------------------------------------------------------
        # . Create the solvated system.
        solution       = SolvateSystemBySuperposition ( solute, solvent, reorientSolute = False )
        solution.label = "Solvated {:s}".format ( system.label )
        self.define_NBModel (_type = 1 , parameters =  None, system = solution)
        #--------------------------------------------------------------------------------------------

        #--------------------------------------------------------------------------------------------
        new_system = self.add_new_system_to_psession (system = solution)
        self.main.main_treeview.add_new_system_to_treeview (new_system)
        ff  =  getattr(new_system.mmModel, 'forceField', "None")

        self.main.bottom_notebook.status_teeview_add_new_item(message = 'New System:  {} ({}) - Force Field:  {}'.format(new_system.label, new_system.e_tag, ff), system = new_system)
        self._add_vismol_object_to_easyhybrid_session (new_system, True) #, name = 'olha o  coco')
        #--------------------------------------------------------------------------------------------

    def clone_system (self, e_id = None, vobject = None, name = 'Unknow', tag = 'UNK', color = [0,1,1]):
        
        if e_id != None:
            system = self.psystem[e_id]
        else:
            system = self.psystem[self.active_id]
        #print(system.label)
        #backup = []
        #backup.append(system.e_treeview_iter)
        #backup.append(system.e_liststore_iter)
        
        #system.e_treeview_iter   = None
        #system.e_liststore_iter  = None

        new_system = copy.deepcopy(system)
        new_system.e_step_counter = 0
        #system.e_treeview_iter   = backup[0]
        #system.e_liststore_iter  = backup[1]
        
        #print('menuitem_clone')

        new_system = self.add_new_system_to_psession (system = new_system,
                                                            name   = name  , 
                                                            tag    = tag   , 
                                                            color  = color )

        self.main.main_treeview.add_new_system_to_treeview (new_system)
        ff  =  getattr(new_system.mmModel, 'forceField', "None")

        self.main.bottom_notebook.status_teeview_add_new_item(message = 'New System:  {} ({}) - Force Field:  {}'.format(new_system.label, new_system.e_tag, ff), system = new_system)
        self._add_vismol_object_to_easyhybrid_session (new_system, True) #, name = 'olha o  coco')
   

    def merge_system (self,  e_id1   = None , 
                             e_id2   = None , 
                             vob_id1 = None ,
                             vob_id2 = None ,
                             name    = 'NoName', 
                             tag     = 'NoTag', 
                             color   = [0,1,1]):
        """ Function doc """
        
        #print (e_id1,e_id2,vob_id1,vob_id2,name,tag,color)
        
        system1 = self.psystem[e_id1]
        vob1 = self.vm_session.vm_objects_dic[vob_id1]
        self.set_psystem_coordinates_from_vobject(vobject   = vob1, 
                                                            system_id = e_id1)
        
        system2 = self.psystem[e_id2]        
        vob2    = self.vm_session.vm_objects_dic[vob_id2]        
        self.set_psystem_coordinates_from_vobject(vobject   = vob2, 
                                                            system_id = e_id2)
        
        system2.Energy()
        system1.Energy()
        
        system = MergeByAtom( system1, system2 )
        self.define_NBModel ( system = system  )
        system.Summary ( )
        system = self.add_new_system_to_psession (system = system,  
                                                        name   = name  , 
                                                        tag    = tag   , 
                                                        color  = color )

        self.main.main_treeview.add_new_system_to_treeview (system)
        ff  =  getattr(system.mmModel, 'forceField', "None")
        self.main.bottom_notebook.status_teeview_add_new_item(message = 'New System:  {} ({}) - Force Field:  {}'.format(system.label, system.e_tag, ff), system =system)
        self._add_vismol_object_to_easyhybrid_session (system, True) #, name = 'olha o  coco')


    def prune_system (self, selection = None, name = 'Pruned System', summary = True, tag = "molsys", color = None):
        """ Function doc """
        p_selection   = Selection.FromIterable ( selection )
        system        = PruneByAtom ( self.psystem[self.active_id], p_selection )
        
        self.define_NBModel(_type = 1, system = system)

        
        
        system.label  = name        
        if summary:
            system.Summary ( )
        
        #print('color', color)
        #print('color', color)
        system = self.add_new_system_to_psession ( 
                                                        system = system,  name = name, tag = tag, color = color )
        
        '''-----------------------------------------------------------------'''
        #system
        self.main.main_treeview.add_new_system_to_treeview (system)
        
        self._add_vismol_object_to_easyhybrid_session (system, True)
        #self.main.refresh_active_system_liststore()

    def add_a_new_item_to_selection_list (self, system_id = None, indexes = [], name  = None):
        """ Function doc """

        
        if name != None:
            pass
        
        else:
            name_default = 'sel_'
            
            loop = True
            while loop:
                name = name_default + str(self.psystem[system_id].e_selections_counter)
                
                if name in self.psystem[system_id].e_selections.keys():
                    
                    name = name_default+ str(self.psystem[system_id].e_selections_counter)
                    
                    self.psystem[system_id].e_selections_counter += 1
                else:
                    loop = False


        self.psystem[system_id].e_selections[name] = indexes    
        #print(self.psystem[system_id].e_selections)    
        if self.main.selection_list_window.visible:
            self.main.selection_list_window.update_window()


    def define_a_new_QCModel (self, system = None, parameters = None, vismol_object = None):
        """ Function doc """
        
        '''Here we have to reload the mmModel original charges.
        This is postulated because multiple associations of QC 
        regions can distort the charge distribution of some residues. 
        (because charge rescale)'''
        
        self.clean_QC_representation_to_vobject()
        
        if vismol_object is None:
            system =  self.psystem[self.active_id]
        
        
        if system:
            pass
        else:
            #system = self.psystem[self.active_id]
            system = self.psystem[vismol_object.e_id]
            
        electronicState = ElectronicState.WithOptions ( charge = parameters['charge'], 
                                                  multiplicity = parameters['multiplicity'], 
                                              isSpinRestricted =  parameters['isSpinRestricted'])
                                              
        system.electronicState = electronicState
        
        
        
        
        
        if parameters['qcengine'] == 'ORCA':
            #print(parameters)
            qcModel = QCModelORCA.WithOptions ( keywords       = [parameters['orca_options']], 
                                                deleteJobFiles = False,
                                                scratch        = parameters['orca_scratch'  ],
                                                 )
            qcModel.randomScratch = parameters['random_scratch']
        
        elif parameters['qcengine'] == 'xTB':
            #print(parameters)
            qcModel = QCModelXTB.WithOptions ( #gfn = 2, keywords = None, vfukui = True,  randomJob = False, parallel = 10)
                                              gfn        = parameters['gfn'       ] ,
                                              parallel   = parameters['parallel'  ] ,
                                              acc        = parameters['acc'       ] ,
                                              fermi_temp = parameters['fermi_temp'] ,
                                              iterations = parameters['iterations'] ,
                                              keywords   = parameters['keywords'  ] ,
                                              scratch    = parameters['scratch'],
                                              #lmo      = parameters['lmo'     ] ,
                                              #json     = parameters['json'    ] ,
                                              #vfukui   = parameters['vfukui'  ] ,
                                              )
        
        elif parameters['qcengine'] == 'MOPAC':
            #print(parameters)
            qcModel = QCModelMOPAC.WithOptions ( 
                                              method   = parameters['hamiltonian'],
                                              keywords =  parameters['keywords'  ],
                                              scratch  = parameters['scratch'    ],
                                              )
        
        elif parameters['qcengine'] == 'DFTB+':
            #print(parameters)
            qcModel = QCModelDFTB.WithOptions ( deleteJobFiles =       parameters["deleteJobFiles"      ],
                                                extendedInput =        parameters["extendedInput"       ],
                                                fermiTemperature =     parameters["fermiTemperature"    ],
                                                gaussianBlurWidth =    parameters["gaussianBlurWidth"   ],
                                                #hamiltonian =          parameters["hamiltonian"         ],
                                                maximumSCCIterations = parameters["maximumSCCIterations"],
                                                randomScratch =        parameters["randomScratch"       ],
                                                sccTolerance =         parameters["sccTolerance"        ],
                                                scratch =              parameters["scratch"             ],
                                                skfPath =              parameters["skfPath"             ],
                                                useSCC =               parameters["useSCC"              ],
                                                 )
                                                 
            '''
            qcModel = QCModelDFTB.WithOptions ( deleteJobFiles = parameters['delete_job_files'],
                                                randomScratch  = parameters['random_scratch']  ,
                                                scratch        = parameters['dftb+_scratch']   ,
                                                skfPath        = parameters['skf_path']        ,
                                                useSCC         = parameters['use_scc']         ,
                                                
                                                ThirdOrderFull       = parameters['ThirdOrderFull'      ],
                                                zeta                 = parameters['zeta'                ],
                                                HubbardDerivs        = parameters['HubbardDerivs'       ],
                                                fermiTemperature     = parameters['fermiTemperature'    ],
                                                gaussianBlurWidth    = parameters['gaussianBlurWidth'   ],
                                                maximumSCCIterations = parameters['maximumSCCIterations'],
                                                sccTolerance         = parameters['sccTolerance'        ],
                                                )
            '''

        else:
            converger = DIISSCFConverger.WithOptions( energyTolerance   = parameters['energyTolerance'] ,
                                                      densityTolerance  = parameters['densityTolerance'] ,
                                                      maximumIterations = parameters['maximumIterations'] )
            qcModel         = QCModelMNDO.WithOptions ( hamiltonian = parameters['method'], converger=converger )
        
        
        
        
        
        #if len(system.e_qc_table) > 0:
        if len(system.e_qc_table) > 0 and len(system.e_qc_table) != len(system.atoms):
            '''This function rescales the remaining charges in the MM part. The 
            sum of the charges in the MM region must be an integer value!'''
            self.check_charge_fragmentation(vismol_object = vismol_object)
            '''----------------------------------------------------------------'''
            try:
                system.DefineQCModel (qcModel, qcSelection = Selection.FromIterable ( system.e_qc_table) )          
            except MMModelError:
                print('\n\n\n MMModelError. Total active MM charge is neither integral nor zero', MMModelError)
                call_message_dialog(text1 = 'MMModelError', text2 = 'Total active MM charge is neither integral nor zero', transient_for =  None)
                return None
            if system.mmModel:
                if parameters['qcengine'] == 'ORCA':
                    system.DefineNBModel (NBModelORCA.WithDefaults ( ))
                
                elif parameters['qcengine'] == 'xTB':
                    system.DefineNBModel (NBModelORCA.WithDefaults ( )) # this is correct, it's the same
                
                elif parameters['qcengine'] == 'MOPAC':
                    system.DefineNBModel (NBModelMOPAC.WithDefaults ( )) # this is temporary
                
                elif parameters['qcengine'] == 'DFTB+':
                    system.DefineNBModel (NBModelDFTB.WithDefaults ( ))
                else:
                    system.DefineNBModel ( NBModelCutOff.WithDefaults ( ) )
            else:
                pass

        else:
            system.DefineQCModel (qcModel)

        system.Summary()
        self.main.refresh_widgets()
        for vismol_object in self.vm_session.vm_objects_dic.values():
            if vismol_object.e_id == system.e_id:
                self._apply_QC_representation_to_vobject (system_id = None, 
                                                          vismol_object = vismol_object)
            else:
                pass
        self.main.refresh_widgets()
        
        if self.main.selection_list_window.visible:
            self.main.selection_list_window.update_window()
    
    
    def check_charge_fragmentation(self, system = None, vismol_object = None,  correction = True):
        """ Function doc """
        #self.psystem[self.active_id]['system_original_charges']
        if system:
            pass
        else:
            system = self.psystem[self.active_id]
        
        mm_residue_table = {}
        qc_residue_table = self.psystem[self.active_id].e_qc_residue_table 
        
        ##print('\n\n\Sum of total charges(MM)', sum(system.mmState.charges))
        '''----------------------------------------------------------------'''
        '''Restoring the original charges before rescheduling a new region.'''
        original_charges = self.psystem[self.active_id].e_charges_backup.copy()
        
        for index, charge in enumerate(original_charges):
            system.mmState.charges[index]   = original_charges[index]
        '''----------------------------------------------------------------'''
        ##print('\n\n\Sum of total charges(MM)', sum(system.mmState.charges))

        qc_charge        = 0.0
        
        if system.mmModel is None:
            return None 
        
        '''Here we are going to arrange the atoms that are not in the QC part, 
        but are in the same residues as those atoms within the QC part.'''  

        #self._check_ref_vobject_in_pdynamo_system()
        for chain in vismol_object.chains:
            for resi, residue in vismol_object.chains[chain].residues.items():
                if resi in qc_residue_table.keys():
                    mm_residue_table[resi] = []
                    
                    for index, atom in residue.atoms.items():
                        index_v = atom.index-1
                        charge  = system.mmState.charges[index_v]
                        resn    = residue.name 
                        atom.charge = system.mmState.charges[index_v]
                        
                        if index_v in qc_residue_table[resi]:
                            qc_charge += atom.charge
                            pass
                        else:
                            mm_residue_table[resi].append(index_v)
        else:
            pass
                
                
                ##print(atom.index, atom.atomicNumber, system.mmState.charges[idx],self.psystem[self.active_id]['vobject'].atoms[idx].resn )
            
        ##print('mm_residue_table',mm_residue_table)
        '''Here we are going to do the charge rescaling of atoms within the MM part 
        but that the residues do not add up to an entire charge.''' 
        
        #print(mm_residue_table)
        #print(qc_residue_table)
        
        for resi in mm_residue_table.keys():
            
            total = 0.0
            for index in mm_residue_table[resi]:
                pcharge = system.mmState.charges[index]
                total += pcharge
            
            rounded  = float(round(total))
            diff     = rounded - total
            size     = len(mm_residue_table[resi])
            
            if size > 0:
                fraction = diff/size
                ##print('residue: ', resi, 'mm_size:', size, 'difference', diff, 'charge fraction = ',fraction)
            
                for index in mm_residue_table[resi]:
                    system.mmState.charges[index] += fraction
                    #total += pcharge
                
                #print('residue: ', resi, '  mm_size:', size, '  difference', diff, '  Correction factor = ',fraction, )
            else:
                pass


    def get_system_name (self, system_id = None):
        """ Function doc """
        if system_id:
            return self.psystem[system_id]['name']
        else:
            return self.psystem[self.active_id]['name']


    def p_selections(self, system_id=None, _centerAtom=None, _radius=None, _method=None):
        """
        Select atoms from a pDynamo system based on a center atom and a radius.

        Parameters
        ----------
        system_id : int or None
            The identifier of the system stored in self.psystem.
            If None, the function will use the currently active system (self.active_id).

        _centerAtom : str or int
            Atom pattern or atom index used as the reference ("center atom")
            for the spatial selection.

        _radius : float
            Cutoff distance (in Angstroms). All atoms within this distance
            from the center atom will be included in the initial selection.

        _method : int
            Defines how the selection will be processed:
            
            - 0 : Select all atoms within the radius, then expand the selection
                  to include the entire molecular component(s) (e.g., whole residues
                  or molecules). Returns the list of atom indices.

            - 1 : Same as method 0, but then invert the selection.
                  Returns all atoms *not* belonging to the selected component(s).

            - 2 : Select only the atoms strictly within the radius.
                  Does not expand to the whole molecular component(s).

        Returns
        -------
        list of int
            The indices of the selected atoms, depending on the chosen method.
        """

        # -------------------------------------------------------------------------
        # 1. Determine which system to use
        # -------------------------------------------------------------------------
        if system_id is not None:
            pass  # If the user provides a system_id, use it directly.
        else:
            system_id = self.active_id  # Otherwise, fall back to the active system.

        # -------------------------------------------------------------------------
        # 2. Get the reference atom from the given pattern/index
        #    This identifies the "center" of the spherical selection.
        # -------------------------------------------------------------------------
        atomref = AtomSelection.FromAtomPattern(self.psystem[system_id], _centerAtom)

        # -------------------------------------------------------------------------
        # 3. Select all atoms within the radius around the reference atom.
        # -------------------------------------------------------------------------
        core = AtomSelection.Within(self.psystem[system_id], atomref, _radius)

        # -------------------------------------------------------------------------
        # 4. Process the selection depending on the method chosen
        # -------------------------------------------------------------------------

        if _method == 0:
            # Expand the selection to include the entire molecular component(s)
            # (e.g., if an atom belongs to a residue, the whole residue is included).
            core = AtomSelection.ByComponent(self.psystem[system_id], core)
            core = list(core)  # Convert to a Python list of indices.

        if _method == 1:
            # Same expansion as in method 0
            core = AtomSelection.ByComponent(self.psystem[system_id], core)
            core = list(core)

            # Now invert the selection: select all atoms NOT in 'core'.
            inverted = []
            for i in range(len(self.psystem[system_id].atoms)):
                if i not in core:
                    inverted.append(i)

            core = inverted  # Replace the selection with the inverted one.

        if _method == 2:
            # Do not expand to components; just use the atoms strictly within the radius.
            core = list(core)

        # -------------------------------------------------------------------------
        # 5. Return the final list of atom indices
        # -------------------------------------------------------------------------
        return core
    
    
    def remove_NBModel (self, system = None ):
        """ Function doc """
        if system:
            system.nbModel = None
            #system.Summary ( )
        else:
            self.psystem[self.active_id].nbModel = None
            #self.psystem[self.active_id].Summary ( )
        return True
    
    
    def remove_PES_data (self, system = False ):
        """ Function doc """
        if system:
            system.e_logfile_data             = {}
        else:
            self.psystem[self.active_id].e_logfile_data             = {}        
        
        msg = 'Log data has been removed.'
        return msg


    def define_NBModel (self, _type = 1 , parameters =  None, system = None):
        """ Function doc """
        
        if _type == 0:
            self.nbModel = NBModelFull.WithDefaults ( )
        
        elif _type == 1:
            self.nbModel = NBModelCutOff.WithDefaults ( )
        
        elif _type == 2:
            self.nbModel = NBModelORCA.WithDefaults ( )
        
        elif _type == 3:
            self.nbModel = NBModelDFTB.WithDefaults ( )
        
        
        try:
            if system:
                system.DefineNBModel ( self.nbModel )
                system.Summary ( )
            else:
                self.psystem[self.active_id].DefineNBModel ( self.nbModel )
                self.psystem[self.active_id].Summary ( )
            return True
        
        except:
            #print('failed to bind nbModel')
            return False
    
    
    def export_pdynamo_system_coordinates (self, folder, filename, system):
        """
            Export to a file the coordinates that are associated with a 
        pDynamo system
        """
        
        Pickle( os.path.join ( folder, filename), system.coordinates3 )
    
    
    def export_system (self,  parameters = {}, coords_from_vobject = True): 
                              
        """  
        Export system model, as pDynamo serization files or Cartesian coordinates.
            0 : 'pkl - pDynamo system'         ,
            1 : 'pkl - pDynamo coordinates'    ,
            2 : 'pdb - Protein Data Bank'      ,
            3 : 'xyz - cartesian coordinates'  ,
            4 : 'mol'                          ,
            5 : 'mol2'                         ,
        
        """
        if 'export_QC_atoms_only' in parameters.keys():
            pass
        else:
            parameters['export_QC_atoms_only'] = False
        
        reverse_idx_2D_xy = {}
        
        if coords_from_vobject:
            vobject  = self.vm_session.vm_objects_dic[parameters['vobject_id']]
            
            
            #. This block refers to a vobject containing a 2D trajectory.
            if getattr ( vobject, 'idx_2D_xy', False):
                max_x = 0
                max_y = 0
                
                #. Checking the dimensions of each coordinate
                for key in vobject.idx_2D_xy.keys():
                    index = vobject.idx_2D_xy[key]
                    if key[0] > max_x:
                        max_x = key[0]
                    if key[1] > max_y:
                        max_y = key[1]
                    reverse_idx_2D_xy[index] = key


                
                #print (max_x, max_y)
                #print (reverse_idx_2D_xy)
                
                #max_x += 1 
                #max_y += 1
                if parameters['last'] == -1:
                    parameters['last'] = max_x
                if parameters['last2'] == -1:
                    parameters['last2'] = max_y

            #. This block refers to a vobject containing a 1D trajectory.
            else:
                if parameters['last'] == -1:
                    parameters['last'] = len(vobject.frames)-1  #.pDynamo is zero base
                else:
                    pass
                    #parameters['last'] =None
        else:
            vobject = None
            
        folder   = parameters['folder'] 
        filename = parameters['filename'] 
        
        

        
        active_id = self.active_id 
        self.active_id = parameters['system_id']
        
        '''
            When 0 or 1, we will export the pDynamo system, and in this 
        case, only the last coordinates are considered. If the coordinate 
        reference is a vobject with two dimensions, this needs to be 
        taken into account when exporting the coordinates.
        '''
        if parameters['format'] == 0 or parameters['format'] == 1:
            
            if coords_from_vobject:
                
                #.Checking is vobject has 2 dimenstions 
                if getattr ( vobject, 'idx_2D_xy', False):
                    frame =  vobject.idx_2D_xy[(parameters['last'],parameters['last2'])]
                else:
                    frame = parameters['last']   
                self.set_psystem_coordinates_from_vobject(vobject       = vobject, 
                                                                    system_id     = parameters['system_id'], 
                                                                    frame         = frame)
            else:
                pass
                
            system = self.psystem[parameters['system_id']]
            
            '''
            The pkl format is not capable of exporting GTK or openGL 
            elements (objects).
            '''
            #backup = []
            #backup.append(system.e_treeview_iter)
            #backup.append(system.e_liststore_iter)
            #backup.append(system.e_logfile_data)
            
            #system.e_treeview_iter   = None
            #system.e_liststore_iter  = None
            
            if parameters['export_QC_atoms_only']:
                system2 = system.PruneToQCRegion()
                if parameters['format'] == 0:
                    ExportSystem ( os.path.join ( folder, filename+'.pkl'), system2 )
                else:
                    YAMLPickle (os.path.join ( folder, filename+'.yaml'),   system2 )
            else:
                if parameters['format'] == 0:
                    ExportSystem ( os.path.join ( folder, filename+'.pkl'), system )
                else:
                    YAMLPickle (os.path.join ( folder, filename+'.yaml'), system )
            
            #system.e_treeview_iter   = backup[0]
            #system.e_liststore_iter  = backup[1]
            #system.e_logfile_data    = backup[2]
        
        
        
        elif parameters['format'] == 2:
            
            #'''   When is Single File     '''
            if parameters['is_single_file']:
                if getattr ( vobject, 'idx_2D_xy', False):
                    frame =  vobject.idx_2D_xy[(parameters['last'] ,parameters['last2'] )]
                else:
                    frame = parameters['last'] 
                
                self.set_psystem_coordinates_from_vobject(vobject   = vobject, 
                                                                    system_id = parameters['system_id'], 
                                                                    frame     = frame)
                
                system   = self.psystem[parameters['system_id']]
                if parameters['export_QC_atoms_only']:
                    system = system.PruneToQCRegion()
                
                
                Pickle( os.path.join ( folder, filename+'.pkl'), 
                        system.coordinates3 )
            
            
            
            #'''   When are Multiple Files   '''
            else:
                #folder   = parameters['folder'] 
                #filename = parameters['filename']
                #os.chdir(folder)
                if os.path.isdir( os.path.join ( folder,filename+".ptGeo")):
                    pass
                else:
                    os.mkdir(os.path.join ( folder,filename+".ptGeo"))
                
                folder = os.path.join ( folder,filename+".ptGeo")
                
                #.Means that this is a 2D trajectory
                if getattr ( vobject, 'idx_2D_xy', False):
                    i = 0
                    for frame_i in range(parameters['first'], parameters['last']+1, parameters['stride']):
                        j = 0
                        for frame_j in range(parameters['first2'], parameters['last2']+1, parameters['stride2']):
                    
                            frame = vobject.idx_2D_xy[(frame_i, frame_j)]
                            self.set_psystem_coordinates_from_vobject(vobject   = vobject, 
                                                                                system_id = parameters['system_id'], 
                                                                                frame     = frame)                
                            system   = self.psystem[parameters['system_id']]
                            
                            if parameters['export_QC_atoms_only']:
                                system = system.PruneToQCRegion()
                            
                            Pickle( os.path.join ( folder, "frame{}_{}.pkl".format(i, j)), 
                                    system.coordinates3 )
                            
                            j+=1
                        
                        i+=1
                
                
                #.Means that is not a 2D trajectory
                else:
                    i = 0
                    for frame in range(parameters['first'], parameters['last']+1, parameters['stride']):
                        
                        self.set_psystem_coordinates_from_vobject(vobject = vobject, 
                                                                                  system_id     = parameters['system_id'], 
                                                                                  frame         = frame)
                        
                        system   = self.psystem[parameters['system_id']]
                        
                        Pickle( os.path.join ( folder, "frame{}.pkl".format(i) ), 
                                system.coordinates3 )
                        
                        i += 1
        
        
        

        else:
            
            #'''   When is Single File     '''
            if parameters['is_single_file']:
                
                self.set_psystem_coordinates_from_vobject(vobject   = vobject, 
                                                                    system_id = parameters['system_id'], 
                                                                    frame     = parameters['last'])
                
                system = self.psystem[parameters['system_id']]
                
                if parameters['format'] == 3:
                    #ExportSystem ( os.path.join ( folder, filename+'.pdb'), system )
                    export_special_PDB(vobject, parameters['last'], os.path.join ( folder, filename+'.pdb'))
                
                if parameters['format'] == 4:
                    ExportSystem ( os.path.join ( folder, filename+'.xyz'), system )
                
                if parameters['format'] == 5:
                    ExportSystem ( os.path.join ( folder, filename+'.mol'), system )
                
                if parameters['format'] == 6:
                    ExportSystem ( os.path.join ( folder, filename+'.mol2'), system )
                
                if parameters['format'] == 7:
                    ExportSystem ( os.path.join ( folder, filename+'.crd'), system )
                    

            
            
            #'''   When are Multiple Files   '''
            else:
                #folder   = parameters['folder'] 
                #filename = parameters['filename']
                #os.chdir(folder)
                if os.path.isdir( os.path.join ( folder,filename+".ptGeo")):
                    pass
                else:
                    os.mkdir(os.path.join ( folder,filename+".ptGeo"))
                
                folder = os.path.join ( folder,filename+".ptGeo")
                
                
                
                #.Means that this is a 2D trajectory
                if getattr ( vobject, 'idx_2D_xy', False):
                    i = 0
                    for frame_i in range(parameters['first'], parameters['last']+1, parameters['stride']):
                        j = 0
                        for frame_j in range(parameters['first2'], parameters['last2']+1, parameters['stride2']):
                    
                            frame = vobject.idx_2D_xy[(frame_i, frame_j)]
                            self.set_psystem_coordinates_from_vobject(vobject   = vobject, 
                                                                                system_id = parameters['system_id'], 
                                                                                frame     = frame)                
                            system   = self.psystem[parameters['system_id']]
                            
                            
                            if parameters['export_QC_atoms_only']:
                               system = system.PruneToQCRegion()
                            
                            #Pickle( os.path.join ( folder, "frame{}_{}.pkl".format(i, j)), 
                            #        system.coordinates3 )
                            if parameters['format'] == 3:
                                #ExportSystem ( os.path.join ( folder, "frame{}_{}.pdb".format(i, j) ), system )
                                export_special_PDB(vobject, frame, os.path.join ( folder, "frame{}_{}.pdb".format(i, j) ))
                            
                            if parameters['format'] == 4:
                                ExportSystem ( os.path.join ( folder, "frame{}_{}.xyz".format(i, j) ), system )
                            
                            if parameters['format'] == 5:
                                ExportSystem ( os.path.join ( folder, "frame{}_{}.mol".format(i, j) ), system )
                            
                            if parameters['format'] == 6:
                                ExportSystem ( os.path.join ( folder, "frame{}_{}.mol2".format(i, j) ), system )

                            j+=1
                        i+=1
                else:
                    i = 0
                    for frame in range(parameters['first'], parameters['last']+1, parameters['stride']):
                        
                        self.set_psystem_coordinates_from_vobject(vobject   = vobject, 
                                                                            system_id = parameters['system_id'], 
                                                                            frame     = frame)

                        system   = self.psystem[parameters['system_id']]
                        
                        if parameters['export_QC_atoms_only']:
                           system = system.PruneToQCRegion()


                        if parameters['format'] == 3:
                            #ExportSystem ( os.path.join ( folder, 'frame{}.pdb'.format( i) ), system )
                            export_special_PDB(vobject, frame, os.path.join ( folder, 'frame{}.pdb'.format( i) ))
                        
                        if parameters['format'] == 4:
                            ExportSystem ( os.path.join ( folder, 'frame{}.xyz'.format( i)), system )
                        
                        if parameters['format'] == 5:
                            ExportSystem ( os.path.join ( folder, 'frame{}.mol'.format( i)), system )
                        
                        if parameters['format'] == 6:
                            ExportSystem ( os.path.join ( folder, 'frame{}.mol2'.format( i)), system )
                        
                        
                        #Pickle( os.path.join ( folder, "frame{}.pkl".format(i) ), 
                        #        system.coordinates3 )
                        
                        i += 1
        


        
        self.active_id  = active_id

    
    def _convert_pdynamo_coords_to_vismol(self, p_coords):
        '''
        This function converts pDynamo coordinates ( which is a 
        list containg all the coords) into a vObject coordinates 
        (which has a different numpy structure). PS: Carlos is 
        the responsible for this mess :D
        '''
        
        frame     = list(p_coords)
        atom_qtty = len(frame)/3
        #print(atom_qtty)
        coords  = np.empty([1, int(atom_qtty), 3], dtype=np.float32)        
        #print (coords)
        atom_id = 0
        
        for i in range(0, len(frame), 3):
            x = np.float32(frame[i  ])
            y = np.float32(frame[i+1])
            z = np.float32(frame[i+2])
            coords[0,atom_id,:] = x, y, z
            atom_id += 1
        #print (coords)
        return coords
        
        
    def get_summary_log (self, system):
        """ Function doc """
        summary_items = system.SummaryItems()
        
        text = ''
        text += '--------------------------------------------------------------------------------\n'
        text += '                           Summary of System {}                          \n'.format(system.label)
        text += '--------------------------------------------------------------------------------\n'
        
        
        text += '------------------------------------- Atoms ------------------------------------\n'
        items = molecule.atoms.SummaryItems()
        new_items = []
        
        for i  in range(len(items)):
            if item[1] == True:
                pass
            else:
                new_items.append(list(item))
        
        #if len()
                
        
        
        text += '{:28}={10:}  {:28}={:10}\n'
 
        
        text += 'Hydrogens                   =        12                                         \n'
        
        
        text += '--------------------------------- Connectivity ---------------------------------\n'
        text += 'Angles                      =        38  Atoms                       =        24\n'
        text += 'Bonds                       =        23  Dihedrals                   =        57\n'
        text += 'Isolates                    =         1  Ring Sets                   =         0\n'
        
        
        text += '----------------------------------- Sequence -----------------------------------\n'
        text += 'Atoms                       =        24  Components                  =         1\n'
        text += 'Entities                    =         1  Linear Polymers             =         0\n'
        text += 'Links                       =         0  Variants                    =         0\n'
        
        
        text += '-------------------------------- AMBER MM Model --------------------------------\n'
        text += '1-4 Interactions            =        57  1-4 Lennard-Jones Form      =     Amber\n'
        text += '1-4 Lennard-Jones Types     =         6  Electrostatic 1-4 Scaling   =     0.833\n'
        text += 'Exclusions                  =       118  Fourier Dihedral Parameters =         3\n'
        text += 'Fourier Dihedral Terms      =        57  Fourier Improper Parameters =         1\n'
        text += 'Fourier Improper Terms      =         1  Harmonic Angle Parameters   =        12\n'
        text += 'Harmonic Angle Terms        =        38  Harmonic Bond Parameters    =         7\n'
        text += 'Harmonic Bond Terms         =        23  LJ Parameters Form          =     Amber\n'
        text += 'LJ Parameters Types         =         6  Number of MM Atom Types     =         7\n'
        text += 'Number of MM Atoms          =        24  Total MM Charge             =     -0.00\n'
        text += '-------------------------------- CutOff NB Model -------------------------------\n'
        text += 'Cell Size                   =     6.750  CutOff Cell Size Factor     =     0.500\n'
        text += 'Damping Cut-Off             =     0.500  Dielectric                  =     1.000\n'
        text += 'Grid Cell/Cell Method       =      True  Inner Cut-Off               =     8.000\n'
        text += 'List CutOff                 =    13.500  Minimum Cell Extent         =         2\n'
        text += 'Minimum Cell Size           =     3.000  Minimum Extent Factor       =     1.500\n'
        text += 'Minimum Points              =       500  Outer Cut-Off               =    12.000\n'
        text += 'Sort Indices                =     False  Use Centering               =      True\n'
        text += '--------------------------------------------------------------------------------\n'


    def check_for_fragmented_charges (self, system_id = None):
        """ Function doc """
        try:
            system = self.psystem[self.active_id]
            itotal  = sum(system.mmState.charges)
            #print('total charge =', total)

            vobj_tmp = self._build_vobject_from_pdynamo_system (system = system)
            n = 0
            for chain in vobj_tmp.chains:
                for resi, residue in vobj_tmp.chains[chain].residues.items():
            
                #for res_index, residue in vobj_tmp.residues.items():
                    n = 0.0
                    res_charge = 0.0
                    for atom_index, atom in residue.atoms.items():
                        res_charge += system.mmState.charges[atom_index]
                        n += 1
                    
                    
                    difference = res_charge - float(round(res_charge))
                    fraction = difference/n

                    res_charge2 = 0.0
                    for atom_index, atom in residue.atoms.items():
                        system.mmState.charges[atom_index] -= fraction
                        res_charge2 += system.mmState.charges[atom_index]
                    
                    n += 1
                
                #print('Initial charge: {}, Differente{} {}'.format(res_charge, difference, res_charge2))

            system.e_charges_backup = list(system.AtomicCharges()).copy()
            ftotal = sum(system.mmState.charges)
            
            
            
            
            return True, 'Inspection of partial charges performed successfully.\n\nNumber of resdues: {} \nInitial sum of charges: {} \nFinal sum of charges: {}'.format(n, itotal, ftotal)
        except:
            return False, 'Inspection of partial charges failed!'
        
        
    def generate_new_empty_vismol_object (self, system_id = None, name = 'new_coordinates' ):
        """  
        This method creates a new vismol object (without coordinates). 
        It is called whenever a new object needs to be given a trajectory.
        
        system_id =  system's index 
        name      =  new vobject's name
        
        returns a vismol_object
        """
        #- - - - - - - - - - - - - - -  Creating a new easyhybrid/vismol object  - - - - - - - - - - - - - - -
        #-----------------------------------------------------------------------------------------------------------------------------
        vismol_object = self._add_vismol_object_to_easyhybrid_session(system = self.psystem[system_id], 
                                                                        name = name)        
         
        self._apply_fixed_representation_to_vobject(system_id = None, vismol_object = vismol_object)
        self._apply_QC_representation_to_vobject   (system_id = None, vismol_object = vismol_object)            
        #-----------------------------------------------------------------------------------------------------------------------------

        #- - - - - - - - - - - - - - -  Cleaning up the residual coordinates  - - - - - - - - - - - - - - -
        #-----------------------------------------------------------------------------------------------------------------------------
        vismol_object.frames  = np.empty([0, len(self.psystem[system_id].atoms), 3], 
                                          dtype=np.float32)
        #-----------------------------------------------------------------------------------------------------------------------------
        return vismol_object
