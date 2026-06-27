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

class EasyHybridImportTrajectory:
    """Handles loading coordinates (PKL-based) into EasyHybrid Vismol objects."""
    
    def __init__(self):
        pass

    # ----------------------------------------------------------------------
    # Helper: load coordinates + optional symmetry from PKL
    # ----------------------------------------------------------------------
    def _load_pkl_coordinates(self, filepath):
        """
        Load pDynamo coordinates stored as PKL files.
        The file may contain:
            - only coordinates3
            - a tuple/list (coordinates3, symmetryParameters)

        Returns
        -------
        (coord, sym)
            coord : Coordinates3
            sym   : symmetry object or None
        """
        with open(filepath, "rb") as f:
            data = pickle.load(f)

        if isinstance(data, (tuple, list)) and len(data) == 2:
            coord, sym = data
        else:
            coord = data
            sym = None

        return coord, sym

    # ----------------------------------------------------------------------
    # Main function to import coordinates
    # ----------------------------------------------------------------------

    def _import_coordinates_from_file (self, parameters):
        """ Function doc """
        if parameters['new_vobj_name']:
            '''it is necessary to create a new object'''
            self.psystem[parameters['system_id']].coordinates3 = ImportCoordinates3( parameters['data_path'] )
            
            vismol_object = self._add_vismol_object_to_easyhybrid_session(system = self.psystem[parameters['system_id']], 
                                                                            name = parameters['new_vobj_name'])        
            
            parameters['vobject_id'] = vismol_object.index
            parameters['vobject']    = vismol_object
            
            #self.main.refresh_active_system_liststore()
            self._apply_fixed_representation_to_vobject(system_id = None, vismol_object = vismol_object)
            self._apply_QC_representation_to_vobject   (system_id = None, vismol_object = vismol_object)

        else:
            '''append the coordinates to an existing object'''
            p_coords = ImportCoordinates3 ( parameters['data_path'] )
            v_coords = self._convert_pDynamo_coords_to_vismol(p_coords)
            #print (parameters['vobject'].frames)
            #coords = np.vstack((coords, f))
            parameters['vobject'].frames = np.vstack((parameters['vobject'].frames, v_coords))
            self._apply_QC_representation_to_vobject(vismol_object = parameters['vobject'])
    
    
    def _import_pkl_coordinates_from_file(self, parameters):
        """
        Import a single frame of coordinates written by pDynamo into EasyHybrid.
        Coordinates come from PKL because ImportCoordinates3 cannot handle
        symmetry operations applied during MD, optimization, etc.
        """

        system_id = parameters['system_id']
        system = self.psystem[system_id]

        # --------------------------------------------------------------
        # Load coordinates + symmetry from PKL
        # --------------------------------------------------------------
        coord, sym = self._load_pkl_coordinates(parameters['data_path'])

        # Apply coordinates to pDynamo system
        system.coordinates3 = coord

        # Apply symmetry if available
        if sym is not None:
            #system.symmetry.a = sym.a
            #system.symmetry.b = sym.b
            #system.symmetry.c = sym.c
            param = {'a'    : sym.a    , 
                     'b'    : sym.b    , 
                     'c'    : sym.c    , 
                     'alpha': sym.alpha, 
                     'beta' : sym.beta ,  
                     'gamma': sym.gamma, 
                    }
            #print(param)
        # ==============================================================
        # CASE 1: Create a NEW Vismol object
        # ==============================================================
        if parameters['new_vobj_name']:

            vismol_object = self._add_vismol_object_to_easyhybrid_session(
                system = system,
                name   = parameters['new_vobj_name']
            )

            # Reset dynamic cell buffer
            vismol_object.cell_coordinates = None

            # If system has symmetry, append the first cell frame
            if system.symmetry:
                #self.append_cell_frame(vismol_object, system)
                self.append_cell_frame( vismol_object, system = None, param = param)
            # Store reference for future frames
            parameters['vobject'] = vismol_object
            parameters['vobject_id'] = vismol_object.index

            # Apply visual representations
            self._apply_fixed_representation_to_vobject(
                system_id=None,
                vismol_object=vismol_object
            )
            self._apply_QC_representation_to_vobject(
                system_id=None,
                vismol_object=vismol_object
            )

            return


        # ==============================================================
        # CASE 2: Append to EXISTING Vismol object
        # ==============================================================
        vismol_object = parameters['vobject']

        # Convert coords and append
        v_coords = self._convert_pdynamo_coords_to_vismol(coord)
        vismol_object.frames = np.vstack((vismol_object.frames, v_coords))

        # Append dynamic cell frame if needed
        if system.symmetry:
            self.append_cell_frame(vismol_object, system = None, param = param)

        # Update QC atom representation only (fixed representation is static)
        self._apply_QC_representation_to_vobject(vismol_object=vismol_object)

    def _import_coordinates_from_pklfolder (self, parameters):
        files = os.listdir( parameters['data_path'])
        pkl_files = []
        for _file in files:
            # Check if the file is a text file
            if _file.endswith('.pkl'):
                pkl_files.append(_file)


        print ('pDynamo pkl folder:' , parameters)
        
        if parameters['new_vobj_name']:
            '''it is necessary to create a new object'''
            #- - - - - - - - - - - - - - -  Creating a new easyhybrid/vismol object  - - - - - - - - - - - - - - -
            #-----------------------------------------------------------------------------------------------------------------------------
            vismol_object = self.generate_new_empty_vismol_object (system_id = parameters['system_id'], name = parameters['new_vobj_name'])
            parameters['vobject_id'] = vismol_object.index
            parameters['vobject']    = vismol_object     

            '''            
            1D trajectories (.ptGeo) must be interpreted by pDynamo's 
            "ImportTrajectory". This is because the generated pkl files 
            may contain information about symmetry and periodicity, and 
            not interpreted directly by the ImportCoordinates3 function.
            '''            
            vismol_object.cell_coordinates = None
            #-----------------------------------------------------------------------------------------------------------------------------
            trajectory = ImportTrajectory ( parameters['data_path'], self.psystem[parameters['system_id']] )
            trajectory.ReadHeader ( )
            # . Loop over the frames in the trajectory.
            
            n = 0 
            while trajectory.RestoreOwnerData ( ):
                #if  n >= parameters['first']:
                p_coords = self.psystem[parameters['system_id']].coordinates3
                v_coords = self._convert_pdynamo_coords_to_vismol(p_coords)
                vismol_object.frames = np.vstack((vismol_object.frames, v_coords))
                system = self.psystem[parameters['system_id']]
                
                # In NPT simulations, the box can vary; therefore, 
                # it is necessary to create a coordinate vector for the box.
                if  system.symmetry:
                    self.append_cell_frame(vismol_object, system)
                
                n+=1
                
                
            trajectory.ReadFooter ( )
            trajectory.Close ( )
            #-----------------------------------------------------------------------------------------------------------------------------
        else:
            trajectory = ImportTrajectory ( parameters['data_path'], self.psystem[parameters['system_id']] )
            trajectory.ReadHeader ( )
            vismol_object = parameters['vobject']
            #vismol_object.cell_coordinates = None
            # . Loop over the frames in the trajectory.
            while trajectory.RestoreOwnerData ( ):
                p_coords = self.psystem[parameters['system_id']].coordinates3
                v_coords = self._convert_pdynamo_coords_to_vismol(p_coords)
                
                vismol_object.frames = np.vstack((vismol_object.frames, v_coords))
                system = self.psystem[parameters['system_id']]
                
                # In NPT simulations, the box can vary; therefore, 
                # it is necessary to create a coordinate vector for the box.
                if  system.symmetry:
                    self.append_cell_frame(vismol_object, system)


            trajectory.ReadFooter ( )
            trajectory.Close ( )

        self._apply_QC_representation_to_vobject(vismol_object = vismol_object)

    def _import_coordinates_from_pklfolder2D (self, parameters):
        """ Function doc """
        files = os.listdir( parameters['data_path'])
        
        pkl_files = []
        for _file in files:
            # Check if the file is a text file
            if _file.endswith('.pkl'):
                pkl_files.append(_file)
        pkl_files.sort()
        print ('pDynamo pkl folder:' , parameters['data_type'])
        print ('Number of pkl files:', len(pkl_files))
        
        self.psystem[parameters['system_id']].coordinates3 = ImportCoordinates3 (os.path.join(parameters['data_path'], 
                                                                                              'frame0_0.pkl' )
                                                                                 )

        vismol_object = self._add_vismol_object_to_easyhybrid_session(system = self.psystem[parameters['system_id']], 
                                                                        name = parameters['new_vobj_name'])        
        
        vismol_object.frames  = np.empty([0, len(self.psystem[parameters['system_id']].atoms), 3], dtype=np.float32)
        
        vismol_object.idx_2D_xy = {}
        vismol_object.idx_2D_f  = {}
        
        #print("\n\n\n")
        #print("parameters['vobject_id'] = vismol_object.e_id")
        #print("parameters['vobject']    = vismol_object")
        parameters['vobject_id'] = vismol_object.index
        parameters['vobject']    = vismol_object
        #print(parameters['vobject_id'])
        #print(parameters['vobject']   )
        
        n = 0
        for _file in pkl_files:
            x_y = _file[5:-4].split('_')
            #print(x_y)
            p_coords = ImportCoordinates3 (os.path.join(parameters['data_path'], _file))
            v_coords = self._convert_pdynamo_coords_to_vismol(p_coords)
            
            vismol_object.idx_2D_xy[(int(x_y[0]), int(x_y[1]))] = n
            vismol_object.idx_2D_f [n] = (int(x_y[0]), int(x_y[1]))
            
            
            vismol_object.frames = np.vstack((vismol_object.frames, v_coords))
            n+=1
        self._apply_QC_representation_to_vobject(vismol_object = vismol_object)

    def _import_normal_modes_data (self, parameters):
        """ Function doc """
        #- - - - - - - - - - - - - - -  Creating a new easyhybrid/vismol object  - - - - - - - - - - - - - - -
        #-----------------------------------------------------------------------------------------------------------------------------
        if parameters['vobject_id'] == -1:
            vismol_object = self._add_vismol_object_to_easyhybrid_session(system = self.psystem[parameters['system_id']], 
                                                                            name = 'normal modes')        
            parameters['vobject_id'] = vismol_object.index
            parameters['vobject']    = vismol_object           
            self._apply_fixed_representation_to_vobject(system_id = None, vismol_object = vismol_object)
            self._apply_QC_representation_to_vobject   (system_id = None, vismol_object = vismol_object)            
        
            #- - - - - - - - - - - - - - -  Cleaning up the residual coordinates  - - - - - - - - - - - - - - -
            #-----------------------------------------------------------------------------------------------------------------------------
            vismol_object.frames  = np.empty([0, len(self.psystem[parameters['system_id']].atoms), 3], 
                                              dtype=np.float32)
            #-----------------------------------------------------------------------------------------------------------------------------
        
        
        else:
            parameters['vobject'] = self.vm_session.vm_objects_dic[parameters['vobject_id']] #vismol_object
            vismol_object         = self.vm_session.vm_objects_dic[parameters['vobject_id']] #vismol_object
        #-----------------------------------------------------------------------------------------------------------------------------





        '''
        1D trajectories (.ptGeo) must be interpreted by pDynamo's 
        "ImportTrajectory". This is because the generated pkl files 
        may contain information about symmetry and periodicity, and 
        not interpreted directly by the ImportCoordinates3 function.
        '''
        modes_dict = {}
        
        x = len(vismol_object.frames)  
        for trajectory in parameters['trajectories']:
            
            log   = open(os.path.join(trajectory,'frequency.log'), 'r')
            line  = log.readline()
            line2 = line.split()
            mode  = line2[0].split('_')
            mode  = int(mode[-1])
            
            frequency = line2[2]
            
            modes_dict[mode] = [frequency, [x, None]]
            #print(trajectory)
            #-----------------------------------------------------------------------------------------------------------------------------
            trajectory = ImportTrajectory (trajectory, self.psystem[parameters['system_id']] )
            trajectory.ReadHeader ( )
            # . Loop over the frames in the trajectory.
            
            while trajectory.RestoreOwnerData ( ):
                p_coords = self.psystem[parameters['system_id']].coordinates3
                v_coords = self._convert_pdynamo_coords_to_vismol(p_coords)
                vismol_object.frames = np.vstack((vismol_object.frames, v_coords))
                x += 1
            
            modes_dict[mode][1][1] = x
            
            trajectory.ReadFooter ( )
            trajectory.Close ( )
            #-----------------------------------------------------------------------------------------------------------------------------
        vismol_object.normal_modes_dict = modes_dict
        #print (modes_dict)
        return modes_dict

    def _import_dcd_file (self, parameters):
        """ Function doc """
        #system self.psystem[parameters['system_id']]
        #trajectory = ImportTrajectory ( trajectoryPath, system )
        print('_import_dcd_file')
        if parameters['new_vobj_name']:
            #- - - - - - - - - - - - - - -  Creating a new easyhybrid/vismol object  - - - - - - - - - - - - - - -
            #-----------------------------------------------------------------------------------------------------------------------------
            vismol_object = self.generate_new_empty_vismol_object (system_id = parameters['system_id'], name = parameters['new_vobj_name'])
            parameters['vobject_id'] = vismol_object.index
            parameters['vobject']    = vismol_object  
            trajectory = ImportTrajectory ( parameters['data_path'], self.psystem[parameters['system_id']] )
            trajectory.ReadHeader ( )
            
            system = self.psystem[parameters['system_id']]
            
            #we have to delete frame zero (which is not part of the trajectory)
            vismol_object.cell_coordinates = None
            
            while trajectory.RestoreOwnerData ( ):
                p_coords = self.psystem[parameters['system_id']].coordinates3
                v_coords = self._convert_pdynamo_coords_to_vismol(p_coords)
                vismol_object.frames = np.vstack((vismol_object.frames, v_coords))
                system = self.psystem[parameters['system_id']]
                
                # In NPT simulations, the box can vary; therefore, 
                # it is necessary to create a coordinate vector for the box.
                if  system.symmetry:
                    # here we are generating the dynamic cell 
                    self.append_cell_frame(vismol_object, system)
                    #a = system.symmetryParameters.a
                    #b = system.symmetryParameters.b
                    #c = system.symmetryParameters.c
                    #alpha = system.symmetryParameters.alpha
                    #beta  = system.symmetryParameters.beta
                    #gamma = system.symmetryParameters.gamma
                    #
                    ## let us create a new frame of the box (unit cell)
                    #vertices = vismol_object._calculate_unit_cell_vertices(a, b, c, alpha, beta, gamma)
                    #cell_frame = np.zeros((1, 8, 3), dtype=np.float32)
                    ##
                    #
                    #for i, vertex in enumerate(vertices, 0):
                    #    #print("Vertex {}: {:7.3f} {:7.3f} {:7.3f}".format(i, vertex[0] ,vertex[1] ,vertex[2]))
                    #    cell_frame[0,i,:] = vertex[0] ,vertex[1] ,vertex[2] 
                    ##print(cell_frame)
                    #
                    #if vismol_object.cell_coordinates is not None:
                    #    vismol_object.cell_coordinates = np.vstack((vismol_object.cell_coordinates, 
                    #                                                cell_frame))
                    #else:
                    #    vismol_object.cell_coordinates = cell_frame
                                                                
            
            trajectory.ReadFooter ( )
            trajectory.Close ( )
        
        else:
            trajectory = ImportTrajectory ( parameters['data_path'], self.psystem[parameters['system_id']] )
            trajectory.ReadHeader ( )
            vismol_object = parameters['vobject']
            # . Loop over the frames in the trajectory.
            while trajectory.RestoreOwnerData ( ):
                p_coords = self.psystem[parameters['system_id']].coordinates3
                v_coords = self._convert_pdynamo_coords_to_vismol(p_coords)
                #vismol_object.frames = np.vstack((vismol_object.frames, v_coords))
                system = self.psystem[parameters['system_id']]
                vismol_object.frames = np.vstack((vismol_object.frames, v_coords))
                if  system.symmetry:
                    self.append_cell_frame(vismol_object, system)
            
            trajectory.ReadFooter ( )
            trajectory.Close ( )        
        self._apply_QC_representation_to_vobject(vismol_object = vismol_object)    
            
    def import_data (self, parameters):
        
        """ Function doc 
        parameters = {
                      'system_id'    : None,
                                     
                      'data_path'    : None,
                      'data_type'    : 0   ,
                      
                      'new_vobj_name': None, 
                      'vobject_id'   : None,
                      'vobject'      : None,

                      'logfile'      : None,
                      
                      'first'        : None,
                      'last'         : None,
                      'stride'       : None,
                     }
        
        """
        
        if parameters['data_type'] in ['pklfile']:#,'pdbfile', 'xyz', 'mol2', 'crd']:
            self._import_pkl_coordinates_from_file (parameters)
        
        elif parameters['data_type'] == 'pklfolder':
            self._import_coordinates_from_pklfolder (parameters)

        elif parameters['data_type'] == 'pklfolder2D':
            self._import_coordinates_from_pklfolder2D (parameters)

        elif parameters['data_type'] in ['pdbfile', 'xyz', 'mol2', 'crd']:
            self._import_coordinates_from_file (parameters)
        
        elif parameters['data_type'] == 'dcd':
            self._import_dcd_file (parameters)
        
        elif parameters['data_type'] == 'charges':
            charges = self.chrg_file_parser (parameters['data_path'])

            for index, chg in enumerate(charges):
                self.psystem[self.active_id].mmState.charges[index] = float(chg)
            
            #if len(charges) and len(self.psystem[self.system_id]):
            #    pass
            
            
        elif parameters['data_type'] == 'normal_modes':
            data = self._import_normal_modes_data (parameters)
            self.main.main_treeview.refresh_number_of_frames()
            return data
        else:
            pass

        
        
        #print (parameters)
        if parameters['logfile']:
            logfile= LogFileReader(parameters['logfile'])
            data   = logfile.get_data()
        
            #print('vobject', parameters['vobject'], parameters['vobject_id'] )
            
            vobject_id = parameters['vobject_id']
            
            if parameters['isAppend']:
                '''
                When two trajectories are added together. Here EasyHybrid 
                will try to concatenate the log data. As they are not 
                mandatory, and can be the result of different routines, 
                they may not fit perfectly with the final trajectory
                '''
                
                if vobject_id in  self.psystem[parameters['system_id']].e_logfile_data.keys():
                    if len(self.psystem[parameters['system_id']].e_logfile_data[vobject_id]) != 0:
                        '''Adding the lists of RC1 - reaction coordinate and Z - energy'''
                        try:
                            self.psystem[parameters['system_id']].e_logfile_data[vobject_id][0]["RC1"] += data["RC1"]
                            self.psystem[parameters['system_id']].e_logfile_data[vobject_id][0]["Z"] += data["Z"]
                        except:
                            print('Error: could not process the log file.')
                    else:
                        '''In this case, the list already exists, but there is nothing inside.'''
                        self.psystem[parameters['system_id']].e_logfile_data[vobject_id].append(data)
                
                else:
                    '''Here there is not yet a list of data associated with the vobject_id. 
                    In this case, the list needs to be created first.'''
                    self.psystem[parameters['system_id']].e_logfile_data[vobject_id] = []
                    self.psystem[parameters['system_id']].e_logfile_data[vobject_id].append(data)
            
            else:
                if vobject_id in  self.psystem[parameters['system_id']].e_logfile_data.keys():
                    self.psystem[parameters['system_id']].e_logfile_data[vobject_id].append(data)
                else:
                    self.psystem[parameters['system_id']].e_logfile_data[vobject_id] = []
                    self.psystem[parameters['system_id']].e_logfile_data[vobject_id].append(data)
            
            self.main.PES_analysis_window.refresh_vobject_liststore ()
            #print ('\n\n\n\n\n\n\n')
            #print (self.psystem[parameters['system_id']].e_logfile_data)
        else:
            pass
        
        
        self.main.main_treeview.refresh_number_of_frames()
        self.main.main_treeview.refresh_trajectory_scalebar()

    def chrg_file_parser (self,_file = None, _type = None):
        """ Function doc 
        
        for index, chg in enumerate(input_files['charges']):
            system.mmState.charges[index] = chg
        
        """
        #
        
        _type = _file.split('.')
        if len(_type) == 2:
            _type = _type[-1]
        else:
            _type = 'unk'
        
        if _type == 'MOL2' or _type == 'mol2':
            atoms, bonds = read_MOL2(_file)
            charges = atoms['charges']
        
        elif _type == 'txt' or _type == 'unk':
            charges = read_SIMPLE_txt (_file)
            

        elif _type == 'aux':
            HEAT_OF_FORMATION, ATOM_CORE, ATOM_X_OPT, GRADIENTS, CHARGES = read_MOPAC_aux (_file)
            charges = CHARGES
        else:
            pass
        #data = open(_file, 'r')
        #for line in data:
        #    print(line)
        
        print (charges)
        
        return charges
   
    # ----------------------------------------------------------------------
    # Helper function: handles dynamic unit cell extraction and storage
    # ----------------------------------------------------------------------
    def append_cell_frame(self, vismol_obj, system = None, param = None):
        """
        Compute the unit-cell vertices for the current frame and append them to the
        vismol object's dynamic cell-coordinate buffer.

        Parameters
        ----------
        vismol_obj : VismolObject
            Visualization object receiving the new cell frame.

        system : pDynamo System
            Provides symmetry parameters for the current simulation step.
        """
        
        if system is not None:
            # Extract symmetry parameters
            a     = system.symmetryParameters.a
            b     = system.symmetryParameters.b
            c     = system.symmetryParameters.c
            alpha = system.symmetryParameters.alpha
            beta  = system.symmetryParameters.beta
            gamma = system.symmetryParameters.gamma

        else:
            a     = param['a'    ]
            b     = param['b'    ]
            c     = param['c'    ]
            alpha = param['alpha']
            beta  = param['beta' ]
            gamma = param['gamma']
            
        # Compute the 8 vertices of the simulation box
        vertices = vismol_obj._calculate_unit_cell_vertices(
            a, b, c, alpha, beta, gamma
        )

        # Create a new frame for the cell (shape = [1, 8, 3])
        cell_frame = np.zeros((1, 8, 3), dtype=np.float32)

        for i, vertex in enumerate(vertices):
            cell_frame[0, i, :] = vertex[0], vertex[1], vertex[2]

        # Append to existing cell coordinates
        if vismol_obj.cell_coordinates is None:
            vismol_obj.cell_coordinates = cell_frame
        else:
            vismol_obj.cell_coordinates = np.vstack(
                (vismol_obj.cell_coordinates, cell_frame)
            )
