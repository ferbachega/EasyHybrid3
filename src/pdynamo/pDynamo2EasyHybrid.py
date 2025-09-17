#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Lembrar de colocar uma header nesse arquivo

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################
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
from pdynamo.p_methods import MolecularDynamics
from pdynamo.p_methods import ChainOfStatesOptimizePath
from pdynamo.p_methods import NormalModes
from pdynamo.p_methods import EnergyCalculation
from pdynamo.p_methods import EnergyRefinement
from pdynamo.p_methods import UmbrellaSampling

from pdynamo.p_methods import WHAMAnalysis
from pdynamo.LogFileWriter import LogFileReader

from gui.windows.setup.windows_and_dialogs import call_message_dialog


'''
#import our core lib
from SimulationsPreset import Simulation 
#---------------------------------------
from vModel import VismolObject
from vModel.MolecularProperties import ATOM_TYPES_BY_ATOMICNUMBER
from vModel.MolecularProperties import COLOR_PALETTE

'''

'''
from easyhybrid.gui import *
from easyhybrid.gui.PES_analisys_window  import  PotentialEnergyAnalysisWindow 
from easyhybrid.gui.PES_analisys_window  import  parse_2D_scan_logfile 
'''

class LoadAndSaveData:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass


    def save_easyhybrid_session (self, filename = 'session.easy', tmp = False):
        """   
        When the interface makes some modification to the session, 
        a temporary file "filename.easy~" is saved. When the session 
        is saved by the user, the interface checks if there is a 
        temporary file and deletes it.
        """
        easyhybrid_session_data = {}
        backup = {}
        
        '''- - - - - - - - - - pDynamo systems - - - - - - - - - - - '''
        #easyhybrid_session_data['psystem'] = self.psystem
        '''- - - - - - - - - - - - - - - - - - - - - - - - - - - - - '''
        self.main.bottom_notebook.get_active_system_text_from_textbuffer()
        easyhybrid_session_data['systems'] = [ ]
        
        for e_id, system in self.psystem.items():
            
            if system == None:
                pass
            
            else:
                data   = {}
    
                backup[e_id] = []
                backup[e_id].append(system.e_treeview_iter)
                backup[e_id].append(system.e_liststore_iter)
    
                system.e_treeview_iter   = None
                system.e_liststore_iter  = None
                
                data['system'] = system
                data['vobjects'] = []
                for key, vobject in self.vm_session.vm_objects_dic.items():
                    if system.e_id == vobject.e_id:
                        #data['frames'] = vobject.frames
                        #data['color_palette'] = vobject.color_palette
                        vobj_data = {}
                        vobj_data['frames']        = vobject.frames
                        vobj_data['color_palette'] = vobject.color_palette
                        vobj_data['name']          = vobject.name
                        vobj_data['active']        = vobject.active
                        
                        
                        if key in system.e_logfile_data.keys():
                            vobj_data['logfile_data'] = system.e_logfile_data[key]
                        
                        if getattr (vobject, 'idx_2D_xy', False):
                            vobj_data['idx_2D_xy'] = vobject.idx_2D_xy
                        
                        if getattr (vobject, 'is_surface', False):
                            vobj_data['is_surface'] = vobject.is_surface
                        
                        data['vobjects'].append(vobj_data)
                            
                easyhybrid_session_data['systems'].append(data)
        
        #---------------------------------------------------------------
        """   
        When the interface makes some modification to the session, 
        a temporary file "filename.easy~" is saved. When the session 
        is saved by the user, the interface checks if there is a 
        temporary file and deletes it.
        """
        if tmp:
            
            #.Saves a temporary file with the same name as the original session + "~"
            if filename:
                tmpfile = filename+'~'
                with open(tmpfile,'wb') as outfile:
                    pickle.dump(easyhybrid_session_data, outfile)
            
            #.Saves a temporary file when there is no previous session name.
            else:
                tmpfile = os.path.join(self.vm_session.vm_config.easyhybrid_tmp, 
                                       self.random_code+'.easy')
                with open(tmpfile,'wb') as outfile:
                    pickle.dump(easyhybrid_session_data, outfile)
                


        else:
            with open(filename,'wb') as outfile:
                pickle.dump(easyhybrid_session_data, outfile)
                self.changed = False
            if os.path.exists(filename+'~'):
                os.remove(filename+'~')
            else:
                pass
        #---------------------------------------------------------------
            
        for e_id, data in backup.items():
            pass
            #self.psystem[e_id]
            #self.psystem[e_id].e_treeview_iter   = data[0]
            #self.psystem[e_id].e_liststore_iter  = data[1]
        
        self.main.session_filename = filename
        
        if tmp:
            pass
        else:
            self.main.bottom_notebook.status_teeview_add_new_item(message = ':  {}  saved'.format(filename), 
                                                               system =  system )
        
        #'''- - - - - - - - - - - - vismol obejcts - - - - - - - - - - - '''
        #vobjects = {}
        #for key, vobject in self.vm_session.vm_objects_dic.items():
        #    parameters = {
        #                  'index'             : vobject.index            ,
        #                  'name'              : vobject.name             ,
        #                  'active'            : vobject.active           ,
        #                  'frames'            : vobject.frames           ,
        #                  'color_palette'     : vobject.color_palette    ,
        #                  'mass_center'       : vobject.mass_center      ,
        #                  'selected_atom_ids' : vobject.selected_atom_ids,
        #                  'index_bonds'       : vobject.index_bonds      ,
        #                                        
        #                  'colors'            : vobject.colors           ,
        #                  'color_indexes'     : vobject.color_indexes    ,
        #                 }
        #    vobjects[key] = parameters
        #easyhybrid_session_data['vobjects'] = vobjects
        #'''- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '''

    def load_easyhybrid_serialization_file (self, filename = None, tmp = False):
        self.main.restart() 
        self.main.bottom_notebook.status_teeview_add_new_item(message = ':  {}  loaded'.format(filename), 
                                                               system =  None )
        if filename is None:
            return None
        with open(filename, "rb") as f:
            # Load the object from the file
            easyhybrid_session_data = pickle.load(f)
        #print(easyhybrid_session_data)
        
        for data  in easyhybrid_session_data['systems']:
            system = data['system']
            name   = system.label
            tag    = system.e_tag
            #print('\n\n\n\n',system, name, tag, data['system'], data['vobjects'] )
            if len(data['vobjects']) == 0:
                pass
            else:
                self.append_system_to_pdynamo_session (system = system, name  = name, tag = tag)
                self.main.main_treeview.add_new_system_to_treeview (system)
                ff  =  getattr(system.mmModel, 'forceField', "None")
                self.main.bottom_notebook.status_teeview_add_new_item(message = 'New System:  {} ({}) - Force Field:  {}'.format(system.label, system.e_tag, ff), system = system)
                
                for vobj  in data['vobjects']:
                    frames = vobj['frames']
                    name   = vobj['name']
                    
                    
                    if 'is_surface' in vobj.keys():
                        pass
                    else:
                        vobj['is_surface'] = False
                    
                    if vobj['is_surface']:
                        print(vobj['is_surface'])
                        pass
                    else:
                        vm_object = self._build_vobject_from_pDynamo_system ( system = system, name = name ) 
                        vm_object.frames = frames
                        vm_object.active = vobj['active']
                        self.vm_session._add_vismol_object(vm_object, show_molecule = True)
                        
                        self.main.main_treeview.add_vismol_object_to_treeview(vm_object)
                        
                        self.main.add_vobject_to_vobject_liststore_dict(vm_object)
                        
                        self._apply_fixed_representation_to_vobject(vismol_object =vm_object)
                        self._apply_QC_representation_to_vobject(vismol_object =vm_object)
                        
                        self.main.refresh_widgets()
                        
                        if 'logfile_data' in vobj.keys():
                            system.e_logfile_data[vm_object.index] = vobj['logfile_data']
                        if 'idx_2D_xy' in vobj.keys():
                            vm_object.idx_2D_xy  = vobj['idx_2D_xy']
        if tmp:
            filename = filename.replace('~', '')
            self.main.session_filename = filename
        else:
            self.main.session_filename = filename

    def save_special_PDB (vObject):
        """ Function doc """
        

class EasyHybridImportTrajectory:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass


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

    def _import_coordinates_from_pklfolder (self, parameters):
        files = os.listdir( parameters['data_path'])
        
        pkl_files = []
        for _file in files:
            # Check if the file is a text file
            if _file.endswith('.pkl'):
                pkl_files.append(_file)

        #print ('pDynamo pkl folder:' , parameters['data_type'])
        #print ('Number of pkl files:', len(pkl_files))
        print ('pDynamo pkl folder:' , parameters)
        
        if parameters['new_vobj_name']:
            '''it is necessary to create a new object'''
            
            #- - - - - - - - - - - - - - -  Creating a new easyhybrid/vismol object  - - - - - - - - - - - - - - -
            #-----------------------------------------------------------------------------------------------------------------------------
            vismol_object = self.generate_new_empty_vismol_object (system_id = parameters['system_id'], name = parameters['new_vobj_name'])
            parameters['vobject_id'] = vismol_object.index
            parameters['vobject']    = vismol_object     
            
            
            '''            
            vismol_object = self._add_vismol_object_to_easyhybrid_session(system = self.psystem[parameters['system_id']], 
                                                                            name = parameters['new_vobj_name'])        
            parameters['vobject_id'] = vismol_object.index
            parameters['vobject']    = vismol_object           
            self._apply_fixed_representation_to_vobject(system_id = None, vismol_object = vismol_object)
            self._apply_QC_representation_to_vobject   (system_id = None, vismol_object = vismol_object)            
            #-----------------------------------------------------------------------------------------------------------------------------
    
            
            
            #- - - - - - - - - - - - - - -  Cleaning up the residual coordinates  - - - - - - - - - - - - - - -
            #-----------------------------------------------------------------------------------------------------------------------------
            vismol_object.frames  = np.empty([0, len(self.psystem[parameters['system_id']].atoms), 3], 
                                              dtype=np.float32)
            #-----------------------------------------------------------------------------------------------------------------------------

            1D trajectories (.ptGeo) must be interpreted by pDynamo's 
            "ImportTrajectory". This is because the generated pkl files 
            may contain information about symmetry and periodicity, and 
            not interpreted directly by the ImportCoordinates3 function.
            '''
            
            #-----------------------------------------------------------------------------------------------------------------------------
            trajectory = ImportTrajectory ( parameters['data_path'], self.psystem[parameters['system_id']] )
            trajectory.ReadHeader ( )
            # . Loop over the frames in the trajectory.
            while trajectory.RestoreOwnerData ( ):
                p_coords = self.psystem[parameters['system_id']].coordinates3
                v_coords = self._convert_pDynamo_coords_to_vismol(p_coords)
                vismol_object.frames = np.vstack((vismol_object.frames, v_coords))
                system = self.psystem[parameters['system_id']]
            #vismol_object = self.interpolate_frames_of_a_vobject (system, vismol_object)

            trajectory.ReadFooter ( )
            trajectory.Close ( )
            #-----------------------------------------------------------------------------------------------------------------------------
            '''
            for frame in range(1,len(pkl_files)):
                p_coords = ImportCoordinates3 (os.path.join(parameters['data_path'], 'frame{}.pkl'.format(frame)))
                v_coords = self._convert_pDynamo_coords_to_vismol(p_coords)
                vismol_object.frames = np.vstack((vismol_object.frames, v_coords))
            '''    
        else:
            
            trajectory = ImportTrajectory ( parameters['data_path'], self.psystem[parameters['system_id']] )
            trajectory.ReadHeader ( )
            vismol_object = parameters['vobject']
            # . Loop over the frames in the trajectory.
            while trajectory.RestoreOwnerData ( ):
                p_coords = self.psystem[parameters['system_id']].coordinates3
                v_coords = self._convert_pDynamo_coords_to_vismol(p_coords)
                
                #vismol_object.frames = np.vstack((vismol_object.frames, v_coords))
                vismol_object.frames = np.vstack((vismol_object.frames, v_coords))

            trajectory.ReadFooter ( )
            trajectory.Close ( )
            
            '''
            for frame in range(0,len(pkl_files)):
                p_coords = ImportCoordinates3 (os.path.join(parameters['data_path'], 'frame{}.pkl'.format(frame)))
                v_coords = self._convert_pDynamo_coords_to_vismol(p_coords)
                parameters['vobject'].frames = np.vstack((parameters['vobject'].frames, v_coords))
            '''
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
            v_coords = self._convert_pDynamo_coords_to_vismol(p_coords)
            
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
                v_coords = self._convert_pDynamo_coords_to_vismol(p_coords)
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
            while trajectory.RestoreOwnerData ( ):
                p_coords = self.psystem[parameters['system_id']].coordinates3
                v_coords = self._convert_pDynamo_coords_to_vismol(p_coords)
                vismol_object.frames = np.vstack((vismol_object.frames, v_coords))
                system = self.psystem[parameters['system_id']]
            #vismol_object = self.interpolate_frames_of_a_vobject (system, vismol_object)

            trajectory.ReadFooter ( )
            trajectory.Close ( )
        
        else:
            trajectory = ImportTrajectory ( parameters['data_path'], self.psystem[parameters['system_id']] )
            trajectory.ReadHeader ( )
            vismol_object = parameters['vobject']
            # . Loop over the frames in the trajectory.
            while trajectory.RestoreOwnerData ( ):
                p_coords = self.psystem[parameters['system_id']].coordinates3
                v_coords = self._convert_pDynamo_coords_to_vismol(p_coords)
                #vismol_object.frames = np.vstack((vismol_object.frames, v_coords))
                vismol_object.frames = np.vstack((vismol_object.frames, v_coords))

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
        
        if parameters['data_type'] in ['pklfile','pdbfile', 'xyz', 'mol2', 'crd']:
            self._import_coordinates_from_file (parameters)
        
        elif parameters['data_type'] == 'pklfolder':
            self._import_coordinates_from_pklfolder (parameters)

        elif parameters['data_type'] == 'pklfolder2D':
            self._import_coordinates_from_pklfolder2D (parameters)

        elif parameters['data_type'] == 'pdbfile':
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
            
            #print ('\n\n\n\n\n\n\n')
            #print (self.psystem[parameters['system_id']].e_logfile_data)
        else:
            pass
        
        
        self.main.main_treeview.refresh_number_of_frames()

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
   
        

class pSimulations:
    """Class responsible for managing and running molecular simulations 
    using multiprocessing and GUI integration with EasyHybrid.
    """

    def __init__(self):
        """Class initializer."""
        self.process_pool = {}   # Holds active processes and their communication queues
        self.psystem = {}        # Dictionary of systems currently loaded
        self.active_id = None    # ID of the currently active system
        self.changed = False     # Tracks whether the session/project has been modified
        self.main = None         # Reference to the main application object
        self.target_process = None

    # ========================================================================
    # RESTRAINTS
    # ========================================================================

    def _apply_restraints(self, parameters):
        """Apply distance restraints to the system if specified in parameters.
        
        parameters['system'].e_restraints_dict[rest_name] = [
            True,                # Active or not
            rest_name,           # Name of the restraint
            'distance',          # Type of restraint
            [atom1, atom2],      # Atoms involved
            distance,            # Target distance
            force_constant       # Force constant
        ]
        """
        # Reset and define a new restraint model
        parameters['system'].DefineRestraintModel(None)
        restraints = RestraintModel()
        parameters['system'].DefineRestraintModel(restraints)

        # Iterate through all defined restraints
        for name, restraint_list in parameters['system'].e_restraints_dict.items():
            if restraint_list[0]:  # Only apply if marked as active
                if restraint_list[2] == 'distance':
                    rmodel = RestraintEnergyModel.Harmonic(restraint_list[4], restraint_list[5])
                    restraint = RestraintDistance.WithOptions(
                        energyModel=rmodel,
                        point1=restraint_list[3][0],
                        point2=restraint_list[3][1]
                    )
                    restraints[name] = restraint

        return parameters

    # ========================================================================
    # PROCESS QUEUE HANDLING
    # ========================================================================

    def _check_queue(self):
        """Check all process queues for messages.

        This function is periodically called by GLib.timeout_add to handle 
        inter-process communication asynchronously without blocking the GTK loop.
        """
        for e_id, (queue, process, result, path, treeiter) in list(self.process_pool.items()):
            try:
                # Read all messages available in the queue
                while not queue.empty():
                    msg = queue.get_nowait()

                    if isinstance(msg, tuple) and msg[0] == "RESULT":
                        # Worker process returned results
                        self._handle_result(msg[1])

                    elif msg == "DONE":
                        # Worker process finished execution
                        self._handle_done(e_id, process, path, treeiter)

                    elif msg == "Running":
                        # Update process status to "Running..."
                        self.main.process_manager_window.set_status(treeiter, "Running...")

            except Exception as e:
                # Prevent GUI crash due to queue handling issues
                print(f"Error checking the queue of process {e_id}: {e}")

        # Keep the GLib timeout active
        return True

    def _handle_result(self, results):
        """Handle 'RESULT' messages from a worker process.

        Updates system coordinates, visualization objects, and job history.
        """
        system = self.psystem[results['e_id']]

        if results.get('new_vobject'):
            # Update coordinates
            system.coordinates3 = results['coords']

            # Generate unique name for visualization object
            name = f"{results['simulation_type']}_{system.e_step_counter}"

            # Add visual object to session
            vobject = self._add_vismol_object_to_easyhybrid_session(system=system, name=name)
            vobject.results = results

            # Save job results in system history
            system.e_job_history[system.e_step_counter] = results
        else:
            # Save job results in system history (no new visual object)
            system.e_job_history[system.e_step_counter] = results

    def _handle_done(self, e_id, process, path, treeiter):
        """Handle 'DONE' messages from a worker process.

        Updates GUI, process pool status, and increments the step counter.
        """
        # Update GUI process manager window
        self.main.process_manager_window.set_status(treeiter, "Finished")
        self.main.process_manager_window.set_time(treeiter, False, True)
        self.main.process_manager_window.set_step_counter(treeiter, self.psystem[e_id].e_step_counter)

        # Mark process as completed in pool
        self.process_pool[e_id][2] = "Finished"

        # Update bottom notebook GUI
        iter_ = self.main.bottom_notebook.status_liststore.get_iter(path)
        value = self.main.bottom_notebook.status_liststore.get_value(iter_, 1)
        self.main.bottom_notebook.status_liststore.set_value(iter_, 1, value.replace('Running...', 'Finished!'))

        # Release process resources
        process.join()

        # Increment step counter
        self.psystem[e_id].e_step_counter += 1

        # Mark project/session as changed
        self.changed = True

    # ========================================================================
    # SIMULATION PROCESS
    # ========================================================================

    def _target_process(self, parameters):
        """Target function executed inside a separate process.

        Runs the selected simulation type and communicates results back 
        to the parent process via a multiprocessing.Queue.
        """
        queue = parameters['queue']
        queue.put("Running")  # Notify parent process that job started

        # Map simulation type to class
        simulation_classes = {
            'Energy_Single_Point': EnergyCalculation,
            'Energy_Refinement': EnergyRefinement,
            'Geometry_Optimization': GeometryOptimization,
            'Molecular_Dynamics': MolecularDynamics,
            'Relaxed_Surface_Scan': RelaxedSurfaceScan,
            'Umbrella_Sampling': UmbrellaSampling,
            'Nudged_Elastic_Band': ChainOfStatesOptimizePath,
            'Normal_Modes': NormalModes,
        }

        sim_type = parameters['simulation_type']
        cls = simulation_classes.get(sim_type)

        if not cls:
            # Unknown simulation type → exit silently
            return

        # Instantiate simulation class
        self.target_process = cls()

        # Special flags for specific simulation types
        if sim_type == 'Energy_Single_Point':
            parameters['energy'] = True
            parameters['new_vobject'] = False
        elif sim_type == 'Energy_Refinement':
            parameters['new_vobject'] = False

        # Avoid passing queue object to child processes inside parameters
        parameters['queue'] = None

        # Run the simulation
        self.target_process.run(parameters)

        # Collect results
        results = {
            'new_vobject': parameters.get('new_vobject', True),
            'energy': parameters.get('energy', False),
            'coords': parameters['system'].coordinates3,
            'e_id': parameters['system'].e_id,
            'simulation_type': sim_type,
            'logfile': parameters['logfile']
        }

        # Send results and completion notification
        queue.put(("RESULT", results))
        queue.put("DONE")

    def run_simulation(self, parameters):
        """Start a new subprocess for a molecular simulation.

        Ensures that only one process per system (e_id) runs at a time.
        """
        # Assign system and apply restraints
        system = self.psystem[self.active_id]
        parameters['system'] = system
        parameters = self._apply_restraints(parameters)

        # Ensure process manager window is visible
        self.main.process_manager_window.OpenWindow()

        # Configure logfile path
        if 'logfile' not in parameters:
            if 'filename' in parameters:
                parameters['logfile'] = os.path.join(parameters['folder'], parameters['filename']) + '.log'
            elif 'trajectory_name' in parameters and parameters['trajectory_name']:
                parameters['logfile'] = os.path.join(parameters['folder'], parameters['trajectory_name'], 'output.log')
            else:
                parameters['logfile'] = os.path.join(parameters['folder'], 'output.log')

        pprint(parameters)

        e_id = system.e_id
        name = system.label

        # Prevent multiple processes for the same system
        if e_id in self.process_pool and self.process_pool[e_id][1].is_alive():
            self.main.simple_dialog.error(
                msg=f"There is already a process underway for the system:\n {e_id} - {name}"
            )
            return False

        # Create queue for inter-process communication
        queue = multiprocessing.Queue()
        parameters['new_vobject'] = True
        parameters['queue'] = queue

        # Create and start subprocess
        process = multiprocessing.Process(
            target=self._target_process,
            args=(parameters,)
        )
        
        try:
            process.start()
            # Add entry to process manager window
            message = f"{parameters['simulation_type']} {system.e_step_counter} - Running..."
            hamiltonian = getattr(system.qcModel, 'hamiltonian', 'unk')
            status = 'Queued'
        except:
            message = f"{parameters['simulation_type']} {system.e_step_counter} - Aborted!"
            hamiltonian = getattr(system.qcModel, 'hamiltonian', 'unk')
            status = 'Aborted!'
        
        treeiter = self.main.process_manager_window.add_new_process(
            system=system,
            _type=parameters['simulation_type'],
            potential=hamiltonian,
            status = status
            )

        # Add entry to bottom notebook
        path = self.main.bottom_notebook.status_teeview_add_new_item(
            message=message,
            system=system
        )

        # Store process information in process pool
        self.process_pool[e_id] = [queue, process, None, path, treeiter]

        # Schedule periodic queue checks
        GLib.timeout_add(200, self._check_queue)
    
            
        
        
        return False


class pAnalysis:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        self.imgPlot = None
        
    def on_mouse_button_press (self, widget, event):
        print (widget, event )
    
    def on_motion (self, widget, event):
        """ Function doc """
        #print(widget, event)

        x_on_plot, y_on_plot, x, y = self.imgPlot.get_i_and_j_from_click_event (event)
        i = y_on_plot
        j = x_on_plot
        
        j_size = len(self.imgPlot.norm_data)
        i_size = len(self.imgPlot.norm_data[0])
        
        if i < 0 or j < 0:
            pass
        else:
            if i < i_size and j < j_size:
                text = 'i {}    |    j {}    |    rc1 {:4.2f}    |    rc2 {:4.2f}    |    E = {:6.3f}'.format(i, j,  
                                                                                         self.imgPlot.dataRC1[j][i], 
                                                                                         self.imgPlot.dataRC2[j][i],  
                                                                                         self.imgPlot.data[j][i])
                print(text)
                    #self.RC_label.set_text(text)






    def run_analysis (self, parameters):
        """ Function doc """
        if parameters['analysis_type'] == 'wham':
            wham_analysis = WHAMAnalysis()
            TrueFalse, results = wham_analysis.run ( parameters )
            
            
            if TrueFalse:
                self.main.simple_dialog.info(msg = results['msg'] )
                
                if results['type'] == 0:
                    from util.easyplot import ImagePlot, XYPlot
                    import random
                    
                    '''                 Histograms                 '''
                    #self.plot = XYPlot(bg_color = [0,0,0])
                    self.plot = XYPlot( )
                    
                    for i , log in enumerate(results['histograms']):
                        
                        X = [] 
                        Y = [] 
                        
                        r = random.random()
                        g = random.random()
                        b = random.random()
                        rgb = [r,g,b,]
                        data = open(log, 'r')
                        for line in data:
                            line2 = line.split()
                            X.append(float(line2[0]) )
                            Y.append(float(line2[1]))
                        self.plot.add ( X = X, Y = Y,
                                        symbol = None, sym_color = [1,1,1], sym_fill = False, 
                                        line = 'solid', line_color = rgb, energy_label = None)
                    
                    #self.plot.Ymax_list= [100]
                    window =  Gtk.Window()
                    window.set_default_size(800, 300)
                    window.move(900, 300)
                    window.set_title('Histograms')
                    window.add(self.plot)
                    window.show_all()
                    
                    X = []
                    Y = []
                    
                    '''                   PMF plot                  '''
                    #self.plot2 = XYPlot(bg_color = [0,0,0])
                    self.plot2 = XYPlot( )
                    print(results['pmf'])
                    data2 = open(results['pmf'], 'r')
                    for line in data2:
                        line2 = line.split()
                        if line2[1] == 'inf':
                            pass
                        else:
                            X.append(float(line2[0]) )
                            Y.append(float(line2[1]))

                    
                    
                    self.plot2.add ( X = X, Y = Y,
                                    symbol = 'dot', sym_color = [0,0,0], sym_fill = False, 
                                    line = 'solid', line_color = [0,0,0], energy_label = None)
                    
                    window2 =  Gtk.Window()
                    window2.set_default_size(800, 300)
                    window2.move(100, 300)
                    window2.set_title(results['pmf'])
                    window2.add(self.plot2)
                    window2.show_all()
                    
                if results['type'] == 1:
                    from util.easyplot import ImagePlot, XYPlot
                    data = open(results['pmf'], 'r')
                    self.imgPlot = ImagePlot()
                    self.imgPlot.connect("button_press_event", self.on_mouse_button_press)
                    self.imgPlot.connect("motion-notify-event", self.on_motion)
                    X =[]
                    Y =[]
                    Z =[]
                    for line in data:
                        line2 = line.split()

                        X.append(float(line2[0]) )
                        Y.append(float(line2[1]))
                        Z.append(float(line2[2]))


                    sizex = X.count(X[0])
                    sizey = Y.count(Y[0])

                    RC1 = []
                    RC2 = []
                    _Z  = []
                    for i in range(sizex):
                        RC1.append(X[sizey*i:sizey*(i+1)])
                        _Z.append( Z[sizey*i:sizey*(i+1)] )
                            #print(i,j, X[sizey*i:sizey*(i+1)])

                    for j in range(len(RC1)):
                        RC2.append(  Y[sizey*j:sizey*(j+1)]  )

                    #print(len(RC1),len(RC1[1]) )
                    #print(len(RC2),len(RC2[1]) )
                    #
                    #print(RC1[1])
                    #print(RC2[1])
                    #print(_Z )
                    data = {
                            'name': 'output.log', 
                            'type': 'self.imgPlot2D', 
                            'RC1': RC1,
                            'RC2': RC2,
                            'Z'  : _Z,
                           }

                    self.imgPlot.show()
                    self.imgPlot.data    =  data['Z']
                    self.imgPlot.dataRC1 =  data['RC1']
                    self.imgPlot.dataRC2 =  data['RC2']
                    self.imgPlot.set_label_mode(mode = 1)

                    #self.imgPlot.set_threshold_color ( _min = 0, _max = 100, cmap = 'jet')
                    window =  Gtk.Window()
                    window.set_default_size(800, 300)
                    window.move(900, 300)
                    window.set_title(results['pmf'])
                    window.add(self.imgPlot)
                    window.show_all()
                    
                    

            else:
                self.main.simple_dialog.error(msg = results['msg'] )
                #dialog = Gtk.MessageDialog(
                #        parent=self.main.window,
                #        flags=Gtk.DialogFlags.MODAL,
                #        type=Gtk.MessageType.ERROR,
                #        buttons=Gtk.ButtonsType.OK,
                #        message_format=msg
                #    )
                #dialog.run()
                #dialog.destroy()

            
            #print (TrueFalse,msg )
            return TrueFalse


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

class Restraints:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
    
    def update_restaint_representation (self, e_id = None):
        """ Function doc """
        try:
            if e_id:
                pass
            else:
                e_id = self.active_id
            
            indexes = []
            for name, restraint in self.psystem[e_id].e_restraints_dict.items():
                _bool = restraint[0]
                name  = restraint[1]
                _type = restraint[2]
                if _type == 'distance':
                    atons = '{} / {}'.format(restraint[3][0],restraint[3][1]) 
                    dist_or_angle = '{:.4f}'.format(restraint[4])
                    force_const   = str(restraint[5])
                    e_id          =  restraint[6] 
                    if _bool:
                        indexes.append(restraint[3][0])
                        indexes.append(restraint[3][1])


            for vobject in self.main.vm_session.vm_objects_dic.values():
                
                if vobject.e_id == e_id:
                    
                    if indexes == []:
                        
                        vobject.representations["restraints"] = None
                    
                    else:
                        
                        if 'restraints' in vobject.representations.keys():
                            
                            if vobject.representations["restraints"] is not None:
                                #print('["restraints"] is not None',indexes)
                                #vobject.representations["restraints"].define_new_indexes_to_vbo(indexes)
                                vobject.representations["restraints"] = DashedLinesRepresentation(vobject, self.vm_session.vm_glcore,
                                                                                          active=True, indexes = indexes)
                            else:
                                #print('else',indexes)
                                vobject.representations["restraints"] = DashedLinesRepresentation(vobject, self.vm_session.vm_glcore,
                                                                                          active=True, indexes = indexes)
                        else:
                            vobject.representations["restraints"] = DashedLinesRepresentation(vobject, self.vm_session.vm_glcore,
                                                                      active=True, indexes = indexes)    
            self.main.vm_session.vm_glcore.queue_draw()
            #print(indexes)
    
        except:
            print('\n\n Failed when trying to represent harmonic potential constraints. This is just a representation error, the potencies are working normally.\n\n')
        
    def add_new_harmonic_restraint (self, parameters, _type = 'distance'):
        """ Function doc """
        #restraints = RestraintModel()
        #parameters['system'].DefineRestraintModel( restraints )
        
        if _type == 'distance':
            
            atom1 = parameters['atom1']
            atom2 = parameters['atom2']
           
            rest_name = str(parameters['system'].e_restraint_counter)
            parameters['system'].e_restraints_dict[rest_name] = [True, 
                                                                 rest_name ,
                                                                 'distance', 
                                                                 [parameters['atom1'],parameters['atom2']], 
                                                                 parameters['distance'],  
                                                                 parameters['force_constant'],
                                                                 parameters['system'].e_id] 
            
            
            parameters['system'].e_restraint_counter += 1
        
        else:
            pass

    def define_harmonic_restraints (self, parameters):
        """ Function doc """
        parameters['equilibriumValue'] = 0.0
        parameters['forceConstant']    = 0.0
        parameters['reference']        = None # <- Clone ( system.coordinates3 )
        parameters['period']           = None
        parameters['selection']        = None
                
        tetherEnergyModel  = RestraintEnergyModel.Harmonic ( 0.0, forceConstant )
        tethers            = RestraintModel ( )        
        tethers["Tethers"] = RestraintMultipleTether.WithOptions ( energyModel = tetherEnergyModel ,
                                                                   reference   = reference         , 
                                                                   selection   = heavies           )
        system.DefineRestraintModel ( tethers )

    def define_distance_harmonic_restraints (self, parameters):
        """ Function doc """
        restraints = RestraintModel()
        #print('\n\n\n\n',parameters, '\n\n\n\n')
        parameters['system'].DefineRestraintModel( restraints )
        
        distance = parameters['dminimum_RC1'] + ( parameters['dincre_RC1'] * float(i) )
        
        rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['force_constant_1'])
        restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom1, point2= atom2)
        restraints["RC1"] = restraint

    def define_position_harmonic_restraints (self, parameters):
        """ Function doc """
        # . parameters['reference']
        # . parameters['selection']
        # . parameters['system']
        # . parameters['force_constant']

        # . Harmonically restrain heavy atoms.
        tethers       = None
        forceConstant = parameters['force_constant']
        
        reference          = Clone ( parameters['reference'] ) # Change it later
        tetherEnergyModel  = RestraintEnergyModel.Harmonic ( 0.0, forceConstant )
        tethers            = RestraintModel ( )
        tethers["Tethers"] = RestraintMultipleTether.WithOptions ( energyModel = tetherEnergyModel ,
                                                                   reference   = reference         , 
                                                                   selection   = parameters['selection'])
        parameters['system'].DefineRestraintModel ( tethers )



    def clear_restraints (self):
        """ Function doc """
        parameters['system'].DefineRestraintModel(None)


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
        
        
    def load_a_new_pDynamo_system_from_dict (self, input_files = {}, system_type = 0, name = None, tag = None, color = None, working_folder = None):
        """ Function doc """
        #print('\n\n\ init - load_a_new_pDynamo_system_from_dict')
        # This commented section prints information about the existing psystem dictionary

        #for index , psystem in self.psystem.items():
        #    if psystem:
        #        print(' in load_a_new_pDynamo_system_from_dict', index, psystem.e_color_palette['C'])
        #    else:
        #        print(' in load_a_new_pDynamo_system_from_dict', index, None)
        
        #psystem = self.generate_pSystem_dictionary()
        
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
        system = self.append_system_to_pdynamo_session (system = system,  name = name, tag = tag, color = color, changed = True )

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
   
    
    def _add_vismol_object_to_easyhybrid_session (self, system, show_molecule=True, name = 'new_coords'):
        """ Function doc """
        # Create a VisMol object from the given pDynamo system
        vm_object = self._build_vobject_from_pDynamo_system ( system = system, name = name ) 
        
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

    
    def append_system_to_pdynamo_session (self, 
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


    def _get_sequence_from_pDynamo_system (self, system = None):
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
        #v_coords = self._convert_pDynamo_coords_to_vismol(p_coords)
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


    def get_coordinates_from_vobject_to_pDynamo_system (self, vobject = None, system_id =  None, frame = -1):
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
    
    
    def _build_vobject_from_pDynamo_system (self                                          , 
                                            system                    = None              , 
                                            name                      = 'a_new_vismol_obj',
                                            e_id                      = None              ,
                                            is_vobject_active         = True              ,
                                            autocenter                = True              ,
                                            #refresh_qc_and_fixed     = True              ,
                                            #add_vobject_to_vm_session = True             ,
                                            frames                    = None
                                           ):
        
        #for i, atom in enumerate(system.atoms):
        #    print(i, atom.atomicNumber, atom.label)
        
        sequence  = self._get_sequence_from_pDynamo_system(system)
        
        #for i, atom in enumerate(system.atoms):
        #    print(i, atom.atomicNumber, atom.label)
        
        atoms     = []     
        atom_qtty = len(system.atoms.items) 
        coords    = np.empty([1, atom_qtty, 3], dtype=np.float32)
        j = 0
        for atom in system.atoms.items:
            xyz = system.coordinates3[atom.index]
            x = np.float32(xyz[0])
            y = np.float32(xyz[1])
            z = np.float32(xyz[2])
            coords[0,j,:] = x, y, z
            
            is_from_mol2 = getattr(system, 'sequence_from_mol2', False)
            #print(atom, atom.label )
            atoms.append(self._get_atom_info_from_pdynamo_atom_obj(sequence =  sequence, atom = atom, is_from_mol2 = is_from_mol2))
            j += 1
            
        vm_object = VismolObject(self.vm_session, 
                                 len(self.vm_session.vm_objects_dic), 
                                 name = name, 
                                 color_palette = system.e_color_palette)
                                 
        vm_object.set_model_matrix(self.vm_session.vm_glcore.model_mat)
        
        # Maybe we have to change it later
        #unique_id = len(self.vm_session.atom_dic_id)

        initial = time.time()
        atom_id = 0
        for _atom in atoms:
            if _atom["chain"] not in vm_object.chains.keys():
                vm_object.chains[_atom["chain"]] = Chain(vm_object, name=_atom["chain"])
            _chain = vm_object.chains[_atom["chain"]]
            
            if _atom["resi"] not in _chain.residues.keys():
                #print(_atom["resn"], _atom["chain"],_atom["resi"] )
                _r = Residue(vm_object, name=_atom["resn"], index=_atom["resi"], chain=_chain)
                #vm_object.residues[_atom["resi"]] = _r
                _chain.residues[_atom["resi"]] = _r
            
            _residue = _chain.residues[_atom["resi"]]
            #print(_residue.name, _residue.index )
            #print(_atom)
            #atom = Atom()
            atom = Atom(vismol_object = vm_object,#
                        name          =_atom["name"] , 
                        index         =_atom["index"],
                        symbol        =_atom["symbol"],
                        residue       =_residue      , 
                        chain         =_chain, 
                        atom_id       = atom_id,
                        occupancy     =_atom["occupancy"], 
                        bfactor       =_atom["bfactor"],
                        charge        =_atom["charge"])
                        

                        
                        
            atom.unique_id = self.vm_session.atom_id_counter
            atom._generate_atom_unique_color_id()
            _residue.atoms[atom_id] = atom
            vm_object.atoms[atom_id] = atom
            atom_id += 1
            #unique_id += 1
            self.vm_session.atom_id_counter += 1
            
        logger.debug("Time used to build the tree: {:>8.5f} secs".format(time.time() - initial))
        

        vm_object.frames = coords
        vm_object.mass_center = np.mean(vm_object.frames[0], axis=0)
        
        '''
            -----------------------------------------------------------
            Extracting chemical bonds from the pdynamo system topology.
        '''
        is_mmState  = getattr(system, 'mmState', None)
        index_bonds = None
        if is_mmState:
            for term in system.mmState.mmTerms:
                if term.label == 'Harmonic Bond':
                    print('Bonds defined from pDynamo system topology.')
                    index_bonds = term.Get12Indices() # some thing like:[1, 0, 2, 1, 3, 2, 4, 3, 5, 0, 5, 4, 6, 0, 7, 6, 8, 7, 9, 8, 10, 1, 10, 9, 11, 10, 12, 9, 13, 6, 14, 13, 15, 14, 16, 15, 17, 16, 18, 13, 18, 17, 19...]
            
            if index_bonds:
                vm_object.define_bonds_from_external(index_bonds = index_bonds)
            else:
                vm_object.find_bonded_and_nonbonded_atoms()
                print('Bonds defined from distance.')
        else:
            vm_object.find_bonded_and_nonbonded_atoms()
            print('Bonds defined from distance.')

        vm_object.e_id = system.e_id
        vm_object._generate_color_vectors(self.vm_session.atom_id_counter)
        vm_object.active = True
        #vm_object.e_treeview_iter_parent_key = None
        '''
            -----------------------------------------------------------
        '''
        
        
        
        '''-----------------------------------------------------'''
        '''Now each pdynamo system object and vismol object has a 
        "treeview_iter" / "liststore_iter" attribute, which is a reference to access 
        the main treeview elements  and the self.main.vobject_liststore_dict[sys_e_id]'''
        vm_object.e_treeview_iter = None
        vm_object.liststore_iter = None
        
        '''When an object is removed it has to be removed from the treeview and 
        vobject_liststore_dict, in addition to the vm_object_dic in the .vm_session.'''
        
        '''-----------------------------------------------------'''




        '''Unit Cell'''
        #'''
        cell_object = None
        if system.symmetry:

            a = system.symmetryParameters.a
            b = system.symmetryParameters.b
            c = system.symmetryParameters.c

            alpha = system.symmetryParameters.alpha 
            beta  = system.symmetryParameters.beta  
            gamma = system.symmetryParameters.gamma
            
            vm_object.set_cell (a, b, c, alpha, beta, gamma, color = [0.7, 0.7, 0.2] )
       
        #'''
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


        new_system = self.append_system_to_pdynamo_session (system = solvent)
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
        new_system = self.append_system_to_pdynamo_session (system = solution)
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
        #system.e_treeview_iter   = backup[0]
        #system.e_liststore_iter  = backup[1]
        
        #print('menuitem_clone')

        new_system = self.append_system_to_pdynamo_session (system = new_system,
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
        self.get_coordinates_from_vobject_to_pDynamo_system(vobject   = vob1, 
                                                            system_id = e_id1)
        
        system2 = self.psystem[e_id2]        
        vob2    = self.vm_session.vm_objects_dic[vob_id2]        
        self.get_coordinates_from_vobject_to_pDynamo_system(vobject   = vob2, 
                                                            system_id = e_id2)
        
        system2.Energy()
        system1.Energy()
        
        system = MergeByAtom( system1, system2 )
        self.define_NBModel ( system = system  )
        system.Summary ( )
        system = self.append_system_to_pdynamo_session (system = system,  
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
        system = self.append_system_to_pdynamo_session ( 
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
        
        
        
        
        
        if parameters['method'] == 'ab initio - ORCA':
            #print(parameters)
            qcModel = QCModelORCA.WithOptions ( keywords       = [parameters['orca_options']], 
                                                deleteJobFiles = False,
                                                scratch        = parameters['orca_scratch'  ],
                                                 )
            qcModel.randomScratch = parameters['random_scratch']
        
        elif parameters['method'] == 'xTB':
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
            
        elif parameters['method'] == 'DFTB+':
            #print(parameters)
            qcModel = QCModelDFTB.WithOptions ( deleteJobFiles =       parameters["deleteJobFiles"      ],
                                                extendedInput =        parameters["extendedInput"       ],
                                                fermiTemperature =     parameters["fermiTemperature"    ],
                                                gaussianBlurWidth =    parameters["gaussianBlurWidth"   ],
                                                hamiltonian =          parameters["hamiltonian"         ],
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
                if parameters['method'] == 'ab initio - ORCA':
                    system.DefineNBModel (NBModelORCA.WithDefaults ( ))
                
                elif parameters['method'] == 'xTB':
                    system.DefineNBModel (NBModelORCA.WithDefaults ( ))
                
                elif parameters['method'] == 'DFTB+':
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


    def p_selections (self, system_id = None, _centerAtom = None, _radius = None, _method = None):
        """ Function doc """
       
        if system_id != None:
            pass
        else:
            system_id = self.active_id

        
        
        atomref = AtomSelection.FromAtomPattern( self.psystem[system_id], _centerAtom )
        
        core    = AtomSelection.Within(self.psystem[system_id],
                                            atomref,
                                            _radius)

        if _method ==0:
            core    = AtomSelection.ByComponent(self.psystem[system_id],core)
            core    = list(core)
        
        if _method == 1:
            core    = AtomSelection.ByComponent(self.psystem[system_id],core)
            core    = list(core)
            #'''******************** invert ? **********************
            inverted = []
            for i in range(0, len(self.psystem[system_id].atoms)):
                if i in core:
                    pass
                else:
                    inverted.append(i)
            
            core =  inverted

        if _method == 2:
            core    = list(core)
        
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
                self.get_coordinates_from_vobject_to_pDynamo_system(vobject       = vobject, 
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
                
                self.get_coordinates_from_vobject_to_pDynamo_system(vobject   = vobject, 
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
                            self.get_coordinates_from_vobject_to_pDynamo_system(vobject   = vobject, 
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
                        
                        self.get_coordinates_from_vobject_to_pDynamo_system(vobject = vobject, 
                                                                                  system_id     = parameters['system_id'], 
                                                                                  frame         = frame)
                        
                        system   = self.psystem[parameters['system_id']]
                        
                        Pickle( os.path.join ( folder, "frame{}.pkl".format(i) ), 
                                system.coordinates3 )
                        
                        i += 1
        
        
        

        else:
            
            #'''   When is Single File     '''
            if parameters['is_single_file']:
                
                self.get_coordinates_from_vobject_to_pDynamo_system(vobject   = vobject, 
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
                            self.get_coordinates_from_vobject_to_pDynamo_system(vobject   = vobject, 
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
                        
                        self.get_coordinates_from_vobject_to_pDynamo_system(vobject   = vobject, 
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

    
    def _convert_pDynamo_coords_to_vismol(self, p_coords):
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

            vobj_tmp = self._build_vobject_from_pDynamo_system (system = system)
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


class Atom:
    """ Class doc """
    
    def __init__(self, vismol_object, name="Xx", index=None, residue=None,
                 chain=None, pos=None, symbol=None, atom_id=None, color=None,
                 vdw_rad=None, cov_rad=None, ball_rad=None,
                 occupancy=0.0, bfactor=0.0, charge=0.0, bonds_indexes=None):
        """ Class initializer """
        self.vm_object   = vismol_object
        self.vm_session  = vismol_object.vm_session
        
        
        self.name     = name
        self.index    = index   # - Remember that the "index" attribute refers to the numbering of atoms (it is not a zero base, it starts at 1 for the first atom)
        self.residue  = residue
        self.chain    = chain
        self.molecule = None

        self.pos = pos     # - coordinates of the first frame
        self.unique_id = None
        
        # self.symbol = self._get_symbol(self.name) if (symbol is None) else symbol
        if symbol is None:
            self.symbol = self._get_symbol()
        else:
            self.symbol = symbol
        self.atom_id = atom_id
        
        if color is None:
            self.color = self._init_color()
        else:
            self.color = color
        
        if vdw_rad is None:
            self.vdw_rad = self._init_vdw_rad()
        else:
            self.vdw_rad = vdw_rad
        
        if cov_rad is None:
            self.cov_rad = self._init_cov_rad()
        else:
            self.cov_rad = cov_rad
        
        if ball_rad is None:
            self.ball_rad = self._init_ball_rad()
        else:
            self.ball_rad = ball_rad
        
        self.color_id = None
        self.occupancy = occupancy
        self.bfactor = bfactor
        self.charge = charge
        if bonds_indexes is None:
            self.bonds_indexes = []
        else:
            self.bonds_indexes = bonds_indexes
        
        self.selected       = False
        self.lines          = True
        self.dots           = False
        self.nonbonded      = False
        self.impostor       = False
        self.ribbons        = False
        self.ribbon_sphere  = False
        self.dynamic        = False
        self.ball_and_stick = False
        self.sticks         = False
        self.stick_spheres  = False
        self.spheres        = False
        self.vdw_spheres    = False
        self.dash           = False
        self.surface        = False
        self.bonds          = []
        self.isfree         = True
        self.labels         = False
        self.label_text     = None
        
    def _get_symbol(self):
        """ Function doc """
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        
        name = self.name.strip()
        if name == "":
            return ""
        
        _n = name
        for char in name:
            if char.isnumeric():
                _n = _n.replace(char, "")
        name = _n
        
        if len(name) >= 3:
            name = name[:2]
        
        # This can't happen since before all numbers have been converted to strings
        # if len(name) == 2:
        #     if name[1].isnumeric():
        #         symbol = name[0]
        
        # The capitalization of hte name can solve all the next ifs
        # name = name.lower().capitalize()
        if name in ATOM_TYPES.keys():
            return name
        else:
            if name[0] == "H":
                if name[1] == "g":
                    symbol =  "Hg"
                elif name[1] =="e":
                    symbol = "He"
                else:
                    symbol =  "H"
            
            elif name[0] == "C":
                if name[1] == "a":
                    symbol = "Ca"
                elif name[1] =="l":
                    symbol = "Cl"
                elif name[1] =="L":
                    symbol = "Cl"
                elif name[1] =="d":
                    symbol = "Cd"
                elif name[1] =="u":
                    symbol = "Cu"
                elif name[1] =="U":
                    symbol = "Cu"
                else:
                    symbol = "C"
            
            elif name[0] == "N":
                if name[1] == "i" or name[1] == "I":
                    symbol = "Ni"
                elif name[1] == "a":
                    symbol = "Na"
                elif name[1] == "e":
                    symbol = "Ne"
                elif name[1] == "b":
                    symbol = "Nb"
                else:
                    symbol = "N"
            
            elif name[0] == "O":
                if name[1] == "s":
                    symbol = "Os"
                else:
                    symbol = "O"
            
            elif name[0] == "S":
                if name[1] == "I":
                    symbol = "Si"
                elif name[1] == "e":
                    symbol = "Se"
                
                elif name[1] == "O":
                    try:
                        if name[2] == "D":
                            symbol = "Na"
                        else:
                            pass
                    except:
                        symbol = "S"
                else:
                    symbol = "S"
            
            elif name[0] == "P":
                if name[1] == "d":
                    symbol = "Pd"
                elif  name[1] == "b":
                    symbol = "Pb"
                elif  name[1] == "o":
                    symbol = "Po"
                else:
                    symbol = "P" 
            
            elif name[0] == "Z":
                if name[1] == "r":
                    symbol = "Zr"
                elif  name[1] == "N":
                    symbol = "Zn"
                else:
                    symbol = "Zn" 
            
            elif name[0] == "F":
                if name[1] == "E":
                    symbol = "Fe"
                elif  name[1] == "e":
                    symbol = "Fe"
                else:
                    symbol = "F" 
            
            elif name[0] == "M":
                if name[1] == "n":
                    symbol = "Mn"
                elif name[1] == "N":
                    symbol = "Mn"
                elif name[1] == "o":
                    symbol = "Mo"
                elif name[1] == "G":
                    symbol = "Mg"
                else:
                    symbol = "X"
            else:
                symbol = "X"
        return symbol
    
    def _init_color(self):
        """ Return the color of an atom in RGB. Note that the returned
            value is in scale of 0 to 1, but you can change this in the
            index. If the atomname does not match any of the names
            given, it returns the default dummy value of atom X.
        """
        try:
            color = self.vm_object.color_palette[self.name]
        except KeyError:
            # #print(self.symbol, "Atom")
            color = self.vm_object.color_palette[self.symbol]
        return np.array(color, dtype=np.float32)
    
    def _generate_atom_unique_color_id(self):
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        """ Function doc """
        r = (self.unique_id & 0x000000FF) >> 0
        g = (self.unique_id & 0x0000FF00) >> 8
        b = (self.unique_id & 0x00FF0000) >> 16
        # pickedID = r + g * 256 + b * 256*256
        self.color_id = np.array([r/255.0, g/255.0, b/255.0], dtype=np.float32)
    
    def _init_vdw_rad(self):
        """ Function doc """
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        try:
            vdw = ATOM_TYPES[self.name][6]
        except KeyError:
            vdw = ATOM_TYPES[self.symbol][6]
        return vdw
    
    def _init_cov_rad(self):
        """ Function doc """
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        try:
            cov = ATOM_TYPES[self.name][5]
        except KeyError:
            cov = ATOM_TYPES[self.symbol][5]
        return cov
    
    def _init_ball_rad(self):
        """ Function doc """
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        try:
            ball = ATOM_TYPES[self.name][6]
        except KeyError:
            ball = ATOM_TYPES[self.symbol][6]
        return ball
    
    def coords(self, frame=None):
        """ 
        frame = int
        
        Returns the coordinates of an atom according to the specified frame. 
        If no frame is specified, the frame set by easyhybrid (probably by the 
        scale bar of the trajectory manipulation window) is used. If the object 
        (vobject) has a smaller number of frames than the one set by the 
        interface, the last frame of the object is used.
        
        return  xyz 
        """
        if frame is None:
            frame  = self.vm_object.vm_session.frame
            #print (frame, len(self.vm_object.frames))
            if len(self.vm_object.frames)-1 <= frame:
                frame = len(self.vm_object.frames)-1
            else:
                pass
        return self.vm_object.frames[frame, self.atom_id]
    
    def get_grid_position(self, gridsize=3, frame=None):
        """ Function doc """
        coords = self.coords(frame)
        gridpos = (int(coords[0]/gridsize), int(coords[1]/gridsize), int(coords[2]/gridsize))
        return gridpos
    
    def get_cov_rad(self):
        """ Function doc """
        return self.cov_rad
    
    def init_color_rgb(self, name):
        """ Return the color of an atom in RGB. Note that the returned
            value is in scale of 0 to 1, but you can change this in the
            index. If the atomname does not match any of the names
            given, it returns the default dummy value of atom X.
        """
        
        try:
            color = color =self.vm_object.color_palette[name]
            #color = ATOM_TYPES[name][1]
        except KeyError:
            if name[0] == "H":# or name in self.hydrogen:
                #color = ATOM_TYPES["H"][1]
                color = self.vm_object.color_palette["H"]
            
            elif name[0] == "C":
                #color = ATOM_TYPES["C"][1]
                color = self.vm_object.color_palette["C"]
            
            elif name[0] == "O":
                #color = ATOM_TYPES["O"][1]
                color = self.vm_object.color_palette["O"]
            
            elif name[0] == "N":
                #color = ATOM_TYPES["N"][1]
                color = self.vm_object.color_palette["N"]
                
            elif name[0] == "S":
                #color = ATOM_TYPES["S"][1]
                color = self.vm_object.color_palette["S"]
            else:
                #color = ATOM_TYPES["X"][1]
                color = self.vm_object.color_palette["X"]
                
        color = [int(color[0]*250), int(color[1]*250), int(color[2]*250)]
        return color
    
    def init_radius(self, name):
        """
        """
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        try:
            rad = ATOM_TYPES[name][6]/5.0
        except KeyError:
            if name[0] == "H" or name in self.hydrogen:
                rad = ATOM_TYPES["H"][6]/5.0
            elif name[0] == "C":
                rad = ATOM_TYPES["C"][6]/5.0
            elif name[0] == "O":
                rad = ATOM_TYPES["O"][6]/5.0
            elif name[0] == "N":
                rad = ATOM_TYPES["N"][6]/5.0
            elif name[0] == "S":
                rad = ATOM_TYPES["S"][6]/5.0
            else:
                rad = 0.30
        return rad


def generate_random_code(length):
    # Define the character set: uppercase letters, lowercase letters, and digits
    characters = string.ascii_letters + string.digits
    # Generate the random code by randomly selecting characters from the set
    random_code = ''.join(random.choice(characters) for _ in range(length))
    return random_code









def export_special_PDB (vobject = None, frame = -1, output = 'temp.pdb'):
    """ Function doc """
    
    
    data = open(output, 'w')
    change_chain = False
    text = ''
    for index in vobject.atoms.keys():
        atom    = vobject.atoms[index]
        name    = atom.name
        
        resi    = atom.residue
        chain   = atom.chain
        chain_name = chain.name
        
        if chain_name =='':
            if change_chain:
                chain_name = 'B'
            else:
                chain_name = 'A'
        #resn    = vobject.chain[chain].
        
        #chain   = atom.chain
        
        symbol  = atom.symbol
        coords  = atom.coords(frame)
        
        resi_index =resi.index
        
        if resi_index > 9999:
            change_chain = True
            resi_index = resi_index-9999
        
        #print (atom, name, resi, chain)
        #           ATOM    120  CA  VAL A  17      36.365  31.982  42.405  1.00 53.96           C  
        #           ATOM    514 H514 SER  34        17.413  42.503  16.453  1.00  1.00      H   
        #           HETATM63640  H1  WAT  1916       7.005  19.149   8.699  0.00  0.00          H   
        #           ATOM    116 S116 CYS  9         23.989  35.368   6.299  1.00  1.00      S   
        #           HETATM63640  H1         WAT  1916       7.005  19.149   8.699  0.00  0.00          H   
        ATOMLINE = "{:<6s}{:5d} {:<4s} {:3s} {:1s}{:>4s}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:<4s}\n".format(
        #HETATM
        'ATOM  ',
        index+1,
        str(name),
        str(resi.name),
        str(chain_name),
        str(resi_index),
        coords[0],
        coords[1],
        coords[2],
        1.0,
        1.0,
        symbol
        )
        #print(ATOMLINE)
        text += ATOMLINE
    data.write(text)



#class pSimulations:
#    """ Class doc """
#    
#    def __init__ (self):
#        """ Class initialiser """
#
#    def _apply_restraints (self, parameters):
#        """ Function doc 
#        
#        parameters['system'].e_restraints_dict[rest_name] = [True, 
#                                                     rest_name ,
#                                                     'distance', 
#                                                     [parameters['atom1'],parameters['atom2']], 
#                                                     parameters['distance'],  
#                                                     parameters['force_constant']] 
#
#        """
#        #--------------------------------------------------------------------------
#        parameters['system'].DefineRestraintModel(None)
#        restraints = RestraintModel()
#        parameters['system'].DefineRestraintModel( restraints )
#        for name, restraint_list in parameters['system'].e_restraints_dict.items():
#
#            if restraint_list[0]:
#                
#                if restraint_list[2] == 'distance':
#                    rmodel    = RestraintEnergyModel.Harmonic(restraint_list[4],restraint_list[5])
#                    restraint = RestraintDistance.WithOptions(energyModel = rmodel, 
#                                                                   point1 = restraint_list[3][0] , 
#                                                                   point2 = restraint_list[3][1] )
#                    restraints[name] = restraint
#        #--------------------------------------------------------------------------
#        return parameters
#    
#    
#    def _check_queue(self):
#        """
#        Check all process queues for messages.
#        This function is periodically called by GLib.timeout_add 
#        to handle inter-process communication asynchronously 
#        without blocking the GTK main loop.
#        """
#        for e_id, (queue, process, result, path, treeiter ) in list(self.process_pool.items()):
#            try:
#                # Read all messages currently available in the queue
#                while not queue.empty():
#                    msg = queue.get_nowait()
#
#                    # Message contains results data
#                    if isinstance(msg, tuple) and msg[0] == "RESULT":
#                        self._handle_result(msg[1])
#
#                    # Process has finished execution
#                    elif msg == "DONE":
#                        self._handle_done(e_id, process, path, treeiter)
#                    
#                    elif msg == "Running":
#                        self.main.process_manager_window.set_status(treeiter, "Running...")
#
#            except Exception as e:
#                # Any exception during queue handling is caught here.
#                # This prevents the GUI loop from breaking.
#                print(f"Error checking the queue of process {e_id}: {e}")
#
#        # Keep the GLib timeout active
#        return True
#
#
#    def _handle_result(self, results):
#        """
#        Process 'RESULT' messages received from a worker process.
#        Updates molecular coordinates and optionally adds a new
#        visual object to the EasyHybrid session.
#        """
#        
#        if results.get('new_vobject'):
#            system = self.psystem[results['e_id']]
#            system.coordinates3 = results['coords']
#
#            # Generate a unique name for the new visual object
#            name = f"{results['simulation_type']}_{system.e_step_counter}"
#
#            # Add the new object into the visualization environment
#            vobject = self._add_vismol_object_to_easyhybrid_session(system=system, name=name)
#            vobject.results = results
#            
#            system.e_job_history[system.e_step_counter] = results
#
#        else:
#            system = self.psystem[results['e_id']]
#            system.e_job_history[system.e_step_counter] = results
#            
#            
#            #log  = LogFile(system)
#            #    path = log.path
#            #    with open(path, "r") as f:
#            #        header = f.read()
#            #else:
#            #    header = '' 
#            #
#            #text = header+text
#            #textwindow = TextWindow(text)
#
#    def _handle_done(self, e_id, process, path, treeiter):
#        """
#        Process 'DONE' messages from a worker process.
#        Updates the GUI to reflect the finished process and 
#        increments the step counter of the corresponding system.
#        """
#        self.main.process_manager_window.set_status (treeiter,"Finished")
#        self.main.process_manager_window.set_time(treeiter, False, True)
#        self.main.process_manager_window.set_step_counter(treeiter, self.psystem[e_id].e_step_counter)
#        
#        # Mark process as completed in the pool
#        self.process_pool[e_id][2] = "Finished"
#        
#        # Update bottom_notebook GUI (replace "Running..." with "Done!")
#        iter_ = self.main.bottom_notebook.status_liststore.get_iter(path)
#        value = self.main.bottom_notebook.status_liststore.get_value(iter_, 1)
#        self.main.bottom_notebook.status_liststore.set_value(iter_, 1, value.replace('Running...', 'Finished!'))
#        # Ensure process resources are released
#        process.join()
#
#        # Increment step counter for the system
#        self.psystem[e_id].e_step_counter += 1
#
#        # Mark the project/session as changed
#        self.changed = True
#
#
#    def _target_process(self, parameters):
#        """
#        Target function executed inside a separate process.
#        Runs the appropriate type of simulation and communicates
#        results back to the parent process via a multiprocessing.Queue.
#        """
#        #queue = self.process_pool[parameters['system'].e_id][0]   
#        queue = parameters['queue']
#        queue.put("Running")  # Notify parent that the job has started
#
#        # Mapping between simulation types and corresponding classes
#        simulation_classes = {
#            'Energy_Single_Point': EnergyCalculation,
#            'Energy_Refinement': EnergyRefinement,
#            'Geometry_Optimization': GeometryOptimization,
#            'Molecular_Dynamics': MolecularDynamics,
#            'Relaxed_Surface_Scan': RelaxedSurfaceScan,
#            'Umbrella_Sampling': UmbrellaSampling,
#            'Nudged_Elastic_Band': ChainOfStatesOptimizePath,
#            'Normal_Modes': NormalModes,
#        }
#
#        sim_type = parameters['simulation_type']
#        cls = simulation_classes.get(sim_type)
#
#        if cls:
#            # Instantiate the simulation class
#            self.target_process = cls()
#
#            # Special flags for specific simulations
#            if sim_type == 'Energy_Single_Point':
#                parameters['energy'] = True
#                parameters['new_vobject'] = False
#            elif sim_type == 'Energy_Refinement':
#                parameters['new_vobject'] = False
#        else:
#            # Unknown simulation type → exit silently
#            return
#        parameters['queue'] =  None
#        # Run the actual simulation
#        self.target_process.run(parameters)
#
#        # Collect results
#        results = {
#            'new_vobject'    : parameters.get('new_vobject', True),
#            'energy'         : parameters.get('energy', False),
#            'coords'         : parameters['system'].coordinates3,
#            'e_id'           : parameters['system'].e_id,
#            'simulation_type': sim_type,
#            'logfile'        : parameters['logfile']
#        }
#
#        # Send results and completion message to parent process
#        queue.put(("RESULT", results))
#        queue.put("DONE")
#
#
#    def run_simulation(self, parameters):
#        """
#        Starts a new subprocess for running a molecular simulation.
#        Ensures that only one process per system (e_id) runs at a time.
#        """
#        # Assign system and apply restraints before running
#        system = self.psystem[self.active_id]
#        parameters['system'] = system
#        parameters = self._apply_restraints(parameters)
#        self.main.process_manager_window.OpenWindow()
#        
#        
#        
#        
#        if 'logfile' in parameters.keys():
#            pass
#        else:
#            if 'filename' in parameters.keys():
#            
#                full_path_trajectory = os.path.join(parameters['folder'], 
#                                                    parameters['filename'])
#                parameters['logfile'] = full_path_trajectory+'.log'
#            
#            else:
#                if 'trajectory_name'in parameters.keys(): 
#                    if parameters['trajectory_name']:
#                        full_path_trajectory  = os.path.join(parameters['folder'], parameters['trajectory_name'],'output.log' )
#                    else:
#                        full_path_trajectory  = os.path.join(parameters['folder'], 'output.log')
#                else:
#                    full_path_trajectory  = os.path.join(parameters['folder'], 'output.log')
#                
#                parameters['logfile'] =  full_path_trajectory 
#        
#        pprint(parameters)
#        
#         
#
#
#        e_id = system.e_id
#        name = system.label
#        # If a process is already running for this system, abort
#        if e_id in self.process_pool and self.process_pool[e_id][1].is_alive():
#            self.main.simple_dialog.error(
#                msg=f"There is already a process underway for the system:\n {e_id} - {name}"
#            )
#            return False
#
#        # Create new queue for inter-process communication
#        queue = multiprocessing.Queue()
#        parameters['new_vobject'] = True
#        parameters['queue']       = queue
#        
#        
#        # Create and start subprocess
#        process = multiprocessing.Process(
#            target=self._target_process,
#            args=(parameters,)
#        )
#        process.start()
#
#        # Add new entry in bottom_notebook to track progress
#        message = f"{parameters['simulation_type']} {system.e_step_counter} - Running..."
#        hamiltonian = getattr(system.qcModel, 'hamiltonian', 'unk')
#        
#        treeiter = self.main.process_manager_window.add_new_process(system    = system,
#                                                                    _type     = parameters['simulation_type'],
#                                                                    potential = hamiltonian,
#                                                                    status    = 'Queued')
#          
#          
#        
#        
#        path = self.main.bottom_notebook.status_teeview_add_new_item(
#            message=message,
#            system=system
#        )
#
#        # Store process info in process_pool:
#        # [queue, process, result_status, GUI_path]
#        self.process_pool[e_id] = [queue, process, None, path, treeiter]
#
#        # Schedule periodic queue checks
#        GLib.timeout_add(200, self._check_queue)
#        #system.e_treeview_iter   = backup[0]
#        #system.e_liststore_iter  = backup[1]
#        return False
#    
#    
#    def cancel_process(self, e_id):
#        """
#        Cancel (terminate) a running process for the given system ID (e_id).
#        Updates the GUI and cleans up the process pool.
#        """
#        if e_id not in self.process_pool:
#            print(f"No active process found for system {e_id}.")
#            return False
#
#        queue, process, result, path = self.process_pool[e_id]
#
#        if process.is_alive():
#            # Forcefully terminate the process
#            process.terminate()
#            process.join()
#
#            # Update GUI status to reflect cancellation
#            iter_ = self.main.bottom_notebook.status_liststore.get_iter(path)
#            value = self.main.bottom_notebook.status_liststore.get_value(iter_, 1)
#            self.main.bottom_notebook.status_liststore.set_value(iter_, 1, value.replace('Running...', 'Canceled!'))
#
#            # Mark the process as canceled
#            self.process_pool[e_id][2] = "Canceled"
#            print(f"Process for system {e_id} was canceled.")
#
#            return True
#        else:
#            print(f"Process for system {e_id} is not running.")
#            return False
#
#
