#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Lembrar de colocar uma header nesse arquivo

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

import glob, math, os, os.path, sys
import pickle

from datetime import date
import time

import numpy as np

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
from vismol.model.molecular_properties import ATOM_TYPES
from vismol.model.molecular_properties import COLOR_PALETTE

from pdynamo.p_methods import GeometryOptimization
from pdynamo.p_methods import RelaxedSurfaceScan
from pdynamo.p_methods import MolecularDynamics
from pdynamo.p_methods import ChainOfStatesOptimizePath
from pdynamo.p_methods import NormalModes
from pdynamo.LogFileWriter import LogFileReader

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


    def save_easyhybrid_session (self, filename = 'session.easy'):
        """ Function doc """
        easyhybrid_session_data = {}
        
        
        '''- - - - - - - - - - pDynamo systems - - - - - - - - - - - '''
        easyhybrid_session_data['psystem'] = self.psystem
        '''- - - - - - - - - - - - - - - - - - - - - - - - - - - - - '''
        
        
        
        '''- - - - - - - - - - - - vismol obejcts - - - - - - - - - - - '''
        vobjects = {}
        for key, vobject in self.vm_session.vm_objects_dic.items():
            parameters = {
                          'index'             : vobject.index            ,
                          'name'              : vobject.name             ,
                          'active'            : vobject.active           ,
                          'frames'            : vobject.frames           ,
                          'color_palette'     : vobject.color_palette    ,
                          'mass_center'       : vobject.mass_center      ,
                          'selected_atom_ids' : vobject.selected_atom_ids,
                          'index_bonds'       : vobject.index_bonds      ,
                                                
                          'colors'            : vobject.colors           ,
                          'color_indexes'     : vobject.color_indexes    ,
                         }
            vobjects[key] = parameters
        easyhybrid_session_data['vobjects'] = vobjects
        '''- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '''

        
        with open(filename,'wb') as outfile:
            pickle.dump(easyhybrid_session_data, outfile)



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

    
    def _import_coordinates_from_pklfolder (self, parameters):
        files = os.listdir( parameters['data_path'])
        
        pkl_files = []
        for _file in files:
            # Check if the file is a text file
            if _file.endswith('.pkl'):
                pkl_files.append(_file)

        print ('pDynamo pkl folder:' , parameters['data_type'])
        print ('Number of pkl files:', len(pkl_files))
        
        
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
            # . Loop over the frames in the trajectory.
            while trajectory.RestoreOwnerData ( ):
                p_coords = self.psystem[parameters['system_id']].coordinates3
                v_coords = self._convert_pDynamo_coords_to_vismol(p_coords)
                vismol_object.frames = np.vstack((vismol_object.frames, v_coords))

            trajectory.ReadFooter ( )
            trajectory.Close ( )
            
            '''
            for frame in range(0,len(pkl_files)):
                p_coords = ImportCoordinates3 (os.path.join(parameters['data_path'], 'frame{}.pkl'.format(frame)))
                v_coords = self._convert_pDynamo_coords_to_vismol(p_coords)
                parameters['vobject'].frames = np.vstack((parameters['vobject'].frames, v_coords))
            '''
            
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
            print(trajectory)
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
        print (modes_dict)
        return modes_dict

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
        
        elif parameters['data_type'] == 'normal_modes':
            data = self._import_normal_modes_data (parameters)
            self.main.main_treeview.refresh_number_of_frames()
            return data
        else:
            pass

        
        
        print (parameters)
        if parameters['logfile']:
            logfile= LogFileReader(parameters['logfile'])
            data   = logfile.get_data()
        
            #print('vobject', parameters['vobject'], parameters['vobject_id'] )
            
            vobject_id = parameters['vobject_id']
            
            
            if vobject_id in  self.psystem[parameters['system_id']].e_logfile_data.keys():
                self.psystem[parameters['system_id']].e_logfile_data[vobject_id].append(data)
            else:
                self.psystem[parameters['system_id']].e_logfile_data[vobject_id] = []
                self.psystem[parameters['system_id']].e_logfile_data[vobject_id].append(data)
            print ('\n\n\n\n\n\n\n')
            print (self.psystem[parameters['system_id']].e_logfile_data)
        else:
            pass
        
        
        self.main.main_treeview.refresh_number_of_frames()



class pSimulations:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass


    def run_simulation (self, parameters):
        """ Function doc """
        parameters['system'] = self.psystem[self.active_id] 
        print(parameters)
        

        if parameters['simulation_type'] == 'Geometry_Optimization':
            self.geometry_optimization_object = GeometryOptimization()
            self.geometry_optimization_object.run(parameters)
            
        
        elif parameters['simulation_type'] == 'Molecular_Dynamics':
            pprint(parameters)
            self.molecular_dynamics_object = MolecularDynamics()
            self.molecular_dynamics_object.run(parameters)
        
        
        elif parameters['simulation_type'] == 'Relaxed_Surface_Scan':
            self.relaxed_surface_scan = RelaxedSurfaceScan()
            self.relaxed_surface_scan.run(parameters)
            

        
        elif parameters['simulation_type'] == 'Nudged_Elastic_Band':
            self.chain_of_states_opt = ChainOfStatesOptimizePath()
            self.chain_of_states_opt.run(parameters)
        
        
        elif parameters['simulation_type'] == 'Normal_Modes':
            self.normal_modes = NormalModes()
            self.normal_modes.run(parameters)
        
        
        else:
            pass
        
        
        
        self.psystem[self.active_id].e_working_folder == parameters['folder']

        vismol_object = self._add_vismol_object_to_easyhybrid_session(system = self.psystem[self.active_id], 
                                                                        name = parameters['simulation_type']+ ' ' + str(parameters['system'].e_step_counter))        
        
        self.main.bottom_notebook.status_teeview_add_new_item(message = 'New VObejct:  {} '.format(parameters['simulation_type']+ ' ' + str(parameters['system'].e_step_counter)))
                                                                                                                          

        self._apply_fixed_representation_to_vobject(system_id = None, vismol_object = vismol_object)
        self._apply_QC_representation_to_vobject   (system_id = None, vismol_object = vismol_object)

        parameters['system'].e_step_counter += 1
        
class ModifyRepInVismol:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass


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

            self.psystem[system_id].e_qc_table = list(self.psystem[system_id].qcState.pureQCAtoms)
            
            
            for atom in vismol_object.atoms.values():
                atom.spheres = False
                atom.sticks = False
                    
                
            selection = self.vm_session.create_new_selection()
            selection.selecting_by_indexes(vismol_object, self.psystem[system_id].e_qc_table, clear=True)
            
            '''
            for _id , vismol_object in self.vm_session.vm_objects_dic.items():
                selection.selecting_by_indexes(vismol_object, self.psystem[system_id].e_qc_table, clear=False)
            
            #'''

            self.vm_session.show_or_hide (rep_type = 'spheres',selection= selection ,  show = True )
            #self.vm_session.show_or_hide (rep_type = 'spheres', show = True )

            if static:

                self.vm_session.show_or_hide(rep_type = 'sticks',selection= selection , show = True )
                #self.vm_session.show_or_hide(rep_type = 'sticks' , show = True )

            else:
                pass
                #self.vm_session.show_or_hide(_type = 'dynamic_bonds' , vobject = vobject, selection_table = self.systems[system_id]['qc_table'] , show = True )

            for atom in vismol_object.atoms.values():
                atom.selected  =  False
                atom.vm_object.selected_atom_ids.discard(atom.atom_id)


class Restraints:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
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

    def define_distance_farmonic_restraints (self, parameters):
        """ Function doc """
        restraints = RestraintModel()
        
        parameters['system'].DefineRestraintModel( restraints )
        
        distance = parameters['dminimum_RC1'] + ( parameters['dincre_RC1'] * float(i) )
        
        rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['force_constant_1'])
        restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom1, point2= atom2)
        restraints["RC1"] = restraint

    def clear_restraints (self):
        """ Function doc """
        parameters['system'].DefineRestraintModel(None)


class pDynamoSession (pSimulations, ModifyRepInVismol, LoadAndSaveData, EasyHybridImportTrajectory):
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
        
        
    def load_a_new_pDynamo_system_from_dict (self, input_files = {}, system_type = 0, name = None, tag = None):
        """ Function doc """

        #psystem = self.generate_pSystem_dictionary()
        
        system = None 
        if system_type == 0:
            system              = ImportSystem       ( input_files['amber_prmtop'] )
            system.coordinates3 = ImportCoordinates3 ( input_files['coordinates'] )
            self.define_NBModel(_type = 1, system = system)                      
        
        elif system_type == 1:
            parameters          = CHARMMParameterFileReader.PathsToParameters (input_files['charmm_par'])
            system              = ImportSystem       ( input_files['charmm_psf'] , isXPLOR = True, parameters = parameters )
            system.coordinates3 = ImportCoordinates3 ( input_files['coordinates'] )
            self.define_NBModel(_type = 1, system = system)        
        
        elif system_type == 2:
            mmModel        = MMModelOPLS.WithParameterSet ( input_files['opls_folder'] )       
            system         = ImportSystem ( input_files['coordinates'])
            system.DefineMMModel ( mmModel )
            self.define_NBModel(_type = 1, system = system)          
        
        elif system_type == 5:
            mmModel        = MMModelDYFF.WithParameterSet ( input_files['opls_folder'] )       
            system         = ImportSystem ( input_files['coordinates'])
            system.DefineMMModel ( mmModel )
            self.define_NBModel(_type = 1, system = system)          
        
        elif system_type == 3 or system_type == 4 :
            system = ImportSystem (input_files['coordinates'])

        else:
            pass
        system.Summary()
        
        
        
        '''-----------------------------------------------------------------'''
        
        name =  self._check_name(name)
        
        system.e_input_files = input_files
        system = self.append_system_to_pdynamo_session ( 
                                                        system = system,  name = name, tag = tag )
        
        '''-----------------------------------------------------------------'''
        
        #system
        self.main.main_treeview.add_new_system_to_treeview (system)
        
        #sys_type = {0: 'AMBER', 1:'CHARMM'}
        ff  =  getattr(system.mmModel, 'forceField', "None")
        self.main.bottom_notebook.status_teeview_add_new_item(message = 'New System:  {} ({}) - Force Field:  {}'.format(system.label, system.e_tag, ff))

        self._add_vismol_object_to_easyhybrid_session (system, True) #, name = 'olha o  coco')
        #self.main.refresh_active_system_liststore()
        #self.main.refresh_system_liststore ()
        
    
    def _add_vismol_object_to_easyhybrid_session (self, system, show_molecule=True, name = 'new_coords'):
        """ Function doc """
        # vismol obejct
        vm_object = self._build_vobject_from_pDynamo_system ( system = system, name = name ) 
        self.vm_session._add_vismol_object(vm_object, show_molecule=True)
        self.main.main_treeview.add_vismol_object_to_treeview(vm_object)
        
        self.main.add_vobject_to_vobject_liststore_dict(vm_object)
        
        self._apply_fixed_representation_to_vobject(vismol_object =vm_object)
        self._apply_QC_representation_to_vobject(vismol_object =vm_object)
        
        self.main.refresh_widgets()
        return vm_object

    
    def append_system_to_pdynamo_session (self, 
                                          system         = None            ,
                                          name           = 'pDynamo system',
                                          tag            = 'molsys'        ,
                                          working_folder = None            ,
                                                        ):
        """ Function doc """
        
        if name:
            system.label = name
        else:
            pass
        
        if tag:
            system.e_tag = tag
        else:
            system.e_tag ='MolSys'
        
        system.e_working_folder           = working_folder
        
        try:
            system.e_charges_backup           = list(system.AtomicCharges()).copy()
        except:
            system.e_charges_backup           = []

        
        system.e_id                      =  self.counter
        self.psystem[system.e_id]        =  system 
        
        self.active_id   = self.counter  
        self.counter    += 1
        
        
        system.e_active                   = False
        system.e_date                     = date.today()               # Time     
        system.e_color_palette            = self.color_palette_counter # will be replaced by a dict
        
        
        '''When the number of available colors runs out, we have to reset the counter'''
        self.color_palette_counter        += 1
        if self.color_palette_counter > len(COLOR_PALETTE.keys())-1:
            self.color_palette_counter = 0
        else:
            pass
        
        #system.e_vobject                  = None            # Vismol object associated with the system -> is the object that will 

        system.e_active                   = False          
        system.e_bonds                    = None           
        system.e_sequence                 = None           

        #system.e_QCMM                     = False
        system.e_qc_table                 = []           
        system.e_qc_residue_table         = {}             #yes, it's a dict           
        system.e_fixed_table              = []             
        system.e_selections               = {}             
        system.e_selections_counter       = 0             

        system.e_vm_objects               = {}                  
        system.e_logfile_data             = {}             # <--- vobject_id : [data0, data1, data2], ...]  obs: each "data" is a dictionary 
        system.e_step_counter             = 0              
        
        '''Now each pdynamo system object and vismol object has a 
        "treeview_iter" attribute, which is a reference to access 
        the treeview elements
        
        the method "add_new_system_to_treeview" 
        
        '''
        system.e_treeview_iter            = None  #the method "add_new_system_to_treeview" modify this attribute 
        system.e_liststore_iter           = None  #the method "add_new_system_to_treeview" modify this attribute
        
        
        
        return system
        
        
        # obtaing the color palette
        '''
        if self.color_palette_counter >= len(COLOR_PALETTE)-1:
            self.color_palette_counter = 0
        else:
            self.color_palette_counter += 1
        '''
        

    def _get_sequence_from_pDynamo_system (self, system = None):
        """ Function doc """
        
        if system:
            sequence = getattr ( system, "sequence", None )
            
            if sequence is None: 
                sequence = Sequence.FromAtoms (system.atoms, componentLabel = "UNK.1" )
            else:
                pass
        else:
            sequence = None
        
        return sequence


    def _get_atom_info_from_pdynamo_atom_obj (self, sequence = None, atom = None):
        """
        To extract information from atom objects, 
        belonging to pdynamo, and to  organize it as a list
        that will be delivered later to build the 
        vismolObj
        
        """
        entityLabel = atom.parent.parent.label
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
        
        ###print(atom.index, atom.label,resName, resSeq, iCode,chainID, segID,  atom.atomicNumber, atom.connections)#, xyz[0], xyz[1], xyz[2] )
        
        index        = atom.index+1
        at_name      = atom.label
        at_resi      = int(resSeq)
        at_resn      = resName
        at_ch        = chainID
        try:
            at_symbol    = ATOM_TYPES_BY_ATOMICNUMBER[atom.atomicNumber] # at.get_symbol(at_name)
        except:
            at_symbol = "O"
        
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
        
        return atom


    def _get_qc_table_from_pDynamo_system(self, system):
        system.e_qc_table = list(system.qcState.pureQCAtoms)


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
        """ Function doc """
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
    
    def get_color_palette (self, system_id = None):
        """ Function doc """
        if system_id:
            system = self.psystem[system_id]
            return COLOR_PALETTE[system.e_color_palette]
        else:
            system = self.psystem[self.active_id]
            return COLOR_PALETTE[system.e_color_palette]
    
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
        sequence  = self._get_sequence_from_pDynamo_system(system)
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

            atoms.append(self._get_atom_info_from_pdynamo_atom_obj(sequence =  sequence, atom   = atom))
            j += 1
            
        vm_object = VismolObject(self.vm_session, 
                                 len(self.vm_session.vm_objects_dic), 
                                 name = name, 
                                 color_palette = COLOR_PALETTE[system.e_color_palette])
                                 
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
                _r = Residue(vm_object, name=_atom["resn"], index=_atom["resi"], chain=_chain)
                vm_object.residues[_atom["resi"]] = _r
                _chain.residues[_atom["resi"]] = _r
            _residue = _chain.residues[_atom["resi"]]
            
            #print(_atom)
            #atom = Atom()
            atom = Atom(vismol_object = vm_object,#
                        name          =_atom["name"] , 
                        index         =_atom["index"],
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
        vm_object.find_bonded_and_nonbonded_atoms()
        
        vm_object.e_id = system.e_id
        vm_object._generate_color_vectors(self.vm_session.atom_id_counter)
        vm_object.active = True
        #vm_object.e_treeview_iter_parent_key = None
        
        '''-----------------------------------------------------'''
        '''Now each pdynamo system object and vismol object has a 
        "treeview_iter" / "liststore_iter" attribute, which is a reference to access 
        the main treeview elements  and the self.main.vobject_liststore_dict[sys_e_id]'''
        vm_object.e_treeview_iter = None
        vm_object.liststore_iter = None
        
        '''When an object is removed it has to be removed from the treeview and 
        vobject_liststore_dict, in addition to the vm_object_dic in the .vm_session.'''
        
        '''-----------------------------------------------------'''
        
        return vm_object


    def delete_system (self, system_e_id = None):
        """ Function doc """
        if system_e_id != None:
            
            if self.psystem[system_e_id].label in self.psystem_name_list:
                index = self.psystem_name_list.index(self.psystem[system_e_id].label)
                self.psystem_name_list.pop(index)
            
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


    def prune_system (self, selection = None, name = 'Pruned System', summary = True, tag = "molsys"):
        """ Function doc """
        p_selection   = Selection.FromIterable ( selection )
        system        = PruneByAtom ( self.psystem[self.active_id], p_selection )
        
        self.define_NBModel(_type = 1, system = system)

        
        
        system.label  = name        
        if summary:
            system.Summary ( )
            
        system = self.append_system_to_pdynamo_session ( 
                                                        system = system,  name = name, tag = tag )
        
        '''-----------------------------------------------------------------'''
        #system
        self.main.main_treeview.add_new_system_to_treeview (system)
        
        self._add_vismol_object_to_easyhybrid_session (system, True)
        self.main.refresh_active_system_liststore()
        self.main.refresh_system_liststore ()


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
        
        if system:
            pass
        else:
            system = self.psystem[self.active_id]
        
        electronicState = ElectronicState.WithOptions ( charge = parameters['charge'], 
                                                  multiplicity = parameters['multiplicity'], 
                                              isSpinRestricted =  parameters['isSpinRestricted'])
                                              
        system.electronicState = electronicState
        
        converger = DIISSCFConverger.WithOptions( energyTolerance   = parameters['energyTolerance'] ,
                                                  densityTolerance  = parameters['densityTolerance'] ,
                                                  maximumIterations = parameters['maximumIterations'] )


        qcModel         = QCModelMNDO.WithOptions ( hamiltonian = parameters['method'], converger=converger )
        
        if len(system.e_qc_table) > 0:
            
            '''This function rescales the remaining charges in the MM part. The 
            sum of the charges in the MM region must be an integer value!'''
            self.check_charge_fragmentation(vismol_object = vismol_object)
            '''----------------------------------------------------------------'''
            system.DefineQCModel (qcModel, qcSelection = Selection.FromIterable ( system.e_qc_table) )          
            
            if system.mmModel:
                system.DefineNBModel ( NBModelCutOff.WithDefaults ( ) )
            else:
                pass
        
        else:
            system.DefineQCModel (qcModel)
                       
        system.Summary()
        self.main.refresh_widgets()
        for vismol_object in self.vm_session.vm_objects_dic.values():
            if vismol_object.e_id == system.e_id:
                self._apply_QC_representation_to_vobject   (system_id = None, vismol_object = vismol_object)
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
        if vismol_object:
            
            
            for resi in vismol_object.residues.keys():
                
                if resi in qc_residue_table.keys():
                    
                    mm_residue_table[resi] = []
                    
                    for index, atom in vismol_object.residues[resi].atoms.items():
                        index_v = atom.index-1
                        charge  = system.mmState.charges[index_v]
                        resn    = vismol_object.residues[resi].name 
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
        '''Here we are going to do a charge rescaling of atoms within the MM part 
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
        

    def export_system (self,  parameters ): 
                              
        """  
        Export system model, as pDynamo serization files or Cartesian coordinates.
            0 : 'pkl - pDynamo system'         ,
            1 : 'pkl - pDynamo coordinates'    ,
            2 : 'pdb - Protein Data Bank'      ,
            3 : 'xyz - cartesian coordinates'  ,
            4 : 'mol'                          ,
            5 : 'mol2'                         ,
        
        """
        vobject  = self.vm_session.vm_objects_dic[parameters['vobject_id']]
        folder   = parameters['folder'] 
        filename = parameters['filename'] 
        
        
        if parameters['last'] == -1:
            parameters['last'] = len(vobject.frames)-1
        
        active_id = self.active_id 
        self.active_id = parameters['system_id']
        
        if parameters['format'] == 0:
            self.get_coordinates_from_vobject_to_pDynamo_system(vobject       = vobject, 
                                                                system_id     = parameters['system_id'], 
                                                                frame         = parameters['last'])
            
            system = self.psystem[parameters['system_id']]
            
            '''
            The pkl format is not capable of exporting GTK or openGL 
            elements (objects).
            '''
            backup = []
            backup.append(system.e_treeview_iter)
            backup.append(system.e_liststore_iter)
            
            system.e_treeview_iter   = None
            system.e_liststore_iter  = None

            ExportSystem ( os.path.join ( folder, filename+'.pkl'), system )

            system.e_treeview_iter   = backup[0]
            system.e_liststore_iter  = backup[1]
        
        
        
        if parameters['format'] == 1:
            
            #'''   When is Single File     '''
            if parameters['is_single_file']:
                self.get_coordinates_from_vobject_to_pDynamo_system(vobject = vobject, 
                                                                          system_id     = parameters['system_id'], 
                                                                          frame         = parameters['last'])
                
                system   = self.psystem[parameters['system_id']]
                
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
                
                
                i = 0
                for frame in range(parameters['first'], parameters['last'], parameters['stride']):
                    
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
                
                if parameters['format'] == 2:
                    ExportSystem ( os.path.join ( folder, filename+'.pdb'), system )
                
                if parameters['format'] == 3:
                    ExportSystem ( os.path.join ( folder, filename+'.xyz'), system )
                
                if parameters['format'] == 4:
                    ExportSystem ( os.path.join ( folder, filename+'.mol'), system )
                
                if parameters['format'] == 5:
                    ExportSystem ( os.path.join ( folder, filename+'.mol2'), system )
                    

            
            
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
                
                
                i = 0
                for frame in range(parameters['first'], parameters['last'], parameters['stride']):
                    
                    self.get_coordinates_from_vobject_to_pDynamo_system(vobject   = vobject, 
                                                                        system_id = parameters['system_id'], 
                                                                        frame     = frame)

                    system   = self.psystem[parameters['system_id']]
                    
                    if parameters['format'] == 2:
                        ExportSystem ( os.path.join ( folder, filename+'.pdb'), system )
                    
                    if parameters['format'] == 3:
                        ExportSystem ( os.path.join ( folder, filename+'.xyz'), system )
                    
                    if parameters['format'] == 4:
                        ExportSystem ( os.path.join ( folder, filename+'.mol'), system )
                    
                    if parameters['format'] == 5:
                        ExportSystem ( os.path.join ( folder, filename+'.mol2'), system )
                    
                    
                    #Pickle( os.path.join ( folder, "frame{}.pkl".format(i) ), 
                    #        system.coordinates3 )
                    
                    i += 1
        


        
        self.active_id  = active_id

    
    def _convert_pDynamo_coords_to_vismol(self, p_coords):
        '''
        This function converts pDynamo coordinates ( which is a l
        ist containg all the coords) into a vObject coordinates 
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
        print (coords)
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
    
    def __init__(self, vismol_object=None, name="Xx", index=None, residue=None,
                 chain=None, pos=None, symbol=None, atom_id=None, color=None,
                 vdw_rad=None, cov_rad=None, ball_rad=None,
                 occupancy=0.0, bfactor=0.0, charge=0.0, bonds_indexes=None):
        """ Class initializer """
        self.vm_object = vismol_object
        self.name = name
        self.index = index   # - Remember that the "index" attribute refers to the numbering of atoms (it is not a zero base, it starts at 1 for the first atom)
        self.residue = residue
        self.chain = chain
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
        self.ball_and_stick = False
        self.sticks         = False
        self.spheres        = False
        self.dash           = False
        self.surface        = False
        self.bonds          = []
        self.isfree         = True
    
    def _get_symbol(self):
        """ Function doc """
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
                elif name[1] =="d":
                    symbol = "Cd"
                elif name[1] =="u":
                    symbol = "Cu"
                else:
                    symbol = "C"
            
            elif name[0] == "N":
                if name[1] == "i":
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
        """ Function doc """
        r = (self.unique_id & 0x000000FF) >> 0
        g = (self.unique_id & 0x0000FF00) >> 8
        b = (self.unique_id & 0x00FF0000) >> 16
        # pickedID = r + g * 256 + b * 256*256
        self.color_id = np.array([r/255.0, g/255.0, b/255.0], dtype=np.float32)
    
    def _init_vdw_rad(self):
        """ Function doc """
        try:
            vdw = ATOM_TYPES[self.name][6]
        except KeyError:
            vdw = ATOM_TYPES[self.symbol][6]
        return vdw
    
    def _init_cov_rad(self):
        """ Function doc """
        try:
            cov = ATOM_TYPES[self.name][5]
        except KeyError:
            cov = ATOM_TYPES[self.symbol][5]
        return cov
    
    def _init_ball_rad(self):
        """ Function doc """
        try:
            ball = ATOM_TYPES[self.name][6]
        except KeyError:
            ball = ATOM_TYPES[self.symbol][6]
        return ball
    
    def coords(self, frame=None):
        """ Function doc """
        if frame is None:
            frame  = self.vm_object.vm_session.frame
        # x = self.vm_object.frames[frame][(self.index-1)*3  ]
        # y = self.vm_object.frames[frame][(self.index-1)*3+1]
        # z = self.vm_object.frames[frame][(self.index-1)*3+2]
        # return np.array([x, y, z])
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





