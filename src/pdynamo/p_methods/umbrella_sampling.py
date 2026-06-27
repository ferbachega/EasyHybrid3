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


#from LogFile import LogFileWriter
# pDynamo
from pBabel                    import *                                     
from pCore                     import *                                     
from pMolecule                 import *                  
from pScientific               import *                                     
from pScientific.Arrays        import *                                     
from pScientific.Geometry3     import *                 
from pSimulation               import *
#*********************************************************************************
import multiprocessing
import copy

from pScientific.RandomNumbers import NormalDeviateGenerator                       , \
                                      RandomNumberGenerator
import json
from pprint import pprint
import os, time, sys

# --- imports entre modulos adicionados na refatoracao ---
from pdynamo.p_methods._common import write_header, get_hamiltonian

class UmbrellaSampling:
    def __init__ (self):
        """ Class initialiser """
        pass

    def write_header (self, parameters, logfile = 'output.log'):
        """ Function doc """
        
        arq = open(os.path.join( parameters['folder'], parameters['traj_folder_name'], logfile), "a")
        text = ""
        
        '''
        if parameters["RC2"] is not None:
            text = text + "\n"
            text = text + "\n--------------------------------------------------------------------------------"
            text = text + "\nTYPE                         EasyHybrid-SCAN2D                                  "
            text = text + "\n--------------------------------------------------------------------------------"

        else:
            text = text + "\n"
            text = text + "\n--------------------------------------------------------------------------------"
            text = text + "\nTYPE                          EasyHybrid-SCAN                                   "
            text = text + "\n--------------------------------------------------------------------------------"
        '''



        '''This part writes the parameters used in the first reaction coordinate'''
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------
        '''
        if parameters['RC1']["rc_type"] == 'simple_distance':
            text = text + "\n"
            text = text + "\n----------------------- Coordinate 1 - Simple-Distance -------------------------"
            text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC1']['ATOMS'][0], parameters['RC1']['ATOM_NAMES'][0] )
            text = text + "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['RC1']['ATOMS'][1], parameters['RC1']['ATOM_NAMES'][1] )
            text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['RC1']['nsteps']  , parameters['RC1']['force_constant'] )
            text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC1']['dminimum'], parameters['OPT_parm']['maximumIterations']         )
            text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC1']['dincre']  , parameters['OPT_parm']['rmsGradientTolerance']           )
            text = text + "\n--------------------------------------------------------------------------------"

        
        elif parameters['RC1']["rc_type"] == 'multiple_distance':
            text = text + "\n"
            text = text + "\n---------------------- Coordinate 1 - multiple-Distance ------------------------"	
            text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC1']['ATOMS'][0]    , parameters['RC1']['ATOM_NAMES'][0] )
            text = text + "\nATOM2*                 =%15i  ATOM NAME2             =%15s"     % (parameters['RC1']['ATOMS'][1]    , parameters['RC1']['ATOM_NAMES'][1] )
            text = text + "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC1']['ATOMS'][2]    , parameters['RC1']['ATOM_NAMES'][2] )
            text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['RC1']['nsteps']      , parameters['RC1']['force_constant'] ) 
            text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC1']['dminimum']    , parameters['OPT_parm']['maximumIterations']         )
            text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC1']['dincre']      , parameters['OPT_parm']['rmsGradientTolerance']           )
            text = text + "\nSigma atom1 - atom3    =%15.5f  Sigma atom3 - atom1    =%15.5f" % (parameters['RC1']['sigma_pk1pk3'], parameters['RC1']['sigma_pk3pk1']   )
            text = text + "\n--------------------------------------------------------------------------------"

        else:
            pass
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------
        '''
        
        
        
        #if parameters["RC2"] is not None :
        #
        #    '''This part writes the parameters used in the second reaction coordinate'''
        #    #-----------------------------------------------------------------------------------------------------------------------------------------------------------
        #    if parameters['RC2']["rc_type"] == 'simple_distance':
        #        text = text + "\n"
        #        text = text + "\n----------------------- Coordinate 2 - Simple-Distance -------------------------"
        #        text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC2']['ATOMS'][0], parameters['RC2']['ATOM_NAMES'][0] )
        #        text = text + "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['RC2']['ATOMS'][1], parameters['RC2']['ATOM_NAMES'][1] )
        #        text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['RC2']['nsteps']  , parameters['RC2']['force_constant'] )
        #        text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC2']['dminimum'], parameters['maximumIterations']         )
        #        text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC2']['dincre']  , parameters['rmsGradientTolerance']           )
        #        text = text + "\n--------------------------------------------------------------------------------"
        #
        #    
        #    elif parameters['RC2']["rc_type"] == 'multiple_distance':
        #        text = text + "\n"
        #        text = text + "\n---------------------- Coordinate 2 - multiple-Distance ------------------------"	
        #        text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC2']['ATOMS'][0]    , parameters['RC2']['ATOM_NAMES'][0] )
        #        text = text + "\nATOM2*                 =%15i  ATOM NAME2             =%15s"     % (parameters['RC2']['ATOMS'][1]    , parameters['RC2']['ATOM_NAMES'][1] )
        #        text = text + "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC2']['ATOMS'][2]    , parameters['RC2']['ATOM_NAMES'][2] )
        #        text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['RC2']['nsteps']      , parameters['RC2']['force_constant'] )
        #        text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC2']['dminimum']    , parameters['maximumIterations']         )
        #        text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC2']['dincre']      , parameters['rmsGradientTolerance']           )
        #        text = text + "\nSigma atom1 - atom3    =%15.5f  Sigma atom3 - atom1    =%15.5f" % (parameters['RC2']['sigma_pk1pk3'], parameters['RC2']['sigma_pk3pk1']   )
        #        text = text + "\n--------------------------------------------------------------------------------"
        #    else:
        #        pass
        #    #-----------------------------------------------------------------------------------------------------------------------------------------------------------
        #    
        #
        #
        #
        #
        #
        #
        #
        ##-----------------------------------------------------------------------------------------------------------------------------------------------------------
        #'''This part writes the header of the frames distances and energy'''
        #if parameters["RC2"] is not None :
        #    text = text + "\n\n--------------------------------------------------------------------------------"
        #    text = text + "\n   Frame i  /  j        RCOORD-1             RCOORD-2                Energy     "
        #    text = text + "\n--------------------------------------------------------------------------------"
        #
        #
        #
        #
        #else:
        #    if parameters['RC1']["rc_type"] == 'simple_distance':
        #        text = text + "\n\n-------------------------------------------------------------"
        #        text = text + "\n           Frame    dist-ATOM1-ATOM2             Energy      "
        #        text = text + "\n-------------------------------------------------------------"
        #    
        #    elif parameters['RC1']["rc_type"] == 'multiple_distance':
        #        text = text + "\n\n--------------------------------------------------------------------------------"
        #        text = text + "\n           Frame     dist-ATOM1-ATOM2      dist-ATOM2-ATOM3         Energy        "
        #        text = text + "\n--------------------------------------------------------------------------------  "
        ##-----------------------------------------------------------------------------------------------------------------------------------------------------------

            
        
        
        
            
            
            #text = text + "\n\n------------------------------------------------------"
            #text = text + "\n       Frame     distance-pK1-pK2         Energy      "
            #text = text + "\n------------------------------------------------------"
        
        arq.write(text)
        return arq

    def write_parameters (self, parameters, logfile = None):
        """ Function doc """
        ##logfile = open( os.path.join(full_path_trajectory, 'output.log'), 'a')
        #json_data = json.dumps(parameters, indent=4)
        #with open(logfile, "a") as file:
        #    file.write(json_data)

        with open(logfile, 'a') as file:
            sys.stdout = file  # Redirect stdout to the file
            pprint(parameters)
            sys.stdout = sys.__stdout__  # Restore stdout

    def run (self, parameters, interface = False):
        """ Function doc """
        
        #this a temorary code- change it later
        
        if 'RC' in parameters['RC1'].keys():
            pass
        else:
            #----------------------------------------------------------------------------------------    
            atom1 = parameters['RC1']['ATOMS'][0]
            atom2 = parameters['RC1']['ATOMS'][1]                   
            atom1_name = parameters['RC1']['ATOM_NAMES'][0]
            atom2_name = parameters['RC1']['ATOM_NAMES'][1]
            #----------------------------------------------------------------------------------------    
            if parameters['RC1']["rc_type"] == 'simple_distance':
                parameters['RC1']['RC'] = [[atom1_name, atom1, atom2_name, atom2, 1.0, 1.0]]
                #distance = system.coordinates3.Distance( atom1, atom2)
                #restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom1, point2= atom2)

            elif parameters['RC1']["rc_type"] == 'multiple_distance':

                atom3      = parameters['RC1']['ATOMS'][2]
                atom3_name = parameters['RC1']['ATOM_NAMES'][2]
                sigma_pk1pk3 = parameters['RC1']['sigma_pk1pk3']
                sigma_pk3pk1 = parameters['RC1']['sigma_pk3pk1']
                parameters['RC1']['RC'] = [[atom1_name, atom1, atom2_name, atom2, sigma_pk1pk3, 1.0]]
                parameters['RC1']['RC'].append([atom2_name, atom2, atom3_name, atom3, sigma_pk3pk1, 1.0])
 
            elif parameters['RC1']["rc_type"] == 'multiple_distance*4atoms':
                atom3      = parameters['RC1']['ATOMS'][2]
                atom3_name = parameters['RC1']['ATOM_NAMES'][2]

                atom4      = parameters['RC1']['ATOMS'][3]
                atom4_name = parameters['RC1']['ATOM_NAMES'][3]

                weight1 =  1.0#parameters['RC1']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
                weight2 = -1.0#parameters['RC1']['sigma_pk3pk1'] #self.sigma_a3_a1[0] 
                parameters['RC1']['RC'] = [[atom1_name, atom1, atom2_name, atom2, weight1, 1.0]]
                parameters['RC1']['RC'].append([atom3_name, atom3, atom4_name, atom4, weight2, 1.0])
                #--------------------------------------------------------------------
            else:
                pass        

        if parameters['RC2'] is not None:
            if 'RC' in parameters['RC2'].keys():
                pass
            else:
                #----------------------------------------------------------------------------------------    
                atom1 = parameters['RC2']['ATOMS'][0]
                atom2 = parameters['RC2']['ATOMS'][1]                   
                atom1_name = parameters['RC2']['ATOM_NAMES'][0]
                atom2_name = parameters['RC2']['ATOM_NAMES'][1]
                #----------------------------------------------------------------------------------------    
                if parameters['RC2']["rc_type"] == 'simple_distance':
                    parameters['RC2']['RC'] = [[atom1_name, atom1, atom2_name, atom2, 1.0, 1.0]]
                    #distance = system.coordinates3.Distance( atom1, atom2)
                    #restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom1, point2= atom2)

                elif parameters['RC2']["rc_type"] == 'multiple_distance':

                    atom3      = parameters['RC2']['ATOMS'][2]
                    atom3_name = parameters['RC2']['ATOM_NAMES'][2]
                    sigma_pk1pk3 = parameters['RC2']['sigma_pk1pk3']
                    sigma_pk3pk1 = parameters['RC2']['sigma_pk3pk1']
                    parameters['RC2']['RC'] = [[atom1_name, atom1, atom2_name, atom2, sigma_pk1pk3, 1.0]]
                    parameters['RC2']['RC'].append([atom2_name, atom2, atom3_name, atom3, sigma_pk3pk1, 1.0])
     
                elif parameters['RC2']["rc_type"] == 'multiple_distance*4atoms':
                    atom3      = parameters['RC2']['ATOMS'][2]
                    atom3_name = parameters['RC2']['ATOM_NAMES'][2]

                    atom4      = parameters['RC2']['ATOMS'][3]
                    atom4_name = parameters['RC2']['ATOM_NAMES'][3]

                    weight1 =  1.0#parameters['RC2']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
                    weight2 = -1.0#parameters['RC2']['sigma_pk3pk1'] #self.sigma_a3_a1[0] 
                    parameters['RC2']['RC'] = [[atom1_name, atom1, atom2_name, atom2, weight1, 1.0]]
                    parameters['RC2']['RC'].append([atom3_name, atom3, atom4_name, atom4, weight2, 1.0])
                    #--------------------------------------------------------------------
                else:
                    pass        
               
        full_path_trajectory = os.path.join(parameters['folder'], 
                                            parameters['traj_folder_name'])

        os.mkdir(
                 full_path_trajectory
                 )
        

        # - - - - - - - - - - - - - Checking trajectory - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        self.logFile2 = TextLogFileWriter.WithOptions ( path = os.path.join(full_path_trajectory, 'output.log') )
        parameters['system'].Summary(log = self.logFile2)
        self.logFile2.Header ( )
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        
        self.write_parameters(parameters = parameters, logfile = os.path.join(full_path_trajectory, 'output.log'))

        #----------------------------------------------------------------------------------------
        full_path_trajectory_opt = os.path.join(parameters['folder'], parameters['traj_folder_name'], 'optmization.ptGeo')
        isExist = os.path.exists(full_path_trajectory_opt)
        if isExist:
            pass
        else:
            os.mkdir(full_path_trajectory_opt)
        parameters['trajectory_path_opt'] = full_path_trajectory_opt
        #----------------------------------------------------------------------------------------
        
        #----------------------------------------------------------------------------------------
        full_path_trajectory_eq = os.path.join(parameters['folder'], parameters['traj_folder_name'], 'equilibration')
        isExist = os.path.exists(full_path_trajectory_eq)
        if isExist:
            pass
        else:
            os.mkdir(full_path_trajectory_eq)
        parameters['trajectory_path_eq'] = full_path_trajectory_eq
        #----------------------------------------------------------------------------------------
        
        
        #----------------------------------------------------------------------------------------
        full_path_trajectory_dc = os.path.join(parameters['folder'], parameters['traj_folder_name'], 'data_collection')
        isExist = os.path.exists(full_path_trajectory_dc)
        if isExist:
            pass
        else:
            os.mkdir(full_path_trajectory_dc)
        parameters['trajectory_path_dc'] = full_path_trajectory_dc
        #----------------------------------------------------------------------------------------
        
        
        #----------------------------------------------------------------------------------------
        full_path_trajectory_dc = os.path.join(parameters['folder'], parameters['traj_folder_name'], 'data_collection_traj')
        isExist = os.path.exists(full_path_trajectory_dc)
        if isExist:
            pass
        else:
            os.mkdir(full_path_trajectory_dc)
        parameters['trajectory_path_dc_traj'] = full_path_trajectory_dc
        #----------------------------------------------------------------------------------------

        #pprint (parameters)
        
        ##-------------------------------------------------------------------------
        #if interface:
        #    e_treeview_iter  = getattr(parameters['system'], 'e_treeview_iter', None)
        #    e_liststore_iter = getattr(parameters['system'], 'e_liststore_iter', None)
        #    parameters['system'].e_treeview_iter   = None
        #    parameters['system'].e_liststore_iter  = None            
        ##-------------------------------------------------------------------------
        
        if parameters['RC2'] is not None:
            #self._run_umbrella_sampling_2D(parameters = parameters, interface = False)
            self._run_umbrella_sampling_2D(parameters = parameters, interface = False)
        else:
            #self._run_umbrella_sampling_1D(parameters)
            self._run_umbrella_sampling_1D(parameters)
        
        ##-------------------------------------------------------------------------
        #if interface:
        #    parameters['system'].e_treeview_iter  = e_treeview_iter 
        #    parameters['system'].e_liststore_iter = e_liststore_iter 
        ##-------------------------------------------------------------------------


        '''
        if parameters['RC2'] is not None:
            self._run_umbrella_sampling_2D(parameters = parameters, interface = False)
        else:
            self._run_umbrella_sampling_1D(parameters)
        '''
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        #self.logFile2.Footer ( )
        self.logFile2.Close()
        self.logFile2 = None
 
    def _run_umbrella_sampling (self, parameters, interface = False):
        """ Function doc """
        
        if parameters['input_type'] == 1: # means = parallel
            #------------------------------------------------------------------
            files = os.listdir( parameters['source_folder'])
            pkl_files = []
            for _file in files:
                # Check file extention
                if _file.endswith('.pkl'):
                    pkl_files.append(_file)
            pkl_files.sort()
    
            
            if interface:
                #-------------------------------------------------------------------------
                e_treeview_iter  = getattr(parameters['system'], 'e_treeview_iter', None)
                e_liststore_iter = getattr(parameters['system'], 'e_liststore_iter', None)
                parameters['system'].e_treeview_iter   = None
                parameters['system'].e_liststore_iter  = None            
                #-------------------------------------------------------------------------
            else:
                pass
            
            #-------------------------------------------------------------------------

            
            if parameters['RC2']:
                
                #-------------------------------------------------------------------------
                joblist = []
                for pkl in pkl_files:
                    #  - - - - adding to the joblist - - - - 
                    system  = parameters['system']
                    
                    i_j = pkl[5:-4] 
                    i_j = i_j.split('_')
                    
                    i = int(i_j[0])
                    j = int(i_j[1])
                                   
                    joblist.append([i, j, pkl, parameters, parameters['system']])
                
                
                p = multiprocessing.Pool(processes=parameters['NmaxThreads'])
                results = p.map(_run_parallel_umbrella_sampling_2D, joblist)
                #-------------------------------------------------------------------------

                if interface:
                    parameters['system'].e_treeview_iter  = e_treeview_iter 
                    parameters['system'].e_liststore_iter = e_liststore_iter 





            else:
                joblist = []
                for pkl in pkl_files:
                    #  - - - - adding to the joblist - - - - 
                    system  = parameters['system']
                    i = int(pkl[5:-4]) 
                    joblist.append([i, pkl, parameters, parameters['system']])
                    
                p = multiprocessing.Pool(processes=parameters['NmaxThreads'])
                results = p.map(_run_parallel_umbrella_sampling_1D, joblist)

                #-------------------------------------------------------------------------
                
                if interface:
                    parameters['system'].e_treeview_iter  = e_treeview_iter 
                    parameters['system'].e_liststore_iter = e_liststore_iter 
            
        
        else:
            self._run_serial_umbrella_sampling_1D (parameters)
    
    def _run_umbrella_sampling_1D (self, parameters):
        """ Function doc """
        '''
        #----------------------------------------------------------------------------------------
        full_path_trajectory_eq = os.path.join(parameters['folder'], parameters['traj_folder_name'], 'equilibration')
        isExist = os.path.exists(full_path_trajectory_eq)
        if isExist:
            pass
        else:
            os.mkdir(full_path_trajectory_eq)
        parameters['trajectory_path_eq'] = full_path_trajectory_eq
        #----------------------------------------------------------------------------------------
        
        
        #----------------------------------------------------------------------------------------
        full_path_trajectory_dc = os.path.join(parameters['folder'], parameters['traj_folder_name'], 'data_collection')
        isExist = os.path.exists(full_path_trajectory_dc)
        if isExist:
            pass
        else:
            os.mkdir(full_path_trajectory_dc)
        parameters['trajectory_path_dc'] = full_path_trajectory_dc
        #----------------------------------------------------------------------------------------
        '''
        
        if parameters['input_type'] == 1: # means = parallel
            #------------------------------------------------------------------
            files = os.listdir( parameters['source_folder'])
            pkl_files = []
            for _file in files:
                # Check file extention
                if _file.endswith('.pkl'):
                    pkl_files.append(_file)
            pkl_files.sort()
            
            ##-------------------------------------------------------------------------
            #e_treeview_iter  = getattr(parameters['system'], 'e_treeview_iter', None)
            #e_liststore_iter = getattr(parameters['system'], 'e_liststore_iter', None)
            #parameters['system'].e_treeview_iter   = None
            #parameters['system'].e_liststore_iter  = None            
            ##-------------------------------------------------------------------------
            
            
            #-------------------------------------------------------------------------
            joblist = []
            for pkl in pkl_files:
                #  - - - - adding to the joblist - - - - 
                system  = parameters['system']
                i = int(pkl[5:-4]) 
                joblist.append([i, pkl, parameters, parameters['system']])
            
            print('input_type', parameters['input_type'],  joblist, parameters['NmaxThreads'])
            p = multiprocessing.Pool(processes=parameters['NmaxThreads'])
            #results = p.map(_run_parallel_umbrella_sampling_1D, joblist)
            results = p.map(_run_advanced_parallel_umbrella_sampling_1D, joblist)
            #-------------------------------------------------------------------------           
        
        else:
            #self._run_serial_umbrella_sampling_1D (parameters)
            self._run_advanced_serial_umbrella_sampling_1D (parameters)
    
    def _run_serial_umbrella_sampling_1D (self, parameters):
        """ Function doc """
        
        #-------------------------------------------------------------------------
        #Setting some local vars to ease the notation in the pDynamo methods
        #----------------------------------
        opt_parameters  = parameters['OPT_parm']
        opt_parameters['trajectory_path_opt'] = parameters['trajectory_path_opt']
        
        md_paramters    = parameters['MD_parm']
        
        atom1 = parameters['RC1']['ATOMS'][0]
        atom2 = parameters['RC1']['ATOMS'][1]                   
        #---------------------------------
        restraints = RestraintModel()
        parameters['system'].DefineRestraintModel( restraints )                     
        #----------------------------------------------------------------------------------------
        
        
        full_path_trajectory_eq = parameters['trajectory_path_eq']
        full_path_trajectory_dc = parameters['trajectory_path_dc']
        '''
        #----------------------------------------------------------------------------------------
        full_path_trajectory_eq = os.path.join(parameters['folder'], parameters['traj_folder_name'], 'equilibration')
        os.mkdir(full_path_trajectory_eq)
        
        full_path_trajectory_dc = os.path.join(parameters['folder'], parameters['traj_folder_name'], 'data_collection')
        os.mkdir(full_path_trajectory_dc)
        #----------------------------------------------------------------------------------------
        '''
        
        arq = self.write_header(parameters)
        data = []
        
        for i in range(parameters['RC1']['nsteps']):       
            distance = parameters['RC1']['dminimum'] + ( parameters['RC1']['dincre'] * float(i) )
            
            
            '''----------------------------------------------------------------------------------------------------------------'''
            if parameters['RC1']["rc_type"] == 'simple_distance':
                #---------------------------------------------------------------------------------------------------------
                rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
                restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom1, point2= atom2)
                restraints["RC1"] = restraint            
                #---------------------------------------------------------------------------------------------------------
                
            elif parameters['RC1']["rc_type"] == 'multiple_distance':
                #--------------------------------------------------------------------
                atom3   = parameters['RC1']['ATOMS'][2]
                weight1 = parameters['RC1']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
                weight2 = parameters['RC1']['sigma_pk3pk1'] #self.sigma_a3_a1[0] 
                
                rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
                restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom2, atom1, weight1 ], [ atom2, atom3, weight2 ] ] )
                restraints["RC1"] = restraint            
                #--------------------------------------------------------------------
            
            elif parameters['RC1']["rc_type"] == 'multiple_distance*4atoms':
                #--------------------------------------------------------------------

                atom_3  = parameters['RC1']['ATOMS'][2]
                atom_4  = parameters['RC1']['ATOMS'][3]
                weight1 =  1.0#parameters['RC1']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
                weight2 = -1.0#parameters['RC1']['sigma_pk3pk1'] #self.sigma_a3_a1[0] 
                
                distance_a1_a2 = system.coordinates3.Distance( atom_1, atom_2)
                distance_a2_a3 = system.coordinates3.Distance( atom_3, atom_4)
                
                distance = (weight1 * distance_a1_a2) - (weight2 * distance_a2_a3*-1)
                
                rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
                restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom_2, atom_1, weight1 ], [ atom_3, atom_4, weight2 ] ] )
                restraints["RC1"] = restraint            
                #--------------------------------------------------------------------
            
            
            
            else:
                pass
                
            
            
            '''             G E O M E T R Y   O P T I M I Z A T I O N            '''
            if parameters['OPT_parm'] is not None:
                _us_geo_opt(parameters['system'], opt_parameters)
                Pickle( os.path.join(parameters['trajectory_path_opt'], 
                                                 "frame{}.pkl".format(i) ), 
                                        parameters['system'].coordinates3 )
            
            
            '''                M O L E C U L A R   D Y N A M I C S               '''
            _us_molecular_dynamics (system          = parameters['system'], 
                                    parameters      = md_paramters, 
                                    path_trajectory = full_path_trajectory_eq, 
                                    mode            = 0, 
                                    i = i, 
                                    j = None)
            
            _us_molecular_dynamics (system          = parameters['system'], 
                                    parameters      = md_paramters, 
                                    path_trajectory = full_path_trajectory_dc, 
                                    mode            = 1, 
                                    i = i, 
                                    j = None)
        #---------------------------------------
        parameters['system'].DefineRestraintModel(None)
        
    def _run_advanced_serial_umbrella_sampling_1D (self, parameters):
        """ Function doc """
        
        #----------------------------------------------------------------------------------------
        #    Setting some local vars to ease the notation in the pDynamo methods
        #----------------------------------------------------------------------------------------
        opt_parameters  = parameters['OPT_parm']
        if opt_parameters:
            opt_parameters['trajectory_path_opt'] = parameters['trajectory_path_opt']
        
        md_paramters    = parameters['MD_parm']
        #----------------------------------------------------------------------------------------
        restraints = RestraintModel()
        parameters['system'].DefineRestraintModel( restraints )                     
        #----------------------------------------------------------------------------------------
        
        
        full_path_trajectory_eq = parameters['trajectory_path_eq']
        full_path_trajectory_dc = parameters['trajectory_path_dc']
        
        arq = self.write_header(parameters)
        data = []
        
        #----------------------------------------------------------------- 
        RC = parameters['RC1']['RC']                                       
        RC1 = []                                                           
        distance = 0.0                                                     
        for rc in RC:                                                      
            dist = parameters['system'].coordinates3.Distance(int(rc[1]),  int(rc[3]))   
            dist = dist*float(rc[4]) #weighted distance                    
            distance += dist                                                                                                                          
            RC1.append([ int(rc[1]),  int(rc[3]), float(rc[4]) ])          
        #----------------------------------------------------------------- 
        
        
        for i in range(parameters['RC1']['nsteps']):       
            
            distance = parameters['RC1']['dminimum'] + ( parameters['RC1']['dincre'] * float(i) )
            rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
            restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = RC1)
            restraints["RC1"] = restraint            
            #--------------------------------------------------------------------
               
            
            
            '''             G E O M E T R Y   O P T I M I Z A T I O N            '''
            if parameters['OPT_parm'] is not None:
                _us_geo_opt(parameters['system'], opt_parameters)
                Pickle( os.path.join(parameters['trajectory_path_opt'], 
                                                 "frame{}.pkl".format(i) ), 
                                        parameters['system'].coordinates3 )
            
            
            '''                M O L E C U L A R   D Y N A M I C S               '''
            _us_molecular_dynamics (system          = parameters['system'], 
                                    parameters      = md_paramters, 
                                    path_trajectory = full_path_trajectory_eq, 
                                    mode            = 0, 
                                    i = i, 
                                    j = None)
            
            _us_molecular_dynamics (system          = parameters['system'], 
                                    parameters      = md_paramters, 
                                    path_trajectory = full_path_trajectory_dc, 
                                    mode            = 1, 
                                    i = i, 
                                    j = None)
        #---------------------------------------
        parameters['system'].DefineRestraintModel(None)
        
    def _run_umbrella_sampling_2D (self, parameters, interface = False):
        """ Function doc """
        '''
        #----------------------------------------------------------------------------------------
        full_path_trajectory_eq = os.path.join(parameters['folder'], parameters['traj_folder_name'], 'equilibration')
        isExist = os.path.exists(full_path_trajectory_eq)
        if isExist:
            pass
        else:
            os.mkdir(full_path_trajectory_eq)
        parameters['trajectory_path_eq'] = full_path_trajectory_eq
        #----------------------------------------------------------------------------------------
        
        
        #----------------------------------------------------------------------------------------
        full_path_trajectory_dc = os.path.join(parameters['folder'], parameters['traj_folder_name'], 'data_collection')
        isExist = os.path.exists(full_path_trajectory_dc)
        if isExist:
            pass
        else:
            os.mkdir(full_path_trajectory_dc)
        parameters['trajectory_path_dc'] = full_path_trajectory_dc
        #----------------------------------------------------------------------------------------
        '''
        
        if parameters['input_type'] == 1: # means = parallel
            #------------------------------------------------------------------
            files = os.listdir( parameters['source_folder'])
            pkl_files = []
            for _file in files:
                # Check file extention
                if _file.endswith('.pkl'):
                    pkl_files.append(_file)
            pkl_files.sort()
            
            
            #-------------------------------------------------------------------#
            #                       BUILDING THE JOBLIST                        #
            #-------------------------------------------------------------------#
            joblist = []
            for pkl in pkl_files:
                #  - - - - adding to the joblist - - - - 
                system  = parameters['system']
                
                i_j = pkl[5:-4] 
                i_j = i_j.split('_')
                
                i = int(i_j[0])
                j = int(i_j[1])
                joblist.append([i, j, pkl, parameters, parameters['system']])
            #-------------------------------------------------------------------------
            
            #print([i, j, pkl, parameters, parameters['system']])
            #'''
            #print('input_type', parameters['input_type'],  joblist, parameters['NmaxThreads'])
            p = multiprocessing.Pool(processes=parameters['NmaxThreads'])
            #results = p.map(_run_parallel_umbrella_sampling_2D, joblist)
            results = p.map(_run_advanced_parallel_umbrella_sampling_2D, joblist)


def _run_advanced_parallel_umbrella_sampling_2D (job):
    """ Function doc """
    i          = job[0]
    j          = job[1]
    pkl        = job[2]
    parameters = job[3]
    system     = job[4]


    hamiltonian = get_hamiltonian (system)
    ''' 
    
    '''
    if hamiltonian in ['DFTB QC Model', 'ORCA QC Model','XTB QC Model', 'external']:
        try:
            os.mkdir(system.qcModel.scratch +'/process_'+str(i)+'_'+str(j))
        except:
            pass
     
        try:
            system.qcState.DeterminePaths(system.qcModel.scratch +'/process_'+str(i)+'_'+str(j))
        except:
            pass


    parameters['system'].coordinates3 = ImportCoordinates3(os.path.join(parameters['source_folder'], pkl ))
    
    opt_parameters  = parameters['OPT_parm']
    md_paramters    = parameters['MD_parm']
    #----------------------------------------------------------------------------------------
    full_path_trajectory_eq = parameters['trajectory_path_eq']
    full_path_trajectory_dc = parameters['trajectory_path_dc']
    #----------------------------------------------------------------------------------------

    #----------------------------------------------------------------------------------------    
    restraints = RestraintModel()
    parameters['system'].DefineRestraintModel( restraints )                     
    #----------------------------------------------------------------------------------------
    
    
    #----------------------------------------------------------------------------------------
    #                                       R C 1
    #----------------------------------------------------------------------------------------
    RC = parameters['RC1']['RC']
    RC1 = []
    distance = 0.0
    for rc in RC:
        dist = system.coordinates3.Distance(int(rc[1]),  int(rc[3]))
        
        dist = dist*float(rc[4]) #weighted distance
        distance += dist
            
        RC1.append([ int(rc[1]),  int(rc[3]), float(rc[4]) ])
    #----------------------------------------------------------------------------------------
    '''reaction coordinate 1  '''
    rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
    restraint         = RestraintMultipleDistance.WithOptions( energyModel  = rmodel, 
                                                                  distances = RC1   )
    restraints["RC1"] = restraint            
    #----------------------------------------------------------------------------------------



    #----------------------------------------------------------------------------------------
    #                                       R C 2
    #----------------------------------------------------------------------------------------
    RC  = parameters['RC2']['RC']
    RC2 = []
    distance2 = 0.0
    for rc in RC:
        dist = system.coordinates3.Distance(int(rc[1]),  int(rc[3]))
        
        dist = dist*float(rc[4]) #weighted distance
        distance2 += dist
            
        RC2.append([ int(rc[1]),  int(rc[3]), float(rc[4]) ])
    #----------------------------------------------------------------------------------------

    
    #----------------------------------------------------------------------------------------
    rmodel            = RestraintEnergyModel.Harmonic(distance2, parameters['RC2']['force_constant'])
    restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = RC2 )
    restraints["RC2"] = restraint            
    #----------------------------------------------------------------------------------------
   
   
   
    '''             G E O M E T R Y   O P T I M I Z A T I O N            '''
    if parameters['OPT_parm'] is not None:
        _us_geo_opt(parameters['system'], opt_parameters)
        Pickle( os.path.join(parameters['trajectory_path_opt'], 
                                         "frame{}_{}.pkl".format(i, j) ), 
                                parameters['system'].coordinates3 )


    '''                M O L E C U L A R   D Y N A M I C S               '''
    _us_molecular_dynamics (system          = parameters['system'], 
                            parameters      = md_paramters, 
                            path_trajectory = full_path_trajectory_eq, 
                            mode            = 0, 
                            i = i, 
                            j = j)

    _us_molecular_dynamics (system          = parameters['system'], 
                            parameters      = md_paramters, 
                            path_trajectory = full_path_trajectory_dc, 
                            mode            = 1, 
                            i = i, 
                            j = j)
    #---------------------------------------
    parameters['system'].DefineRestraintModel(None)


def _run_advanced_parallel_umbrella_sampling_1D (job):
    """ Function doc """
    #print('hello')
    i          = job[0]
    pkl        = job[1]
    parameters = job[2]
    system     = job[3]

    
    hamiltonian = get_hamiltonian (system)
    ''' 
    '''
    if hamiltonian in ['DFTB QC Model', 'ORCA QC Model', 'XTB QC Model','external']:
        try:
            os.mkdir(system.qcModel.scratch +'/process_'+str(i))
        except:
            pass
     
        try:
            system.qcState.DeterminePaths(system.qcModel.scratch +'/process_'+str(i))
        except:
            pass

    
    parameters['system'].coordinates3 = ImportCoordinates3(os.path.join(parameters['source_folder'], pkl ))
    #-------------------------------------------------------------------------
    #Setting some local vars to ease the notation in the pDynamo methods
    #----------------------------------
    opt_parameters  = parameters['OPT_parm']
    md_paramters    = parameters['MD_parm']
    
    
    #----------------------------------------------------------------------------------------
    full_path_trajectory_eq = parameters['trajectory_path_eq']
    full_path_trajectory_dc = parameters['trajectory_path_dc']
    #----------------------------------------------------------------------------------------    
    restraints = RestraintModel()
    parameters['system'].DefineRestraintModel( restraints )                     
    #----------------------------------------------------------------------------------------

    #----------------------------------------------------------------------------------------
    RC = parameters['RC1']['RC']
    RC1 = []
    distance = 0.0
    #print('here:', RC)
    
    pprint(parameters)
    
    for rc in parameters['RC1']['RC']:
        #print('\n\n\nlen(rc)',len(rc), rc)
        #print(a1, a2)
        a1 = int(rc[1])
        a2 = int(rc[3])
        dist = system.coordinates3.Distance( a1,  a2 )
        
        dist = dist*float(rc[4]) #weighted distance
        distance += dist
            
        RC1.append([ int(rc[1]),  int(rc[3]), float(rc[4]) ])
    #----------------------------------------------------------------------------------------
    #distance_a1_a2 = system.coordinates3.Distance( atom1, atom2)
    #distance_a2_a3 = system.coordinates3.Distance( atom2, atom3)
    #distance = (weight1 * distance_a1_a2) - (weight2 * distance_a2_a3*-1)
    #----------------------------------------------------------------------------------------


    rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
    restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, 
                                                                 distances = RC1 )
    restraints["RC1"] = restraint            





    '''             G E O M E T R Y   O P T I M I Z A T I O N            '''
    if parameters['OPT_parm'] is not None:
        _us_geo_opt(parameters['system'], opt_parameters)
        
        if parameters['trajectory_path_opt']:
            Pickle( os.path.join(parameters['trajectory_path_opt'], 
                                 "frame{}.pkl".format(i) ), 
                        parameters['system'].coordinates3 ) 


    '''                M O L E C U L A R   D Y N A M I C S               '''
    _us_molecular_dynamics (system          = parameters['system'], 
                            parameters      = md_paramters, 
                            path_trajectory = full_path_trajectory_eq, 
                            mode            = 0, 
                            i = i, 
                            j = None)

    _us_molecular_dynamics (system          = parameters['system'], 
                            parameters      = md_paramters, 
                            path_trajectory = full_path_trajectory_dc, 
                            mode            = 1, 
                            i = i, 
                            j = None)
    #---------------------------------------
    parameters['system'].DefineRestraintModel(None)


def _run_parallel_umbrella_sampling_2D (job):
    """ Function doc """
    i          = job[0]
    j          = job[1]
    pkl        = job[2]
    parameters = job[3]
    system     = job[4]


    hamiltonian = get_hamiltonian (system)
    ''' 
    
    '''
    if hamiltonian in ['DFTB QC Model', 'ORCA QC Model','XTB QC Model', 'external']:
        try:
            os.mkdir(system.qcModel.scratch +'/process_'+str(i)+'_'+str(j))
        except:
            pass
     
        try:
            system.qcState.DeterminePaths(system.qcModel.scratch +'/process_'+str(i)+'_'+str(j))
        except:
            pass


    parameters['system'].coordinates3 = ImportCoordinates3(os.path.join(parameters['source_folder'], pkl ))
    
    opt_parameters  = parameters['OPT_parm']
    md_paramters    = parameters['MD_parm']
    #----------------------------------------------------------------------------------------
    full_path_trajectory_eq = parameters['trajectory_path_eq']
    full_path_trajectory_dc = parameters['trajectory_path_dc']
    #----------------------------------------------------------------------------------------

    
    
    #----------------------------------------------------------------------------------------    
    atom_RC1_1 = parameters['RC1']['ATOMS'][0]
    atom_RC1_2 = parameters['RC1']['ATOMS'][1]                   

    atom_RC2_1 = parameters['RC2']['ATOMS'][0]
    atom_RC2_2 = parameters['RC2']['ATOMS'][1]  
    
    #---------------------------------
    restraints = RestraintModel()
    parameters['system'].DefineRestraintModel( restraints )                     
    #----------------------------------------------------------------------------------------

    '''----------------------------------------------------------------------------------------------------------------'''
    if parameters['RC1']["rc_type"] == 'simple_distance':
        distance1 = system.coordinates3.Distance( atom_RC1_1, atom_RC1_2)
        #---------------------------------------------------------------------------------------------------------
        rmodel            = RestraintEnergyModel.Harmonic(distance1, parameters['RC1']['force_constant'])
        restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom_RC1_1, point2= atom_RC1_2)
        restraints["RC1"] = restraint            
        #---------------------------------------------------------------------------------------------------------
    elif parameters['RC1']["rc_type"] == 'multiple_distance':
        #--------------------------------------------------------------------

        atom_RC1_3 = parameters['RC1']['ATOMS'][2]
        weight1 = parameters['RC1']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
        weight2 = parameters['RC1']['sigma_pk3pk1'] #self.sigma_a3_a1[0] 
        
        distance_a1_a2 = system.coordinates3.Distance( atom_RC1_1, atom_RC1_2)
        distance_a2_a3 = system.coordinates3.Distance( atom_RC1_2, atom_RC1_3)
        
        distance = (weight1 * distance_a1_a2) - (weight2 * distance_a2_a3*-1)
        
        rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
        restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ 
                                                                                                      [ atom_RC1_2, atom_RC1_1, weight1 ], 
                                                                                                      [ atom_RC1_2, atom_RC1_3, weight2 ] 
                                                                                                    ] )
        restraints["RC1"] = restraint            
        #--------------------------------------------------------------------
    elif parameters['RC1']["rc_type"] == 'multiple_distance*4atoms':
        #--------------------------------------------------------------------
        atom_RC1_3   = parameters['RC1']['ATOMS'][2]
        atom_RC1_4   = parameters['RC1']['ATOMS'][3]
        weight1 =  1.0#parameters['RC1']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
        weight2 = -1.0#parameters['RC1']['sigma_pk3pk1'] #self.sigma_a3_a1[0] 
        
        distance_a1_a2 = system.coordinates3.Distance( atom_RC1_1, atom_RC1_2)
        distance_a2_a3 = system.coordinates3.Distance( atom_RC1_3, atom_RC1_4)
        
        distance = (weight1 * distance_a1_a2) - (weight2 * distance_a2_a3*-1)
        
        rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
        restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom_RC1_2, atom_RC1_1, weight1 ], [ atom_RC1_3, atom_RC1_4, weight2 ] ] )
        restraints["RC1"] = restraint            
        #--------------------------------------------------------------------
    else:
        pass


    '''reaction coordinate 2 - ONLY at 0 first'''
    '''----------------------------------------------------------------------------------------------------------------'''
    #distance2 = parameters['RC2']['dminimum'] + ( parameters['RC2']['dincre'] * float(j) )
    
    if parameters['RC2']["rc_type"] == 'simple_distance':
        distance2 = system.coordinates3.Distance( atom_RC2_1, atom_RC2_2)
        #---------------------------------------------------------------------------------------------------------
        rmodel            = RestraintEnergyModel.Harmonic(distance2, parameters['RC2']['force_constant'])
        restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom_RC2_1, point2= atom_RC2_2)
        restraints["RC2"] = restraint            
        #---------------------------------------------------------------------------------------------------------
    
    elif parameters['RC2']["rc_type"] == 'multiple_distance':
        #--------------------------------------------------------------------
        atom_RC2_3 = parameters['RC2']['ATOMS'][2]
        weight1 = parameters['RC2']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
        weight2 = parameters['RC2']['sigma_pk3pk1'] #self.sigma_a3_a1[0] 
        
        distance_a1_a2 = system.coordinates3.Distance( atom_RC2_1, atom_RC2_2)
        distance_a2_a3 = system.coordinates3.Distance( atom_RC2_2, atom_RC2_3)
        
        distance2 = (weight1 * distance_a1_a2) - (weight2 * distance_a2_a3*-1)
        
        
        rmodel            = RestraintEnergyModel.Harmonic(distance2, parameters['RC2']['force_constant'])
        restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ 
                                                                                                      [ atom_RC2_2, atom_RC2_1, weight1 ], 
                                                                                                      [ atom_RC2_2, atom_RC2_3, weight2 ] 
                                                                                                      ] )
        restraints["RC2"] = restraint            
        #--------------------------------------------------------------------
   
    elif parameters['RC2']["rc_type"] == 'multiple_distance*4atoms':
        #--------------------------------------------------------------------

        atom_RC2_3   = parameters['RC2']['ATOMS'][2]
        atom_RC2_4   = parameters['RC2']['ATOMS'][3]
        weight1 =  1.0#parameters['RC1']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
        weight2 = -1.0#parameters['RC1']['sigma_pk3pk1'] #self.sigma_a3_a1[0] 
        
        distance_a1_a2 = system.coordinates3.Distance( atom_RC2_1, atom_RC2_2)
        distance_a2_a3 = system.coordinates3.Distance( atom_RC2_3, atom_RC2_4)
        
        distance = (weight1 * distance_a1_a2) - (weight2 * distance_a2_a3*-1)
        
        rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC2']['force_constant'])
        restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom_RC2_1, atom_RC2_2, weight1 ], [ atom_RC2_3, atom_RC2_4, weight2 ] ] )
        restraints["RC2"] = restraint            
        #--------------------------------------------------------------------
   
   
   
   
   
   
   
    '''             G E O M E T R Y   O P T I M I Z A T I O N            '''
    if parameters['OPT_parm'] is not None:
        _us_geo_opt(parameters['system'], opt_parameters)
        Pickle( os.path.join(parameters['trajectory_path_opt'], 
                                         "frame{}_{}.pkl".format(i, j) ), 
                                parameters['system'].coordinates3 )


    '''                M O L E C U L A R   D Y N A M I C S               '''
    _us_molecular_dynamics (system          = parameters['system'], 
                            parameters      = md_paramters, 
                            path_trajectory = full_path_trajectory_eq, 
                            mode            = 0, 
                            i = i, 
                            j = j)

    _us_molecular_dynamics (system          = parameters['system'], 
                            parameters      = md_paramters, 
                            path_trajectory = full_path_trajectory_dc, 
                            mode            = 1, 
                            i = i, 
                            j = j)
    #---------------------------------------
    parameters['system'].DefineRestraintModel(None)


def _run_parallel_umbrella_sampling_1D (job):
    """ Function doc """
    #print('hello')
    i          = job[0]
    pkl        = job[1]
    parameters = job[2]
    system     = job[3]

    
    hamiltonian = get_hamiltonian (system)
    ''' 
    '''
    if hamiltonian in ['DFTB QC Model', 'ORCA QC Model', 'XTB QC Model','external']:
        try:
            os.mkdir(system.qcModel.scratch +'/process_'+str(i))
        except:
            pass
     
        try:
            system.qcState.DeterminePaths(system.qcModel.scratch +'/process_'+str(i))
        except:
            pass

    
    parameters['system'].coordinates3 = ImportCoordinates3(os.path.join(parameters['source_folder'], pkl ))
    #-------------------------------------------------------------------------
    #Setting some local vars to ease the notation in the pDynamo methods
    #----------------------------------
    opt_parameters  = parameters['OPT_parm']
    md_paramters    = parameters['MD_parm']
    
    
    #----------------------------------------------------------------------------------------
    full_path_trajectory_eq = parameters['trajectory_path_eq']
    full_path_trajectory_dc = parameters['trajectory_path_dc']
    #----------------------------------------------------------------------------------------    
    
    
    #----------------------------------------------------------------------------------------    
    atom1 = parameters['RC1']['ATOMS'][0]
    atom2 = parameters['RC1']['ATOMS'][1]                   
    #---------------------------------
    restraints = RestraintModel()
    parameters['system'].DefineRestraintModel( restraints )                     
    #----------------------------------------------------------------------------------------


    '''----------------------------------------------------------------------------------------------------------------'''
    if parameters['RC1']["rc_type"] == 'simple_distance':
        distance = system.coordinates3.Distance( atom1, atom2)
        #---------------------------------------------------------------------------------------------------------
        rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
        restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom1, point2= atom2)
        restraints["RC1"] = restraint            
        #---------------------------------------------------------------------------------------------------------
    elif parameters['RC1']["rc_type"] == 'multiple_distance':
        #--------------------------------------------------------------------

        atom3   = parameters['RC1']['ATOMS'][2]
        weight1 = parameters['RC1']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
        weight2 = parameters['RC1']['sigma_pk3pk1'] #self.sigma_a3_a1[0] 
        
        distance_a1_a2 = system.coordinates3.Distance( atom1, atom2)
        distance_a2_a3 = system.coordinates3.Distance( atom2, atom3)
        
        distance = (weight1 * distance_a1_a2) - (weight2 * distance_a2_a3*-1)
        
        rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
        restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom2, atom1, weight1 ], [ atom2, atom3, weight2 ] ] )
        restraints["RC1"] = restraint            
        #--------------------------------------------------------------------
    
    elif parameters['RC1']["rc_type"] == 'multiple_distance*4atoms':
        #--------------------------------------------------------------------

        atom3   = parameters['RC1']['ATOMS'][2]
        atom4   = parameters['RC1']['ATOMS'][3]
        weight1 =  1.0#parameters['RC1']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
        weight2 = -1.0#parameters['RC1']['sigma_pk3pk1'] #self.sigma_a3_a1[0] 
        
        distance_a1_a2 = system.coordinates3.Distance( atom1, atom2)
        distance_a2_a3 = system.coordinates3.Distance( atom3, atom4)
        
        distance = (weight1 * distance_a1_a2) - (weight2 * distance_a2_a3*-1)
        
        rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
        restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom2, atom1, weight1 ], [ atom3, atom4, weight2 ] ] )
        restraints["RC1"] = restraint            
        #--------------------------------------------------------------------
    else:
        pass




    '''             G E O M E T R Y   O P T I M I Z A T I O N            '''
    if parameters['OPT_parm'] is not None:
        _us_geo_opt(parameters['system'], opt_parameters)
        
        if parameters['trajectory_path_opt']:
            Pickle( os.path.join(parameters['trajectory_path_opt'], 
                                 "frame{}.pkl".format(i) ), 
                        parameters['system'].coordinates3 ) 
        #opt_parameters['trajectory_path_opt']


    '''                M O L E C U L A R   D Y N A M I C S               '''
    _us_molecular_dynamics (system          = parameters['system'], 
                            parameters      = md_paramters, 
                            path_trajectory = full_path_trajectory_eq, 
                            mode            = 0, 
                            i = i, 
                            j = None)

    _us_molecular_dynamics (system          = parameters['system'], 
                            parameters      = md_paramters, 
                            path_trajectory = full_path_trajectory_dc, 
                            mode            = 1, 
                            i = i, 
                            j = None)
    #---------------------------------------
    parameters['system'].DefineRestraintModel(None)


def _us_geo_opt (system, parameters):
    """ Function doc """
    parameters['logFrequency'] = 50
    
    if parameters['optimizer'] == 'ConjugatedGradient':
        #-------------------------------------------------------------------------------------------------------------
        ConjugateGradientMinimize_SystemGeometry(system                     ,                
                                                 logFrequency           = parameters['logFrequency'],
                                                 maximumIterations      = parameters['maximumIterations'],
                                                 rmsGradientTolerance   = parameters['rmsGradientTolerance'])
        #-------------------------------------------------------------------------------------------------------------
    
    elif parameters['optimizer'] == 'SteepestDescent':
        SteepestDescentMinimize_SystemGeometry( system                                            ,               
                                            logFrequency            = parameters['logFrequency'] ,
                                            maximumIterations       = parameters['maximumIterations'] ,
                                            rmsGradientTolerance    = parameters['rmsGradientTolerance']   )
    
    elif parameters['optimizer'] == 'LFBGS':
        LBFGSMinimize_SystemGeometry(parameters['system']                               ,                
                                     logFrequency         = parameters['logFrequency'] ,
                                     maximumIterations    = parameters['maximumIterations'] ,
                                     rmsGradientTolerance = parameters['rmsGradientTolerance']   )
    else:
        pass


def _us_molecular_dynamics (system, parameters, path_trajectory, mode = 0, i = None, j = None):
    """ 
    
    mode =  0 # equilibration
    mode =  1 # data collection
    
    i = int  - 1a  RC index
    j = int  - 2sd RC index  
    Function doc """
    

    #-----------------------------------------------------------------------------------------------------
    # . Equilibration.
    if mode == 0:        
        if j == None:
            full_path_trajectory = os.path.join(path_trajectory, 'window{:d}.ptGeo'.format ( i ))
        else:
            full_path_trajectory = os.path.join(path_trajectory, 'window{:d}_{:d}.ptGeo'.format ( i, j ))
        parameters['steps']  = parameters['steps_eq'] 
        
        
        trajectory = ExportTrajectory(full_path_trajectory,  system , log=None )
        trajectories = [(trajectory, parameters['trajectory_frequency'])]
    #-----------------------------------------------------------------------------------------------------
        
    
    
    
    #-----------------------------------------------------------------------------------------------------
    # . Data Collection
    elif mode == 1:
        
        if j == None:
            full_path_trajectory = os.path.join(path_trajectory, 'window{:d}.ptRes'.format ( i ))
        else:
            full_path_trajectory = os.path.join(path_trajectory, 'window{:d}_{:d}.ptRes'.format ( i, j ))
        parameters['steps']  = parameters['steps_dc']
        
        path_trajectory2 =  path_trajectory + '_traj'
        if j == None:
            full_path_trajectory2 = os.path.join(path_trajectory2, 'window{:d}.ptGeo'.format ( i ))
        else:
            full_path_trajectory2 = os.path.join(path_trajectory2, 'window{:d}_{:d}.ptGeo'.format ( i, j ))
        #parameters['steps']  = parameters['steps_dc']
        
        
        trajectory1 = ExportTrajectory(full_path_trajectory,   system , log=None )
        trajectory2 = ExportTrajectory(full_path_trajectory2,  system , log=None )
        trajectories = [(trajectory1, parameters['trajectory_frequency_dc_ptRes']),
                        (trajectory2, parameters['trajectory_frequency_dc_ptGeo'])]
    
    
    #-----------------------------------------------------------------------------------------------------
    else:
        pass
    
    
    #trajectoy instances
    #trajectory = ExportTrajectory(full_path_trajectory,  system , log=None )
    
    
    
    
        
    if parameters['integrator'] == 'Verlet':
        #_us_velocity_verlet_dynamics(system, trajectory, parameters)
        _us_velocity_verlet_dynamics(system, trajectories, parameters)
    
    elif parameters['integrator'] == 'LeapFrog':
        #_us_leap_frog_dynamics (system, trajectory, parameters)
        _us_leap_frog_dynamics (system, trajectories, parameters)
    
    elif parameters['integrator'] == 'Langevin':
        #_us_langevin_dynamics (system, trajectory, parameters)
        _us_langevin_dynamics (system, trajectories, parameters)
    
    else:
        pass


def _us_leap_frog_dynamics (system, trajectories, parameters ):
    """ Function doc """
    if parameters['temperatureControl']:
        LeapFrogDynamics_SystemGeometry ( parameters['system']                       ,
                                          logFrequency           = parameters['logFrequency'] ,
                                          #normalDeviateGenerator = normalDeviateGenerator ,
                                          pressure               = parameters['pressure'],
                                          pressureControl        = parameters['pressureControl'],
                                          pressureCoupling       = parameters['pressureCoupling'],
                                          steps                  = parameters['steps'],
                                          temperature            = parameters['temperatureStart'],
                                          temperatureControl     = parameters['temperatureControl'],
                                          temperatureCoupling    = parameters['temperatureCoupling'],
                                          timeStep               = parameters['timeStep'],
                                          trajectories           = trajectories, #[(trajectory, parameters['trajectory_frequency'])],
                                          )
                    

    else:
        LeapFrogDynamics_SystemGeometry ( system                                               ,
                                          logFrequency           =  parameters['logFrequency'] ,
                                          #normalDeviateGenerator = normalDeviateGenerator ,
                                          steps                  =  parameters['steps'] ,
                                          temperature            =  parameters['temperatureStart'] ,
                                          temperatureControl     =  parameters['temperatureControl'] ,
                                          temperatureCoupling    =  parameters['temperatureCoupling'] ,
                                          timeStep               =  parameters['timeStep'] ,
                                          trajectories           = trajectories,#[(trajectory,  parameters['trajectory_frequency'])]
                                          
                                          )


def _us_velocity_verlet_dynamics (system, trajectories, parameters):
    """ Function doc """
    # . Define a normal deviate generator in a given state.
    #normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( parameters['seed'] ) )        

    VelocityVerletDynamics_SystemGeometry ( system                                                     ,
                                            logFrequency              = parameters['logFrequency']                  ,
                                            #normalDeviateGenerator    = normalDeviateGenerator                      ,
                                            steps                     = parameters['steps']                         ,
                                            timeStep                  = parameters['timeStep']                      ,
                                            temperatureScaleFrequency = parameters['temperatureScaleFrequency']     ,
                                            temperatureScaleOption    = "constant"                                  ,
                                            temperatureStart          = parameters['temperatureStart']              ,
                                            trajectories              = trajectories,#[(trajectory, parameters['trajectory_frequency'])],
                                            #log                       = self.logFile2
                                            )


def _us_langevin_dynamics (system, trajectory, parameters):
    
    """ Function doc """
    # . Define a normal deviate generator in a given state.
    #normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( parameters['seed'] ) )        
    LangevinDynamics_SystemGeometry ( system                                                             ,
                                      collisionFrequency     = parameters['collisionFrequency']                       ,
                                      logFrequency           = parameters['logFrequency']                             ,
                                      #normalDeviateGenerator = normalDeviateGenerator                                 ,
                                      steps                  = parameters['steps']                                    ,
                                      temperature            = parameters['temperatureStart']                         ,
                                      timeStep               = parameters['timeStep']                                 ,
                                      trajectories           = trajectories, #[(trajectory, parameters['trajectory_frequency'])] #, (traj2, int)],
                                      #log                    = self.logFile2
                                      )
