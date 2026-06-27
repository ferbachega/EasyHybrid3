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
from pdynamo.p_methods._common import backup_orca_files, write_header, get_hamiltonian, plot_data

class AdvancedRelaxedSurfaceScan:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass
    
    def run (self, parameters, interface = False):
        """ Function doc """
        full_path_trajectory =os.path.join(parameters['folder'], 
                              parameters['traj_folder_name']+".ptGeo")
        os.mkdir(
                 full_path_trajectory
                 )
        
        # - - - - - - - - - - - - - Checking trajectory - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        self.logFile2 = TextLogFileWriter.WithOptions ( path = os.path.join(full_path_trajectory, 'output.log') )
        parameters['system'].Summary(log = self.logFile2)
        self.logFile2.Header ( )
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        #if parameters['_is_ts_centered']: 
        #    advanced = True
        
        if parameters['RC2'] is not None:
            #pass
            #print('\n\n\nself._run_scan_2D(parameters = parameters, interface = False)')
            #
            if parameters['_is_ts_centered']:
                pass
            #    parameters['RC1']['up_nsteps'  ] = parameters['RC1']['nsteps']
            #    parameters['RC1']['down_nsteps'] = parameters['RC1']['nsteps_back']
            #    
            #    parameters['RC2']['left_nsteps' ] = parameters['RC2']['nsteps']
            #    parameters['RC2']['right_nsteps'] = parameters['RC2']['nsteps_back']
            #    
            #    
            #    _run_advanded_scan_2D (parameters = parameters, interface = False)
            #
            else:
                self._run_scan_2D(parameters = parameters, interface = False)
        else:
            self._run_scan_1D(parameters)
        
        
        
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        self.logFile2.Footer ( )
        self.logFile2.Close()
        self.logFile2 = None
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def write_header (self, parameters, logfile = 'output.log'):
        """ Function doc """
        arq = open(os.path.join( parameters['folder'], parameters['traj_folder_name']+".ptGeo", logfile), "a")
        text = ""
        
        
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


        text = text + "\n--------------------------------- Coordinate 1 ---------------------------------"
        for rc1 in parameters["RC1"]['RC']:
            text = text + "\nATOM                   =%15i  ATOM NAME              =%15s"     % ( int(rc1[1]), str(rc1[0]) )
            text = text + "\nATOM                   =%15i  ATOM NAME              =%15s"     % ( int(rc1[3]), str(rc1[2]) )
            text = text + "\nSigma                  =%15.5f"                                 % ( float(rc1[4]) )
        text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['RC1']['nsteps']  , parameters['RC1']['force_constant'] )
        text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC1']['dminimum'], parameters['maximumIterations']         )
        text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC1']['dincre']  , parameters['rmsGradientTolerance']           )
        text = text + "\n--------------------------------------------------------------------------------"
        
        if parameters["RC2"] is not None :
            text = text + "\n--------------------------------- Coordinate 2 ---------------------------------"
            for rc2 in parameters["RC2"]['RC']:#int(rc[1]),  int(rc[3]), float(rc[4])
                text = text + "\nATOM                   =%15i  ATOM NAME              =%15s"     % ( int(rc2[1]), str(rc2[0]) )
                text = text + "\nATOM                   =%15i  ATOM NAME              =%15s"     % ( int(rc2[3]), str(rc2[2]) )
                text = text + "\nSigma                  =%15.5f"                                 % ( float(rc2[4]) )
            text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['RC1']['nsteps']  , parameters['RC1']['force_constant'] )
            text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC1']['dminimum'], parameters['maximumIterations']         )
            text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC1']['dincre']  , parameters['rmsGradientTolerance']           )
            text = text + "\n--------------------------------------------------------------------------------"


        arq.write(text)

        return arq

    def _run_scan_1D (self, parameters):
        """ Function doc """
        #-------------------------------------------------------------------------
        #Setting some local vars to ease the notation in the pDynamo methods
        #----------------------------------              
        #---------------------------------
        restraints = RestraintModel()
        parameters['system'].DefineRestraintModel( restraints )                     
        #----------------------------------------------------------------------------------------
        RC = parameters['RC1']['RC']
        RC1 = []
        #dminimum = 0.0
        for rc in RC:
            dist = float(rc[5])*float(rc[4])
            #dminimum += dist
            RC1.append([ int(rc[1]),  int(rc[3]), float(rc[4])])
        
        #parameters['RC1']['dminimum'] = dminimum
        
        arq  = self.write_header(parameters)
        data = []
        
        for i in range(parameters['RC1']['nsteps']):       
            distance = parameters['RC1']['dminimum'] + ( parameters['RC1']['dincre'] * float(i) )

            '''----------------------------------------------------------------------------------------------------------------'''
            #if parameters['RC1']["rc_type"] == 'multiple_distance':
                #--------------------------------------------------------------------
            rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
            #restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = [ [ atom2, atom1, weight1 ], [ atom2, atom3, weight2 ] ] )
            restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = RC1 )
            restraints["RC1"] = restraint            
            
            #--------------------------------------------------------------------
        
            #-------------------------------------------------------------------------------------------------------------
            ConjugateGradientMinimize_SystemGeometry(parameters['system']                                ,                
                                                     logFrequency           = parameters['logFrequency'],
                                                     maximumIterations      = parameters['maximumIterations'],
                                                     rmsGradientTolerance   = parameters['rmsGradientTolerance'])
            #-------------------------------------------------------------------------------------------------------------
        
            #distance1 = parameters['system'].coordinates3.Distance( atom1 , atom2  )
            #distance2 = parameters['system'].coordinates3.Distance( atom2 , atom3  )
            energy   = parameters['system'].Energy(log=None)
            print(distance, energy)
            #data.append([i, distance1, distance2, energy])             
            data.append([i, distance, energy])             
            
            #text = "\nDATA %9i      %13.12f        %13.12f"% (int(i), float(distance), float(energy))
            text = "\nDATA {:>9d}   {:>18.12f}   {:>18.12f}".format(int(i), float(distance), float(energy))
            arq.write(text)
            
            
            Pickle( os.path.join( parameters['folder'], 
                                  parameters['traj_folder_name']+".ptGeo", 
                                  "frame{}.pkl".format(i) ), 
                    parameters['system'].coordinates3 ) 
            
            backup_orca_files(system        = parameters['system'], 
                              output_folder = os.path.join(parameters['folder'],parameters['traj_folder_name']+".ptGeo") , 
                              output_name   = "frame{}".format(i))
                

        parameters['system'].DefineRestraintModel(None)
        pprint(data)
        
    def _run_scan_2D (self, parameters = None, interface = False):
        """ Function doc """
        #-------------------------------------------------------------------------
        #Setting some local vars to ease the notation in the pDynamo methods
        #----------------------------------
        restraints = RestraintModel()
        parameters['system'].DefineRestraintModel( restraints )                     
        #----------------------------------------------------------------------------------------
 
        #----------------------------------------------------------------------------------------
        RC = parameters['RC1']['RC']
        RC1 = []
        #dminimum = 0.0
        for rc in RC:
            dist = float(rc[5])*float(rc[4])
            #dminimum += dist
            RC1.append([ int(rc[1]),  int(rc[3]), float(rc[4])])
        #parameters['RC1']['dminimum'] = dminimum
        #----------------------------------------------------------------------------------------
        
        
        #----------------------------------------------------------------------------------------
        RC = parameters['RC2']['RC']
        RC2 = []
        #dminimum2 = 0.0
        for rc in RC:
            dist = float(rc[5])*float(rc[4])
            #dminimum2 += dist
            RC2.append([ int(rc[1]),  int(rc[3]), float(rc[4])])
        #parameters['RC2']['dminimum'] = dminimum2
        #----------------------------------------------------------------------------------------
        
        
        #----------------------------------------------------------------------------------------
        arq = self.write_header(parameters)
        data = {} 
        #----------------------------------------------------------------------------------------
        '''We have to first run a sequence (first line) of optimizations sequentially. 
        This will generate frames (0,0) to (n, 0) serially, where n is the number of 
        steps from coordinate x)
        '''
        #try:
        joblist = []
        j = 0
        for i in range(parameters['RC1']['nsteps']):       
            distance = parameters['RC1']['dminimum'] + ( parameters['RC1']['dincre'] * float(i) )
            
            '''reaction coordinate 1 - starts at 0 and goes to nstpes'''
            '''----------------------------------------------------------------------------------------------------------------'''
            rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
            restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = RC1 )
            restraints["RC1"] = restraint 


            '''reaction coordinate 2 - ONLY at 0 first'''
            '''----------------------------------------------------------------------------------------------------------------'''
            distance2 = parameters['RC2']['dminimum'] + ( parameters['RC2']['dincre'] * float(j) )
            #--------------------------------------------------------------------
            rmodel            = RestraintEnergyModel.Harmonic(distance2, parameters['RC2']['force_constant'])
            restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = RC2 )
            restraints["RC2"] = restraint            
            #--------------------------------------------------------------------
            #ConjugateGradientMinimize_SystemGeometry(parameters['system']                              ,                
            #                                         logFrequency           = parameters['logFrequency'],
            #                                         maximumIterations      = parameters['maximumIterations'],
            #                                         rmsGradientTolerance   = parameters['rmsGradientTolerance'], 
            #                                         log= None) 
            #energy = parameters['system'].Energy(log=None)
            try:
                #-------------------------------------------------------------------------------------------------------------
                ConjugateGradientMinimize_SystemGeometry(parameters['system']                                ,                
                                                         logFrequency           = parameters['logFrequency'],
                                                         maximumIterations      = parameters['maximumIterations'],
                                                         rmsGradientTolerance   = parameters['rmsGradientTolerance']  ,
                                                         log                    = None)
            
                energy = parameters['system'].Energy(log=None)
                opt_convergency = ''
                #-------------------------------------------------------------------------------------------------------------
            except:
                energy   = 0
                opt_convergency = 'Ops! Error in geometry optimization.'

            #--------------------------------------------------------------------------------------
            #                     Calculating the reaction coordinate 1
            #--------------------------------------------------------------------------------------
            '''
            distance1 = parameters['system'].coordinates3.Distance( atom_RC1_1 , atom_RC1_2  )
            
            if parameters['RC1']["rc_type"] == 'multiple_distance':
                distance2 = parameters['system'].coordinates3.Distance( atom_RC1_2,  atom_RC1_3  )

            elif parameters['RC1']["rc_type"] == 'multiple_distance*4atoms':
                distance2 = parameters['system'].coordinates3.Distance( atom_RC1_3,  atom_RC1_4  )
            else:
                distance2 = 0
            
            RC1_d1_minus_d2 = distance1 - distance2
            #--------------------------------------------------------------------------------------
            
            
            
            #--------------------------------------------------------------------------------------
            #                     Calculating the reaction coordinate 2
            #--------------------------------------------------------------------------------------
            distance1 = parameters['system'].coordinates3.Distance( atom_RC2_1 , atom_RC2_2  )
            
            if parameters['RC2']["rc_type"] == 'multiple_distance':
                distance2 = parameters['system'].coordinates3.Distance( atom_RC2_2,  atom_RC2_3  )
            
            elif parameters['RC2']["rc_type"] == 'multiple_distance*4atoms':
                distance2 = parameters['system'].coordinates3.Distance( atom_RC2_3,  atom_RC2_4  )
            
            else:
                distance2 = 0
            
            RC2_d1_minus_d2 = distance1 - distance2
            #--------------------------------------------------------------------------------------
            #'''
            
            
            #--------------------------------------------------------------------------------------
            #                          Energy and data dictionary
            #--------------------------------------------------------------------------------------
            data[(i,j)] = [distance, distance2, energy]    
            #--------------------------------------------------------------------------------------
            
            
            
            #--------------------------------------------------------------------------------------
            #                             Exporting Coordinates
            #--------------------------------------------------------------------------------------
            pkl = os.path.join( parameters['folder'], 
                                  parameters['traj_folder_name']+".ptGeo", 
                                  "frame{}_{}.pkl".format(i,j) )
                                   
            Pickle( pkl, parameters['system'].coordinates3 )
            
            backup_orca_files(system        = parameters['system'], 
                              output_folder = os.path.join(parameters['folder'],parameters['traj_folder_name']+".ptGeo") , 
                              output_name   = "frame{}_{}".format(i,j))
            
            #--------------------------------------------------------------------------------------

            new_parameters = {
                             
                             'RC1'          : parameters['RC1'],
                             'RC2'          : parameters['RC2'],
                             'dminimum_RC1' : parameters['RC1']['dminimum'],
                             'dminimum_RC2' : parameters['RC2']['dminimum'],

                             'force_constant_1': parameters['RC1']['force_constant'],
                             'nsteps_RC1'      : parameters['RC1']['nsteps'],
                             'dincre_RC1'      : parameters['RC1']['dincre'],

                             'force_constant_2': parameters['RC2']['force_constant'],
                             'nsteps_RC2'      : parameters['RC2']['nsteps'],
                             'dincre_RC2'      : parameters['RC2']['dincre'],
                             
                             'logFrequency'        : parameters['logFrequency'],
                             'maximumIterations'   : parameters['maximumIterations'],
                             'rmsGradientTolerance': parameters['rmsGradientTolerance'],

                             'folder'          : parameters['folder']    , 
                             'traj_folder_name': parameters['traj_folder_name']
                             }
            
            #print(i, j, 'Energy:', energy, opt_convergency)
            print(i, j, 'Energy:', energy)
            #  - - - - adding to the joblist - - - - 
            joblist.append([i, pkl, new_parameters, parameters['system']])
            #--------------------------------------------------------------------------------------


        #print(joblist)
        '''
        results = []
        for job in joblist:
            result = _run_advanced_second_coordinate_in_parallel(job)
            results.append(result)
        '''
        p = multiprocessing.Pool(processes=parameters['NmaxThreads'])
        #p = multiprocessing.Pool(processes=8)
        results = p.map(_run_advanced_second_coordinate_in_parallel, joblist)
        #results = p.map(_run_advanced_second_coordinate_in_parallel, joblist)
        
        

        
        for partial_data_dict in  results:
            for key, partial_data in partial_data_dict.items():
                data[key] = partial_data

        #--------------------------------------------------------------------------------------
        #                             Writing the log information 
        #--------------------------------------------------------------------------------------
        for i in range(parameters['RC1']['nsteps']):
            for j in range(parameters['RC2']['nsteps']):
        
                #text = "\nDATA  %4i  %4i     %13.12f       %13.12f       %13.12f"% (int(i), int(j),  float(data[(i,j)][0]), float(data[(i,j)][1]), float(data[(i,j)][2]))
                text = "\nDATA  {:>4d} {:>4d}   {:>18.12f}   {:>18.12f}   {:>18.12f}".format(int(i), int(j),  float(data[(i,j)][0]), float(data[(i,j)][1]), float(data[(i,j)][2]))
                arq.write(text)


class RelaxedSurfaceScan:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass

    def run (self, parameters, interface = False):
        """ Function doc """
        full_path_trajectory =os.path.join(parameters['folder'], 
                              parameters['traj_folder_name']+".ptGeo")
        os.mkdir(
                 full_path_trajectory
                 )
        
        # - - - - - - - - - - - - - Checking trajectory - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        self.logFile2 = TextLogFileWriter.WithOptions ( path = os.path.join(full_path_trajectory, 'output.log') )
        parameters['system'].Summary(log = self.logFile2)
        self.logFile2.Header ( )
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        #if parameters['_is_ts_centered']: 
        #    advanced = True
        
        if parameters['RC2'] is not None:
            print('\n\n\nself._run_scan_2D(parameters = parameters, interface = False)')
            
            if parameters['_is_ts_centered']:
                parameters['RC1']['up_nsteps'  ] = parameters['RC1']['nsteps']
                parameters['RC1']['down_nsteps'] = parameters['RC1']['nsteps_back']
                
                parameters['RC2']['left_nsteps' ] = parameters['RC2']['nsteps']
                parameters['RC2']['right_nsteps'] = parameters['RC2']['nsteps_back']
                
                
                _run_advanded_scan_2D (parameters = parameters, interface = False)
            
            else:
                self._run_scan_2D(parameters = parameters, interface = False)
        else:
            self._run_scan_1D(parameters)
        
        
        
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        self.logFile2.Footer ( )
        self.logFile2.Close()
        self.logFile2 = None
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
            
        
    def write_header (self, parameters, logfile = 'output.log'):
        """ Function doc """
        
        arq = open(os.path.join( parameters['folder'], parameters['traj_folder_name']+".ptGeo", logfile), "a")
        text = ""

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




        '''This part writes the parameters used in the first reaction coordinate'''
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------
        if parameters['RC1']["rc_type"] == 'simple_distance':
            #text = text + "\n"
            #text = text + "\n----------------------- Coordinate 1 - Simple-Distance -------------------------"
            #text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['ATOMS_RC1'][0], parameters['ATOMS_RC1_NAMES'][0] )
            #text = text + "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['ATOMS_RC1'][1], parameters['ATOMS_RC1_NAMES'][1] )
            #text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['nsteps_RC1']  , parameters['force_constant_1']   )
            #text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['dminimum_RC1'], parameters['maximumIterations']      )
            #text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['dincre_RC1']  , parameters['rmsGradientTolerance']        )
            #text = text + "\n--------------------------------------------------------------------------------"
            text = text + "\n"
            text = text + "\n----------------------- Coordinate 1 - Simple-Distance -------------------------"
            text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC1']['ATOMS'][0], parameters['RC1']['ATOM_NAMES'][0] )
            text = text + "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['RC1']['ATOMS'][1], parameters['RC1']['ATOM_NAMES'][1] )
            text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['RC1']['nsteps']  , parameters['RC1']['force_constant'] )
            text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC1']['dminimum'], parameters['maximumIterations']         )
            text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC1']['dincre']  , parameters['rmsGradientTolerance']           )
            text = text + "\n--------------------------------------------------------------------------------"

        
        elif parameters['RC1']["rc_type"] == 'multiple_distance*4atoms':
            text = text + "\n"
            text = text + "\n---------------------- Coordinate 1 - multiple-Distance ------------------------"	
            text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC1']['ATOMS'][0]    , parameters['RC1']['ATOM_NAMES'][0] )
            text = text + "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['RC1']['ATOMS'][1]    , parameters['RC1']['ATOM_NAMES'][1] )
            text = text + "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC1']['ATOMS'][2]    , parameters['RC1']['ATOM_NAMES'][2] )
            text = text + "\nATOM4                  =%15i  ATOM NAME4             =%15s"     % (parameters['RC1']['ATOMS'][3]    , parameters['RC1']['ATOM_NAMES'][3] )
            text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['RC1']['nsteps']      , parameters['RC1']['force_constant'] ) 
            text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC1']['dminimum']    , parameters['maximumIterations']         )
            text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC1']['dincre']      , parameters['rmsGradientTolerance']           )
            #text = text + "\nSigma atom1 - atom3    =%15.5f  Sigma atom3 - atom1    =%15.5f" % (parameters['RC1']['sigma_pk1pk3'], parameters['RC1']['sigma_pk3pk1']   )
            text = text + "\n--------------------------------------------------------------------------------"

        elif parameters['RC1']["rc_type"] == 'multiple_distance':
            text = text + "\n"
            text = text + "\n---------------------- Coordinate 1 - multiple-Distance ------------------------"	
            text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC1']['ATOMS'][0]    , parameters['RC1']['ATOM_NAMES'][0] )
            text = text + "\nATOM2*                 =%15i  ATOM NAME2             =%15s"     % (parameters['RC1']['ATOMS'][1]    , parameters['RC1']['ATOM_NAMES'][1] )
            text = text + "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC1']['ATOMS'][2]    , parameters['RC1']['ATOM_NAMES'][2] )
            text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['RC1']['nsteps']      , parameters['RC1']['force_constant'] ) 
            text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC1']['dminimum']    , parameters['maximumIterations']         )
            text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC1']['dincre']      , parameters['rmsGradientTolerance']           )
            text = text + "\nSigma atom1 - atom3    =%15.5f  Sigma atom3 - atom1    =%15.5f" % (parameters['RC1']['sigma_pk1pk3'], parameters['RC1']['sigma_pk3pk1']   )
            text = text + "\n--------------------------------------------------------------------------------"

        else:
            pass
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------
        
        
        
        
        if parameters["RC2"] is not None :

            '''This part writes the parameters used in the second reaction coordinate'''
            #-----------------------------------------------------------------------------------------------------------------------------------------------------------
            if parameters['RC2']["rc_type"] == 'simple_distance':
                text = text + "\n"
                text = text + "\n----------------------- Coordinate 2 - Simple-Distance -------------------------"
                text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC2']['ATOMS'][0], parameters['RC2']['ATOM_NAMES'][0] )
                text = text + "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['RC2']['ATOMS'][1], parameters['RC2']['ATOM_NAMES'][1] )
                text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['RC2']['nsteps']  , parameters['RC2']['force_constant'] )
                text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC2']['dminimum'], parameters['maximumIterations']         )
                text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC2']['dincre']  , parameters['rmsGradientTolerance']           )
                text = text + "\n--------------------------------------------------------------------------------"

            
            elif parameters['RC2']["rc_type"] == 'multiple_distance*4atoms':
                text = text + "\n"
                text = text + "\n---------------------- Coordinate 2 - multiple-Distance ------------------------"	
                text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC2']['ATOMS'][0]    , parameters['RC2']['ATOM_NAMES'][0] )
                text = text + "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['RC2']['ATOMS'][1]    , parameters['RC2']['ATOM_NAMES'][1] )
                text = text + "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC2']['ATOMS'][2]    , parameters['RC2']['ATOM_NAMES'][2] )
                text = text + "\nATOM4                  =%15i  ATOM NAME4             =%15s"     % (parameters['RC2']['ATOMS'][3]    , parameters['RC2']['ATOM_NAMES'][3] )
                text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['RC2']['nsteps']      , parameters['RC2']['force_constant'] )
                text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC2']['dminimum']    , parameters['maximumIterations']         )
                text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC2']['dincre']      , parameters['rmsGradientTolerance']           )
                text = text + "\n--------------------------------------------------------------------------------"
            
            elif parameters['RC2']["rc_type"] == 'multiple_distance':
                text = text + "\n"
                text = text + "\n---------------------- Coordinate 2 - multiple-Distance ------------------------"	
                text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC2']['ATOMS'][0]    , parameters['RC2']['ATOM_NAMES'][0] )
                text = text + "\nATOM2*                 =%15i  ATOM NAME2             =%15s"     % (parameters['RC2']['ATOMS'][1]    , parameters['RC2']['ATOM_NAMES'][1] )
                text = text + "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC2']['ATOMS'][2]    , parameters['RC2']['ATOM_NAMES'][2] )
                text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['RC2']['nsteps']      , parameters['RC2']['force_constant'] )
                text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC2']['dminimum']    , parameters['maximumIterations']         )
                text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC2']['dincre']      , parameters['rmsGradientTolerance']           )
                text = text + "\nSigma atom1 - atom3    =%15.5f  Sigma atom3 - atom1    =%15.5f" % (parameters['RC2']['sigma_pk1pk3'], parameters['RC2']['sigma_pk3pk1']   )
                text = text + "\n--------------------------------------------------------------------------------"
            else:
                pass
            #-----------------------------------------------------------------------------------------------------------------------------------------------------------
            
        
        
        
        
        
        
        
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------
        '''This part writes the header of the frames distances and energy'''
        if parameters["RC2"] is not None :
            text = text + "\n\n--------------------------------------------------------------------------------"
            text = text + "\n   Frame i  /  j        RCOORD-1             RCOORD-2                Energy     "
            text = text + "\n--------------------------------------------------------------------------------"




        else:
            if parameters['RC1']["rc_type"] == 'simple_distance':
                text = text + "\n\n-------------------------------------------------------------"
                text = text + "\n           Frame    dist-ATOM1-ATOM2             Energy      "
                text = text + "\n-------------------------------------------------------------"
            
            elif parameters['RC1']["rc_type"] == 'multiple_distance':
                text = text + "\n\n--------------------------------------------------------------------------------"
                text = text + "\n           Frame     dist-ATOM1-ATOM2      dist-ATOM2-ATOM3         Energy        "
                text = text + "\n--------------------------------------------------------------------------------  "
        
            elif parameters['RC1']["rc_type"] == 'multiple_distance*4atoms':
                text = text + "\n\n--------------------------------------------------------------------------------"
                text = text + "\n           Frame     dist-ATOM1-ATOM2      dist-ATOM3-ATOM4         Energy        "
                text = text + "\n--------------------------------------------------------------------------------  "
        
        
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------

            
        
        
        
            
            
            #text = text + "\n\n------------------------------------------------------"
            #text = text + "\n       Frame     distance-pK1-pK2         Energy      "
            #text = text + "\n------------------------------------------------------"
        
        arq.write(text)
        return arq

    
    def _run_scan_1D (self, parameters):
        """ Function doc """
        #-------------------------------------------------------------------------
        #Setting some local vars to ease the notation in the pDynamo methods
        #----------------------------------
        atom1 = parameters['RC1']['ATOMS'][0]
        atom2 = parameters['RC1']['ATOMS'][1]                   
        #---------------------------------
        restraints = RestraintModel()
        parameters['system'].DefineRestraintModel( restraints )                     
        #----------------------------------------------------------------------------------------
        
        arq = self.write_header(parameters)
        data = []
        
        for i in range(parameters['RC1']['nsteps']):       
            distance = parameters['RC1']['dminimum'] + ( parameters['RC1']['dincre'] * float(i) )
            #print(parameters["rc_type_1"])
            
            
            '''----------------------------------------------------------------------------------------------------------------'''
            if parameters['RC1']["rc_type"] == 'simple_distance':
                #---------------------------------------------------------------------------------------------------------
                rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
                restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom1, point2= atom2)
                restraints["RC1"] = restraint            
                #---------------------------------------------------------------------------------------------------------
            
                #-------------------------------------------------------------------------------------------------------------
                ConjugateGradientMinimize_SystemGeometry(parameters['system']                                ,                
                                                         logFrequency           = parameters['logFrequency'],
                                                         maximumIterations      = parameters['maximumIterations'],
                                                         rmsGradientTolerance   = parameters['rmsGradientTolerance'])
                #-------------------------------------------------------------------------------------------------------------
               
                distance = parameters['system'].coordinates3.Distance( atom1 , atom2  )
                energy   = parameters['system'].Energy(log=None)
                data.append([i, distance, energy])             
                
                text = "\nDATA %9i       %13.12f        %13.12f"% (int(i), float(distance), float(energy))
                arq.write(text)
                Pickle( os.path.join( parameters['folder'], 
                                      parameters['traj_folder_name']+".ptGeo", 
                                      "frame{}.pkl".format(i) ), 
                        parameters['system'].coordinates3 ) 
                
                backup_orca_files(system        = parameters['system'], 
                                  output_folder = os.path.join(parameters['folder'],parameters['traj_folder_name']+".ptGeo") , 
                                  output_name   = "frame{}".format(i))
        
                '''----------------------------------------------------------------------------------------------------------------'''

            
                '''----------------------------------------------------------------------------------------------------------------'''
            
            elif parameters['RC1']["rc_type"] == 'multiple_distance*4atoms':
                atom3   = parameters['RC1']['ATOMS'][2]
                atom4   = parameters['RC1']['ATOMS'][3]
                weight1 =  1.0#parameters['RC1']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
                weight2 = -1.0#parameters['RC1']['sigma_pk3pk1'] #self.sigma_a3_a1[0]
                rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
                restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom2, atom1, weight1 ], [ atom3, atom4, weight2 ] ] )
                restraints["RC1"] = restraint            
                
                #-------------------------------------------------------------------------------------------------------------
                ConjugateGradientMinimize_SystemGeometry(parameters['system']                                ,                
                                                         logFrequency           = parameters['logFrequency'],
                                                         maximumIterations      = parameters['maximumIterations'],
                                                         rmsGradientTolerance   = parameters['rmsGradientTolerance'])
                #-------------------------------------------------------------------------------------------------------------
                distance1 = parameters['system'].coordinates3.Distance( atom1 , atom2  )
                distance2 = parameters['system'].coordinates3.Distance( atom3 , atom4  )
                energy   = parameters['system'].Energy(log=None)
                data.append([i, distance1, distance2, energy])             
                
                text = "\nDATA %9i       %13.12f        %13.12f        %13.12f"% (int(i), float(distance1), float(distance2), float(energy))
                arq.write(text)
                
                Pickle( os.path.join( parameters['folder'], 
                                      parameters['traj_folder_name']+".ptGeo", 
                                      "frame{}.pkl".format(i) ), 
                        parameters['system'].coordinates3 ) 
                
                backup_orca_files(system        = parameters['system'], 
                                  output_folder = os.path.join(parameters['folder'],parameters['traj_folder_name']+".ptGeo") , 
                                  output_name   = "frame{}".format(i))
            
            
            elif parameters['RC1']["rc_type"] == 'multiple_distance':
                #--------------------------------------------------------------------
                atom3   = parameters['RC1']['ATOMS'][2]
                weight1 = parameters['RC1']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
                weight2 = parameters['RC1']['sigma_pk3pk1'] #self.sigma_a3_a1[0] 
                
                rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
                restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom2, atom1, weight1 ], [ atom2, atom3, weight2 ] ] )
                restraints["RC1"] = restraint            
                #--------------------------------------------------------------------
            
                #-------------------------------------------------------------------------------------------------------------
                ConjugateGradientMinimize_SystemGeometry(parameters['system']                                ,                
                                                         logFrequency           = parameters['logFrequency'],
                                                         maximumIterations      = parameters['maximumIterations'],
                                                         rmsGradientTolerance   = parameters['rmsGradientTolerance'])
                #-------------------------------------------------------------------------------------------------------------
            
                distance1 = parameters['system'].coordinates3.Distance( atom1 , atom2  )
                distance2 = parameters['system'].coordinates3.Distance( atom2 , atom3  )
                energy   = parameters['system'].Energy(log=None)
                data.append([i, distance1, distance2, energy])             
                
                text = "\nDATA %9i       %13.12f        %13.12f        %13.12f"% (int(i), float(distance1), float(distance2), float(energy))
                arq.write(text)
                
                
                Pickle( os.path.join( parameters['folder'], 
                                      parameters['traj_folder_name']+".ptGeo", 
                                      "frame{}.pkl".format(i) ), 
                        parameters['system'].coordinates3 ) 
                
                backup_orca_files(system        = parameters['system'], 
                                  output_folder = os.path.join(parameters['folder'],parameters['traj_folder_name']+".ptGeo") , 
                                  output_name   = "frame{}".format(i))
                
                
                '''----------------------------------------------------------------------------------------------------------------'''

            elif parameters['RC1']["rc_type"] == 'dihedral':
                atom3   = parameters['RC1']['ATOMS'][2]
                atom4   = parameters['RC1']['ATOMS'][3]
                
                angle = parameters['RC1']['dminimum'] + ( parameters['RC1']['dincre'] * float(i) )
                #print('\n\n\n\n\n\n\n\n',distance, parameters['system'].coordinates3.Dihedral (atom1, atom2, atom3, atom4), '\n\n\n\n\n\n\n\n')
                rmodel    = RestraintEnergyModel.Harmonic(angle, parameters['RC1']['force_constant'], period = 360.0 )
                restraint = RestraintDihedral.WithOptions ( energyModel = rmodel     ,
                                                            point1      = atom1 ,
                                                            point2      = atom2 ,
                                                            point3      = atom3 ,
                                                            point4      = atom4 )                
                restraints["RC1"] = restraint            
                #--------------------------------------------------------------------
            
                #-------------------------------------------------------------------------------------------------------------
                ConjugateGradientMinimize_SystemGeometry(parameters['system']                                ,                
                                                         logFrequency           = parameters['logFrequency'],
                                                         maximumIterations      = parameters['maximumIterations'],
                                                         rmsGradientTolerance   = parameters['rmsGradientTolerance'])
                #-------------------------------------------------------------------------------------------------------------
                angle0 = parameters['system'].coordinates3.Dihedral (atom1, atom2, atom3, atom4)
                #distance1 = parameters['system'].coordinates3.Distance( atom1 , atom2  )
                #distance2 = parameters['system'].coordinates3.Distance( atom2 , atom3  )
                energy   = parameters['system'].Energy(log=None)
                #data.append([i, distance1, distance2, energy])             
                
                text = "\nDATA %9i       %13.12f        %13.12f"% (int(i), float(angle0), float(energy))
                arq.write(text)

                Pickle( os.path.join( parameters['folder'], 
                                      parameters['traj_folder_name']+".ptGeo", 
                                      "frame{}.pkl".format(i) ), 
                        parameters['system'].coordinates3 ) 
                
                backup_orca_files(system        = parameters['system'], 
                                  output_folder = os.path.join(parameters['folder'],parameters['traj_folder_name']+".ptGeo") , 
                                  output_name   = "frame{}".format(i))



            else:
                pass
        #---------------------------------------
        parameters['system'].DefineRestraintModel(None)
        pprint(data)
        
        
    def _run_scan_2D (self, parameters = None, interface = False):
        """ Function doc """
        #-------------------------------------------------------------------------
        #Setting some local vars to ease the notation in the pDynamo methods
        #----------------------------------
        atom_RC1_1 = parameters['RC1']['ATOMS'][0]
        atom_RC1_2 = parameters['RC1']['ATOMS'][1]                   
        
        atom_RC2_1 = parameters['RC2']['ATOMS'][0]
        atom_RC2_2 = parameters['RC2']['ATOMS'][1]                   
        
        print(atom_RC1_1, atom_RC1_2, atom_RC2_1, atom_RC2_2)
        print(parameters['RC1']['ATOMS'], parameters['RC2']['ATOMS'])
        
        #return True
        #---------------------------------
        restraints = RestraintModel()
        parameters['system'].DefineRestraintModel( restraints )                     
        #----------------------------------------------------------------------------------------
        
        arq = self.write_header(parameters)
        data = {}
        #print(parameters["rc_type_1"])
        
        
        '''This part of the program is crucial because the Clone function 
        is unable to handle data that has been allocated using the GTK library. 
        When a system is generalized using EasyHybrid, it receives information 
        from the treeview, enabling easy interconnection of the data.'''
        #if interface:
        backup = []
        try:
            backup.append(parameters['system'].e_treeview_iter)
            backup.append(parameters['system'].e_liststore_iter)
            parameters['system'].e_treeview_iter   = None
            parameters['system'].e_liststore_iter  = None
        except:
            pass
       
        
        
        '''We have to first run a sequence (first line) of optimizations sequentially. 
        This will generate frames (0,0) to (n, 0) serially, where n is the number of 
        steps from coordinate x)
        '''
        #try:
        joblist = []
        j = 0
        for i in range(parameters['RC1']['nsteps']):       
            distance = parameters['RC1']['dminimum'] + ( parameters['RC1']['dincre'] * float(i) )
            
            '''reaction coordinate 1 - starts at 0 and goes to nstpes'''
            '''----------------------------------------------------------------------------------------------------------------'''
            if parameters['RC1']["rc_type"] == 'simple_distance':
                #---------------------------------------------------------------------------------------------------------
                rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
                restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom_RC1_1, point2= atom_RC1_2)
                restraints["RC1"] = restraint            
                #---------------------------------------------------------------------------------------------------------
            
            elif parameters['RC1']["rc_type"] == 'multiple_distance*4atoms':
                atom_RC1_3   = parameters['RC1']['ATOMS'][2]
                atom_RC1_4   = parameters['RC1']['ATOMS'][3]
                weight1 =  1.0#parameters['RC1']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
                weight2 = -1.0#parameters['RC1']['sigma_pk3pk1'] #self.sigma_a3_a1[0]
                rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
                restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom_RC1_1, atom_RC1_2, weight1 ], 
                                                                                                              [ atom_RC1_3, atom_RC1_4, weight2 ] 
                                                                                                              ] )
                restraints["RC1"] = restraint  
            
            elif parameters['RC1']["rc_type"] == 'multiple_distance':
                #--------------------------------------------------------------------
                atom_RC1_3   = parameters['RC1']['ATOMS'][2]
                weight1 = parameters['RC1']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
                weight2 = parameters['RC1']['sigma_pk3pk1'] #self.sigma_a3_a1[0] 
                
                rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
                restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom_RC1_2, atom_RC1_1, weight1 ], 
                                                                                                              [ atom_RC1_2, atom_RC1_3, weight2 ] 
                                                                                                            ] )
                restraints["RC1"] = restraint            
                #--------------------------------------------------------------------
            else:
                pass
                


            '''reaction coordinate 2 - ONLY at 0 first'''
            '''----------------------------------------------------------------------------------------------------------------'''
            distance2 = parameters['RC2']['dminimum'] + ( parameters['RC2']['dincre'] * float(j) )
            
            if parameters['RC2']["rc_type"] == 'simple_distance':
                #---------------------------------------------------------------------------------------------------------
                rmodel            = RestraintEnergyModel.Harmonic(distance2, parameters['RC2']['force_constant'])
                restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom_RC2_1, point2= atom_RC2_2)
                restraints["RC2"] = restraint            
                #---------------------------------------------------------------------------------------------------------
            
            
            elif parameters['RC2']["rc_type"] == 'multiple_distance*4atoms':
                atom_RC2_3   = parameters['RC2']['ATOMS'][2]
                atom_RC2_4   = parameters['RC2']['ATOMS'][3]
                weight1 =  1.0#parameters['RC1']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
                weight2 = -1.0#parameters['RC1']['sigma_pk3pk1'] #self.sigma_a3_a1[0]
                rmodel            = RestraintEnergyModel.Harmonic(distance2, parameters['RC2']['force_constant'])
                restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom_RC2_1, atom_RC2_2, weight1 ], 
                                                                                                              [ atom_RC2_3, atom_RC2_4, weight2 ] 
                                                                                                              ] )
                restraints["RC2"] = restraint  
            
            
            
            elif parameters['RC2']["rc_type"] == 'multiple_distance':
                #--------------------------------------------------------------------
                atom_RC2_3 = parameters['RC2']['ATOMS'][2]
                weight1 = parameters['RC2']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
                weight2 = parameters['RC2']['sigma_pk3pk1'] #self.sigma_a3_a1[0] 
                
                rmodel            = RestraintEnergyModel.Harmonic(distance2, parameters['RC2']['force_constant'])
                restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ 
                                                                                                              [ atom_RC2_2, atom_RC2_1, weight1 ], 
                                                                                                              [ atom_RC2_2, atom_RC2_3, weight2 ] 
                                                                                                              ] )
                restraints["RC2"] = restraint            
                #--------------------------------------------------------------------
           
            
            
            
            try:
                #-------------------------------------------------------------------------------------------------------------
                ConjugateGradientMinimize_SystemGeometry(parameters['system']                                ,                
                                                         logFrequency           = parameters['logFrequency'],
                                                         maximumIterations      = parameters['maximumIterations'],
                                                         rmsGradientTolerance   = parameters['rmsGradientTolerance']  ,
                                                         log                    = None)

                energy   = parameters['system'].Energy(log=None)
                opt_convergency = ''
                #-------------------------------------------------------------------------------------------------------------
            except:
                energy   = 0
                opt_convergency = 'Ops! Error in geometry optimization.'

            #--------------------------------------------------------------------------------------
            #                     Calculating the reaction coordinate 1
            #--------------------------------------------------------------------------------------
            distance1 = parameters['system'].coordinates3.Distance( atom_RC1_1 , atom_RC1_2  )
            
            if parameters['RC1']["rc_type"] == 'multiple_distance':
                distance2 = parameters['system'].coordinates3.Distance( atom_RC1_2,  atom_RC1_3  )

            elif parameters['RC1']["rc_type"] == 'multiple_distance*4atoms':
                distance2 = parameters['system'].coordinates3.Distance( atom_RC1_3,  atom_RC1_4  )
            
            else:
                distance2 = 0
            
            RC1_d1_minus_d2 = distance1 - distance2
            #--------------------------------------------------------------------------------------
            
            
            
            #--------------------------------------------------------------------------------------
            #                     Calculating the reaction coordinate 2
            #--------------------------------------------------------------------------------------
            distance1 = parameters['system'].coordinates3.Distance( atom_RC2_1 , atom_RC2_2  )
            
            if parameters['RC2']["rc_type"] == 'multiple_distance':
                distance2 = parameters['system'].coordinates3.Distance( atom_RC2_2,  atom_RC2_3  )
            
            elif parameters['RC2']["rc_type"] == 'multiple_distance*4atoms':
                distance2 = parameters['system'].coordinates3.Distance( atom_RC2_3,  atom_RC2_4  )
            
            else:
                distance2 = 0
            
            RC2_d1_minus_d2 = distance1 - distance2
            #--------------------------------------------------------------------------------------
            
            
            
            #--------------------------------------------------------------------------------------
            #                          Energy and data dictionary
            #--------------------------------------------------------------------------------------
            data[(i,j)] = [RC1_d1_minus_d2, RC2_d1_minus_d2, energy]    
            #--------------------------------------------------------------------------------------
            
            
            
            #--------------------------------------------------------------------------------------
            #                             Exporting Coordinates
            #--------------------------------------------------------------------------------------
            pkl = os.path.join( parameters['folder'], 
                                  parameters['traj_folder_name']+".ptGeo", 
                                  "frame{}_{}.pkl".format(i,j) )
                                   
            Pickle( pkl, parameters['system'].coordinates3 )
            
            backup_orca_files(system        = parameters['system'], 
                              output_folder = os.path.join(parameters['folder'],parameters['traj_folder_name']+".ptGeo") , 
                              output_name   = "frame{}_{}".format(i,j))
            
            #--------------------------------------------------------------------------------------
            
            

            #system = parameters['system']
            #backup = []
            #backup.append(system.e_treeview_iter)
            #backup.append(system.e_liststore_iter)
            
            #For some reason the previous dictionary had problems when submitted to the pool
            new_parameters = {
                             
                             'ATOMS_RC1'       : parameters['RC1']['ATOMS'],
                             'ATOMS_RC2'       : parameters['RC2']['ATOMS'],
                             'rc_type_1'       : parameters['RC1']['rc_type'],
                             'rc_type_2'       : parameters['RC2']['rc_type'],
                             'dminimum_RC1'    : parameters['RC1']['dminimum'],
                             'dminimum_RC2'    : parameters['RC2']['dminimum'],
                             
                             'sigma_pk1pk3_rc1': parameters['RC1']['sigma_pk1pk3'],
                             'sigma_pk3pk1_rc1': parameters['RC1']['sigma_pk3pk1'],
                             'force_constant_1': parameters['RC1']['force_constant'],
                             'nsteps_RC1'      : parameters['RC1']['nsteps'],
                             'dincre_RC1'      : parameters['RC1']['dincre'],
                             
                             
                             'sigma_pk1pk3_rc2': parameters['RC2']['sigma_pk1pk3'],
                             'sigma_pk3pk1_rc2': parameters['RC2']['sigma_pk3pk1'],
                             'force_constant_2': parameters['RC2']['force_constant'],
                             'nsteps_RC2'      : parameters['RC2']['nsteps'],
                             'dincre_RC2'      : parameters['RC2']['dincre'],
                             
                             'logFrequency'   : parameters['logFrequency'],
                             'maximumIterations'   : parameters['maximumIterations'],
                             'rmsGradientTolerance'     : parameters['rmsGradientTolerance'],

                             'folder'          : parameters['folder']    , 
                             'traj_folder_name': parameters['traj_folder_name']
                             }
            
            print(i, j, 'Energy:', energy, opt_convergency)
            #  - - - - adding to the joblist - - - - 
            joblist.append([i, pkl, new_parameters, parameters['system']])
            #--------------------------------------------------------------------------------------


        #print(joblist)
        p = multiprocessing.Pool(processes=parameters['NmaxThreads'])
        #p = multiprocessing.Pool(processes=8)
        #results = p.map(self._run_second_coordinate_in_parallel, joblist)
        results = p.map(_run_second_coordinate_in_parallel, joblist)
        
        

        
        for partial_data_dict in  results:
            for key, partial_data in partial_data_dict.items():
                data[key] = partial_data

        #--------------------------------------------------------------------------------------
        #                             Writing the log information 
        #--------------------------------------------------------------------------------------
        for i in range(parameters['RC1']['nsteps']):
            for j in range(parameters['RC2']['nsteps']):
        
                text = "\nDATA  %4i  %4i     %13.12f       %13.12f       %13.12f"% (int(i), int(j),  float(data[(i,j)][0]), float(data[(i,j)][1]), float(data[(i,j)][2]))
                arq.write(text)
        #--------------------------------------------------------------------------------------
        

        #if interface:
        try:
            parameters['system'].e_treeview_iter   = backup[0]
            parameters['system'].e_liststore_iter  = backup[1]
        except:
            pass
    
    
    '''
    def _run_scan_2D_old (self, parameters):
        """ Function doc """
        #-------------------------------------------------------------------------
        #Setting some local vars to ease the notation in the pDynamo methods
        #----------------------------------
        atom_RC1_1 = parameters['ATOMS_RC1'][0]
        atom_RC1_2 = parameters['ATOMS_RC1'][1]                   
        
        atom_RC2_1 = parameters['ATOMS_RC2'][0]
        atom_RC2_2 = parameters['ATOMS_RC2'][1]                   
        #---------------------------------
        restraints = RestraintModel()
        parameters['system'].DefineRestraintModel( restraints )                     
        #----------------------------------------------------------------------------------------
        
        arq = self.write_header(parameters)
        data = {}
        #print(parameters["rc_type_1"])

        for i in range(parameters['nsteps_RC1']):       
            distance = parameters['dminimum_RC1'] + ( parameters['dincre_RC1'] * float(i) )

            if parameters["rc_type_1"] == 'simple_distance':
                #---------------------------------------------------------------------------------------------------------
                rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['force_constant_1'])
                restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom_RC1_1, point2= atom_RC1_2)
                restraints["RC1"] = restraint            
                #---------------------------------------------------------------------------------------------------------
            
            elif parameters["rc_type_1"] == 'multiple_distance':
                #--------------------------------------------------------------------
                atom_RC1_3   = parameters['ATOMS_RC1'][2]
                weight1 = parameters['sigma_pk1pk3_rc1'] #self.sigma_a1_a3[0]
                weight2 = parameters['sigma_pk3pk1_rc1'] #self.sigma_a3_a1[0] 
                
                rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['force_constant_1'])
                restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom_RC1_2, atom_RC1_1, weight1 ], 
                                                                                                              [ atom_RC1_2, atom_RC1_3, weight2 ] 
                                                                                                            ] )
                restraints["RC1"] = restraint            
                #--------------------------------------------------------------------
            else:
                pass
                

            for j in range(parameters['nsteps_RC2']):       
                distance2 = parameters['dminimum_RC2'] + ( parameters['dincre_RC2'] * float(j) )
                
                if parameters["rc_type_2"] == 'simple_distance':
                    #---------------------------------------------------------------------------------------------------------
                    rmodel            = RestraintEnergyModel.Harmonic(distance2, parameters['force_constant_2'])
                    restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom_RC2_1, point2= atom_RC2_2)
                    restraints["RC2"] = restraint            
                    #---------------------------------------------------------------------------------------------------------
                
                elif parameters["rc_type_2"] == 'multiple_distance':
                    #--------------------------------------------------------------------
                    atom_RC2_3 = parameters['ATOMS_RC2'][2]
                    weight1 = parameters['sigma_pk1pk3_rc2'] #self.sigma_a1_a3[0]
                    weight2 = parameters['sigma_pk3pk1_rc2'] #self.sigma_a3_a1[0] 
                    
                    rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['force_constant_2'])
                    restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ 
                                                                                                                  [ atom_RC2_2, atom_RC2_1, weight1 ], 
                                                                                                                  [ atom_RC2_2, atom_RC2_3, weight2 ] 
                                                                                                                  ] )
                    restraints["RC2"] = restraint            
                    #--------------------------------------------------------------------
               
                
                #-------------------------------------------------------------------------------------------------------------
                ConjugateGradientMinimize_SystemGeometry(parameters['system']                                ,                
                                                         logFrequency           = parameters['logFrequency'],
                                                         maximumIterations      = parameters['maximumIterations'],
                                                         rmsGradientTolerance   = parameters['rmsGradientTolerance'])
                #-------------------------------------------------------------------------------------------------------------

                
                #--------------------------------------------------------------------------------------
                #                     Calculating the reaction coordinate 1
                #--------------------------------------------------------------------------------------
                distance1 = parameters['system'].coordinates3.Distance( atom_RC1_1 , atom_RC1_2  )
                if parameters["rc_type_1"] == 'multiple_distance':
                    distance2 = parameters['system'].coordinates3.Distance( atom_RC1_2,  atom_RC1_3  )
                else:
                    distance2 = 0
                RC1_d1_minus_d2 = distance1 - distance2
                #--------------------------------------------------------------------------------------
                
                
                
                #--------------------------------------------------------------------------------------
                #                     Calculating the reaction coordinate 2
                #--------------------------------------------------------------------------------------
                distance1 = parameters['system'].coordinates3.Distance( atom_RC2_1 , atom_RC2_2  )
                if parameters["rc_type_2"] == 'multiple_distance':
                    distance2 = parameters['system'].coordinates3.Distance( atom_RC2_2,  atom_RC2_3  )
                else:
                    distance2 = 0
                RC2_d1_minus_d2 = distance1 - distance2
                #--------------------------------------------------------------------------------------
                
                
                
                #--------------------------------------------------------------------------------------
                #                          Energy and data dictionary
                #--------------------------------------------------------------------------------------
                energy   = parameters['system'].Energy(log=None)
                data[(i,j)] = [RC1_d1_minus_d2, RC2_d1_minus_d2, energy]    
                #--------------------------------------------------------------------------------------
                
                
                
                #--------------------------------------------------------------------------------------
                #                             Exporting Coordinates
                #--------------------------------------------------------------------------------------
                Pickle( os.path.join( parameters['folder'], 
                                      parameters['traj_folder_name']+".ptGeo", 
                                      "frame{}_{}.pkl".format(i,j) ), 
                                      parameters['system'].coordinates3 )
                #--------------------------------------------------------------------------------------
                backup_orca_files(system        = parameters['system'], 
                                  output_folder = os.path.join(parameters['folder'],parameters['traj_folder_name']+".ptGeo") , 
                                  output_name   = "frame{}_{}".format(i,j))
        
        
        #--------------------------------------------------------------------------------------
        #                             Writing the log information 
        #--------------------------------------------------------------------------------------
        for i in range(parameters['nsteps_RC1']):
            for j in range(parameters['nsteps_RC2']):
        
                text = "\nDATA  %4i  %4i     %13.12f       %13.12f       %13.12f"% (int(i), int(j),  float(data[(i,j)][0]), float(data[(i,j)][1]), float(data[(i,j)][2]))
                arq.write(text)
        #--------------------------------------------------------------------------------------
    '''


def _run_advanced_second_coordinate_in_parallel (job):
    """ Function doc """
    #print('hello')
    i          = job[0]
    pkl        = job[1]
    parameters = job[2]
    system     = job[3]
    
    #----------------------------------------------------------------------------------------
    RC = parameters['RC1']['RC']
    RC1 = []
    dminimum = 0.0
    for rc in RC:
        dist = float(rc[5])*float(rc[4])
        dminimum += dist
        RC1.append([ int(rc[1]),  int(rc[3]), float(rc[4])])
    parameters['RC1']['dminimum'] = dminimum
    #----------------------------------------------------------------------------------------
    
    
    #----------------------------------------------------------------------------------------
    RC = parameters['RC2']['RC']
    RC2 = []
    dminimum2 = 0.0
    for rc in RC:
        dist = float(rc[5])*float(rc[4])
        dminimum2 += dist
        RC2.append([ int(rc[1]),  int(rc[3]), float(rc[4])])
    parameters['RC2']['dminimum'] = dminimum2
    #----------------------------------------------------------------------------------------
    
    hamiltonian = get_hamiltonian (system)

    if hamiltonian in ['DFTB QC Model', 'ORCA QC Model', 'XTB QC Model', 'external']: #
        try:
            os.mkdir(system.qcModel.scratch +'/process_'+str(i))
        except:
            pass
     
        try:
            system.qcState.DeterminePaths(system.qcModel.scratch +'/process_'+str(i))
        except:
            pass


    restraints = RestraintModel()
    system.DefineRestraintModel( restraints )                     
    #----------------------------------------------------------------------------------------
    system.coordinates3 = ImportCoordinates3(pkl )
    #arq = self.write_header(parameters)
    data = {}
    
    
    for j in range(parameters['nsteps_RC2']):       
        
        distance = parameters['dminimum_RC1'] + ( parameters['dincre_RC1'] * float(i) )
        
        '''reaction coordinate 1 - starts at 0 and goes to nstpes'''
        '''----------------------------------------------------------------------------------------------------------------'''
        rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
        restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = RC1 )
        restraints["RC1"] = restraint 
        '''----------------------------------------------------------------------------------------------------------------'''

        

        '''reaction coordinate 2 '''
        distance2 = parameters['dminimum_RC2'] + ( parameters['dincre_RC2'] * float(j) )
        rmodel            = RestraintEnergyModel.Harmonic(distance2, parameters['RC2']['force_constant'])
        restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances = RC2 )
        restraints["RC2"] = restraint         

        #ConjugateGradientMinimize_SystemGeometry(system                                ,                
        #                                         logFrequency           = parameters['logFrequency'],
        #                                         maximumIterations      = parameters['maximumIterations'],
        #                                         rmsGradientTolerance   = parameters['rmsGradientTolerance'], 
        #                                         log= None) 
        #energy = system.Energy(log=None)
        try:
            #-------------------------------------------------------------------------------------------------------------
            ConjugateGradientMinimize_SystemGeometry(system                                ,                
                                                     logFrequency           = parameters['logFrequency'],
                                                     maximumIterations      = parameters['maximumIterations'],
                                                     rmsGradientTolerance   = parameters['rmsGradientTolerance'], 
                                                     log                    =                      None)
            #-------------------------------------------------------------------------------------------------------------
            energy   = system.Energy(log=None)
            opt_convergency = ''
            #-------------------------------------------------------------------------------------------------------------
        except:
            energy   = 0
            opt_convergency = 'Ops! Error in geometry optimization.'
        

        #--------------------------------------------------------------------------------------
        #                          Energy and data dictionary
        #--------------------------------------------------------------------------------------
        data[(i,j)] = [distance, distance2, energy]    
        #--------------------------------------------------------------------------------------
        
        
        
        #--------------------------------------------------------------------------------------
        #                             Exporting Coordinates
        #--------------------------------------------------------------------------------------
        Pickle( os.path.join( parameters['folder'], 
                              parameters['traj_folder_name']+".ptGeo", 
                              "frame{}_{}.pkl".format(i,j) ), 
                              
                              system.coordinates3 )
                              
        backup_orca_files(system        = system, 
                          output_folder = os.path.join(parameters['folder'],parameters['traj_folder_name']+".ptGeo") , 
                          output_name   = "frame{}_{}".format(i,j))
        #--------------------------------------------------------------------------------------
        #print(i, j, 'Energy:', energy, opt_convergency)
        print(i, j, 'Energy:', energy)
    return data


def _run_second_coordinate_in_parallel (job):
    """ Function doc """
    #print('hello')
    i          = job[0]
    pkl        = job[1]
    parameters = job[2]
    system     = job[3]
    
    
    hamiltonian = get_hamiltonian (system)
    ''' 
    
    '''
    if hamiltonian in ['DFTB QC Model', 'ORCA QC Model', 'XTB QC Model', 'external']: #
        try:
            os.mkdir(system.qcModel.scratch +'/process_'+str(i))
        except:
            pass
     
        try:
            system.qcState.DeterminePaths(system.qcModel.scratch +'/process_'+str(i))
        except:
            pass

    
    
    
    atom_RC1_1 = parameters['ATOMS_RC1'][0]
    atom_RC1_2 = parameters['ATOMS_RC1'][1]                   
    
    atom_RC2_1 = parameters['ATOMS_RC2'][0]
    atom_RC2_2 = parameters['ATOMS_RC2'][1]  
    
    restraints = RestraintModel()
    system.DefineRestraintModel( restraints )                     
    #----------------------------------------------------------------------------------------
    system.coordinates3 = ImportCoordinates3(pkl )
    #arq = self.write_header(parameters)
    data = {}
    
    
    for j in range(parameters['nsteps_RC2']):       
        
        distance = parameters['dminimum_RC1'] + ( parameters['dincre_RC1'] * float(i) )
        
        '''reaction coordinate 1 - starts at 0 and goes to nstpes'''
        '''----------------------------------------------------------------------------------------------------------------'''
        if parameters["rc_type_1"] == 'simple_distance':
            #---------------------------------------------------------------------------------------------------------
            rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['force_constant_1'])
            restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom_RC1_1, point2= atom_RC1_2)
            restraints["RC1"] = restraint            
            #---------------------------------------------------------------------------------------------------------
        
        elif parameters["rc_type_1"] == 'multiple_distance*4atoms':
            atom_RC1_3   = parameters['ATOMS_RC1'][2]
            atom_RC1_4   = parameters['ATOMS_RC1'][3]
            weight1 =  1.0 #parameters['sigma_pk1pk3_rc1'] #self.sigma_a1_a3[0]
            weight2 = -1.0 #parameters['sigma_pk3pk1_rc1'] #self.sigma_a3_a1[0] 
            
            rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['force_constant_1'])
            restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom_RC1_2, atom_RC1_1, weight1 ], 
                                                                                                          [ atom_RC1_3, atom_RC1_4, weight2 ] 
                                                                                                        ] )
            restraints["RC1"] = restraint 
        
        elif parameters["rc_type_1"] == 'multiple_distance':
            #--------------------------------------------------------------------
            atom_RC1_3   = parameters['ATOMS_RC1'][2]
            weight1 = parameters['sigma_pk1pk3_rc1'] #self.sigma_a1_a3[0]
            weight2 = parameters['sigma_pk3pk1_rc1'] #self.sigma_a3_a1[0] 
            
            rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['force_constant_1'])
            restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ atom_RC1_2, atom_RC1_1, weight1 ], 
                                                                                                          [ atom_RC1_2, atom_RC1_3, weight2 ] 
                                                                                                        ] )
            restraints["RC1"] = restraint            
            #--------------------------------------------------------------------
        else:
            pass
        '''----------------------------------------------------------------------------------------------------------------'''

        

        
        distance2 = parameters['dminimum_RC2'] + ( parameters['dincre_RC2'] * float(j) )
        
        if parameters["rc_type_2"] == 'simple_distance':
            #---------------------------------------------------------------------------------------------------------
            rmodel            = RestraintEnergyModel.Harmonic(distance2, parameters['force_constant_2'])
            restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom_RC2_1, point2= atom_RC2_2)
            restraints["RC2"] = restraint            
            #---------------------------------------------------------------------------------------------------------
        
        elif parameters["rc_type_2"] == 'multiple_distance*4atoms':
            #--------------------------------------------------------------------
            atom_RC2_3 = parameters['ATOMS_RC2'][2]
            atom_RC2_4 = parameters['ATOMS_RC2'][3]
            weight1 =   1.0 #parameters['sigma_pk1pk3_rc2'] #self.sigma_a1_a3[0]
            weight2 =  -1.0 #parameters['sigma_pk3pk1_rc2'] #self.sigma_a3_a1[0] 
            
            rmodel            = RestraintEnergyModel.Harmonic(distance2, parameters['force_constant_2'])
            restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ 
                                                                                                          [ atom_RC2_2, atom_RC2_1, weight1 ], 
                                                                                                          [ atom_RC2_3, atom_RC2_4, weight2 ] 
                                                                                                          ] )
            restraints["RC2"] = restraint            
            #--------------------------------------------------------------------
       
        elif parameters["rc_type_2"] == 'multiple_distance':
            #--------------------------------------------------------------------
            atom_RC2_3 = parameters['ATOMS_RC2'][2]
            weight1 = parameters['sigma_pk1pk3_rc2'] #self.sigma_a1_a3[0]
            weight2 = parameters['sigma_pk3pk1_rc2'] #self.sigma_a3_a1[0] 
            
            rmodel            = RestraintEnergyModel.Harmonic(distance2, parameters['force_constant_2'])
            restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ 
                                                                                                          [ atom_RC2_2, atom_RC2_1, weight1 ], 
                                                                                                          [ atom_RC2_2, atom_RC2_3, weight2 ] 
                                                                                                          ] )
            restraints["RC2"] = restraint            
            #--------------------------------------------------------------------
       
        else:
            pass
        
        ##-------------------------------------------------------------------------------------------------------------
        #ConjugateGradientMinimize_SystemGeometry(system                                ,                
        #                                         logFrequency           = parameters['logFrequency'],
        #                                         maximumIterations      = parameters['maximumIterations'],
        #                                         rmsGradientTolerance   = parameters['rmsGradientTolerance'], 
        #                                         log                    =                      None)
        ##-------------------------------------------------------------------------------------------------------------
        try:
            #-------------------------------------------------------------------------------------------------------------
            ConjugateGradientMinimize_SystemGeometry(system                                ,                
                                                     logFrequency           = parameters['logFrequency'],
                                                     maximumIterations      = parameters['maximumIterations'],
                                                     rmsGradientTolerance   = parameters['rmsGradientTolerance'], 
                                                     log                    =                      None)
            #-------------------------------------------------------------------------------------------------------------
            energy   = system.Energy(log=None)
            opt_convergency = ''
            #-------------------------------------------------------------------------------------------------------------
        except:
            energy   = 0
            opt_convergency = 'Ops! Error in geometry optimization.'
        
        #--------------------------------------------------------------------------------------
        #                     Calculating the reaction coordinate 1
        #--------------------------------------------------------------------------------------
        distance1 = system.coordinates3.Distance( atom_RC1_1 , atom_RC1_2  )
        if parameters["rc_type_1"] == 'multiple_distance':
            distance2 = system.coordinates3.Distance( atom_RC1_2,  atom_RC1_3  )
        else:
            distance2 = 0
        RC1_d1_minus_d2 = distance1 - distance2
        #--------------------------------------------------------------------------------------
        
        
        
        #--------------------------------------------------------------------------------------
        #                     Calculating the reaction coordinate 2
        #--------------------------------------------------------------------------------------
        distance1 = system.coordinates3.Distance( atom_RC2_1 , atom_RC2_2  )
        if parameters["rc_type_2"] == 'multiple_distance':
            distance2 = system.coordinates3.Distance( atom_RC2_2,  atom_RC2_3  )
        else:
            distance2 = 0
        RC2_d1_minus_d2 = distance1 - distance2
        #--------------------------------------------------------------------------------------
        
        
        
        #--------------------------------------------------------------------------------------
        #                          Energy and data dictionary
        #--------------------------------------------------------------------------------------
        #energy   = system.Energy(log=None)
        data[(i,j)] = [RC1_d1_minus_d2, RC2_d1_minus_d2, energy]    
        #--------------------------------------------------------------------------------------
        
        
        
        #--------------------------------------------------------------------------------------
        #                             Exporting Coordinates
        #--------------------------------------------------------------------------------------
        Pickle( os.path.join( parameters['folder'], 
                              parameters['traj_folder_name']+".ptGeo", 
                              "frame{}_{}.pkl".format(i,j) ), 
                              
                              system.coordinates3 )
                              
        backup_orca_files(system        = system, 
                          output_folder = os.path.join(parameters['folder'],parameters['traj_folder_name']+".ptGeo") , 
                          output_name   = "frame{}_{}".format(i,j))
        #--------------------------------------------------------------------------------------
        print(i, j, 'Energy:', energy, opt_convergency)
    return data


def _run_advanded_scan_1D ( parameters = None, interface = False):
    '''under construction'''
    pass


def _run_advanded_scan_2D ( parameters = None, interface = False):
    """ Function doc """
    #-------------------------------------------------------------------------
    #Setting some local vars to ease the notation in the pDynamo methods
    #----------------------------------
                

    n_south = parameters['RC1']['up_nsteps'  ]
    n_north = parameters['RC1']['down_nsteps']
              
    n_left  = parameters['RC2']['left_nsteps' ]
    n_right = parameters['RC2']['right_nsteps']
    
    
    
    full_path_trajectory =os.path.join(parameters['folder'], 
                          parameters['traj_folder_name']+".ptGeo")
    #os.mkdir(
    #         full_path_trajectory
    #         )
    
    # - - - - - - - - - - - - - Checking trajectory - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    logFile2 = TextLogFileWriter.WithOptions ( path = os.path.join(full_path_trajectory, 'output.log') )
    parameters['system'].Summary(log = logFile2)
    logFile2.Header ( )
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    arq = write_header(parameters)


    '''the total number of steps in both directions'''
    i_size = n_south + n_north
    j_size = n_left + n_right

    #print(parameters['RC1']['ATOMS'][0], parameters['RC1']['ATOMS'][1], parameters['RC2']['ATOMS'][0], parameters['RC2']['ATOMS'][1])
    #print(parameters['RC1']['ATOMS'], parameters['RC2']['ATOMS'])

    '''matrix_i_j is the energy matrix (i x j)-  starts with only zeros'''
    matrix_i_j = []#i_size*[[0.0]*j_size]
    for i in range(i_size):
        matrix_i_j.append([])
        for j in range(j_size):
            matrix_i_j[i].append(0.0)


    #return True
    #----------------------------------------------------------------------------------------
    restraints = RestraintModel()
    parameters['system'].DefineRestraintModel( restraints )                     
    #----------------------------------------------------------------------------------------
    
    #arq = self.write_header(parameters)
    
    data = {}
    
    
    '''This part of the program is crucial because the Clone function 
    is unable to handle data that has been allocated using the GTK library. 
    When a system is generalized using EasyHybrid, it receives information 
    from the treeview, enabling easy interconnection of the data.'''
    #if interface:
    backup = []
    try:
        backup.append(parameters['system'].e_treeview_iter)
        backup.append(parameters['system'].e_liststore_iter)
        parameters['system'].e_treeview_iter   = None
        parameters['system'].e_liststore_iter  = None
    except:
        pass

    '''
    o------o------o------o------o------o
    |      |      |      |      |      |
    |      |      |      |      |      |
    o------o------o------o------o------o
    |      |      |      |      |      |
    |      |      |      |      |      |
    o------o------o------o------o------o
    |      | <--  | init |  --> |      | scaning from mid to left
    |      | <--  |  pos |  --> |      | scaning from mid to right
    o------o------o------o------o------o
    |      |      |      |      |      |
    |      |      |      |      |      |
    o------o------o------o------o------o
    |      |      |      |      |      |
    |      |      |      |      |      |
    o------o------o------o------o------o
    '''
    i = 0
    pkl = 'pkl'
    new_parameters = parameters.copy()
    new_parameters['is_right'] = True
    joblist = []
    joblist.append([i, pkl, new_parameters, parameters['system']])

    new_parameters2 =  new_parameters.copy()
    new_parameters2['is_right'] = False
    joblist.append([i, pkl, new_parameters2, parameters['system']])

    #pprint(joblist)

    p = multiprocessing.Pool(processes = parameters['NmaxThreads'])
    results = p.map(_run_left_right_line_in_parallel, joblist)
    #pprint(results)

    
    
    '''
    o------o------o------o------o------o
    |      |      |      |      |      |
    |      |      |      |      |      |
    o------o------o------o------o------o
    |  ^   |   ^  |  ^   |  ^   |  ^   |
    |  |   |   |  |  |   |  |   |  |   |
    o------o------o------o------o------o
    |  XX  |  XX  |  XX  |  XX  |  XX  | scanning from mid to north
    |  XX  |  XX  |  XX  |  XX  |  XX  | scanning from mid to south
    o------o------o------o------o------o
    |   |  |  |   |  |   |  |   |  |   |
    |   v  |  v   |  v   |  v   |  v   |
    o------o------o------o------o------o
    |      |      |      |      |      |
    |      |      |      |      |      |
    o------o------o------o------o------o
    '''
        
    joblist = []
    for result in results:
        for key in result.keys():
            #print(key)
            i = key[0]
            j = key[1]
            new_parameters = parameters.copy()
            new_parameters['j'] = j
            pkl = 'frame{}_{}.pkl'.format(i, j)
            joblist.append([i, pkl, new_parameters, parameters['system']])
    
    #'''
    p = multiprocessing.Pool(processes = parameters['NmaxThreads'])
    results2 = p.map(_run_up_down_column_in_parallel, joblist)

    #pprint(results2)

    data = {}
    for result in results:
        for key in result.keys():
            #print(key)
            i = key[0]
            j = key[1]
            matrix_i_j[i][j] = result[key][2]
            data[key] = result[key]
    for result in results2:
        for key in result.keys():
            #print(key)
            i = key[0]
            j = key[1]
            matrix_i_j[i][j] = result[key][2]
            data[key] = result[key]


    #pprint (matrix_i_j)
    #plot_data (matrix_i_j)
    

    #--------------------------------------------------------------------------------------
    #                             Writing the log information 
    #--------------------------------------------------------------------------------------
    for i in range( i_size ):
        for j in range(j_size):

            text = "\nDATA  %4i  %4i     %13.12f       %13.12f       %13.12f"% (int(i), int(j),  float(data[(i,j)][0]), float(data[(i,j)][1]), float(data[(i,j)][2]))
            arq.write(text)


def define_RC1_restraint (parameters, distance):
    """ Function doc """
    
    #distance = parameters['RC1']['dminimum'] + ( parameters['RC1']['dincre'] * float(i) )
    '''reaction coordinate 1 - starts at 0 and goes to nstpes'''
    '''----------------------------------------------------------------------------------------------------------------'''
    if parameters['RC1']["rc_type"] == 'simple_distance':
        #---------------------------------------------------------------------------------------------------------
        rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
        restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= parameters['RC1']['ATOMS'][0], point2= parameters['RC1']['ATOMS'][1])
        #restraints["RC1"] = restraint            
        #---------------------------------------------------------------------------------------------------------
    elif parameters['RC1']["rc_type"] == 'multiple_distance*4atoms':
        atom_RC1_3   = parameters['RC1']['ATOMS'][2]
        atom_RC1_4   = parameters['RC1']['ATOMS'][3]
        weight1 =  1.0#parameters['RC1']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
        weight2 = -1.0#parameters['RC1']['sigma_pk3pk1'] #self.sigma_a3_a1[0]
        rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
        restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ parameters['RC1']['ATOMS'][0], parameters['RC1']['ATOMS'][1], weight1 ], 
                                                                                                      [ parameters['RC1']['ATOMS'][2], parameters['RC1']['ATOMS'][3], weight2 ] 
                                                                                                      ] )
        #restraints["RC1"] = restraint  
    
    elif parameters['RC1']["rc_type"] == 'multiple_distance':
        #--------------------------------------------------------------------
        atom_RC1_3   = parameters['RC1']['ATOMS'][2]
        weight1 = parameters['RC1']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
        weight2 = parameters['RC1']['sigma_pk3pk1'] #self.sigma_a3_a1[0] 
        
        rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['RC1']['force_constant'])
        restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ parameters['RC1']['ATOMS'][1], parameters['RC1']['ATOMS'][0], weight1 ], 
                                                                                                      [ parameters['RC1']['ATOMS'][1], parameters['RC1']['ATOMS'][2], weight2 ] 
                                                                                                    ] )
        #restraints["RC1"] = restraint            
        #--------------------------------------------------------------------
    else:
        restraint = None
    
    return restraint


def define_RC2_restraint (parameters, distance2):

    if parameters['RC2']["rc_type"] == 'simple_distance':
        #---------------------------------------------------------------------------------------------------------
        rmodel            = RestraintEnergyModel.Harmonic(distance2, parameters['RC2']['force_constant'])
        restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= parameters['RC2']['ATOMS'][0], point2= parameters['RC2']['ATOMS'][1])
        #restraints["RC2"] = restraint            
        #---------------------------------------------------------------------------------------------------------
                
    elif parameters['RC2']["rc_type"] == 'multiple_distance*4atoms':
        atom_RC2_3   = parameters['RC2']['ATOMS'][2]
        atom_RC2_4   = parameters['RC2']['ATOMS'][3]
        weight1 =  1.0 #parameters['RC1']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
        weight2 = -1.0 #parameters['RC1']['sigma_pk3pk1'] #self.sigma_a3_a1[0]
        rmodel            = RestraintEnergyModel.Harmonic(distance2, parameters['RC2']['force_constant'])
        restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ [ parameters['RC2']['ATOMS'][0], parameters['RC2']['ATOMS'][1], weight1 ], 
                                                                                                      [ parameters['RC2']['ATOMS'][2], parameters['RC2']['ATOMS'][3], weight2 ] 
                                                                                                      ] )
        #restraints["RC2"] = restraint  



    elif parameters['RC2']["rc_type"] == 'multiple_distance':
        #--------------------------------------------------------------------
        atom_RC2_3 = parameters['RC2']['ATOMS'][2]
        weight1 = parameters['RC2']['sigma_pk1pk3'] #self.sigma_a1_a3[0]
        weight2 = parameters['RC2']['sigma_pk3pk1'] #self.sigma_a3_a1[0] 
        
        rmodel            = RestraintEnergyModel.Harmonic(distance2, parameters['RC2']['force_constant'])
        restraint         = RestraintMultipleDistance.WithOptions( energyModel = rmodel, distances= [ 
                                                                                                      [ parameters['RC2']['ATOMS'][1], parameters['RC2']['ATOMS'][0], weight1 ], 
                                                                                                      [ parameters['RC2']['ATOMS'][1], parameters['RC2']['ATOMS'][2], weight2 ] 
                                                                                                      ] )
         #--------------------------------------------------------------------
    else:
        restraint = None
    
    return restraint


def get_RC_d1_minus_d2 (parameters):
    """ Function doc """
    #--------------------------------------------------------------------------------------
    #                     Calculating the reaction coordinate 1
    #--------------------------------------------------------------------------------------
    distance1 = parameters['system'].coordinates3.Distance( parameters['RC1']['ATOMS'][0] , parameters['RC1']['ATOMS'][1]  )
    if parameters['RC1']["rc_type"] == 'multiple_distance':
        distance2 = parameters['system'].coordinates3.Distance( parameters['RC1']['ATOMS'][1],  parameters['RC1']['ATOMS'][2]  )

    elif parameters['RC1']["rc_type"] == 'multiple_distance*4atoms':
        distance2 = parameters['system'].coordinates3.Distance( parameters['RC1']['ATOMS'][2],  parameters['RC1']['ATOMS'][3]  )

    else:
        distance2 = 0
    
    RC1_d1_minus_d2 = distance1 - distance2
    
    
    #--------------------------------------------------------------------------------------
    #                     Calculating the reaction coordinate 2
    #--------------------------------------------------------------------------------------
    distance1 = parameters['system'].coordinates3.Distance( parameters['RC2']['ATOMS'][0] , parameters['RC2']['ATOMS'][1]  )
    
    if parameters['RC2']["rc_type"] == 'multiple_distance':
        distance2 = parameters['system'].coordinates3.Distance( parameters['RC2']['ATOMS'][1],  parameters['RC2']['ATOMS'][2]  )
    
    elif parameters['RC2']["rc_type"] == 'multiple_distance*4atoms':
        distance2 = parameters['system'].coordinates3.Distance( parameters['RC2']['ATOMS'][2],  parameters['RC2']['ATOMS'][3]  )
    
    else:
        distance2 = 0
    
    RC2_d1_minus_d2 = distance1 - distance2

    return RC1_d1_minus_d2, RC2_d1_minus_d2


def run_scan_energy (parameters):
    """ Function doc """
    

    try:
        #-------------------------------------------------------------------------------------------------------------
        ConjugateGradientMinimize_SystemGeometry(parameters['system']                                ,                
                                                 logFrequency           = parameters['logFrequency'],
                                                 maximumIterations      = parameters['maximumIterations'],
                                                 rmsGradientTolerance   = parameters['rmsGradientTolerance']  ,
                                                 log                    = None)
        energy   = parameters['system'].Energy(log=None)
        opt_convergency = ''
        #-------------------------------------------------------------------------------------------------------------
    except:
        energy   = 0
        opt_convergency = 'Ops! Error in geometry optimization.'
        #-------------------------------------------------------------------------------------------------------------
    
    return  energy


def _run_left_right_line_in_parallel (job):
    """ Function doc """

    i          = job[0]
    
    pkl        = job[1]
    parameters = job[2]
    system     = job[3]
    restraints = RestraintModel()
    parameters['system'].DefineRestraintModel( restraints )
    
    n_north = parameters['RC1']['up_nsteps'  ]
    n_south = parameters['RC1']['down_nsteps']
    
    n_right = parameters['RC2']['left_nsteps' ]
    n_left  = parameters['RC2']['right_nsteps']
                        

    matrix_index_i = n_south + i   
    
    data = {}
     
    if parameters['is_right']:
        for j in range(0, n_right):
            
            distance = parameters['RC1']['dminimum'] + ( parameters['RC1']['dincre'] * float(i) )
            restraints["RC1"] = define_RC1_restraint (parameters, distance)
            
            distance2 = parameters['RC2']['dminimum'] + ( parameters['RC2']['dincre'] * float(j) )
            restraints["RC2"] = define_RC2_restraint(parameters, distance2)
            

            energy   =  run_scan_energy (parameters)
            
            
            RC1_d1_minus_d2, RC2_d1_minus_d2 = get_RC_d1_minus_d2(parameters)
            
            matrix_index_j =  n_left + j 

            E = energy 
            
            #print (matrix_index_i, matrix_index_j, E )            
            data[(matrix_index_i, matrix_index_j)] = [RC1_d1_minus_d2, RC2_d1_minus_d2, E]
    


            pkl = os.path.join( parameters['folder'], 
                                parameters['traj_folder_name']+".ptGeo", 
                                "frame{}_{}.pkl".format(matrix_index_i, matrix_index_j) )
            
            Pickle( pkl, parameters['system'].coordinates3 )
            
            #backup_orca_files(system        = parameters['system'], 
            #                  output_folder = os.path.join(parameters['folder'],parameters['traj_folder_name']+".ptGeo") , 
            #                  output_name   = "frame{}_{}".format(i,j))
            
   
    #'''
    else:
        for j in range(n_left):
            
            distance = parameters['RC1']['dminimum'] + ( parameters['RC1']['dincre'] * float(i) )
            restraints["RC1"] = define_RC1_restraint (parameters, distance)
            
            distance2 = parameters['RC2']['dminimum'] + ( parameters['RC2']['dincre'] * float(-(j+1)) )
            restraints["RC2"] = define_RC2_restraint(parameters, distance2)
            
            matrix_index_j = n_left -j-1
            #RC1_d1_minus_d2 = parameters['dminimum_RC1'] + ( parameters['dincre_RC1'] * float(i) )     #r_i + r_i_incr*+i
            #RC2_d1_minus_d2 = parameters['dminimum_RC2'] + ( parameters['dincre_RC2'] * float(-(j+1)) )#r_j + r_j_incr*j
            
            energy   =  run_scan_energy (parameters)
            
            RC1_d1_minus_d2, RC2_d1_minus_d2 = get_RC_d1_minus_d2(parameters)
            
            E = energy 
            #print ( matrix_index_i, matrix_index_j, E )            
            data[(matrix_index_i, matrix_index_j)] = [RC1_d1_minus_d2, RC2_d1_minus_d2, E]

            pkl = os.path.join( parameters['folder'], 
                                parameters['traj_folder_name']+".ptGeo", 
                                "frame{}_{}.pkl".format(matrix_index_i, matrix_index_j) )
            
            Pickle( pkl, parameters['system'].coordinates3 )
            
    #'''
    return data


def _run_up_down_column_in_parallel (job):
    
    i          = job[0]
    pkl        = job[1]
    parameters = job[2]
    system     = job[3]

    n_north = parameters['RC1']['up_nsteps'  ]
    n_south = parameters['RC1']['down_nsteps']
    n_right = parameters['RC2']['left_nsteps' ]
    n_left  = parameters['RC2']['right_nsteps']
    
    restraints = RestraintModel()
    
    pkl = os.path.join( parameters['folder'], 
                        parameters['traj_folder_name']+".ptGeo",pkl) 
    
    parameters['system'].DefineRestraintModel( restraints )
    parameters['system'].coordinates3 = ImportCoordinates3(pkl )
    
    j = parameters ['j']

    matrix_index_i = n_south + i   
    
    data = {}
    
    #print (j, j -n_left)
    matrix_index_j = j
    
    j = j -n_left
    
    for i in range(0, n_north):
        matrix_index_i = n_south + i  
        
        
        #RC1_d1_minus_d2 = parameters['dminimum_RC1'] + ( parameters['dincre_RC1'] * float(i) ) 
        #RC2_d1_minus_d2 = parameters['dminimum_RC2'] + ( parameters['dincre_RC2'] * float((j)) )
        
        distance = parameters['RC1']['dminimum'] + ( parameters['RC1']['dincre'] * float(i) )
        restraints["RC1"] = define_RC1_restraint (parameters, distance)
        
        distance2 = parameters['RC2']['dminimum'] + ( parameters['RC2']['dincre'] * float(j) )
        restraints["RC2"] = define_RC2_restraint(parameters, distance2)
        
        energy   =  run_scan_energy (parameters)
        RC1_d1_minus_d2, RC2_d1_minus_d2 = get_RC_d1_minus_d2(parameters)
        E = energy 


        pkl = os.path.join( parameters['folder'], 
                            parameters['traj_folder_name']+".ptGeo", 
                            "frame{}_{}.pkl".format(matrix_index_i, matrix_index_j) )
        
        Pickle( pkl, parameters['system'].coordinates3 )


        #print (matrix_index_i, matrix_index_j, E )            
        data[(matrix_index_i, matrix_index_j)] = [RC1_d1_minus_d2, RC2_d1_minus_d2, E]





    parameters['system'].coordinates3 = ImportCoordinates3(pkl )
    for i in range(n_south):
        matrix_index_i = n_south -i -1
        
        #RC1_d1_minus_d2 = parameters['dminimum_RC1'] + ( parameters['dincre_RC1'] * float(-i-1) ) 
        #RC2_d1_minus_d2 = parameters['dminimum_RC2'] + ( parameters['dincre_RC2'] * float((j)) )
        
        distance = parameters['RC1']['dminimum'] + ( parameters['RC1']['dincre'] * float(-i-1) )
        restraints["RC1"] = define_RC1_restraint (parameters, distance)
        
        distance2 = parameters['RC2']['dminimum'] + ( parameters['RC2']['dincre'] * float(j) )
        restraints["RC2"] = define_RC2_restraint(parameters, distance2)
        
        #matrix_index_j = n_left -j-1
        energy         =  run_scan_energy (parameters)
        RC1_d1_minus_d2, RC2_d1_minus_d2 = get_RC_d1_minus_d2(parameters)
        
        E = energy 
        
        pkl = os.path.join( parameters['folder'], 
                            parameters['traj_folder_name']+".ptGeo", 
                            "frame{}_{}.pkl".format(matrix_index_i, matrix_index_j) )
        
        Pickle( pkl, parameters['system'].coordinates3 )

        #print ( matrix_index_i, matrix_index_j, E )            
        data[(matrix_index_i, matrix_index_j)] = [RC1_d1_minus_d2, RC2_d1_minus_d2, E]
    
    return data
