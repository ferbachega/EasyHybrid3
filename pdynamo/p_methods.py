#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  geometry_optmization.py
#  
#  Copyright 2022 Fernando Bachega <ferbachega@ufcpa.edu.br>
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


def backup_orca_files (system, output_folder = None, output_name = None):
    """ Function doc """
    
    
    if system.qcModel:
        items = system.qcModel.SummaryItems()
        '''Checking if ORCA is being used'''
        if items[0][0] == 'ORCA QC Model':
            
            scratch =  system.qcModel.scratch
            _time = time.asctime()
            
            
            if output_folder is None:
                folder = scratch
            else:
                folder = output_folder
            
            
            
            if output_name is None:
                output_name = 'orcaJob'+_time 
            else:
                pass
            
         
            
            print ('\nChecking for ORCA files at: ', scratch)

            try:
                os.rename(os.path.join(scratch, 'orcaJob.log'),   os.path.join(folder , output_name+'.orca.log'))
                print('Renaming logfile to:', os.path.join(folder , output_name+'.orca.log'))
            except:
                print('File {} not found at:'.format(output_name+'.orca.log'), scratch)
            
            try:
                os.rename(os.path.join(scratch, 'orcaJob.gbw'), os.path.join(folder , output_name+'.orca.gbw'))
                print('Renaming gbwfile to:', os.path.join(folder , output_name+'.orca.gbw'))
            except:
                print('File {} not found at:'.format(output_name+'.orca.gbw'), scratch)
            print ('\n')

            
        

HEADER ='''
----------------------------------------------------------------------------- 
                                EasyHybrid 3.0                                
                   - A pDynamo3 graphical user interface -                    
----------------------------------------------------------------------------- 

''' 
class LogFile:
    """ Class doc """
    
    def __init__ (self, system):
        """ Class initialiser """
        self.path = os.path.join(os.environ.get('PDYNAMO3_SCRATCH'), 'summary_temp.log')
        self.logFile2 = TextLogFileWriter.WithOptions ( path = self.path )
        system.Summary(log = self.logFile2)
        self.logFile2.Close()
        #return path

class EnergyCalculation:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass
    
    def run (self, parameters):
        """ Function doc """
        full_path_file = os.path.join(parameters['folder'])
        self.logFile2 = TextLogFileWriter.WithOptions ( path = os.path.join(full_path_file, parameters['filename']+'.log') )
        
        parameters['system'].Summary(log = self.logFile2)
        #try:
        energy = parameters['system'].Energy(log = self.logFile2)
        #except :
        #    msg = 'Error!'
        #    return False, msg
        
        backup_orca_files(system        = parameters['system']         , 
                          output_folder = parameters['folder']         , 
                          output_name   = parameters['filename'])
        
        self.logFile2.Footer ( )
        self.logFile2.Close()
        self.logFile2 = None
        
        return energy, 'Energy: '+str(energy)
        # - - - - - -


class EnergyRefinement:
    
    def __init__ (self):
        """ Class initialiser """
        pass
    
    def run (self, parameters, interface = False):
        
        full_path_trajectory = os.path.join(parameters['folder'], 
                               parameters['filename'])
        os.mkdir(
                 full_path_trajectory
                 )
        
        # - - - - - - - - - - - - - Checking trajectory - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        self.logFile2 = TextLogFileWriter.WithOptions ( path = os.path.join(full_path_trajectory, 'output.log') )
        parameters['system'].Summary(log = self.logFile2)
        self.logFile2.Header ( )
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        self.logFile2.Footer ( )
        self.logFile2.Close()
        self.logFile2 = None
        
        logfile = self.write_header (parameters = parameters,
                                     logfile    = os.path.join(full_path_trajectory, 'output.log') )
        
        if parameters['traj_type'] == 'pklfolder':
            data_energy = {}
            energy_list = []
            text = ""
            for frame in parameters["trajectory"]:
                x = frame[5:-4]
                parameters['system'].coordinates3 = ImportCoordinates3 (os.path.join(parameters['data_path'],frame)) 
                energy = parameters['system'].Energy()
                energy_list.append(energy)
                
                if parameters['RC1']['rc_type'] == 'simple_distance':
                    atom1 = parameters['RC1']['ATOMS'][0]
                    atom2 = parameters['RC1']['ATOMS'][1]
                    dist =  parameters['system'].coordinates3.Distance (atom1, atom2) 
                    data_energy[int(x)] = [energy, dist ]
                
                elif parameters['rc_type_1'] == 'multiple_distance':
                    atom1 = parameters['RC1']['ATOMS'][0]
                    atom2 = parameters['RC1']['ATOMS'][1]
                    atom3 = parameters['RC1']['ATOMS'][2]
                    dist1 =  parameters['system'].coordinates3.Distance (atom1, atom2) 
                    dist2 =  parameters['system'].coordinates3.Distance (atom2, atom3) 
                    dist  =  dist1 - dist2 
                    data_energy[int(x)] = [energy, dist1, dist2]
                
                print ('frame:', x, dist,  energy)
                backup_orca_files(system        = parameters['system']         , 
                                  output_folder = full_path_trajectory         , 
                                  output_name   = 'frame'+x)
            
            lowest_energy = min(energy_list)
            keys = data_energy.keys()
            for i in range(0, len(keys)):

                
                if parameters['RC1']['rc_type'] == 'simple_distance':
                    distance = data_energy[i][1]
                    energy   = data_energy[i][0]
                    text += "\nDATA %9i       %13.12f        %13.12f"% (int(i), float(distance), float(energy))
                
                elif parameters['RC1']['rc_type'] == 'multiple_distance':
                    distance1 = data_energy[i][1]
                    distance2 = data_energy[i][2]
                    energy    = data_energy[i][0]
                    text += "\nDATA %9i       %13.12f        %13.12f        %13.12f"% (int(i), float(distance1), float(distance2), float(energy))
                
                else:
                    pass
        
        
            logfile.write(text)
        
        
        elif parameters['traj_type'] == 'pklfolder2D': 
            max_i = 0
            max_j = 0
            i_list = []
            j_list = []
            data   = {}
            energy_list = []
            
            text = ""
            for frame in parameters["trajectory"]:
                i_j = frame[5:-4].split('_')
                parameters['system'].coordinates3 = ImportCoordinates3 (os.path.join(parameters['data_path'],frame)) 
                energy = parameters['system'].Energy()
                i = int(i_j[0])
                j = int(i_j[1])
                i_list.append(i)
                j_list.append(j)
                energy_list.append(energy)
                
                
                
                if parameters['RC1']['rc_type'] == 'simple_distance':
                    atom1 = parameters['RC1']['ATOMS'][0]
                    atom2 = parameters['RC1']['ATOMS'][1]
                    
                    dist1 =  parameters['system'].coordinates3.Distance (atom1, atom2) 
 
                
                elif parameters['RC1']['rc_type'] == 'multiple_distance':
                    atom1 = parameters['RC1']['ATOMS'][0]
                    atom2 = parameters['RC1']['ATOMS'][1]
                    atom3 = parameters['RC1']['ATOMS'][2]
                    dist1 =  parameters['system'].coordinates3.Distance (atom1, atom2) 
                    dist2 =  parameters['system'].coordinates3.Distance (atom2, atom3) 
                    
                    dist1 =  dist1 - dist2 
                 
                
                if parameters['RC2']['rc_type'] == 'simple_distance':
                    atom1 = parameters['RC2']['ATOMS'][0]
                    atom2 = parameters['RC2']['ATOMS'][1]
                    dist2 = parameters['system'].coordinates3.Distance (atom1, atom2) 
                 
                elif parameters['RC2']['rc_type'] == 'multiple_distance':
                    atom1 = parameters['RC2']['ATOMS'][0]
                    atom2 = parameters['RC2']['ATOMS'][1]
                    atom3 = parameters['RC2']['ATOMS'][2]
                    dist1 =  parameters['system'].coordinates3.Distance (atom1, atom2) 
                    dist2 =  parameters['system'].coordinates3.Distance (atom2, atom3) 
                    
                    dist2 = dist1 - dist2
                 
                
                data[(i,j)] = [dist1, dist2, energy]
                print ('frame:', i_j, energy)

            max_i = max(i_list)
            max_j = max(j_list)
            
            for i in range(0, max_i+1):
                for j in range(0, max_j+1):

                    text = "\nDATA  %4i  %4i     %13.12f       %13.12f       %13.12f"% (int(i), int(j),  float(data[(i,j)][0]), float(data[(i,j)][1]), float(data[(i,j)][2]))
                    logfile.write(text)

    def write_header (self, parameters, logfile = 'output.log'):
        """ Function doc """
        
        arq = open(logfile, "a")
        text = ""

        if parameters['RC2'] is not None :
            text = text + "\n"
            text = text + "\n--------------------------------------------------------------------------------"
            text = text + "\nTYPE                 EasyHybrid Energy Refinement 2D                            "
            text = text + "\n--------------------------------------------------------------------------------"

        else:
            text = text + "\n"
            text = text + "\n--------------------------------------------------------------------------------"
            text = text + "\nTYPE                   EasyHybrid Energy Refinement                             "
            text = text + "\n--------------------------------------------------------------------------------"




        '''This part writes the parameters used in the first reaction coordinate'''
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------
        if parameters['RC1']["rc_type"] == 'simple_distance':
            text = text + "\n"
            text = text + "\n----------------------- Coordinate 1 - Simple-Distance -------------------------"
            text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC1']['ATOMS'][0], parameters['RC1']['ATOM_NAMES'][0] )
            text = text + "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['RC1']['ATOMS'][1], parameters['RC1']['ATOM_NAMES'][1] )
            text = text + "\n--------------------------------------------------------------------------------"

        
        elif parameters['RC1']["rc_type"] == 'multiple_distance':
            text = text + "\n"
            text = text + "\n---------------------- Coordinate 1 - multiple-Distance ------------------------"	
            text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC1']['ATOMS'][0]    , parameters['RC1']['ATOM_NAMES'][0] )
            text = text + "\nATOM2*                 =%15i  ATOM NAME2             =%15s"     % (parameters['RC1']['ATOMS'][1]    , parameters['RC1']['ATOM_NAMES'][1] )
            text = text + "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC1']['ATOMS'][2]    , parameters['RC1']['ATOM_NAMES'][2] )
            text = text + "\n--------------------------------------------------------------------------------"
        else:
            pass
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------
        
        
        
        
        if parameters['RC2'] is not None :

            '''This part writes the parameters used in the second reaction coordinate'''
            #-----------------------------------------------------------------------------------------------------------------------------------------------------------
            if parameters['RC2']["rc_type"] == 'simple_distance':
                text = text + "\n"
                text = text + "\n----------------------- Coordinate 2 - Simple-Distance -------------------------"
                text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC2']['ATOMS'][0], parameters['RC2']['ATOM_NAMES'][0] )
                text = text + "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['RC2']['ATOMS'][1], parameters['RC2']['ATOM_NAMES'][1] )
                text = text + "\n--------------------------------------------------------------------------------"

            
            elif parameters['RC2']["rc_type"] == 'multiple_distance':
                text = text + "\n"
                text = text + "\n---------------------- Coordinate 2 - multiple-Distance ------------------------"	
                text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC2']['ATOMS'][0]    , parameters['RC2']['ATOM_NAMES'][0] )
                text = text + "\nATOM2*                 =%15i  ATOM NAME2             =%15s"     % (parameters['RC2']['ATOMS'][1]    , parameters['RC2']['ATOM_NAMES'][1] )
                text = text + "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC2']['ATOMS'][2]    , parameters['RC2']['ATOM_NAMES'][2] )
                text = text + "\n--------------------------------------------------------------------------------"
            else:
                pass
            #-----------------------------------------------------------------------------------------------------------------------------------------------------------
            
        
        
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------
        '''This part writes the header of the frames distances and energy'''
        if parameters['RC2'] is not None :
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
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------


        
        arq.write(text)
        return arq

    







class GeometryOptimization:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass
        self.trajectory = None
        self.logFile2   = None
        
    def run (self, parameters):
        """ Function doc """
        
        # - - - - - - - - - - - - - Checking trajectory - - - - - - - - - - - - - -
        '''If a trajectory is provided, it holds for all optimization functions'''
        if parameters['trajectory_name']:
            full_path_trajectory = os.path.join(parameters['folder'], parameters['trajectory_name'])
            self.trajectory = ExportTrajectory(full_path_trajectory, parameters['system'], log=None )
            
            self.logFile2 = TextLogFileWriter.WithOptions ( path = os.path.join(full_path_trajectory, 'output.log') )
            parameters['system'].Summary(log = self.logFile2)
            self.logFile2.Header ( )
        else:
            pass
            parameters['trajectory_name'] = None
            parameters['folder'] = None
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        
        if parameters['optimizer'] == 'ConjugatedGradient':
            self._run_conjugated_gradient(parameters)
        
        elif parameters['optimizer'] == 'SteepestDescent':
            self._run_steepest_descent(parameters)
        
        elif parameters['optimizer'] == 'LFBGS':
            self._run_LBFGS(parameters)
        
        elif parameters['optimizer'] == 'QuasiNewton':
            self._run_quasi_newton(parameters)
        
        elif parameters['optimizer'] == 'FIRE':
            self._run_FIRE(parameters)
        
        else:
            print('Geometry Optimization method not found!' )
        
        
        '''ORCA backup files'''
        backup_orca_files(system = parameters['system'], output_folder = parameters['folder'], output_name = parameters['trajectory_name'])
        
        
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if parameters['trajectory_name']:
            self.logFile2.Footer ( )
            self.logFile2.Close()
            self.logFile2 = None
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    def _run_conjugated_gradient(self, parameters):
        """ Function doc """
        
        
        
        if parameters['trajectory_name'] == None:
            ConjugateGradientMinimize_SystemGeometry(parameters['system']                                ,                
                                                     logFrequency           = parameters['log_frequency'],
                                                     maximumIterations      = parameters['maxIterations'],
                                                     rmsGradientTolerance   = parameters['rmsGradient'])
        else:            
            
            ConjugateGradientMinimize_SystemGeometry(parameters['system']                                                ,                
                                                     logFrequency           = parameters['log_frequency']                ,
                                                     trajectories           = [(self.trajectory, parameters['save_frequency'])],
                                                     maximumIterations      = parameters['maxIterations']                ,
                                                     rmsGradientTolerance   = parameters['rmsGradient']                  ,
                                                     log                    = self.logFile2                                   )

    def _run_steepest_descent(self, parameters):

        if parameters['trajectory_name'] == None:
            SteepestDescentMinimize_SystemGeometry(parameters['system']                               ,               
                                                logFrequency            = parameters['log_frequency'] ,
                                                maximumIterations       = parameters['maxIterations'] ,
                                                rmsGradientTolerance    = parameters['rmsGradient']   )
        else:
            SteepestDescentMinimize_SystemGeometry(parameters['system']                                               ,               
                                                logFrequency            = parameters['log_frequency']                 ,
                                                trajectories            = [(self.trajectory, parameters['save_frequency'])],
                                                maximumIterations       = parameters['maxIterations']                 ,
                                                rmsGradientTolerance    = parameters['rmsGradient']                   ,
                                                log                     = self.logFile2                                   )

    def _run_LBFGS(self, parameters):
        
        if parameters['trajectory_name'] == None:
            LBFGSMinimize_SystemGeometry(parameters['system']                               ,                
                                         logFrequency         = parameters['log_frequency'] ,
                                         maximumIterations    = parameters['maxIterations'] ,
                                         rmsGradientTolerance = parameters['rmsGradient']   )
        else:

            LBFGSMinimize_SystemGeometry(parameters['system']                                                ,                
                                         logFrequency         = parameters['log_frequency']                  ,
                                         trajectories         = [(self.trajectory, parameters['save_frequency'])] ,
                                         maximumIterations    = parameters['maxIterations']                  ,
                                         rmsGradientTolerance = parameters['rmsGradient']                    ,
                                         log                  = self.logFile2                                   )  

    def _run_quasi_newton(self, parameters):
        '''
        Class method to apply the Quaisi-Newton minimizer
        '''        
        if parameters['trajectory_name'] == None:
            QuasiNewtonMinimize_SystemGeometry( parameters['system']                               ,               
                                                logFrequency         = parameters['log_frequency'] ,
                                                maximumIterations    = parameters['maxIterations'] ,
                                                rmsGradientTolerance = parameters['rmsGradient']   )
        else:
            QuasiNewtonMinimize_SystemGeometry( parameters['system']                                                ,
                                                logFrequency         = parameters['log_frequency']                  ,
                                                trajectories         = [(self.trajectory, parameters['save_frequency'])] ,
                                                maximumIterations    = parameters['maxIterations']                  ,
                                                rmsGradientTolerance = parameters['rmsGradient']                    ,
                                                log                  = self.logFile2                                   )
    
    def _run_FIRE(self, parameters):
        '''
        '''
        if parameters['trajectory_name'] == None:
            FIREMinimize_SystemGeometry( parameters['system']                               ,          
                                         logFrequency         = parameters['log_frequency'] ,
                                         maximumIterations    = parameters['maxIterations'] ,
                                         rmsGradientTolerance = parameters['rmsGradient']   )
        else:

            FIREMinimize_SystemGeometry( parameters['system']                                                ,                
                                         logFrequency         = parameters['log_frequency']                  ,
                                         trajectories         = [(self.trajectory, parameters['save_frequency'])] ,
                                         maximumIterations    = parameters['maxIterations']                  ,
                                         rmsGradientTolerance = parameters['rmsGradient']                    ,
                                         log                  = self.logFile2                                   )



class MolecularDynamics:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass

    def run (self, parameters = None):
        """ Function doc """
        
        # - - - - - - - - - - - - - Checking trajectory - - - - - - - - - - - - - -
        '''If a trajectory is provided, it holds for all optimization functions'''
        if parameters['trajectory_name']:
            full_path_trajectory = os.path.join(parameters['folder'], parameters['trajectory_name'] + '.ptGeo')
            #print('\n\n\n')
            #print('full_path_trajectory:', full_path_trajectory)
            self.trajectory = ExportTrajectory(full_path_trajectory, parameters['system'], log=None )
            
            self.logFile2 = TextLogFileWriter.WithOptions ( path = os.path.join(full_path_trajectory, 'output.log') )
            parameters['system'].Summary(log = self.logFile2)
            self.logFile2.Header ( )
        else:
            pass
            parameters['trajectory_name'] = None
            parameters['folder'] = None
        
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        
        if parameters['integrator'] == 'Verlet':
            self._velocity_verlet_dynamics (parameters)
        
        elif parameters['integrator'] == 'LeapFrog':
            self._leap_frog_dynamics(parameters)
        
        elif parameters['integrator'] == 'Langevin':
            self._langevin_dynamics(parameters)
        
        else:
            pass    
        
        backup_orca_files(system = parameters['system'], output_folder = parameters['folder'], output_name = parameters['trajectory_name'])
        
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if parameters['trajectory_name']:
            self.logFile2.Footer ( )
            self.logFile2.Close()
            self.logFile2 = None
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


    
    def _velocity_verlet_dynamics (self, parameters = None):
        """ Function doc """
        # . Define a normal deviate generator in a given state.
        normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( parameters['seed'] ) )        
    
        
        if parameters['temperatureScaleOption'] == 'constant':
            # . Equilibration. parameters['save_frequency'])]
            
            if parameters['trajectory_name'  ]:               
                
                #logFile2 = TextLogFileWriter.WithOptions ( path = os.path.join(full_path_trajectory, 'output.log') )
                #parameters['system'].Summary(log = logFile2)
                #logFile2.Header ( )
                
                VelocityVerletDynamics_SystemGeometry ( parameters['system']                                                    ,
                                                        logFrequency              = parameters['logFrequency']                  ,
                                                        normalDeviateGenerator    = normalDeviateGenerator                      ,
                                                        steps                     = parameters['steps']                         ,
                                                        timeStep                  = parameters['timeStep']                      ,
                                                        temperatureScaleFrequency = parameters['temperatureScaleFrequency']     ,
                                                        temperatureScaleOption    = parameters['temperatureScaleOption']        ,
                                                        temperatureStart          = parameters['temperatureStart']              ,
                                                        trajectories              = [(self.trajectory, parameters['trajectory_frequency'])],
                                                        log                       = self.logFile2
                                                        )
                #logFile2.Footer ( )
                #logFile2.Close()
            
            else:
                VelocityVerletDynamics_SystemGeometry ( parameters['system']                                               ,
                                                        logFrequency              = parameters['logFrequency']             ,
                                                        normalDeviateGenerator    = normalDeviateGenerator                 ,
                                                        steps                     = parameters['steps']                    ,
                                                        timeStep                  = parameters['timeStep']                 ,
                                                        temperatureScaleFrequency = parameters['temperatureScaleFrequency'],
                                                        temperatureScaleOption    = parameters['temperatureScaleOption']   ,
                                                        temperatureStart          = parameters['temperatureStart']         )


        else:
            # . Heating.
            if parameters['trajectory_name'  ]:


                
                VelocityVerletDynamics_SystemGeometry ( parameters['system']                                               ,
                                                        logFrequency              = parameters['logFrequency']             ,
                                                        normalDeviateGenerator    = normalDeviateGenerator                 ,
                                                        steps                     = parameters['steps']                    ,
                                                        timeStep                  = parameters['timeStep']                 ,
                                                        temperatureScaleFrequency = parameters['temperatureScaleFrequency'],
                                                        temperatureScaleOption    = parameters['temperatureScaleOption']   ,
                                                        temperatureStart          = parameters['temperatureStart']         ,
                                                        temperatureStop           = parameters['temperatureStop']          ,
                                                        trajectories              = [(self.trajectory, parameters['trajectory_frequency'])],
                                                        log                       = self.logFile2
                                                        )                                                    
                


            else:
                VelocityVerletDynamics_SystemGeometry ( parameters['system']                                               ,
                                                        logFrequency              = parameters['logFrequency']             ,
                                                        normalDeviateGenerator    = normalDeviateGenerator                 ,
                                                        steps                     = parameters['steps']                    ,
                                                        timeStep                  = parameters['timeStep']                 ,
                                                        temperatureScaleFrequency = parameters['temperatureScaleFrequency'],
                                                        temperatureScaleOption    = parameters['temperatureScaleOption']   ,
                                                        temperatureStart          = parameters['temperatureStart']         ,
                                                        temperatureStop           = parameters['temperatureStop']          ,
                                                         )



    def _leap_frog_dynamics (self, parameters = None):
        """ Function doc """
        # . Define a normal deviate generator in a given state.
        normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( parameters['seed'] ) )        
        '''
        #normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( 491831 ) )
        if parameters['trajectory_name'  ]:
            full_path_trajectory = os.path.join(parameters['folder'], parameters['trajectory_name'] +".ptGeo")
            print(full_path_trajectory)
            trajectory = ExportTrajectory(full_path_trajectory, parameters['system'], log=None )        
        # . Data-collection.
        #'''
        if parameters['pressureControl']:
            if parameters['trajectory_name']:
                LeapFrogDynamics_SystemGeometry ( parameters['system']                       ,
                                                  logFrequency           = parameters['logFrequency'] ,
                                                  normalDeviateGenerator = normalDeviateGenerator ,
                                                  pressure               = parameters['pressure'],
                                                  pressureControl        = parameters['pressureControl'],
                                                  pressureCoupling       = parameters['pressureCoupling'],
                                                  steps                  = parameters['steps'],
                                                  temperature            = parameters['temperatureStart'],
                                                  temperatureControl     = parameters['temperatureControl'],
                                                  temperatureCoupling    = parameters['temperatureCoupling'],
                                                  timeStep               = parameters['timeStep'],
                                                  trajectories           = [(self.trajectory, parameters['trajectory_frequency'])],
                                                  log                    = self.logFile2
                                                  )

            else:
                # . Equilibration.
                LeapFrogDynamics_SystemGeometry ( parameters['system']                       ,
                                                  logFrequency           = parameters['logFrequency'] ,
                                                  normalDeviateGenerator = normalDeviateGenerator ,
                                                  pressure               = parameters['pressure'],
                                                  pressureControl        = parameters['pressureControl'],
                                                  pressureCoupling       = parameters['pressureCoupling'],
                                                  steps                  = parameters['steps'],
                                                  temperature            = parameters['temperatureStart'],
                                                  temperatureControl     = parameters['temperatureControl'],
                                                  temperatureCoupling    = parameters['temperatureCoupling'],
                                                  timeStep               = parameters['timeStep'])
        
        else:
            if parameters['trajectory_name']:
                LeapFrogDynamics_SystemGeometry ( parameters['system']                                                            ,
                                                  logFrequency           = parameters['logFrequency']                             ,
                                                  normalDeviateGenerator = normalDeviateGenerator                                 ,
                                                  steps                  = parameters['steps']                                    ,
                                                  temperature            = parameters['temperatureStart']                         ,
                                                  temperatureControl     = parameters['temperatureControl']                       ,
                                                  temperatureCoupling    = parameters['temperatureCoupling']                      ,
                                                  timeStep               = parameters['timeStep']                                 ,
                                                  trajectories           = [(self.trajectory, parameters['trajectory_frequency'])],
                                                  log                    = self.logFile2
                                                  )
            else:
                LeapFrogDynamics_SystemGeometry ( parameters['system']                                          ,
                                                  logFrequency           = parameters['logFrequency']           ,
                                                  normalDeviateGenerator = normalDeviateGenerator               ,
                                                  steps                  = parameters['steps']                  ,
                                                  temperature            = parameters['temperatureStart']       ,
                                                  temperatureControl     = parameters['temperatureControl']     ,
                                                  temperatureCoupling    = parameters['temperatureCoupling']    ,
                                                  timeStep               = parameters['timeStep']               )


    def _langevin_dynamics (self, parameters = None):
        
        """ Function doc """
        
        # . Define a normal deviate generator in a given state.
        normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( parameters['seed'] ) )        
        '''
        if parameters['trajectory_name'  ]:
            full_path_trajectory = os.path.join(parameters['folder'], parameters['trajectory_name'] +".ptGeo")
            print(full_path_trajectory)
            trajectory = ExportTrajectory(full_path_trajectory, parameters['system'], log=None )        
        '''
        
        if parameters['trajectory_name']:
            LangevinDynamics_SystemGeometry ( parameters['system']                                                            ,
                                              collisionFrequency     = parameters['collisionFrequency']                       ,
                                              logFrequency           = parameters['logFrequency']                             ,
                                              normalDeviateGenerator = normalDeviateGenerator                                 ,
                                              steps                  = parameters['steps']                                    ,
                                              temperature            = parameters['temperatureStart']                         ,
                                              timeStep               = parameters['timeStep']                                 ,
                                              trajectories           = [(self.trajectory, parameters['trajectory_frequency'])],
                                              log                    = self.logFile2
                                              )

        else:
            # . Dynamics.
            LangevinDynamics_SystemGeometry ( parameters['system']                                      ,
                                              collisionFrequency     = parameters['collisionFrequency'] ,
                                              logFrequency           = parameters['logFrequency']       ,
                                              normalDeviateGenerator = normalDeviateGenerator           ,
                                              steps                  = parameters['steps']              ,
                                              temperature            = parameters['temperatureStart']   ,
                                              timeStep               = parameters['timeStep']           )
            

            
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


        
        if parameters['RC2'] is not None:
            print('\n\n\nself._run_scan_2D(parameters = parameters, interface = False)')
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
            #text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['dminimum_RC1'], parameters['maxIterations']      )
            #text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['dincre_RC1']  , parameters['rmsGradient']        )
            #text = text + "\n--------------------------------------------------------------------------------"
            text = text + "\n"
            text = text + "\n----------------------- Coordinate 1 - Simple-Distance -------------------------"
            text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC1']['ATOMS'][0], parameters['RC1']['ATOM_NAMES'][0] )
            text = text + "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['RC1']['ATOMS'][1], parameters['RC1']['ATOM_NAMES'][1] )
            text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['RC1']['nsteps']  , parameters['RC1']['force_constant'] )
            text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC1']['dminimum'], parameters['maxIterations']         )
            text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC1']['dincre']  , parameters['rmsGradient']           )
            text = text + "\n--------------------------------------------------------------------------------"

        
        elif parameters['RC1']["rc_type"] == 'multiple_distance':
            text = text + "\n"
            text = text + "\n---------------------- Coordinate 1 - multiple-Distance ------------------------"	
            text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC1']['ATOMS'][0]    , parameters['RC1']['ATOM_NAMES'][0] )
            text = text + "\nATOM2*                 =%15i  ATOM NAME2             =%15s"     % (parameters['RC1']['ATOMS'][1]    , parameters['RC1']['ATOM_NAMES'][1] )
            text = text + "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC1']['ATOMS'][2]    , parameters['RC1']['ATOM_NAMES'][2] )
            text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['RC1']['nsteps']      , parameters['RC1']['force_constant'] ) 
            text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC1']['dminimum']    , parameters['maxIterations']         )
            text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC1']['dincre']      , parameters['rmsGradient']           )
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
                text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC2']['dminimum'], parameters['maxIterations']         )
                text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC2']['dincre']  , parameters['rmsGradient']           )
                text = text + "\n--------------------------------------------------------------------------------"

            
            elif parameters['RC2']["rc_type"] == 'multiple_distance':
                text = text + "\n"
                text = text + "\n---------------------- Coordinate 2 - multiple-Distance ------------------------"	
                text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC2']['ATOMS'][0]    , parameters['RC2']['ATOM_NAMES'][0] )
                text = text + "\nATOM2*                 =%15i  ATOM NAME2             =%15s"     % (parameters['RC2']['ATOMS'][1]    , parameters['RC2']['ATOM_NAMES'][1] )
                text = text + "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC2']['ATOMS'][2]    , parameters['RC2']['ATOM_NAMES'][2] )
                text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['RC2']['nsteps']      , parameters['RC2']['force_constant'] )
                text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC2']['dminimum']    , parameters['maxIterations']         )
                text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC2']['dincre']      , parameters['rmsGradient']           )
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
                                                         logFrequency           = parameters['log_frequency'],
                                                         maximumIterations      = parameters['maxIterations'],
                                                         rmsGradientTolerance   = parameters['rmsGradient'])
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
                                                         logFrequency           = parameters['log_frequency'],
                                                         maximumIterations      = parameters['maxIterations'],
                                                         rmsGradientTolerance   = parameters['rmsGradient'])
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
                                                         logFrequency           = parameters['log_frequency'],
                                                         maximumIterations      = parameters['maxIterations'],
                                                         rmsGradientTolerance   = parameters['rmsGradient'])
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
        backup.append(parameters['system'].e_treeview_iter)
        backup.append(parameters['system'].e_liststore_iter)
        parameters['system'].e_treeview_iter   = None
        parameters['system'].e_liststore_iter  = None
        
       
        
        
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
                                                         logFrequency           = parameters['log_frequency'],
                                                         maximumIterations      = parameters['maxIterations'],
                                                         rmsGradientTolerance   = parameters['rmsGradient']  ,
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
                             
                             'log_frequency'   : parameters['log_frequency'],
                             'maxIterations'   : parameters['maxIterations'],
                             'rmsGradient'     : parameters['rmsGradient'],

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
        parameters['system'].e_treeview_iter   = backup[0]
        parameters['system'].e_liststore_iter  = backup[1]

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
                                                         logFrequency           = parameters['log_frequency'],
                                                         maximumIterations      = parameters['maxIterations'],
                                                         rmsGradientTolerance   = parameters['rmsGradient'])
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
            text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC1']['dminimum'], parameters['OPT_parm']['maxIterations']         )
            text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC1']['dincre']  , parameters['OPT_parm']['rmsGradient']           )
            text = text + "\n--------------------------------------------------------------------------------"

        
        elif parameters['RC1']["rc_type"] == 'multiple_distance':
            text = text + "\n"
            text = text + "\n---------------------- Coordinate 1 - multiple-Distance ------------------------"	
            text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC1']['ATOMS'][0]    , parameters['RC1']['ATOM_NAMES'][0] )
            text = text + "\nATOM2*                 =%15i  ATOM NAME2             =%15s"     % (parameters['RC1']['ATOMS'][1]    , parameters['RC1']['ATOM_NAMES'][1] )
            text = text + "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC1']['ATOMS'][2]    , parameters['RC1']['ATOM_NAMES'][2] )
            text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['RC1']['nsteps']      , parameters['RC1']['force_constant'] ) 
            text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC1']['dminimum']    , parameters['OPT_parm']['maxIterations']         )
            text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC1']['dincre']      , parameters['OPT_parm']['rmsGradient']           )
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
        #        text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC2']['dminimum'], parameters['maxIterations']         )
        #        text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC2']['dincre']  , parameters['rmsGradient']           )
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
        #        text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['RC2']['dminimum']    , parameters['maxIterations']         )
        #        text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['RC2']['dincre']      , parameters['rmsGradient']           )
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

        pprint (parameters)
        
        #-------------------------------------------------------------------------
        if interface:
            e_treeview_iter  = getattr(parameters['system'], 'e_treeview_iter', None)
            e_liststore_iter = getattr(parameters['system'], 'e_liststore_iter', None)
            parameters['system'].e_treeview_iter   = None
            parameters['system'].e_liststore_iter  = None            
        #-------------------------------------------------------------------------
        
        if parameters['RC2'] is not None:
            self._run_umbrella_sampling_2D(parameters = parameters, interface = False)
        else:
            self._run_umbrella_sampling_1D(parameters)
        
        #-------------------------------------------------------------------------
        if interface:
            parameters['system'].e_treeview_iter  = e_treeview_iter 
            parameters['system'].e_liststore_iter = e_liststore_iter 
        #-------------------------------------------------------------------------


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
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            
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
            results = p.map(_run_parallel_umbrella_sampling_1D, joblist)
            #-------------------------------------------------------------------------


            #parameters['system'].e_treeview_iter  = e_treeview_iter 
            #parameters['system'].e_liststore_iter = e_liststore_iter 
            
        
        else:
            self._run_serial_umbrella_sampling_1D (parameters)
    

    def _run_serial_umbrella_sampling_1D (self, parameters):
        """ Function doc """
        
        #-------------------------------------------------------------------------
        #Setting some local vars to ease the notation in the pDynamo methods
        #----------------------------------
        opt_parameters  = parameters['OPT_parm']
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
            else:
                pass
                
            
            
            '''             G E O M E T R Y   O P T I M I Z A T I O N            '''
            if parameters['OPT_parm'] is not None:
                _us_geo_opt(parameters['system'], opt_parameters)
            
            
            
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
            results = p.map(_run_parallel_umbrella_sampling_2D, joblist)
            #-------------------------------------------------------------------------
            #'''



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
    if hamiltonian in ['DFTB QC Model', 'ORCA QC Model', 'external']:
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
   
    '''             G E O M E T R Y   O P T I M I Z A T I O N            '''
    if parameters['OPT_parm'] is not None:
        _us_geo_opt(parameters['system'], opt_parameters)



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
    if hamiltonian in ['DFTB QC Model', 'ORCA QC Model', 'external']:
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
    else:
        pass




    '''             G E O M E T R Y   O P T I M I Z A T I O N            '''
    if parameters['OPT_parm'] is not None:
        _us_geo_opt(parameters['system'], opt_parameters)



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
    parameters['log_frequency'] = 50
    
    if parameters['optimizer'] == 'ConjugatedGradient':
        #-------------------------------------------------------------------------------------------------------------
        ConjugateGradientMinimize_SystemGeometry(system                     ,                
                                                 logFrequency           = parameters['log_frequency'],
                                                 maximumIterations      = parameters['maxIterations'],
                                                 rmsGradientTolerance   = parameters['rmsGradient'])
        #-------------------------------------------------------------------------------------------------------------
    
    elif parameters['optimizer'] == 'SteepestDescent':
        SteepestDescentMinimize_SystemGeometry( system                                            ,               
                                            logFrequency            = parameters['log_frequency'] ,
                                            maximumIterations       = parameters['maxIterations'] ,
                                            rmsGradientTolerance    = parameters['rmsGradient']   )
    
    elif parameters['optimizer'] == 'LFBGS':
        LBFGSMinimize_SystemGeometry(parameters['system']                               ,                
                                     logFrequency         = parameters['log_frequency'] ,
                                     maximumIterations    = parameters['maxIterations'] ,
                                     rmsGradientTolerance = parameters['rmsGradient']   )
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
    #-----------------------------------------------------------------------------------------------------
    
    
    
    
    #-----------------------------------------------------------------------------------------------------
    # . Data Collection
    elif mode == 1:
        if j == None:
            full_path_trajectory = os.path.join(path_trajectory, 'window{:d}.ptRes'.format ( i ))
        else:
            full_path_trajectory = os.path.join(path_trajectory, 'window{:d}_{:d}.ptRes'.format ( i, j ))
        parameters['steps']  = parameters['steps_dc']
    #-----------------------------------------------------------------------------------------------------
    else:
        pass
    
    
    #trajectoy instances
    trajectory = ExportTrajectory(full_path_trajectory,  system , log=None )
    
    
    
    
        
    if parameters['integrator'] == 'Verlet':
        _us_velocity_verlet_dynamics(system, trajectory, parameters)
    
    elif parameters['integrator'] == 'LeapFrog':
        _us_leap_frog_dynamics (system, trajectory, parameters)
    
    elif parameters['integrator'] == 'Langevin':
        _us_langevin_dynamics (system, trajectory, parameters)
    
    else:
        pass
        
        
def _us_leap_frog_dynamics (system, trajectory, parameters ):
    """ Function doc """
    if parameters['pressureControl']:
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
                                          trajectories           = [(trajectory, parameters['trajectory_frequency'])],
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
                                          trajectories           = [(trajectory,  parameters['trajectory_frequency'])]
                                          
                                          )
                                          

def _us_velocity_verlet_dynamics (system, trajectory, parameters):
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
                                            trajectories              = [(trajectory, parameters['trajectory_frequency'])],
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
                                      trajectories           = [(trajectory, parameters['trajectory_frequency'])],
                                      #log                    = self.logFile2
                                      )



class WHAMAnalysis:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass

    def run (self, parameters, interface = False):
        """ Function doc """
        
        try:
            self.logFile2 = TextLogFileWriter.WithOptions ( path = os.path.join(parameters['folder'], parameters['logfile']+'.log') )
        except TypeError:
            return False, 'Error generating log file! \n\nPlease, check your working folder and output logfile name'
        
        
        
        # . Calculate the PMF.
        if parameters['system']:
            parameters['system'].Summary(log = self.logFile2)
        else:
            pass

        
        
        try:
            state = WHAM_ConjugateGradientMinimize ( parameters['fileNames']                                  ,
                                                     bins                 = parameters['bins'                 ],
                                                     logFrequency         = parameters['logFrequency'         ],
                                                     maximumIterations    = parameters['maximumIterations'    ],
                                                     rmsGradientTolerance = parameters['rmsGradientTolerance' ],
                                                     temperature          = parameters['temperature'          ],
                                                     
                                                     log = self.logFile2
                                                     )
        except TypeError:
            return False, 'Error on pDynamo WHAM_ConjugateGradientMinimize! \n\nPlease, check your input data.'
        
        
        
        try:
            # . Write the PMF to a file.
            histogram = state["Histogram"]
            pmf       = state["PMF"      ]
            
            PMF_file = os.path.join(parameters['folder'], parameters['logfile']+'_pmf.log')
            histogram.ToTextFileWithData (  PMF_file , [ pmf ], format = "{:20.3f} {:20.3f}\n" )
        
        except TypeError:
            return False, 'Error on Write the PMF to a file!\n\n Please, check your input data.'

        msg = 'WHAM calculation performed successfully! \n\nPlease check:'
        return True, msg+PMF_file


'''
try:
    num1 = int(input("Enter a number: "))
    num2 = int(input("Enter another number: "))
    result = num1 / num2
    print("Result:", result)
except ValueError:
    print("Invalid input. Please enter a valid integer.")
except ZeroDivisionError:
    print("Cannot divide by zero.")
'''














class ChainOfStatesOptimizePath:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass                
    
    def run (self, parameters):         
          
          
        '''-----------------------------------------------------------------------------------------------------'''        
        trajectoryPath = os.path.join ( parameters['folder'], parameters['trajectory_name'] +".ptGeo" )
        self.trajectory = ExportTrajectory(trajectoryPath, parameters['system'], log=None )

        self.logFile2 = TextLogFileWriter.WithOptions ( path = os.path.join(trajectoryPath, 'output.log') )
        parameters['system'].Summary(log = self.logFile2)
        self.logFile2.Header ( )
        '''-----------------------------------------------------------------------------------------------------'''        
        #logFile.Header ( )
        system = parameters['system']

        # . Assign the reactant and product coordinates.
        reactants = parameters['reac_coordinates']#ImportCoordinates3 ( os.path.join ( xyzPath, "cyclohexane_chair.xyz"     ) )
        products  = parameters['prod_coordinates']#ImportCoordinates3 ( os.path.join ( xyzPath, "cyclohexane_twistboat.xyz" ) )

        # . Get an initial path.

        GrowingStringInitialPath ( system                         , 
                                   parameters['number_of_structures'], 
                                   reactants                      , 
                                   products                       , 
                                   trajectoryPath                 ,
                                   log  = self.logFile2)

        # . Optimization.
        trajectory = ExportTrajectory ( trajectoryPath, system, append = True )
        
        ChainOfStatesOptimizePath_SystemGeometry ( system                                                   ,
                                                   trajectory                                               ,
                                                   logFrequency         = parameters['log_frequency']       ,
                                                   maximumIterations    = parameters['maximumIterations']   ,
                                                   rmsGradientTolerance = parameters['rmsGradientTolerance'],
                                                   
                                                   springForceConstant           = parameters['spring_force_constant']          ,
                                                   splineRedistributionTolerance = parameters['spline_redistribution_tolerance'],
                                                   
                                                   forceSplineRedistributionCheckPerIteration = parameters['force_spline_redistribution_check_per_iteration'],
                                                   fixedTerminalImages                        = parameters['fixed_terminal_images']                          ,
                                                   
                                                   
                                                   log                                        = self.logFile2                     )



        # . Footer.
        self.logFile2.Footer ( )
        self.logFile2.Close ( )
        self.logFile2 = None



class NormalModes:
    def __init__ (self):
        """ Class initialiser """
        pass                
    
    def run (self, parameters):         
        
        trajectoryPath  = os.path.join ( parameters['folder'], parameters['trajectory_name'])                            
        
        #
        isExist = os.path.exists(trajectoryPath)
        if isExist:
            pass
        else:
            os.mkdir(trajectoryPath)
        #
        
        
        # . Logfile
        self.logFile2 = TextLogFileWriter.WithOptions ( path = os.path.join(trajectoryPath, parameters['trajectory_name']+'.log') )
        parameters['system'].Summary(log = self.logFile2)
        self.logFile2.Header ( )
        
        # . Calculate the normal modes.
        NormalModes_SystemGeometry ( parameters['system'], modify = ModifyOption.Project, log = self.logFile2 )
        
        # . Exporting trajectories
        mode = 0
        for frequency in  list(parameters['system'].scratch.nmState.frequencies):
            

            self.trajectory = ExportTrajectory(os.path.join(trajectoryPath              ,
                                               'mode'+str(mode)+'.ptGeo')               , 
                                               parameters['system']                     , 
                                               log = None                               )
            
            
            NormalModesTrajectory_SystemGeometry ( parameters['system']                    ,
                                                   self.trajectory                         ,
                                                   mode        = mode                      ,
                                                   cycles      = parameters['cycles']      ,
                                                   frames      = parameters['frames']      ,
                                                   temperature = parameters['temperature'] )
        
            #-------------------------------------------------------------------------------------
            log = os.path.join(trajectoryPath,'mode'+str(mode)+'.ptGeo', 'frequency.log')
            _file = open(log, 'w')
            
            line = ['mode_'+str(mode) +' = '+str(frequency)+ '  (cm^-1)']
            _file.writelines(line)
            _file.close()
            #-------------------------------------------------------------------------------------

        
            mode += 1
        
        self.logFile2.Footer ( )
        self.logFile2.Close()
        self.logFile2 = None





def func (job):
    """ Function doc """
    print(job)


def _run_second_coordinate_in_parallel ( job):
    """ Function doc """
    #print('hello')
    i          = job[0]
    pkl        = job[1]
    parameters = job[2]
    system     = job[3]
    
    
    hamiltonian = get_hamiltonian (system)
    ''' 
    
    '''
    if hamiltonian in ['DFTB QC Model', 'ORCA QC Model', 'external']:
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
        #                                         logFrequency           = parameters['log_frequency'],
        #                                         maximumIterations      = parameters['maxIterations'],
        #                                         rmsGradientTolerance   = parameters['rmsGradient'], 
        #                                         log                    =                      None)
        ##-------------------------------------------------------------------------------------------------------------
        try:
            #-------------------------------------------------------------------------------------------------------------
            ConjugateGradientMinimize_SystemGeometry(system                                ,                
                                                     logFrequency           = parameters['log_frequency'],
                                                     maximumIterations      = parameters['maxIterations'],
                                                     rmsGradientTolerance   = parameters['rmsGradient'], 
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
        

        
    

def get_hamiltonian (system):
    """ Function doc """
    
    hamiltonian   = getattr(system.qcModel, 'hamiltonian', False)
    
    if hamiltonian:
        pass
    
    else:
        try:
            itens = system.qcModel.SummaryItems()
            #print(itens)
            hamiltonian = itens[0][0]
        except:
            hamiltonian = 'external'
    
    return hamiltonian


