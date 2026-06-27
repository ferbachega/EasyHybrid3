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
from pdynamo.p_methods._common import backup_orca_files

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
            
            parameters['logfile'] = os.path.join(full_path_trajectory, 'output.log')  
            parameters['system'].Summary(log = self.logFile2)
            
            self.logFile2.Header ( )
        else:
            pass
            self.logFile2 = TextLogFileWriter.WithOptions ( path = parameters['logfile'] )
            parameters['system'].Summary(log = self.logFile2)
            self.logFile2.Header ( )
            
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
                                                        temperatureStart          = parameters['temperatureStart'],
                                                        log                       = self.logFile2
                                                        )

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
                                                        temperatureStop           = parameters['temperature_stop']          ,
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
                                                        temperatureStop           = parameters['temperature_stop'],
                                                        log                       = self.logFile2
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
        if parameters['temperatureControl']:
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
                                                  timeStep               = parameters['timeStep'],
                                                  log                    = self.logFile2
                                                  )
        
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
                                                  timeStep               = parameters['timeStep'],
                                                  log                    = self.logFile2
                                                  )


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
                                              timeStep               = parameters['timeStep'],
                                              log                    = self.logFile2
                                              )
