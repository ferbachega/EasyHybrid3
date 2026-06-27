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
            self.logFile2 = TextLogFileWriter.WithOptions ( path = parameters['logfile'] )
            parameters['system'].Summary(log = self.logFile2)
            self.logFile2.Header ( )
            
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
                                                     logFrequency           = parameters['logFrequency'],
                                                     maximumIterations      = parameters['maximumIterations'],
                                                     rmsGradientTolerance   = parameters['rmsGradientTolerance'],
                                                     log                    = self.logFile2 )
        else:            
            
            ConjugateGradientMinimize_SystemGeometry(parameters['system']                                                ,                
                                                     logFrequency           = parameters['logFrequency']                ,
                                                     trajectories           = [(self.trajectory, parameters['save_frequency'])],
                                                     maximumIterations      = parameters['maximumIterations']                ,
                                                     rmsGradientTolerance   = parameters['rmsGradientTolerance']                  ,
                                                     log                    = self.logFile2                                   )

    def _run_steepest_descent(self, parameters):

        if parameters['trajectory_name'] == None:
            SteepestDescentMinimize_SystemGeometry(parameters['system']                               ,               
                                                logFrequency            = parameters['logFrequency'] ,
                                                maximumIterations       = parameters['maximumIterations'] ,
                                                rmsGradientTolerance    = parameters['rmsGradientTolerance'],
                                                log                    = self.logFile2   )
        else:
            SteepestDescentMinimize_SystemGeometry(parameters['system']                                               ,               
                                                logFrequency            = parameters['logFrequency']                 ,
                                                trajectories            = [(self.trajectory, parameters['save_frequency'])],
                                                maximumIterations       = parameters['maximumIterations']                 ,
                                                rmsGradientTolerance    = parameters['rmsGradientTolerance']                   ,
                                                log                     = self.logFile2                                   )

    def _run_LBFGS(self, parameters):
        
        if parameters['trajectory_name'] == None:
            LBFGSMinimize_SystemGeometry(parameters['system']                               ,                
                                         logFrequency         = parameters['logFrequency'] ,
                                         maximumIterations    = parameters['maximumIterations'] ,
                                         rmsGradientTolerance = parameters['rmsGradientTolerance'],
                                         log                  = self.logFile2)
        else:

            LBFGSMinimize_SystemGeometry(parameters['system']                                                ,                
                                         logFrequency         = parameters['logFrequency']                  ,
                                         trajectories         = [(self.trajectory, parameters['save_frequency'])] ,
                                         maximumIterations    = parameters['maximumIterations']                  ,
                                         rmsGradientTolerance = parameters['rmsGradientTolerance']                    ,
                                         log                  = self.logFile2                                   )  

    def _run_quasi_newton(self, parameters):
        '''
        Class method to apply the Quaisi-Newton minimizer
        '''        
        if parameters['trajectory_name'] == None:
            QuasiNewtonMinimize_SystemGeometry( parameters['system']                               ,               
                                                logFrequency         = parameters['logFrequency'] ,
                                                maximumIterations    = parameters['maximumIterations'] ,
                                                rmsGradientTolerance = parameters['rmsGradientTolerance'],
                                                log                  = self.logFile2   )
        else:
            QuasiNewtonMinimize_SystemGeometry( parameters['system']                                                ,
                                                logFrequency         = parameters['logFrequency']                  ,
                                                trajectories         = [(self.trajectory, parameters['save_frequency'])] ,
                                                maximumIterations    = parameters['maximumIterations']                  ,
                                                rmsGradientTolerance = parameters['rmsGradientTolerance']                    ,
                                                log                  = self.logFile2                                   )
    
    def _run_FIRE(self, parameters):
        '''
        '''
        if parameters['trajectory_name'] == None:
            FIREMinimize_SystemGeometry( parameters['system']                               ,          
                                         logFrequency         = parameters['logFrequency'] ,
                                         maximumIterations    = parameters['maximumIterations'] ,
                                         rmsGradientTolerance = parameters['rmsGradientTolerance'],
                                         log                  = self.logFile2   )
        else:

            FIREMinimize_SystemGeometry( parameters['system']                                                ,                
                                         logFrequency         = parameters['logFrequency']                  ,
                                         trajectories         = [(self.trajectory, parameters['save_frequency'])] ,
                                         maximumIterations    = parameters['maximumIterations']                  ,
                                         rmsGradientTolerance = parameters['rmsGradientTolerance']                    ,
                                         log                  = self.logFile2                                   )
