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
        parameters['logfile'] = os.path.join(trajectoryPath, 'output.log')
        parameters['system'].Summary(log = self.logFile2)
        self.logFile2.Header ( )
        '''-----------------------------------------------------------------------------------------------------'''        
        #logFile.Header ( )
        system = parameters['system']

        # . Assign the reactant and product coordinates.
        if parameters['external_coords']:
            reactants = ImportCoordinates3 (parameters['reac_coordinates'])
            products  = ImportCoordinates3 (parameters['prod_coordinates'])
        else:
            reactants = parameters['reac_coordinates']
            products  = parameters['prod_coordinates']

        # . Get an initial path.
        GrowingStringInitialPath ( system                         , 
                                   parameters['number_of_structures'], 
                                   reactants                      , 
                                   products                       , 
                                   trajectoryPath                 ,
                                   log  = self.logFile2)

        # . Optimization.
        trajectory = ExportTrajectory ( trajectoryPath, system, append = True )
        
        
        #poolfactory
        if parameters["poolFactory"] == None or parameters["poolFactory"] <=1:
            parameters["poolFactory"] = None
        else:
            mProc = parameters["poolFactory"]
            #parameters["poolFactory"] = SGOFProcessPoolFactory ( maximumProcesses = 6, poolType = "Multiprocessing" )
            parameters["poolFactory"] = SGOFProcessPoolFactory ( maximumProcesses = mProc, poolType = "Multiprocessing" )
                                                                 
        ChainOfStatesOptimizePath_SystemGeometry ( system                                                   ,
                                                   trajectory                                               ,
                                                   
                                                   logFrequency                  = parameters['logFrequency']       ,
                                                   maximumIterations             = parameters['maximumIterations']   ,
                                                   rmsGradientTolerance          = parameters['rmsGradientTolerance'],
                                                   
                                                   springForceConstant           = parameters['springForceConstant']          ,
                                                   splineRedistributionTolerance = parameters['splineRedistributionTolerance'],
                                                   useSplineRedistribution       = parameters["useSplineRedistribution"        ],

                                                   fixedTerminalImages                 = parameters['fixedTerminalImages']                          ,
                                                   
                                                   forceOneSingleImageOptimization     = parameters['forceOneSingleImageOptimization'],            
                                                   forceSingleImageOptimizations       = parameters['forceSingleImageOptimizations'],    
                                                   forceSplineRedistributionCheckPerIteration = parameters['forceSplineRedistributionCheckPerIteration'],
                                                   
                                                   freezeRMSGradientTolerance = parameters['freezeRMSGradientTolerance'],
                                                   
                                                   optimizer                 = parameters["optimizer"                 ] ,
                                                   poolFactory               = parameters["poolFactory"               ] ,
                                                   rmsGradientToleranceScale = parameters["rmsGradientToleranceScale" ] ,
                                                   
                                                   
                                                   
                                                   log                                        = self.logFile2                     )



        # . Footer.
        self.logFile2.Footer ( )
        self.logFile2.Close ( )
        self.logFile2 = None
