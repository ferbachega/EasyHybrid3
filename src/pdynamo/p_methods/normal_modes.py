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
