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

from pprint import pprint
import os


HEADER ='''
----------------------------------------------------------------------------- 
                                EasyHybrid 3.0                                
                   - A pDynamo3 graphical user interface -                    
----------------------------------------------------------------------------- 

''' 
class LogFile:
    """ Class doc """
    
    def __init__ (self, parameters):
        """ Class initialiser """
        pass


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
            print('\n\n\n')
            print('full_path_trajectory:', full_path_trajectory)
            self.trajectory = ExportTrajectory(full_path_trajectory, parameters['system'], log=None )
            
            self.logFile2 = TextLogFileWriter.WithOptions ( path = os.path.join(full_path_trajectory, 'output.log') )
            parameters['system'].Summary(log = self.logFile2)
            self.logFile2.Header ( )
        else:
            pass
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        
        if parameters['integrator'] == 'Verlet':
            self._velocity_verlet_dynamics (parameters)
        
        elif parameters['integrator'] == 'LeapFrog':
            self._leap_frog_dynamics(parameters)
        
        elif parameters['integrator'] == 'Langevin':
            self._langevin_dynamics(parameters)
        
        else:
            pass    
        

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

    def run (self, parameters):
        """ Function doc """
        full_path_trajectory =os.path.join(parameters['folder'], 
                              parameters['traj_folder_name']+".ptGeo")
        os.mkdir(
                 full_path_trajectory
                 )
        
        # - - - - - - - - - - - - - Checking trajectory - - - - - - - - - - - - - -
        self.logFile2 = TextLogFileWriter.WithOptions ( path = os.path.join(full_path_trajectory, 'output.log') )
        parameters['system'].Summary(log = self.logFile2)
        self.logFile2.Header ( )
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


        
        if parameters['second_coordinate']:
            self._run_scan_2D(parameters)
            
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

        if parameters["second_coordinate"] :
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
        if parameters["rc_type_1"] == 'simple_distance':
            text = text + "\n"
            text = text + "\n----------------------- Coordinate 1 - Simple-Distance -------------------------"
            text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['ATOMS_RC1'][0], parameters['ATOMS_RC1_NAMES'][0] )
            text = text + "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['ATOMS_RC1'][1], parameters['ATOMS_RC1_NAMES'][1] )
            text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['nsteps_RC1']  , parameters['force_constant_1']   )
            text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['dminimum_RC1'], parameters['maxIterations']      )
            text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['dincre_RC1']  , parameters['rmsGradient']        )
            text = text + "\n--------------------------------------------------------------------------------"

        
        elif parameters["rc_type_1"] == 'multiple_distance':
            text = text + "\n"
            text = text + "\n---------------------- Coordinate 1 - multiple-Distance ------------------------"	
            text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['ATOMS_RC1'][0]    , parameters['ATOMS_RC1_NAMES'][0] )
            text = text + "\nATOM2*                 =%15i  ATOM NAME2             =%15s"     % (parameters['ATOMS_RC1'][1]    , parameters['ATOMS_RC1_NAMES'][1] )
            text = text + "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['ATOMS_RC1'][2]    , parameters['ATOMS_RC1_NAMES'][2] )
            text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['nsteps_RC1']      , parameters['force_constant_1']   ) 
            text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['dminimum_RC1']    , parameters['maxIterations']      )
            text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['dincre_RC1']      , parameters['rmsGradient']        )
            text = text + "\nSigma atom1 - atom3    =%15.5f  Sigma atom3 - atom1    =%15.5f" % (parameters['sigma_pk1pk3_rc1'], parameters['sigma_pk3pk1_rc1']      )
            text = text + "\n--------------------------------------------------------------------------------"

        else:
            pass
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------
        
        
        
        
        if parameters["second_coordinate"] :

            '''This part writes the parameters used in the second reaction coordinate'''
            #-----------------------------------------------------------------------------------------------------------------------------------------------------------
            if parameters["rc_type_2"] == 'simple_distance':
                text = text + "\n"
                text = text + "\n----------------------- Coordinate 2 - Simple-Distance -------------------------"
                text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['ATOMS_RC2'][0], parameters['ATOMS_RC2_NAMES'][0] )
                text = text + "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['ATOMS_RC2'][1], parameters['ATOMS_RC2_NAMES'][1] )
                text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['nsteps_RC2']  , parameters['force_constant_1']   )
                text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['dminimum_RC2'], parameters['maxIterations']      )
                text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['dincre_RC2']  , parameters['rmsGradient']        )
                text = text + "\n--------------------------------------------------------------------------------"

            
            elif parameters["rc_type_2"] == 'multiple_distance':
                text = text + "\n"
                text = text + "\n---------------------- Coordinate 2 - multiple-Distance ------------------------"	
                text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['ATOMS_RC2'][0]    , parameters['ATOMS_RC2_NAMES'][0] )
                text = text + "\nATOM2*                 =%15i  ATOM NAME2             =%15s"     % (parameters['ATOMS_RC2'][1]    , parameters['ATOMS_RC2_NAMES'][1] )
                text = text + "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['ATOMS_RC2'][2]    , parameters['ATOMS_RC2_NAMES'][2] )
                text = text + "\nNUMBER OF STEPS        =%15i  FORCE CONSTANT         =%15i"     % (parameters['nsteps_RC2']      , parameters['force_constant_2']   ) 
                text = text + "\nDMINIMUM               =%15.5f  MAX INTERACTIONS       =%15i"   % (parameters['dminimum_RC2']    , parameters['maxIterations']      )
                text = text + "\nSTEP SIZE              =%15.7f  RMS GRAD               =%15.7f" % (parameters['dincre_RC2']      , parameters['rmsGradient']        )
                text = text + "\nSigma atom1 - atom3    =%15.5f  Sigma atom3 - atom1    =%15.5f" % (parameters['sigma_pk1pk3_rc2'], parameters['sigma_pk3pk1_rc2']      )
                text = text + "\n--------------------------------------------------------------------------------"
            else:
                pass
            #-----------------------------------------------------------------------------------------------------------------------------------------------------------
            
        
        
        
        
        
        
        
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------
        '''This part writes the header of the frames distances and energy'''
        if parameters["second_coordinate"] :
            text = text + "\n\n--------------------------------------------------------------------------------"
            text = text + "\n   Frame i  /  j        RCOORD-1             RCOORD-2                Energy     "
            text = text + "\n--------------------------------------------------------------------------------"




        else:
            if parameters["rc_type_1"] == 'simple_distance':
                text = text + "\n\n-------------------------------------------------------------"
                text = text + "\n           Frame    dist-ATOM1-ATOM2             Energy      "
                text = text + "\n-------------------------------------------------------------"
            
            elif parameters["rc_type_1"] == 'multiple_distance':
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
        atom1 = parameters['ATOMS_RC1'][0]
        atom2 = parameters['ATOMS_RC1'][1]                   
        #---------------------------------
        restraints = RestraintModel()
        parameters['system'].DefineRestraintModel( restraints )                     
        #----------------------------------------------------------------------------------------
        
        arq = self.write_header(parameters)
        data = []
        
        for i in range(parameters['nsteps_RC1']):       
            distance = parameters['dminimum_RC1'] + ( parameters['dincre_RC1'] * float(i) )
            print(parameters["rc_type_1"])
            
            
            '''----------------------------------------------------------------------------------------------------------------'''
            if parameters["rc_type_1"] == 'simple_distance':
                #---------------------------------------------------------------------------------------------------------
                rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['force_constant_1'])
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
                '''----------------------------------------------------------------------------------------------------------------'''

            
                '''----------------------------------------------------------------------------------------------------------------'''
            elif parameters["rc_type_1"] == 'multiple_distance':
                #--------------------------------------------------------------------
                atom3   = parameters['ATOMS_RC1'][2]
                weight1 = parameters['sigma_pk1pk3_rc1'] #self.sigma_a1_a3[0]
                weight2 = parameters['sigma_pk3pk1_rc1'] #self.sigma_a3_a1[0] 
                
                rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['force_constant_1'])
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
                '''----------------------------------------------------------------------------------------------------------------'''

            else:
                pass
        #---------------------------------------
        parameters['system'].DefineRestraintModel(None)
        pprint(data)
        
        
    def _run_scan_2D (self, parameters):
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
        print(parameters["rc_type_1"])
        
        
        '''This part of the program is crucial because the Clone function 
        is unable to handle data that has been allocated using the GTK library. 
        When a system is generalized using EasyHybrid, it receives information 
        from the treeview, enabling easy interconnection of the data.'''
        backup = []
        backup.append(parameters['system'].e_treeview_iter)
        backup.append(parameters['system'].e_liststore_iter)
        parameters['system'].e_treeview_iter   = None
        parameters['system'].e_liststore_iter  = None
        
       
        
        
        '''We have to first run a sequence (first line) of optimizations sequentially. 
        This will generate frames (0,0) to (n, 0) serially, where n is the number of 
        steps from coordinate x)
        '''
        try:
            joblist = []
            j = 0
            for i in range(parameters['nsteps_RC1']):       
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
                    


                '''reaction coordinate 2 - ONLY at 0 first'''
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
                data[(i,j)] = [RC1_d1_minus_d2, RC2_d1_minus_d2, energy]    
                #--------------------------------------------------------------------------------------
                
                
                
                #--------------------------------------------------------------------------------------
                #                             Exporting Coordinates
                #--------------------------------------------------------------------------------------
                pkl = os.path.join( parameters['folder'], 
                                      parameters['traj_folder_name']+".ptGeo", 
                                      "frame{}_{}.pkl".format(i,j) )
                                       
                Pickle( pkl, parameters['system'].coordinates3 )
                #--------------------------------------------------------------------------------------
                
                

                #system = parameters['system']
                #backup = []
                #backup.append(system.e_treeview_iter)
                #backup.append(system.e_liststore_iter)
                
                #For some reason the previous dictionary had problems when submitted to the pool
                new_parameters = {
                                 
                                 'ATOMS_RC1'       : parameters['ATOMS_RC1'],
                                 'ATOMS_RC2'       : parameters['ATOMS_RC2'],
                                 'rc_type_1'       : parameters['rc_type_1'],
                                 'rc_type_2'       : parameters['rc_type_2'],
                                 'dminimum_RC1'    : parameters['dminimum_RC1'],
                                 'dminimum_RC2'    : parameters['dminimum_RC2'],
                                 
                                 'sigma_pk1pk3_rc1': parameters['sigma_pk1pk3_rc1'],
                                 'sigma_pk3pk1_rc1': parameters['sigma_pk3pk1_rc1'],
                                 'force_constant_1': parameters['force_constant_1'],
                                 'nsteps_RC1'      : parameters['nsteps_RC1'],
                                 'dincre_RC1'      : parameters['dincre_RC1'],
                                 
                                 
                                 'sigma_pk1pk3_rc2': parameters['sigma_pk1pk3_rc2'],
                                 'sigma_pk3pk1_rc2': parameters['sigma_pk3pk1_rc2'],
                                 'force_constant_2': parameters['force_constant_2'],
                                 'nsteps_RC2'      : parameters['nsteps_RC2'],
                                 'dincre_RC2'      : parameters['dincre_RC2'],
                                 
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
            for i in range(parameters['nsteps_RC1']):
                for j in range(parameters['nsteps_RC2']):
            
                    text = "\nDATA  %4i  %4i     %13.12f       %13.12f       %13.12f"% (int(i), int(j),  float(data[(i,j)][0]), float(data[(i,j)][1]), float(data[(i,j)][2]))
                    arq.write(text)
            #--------------------------------------------------------------------------------------
            parameters['system'].e_treeview_iter   = backup[0]
            parameters['system'].e_liststore_iter  = backup[1]

        except:
            parameters['system'].e_treeview_iter   = backup[0]
            parameters['system'].e_liststore_iter  = backup[1]


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
        print(parameters["rc_type_1"])

        for i in range(parameters['nsteps_RC1']):       
            distance = parameters['dminimum_RC1'] + ( parameters['dincre_RC1'] * float(i) )

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
        
        
        
        #--------------------------------------------------------------------------------------
        #                             Writing the log information 
        #--------------------------------------------------------------------------------------
        for i in range(parameters['nsteps_RC1']):
            for j in range(parameters['nsteps_RC2']):
        
                text = "\nDATA  %4i  %4i     %13.12f       %13.12f       %13.12f"% (int(i), int(j),  float(data[(i,j)][0]), float(data[(i,j)][1]), float(data[(i,j)][2]))
                arq.write(text)
        #--------------------------------------------------------------------------------------
    


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
        #--------------------------------------------------------------------------------------
        print(i, j, 'Energy:', energy, opt_convergency)
    return data
        

        
    





