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


def write_header (parameters, logfile = 'output.log'):
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
    #-------------------------------------------------------------------------------------------------------------------------------------------------------
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
        #text = text + "\nSigma atom1 - atom3    =%15.5f  Sigma atom3 - atom1    =%15.5f" % (parameters['RC1']['sigma_pk1pk3'], parameters['RC1']['sigma_pk3
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
    #-------------------------------------------------------------------------------------------------------------------------------------------------------
    
    
    
    
    if parameters["RC2"] is not None :

        '''This part writes the parameters used in the second reaction coordinate'''
        #---------------------------------------------------------------------------------------------------------------------------------------------------
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
        #---------------------------------------------------------------------------------------------------------------------------------------------------
        
    
    
    
    
    
    
    
    #-------------------------------------------------------------------------------------------------------------------------------------------------------
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
    
    
    #-------------------------------------------------------------------------------------------------------------------------------------------------------
        
        #text = text + "\n\n------------------------------------------------------"
        #text = text + "\n       Frame     distance-pK1-pK2         Energy      "
        #text = text + "\n------------------------------------------------------"
    print(text)
    arq.write(text)
    return arq


def plot_data (matrix_i_j):
    """ Function doc """
    import matplotlib.pyplot as plt
    import numpy as np

    plt.style.use('_mpl-gallery-nogrid')

    # make data with uneven sampling in x
    #x = [-3, -2, -1.6, -1.2, -.8, -.5, -.2, .1, .3, .5, .8, 1.1, 1.5, 1.9, 2.3, 3]
    #X, Y = np.meshgrid(x, np.linspace(-3, 3, 128))
    #Z = (1 - X/2 + X**5 + Y**3)# * np.exp(-X**2 - Y**2)

    # plot
    fig, ax = plt.subplots()

    ax.pcolormesh( matrix_i_j)#, vmin=-0.5, vmax=1.0)

    plt.show()


def func (job):
    """ Function doc """
    print(job)
