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
from util.debug import dprint
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
            
         
            
            dprint ('\nChecking for ORCA files at: ', scratch)

            try:
                os.rename(os.path.join(scratch, 'orcaJob.log'),   os.path.join(folder , output_name+'.orca.log'))
                dprint('Renaming logfile to:', os.path.join(folder , output_name+'.orca.log'))
            except:
                dprint('File {} not found at:'.format(output_name+'.orca.log'), scratch)
            
            try:
                os.rename(os.path.join(scratch, 'orcaJob.gbw'), os.path.join(folder , output_name+'.orca.gbw'))
                dprint('Renaming gbwfile to:', os.path.join(folder , output_name+'.orca.gbw'))
            except:
                dprint('File {} not found at:'.format(output_name+'.orca.gbw'), scratch)
            dprint ('\n')


# . Job files known to be written by pDynamo3's XTB QC model (see
#   QCModelXTBState.DeterminePaths in util/extras/QCModelXTB.py), plus a few
#   extra files xTB itself can write depending on the keywords used (e.g.
#   "molden.input" when "--molden" is passed -- this one is NOT prefixed
#   with the job name, xTB always writes it under that fixed name in the
#   working/scratch directory). The keys are what the "Backup Files"
#   checkboxes in the xTB setup window
#   (gui/windows/setup/windows_and_dialogs/qc_setup/setup_xtb.py) refer to.
#
# . Each entry is ( source_filename_or_None, source_extension_or_None, description ).
#   If "source_filename" is given, it is used as-is (fixed name, independent
#   of the job name). Otherwise the source file is "<_XTB_JOB_NAME>.<source_extension>".
_XTB_JOB_NAME       = 'XTBJob'
_XTB_BACKUP_FILES = {
    'log'    : { 'extension' : 'log'    , 'description' : 'Output'                 } ,
    'inp'    : { 'extension' : 'inp'    , 'description' : 'Input'                  } ,
    'engrad' : { 'extension' : 'engrad' , 'description' : 'Gradient'               } ,
    'coord'  : { 'extension' : 'coord'  , 'description' : 'Coordinates'            } ,
    'pc'     : { 'extension' : 'pc'     , 'description' : 'Point charges'          } ,
    'pcgrad' : { 'extension' : 'pcgrad' , 'description' : 'Point-charge gradient'  } ,
    # . Fixed name, only written when "--molden" is among the xTB keywords.
    'molden' : { 'source_name' : 'molden.input', 'extension' : 'molden.input', 'description' : 'Molden orbitals' } ,
}


def backup_xtb_files (system, output_folder = None, output_name = None, files = None):
    """
    Copies the requested xTB job files from the scratch folder to a
    permanent output folder, analogous to backup_orca_files() above, but
    letting the user choose *which* files to keep (see the "Backup Files"
    option added to the xTB setup window). Only files that were actually
    requested (and that xTB actually wrote -- e.g. "molden.input" is only
    written when "--molden" is among the xTB keywords) are copied; anything
    else is silently skipped, exactly like the ORCA counterpart.

    "files" is a list of keys from _XTB_BACKUP_FILES (e.g. ['log',
    'molden']). If not given, falls back to the "e_xtb_backup_files"
    attribute that EasyHybrid attaches to the system when the xTB QC model
    is defined (see session.define_a_new_QCModel), or just ['log'] if that
    is not set either.
    """
    
    
    if not system.qcModel:
        return

    items = system.qcModel.SummaryItems()
    '''Checking if xTB is being used'''
    if items[0][0] != 'XTB QC Model':
        return

    if files is None:
        files = getattr ( system, 'e_xtb_backup_files', [ 'log' ] )
    if not files:
        return

    # The xTB job files now live in a unique random subfolder created per run
    # (see QCModelXTB.DeterminePaths/_resolve_scratch), NOT directly in the
    # configured scratch base. Read the real folder via activeScratch; fall back
    # to .scratch for older models that don't expose it.
    
    scratch = os.path.dirname(system.qcState.paths["Coord"])
    
    #'''
    #scratch = getattr ( system.qcModel, 'activeScratch', None ) or system.qcModel.scratch
    
    #print('scratch', scratch)
    #print('paths', system.qcState.paths)
    _time = time.asctime()
    #'''   
    
    if output_folder is None:
        folder = scratch
    else:
        folder = output_folder

    if output_name is None:
        output_name = 'XTBJob' + _time

    dprint ('\nChecking for xTB files at: ', scratch)
    
    #print('backup_xtb_files', system, output_folder  , output_name  , files )
    
    for key in files:
        file_info = _XTB_BACKUP_FILES.get ( key )
        if file_info is None:
            dprint ('Unknown xTB backup file type requested: {}'.format ( key ))
            continue

        description = file_info['description']
        extension   = file_info['extension']
        source_name = file_info.get ( 'source_name', '{}.{}'.format ( _XTB_JOB_NAME, extension ) )
        if extension == 'log':
            extension='out'
        
        destination_name = '{}.xtb.{}'.format ( output_name, extension )

        try:
            os.rename ( os.path.join ( scratch, source_name ), os.path.join ( folder, destination_name ) )
            dprint ('Renaming {} ({}) file to:'.format ( description, source_name ), os.path.join ( folder, destination_name ))
        except Exception:
            dprint ('File {} not found at:'.format ( source_name ), scratch)

    dprint ('\n')


def backup_qc_files (system, output_folder = None, output_name = None, xtb_files = None):
    """
    Single entry point for backing up external-QC-program files, used by ALL
    calculation types (single point, geometry optimization, surface scan,
    molecular dynamics, ...). It figures out which external QC engine (if any)
    the system currently uses and dispatches to the matching backup function.

    Safe to call unconditionally -- it is a no-op for systems without a QC
    model, or with a QC model whose engine has no backup routine registered
    (e.g. MOPAC, DFTB+, MNDO).

    EXTENSIBILITY: to support a new QC program, write its backup function
    (following backup_orca_files / backup_xtb_files) and add ONE entry to
    _QC_BACKUP_DISPATCH below, keyed by the engine's SummaryItems label. Every
    calculation type then supports it automatically, because they all funnel
    through here.
    """
    #print('aqui')
    
    if not system.qcModel:
        return

    engine = system.qcModel.SummaryItems()[0][0]

    backup_fn = _QC_BACKUP_DISPATCH.get ( engine )
    
    if backup_fn is None:
        return   # engine with no backup routine -- nothing to do
    '''
    print( 'backup_qc_files')
    print( system        ,
           output_folder ,
           output_name    ,
           xtb_files       )
    '''
    
    
    backup_fn ( system        = system        ,
                output_folder = output_folder ,
                output_name   = output_name    ,
                xtb_files     = xtb_files       )


def _backup_orca_dispatch ( system, output_folder, output_name, xtb_files = None ):
    """Adapter so ORCA's backup matches the common dispatch signature."""
    backup_orca_files ( system = system, output_folder = output_folder, output_name = output_name )


def _backup_xtb_dispatch ( system, output_folder, output_name, xtb_files = None ):
    """Adapter so xTB's backup matches the common dispatch signature."""
    backup_xtb_files ( system = system, output_folder = output_folder,
                       output_name = output_name, files = xtb_files )


# Registry of engine label -> backup function. Add new QC programs here; every
# calculation type picks them up automatically via backup_qc_files().
_QC_BACKUP_DISPATCH = {
    'ORCA QC Model' : _backup_orca_dispatch ,
    'XTB QC Model'  : _backup_xtb_dispatch  ,
    # 'DFTB+ QC Model' : _backup_dftbplus_dispatch ,   # <- future engines go here
    # 'MOPAC QC Model' : _backup_mopac_dispatch    ,
}


def _backup_qc_files_unused_old ( system, output_folder = None, output_name = None, xtb_files = None ):
    """ retained for reference """
    if not system.qcModel:
        return

    engine = system.qcModel.SummaryItems()[0][0]

    if engine == 'ORCA QC Model':
        backup_orca_files ( system = system, output_folder = output_folder, output_name = output_name )
    elif engine == 'XTB QC Model':
        backup_xtb_files  ( system = system, output_folder = output_folder, output_name = output_name, files = xtb_files )


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
    dprint(text)
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
    dprint(job)
