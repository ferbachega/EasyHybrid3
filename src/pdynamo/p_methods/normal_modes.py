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
from pScientific.Symmetry      import Find3DGraphPointGroup
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
        # [EN] BUG FIXED (reported by the user: the Process Manager's
        # "open log" didn't work for Normal Modes -- this method used
        # the CORRECT parameter key ('trajectory_name', matching what
        # simulations_mixin._configure_logfile() expects) but never
        # corrected parameters['logfile'] to match its OWN actual
        # filename scheme ('{trajectory_name}.log' inside the trajectory
        # folder, not the generic 'output.log' _configure_logfile()
        # assumes, nor PDYNAMO3_SCRATCH as the base folder). Same self-
        # correcting pattern molecular_dynamics.py/chain_of_states.py
        # already use.
        parameters['logfile'] = os.path.join(trajectoryPath, parameters['trajectory_name']+'.log')
        parameters['system'].Summary(log = self.logFile2)
        self.logFile2.Header ( )
        
        # . Calculate the normal modes.
        NormalModes_SystemGeometry ( parameters['system'], modify = ModifyOption.Project, log = self.logFile2 )

        # . Identify the vibrational symmetry (irreducible representation)
        # of each mode -- best-effort: point-group/character-table
        # symmetry perception can fail or be ambiguous for some
        # geometries (e.g. numerical noise blurring an exact symmetry),
        # so a failure here degrades to "?" per mode rather than
        # aborting the whole Normal Modes job. Symmetry perception uses
        # a CLONE of the coordinates -- Find3DGraphPointGroup() re-centers
        # whatever Coordinates3 it is given on the (mass-)weighted
        # center, and that must not silently translate the system's own
        # live coordinates that NormalModesTrajectory_SystemGeometry()
        # below still uses as the reference structure for every mode's
        # oscillation trajectory.
        try:
            coordinatesForSymmetry = Clone ( parameters['system'].coordinates3 )
            masses   = Array.FromIterable ( [ atom.mass for atom in parameters['system'].atoms ] )
            numbers  = [ atom.atomicNumber for atom in parameters['system'].atoms ]
            pgReport = Find3DGraphPointGroup ( numbers, coordinatesForSymmetry, doCharacterSymmetryOperations = True, log = None, weights = masses )
            self.logFile2.Paragraph ( "Point group: {:s}".format ( pgReport["Point Group"].label ) )
            modeIRs  = NormalModes_IrreducibleRepresentations ( parameters['system'], results = pgReport, log = self.logFile2 )
        except Exception as error:
            self.logFile2.Paragraph ( "Vibrational symmetry perception failed: {:s}".format ( str ( error ) ) )
            modeIRs = None

        # . Infrared intensities -- opt-in (parameters['compute_ir_intensities'])
        # because NormalModes_InfraredIntensities() does 2 extra energy +
        # dipole-moment evaluations per Cartesian coordinate per FREE
        # atom (one to a side by finite difference), i.e. roughly 6N
        # extra single points for N free atoms -- cheap for the small
        # semiempirical systems this was validated on, but potentially
        # very expensive for a large QC region or a DFT model, so it is
        # never computed unless the user explicitly asks for it in the
        # Run Normal Modes window.
        intensities = None
        if parameters.get ( 'compute_ir_intensities', False ):
            # [EN] System.DipoleMoment() (pDynamo3, pMolecule/System.py)
            # only keeps ONE model's result -- it loops over ("mmModel",
            # "qcModel") and overwrites "dipole" each time a model is
            # present, rather than summing them (a source comment there
            # attributes this to a deliberate past fix: "breaks dipoles
            # from QCModelDFTB/ORCA... not required for QCModelMNDO").
            # For a real QM/MM system (both models defined) this means
            # only the QC region's own dipole survives -- the MM
            # region's point-charge contribution is silently dropped,
            # never summed in. Every IR intensity computed below is
            # therefore only the QC region's contribution for such a
            # system, not the true whole-system intensity -- flag this
            # loudly in the log rather than let it look like a complete
            # result. Applies regardless of which QC engine is active
            # (MNDO, XTB, ORCA, ...); this is a limitation of System.
            # DipoleMoment() itself, not of any one QC model.
            if ( parameters['system'].mmModel is not None ) and ( parameters['system'].qcModel is not None ):
                self.logFile2.Paragraph (
                    "WARNING: this is a QM/MM system (both an MM and a QC model are defined). "
                    "System.DipoleMoment() only returns the QC region's own dipole moment -- "
                    "the MM region's point-charge contribution is NOT included/summed. "
                    "The infrared intensities below are therefore incomplete for this hybrid "
                    "system (missing the MM region's direct contribution), not the true "
                    "whole-system spectrum."
                )
            try:
                irState     = NormalModes_InfraredIntensities ( parameters['system'], log = self.logFile2 )
                intensities = list ( irState.intensities )
            except Exception as error:
                self.logFile2.Paragraph ( "Infrared intensity calculation failed: {:s}".format ( str ( error ) ) )
                intensities = None

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

            symmetry  = modeIRs[mode]     if modeIRs     is not None else '?'
            intensity = intensities[mode] if intensities is not None else '?'
            line = ['mode_'+str(mode) +' = '+str(frequency)+ '  (cm^-1)  symmetry='+str(symmetry)+'  intensity='+str(intensity)]
            _file.writelines(line)
            _file.close()
            #-------------------------------------------------------------------------------------

        
            mode += 1
        
        self.logFile2.Footer ( )
        self.logFile2.Close()
        self.logFile2 = None
