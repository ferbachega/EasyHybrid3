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

        results = {}
        
        options = dict ( bins                 = parameters['bins'                 ],
                         logFrequency         = parameters['logFrequency'         ],
                         maximumIterations    = parameters['maximumIterations'    ],
                         rmsGradientTolerance = parameters['rmsGradientTolerance' ],
                         temperature          = parameters['temperature'          ],
                         log                  = self.logFile2 )
        resamples = int ( parameters.get ( 'resamples', 0 ) or 0 )
        try:
            if resamples > 1:
                # . 2026-10-23: errors by bootstrap resampling (standard error and 95% confidence interval)
                state = WHAM_Bootstrapping ( parameters['file_names'], resamples = resamples, **options )
            else:
                state = WHAM_ConjugateGradientMinimize ( parameters['file_names'], **options )
        except TypeError:
            return False, 'Error on pDynamo WHAM_ConjugateGradientMinimize! \n\nPlease, check your input data.'

        if resamples > 1 and not state.get ( "Converged", True ):
            return False, 'The WHAM equations did not converge, so no bootstrapping was done.\n\nIncrease the maximum number of iterations or the gradient tolerance.'

        #try:
        # . Write the PMF to a file.
        histogram = state["Histogram"]
        pmf       = state["PMF"      ]
        
        PMF_file = os.path.join(parameters['folder'], parameters['logfile']+'_pmf.log')
        if parameters['type'] == 1:
            histogram.ToTextFileWithData (  PMF_file , [ pmf ], format = "{:20.3f} {:20.3f} {:20.3f}\n" )
        else:
            histogram.ToTextFileWithData (  PMF_file , [ pmf ], format = "{:20.3f} {:20.3f}\n" )
        results['pmf'] = PMF_file
        if "Standard Error" in state:
            # . a separate file, so the PMF file keeps its format: bin(s), PMF, standard error, lower and upper 95% limits
            errors_file = os.path.join(parameters['folder'], parameters['logfile']+'_pmf_errors.log')
            columns = len ( histogram.dimensions ) + 4
            histogram.ToTextFileWithData ( errors_file, [ pmf, state["Standard Error"], state["Lower 95% Confidence Interval"],
                                                          state["Upper 95% Confidence Interval"] ],
                                           format = " ".join ( ["{:20.3f}"] * columns ) + "\n" )
            results['pmf_errors'] = errors_file
            results['resamples']  = state.get ( "Number Of Resamples", resamples )
        #except TypeError:
        #    return False, 'Error on Write the PMF to a file!\n\n Please, check your input data.'


        results['type'] = parameters['type']
        results['histograms'] = []
        
        #try:
        files = parameters['file_names']
        try:
            for _file in files:
                dprint (_file)
                basename = os.path.basename(_file)
                output_file = os.path.join(parameters['folder'], parameters['logfile']+'_'+basename+'_histogram.dat')
                # . Histogram the trajectory data.
                handler   = SystemRestraintTrajectoryDataHandler.FromTrajectoryPaths ( [_file] )
                histogram = handler.HistogramData ( [ 360 ] )
                counts    = [ float ( count ) for count    in histogram.counts ]
                histogram.ToTextFileWithData ( output_file , [ counts ], format = "{:20.3f} {:20.3f}\n" )
                results['histograms'].append(output_file)
        except:
            pass


        msg = 'WHAM calculation performed successfully! \n\nPlease check:' + PMF_file
        if not state.get ( "Converged", True ):
            msg += ('\n\nWARNING: the WHAM equations did not converge in {} iterations - the PMF may be '
                    'wrong. Increase the maximum number of iterations.'.format(parameters['maximumIterations']))
        if results.get('pmf_errors'):
            msg += '\n\nBootstrap errors ({} resamples): {}'.format(results['resamples'], results['pmf_errors'])
        results['msg'] = msg

        return True, results
