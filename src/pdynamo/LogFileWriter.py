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

from util.debug import dprint
from pCore import *
from datetime import datetime
from timeit import default_timer as timer
import os
import re
#*************************************************************
from pprint import pprint
HEADER = '''
#-----------------------------------------------------------------------------#
#                                                                             #
#                                EasyHybrid 3.0                               #
#                   - A pDynamo3 graphical user interface -                   #
#                                                                             #
#-----------------------------------------------------------------------------#
#                                                                             #
#             visit: https://sites.google.com/site/EasyHybrid/                #
#                                                                             #
#                                                                             #
#   EasyHybrid team:                                                          #
#   - Fernando Bachega                                                        #
#   - Igor Barden                                                             #
#   - Luis Fernando S M Timmers                                               #
#   - Martin Field                                                            #
#   - Troy Wymore                                                             #
#                                                                             #
#   Cite this work as:                                                        #
#   J.F.R. Bachega, L.F.S.M. Timmers, L.Assirati, L.B. Bachega, M.J. Field,   #
#   T. Wymore. J. Comput. Chem. 2013, 34, 2190-2196. DOI: 10.1002/jcc.23346   #
#                                                                             #
#-----------------------------------------------------------------------------#
'''


class LogFileWriter:
    '''
    Class to create and handle Logfiles of pDynamo
    '''
    #-----------------------------------------------------------------
    def __init__(self, psystem = None):
        '''
        Class constructor.
        Opens the file and initialize the text variable.
        '''
        self.psystem = psystem
        self.start = timer()
        self.end   = 0 

        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

        #self.filePath   = _filePath
        self.text       = "Log File for Simulation project on pDynamo make by EasyHybrid3.0!\n"
        
        
        self.create_header_text()
        #self.text       += "Starting at: " + dt_string + "\n"
        #self.separator()

        #self.fileObj    = open( self.filePath,"w")

    #===================================================================    
    def add_simulation_parameters_text (self, parameters):
        """ Function doc """
        if parameters['simulation_type'] ==  'Relaxed_Surface_Scan':
            if "ATOMS_RC2" in parameters:
                simulation_type = parameters['simulation_type']+"_2D" 
            else:
                simulation_type = parameters['simulation_type']
            
            self.text += "\n---------------- Simulation Setup -------------------"								#				
            self.text += "\nType                   =   %20s"     % ( simulation_type) 
            self.text += "\nSystem                 =   %20s"     % ( parameters['system_name']) 
            self.text += "\nInitial Coortinates    =   %20s"     % ( parameters['initial_coordinates'])            
            self.text += "\nOptimizer              =   %20s"     % ( parameters['optimizer']) 
            self.text += "\nrmsGradient            =   %20f"     % ( parameters['rmsGradient']) 
            self.text += "\nMaximum Iterations     =   %20i"     % ( parameters['maxIterations']) 
            self.text += "\nMaximum Threads        =   %20i"     % ( parameters['nprocs']) 
            self.text += "\n-----------------------------------------------------"   
            self.text += "\n" 
            
            
            if not "rc_type_1" in parameters: parameters["rc_type_1"] = "None"

            if parameters['rc_type_1'] == 'Distance':                
                
                if len(parameters['ATOMS_RC1']) == 2:  # simple distance                    
                    self.text += "\n----------------------- Coordinate 1 - Simple-Distance -------------------------"								#				
                    self.text += "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % ( parameters['ATOMS_RC1'][0] , parameters['ATOMS_RC1_NAMES'][0] ) 
                    self.text += "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % ( parameters['ATOMS_RC1'][1] , parameters['ATOMS_RC1_NAMES'][1] ) 
                    self.text += "\nSTEPS                  =%15i  FORCE CONSTANT         =%15.3f"  %  ( parameters['nsteps_RC1']   , parameters['force_constant_1']   )  
                    self.text += "\nDMINIMUM               =%15.5f  DINCREMENT             =%15.5f" % ( parameters['dminimum_RC1'] , parameters['dincre_RC1']         )  
                    self.text += "\n--------------------------------------------------------------------------------"                              #
                                                                                                                                                #
                else:                                                                                                    #                                                                                                                                          #
                    self.text += "\n--------------------- Coordinate 1 - Multiple-Distance -------------------------"								#				
                    self.text += "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % ( parameters['ATOMS_RC1'][0]     , parameters['ATOMS_RC1_NAMES'][0]) 
                    self.text += "\nATOM2*                 =%15i  ATOM NAME2             =%15s"     % ( parameters['ATOMS_RC1'][1]     , parameters['ATOMS_RC1_NAMES'][1]) 
                    self.text += "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % ( parameters['ATOMS_RC1'][2]     , parameters['ATOMS_RC1_NAMES'][2]) 
                    self.text += "\nSIGMA ATOM1/ATOM3      =%15.5f  SIGMA ATOM3/ATOM1      =%15.5f" % ( parameters['sigma_pk1pk3_rc1'] , parameters['sigma_pk3pk1_rc1']  ) 
                    self.text += "\nSTEPS                  =%15i  FORCE CONSTANT         =%15.3f"  %  ( parameters['nsteps_RC1']       , parameters['force_constant_1']  ) 	
                    self.text += "\nDMINIMUM               =%15.5f  DINCREMENT             =%15.5f" % ( parameters['dminimum_RC1']     , parameters['dincre_RC1']        ) 	
                    self.text += "\n--------------------------------------------------------------------------------"                              #
                #----------------------------------------------------------------------------------------------------------------------------------#
                self.text += "\n"                  

            if not "rc_type_2" in parameters: parameters["rc_type_2"] = "None"

            if parameters['rc_type_2'] == 'Distance':                
                if parameters['ATOMS_RC2']:
                
                    if len(parameters['ATOMS_RC2']) == 2:  # simple distance
                        
                        self.text += "\n----------------------- Coordinate 2 - Simple-Distance -------------------------"								#				
                        self.text += "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % ( parameters['ATOMS_RC2'][0]  , parameters['ATOMS_RC2_NAMES'][0]) 
                        self.text += "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % ( parameters['ATOMS_RC2'][1]  , parameters['ATOMS_RC2_NAMES'][1]) 
                        self.text += "\nSTEPS                  =%15i  FORCE CONSTANT         =%15.3f"  %  ( parameters['nsteps_RC2']    , parameters['force_constant_2']  ) 	
                        self.text += "\nDMINIMUM               =%15.5f  DINCREMENT             =%15.5f" % ( parameters['dminimum_RC2']  , parameters['dincre_RC2']        ) 	
                        self.text += "\n--------------------------------------------------------------------------------"                              #
                                                                                                                                                    #
                    else:                                                                                                    #                                                                                                                                          #
                        self.text += "\n--------------------- Coordinate 2 - Multiple-Distance -------------------------"								#				
                        self.text += "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % ( parameters['ATOMS_RC2'][0]     , parameters['ATOMS_RC2_NAMES'][0] )
                        self.text += "\nATOM2*                 =%15i  ATOM NAME2             =%15s"     % ( parameters['ATOMS_RC2'][1]     , parameters['ATOMS_RC2_NAMES'][1] )
                        self.text += "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % ( parameters['ATOMS_RC2'][2]     , parameters['ATOMS_RC2_NAMES'][2] )
                        self.text += "\nSIGMA ATOM1/ATOM3      =%15.5f  SIGMA ATOM3/ATOM1      =%15.5f" % ( parameters['sigma_pk1pk3_rc2'] , parameters['sigma_pk3pk1_rc2']   )
                        self.text += "\nSTEPS                  =%15i  FORCE CONSTANT         =%15.3f"  %  ( parameters['nsteps_RC2']       , parameters['force_constant_2']   )
                        self.text += "\nDMINIMUM               =%15.5f  DINCREMENT             =%15.5f" % ( parameters['dminimum_RC2']     , parameters['dincre_RC2']         )
                        self.text += "\n--------------------------------------------------------------------------------"                              #
                    self.text += "\n"                  
                   
                #----------------------------------------------------------------------------------------------------------------------------------#
    #===============================================            
    def add_pdynamo_summary_text (self):
        """ Function doc """
        pass

    #===============================================    
    def create_header_text (self, paramters = None):
        """ Function doc """
        self.text += HEADER

    
    def save_logfile (self, filename = 'logfile', path = None ):
        """ Function doc """
        filename = filename+".log"
        logfile = open( os.path.join(path, filename), "w" ) 
        dprint(os.path.join(path, filename))
        logfile.write(self.text)
        logfile.close()
    
    
    def add_text_Line(self, text_line):
        '''
        Insert lines in the text container.
        '''
        self.text += text_line 
        self.text +="\n"

    #======================================================================
    def separator(self):
        '''
        Include a separator in the log Text.
        '''
        self.text += "===================================================\n"

    #======================================================================
    def close(self):
        '''
        Write and close the file object.
        '''
        #--------------------------------------------------------
        self.end = timer()
        cputime = self.end - self.start
        dprint("Cpu time: " + str(cputime) )
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
        #--------------------------------------------------------
        self.separator()
        self.text += "Finishing at: " + dt_string + "\n"
        self.text += "Elapsed time: " + str(cputime) + "\n"
        self.separator()
        #--------------------------------------------------------
        self.fileObj.write(self.text)
        self.fileObj.close()
    
    #---------------------------------------------------------
    #def get_log(self):
    #    '''
    #    Class object to return a TextLogFileWriter pDynamo instance to use in individual methods
    #    '''
    #    logObj = TextLogFileWriter.WithOptions(self.filePath)
    #    return(logObj)
  
#*******************************************************************************
class LogFileReader:
    """Class for reading and parsing EasyHybrid log files."""

    def __init__(self, logfile):
        """Initialize LogFileReader with the path to the logfile."""

        self.type = None  # Type of log file (to be determined)

        # Extract base name (file name only) and directory path of the logfile
        self.logfile  = logfile
        self.basename = os.path.basename(logfile)
        self.dirname  = os.path.dirname(logfile)

        # Open file and read all lines
        data = open(logfile, "r")
        self.data = data.readlines()
        # Suggestion: use a context manager (with open(...)) instead of manually closing

        # Identify the log type from the file content
        self.get_logtype()
        # self.get_data()  # Not needed here, parsing can be done on demand

        data.close()

    #===================================================================
    def get_data(self):
        """Parse the log file depending on its type and return structured data."""

        if self.type == 'EasyHybrid-SCAN2D':
            '''
            Example of EasyHybrid-SCAN2D log file block:
            (Coordinates and scan parameters for 2D scans)
             
            ----------------------- Coordinate 1 - Simple-Distance -------------------------
            ATOM1                  =           3887  ATOM NAME1             =             O1
            ATOM2                  =           3890  ATOM NAME2             =              N
            NUMBER OF STEPS        =             25  FORCE CONSTANT         =           4000
            DMINIMUM               =        2.30677  MAX INTERACTIONS       =           6000
            STEP SIZE              =     -0.0500000  RMS GRAD               =      0.1000000
            --------------------------------------------------------------------------------

            ----------------------- Coordinate 2 - Simple-Distance -------------------------
            ATOM1                  =           3841  ATOM NAME1             =             Fe
            ATOM2                  =           3887  ATOM NAME2             =             O1
            NUMBER OF STEPS        =            250  FORCE CONSTANT         =           4000
            DMINIMUM               =        1.53720  MAX INTERACTIONS       =           6000
            STEP SIZE              =      0.0500000  RMS GRAD               =      0.1000000
            --------------------------------------------------------------------------------
            '''
            rc1_atoms = []  # Atoms involved in reaction coordinate 1
            rc2_atoms = []  # Atoms involved in reaction coordinate 2

            datalines = []
            _is_rc1 = True
            for line in self.data:
                if "DATA" in line:
                    line2 = line.split()
                    if line2[0] == 'DATA':
                        datalines.append(line2[1:])

                if 'Coordinate 2' in line:
                    _is_rc1 = False

                if 'ATOM' in line:
                    line2 = line.split()
                    # Suggestion: improve robustness (check len(line2) before indexing)
                    if line2[1] == '=':
                        if _is_rc1:
                            rc1_atoms.append(line2[2])
                        else:
                            rc2_atoms.append(line2[2])

            # Get grid dimensions from the last "DATA" line
            lastline = datalines[-1]
            x_size = int(lastline[0])
            y_size = int(lastline[1])

            rows = y_size + 1
            cols = x_size + 1

            # Initialize matrices for energy (Z) and reaction coordinates
            Z   = [[0]*cols for _ in range(rows)]
            RC1 = [[0]*cols for _ in range(rows)]
            RC2 = [[0]*cols for _ in range(rows)]

            # Fill matrices with parsed values
            for line in datalines[0:]:
                x = int(line[0])
                y = int(line[1])
                Z[y][x]   = float(line[-1])
                RC1[y][x] = float(line[-3])
                RC2[y][x] = float(line[-2])

            data = {
                'name': self.basename,
                'type': "plot2D",
                'RC1': RC1,
                'RC2': RC2,
                'Z': Z,
                'RC1_indexes': rc1_atoms,
                'RC2_indexes': rc2_atoms
            }
            return data

        elif self.type == 'EasyHybrid-SCAN':
            '''
            Example of EasyHybrid-SCAN log file block:
            (Coordinates and scan parameters for 1D scans)
            ---------------------- Coordinate 1 - multiple-Distance ------------------------
            ATOM1                  =           3841  ATOM NAME1             =             Fe
            ATOM2*                 =           3887  ATOM NAME2             =             O1
            ATOM3                  =           3890  ATOM NAME3             =              N
            NUMBER OF STEPS        =             30  FORCE CONSTANT         =           4000
            DMINIMUM               =        0.60039  MAX INTERACTIONS       =           6000
            STEP SIZE              =      0.0500000  RMS GRAD               =      0.1000000
            Sigma atom1 - atom3    =        0.79949  Sigma atom3 - atom1    =       -0.20051
            --------------------------------------------------------------------------------
            '''
            datalines = []
            rc_atoms  = []
            n = 0
            for line in self.data:
                if "DATA" in line:
                    line2 = line.split()
                    if line2[0] == "DATA":
                        datalines.append(line2[1:])
                        n += 1

                if 'ATOM' in line:
                    line2 = line.split()
                    if line2[1] == '=':
                        rc_atoms.append(line2[2])

            lastline = datalines[-1]
            x_size   = n  # Number of data points in the scan

            Z   = []
            RC1 = []
            RC2 = []  # Declared but unused, consider removing

            for line in datalines:
                Z.append(float(line[-1]))
                RC1.append(float(line[1]))

            data = {
                'name': self.basename,
                'type': "plot1D",
                'RC1': RC1,
                'Z': Z,
                'RC1_indexes': rc_atoms
            }
            return data

        elif self.type == 'Chain-Of-States':
            '''
            Chain-of-states optimization output parsing.
            Extracts energy (Z) and reaction coordinate (RC1).
            '''
            Z   = []
            RC1 = []

            data2   = self.data
            counter = 0
            start   = 0

            # Find where the "Path Summary" starts
            for line in data2:
                if "Path Summary" in line:
                    start = counter
                counter += 1

            # Extract values after "Path Summary"
            for line in data2[start:]:
                line2 = line.split()
                if len(line2) >= 9:
                    try:
                        Z.append(float(line2[1]))
                        RC1.append(float(line2[-1]))
                    except ValueError:
                        dprint('Logfile parsing: line is not valid data')

            data = {
                'name': self.basename,
                'type': "plot1D",
                'RC1': RC1,
                'Z': Z
            }
            return data

        elif self.type == 'Conjugate-Peak-Refinement':
            '''
            addOns.pyCPR output parsing (see pdynamo/p_methods/
            conjugate_peak_refinement.py and util/cpr_log_parser.py,
            which this reuses instead of re-implementing the same table
            regexes here). Plots the LAST "Path Summary after CPRrun N"
            table -- the final, most-refined reaction path -- as a
            plot1D: RC1 = image index along the path, Z = energy. Same
            shape as 'Chain-Of-States' so PES_analysis_window.py needs
            no changes to display it.
            '''
            from util.cpr_log_parser import parse_cpr_log

            parsed = parse_cpr_log(self.logfile)
            finalPath = parsed['final_path']
            if finalPath is None:
                return None

            data = {
                'name': self.basename,
                'type': "plot1D",
                'RC1' : [float(image) for image in finalPath['image']],
                'Z'   : finalPath['energy'],
            }
            return data

        else:
            # Suggestion: raise an exception instead of silently passing
            return None

    #==============================================================
    def get_logtype(self):
        """Determine the type of EasyHybrid log file based on its content."""

        for line in self.data:
            if "EasyHybrid-SCAN2D" in line.split():
                self.type = 'EasyHybrid-SCAN2D'
            elif "EasyHybrid-SCAN" in line.split():
                self.type = 'EasyHybrid-SCAN'
            elif 'EasyHybrid Energy Refinement 2D' in line:
                # p_methods/energy.py's EnergyRefinement.write_header() -- re-scores
                # an existing scan's frames with a different (e.g. higher-level) method.
                # Same TYPE/Coordinate/DATA line layout as EasyHybrid-SCAN2D, so it is
                # parsed by that same branch below.
                self.type = 'EasyHybrid-SCAN2D'
            elif 'EasyHybrid Energy Refinement' in line:
                self.type = 'EasyHybrid-SCAN'
            elif 'Summary of Chain-Of-States Optimizer' in line:
                self.type = 'Chain-Of-States'
            elif 'Python-based Conjugate Peak Refinement' in line:
                self.type = 'Conjugate-Peak-Refinement'
            # Suggestion: consider using `break` after finding the type for efficiency



  
#*******************************************************************************
class LogFileReader_old:
    """ Class doc """
    
    def __init__ (self, logfile):
        """ Class initialiser """
        
        self.type = None
        
        self.basename = os.path.basename(logfile)
        self.dirname  = os.path.dirname(logfile)
        
        data = open(logfile, "r")
        self.data = data.readlines()
        #print(self.data)
        
        #---------------------------------------------------------------
        self.get_logtype() 
        #self.get_data() 
        data.close()
    
    #===================================================================
    def get_data (self):
        """ Function doc """
        #print(self.type)
        if self.type == 'EasyHybrid-SCAN2D':
            '''
            ----------------------- Coordinate 1 - Simple-Distance -------------------------
            ATOM1                  =           3887  ATOM NAME1             =             O1
            ATOM2                  =           3890  ATOM NAME2             =              N
            NUMBER OF STEPS        =             25  FORCE CONSTANT         =           4000
            DMINIMUM               =        2.30677  MAX INTERACTIONS       =           6000
            STEP SIZE              =     -0.0500000  RMS GRAD               =      0.1000000
            --------------------------------------------------------------------------------

            ----------------------- Coordinate 2 - Simple-Distance -------------------------
            ATOM1                  =           3841  ATOM NAME1             =             Fe
            ATOM2                  =           3887  ATOM NAME2             =             O1
            NUMBER OF STEPS        =            250  FORCE CONSTANT         =           4000
            DMINIMUM               =        1.53720  MAX INTERACTIONS       =           6000
            STEP SIZE              =      0.0500000  RMS GRAD               =      0.1000000
            --------------------------------------------------------------------------------
            '''
            rc1_atoms = []
            rc2_atoms = []
            
            datalines = []
            _is_rc1 = True
            for line in self.data:
                if "DATA" in line:
                    line2 = line.split()
                    if line2[0] == 'DATA':
                        datalines.append(line2[1:])
                
                if 'Coordinate 2' in line:
                    _is_rc1 = False
                
                if 'ATOM' in line:
                    line2 = line.split()
                    
                    if line2[1] == '=':
                        if _is_rc1:
                            rc1_atoms.append(line2[2])
                        else:
                            rc2_atoms.append(line2[2])

            #print(datalines)
            
            lastline = datalines[-1]
            x_size = int(lastline[0])
            y_size = int(lastline[1])

            
            rows = y_size+1
            cols = x_size+1
             
            Z       = [[0]*cols for _ in range(rows)]
            RC1     = [[0]*cols for _ in range(rows)]
            RC2     = [[0]*cols for _ in range(rows)]
            
            for line in datalines[0:]:
                #line2 = line.split()
                x = int(line[0])
                y = int(line[1])
                #print(x,y, line2[-1])
                Z[y][x]       = float(line[-1]) 
                RC1[y][x]     = float(line[-3]) 
                RC2[y][x]     = float(line[-2]) 
            

            data = {
                   'name': self.basename,
                   'type': "plot2D",
                   'RC1' : RC1,
                   'RC2' : RC2,
                   'Z'   : Z
                   }
            
            data['RC1_indexes'] = rc1_atoms
            data['RC2_indexes'] = rc2_atoms
            
            #print (data)
            return data       
        
        
        elif self.type == 'EasyHybrid-SCAN':
            '''
            ---------------------- Coordinate 1 - multiple-Distance ------------------------
            ATOM1                  =           3841  ATOM NAME1             =             Fe
            ATOM2*                 =           3887  ATOM NAME2             =             O1
            ATOM3                  =           3890  ATOM NAME3             =              N
            NUMBER OF STEPS        =             30  FORCE CONSTANT         =           4000
            DMINIMUM               =        0.60039  MAX INTERACTIONS       =           6000
            STEP SIZE              =      0.0500000  RMS GRAD               =      0.1000000
            Sigma atom1 - atom3    =        0.79949  Sigma atom3 - atom1    =       -0.20051
            --------------------------------------------------------------------------------
            '''
            
            
            datalines = []
            rc_atoms  = [] 
            n = 0
            for line in self.data:
                if "DATA" in line:
                    #print(line)
                    
                    line2 = line.split()
                    if line2[0] == "DATA":
                        datalines.append(line2[1:])
                        n += 1
                
                if 'ATOM' in line:
                    line2 = line.split()
                    
                    if line2[1] == '=':
                        rc_atoms.append(line2[2])
                    
                    
            lastline = datalines[-1]
            x_size   = n
             
            Z       = []
            RC1     = []
            RC2     = []
            
            for line in datalines:
                dprint(line)
                #if line[0] == 'DATA':
                Z.append(float(line[-1])) 
                RC1.append(float(line[1]))
 
                    
            data = {
                   'name': self.basename,
                   'type': "plot1D",
                   'RC1' : RC1,
                   'Z'   : Z
                   }
            data['RC1_indexes'] = rc_atoms
            #print(data)
            #print(rc_atoms)
            return data        
        
        
        elif self.type == 'Chain-Of-States':
            Z       = []
            RC1     = []
            RC2     = []
            
            datalines = []
            
            data2   = self.data
            counter = 0
            start   = 0
            
            
            
            for line in data2:
                if "Path Summary" in line:
                    start = counter
                counter += 1
            
            
            for line in data2[start:]: 
                line2  = line.split()
                if len(line2) >= 9:
                    try:
                        Z.append(float(line2[1]))
                        RC1.append(float(line2[-1]))
                    except:
                        dprint('Logfile parsing. Line is not a valid data')
            data = {
                   'name': self.basename,
                   'type': "plot1D",
                   'RC1' : RC1,
                   'Z'   : Z
                   }
            #print(data)
            return data  
        
        
        
        else:
            pass
        
    #==============================================================
    def get_logtype (self):
        """ Function doc """
        
        for line in self.data:
            #print(line)
            if "EasyHybrid-SCAN2D" in line.split():
                #line2 = line.split()
                #print(line)
                self.type = 'EasyHybrid-SCAN2D'
            
            elif "EasyHybrid-SCAN" in line.split():
                #line2 = line.split()
                #print(line)
                self.type = 'EasyHybrid-SCAN'
            elif 'Summary of Chain-Of-States Optimizer' in line:
                #print(line)
                self.type = 'Chain-Of-States'
            
            else:
                pass
#==================================================================


# Matches an "ATOMn = idx  ATOM NAMEn = name" pair on a single line, in ANY
# of the 3 shapes EasyHybrid's own scan/refinement backends write it in
# (see detect_reaction_coordinates_from_log below for which):
#   "ATOM1                  =           3841  ATOM NAME1             =             Fe"
#   "ATOM2*                 =           3887  ATOM NAME2             =             O1"   (multiple_distance's middle atom)
#   "ATOM                   =           3841  ATOM NAME              =             Fe"   (advanced/weighted-list, unnumbered)
_ATOM_LINE_RE = re.compile(r'ATOM\d*\*?\s*=\s*(\S+)\s+ATOM NAME\d*\s*=\s*(\S+)')

# Matches the per-pair weight line the advanced/weighted-list format writes
# right after each pair's two ATOM lines -- "Sigma" (surface_scan.py's
# AdvancedRelaxedSurfaceScan) or "WEIGHT" (energy.py's EnergyRefinement).
# Deliberately does NOT match RelaxedSurfaceScan's own informational
# "Sigma atom1 - atom3    =...Sigma atom3 - atom1    =..." line (that text
# sits between "Sigma" and "=", so \s*= alone doesn't match it) -- that line
# only ever appears for the FIXED-SHAPE 'multiple_distance' RC, which must
# NOT be misdetected as advanced.
_WEIGHT_LINE_RE = re.compile(r'(?:Sigma|WEIGHT)\s*=\s*(\S+)')


def detect_reaction_coordinates_from_log(log_path):
    """ Scans an EasyHybrid-written 'output.log' (from a Relaxed Surface
    Scan, Advanced Relaxed Surface Scan, or Energy Refinement run -- see
    p_methods/surface_scan.py and p_methods/energy.py) for the reaction
    coordinate(s) it used, so a window that consumes an EXISTING
    trajectory (e.g. Umbrella Sampling's "From Trajectory" input mode)
    can auto-fill its own RC boxes instead of making the user retype the
    same atoms that produced that trajectory in the first place.

    Returns {'RC1': dict_or_None, 'RC2': dict_or_None}, or None if the
    file doesn't exist / isn't a recognizable EasyHybrid scan log.

    Each RC dict is either:
        - fixed-shape: {'rc_type': 'simple_distance'|'multiple_distance'|
                        'multiple_distance*4atoms',
                        'ATOMS': [int, ...], 'ATOM_NAMES': [str, ...]}
          -- directly usable with ReactionCoordinateBox.set_rc_data().
        - advanced:    {'rc_type': 'advanced',
                        'RC': [[name1, idx1, name2, idx2, weight, '0.0'], ...]}
          -- same row shape AdvancedReactionCoordinateBox's own treeview
          uses. '0.0' is a placeholder for the display-only "dist" column:
          it is never read back by compute_reaction_coordinate() or
          RestraintMultipleDistance, both of which always recompute the
          live distance from the current coordinates, not this stored
          value.
    """
    if not os.path.isfile(log_path):
        return None

    try:
        with open(log_path, 'r', errors='replace') as f:
            text = f.read()
    except Exception:
        return None

    if 'Coordinate 1' not in text:
        return None

    def _parse_section(section_text):
        atom_pairs = _ATOM_LINE_RE.findall(section_text)    # [(idx, name), ...]
        weights    = _WEIGHT_LINE_RE.findall(section_text)  # [weight, ...]

        if not atom_pairs:
            return None

        if weights:
            # Advanced (weighted distance list): each weight consumes the
            # next 2 atoms in order.
            rows = []
            for k, w in enumerate(weights):
                if 2 * k + 1 >= len(atom_pairs):
                    break
                idx1, name1 = atom_pairs[2 * k]
                idx2, name2 = atom_pairs[2 * k + 1]
                rows.append([name1, idx1, name2, idx2, w, '0.0'])
            if not rows:
                return None
            return {'rc_type': 'advanced', 'RC': rows}

        # Fixed-shape: 2/3/4 atoms -> simple_distance/multiple_distance/
        # multiple_distance*4atoms (matches how the writers themselves
        # decide which shape to use -- see surface_scan.py/energy.py).
        n = len(atom_pairs)
        if n == 2:
            rc_type = 'simple_distance'
        elif n == 3:
            rc_type = 'multiple_distance'
        elif n >= 4:
            rc_type = 'multiple_distance*4atoms'
            atom_pairs = atom_pairs[:4]
        else:
            return None

        atoms      = [int(idx) for idx, _name in atom_pairs]
        atom_names = [name for _idx, name in atom_pairs]
        return {'rc_type': rc_type, 'ATOMS': atoms, 'ATOM_NAMES': atom_names}

    idx_c1 = text.index('Coordinate 1')
    if 'Coordinate 2' in text:
        idx_c2 = text.index('Coordinate 2')
        rc1_section = text[idx_c1:idx_c2]
        rc2_section = text[idx_c2:]
    else:
        rc1_section = text[idx_c1:]
        rc2_section = None

    rc1 = _parse_section(rc1_section)
    rc2 = _parse_section(rc2_section) if rc2_section else None

    if rc1 is None and rc2 is None:
        return None

    return {'RC1': rc1, 'RC2': rc2}









