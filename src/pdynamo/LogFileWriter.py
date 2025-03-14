#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#FILE = LofFile.py

##############################################################
#-----------------...EasyHybrid 3.0...-----------------------#
#-----------Credits and other information here---------------#
##############################################################

from pCore import *
from datetime import datetime
from timeit import default_timer as timer
import os
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
        print(os.path.join(path, filename))
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
        print("Cpu time: " + str(cputime) )
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
            datalines = []
            
            for line in self.data:
                if "DATA" in line:
                    line2 = line.split()
                    datalines.append(line2[1:])
            print(datalines)
            
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
            #print(data)
            return data       
        
        
        elif self.type == 'EasyHybrid-SCAN':
            datalines = []
            
            n = 0
            for line in self.data:
                if "DATA" in line:
                    line2 = line.split()
                    datalines.append(line2[1:])
                    n += 1
            lastline = datalines[-1]
            x_size   = n
             
            Z       = []
            RC1     = []
            RC2     = []
            
            for line in datalines:
                Z.append(float(line[-1])) 
                RC1.append(float(line[1])) 
            data = {
                   'name': self.basename,
                   'type': "plot1D",
                   'RC1' : RC1,
                   'Z'   : Z
                   }
            print(data)
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
                        print('Logfile parsing. Line is not a valid data')
            data = {
                   'name': self.basename,
                   'type': "plot1D",
                   'RC1' : RC1,
                   'Z'   : Z
                   }
            print(data)
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
                print(line)
                self.type = 'EasyHybrid-SCAN2D'
            
            elif "EasyHybrid-SCAN" in line.split():
                #line2 = line.split()
                print(line)
                self.type = 'EasyHybrid-SCAN'
            elif 'Summary of Chain-Of-States Optimizer' in line:
                print(line)
                self.type = 'Chain-Of-States'
            
            else:
                pass
#==================================================================













