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
from pdynamo.p_methods._common import backup_orca_files, write_header

class LogFile:
    """ Class doc """
    
    def __init__ (self, system):
        """ Class initialiser """
        self.path = os.path.join(os.environ.get('PDYNAMO3_SCRATCH'), 'summary_temp.log')
        self.logFile2 = TextLogFileWriter.WithOptions ( path = self.path )
        system.Summary(log = self.logFile2)
        self.logFile2.Close()


class EnergyCalculation:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass
    
    def run (self, parameters):
        """ Function doc """
        full_path_file = os.path.join(parameters['folder'])
        self.logFile2  = TextLogFileWriter.WithOptions ( path = os.path.join(full_path_file, parameters['filename']+'.log') )
        
        parameters['system'].Summary(log = self.logFile2)
        #try:
        energy = parameters['system'].Energy(log = self.logFile2)
        
        '''
        #---------------------------------------------------------------
        print_MM   = True # print MM charges
        print_QC   = True # print QC charges (does not inclue boundary atoms)
        print_QCMM = True # print QC/MM 
        try:
            charges = list(parameters['system'].AtomicCharges())
            #except:
            #    charges = []
            MM_charges = list(parameters['system'].mmState.charges)
            
             
            #atomTypes  = parameters['system'].mmState.atomTypes
            if getattr( parameters['system'], 'qcState', False):
                qcAtoms = list(parameters['system'].qcState.qcAtoms)
                e_qc_table = list(parameters['system'].qcState.pureQCAtoms)
            else:
                qcAtoms =[]
                e_qc_table = []
            print(qcAtoms, e_qc_table)
            print(len(e_qc_table), len(qcAtoms), len(parameters['system'].atoms), len(charges)) 
            print('---------------------------------------------')
            print(' Index  Atom  Type   Charge QC/MM  Charge MM ')
            print('---------------------------------------------')
    
            for index, atom in enumerate(parameters['system'].atoms):
                #{:<3s} {:20.10f} 1 {:20.10f}
                if index in qcAtoms:
                    if index in e_qc_table:
                        if print_QC: 
                            line = ' {:<3d}    {:<5s} QC   {:10.6f}  {:10.6f}'.format(index ,
                                                                                    atom.label, 
                                                                                    charges[index], 
                                                                                    MM_charges[index])
                    else:
                        if print_QCMM:
                            line = ' {:<3d}    {:<5s} QC*  {:10.6f}  {:10.6f}'.format(index ,atom.label, charges[index], MM_charges[index])
    
                else:
                    if print_MM:
                        line = ' {:<3d}    {:<5s} MM   {:10.6f}  {:10.6f}'.format(index ,atom.label, charges[index], MM_charges[index])
                print(line)
            print('---------------------------------------------')
        except:
            print('Ops!')
            #print(index ,atom.label, charges[index])
        #---------------------------------------------------------------
        #'''
        
        #except :
        #    msg = 'Error!'
        #    return False, msg
        
        backup_orca_files(system        = parameters['system']         , 
                          output_folder = parameters['folder']         , 
                          output_name   = parameters['filename'])
        
        self.logFile2.Footer ( )
        self.logFile2.Close()
        self.logFile2 = None
        
        return energy, 'Energy: '+str(energy)


class EnergyRefinement:
    
    def __init__ (self):
        """ Class initialiser """
        pass
    
    def run (self, parameters, interface = False):
        
        full_path_trajectory = os.path.join(parameters['folder'], 
                               parameters['filename'])
        os.mkdir(
                 full_path_trajectory
                 )
        
        # - - - - - - - - - - - - - Checking trajectory - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        self.logFile2 = TextLogFileWriter.WithOptions ( path = os.path.join(full_path_trajectory, 'output.log') )
        parameters['system'].Summary(log = self.logFile2)
        self.logFile2.Header ( )
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        self.logFile2.Footer ( )
        self.logFile2.Close()
        self.logFile2 = None
        
        logfile = self.write_header (parameters = parameters,
                                     logfile    = os.path.join(full_path_trajectory, 'output.log') )
        
        
        
        #---------------------------------------------------------------
        backup = []
        try:
            backup.append(parameters['system'].e_treeview_iter)
            backup.append(parameters['system'].e_liststore_iter)
            parameters['system'].e_treeview_iter   = None
            parameters['system'].e_liststore_iter  = None
        except:
            pass
        
        sys  =  Clone(parameters['system'])
        
        
        try:
            parameters['system'].e_treeview_iter   = backup[0]
            parameters['system'].e_liststore_iter  = backup[1]
        except:
            pass
        #---------------------------------------------------------------
        
        parameters['system'] = sys
        
        #-------------------------------------------------------------------------------
        original_charges = None
        if parameters['ignore_mm_charges']:
            print('Adjusting electrical charges in the MM region to zero.')
            
            original_charges = parameters['system'].e_charges_backup.copy()
            
            for index, charge in enumerate(original_charges):
                parameters['system'].mmState.charges[index] = 0.0
            print(parameters['system'].mmState.charges)
            
            for i, charge in enumerate(parameters['system'].mmState.charges):
                parameters['system'].mmState.charges[i] = 0.0
            print(list(parameters['system'].mmState.charges))
        else:
            pass
        #-------------------------------------------------------------------------------
        
        
        if parameters['traj_type'] == 'pklfolder':
            data_energy = {}
            energy_list = []
            text = ""
            for frame in parameters["trajectory"]:
                x = frame[5:-4]
                parameters['system'].coordinates3 = ImportCoordinates3 (os.path.join(parameters['data_path'],frame)) 
                energy = parameters['system'].Energy()
                energy_list.append(energy)
                
                if parameters['RC1']['rc_type'] == 'simple_distance':
                    atom1 = parameters['RC1']['ATOMS'][0]
                    atom2 = parameters['RC1']['ATOMS'][1]
                    dist =  parameters['system'].coordinates3.Distance (atom1, atom2) 
                    data_energy[int(x)] = [energy, dist ]
                
                elif parameters['RC1']['rc_type'] == 'multiple_distance':
                    atom1 = parameters['RC1']['ATOMS'][0]
                    atom2 = parameters['RC1']['ATOMS'][1]
                    atom3 = parameters['RC1']['ATOMS'][2]
                    dist1 =  parameters['system'].coordinates3.Distance (atom1, atom2) 
                    dist2 =  parameters['system'].coordinates3.Distance (atom2, atom3) 
                    dist  =  dist1 - dist2 
                    data_energy[int(x)] = [energy, dist1, dist2]
                
                elif parameters['RC1']['rc_type'] == 'multiple_distance*4atoms':
                    atom1 = parameters['RC1']['ATOMS'][0]
                    atom2 = parameters['RC1']['ATOMS'][1]
                    atom3 = parameters['RC1']['ATOMS'][2]
                    atom4 = parameters['RC1']['ATOMS'][3]
                    dist1 =  parameters['system'].coordinates3.Distance (atom1, atom2) 
                    dist2 =  parameters['system'].coordinates3.Distance (atom3, atom4) 
                    dist  =  dist1 - dist2 
                    data_energy[int(x)] = [energy, dist1, dist2]
                
                print ('frame:', x, dist,  energy)
                backup_orca_files(system        = parameters['system']         , 
                                  output_folder = full_path_trajectory         , 
                                  output_name   = 'frame'+x)
            
            lowest_energy = min(energy_list)
            keys = data_energy.keys()
            
            for i in range(0, len(keys)):

                
                if parameters['RC1']['rc_type'] == 'simple_distance':
                    distance = data_energy[i][1]
                    energy   = data_energy[i][0]
                    text += "\nDATA %9i       %13.12f        %13.12f"% (int(i), float(distance), float(energy))
                
                elif parameters['RC1']['rc_type'] == 'multiple_distance':
                    distance1 = data_energy[i][1]
                    distance2 = data_energy[i][2]
                    energy    = data_energy[i][0]
                    text += "\nDATA %9i       %13.12f        %13.12f        %13.12f"% (int(i), float(distance1), float(distance2), float(energy))
                
                elif parameters['RC1']['rc_type'] == 'multiple_distance*4atoms':
                    distance1 = data_energy[i][1]
                    distance2 = data_energy[i][2]
                    energy    = data_energy[i][0]
                    text += "\nDATA %9i       %13.12f        %13.12f        %13.12f"% (int(i), float(distance1), float(distance2), float(energy))
                
                else:
                    pass
        
        
            logfile.write(text)
        
        
        elif parameters['traj_type'] == 'pklfolder2D': 
            max_i = 0
            max_j = 0
            i_list = []
            j_list = []
            data   = {}
            energy_list = []
            
            text = ""
            for frame in parameters["trajectory"]:
                i_j = frame[5:-4].split('_')
                parameters['system'].coordinates3 = ImportCoordinates3 (os.path.join(parameters['data_path'],frame)) 
                energy = parameters['system'].Energy()
                i = int(i_j[0])
                j = int(i_j[1])
                i_list.append(i)
                j_list.append(j)
                energy_list.append(energy)
                
                
                
                if parameters['RC1']['rc_type'] == 'simple_distance':
                    atom1 = parameters['RC1']['ATOMS'][0]
                    atom2 = parameters['RC1']['ATOMS'][1]
                    
                    dist1 =  parameters['system'].coordinates3.Distance (atom1, atom2) 
 
                
                elif parameters['RC1']['rc_type'] == 'multiple_distance':
                    atom1 = parameters['RC1']['ATOMS'][0]
                    atom2 = parameters['RC1']['ATOMS'][1]
                    atom3 = parameters['RC1']['ATOMS'][2]
                    dist_RC1_1 =  parameters['system'].coordinates3.Distance (atom1, atom2) 
                    dist_RC1_2 =  parameters['system'].coordinates3.Distance (atom2, atom3) 
                    
                    dist1 =  didist_RC1_1 - dist_RC1_2
                 
                elif parameters['RC1']['rc_type'] == 'multiple_distance*4atoms':
                    atom1 = parameters['RC1']['ATOMS'][0]
                    atom2 = parameters['RC1']['ATOMS'][1]
                    atom3 = parameters['RC1']['ATOMS'][2]
                    atom4 = parameters['RC1']['ATOMS'][3]
                    dist_RC1_1 =  parameters['system'].coordinates3.Distance (atom1, atom2) 
                    dist_RC1_2 =  parameters['system'].coordinates3.Distance (atom3, atom4) 
                    dist1  =  dist_RC1_1 - dist_RC1_2 
                 
                
                if parameters['RC2']['rc_type'] == 'simple_distance':
                    atom1 = parameters['RC2']['ATOMS'][0]
                    atom2 = parameters['RC2']['ATOMS'][1]
                    dist2 = parameters['system'].coordinates3.Distance (atom1, atom2) 
                 
                elif parameters['RC2']['rc_type'] == 'multiple_distance':
                    atom1 = parameters['RC2']['ATOMS'][0]
                    atom2 = parameters['RC2']['ATOMS'][1]
                    atom3 = parameters['RC2']['ATOMS'][2]
                    dist_RC2_1 =  parameters['system'].coordinates3.Distance (atom1, atom2) 
                    dist_RC2_2 =  parameters['system'].coordinates3.Distance (atom2, atom3) 
                    
                    dist2 = dist_RC2_1 - dist_RC2_2
                 
                elif parameters['RC2']['rc_type'] == 'multiple_distance*4atoms':
                    atom1 = parameters['RC2']['ATOMS'][0]
                    atom2 = parameters['RC2']['ATOMS'][1]
                    atom3 = parameters['RC2']['ATOMS'][2]
                    atom4 = parameters['RC2']['ATOMS'][3]
                    dist_RC2_1 =  parameters['system'].coordinates3.Distance (atom1, atom2) 
                    dist_RC2_2 =  parameters['system'].coordinates3.Distance (atom3, atom4) 
                    dist2 = dist_RC2_1 - dist_RC2_2
            
                 
                 
                
                data[(i,j)] = [dist1, dist2, energy]
                print ('frame:', i_j, energy)

            max_i = max(i_list)
            max_j = max(j_list)
            
            for i in range(0, max_i+1):
                for j in range(0, max_j+1):

                    text = "\nDATA  %4i  %4i     %13.12f       %13.12f       %13.12f"% (int(i), int(j),  float(data[(i,j)][0]), float(data[(i,j)][1]), float(data[(i,j)][2]))
                    logfile.write(text)

        # not necessary anymore because system now is a clone
        #-------------------------------------------------------------------------------
        #if parameters['ignore_mm_charges']:
        #    print('Restoring electrical charges in the MM.')
        #    for index, charge in enumerate(original_charges):
        #        parameters['system'].mmState.charges[index]   = original_charges[index]
        #               
        #    for i, charge in enumerate(parameters['system'].mmState.charges):
        #        parameters['system'].mmState.charges[i] =  original_charges[i]
        #    print(list(parameters['system'].mmState.charges))
        #
        #
        #
        #else:
        #    pass
        #-------------------------------------------------------------------------------
        
        

    def write_header (self, parameters, logfile = 'output.log'):
        """ Function doc """
        
        arq = open(logfile, "a")
        text = ""

        if parameters['RC2'] is not None :
            text = text + "\n"
            text = text + "\n--------------------------------------------------------------------------------"
            text = text + "\nTYPE                 EasyHybrid Energy Refinement 2D                            "
            text = text + "\n--------------------------------------------------------------------------------"

        else:
            text = text + "\n"
            text = text + "\n--------------------------------------------------------------------------------"
            text = text + "\nTYPE                   EasyHybrid Energy Refinement                             "
            text = text + "\n--------------------------------------------------------------------------------"




        '''This part writes the parameters used in the first reaction coordinate'''
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------
        if parameters['RC1']["rc_type"] == 'simple_distance':
            text = text + "\n"
            text = text + "\n----------------------- Coordinate 1 - Simple-Distance -------------------------"
            text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC1']['ATOMS'][0], parameters['RC1']['ATOM_NAMES'][0] )
            text = text + "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['RC1']['ATOMS'][1], parameters['RC1']['ATOM_NAMES'][1] )
            text = text + "\n--------------------------------------------------------------------------------"

        
        elif parameters['RC1']["rc_type"] == 'multiple_distance':
            text = text + "\n"
            text = text + "\n---------------------- Coordinate 1 - multiple-Distance ------------------------"	
            text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC1']['ATOMS'][0]    , parameters['RC1']['ATOM_NAMES'][0] )
            text = text + "\nATOM2*                 =%15i  ATOM NAME2             =%15s"     % (parameters['RC1']['ATOMS'][1]    , parameters['RC1']['ATOM_NAMES'][1] )
            text = text + "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC1']['ATOMS'][2]    , parameters['RC1']['ATOM_NAMES'][2] )
            text = text + "\n--------------------------------------------------------------------------------"
        
        
        elif parameters['RC1']["rc_type"] == 'multiple_distance*4atoms':
            text = text + "\n"
            text = text + "\n---------------------- Coordinate 1 - multiple-Distance ------------------------"	
            text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC1']['ATOMS'][0]    , parameters['RC1']['ATOM_NAMES'][0] )
            text = text + "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['RC1']['ATOMS'][1]    , parameters['RC1']['ATOM_NAMES'][1] )
            text = text + "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC1']['ATOMS'][2]    , parameters['RC1']['ATOM_NAMES'][2] )
            text = text + "\nATOM4                  =%15i  ATOM NAME4             =%15s"     % (parameters['RC1']['ATOMS'][3]    , parameters['RC1']['ATOM_NAMES'][3] )
            text = text + "\n--------------------------------------------------------------------------------"
        
        
        
        else:
            pass
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------
        
        
        
        
        if parameters['RC2'] is not None :

            '''This part writes the parameters used in the second reaction coordinate'''
            #-----------------------------------------------------------------------------------------------------------------------------------------------------------
            if parameters['RC2']["rc_type"] == 'simple_distance':
                text = text + "\n"
                text = text + "\n----------------------- Coordinate 2 - Simple-Distance -------------------------"
                text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC2']['ATOMS'][0], parameters['RC2']['ATOM_NAMES'][0] )
                text = text + "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['RC2']['ATOMS'][1], parameters['RC2']['ATOM_NAMES'][1] )
                text = text + "\n--------------------------------------------------------------------------------"

            
            elif parameters['RC2']["rc_type"] == 'multiple_distance':
                text = text + "\n"
                text = text + "\n---------------------- Coordinate 2 - multiple-Distance ------------------------"	
                text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC2']['ATOMS'][0]    , parameters['RC2']['ATOM_NAMES'][0] )
                text = text + "\nATOM2*                 =%15i  ATOM NAME2             =%15s"     % (parameters['RC2']['ATOMS'][1]    , parameters['RC2']['ATOM_NAMES'][1] )
                text = text + "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC2']['ATOMS'][2]    , parameters['RC2']['ATOM_NAMES'][2] )
                text = text + "\n--------------------------------------------------------------------------------"
            
            elif parameters['RC2']["rc_type"] == 'multiple_distance*4atoms':
                text = text + "\n"
                text = text + "\n---------------------- Coordinate 1 - multiple-Distance ------------------------"	
                text = text + "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC2']['ATOMS'][0]    , parameters['RC2']['ATOM_NAMES'][0] )
                text = text + "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['RC2']['ATOMS'][1]    , parameters['RC2']['ATOM_NAMES'][1] )
                text = text + "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC2']['ATOMS'][2]    , parameters['RC2']['ATOM_NAMES'][2] )
                text = text + "\nATOM4                  =%15i  ATOM NAME4             =%15s"     % (parameters['RC2']['ATOMS'][3]    , parameters['RC2']['ATOM_NAMES'][3] )
                text = text + "\n--------------------------------------------------------------------------------"
            
                
            else:
                pass
            #-----------------------------------------------------------------------------------------------------------------------------------------------------------
            
        
        
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------
        '''This part writes the header of the frames distances and energy'''
        if parameters['RC2'] is not None :
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
        #-----------------------------------------------------------------------------------------------------------------------------------------------------------


        
        arq.write(text)
        return arq
