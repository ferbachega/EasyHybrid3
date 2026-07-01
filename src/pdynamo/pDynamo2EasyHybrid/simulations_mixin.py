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
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, GLib
import multiprocessing

import glob, math, os, os.path, sys, shutil
import pickle
import threading
from util.file_parser import read_MOL2  
from util.file_parser import read_SIMPLE_txt  
from util.file_parser import read_MOPAC_aux  

from datetime import date
import time

import numpy as np
import copy

import random
import string

from pprint import pprint

#VISMOL_HOME = os.environ.get('VISMOL_HOME')

#path fo the core python files on your machine
#sys.path.append(os.path.join(VISMOL_HOME,"easyhybrid/pDynamoMethods") )
#sys.path.append(os.path.join(VISMOL_HOME,"easyhybrid/gui"))

#from LogFile import LogFileReader

#from gEngine.vismol_object import EVismolObject

#from vismol.model.atom import Atom
from vismol.model.residue import Residue
from vismol.model.chain import Chain
from vismol.core.vismol_object import VismolObject
from vismol.model.atom import Atom
#print ('\n\n\n\n\n\nATOM',Atom,'\n\n\n\n\n\nATOM')
from logging import getLogger
logger = getLogger(__name__)

#---------------------------------------
from pBabel                    import*                                     
from pCore                     import*  
#---------------------------------------
from pMolecule                 import*                              
from pMolecule.MMModel         import*
from pMolecule.NBModel         import*                                     
from pMolecule.QCModel         import*
#---------------------------------------
from pScientific               import*                                     
from pScientific.Arrays        import*                                     
from pScientific.Geometry3     import*                                     
from pScientific.RandomNumbers import*                                     
from pScientific.Statistics    import*
from pScientific.Symmetry      import*
#---------------------------------------                              
from pSimulation               import*
#---------------------------------------


import numpy as np
#from vismol.model.molecular_properties import ATOM_TYPES
from vismol.libgl.representations import DashedLinesRepresentation

from util.colorpalette import CUSTOM_COLOR_PALETTE

from pdynamo.p_methods import GeometryOptimization
from pdynamo.p_methods import RelaxedSurfaceScan
from pdynamo.p_methods import AdvancedRelaxedSurfaceScan
from pdynamo.p_methods import MolecularDynamics
from pdynamo.p_methods import ChainOfStatesOptimizePath
from pdynamo.p_methods import NormalModes
from pdynamo.p_methods import EnergyCalculation
from pdynamo.p_methods import EnergyRefinement
from pdynamo.p_methods import UmbrellaSampling

from pdynamo.p_methods import WHAMAnalysis
from pdynamo.LogFileWriter import LogFileReader

from gui.windows.setup.windows_and_dialogs import call_message_dialog

class pSimulations:
    """Class responsible for managing and running molecular simulations 
    using multiprocessing and GUI integration with EasyHybrid.
    """
    #MSG_RESULT  = "RESULT"
    #MSG_ERROR   = "ERROR"
    #MSG_DONE    = "DONE"
    #MSG_RUNNING = "RUNNING"
    
    def __init__(self):
        """Class initializer."""
        self.process_pool = {}   # Holds active processes and their communication queues
        self.psystem = {}        # Dictionary of systems currently loaded
        self.active_id = None    # ID of the currently active system
        self.changed = False     # Tracks whether the session/project has been modified
        self.main = None         # Reference to the main application object
        self.target_process = None

        # --- constants for messages
        #self.MSG_RESULT  = "RESULT"
        #self.MSG_ERROR   = "ERROR"
        #self.MSG_DONE    = "DONE"
        #self.MSG_RUNNING = "RUNNING"



    # ========================================================================
    # RESTRAINTS and SETUP
    # ========================================================================

    def _apply_restraints(self, parameters):
        """Apply distance restraints to the system if specified in parameters.
        
        parameters['system'].e_restraints_dict[rest_name] = [
            True,                # Active or not
            rest_name,           # Name of the restraint
            'distance',          # Type of restraint
            [atom1, atom2],      # Atoms involved
            distance,            # Target distance
            force_constant       # Force constant
        ]
        """
        # Reset and define a new restraint model
        parameters['system'].DefineRestraintModel(None)
        restraints = RestraintModel()
        parameters['system'].DefineRestraintModel(restraints)

        # Iterate through all defined restraints
        for name, restraint_list in parameters['system'].e_restraints_dict.items():
            if restraint_list[0]:  # Only apply if marked as active
                if restraint_list[2] == 'distance':
                    rmodel = RestraintEnergyModel.Harmonic(restraint_list[4], restraint_list[5])
                    restraint = RestraintDistance.WithOptions(
                        energyModel=rmodel,
                        point1=restraint_list[3][0],
                        point2=restraint_list[3][1]
                    )
                    restraints[name] = restraint
                
                
                if restraint_list[2] == 'position':
                    '''
                    e_restraints_dict[rest_name] = [
                       0 True,                         # active by default
                       1 rest_name,                    # name
                       2 'position',                   # type
                       3 parameters['atomlist']      , # atom list
                       4 parameters['reference']     , # target reference: ( system.coordinates3 ) Dynamic?
                       5 parameters['force_constant'], # force constant
                       6 parameters['system'].e_id     # system id
                       ]    
                        
                    '''
                    #Selection.FromIterable (freeAtoms)
                    print("restraint_list[3]",restraint_list[3])
                    selection          = Selection.FromIterable (restraint_list[3])
                    reference          = Clone (parameters['system'].coordinates3 )
                    tetherEnergyModel  = RestraintEnergyModel.Harmonic ( 0.0, float(restraint_list[5]) )
                    restraint          = RestraintMultipleTether.WithOptions ( energyModel = tetherEnergyModel ,
                                                                               reference   = reference         , 
                                                                               selection   = selection        )

                    restraints[name] = restraint

        return parameters

    def _configure_logfile(self, parameters):
        system    = parameters['system']
        if 'logfile' in parameters:
            system.e_job_history[system.e_step_counter]['logfile'] = parameters['logfile']
            return parameters
        
        #---------------------------------------------------------------
        '''
        here we are setting a temporary name for the log files that 
        will be saved in the scratch directory
        '''
        folder    = os.environ.get('PDYNAMO3_SCRATCH')
        e_id      = str(parameters['system'].e_id)
        key       = system.e_job_history[system.e_step_counter]['key6'] # this is the new key6 - to the object that will be created
        tmp_fname = '{}_{}.log'.format(e_id,key)
        #---------------------------------------------------------------

        
        if 'filename' in parameters:
            parameters['logfile'] = os.path.join(folder, parameters['filename'] + '.log')
            system.e_job_history[system.e_step_counter]['logfile'] = parameters['logfile']
        
        elif parameters.get('trajectory_name'):
            parameters['logfile'] = os.path.join(folder, parameters['trajectory_name'], 'output.log')
            system.e_job_history[system.e_step_counter]['logfile'] = parameters['logfile']
        else:
            parameters['logfile'] = os.path.join(folder, tmp_fname)
            system.e_job_history[system.e_step_counter]['logfile'] = parameters['logfile']
        
        print ('_configure_logfile:', parameters)
        return parameters

    # ========================================================================
    # PROCESS QUEUE HANDLING
    # ========================================================================

    def _check_queue(self):
        """Check all process queues for messages.

        This function is periodically called by GLib.timeout_add to handle 
        inter-process communication asynchronously without blocking the GTK loop.
        """
        for e_id, (queue, process, result, path, treeiter) in list(self.process_pool.items()):
            #print(e_id, (queue, process, result, path, treeiter))
            try:
                # Read all messages available in the queue
                while not queue.empty():
                    msg = queue.get_nowait()
                    #print(msg)
                    self._dispatch_message(msg, e_id, process, path, treeiter)
                    
                    '''
                    if isinstance(msg, tuple) and msg[0] == "RESULT":
                        # Worker process returned results
                        self._handle_result(msg[1])
                    
                    elif isinstance(msg, tuple) and msg[0] == "ERROR":
                        # Worker process returned results
                        results = msg[1]
                        system = self.psystem[results['e_id']]
                        # Update GUI process manager window
                        self.main.process_manager_window.set_status(treeiter, "Error!")
                        self.main.process_manager_window.set_time(treeiter, False, True)
                        self.main.process_manager_window.set_step_counter(treeiter, system.e_step_counter)

                        # Mark process as completed in pool
                        self.process_pool[system.e_id][2] = "Error"

                        # Update bottom notebook GUI
                        iter_ = self.main.bottom_notebook.status_liststore.get_iter(path)
                        value = self.main.bottom_notebook.status_liststore.get_value(iter_, 1)
                        self.main.bottom_notebook.status_liststore.set_value(iter_, 1, value.replace('Running...', 'Error!'))
                        
                        system.e_job_history[system.e_step_counter] = results

                    elif msg == "DONE":
                        # Worker process finished execution
                        self._handle_done(e_id, process, path, treeiter)

                    elif msg == "Running":
                        # Update process status to "Running..."
                        self.main.process_manager_window.set_status(treeiter, "Running...")
                    '''


            except Exception as e:
                # Prevent GUI crash due to queue handling issues
                print(f"Error checking the queue of process {e_id}: {e}")

        # Keep the GLib timeout active
        return True

    def _dispatch_message(self, msg, e_id, process, path, treeiter):
        """ Function doc """
        if isinstance(msg, tuple) and msg[0] ==  self.MSG_RESULT:
            self._handle_result(msg[1])
        elif isinstance(msg, tuple) and msg[0] ==  self.MSG_ERROR:
            self._handle_error(msg[1], path, treeiter)
        elif msg ==  self.MSG_DONE:
            self._handle_done(e_id, process, path, treeiter)
        elif msg ==  self.MSG_RUNNING:
            self.main.process_manager_window.set_status(treeiter, "Running...")
        
    def _handle_error(self, results, path, treeiter):
        """Handle 'ERROR' messages from worker processes."""
        system = self.psystem[results['e_id']]
        step_counter = system.e_step_counter
        
        self.main.process_manager_window.set_status(treeiter, "Error!")
        self.main.process_manager_window.set_time(treeiter, False, True)
        self.main.process_manager_window.set_step_counter(treeiter, system.e_step_counter)
        system.e_job_history[step_counter]['status'] = "Error!"
        
        self.process_pool[system.e_id][2] = "Error"

        iter_ = self.main.bottom_notebook.status_liststore.get_iter(path)
        value = self.main.bottom_notebook.status_liststore.get_value(iter_, 1)
        self.main.bottom_notebook.status_liststore.set_value(iter_, 1, value.replace('Running...', 'Error!'))

        system.e_job_history[system.e_step_counter].update(results)
        
    def _handle_result(self, results):
        """Handle 'RESULT' messages from a worker process.

        Updates system coordinates, visualization objects, and job history.
        """
        e_id  = results['e_id']
        system = self.psystem[e_id]

        if results.get('new_vobject'):
            # Update coordinates
            system.coordinates3 = results['coords']

            # Generate unique name for visualization object
            name = f"{results['simulation_type']}_{system.e_step_counter}"

            # Add visual object to session
            key6 = system.e_job_history[system.e_step_counter]['key6']
            vobject = self._add_vismol_object_to_easyhybrid_session(system=system, name=name, key6 = key6)
            vobject.results = results

            # Save job results in system history
            #system.e_job_history[system.e_step_counter] = results
            system.e_job_history[system.e_step_counter].update(results)
        else:
            # Save job results in system history (no new visual object)
            #system.e_job_history[system.e_step_counter] = results
            system.e_job_history[system.e_step_counter].update(results)
            #self.psystem[e_id].e_step_counter += 1
    
    def _handle_done(self, e_id, process, path, treeiter):
        """Handle 'DONE' messages from a worker process.

        Updates GUI, process pool status, and increments the step counter.
        """
        step_counter = self.psystem[e_id].e_step_counter
        # Update GUI process manager window
        self.main.process_manager_window.set_status(treeiter, "Finished")
        self.main.process_manager_window.set_time(treeiter, False, True)
        self.main.process_manager_window.set_step_counter(treeiter, self.psystem[e_id].e_step_counter)
        self.psystem[e_id].e_job_history[step_counter]['status'] = "Finished"
        # Mark process as completed in pool
        self.process_pool[e_id][2] = "Finished"

        # Update bottom notebook GUI
        iter_ = self.main.bottom_notebook.status_liststore.get_iter(path)
        value = self.main.bottom_notebook.status_liststore.get_value(iter_, 1)
        self.main.bottom_notebook.status_liststore.set_value(iter_, 1, value.replace('Running...', 'Finished!'))

        # Release process resources
        if process.is_alive():
            pass
        else:
            process.join()

        # Increment step counter
        self.psystem[e_id].e_step_counter += 1

        # Mark project/session as changed
        self.changed = True

    # ========================================================================
    # SIMULATION PROCESS
    # ========================================================================

    def _target_process(self, parameters):
        """Target function executed inside a separate process.

        Runs the selected simulation type and communicates results back 
        to the parent process via a multiprocessing.Queue.
        """
        queue = parameters['queue']
        queue.put( self.MSG_RUNNING)  # Notify parent process that job started

        # Map simulation type to class
        simulation_classes = {
            'Energy_Single_Point': EnergyCalculation,
            'Energy_Refinement': EnergyRefinement,
            'Geometry_Optimization': GeometryOptimization,
            'Molecular_Dynamics': MolecularDynamics,
            'Relaxed_Surface_Scan': RelaxedSurfaceScan,
            'Advanced_Relaxed_Surface_Scan':AdvancedRelaxedSurfaceScan,
            'Umbrella_Sampling': UmbrellaSampling,
            'Nudged_Elastic_Band': ChainOfStatesOptimizePath,
            'Normal_Modes': NormalModes,
            }

        sim_type = parameters['simulation_type']
        cls = simulation_classes.get(sim_type)

        if not cls:
            # Unknown simulation type → exit silently
            return

        # Instantiate simulation class
        self.target_process = cls()
        

        
        # Special flags for specific simulation types
        if sim_type == 'Energy_Single_Point':
            parameters['energy'] = True
            parameters['new_vobject'] = False
        elif sim_type == 'Energy_Refinement':
            parameters['new_vobject'] = False

        # Avoid passing queue object to child processes inside parameters
        parameters['queue'] = None
        
        #backup_parameters = parameters
        #backup_parameters = copy.deepcopy(parameters)
        #backup_parameters['system'] = None
        #backup_parameters['step_counter'] = parameters['system'].e_step_counter

        #self.target_process.run(parameters)
        
        try:
            self.target_process.run(parameters)
            results = {
                'new_vobject': parameters.get('new_vobject', True),
                'energy': parameters.get('energy', False),
                'coords': parameters['system'].coordinates3,
                'e_id': parameters['system'].e_id,
                'simulation_type': sim_type,
                'logfile': parameters['logfile'],
                'error': None,
                #'backup_parameters': backup_parameters,
                'step_counter':parameters['system'].e_step_counter
            }
            queue.put((self.MSG_RESULT, results))
            queue.put(self.MSG_DONE)
            
        except Exception as exc:
            print(f"Error {sim_type}: {exc}")
            results = {
                'new_vobject': False,
                'energy': False,
                'coords': None,
                'e_id': parameters['system'].e_id,
                'simulation_type': sim_type,
                'logfile': parameters['logfile'],
                'error': exc,
                #'backup_parameters': backup_parameters,
                'step_counter':parameters['system'].e_step_counter
            }
            queue.put((self.MSG_ERROR, results))
        
    def run_simulation(self, parameters):
        """
        Start a new subprocess for a molecular simulation.

        This method:
          1. Ensures restraints are applied to the active system.
          2. Configures the log file path if not provided.
          3. Guarantees that only one process per system (e_id) runs at a time.
          4. Creates a subprocess and registers it in the process manager.
 
        """
        # --- constants for messages
        self.MSG_RESULT  = "RESULT"
        self.MSG_ERROR   = "ERROR"
        self.MSG_DONE    = "DONE"
        self.MSG_RUNNING = "RUNNING"
        
        # Validate active system
        if self.active_id not in self.psystem:
            self.main.simple_dialog.error(
                msg="No active system selected. Please load or select a system first."
            )
            return False
        
        
        # Assign system and apply restraints
        system = self.psystem[self.active_id]
        key6 = self.vm_session.gen_random_tag_string(length=6)
        #---------------------------------------------------------------
        # Backup parameters
        backup_parameters = copy.deepcopy(parameters)
        backup_parameters['system'] = None
        #---------------------------------------------------------------
        
        system.e_job_history[system.e_step_counter] = {
                'backup_parameters': backup_parameters,
                'new_vobject'      : None,
                'energy'           : None,
                'e_id'             : system.e_id,
                'simulation_type'  : parameters['simulation_type'],
                'error'            : None,
                'step_counter'     : system.e_step_counter,
                'started'          : None,
                'ended'            : None,
                'status'           : 'Queued',
                'potential'        :'UNK',
                'key6'             : key6
            }
        
        
        parameters['system'] = system
        
        # Apply restraints if needed
        parameters = self._apply_restraints(parameters)

        # Ensure process manager window is visible
        self.main.process_manager_window.open_window()
        
        
        parameters = self._configure_logfile(parameters)     
        pprint(parameters)

        e_id, name = system.e_id, system.label

        # Prevent multiple processes for the same system
        if e_id in self.process_pool and self.process_pool[e_id][1].is_alive():
            self.main.simple_dialog.error(
                msg=f"There is already a process underway for the system:\n {e_id} - {name}"
            )
            return False

        # Create queue for inter-process communication
        queue = multiprocessing.Queue()
        parameters.update({
                            'new_vobject': True,
                            'queue': queue
                            })

        
        # Create and start subprocess
        process = multiprocessing.Process(
            target=self._target_process,
            args=(parameters,)
        )
        
        
        hamiltonian = self.get_hamiltonian_type(e_id)
        system.e_job_history[system.e_step_counter]['potential'] = hamiltonian
        #hamiltonian = getattr(system.qcModel, 'hamiltonian', 'unk')
        status = 'Queued'
        
        
        process.start()
        message = f"{parameters['simulation_type']} {system.e_step_counter} - Running..."
        
        #try:
        #    process.start()
        #    # Add entry to process manager window
        #    message = f"{parameters['simulation_type']} {system.e_step_counter} - Running..."
        #    hamiltonian = getattr(system.qcModel, 'hamiltonian', 'unk')
        #    status = 'Queued'
        #except:
        #    message = f"{parameters['simulation_type']} {system.e_step_counter} - Aborted!"
        #    hamiltonian = getattr(system.qcModel, 'hamiltonian', 'unk')
        #    status = 'Aborted!'
        
        treeiter = self.main.process_manager_window.add_new_process(
            system=system,
            _type=parameters['simulation_type'],
            potential=hamiltonian,
            status = status
            )

        # Add entry to bottom notebook
        path = self.main.bottom_notebook.status_teeview_add_new_item(
            message=message,
            system=system
        )

        # Store process information in process pool
        self.process_pool[e_id] = [queue, process, None, path, treeiter]


        if self.GLib_monitor:
            pass
        else:
            # Schedule periodic queue checks
            self.GLib_monitor= GLib.timeout_add(200, self._check_queue)
        return False
