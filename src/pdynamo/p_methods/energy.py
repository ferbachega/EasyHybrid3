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
import pickle
import os

from pScientific.RandomNumbers import NormalDeviateGenerator, RandomNumberGenerator

# --- imports entre modulos adicionados na refatoracao ---
from pdynamo.p_methods._common import backup_orca_files, write_header


# =====================================================================================
#   Helpers  -  used to be duplicated ~6x across the file (1D/2D x pklfolder/vobject).
#   Refactoring them out here means a fix/change only needs to happen in one place,
#   and it also fixes a bug that existed in two of the six copies (see NOTE below).
# =====================================================================================

def compute_reaction_coordinate(coordinates3, rc):
    """ Computes the value(s) of a reaction coordinate for the current frame.

    Returns a tuple (d1, d2):
        - simple_distance            -> (distance, None)
        - multiple_distance          -> (dist(atom1,atom2), dist(atom2,atom3))
        - multiple_distance*4atoms   -> (dist(atom1,atom2), dist(atom3,atom4))

    NOTE: the original code had two copies (2D pklfolder / 2D vobject branches)
    with a typo: `dist1 = didist_RC1_1 - dist_RC1_2` (undefined name
    `didist_RC1_1`), which would raise NameError as soon as anyone used
    'multiple_distance' with a 2D scan. Centralising the logic here removes
    that broken code path entirely.
    """
    rc_type = rc['rc_type']
    atoms   = rc['ATOMS']

    if rc_type == 'simple_distance':
        return coordinates3.Distance(atoms[0], atoms[1]), None

    elif rc_type == 'multiple_distance':
        d1 = coordinates3.Distance(atoms[0], atoms[1])
        d2 = coordinates3.Distance(atoms[1], atoms[2])
        return d1, d2

    elif rc_type == 'multiple_distance*4atoms':
        d1 = coordinates3.Distance(atoms[0], atoms[1])
        d2 = coordinates3.Distance(atoms[2], atoms[3])
        return d1, d2

    return None, None


def _apply_vobject_frame(system, frame):
    """ Copies an in-memory (vobject) frame's xyz array into system.coordinates3. """
    coordinates3 = system.coordinates3
    for idx, xyz in enumerate(frame):
        coordinates3[idx][0] = xyz[0]
        coordinates3[idx][1] = xyz[1]
        coordinates3[idx][2] = xyz[2]


# =====================================================================================
#   Multiprocessing workers
#   ------------------------------------------------------------------------------
#   ASSUMPTION (please validate on your machine): pDynamo QC/MM system objects can
#   be pickled. This seems reasonable given that the existing code already reads/
#   writes '.pkl' trajectory folders, but a full System object (with QC/MM state)
#   is heavier than a coordinates3 object, so this is worth testing on a small
#   case first. If pickling the system fails, `run()` automatically falls back
#   to the original sequential behaviour and prints a warning - it will not
#   silently produce wrong results.
#
#   Each worker process clones the system ONCE (in the pool initializer) and
#   reuses that clone for every frame it processes, so we only pay the pickling
#   cost N_workers times, not once per frame.
# =====================================================================================

_worker_system    = None
_worker_data_path = None


def _pool_initializer(system_pickle, data_path):
    global _worker_system, _worker_data_path
    _worker_system    = pickle.loads(system_pickle)
    _worker_data_path = data_path


def _pool_task_1d(args):
    frame_id, frame, rc1, from_file, full_path_trajectory = args
    system = _worker_system

    if from_file:
        system.coordinates3 = ImportCoordinates3(os.path.join(_worker_data_path, frame))
    else:
        _apply_vobject_frame(system, frame)

    energy  = system.Energy()
    d1, d2  = compute_reaction_coordinate(system.coordinates3, rc1)

    if from_file:
        backup_orca_files(system        = system,
                           output_folder = full_path_trajectory,
                           output_name   = 'frame' + str(frame_id))

    return frame_id, energy, d1, d2


def _pool_task_2d(args):
    key, frame, rc1, rc2, from_file = args
    system = _worker_system

    if from_file:
        system.coordinates3 = ImportCoordinates3(os.path.join(_worker_data_path, frame))
    else:
        _apply_vobject_frame(system, frame)

    energy   = system.Energy()
    d1, _rc1 = compute_reaction_coordinate(system.coordinates3, rc1)
    d2, _rc2 = compute_reaction_coordinate(system.coordinates3, rc2)

    return key, d1, d2, energy


class LogFile:
    """ Class doc """

    def __init__ (self, system):
        """ Class initialiser """
        self.path     = os.path.join(os.environ.get('PDYNAMO3_SCRATCH'), 'summary_temp.log')
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
        energy = parameters['system'].Energy(log = self.logFile2)

        backup_orca_files(system        = parameters['system'],
                          output_folder = parameters['folder'],
                          output_name   = parameters['filename'])

        self.logFile2.Footer ( )
        self.logFile2.Close()
        self.logFile2 = None

        return energy, 'Energy: '+str(energy)


class EnergyRefinement:

    def __init__ (self):
        """ Class initialiser """
        pass

    # ---------------------------------------------------------------------------
    #  Main entry point
    # ---------------------------------------------------------------------------
    def run (self, parameters, interface = False):

        full_path_trajectory = os.path.join(parameters['folder'], parameters['filename'])
        os.mkdir(full_path_trajectory)

        # - - - - - - - - - - - - - Checking trajectory - - - - - - - - - - - - - - - - - -
        self.logFile2 = TextLogFileWriter.WithOptions ( path = os.path.join(full_path_trajectory, 'output.log') )
        parameters['system'].Summary(log = self.logFile2)
        self.logFile2.Header ( )
        self.logFile2.Footer ( )
        self.logFile2.Close()
        self.logFile2 = None

        logfile = self.write_header (parameters = parameters,
                                     logfile    = os.path.join(full_path_trajectory, 'output.log') )

        # ---------------------------------------------------------------
        # Clone the system, keeping GUI-only tree/list-store iterators out
        # of the clone (Clone() has no reason to know about GTK objects).
        # ---------------------------------------------------------------
        backup = []
        try:
            backup.append(parameters['system'].e_treeview_iter)
            backup.append(parameters['system'].e_liststore_iter)
            parameters['system'].e_treeview_iter   = None
            parameters['system'].e_liststore_iter  = None
        except AttributeError:
            pass

        sys_clone = Clone(parameters['system'])

        try:
            parameters['system'].e_treeview_iter   = backup[0]
            parameters['system'].e_liststore_iter  = backup[1]
        except IndexError:
            pass
        # ---------------------------------------------------------------

        parameters['system'] = sys_clone

        # -------------------------------------------------------------------------------
        if parameters['ignore_mm_charges']:
            print('Adjusting electrical charges in the MM region to zero.')
            for i in range(len(parameters['system'].mmState.charges)):
                parameters['system'].mmState.charges[i] = 0.0
        # -------------------------------------------------------------------------------

        n_workers = max(1, int(parameters.get('NmaxThreads', 1) or 1))

        if parameters['is_2D_xy']:
            results = self._run_2d(parameters, full_path_trajectory, n_workers)
            self._write_2d(results, logfile)
        else:
            results = self._run_1d(parameters, full_path_trajectory, n_workers)
            self._write_1d(results, parameters, logfile)

        logfile.close()

    # ---------------------------------------------------------------------------
    #  1D scans  (covers old 'pklfolder' branch + old 'vobject' / not is_2D_xy branch)
    # ---------------------------------------------------------------------------
    def _run_1d (self, parameters, full_path_trajectory, n_workers):
        system    = parameters['system']
        rc1       = parameters['RC1']
        from_file = parameters['traj_type'] in ('pklfolder', 'pklfolder2D')

        tasks = []
        if from_file:
            for frame_name in parameters['trajectory']:
                frame_id = int(frame_name[5:-4])
                tasks.append((frame_id, frame_name))
        else:
            for frame_id, frame in enumerate(parameters['trajectory']):
                tasks.append((frame_id, frame))

        if n_workers > 1:
            try:
                return self._run_parallel_1d(system, tasks, rc1, parameters,
                                              full_path_trajectory, from_file, n_workers)
            except Exception as exc:
                print('[EnergyRefinement] Parallel execution failed ({}); '
                      'falling back to sequential.'.format(exc))

        # --- sequential path (also used as fallback) ---
        results = {}
        for frame_id, frame in tasks:
            if from_file:
                system.coordinates3 = ImportCoordinates3(os.path.join(parameters['data_path'], frame))
            else:
                _apply_vobject_frame(system, frame)

            energy = system.Energy()
            d1, d2 = compute_reaction_coordinate(system.coordinates3, rc1)
            results[frame_id] = (energy, d1, d2)

            print('frame:', frame_id, d1, energy)

            if from_file:
                backup_orca_files(system        = system,
                                  output_folder = full_path_trajectory,
                                  output_name   = 'frame' + str(frame_id))

        return results

    def _run_parallel_1d (self, system, tasks, rc1, parameters, full_path_trajectory, from_file, n_workers):
        system_pickle = pickle.dumps(system)
        data_path      = parameters.get('data_path')

        pool_args = [(frame_id, frame, rc1, from_file, full_path_trajectory)
                     for frame_id, frame in tasks]

        results = {}
        with multiprocessing.Pool(processes  = n_workers,
                                   initializer = _pool_initializer,
                                   initargs    = (system_pickle, data_path)) as pool:
            for frame_id, energy, d1, d2 in pool.imap_unordered(_pool_task_1d, pool_args):
                results[frame_id] = (energy, d1, d2)
                print('frame:', frame_id, d1, energy)

        return results

    def _write_1d (self, results, parameters, logfile):
        rc_type = parameters['RC1']['rc_type']
        lines = []
        for frame_id in sorted(results):
            energy, d1, d2 = results[frame_id]
            if rc_type == 'simple_distance':
                lines.append("\nDATA %9i       %13.12f        %13.12f"
                             % (frame_id, float(d1), float(energy)))
            elif rc_type in ('multiple_distance', 'multiple_distance*4atoms'):
                lines.append("\nDATA %9i       %13.12f        %13.12f        %13.12f"
                             % (frame_id, float(d1), float(d2), float(energy)))

        logfile.write(''.join(lines))

    # ---------------------------------------------------------------------------
    #  2D scans  (covers old 'pklfolder2D' branch + old 'vobject' / is_2D_xy branch)
    # ---------------------------------------------------------------------------
    def _run_2d (self, parameters, full_path_trajectory, n_workers):
        system    = parameters['system']
        rc1       = parameters['RC1']
        rc2       = parameters['RC2']
        from_file = parameters['traj_type'] == 'pklfolder2D'

        tasks = []
        if from_file:
            for frame_name in parameters['trajectory']:
                i_str, j_str = frame_name[5:-4].split('_')
                tasks.append(((int(i_str), int(j_str)), frame_name))
        else:
            for i_j, frame_idx in parameters['idx_2D_xy'].items():
                key = (int(i_j[0]), int(i_j[1]))
                tasks.append((key, parameters['trajectory'][frame_idx]))

        if n_workers > 1:
            try:
                return self._run_parallel_2d(system, tasks, rc1, rc2, parameters,
                                              full_path_trajectory, from_file, n_workers)
            except Exception as exc:
                print('[EnergyRefinement] Parallel execution failed ({}); '
                      'falling back to sequential.'.format(exc))

        # --- sequential path (also used as fallback) ---
        results = {}
        for (i, j), frame in tasks:
            if from_file:
                system.coordinates3 = ImportCoordinates3(os.path.join(parameters['data_path'], frame))
            else:
                _apply_vobject_frame(system, frame)

            energy = system.Energy()
            d1, _  = compute_reaction_coordinate(system.coordinates3, rc1)
            d2, _  = compute_reaction_coordinate(system.coordinates3, rc2)
            results[(i, j)] = (d1, d2, energy)

            print('frame:', (i, j), energy)

        return results

    def _run_parallel_2d (self, system, tasks, rc1, rc2, parameters, full_path_trajectory, from_file, n_workers):
        system_pickle = pickle.dumps(system)
        data_path      = parameters.get('data_path')

        pool_args = [(key, frame, rc1, rc2, from_file) for key, frame in tasks]

        results = {}
        with multiprocessing.Pool(processes  = n_workers,
                                   initializer = _pool_initializer,
                                   initargs    = (system_pickle, data_path)) as pool:
            for key, d1, d2, energy in pool.imap_unordered(_pool_task_2d, pool_args):
                results[key] = (d1, d2, energy)
                print('frame:', key, energy)

        return results

    def _write_2d (self, results, logfile):
        if not results:
            return

        max_i = max(key[0] for key in results)
        max_j = max(key[1] for key in results)

        lines = []
        for i in range(max_i + 1):
            for j in range(max_j + 1):
                if (i, j) in results:
                    d1, d2, energy = results[(i, j)]
                    lines.append("\nDATA  %4i  %4i     %13.12f       %13.12f       %13.12f"
                                 % (i, j, float(d1), float(d2), float(energy)))

        logfile.write(''.join(lines))

    # ---------------------------------------------------------------------------
    #  Header
    # ---------------------------------------------------------------------------
    def write_header (self, parameters, logfile = 'output.log'):
        """ Function doc """

        arq  = open(logfile, "a")
        text = ""

        if parameters['RC2'] is not None :
            text += "\n"
            text += "\n--------------------------------------------------------------------------------"
            text += "\nTYPE                 EasyHybrid Energy Refinement 2D                            "
            text += "\n--------------------------------------------------------------------------------"
        else:
            text += "\n"
            text += "\n--------------------------------------------------------------------------------"
            text += "\nTYPE                   EasyHybrid Energy Refinement                             "
            text += "\n--------------------------------------------------------------------------------"

        # ---- Coordinate 1 -------------------------------------------------------------
        if parameters['RC1']["rc_type"] == 'simple_distance':
            text += "\n"
            text += "\n----------------------- Coordinate 1 - Simple-Distance -------------------------"
            text += "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC1']['ATOMS'][0], parameters['RC1']['ATOM_NAMES'][0] )
            text += "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['RC1']['ATOMS'][1], parameters['RC1']['ATOM_NAMES'][1] )
            text += "\n--------------------------------------------------------------------------------"

        elif parameters['RC1']["rc_type"] == 'multiple_distance':
            text += "\n"
            text += "\n---------------------- Coordinate 1 - multiple-Distance ------------------------"
            text += "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC1']['ATOMS'][0]    , parameters['RC1']['ATOM_NAMES'][0] )
            text += "\nATOM2*                 =%15i  ATOM NAME2             =%15s"     % (parameters['RC1']['ATOMS'][1]    , parameters['RC1']['ATOM_NAMES'][1] )
            text += "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC1']['ATOMS'][2]    , parameters['RC1']['ATOM_NAMES'][2] )
            text += "\n--------------------------------------------------------------------------------"

        elif parameters['RC1']["rc_type"] == 'multiple_distance*4atoms':
            text += "\n"
            text += "\n---------------------- Coordinate 1 - multiple-Distance ------------------------"
            text += "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC1']['ATOMS'][0]    , parameters['RC1']['ATOM_NAMES'][0] )
            text += "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['RC1']['ATOMS'][1]    , parameters['RC1']['ATOM_NAMES'][1] )
            text += "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC1']['ATOMS'][2]    , parameters['RC1']['ATOM_NAMES'][2] )
            text += "\nATOM4                  =%15i  ATOM NAME4             =%15s"     % (parameters['RC1']['ATOMS'][3]    , parameters['RC1']['ATOM_NAMES'][3] )
            text += "\n--------------------------------------------------------------------------------"

        # ---- Coordinate 2 ---------------------------------------------------------------
        if parameters['RC2'] is not None :
            if parameters['RC2']["rc_type"] == 'simple_distance':
                text += "\n"
                text += "\n----------------------- Coordinate 2 - Simple-Distance -------------------------"
                text += "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC2']['ATOMS'][0], parameters['RC2']['ATOM_NAMES'][0] )
                text += "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['RC2']['ATOMS'][1], parameters['RC2']['ATOM_NAMES'][1] )
                text += "\n--------------------------------------------------------------------------------"

            elif parameters['RC2']["rc_type"] == 'multiple_distance':
                text += "\n"
                text += "\n---------------------- Coordinate 2 - multiple-Distance ------------------------"
                text += "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC2']['ATOMS'][0]    , parameters['RC2']['ATOM_NAMES'][0] )
                text += "\nATOM2*                 =%15i  ATOM NAME2             =%15s"     % (parameters['RC2']['ATOMS'][1]    , parameters['RC2']['ATOM_NAMES'][1] )
                text += "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC2']['ATOMS'][2]    , parameters['RC2']['ATOM_NAMES'][2] )
                text += "\n--------------------------------------------------------------------------------"

            elif parameters['RC2']["rc_type"] == 'multiple_distance*4atoms':
                text += "\n"
                text += "\n---------------------- Coordinate 1 - multiple-Distance ------------------------"
                text += "\nATOM1                  =%15i  ATOM NAME1             =%15s"     % (parameters['RC2']['ATOMS'][0]    , parameters['RC2']['ATOM_NAMES'][0] )
                text += "\nATOM2                  =%15i  ATOM NAME2             =%15s"     % (parameters['RC2']['ATOMS'][1]    , parameters['RC2']['ATOM_NAMES'][1] )
                text += "\nATOM3                  =%15i  ATOM NAME3             =%15s"     % (parameters['RC2']['ATOMS'][2]    , parameters['RC2']['ATOM_NAMES'][2] )
                text += "\nATOM4                  =%15i  ATOM NAME4             =%15s"     % (parameters['RC2']['ATOMS'][3]    , parameters['RC2']['ATOM_NAMES'][3] )
                text += "\n--------------------------------------------------------------------------------"

        # ---- Data table header -----------------------------------------------------------
        if parameters['RC2'] is not None :
            text += "\n\n--------------------------------------------------------------------------------"
            text += "\n   Frame i  /  j        RCOORD-1             RCOORD-2                Energy     "
            text += "\n--------------------------------------------------------------------------------"
        else:
            if parameters['RC1']["rc_type"] == 'simple_distance':
                text += "\n\n-------------------------------------------------------------"
                text += "\n           Frame    dist-ATOM1-ATOM2             Energy      "
                text += "\n-------------------------------------------------------------"

            elif parameters['RC1']["rc_type"] == 'multiple_distance':
                text += "\n\n--------------------------------------------------------------------------------"
                text += "\n           Frame     dist-ATOM1-ATOM2      dist-ATOM2-ATOM3         Energy        "
                text += "\n--------------------------------------------------------------------------------  "

            elif parameters['RC1']["rc_type"] == 'multiple_distance*4atoms':
                text += "\n\n--------------------------------------------------------------------------------"
                text += "\n           Frame     dist-ATOM1-ATOM2      dist-ATOM3-ATOM4         Energy        "
                text += "\n--------------------------------------------------------------------------------  "

        arq.write(text)
        return arq
