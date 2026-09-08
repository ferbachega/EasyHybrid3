#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Conjugate Peak Refinement (CPR) log parser
#
#  Copyright 2022-2026 Fernando Bachega
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
#      Parses the "output.log" written by ConjugatePeakRefinementRun (see
#      pdynamo/p_methods/conjugate_peak_refinement.py) into structured
#      data -- pulls out every "Path Summary before/after CPRrun N" table
#      (the reaction-path energy profile as CPR refines it), the exit
#      message, CPU time, and simple run counters (number of CPRrun
#      cycles, "futile loop" retries, unrefinable points in the final
#      path). Meant to feed a future plotting window the same way
#      util/md_analysis.py feeds RMSF_analysis_window.py/RDF_analysis_
#      window.py -- X = image index, Y = energy, for the final path
#      profile; X = CPRrun, Y = saddle energy, for a convergence plot.
#
import re


_ROW_RE = re.compile(
    r'^\s*(?P<image>\d+)\s+'
    r'(?P<energy>[+-]?\d+\.\d+)\s+'
    r'(?P<rel_energy>[+-]?\d+\.\d+)\s+'
    r'(?P<type>\S+)\s+'
    r'(?P<rms_gradient>n/a|[+-]?\d+\.\d+)\s*$'
)
_SEPARATOR_RE = re.compile(r'^-{10,}\s*$')
_TITLE_RE     = re.compile(r'Path Summary (?P<when>before|after) CPRrun (?P<run>\d+)')
_SYSTEM_RE    = re.compile(r'Summary of System "(?P<name>[^"]*)"')
_CPU_TIME_RE  = re.compile(r'CPU Time:\s*([\d.]+)s')


def _parse_table_after(lines, title_index):
    """ Given the line index of a "Path Summary ... CPRrun N" title, find
        and parse the data rows of the table that follows it. Layout is
        fixed (see CPRRefinement.PathSummary() in addOns/pyCPR/
        ConjugatePeakRefinement.py): a separator line, the title, another
        separator, the column headings, another separator, the data
        rows, then a closing separator.
    """
    n = len(lines)
    i = title_index + 1
    # . Skip to the heading line ("  Image  Energy  ...").
    while i < n and not lines[i].strip().startswith('Image'):
        if i - title_index > 5:  # . Defensive: heading should be within a few lines.
            return None
        i += 1
    if i >= n:
        return None
    i += 1  # . Skip the heading line itself.
    # . Skip the separator line under the heading.
    while i < n and not _SEPARATOR_RE.match(lines[i]):
        i += 1
    i += 1
    table = {'image': [], 'energy': [], 'rel_energy': [], 'type': [], 'rms_gradient': []}
    while i < n and not _SEPARATOR_RE.match(lines[i]):
        match = _ROW_RE.match(lines[i])
        if match:
            table['image'].append(int(match.group('image')))
            table['energy'].append(float(match.group('energy')))
            table['rel_energy'].append(float(match.group('rel_energy')))
            table['type'].append(match.group('type'))
            rmsGrad = match.group('rms_gradient')
            table['rms_gradient'].append(None if rmsGrad == 'n/a' else float(rmsGrad))
        i += 1
    if len(table['image']) == 0:
        return None
    return table


def parse_cpr_log(log_path):
    """ Parses a Conjugate Peak Refinement output.log file.

        Returns a dict:
            log_file             : the path passed in
            system_name          : str or None
            converged             : bool (True if the exit message says "finished successfully")
            exit_message          : str or None
            cpu_time_seconds      : float or None
            number_of_cprruns     : int (highest CPRrun index seen)
            number_of_futile_loops: int
            initial_path          : table dict (the "before CPRrun 0" path), or None
            final_path            : table dict (the last "after CPRrun N" path), or None
            path_history          : list of {'run': int, **table dict}, one per "after CPRrun N" block, in order
            saddle_images         : list of (run, image_index, energy, rel_energy, rms_gradient) from final_path
    """
    with open(log_path, 'r') as handle:
        lines = handle.readlines()

    result = {
        'log_file'              : log_path,
        'system_name'           : None,
        'converged'             : False,
        'exit_message'          : None,
        'cpu_time_seconds'      : None,
        'number_of_cprruns'     : 0,
        'number_of_futile_loops': 0,
        'initial_path'          : None,
        'final_path'            : None,
        'path_history'          : [],
        'saddle_images'         : [],
    }

    for index, line in enumerate(lines):
        systemMatch = _SYSTEM_RE.search(line)
        if systemMatch and result['system_name'] is None:
            result['system_name'] = systemMatch.group('name')

        if 'Futile loop detected' in line:
            result['number_of_futile_loops'] += 1

        cpuMatch = _CPU_TIME_RE.search(line)
        if cpuMatch:
            result['cpu_time_seconds'] = float(cpuMatch.group(1))

        if line.strip() == 'CPR terminated':
            for nextLine in lines[index + 1:]:
                if nextLine.strip():
                    result['exit_message'] = nextLine.strip()
                    break

        titleMatch = _TITLE_RE.search(line)
        if titleMatch is None:
            continue

        when = titleMatch.group('when')
        run  = int(titleMatch.group('run'))
        table = _parse_table_after(lines, index)
        if table is None:
            continue

        if when == 'before' and run == 0:
            result['initial_path'] = table
        elif when == 'after':
            result['number_of_cprruns'] = max(result['number_of_cprruns'], run)
            entry = {'run': run}
            entry.update(table)
            result['path_history'].append(entry)
            result['final_path'] = table

    if result['exit_message'] is not None:
        result['converged'] = 'finished successfully' in result['exit_message']

    if result['final_path'] is not None:
        run = result['path_history'][-1]['run'] if result['path_history'] else None
        table = result['final_path']
        for image, energy, relEnergy, imageType, rmsGradient in zip(
                table['image'], table['energy'], table['rel_energy'], table['type'], table['rms_gradient']):
            if imageType == 'saddle':
                result['saddle_images'].append((run, image, energy, relEnergy, rmsGradient))

    return result
