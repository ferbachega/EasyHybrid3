#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: QC engine (XTB/ORCA/DFTB+) scratch folder and executable path checks
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
#      For systems whose QC model is XTB, ORCA or DFTB+, checks whether the
#      configured scratch (temporary files) folder, the external program's
#      executable, and (DFTB+ only) the Slater-Koster parameter folder
#      actually exist on THIS machine -- none of this is verified anywhere
#      else in the app before an actual QC calculation is attempted
#      (pDynamo3's own qcModel.command is a lazily-evaluated property, only
#      ever accessed -- and only then raising NotInstalledError -- when a
#      calculation runs). A system pickled on one machine (e.g. inside a
#      .easy file) and reopened on another can easily have any of these
#      paths pointing somewhere that no longer exists.
#
#      Deliberately GTK-free (no gi.repository imports) -- callers that
#      want to show a report to the user build their own widgets from the
#      plain dicts/strings this module returns (see io_data.py's
#      load_easyhybrid_serialization_file for the live example).
#
#      Real crash this was built to catch ahead of time (both the scratch/
#      executable check and, later, the skfPath one): a system's DFTB+
#      qcModel.skfPath named a specific parameter set folder (e.g.
#      ".../examples/dftbPlus/data/3ob-3-1") that didn't exist on the
#      machine reopening it -- only ".../data/skf" did -- and only failed
#      once an Energy calculation actually ran pMolecule.QCModel.
#      QCModelError: "Error executing program.".
#

import os
import glob

# QCModelXTB/QCModelORCA/QCModelDFTB (pMolecule.QCModel) all read their
# (lazily-evaluated, cached) '.command' property from one of these env
# vars if no command was ever explicitly set on the instance -- see each
# class's own 'command' @property in pDynamo3's pMolecule/QCModel/*.py.
# MNDO/MOPAC and any other QC model type have no external executable/
# scratch-folder concept in the same sense and are not covered here.
_QC_ENGINE_INFO = {
    'QCModelXTB' : {'label': 'XTB'  , 'command_env_var': 'PDYNAMO3_XTBCOMMAND' },
    'QCModelORCA': {'label': 'ORCA' , 'command_env_var': 'PDYNAMO3_ORCACOMMAND'},
    'QCModelDFTB': {'label': 'DFTB+', 'command_env_var': 'PDYNAMO3_DFTBCOMMAND'},
}


def check_qc_engine_paths (system):
    """Checks a System's QC model (if it is XTB, ORCA or DFTB+) for: a)
    whether its configured scratch/temp-files folder exists, b) whether
    its executable's path (resolved the same way pDynamo3's own
    qcModel.command property would) exists and is executable, and c) for
    DFTB+ only, whether its Slater-Koster parameter folder (skfPath)
    exists and actually contains at least one ".skf" file.

    Pure/read-only -- never touches the system. See fix_qc_engine_paths()
    for the auto-correcting counterpart.

    Returns None if the system has no QC model, or its QC model isn't one
    of the three engines this covers (nothing to check). Otherwise returns
    a report dict:
        {
            'engine'                : 'XTB' | 'ORCA' | 'DFTB+',
            'qc_model_class'        : the exact QCModel class name,
            'scratch_path'          : str or None,
            'scratch_exists'        : bool,
            'scratch_parent_exists' : bool,  # informational: pDynamo3's own
                                              # DeterminePaths() auto-mkdir's
                                              # a missing scratch dir, but
                                              # via os.mkdir (not makedirs),
                                              # so it still fails if the
                                              # PARENT doesn't exist either.
            'executable_env_var'    : the env var name this engine reads,
            'executable_path'       : str or None (whatever that env var
                                       is currently set to, if anything),
            'executable_exists'     : bool,   # os.path.isfile
            'executable_executable' : bool,   # also os.access(X_OK)
            'skf_path'              : str or None (DFTB+ only, else None),
            'skf_path_exists'       : bool (DFTB+ only, else True -- so it
                                       never drags down 'ok' for XTB/ORCA),
            'skf_path_has_skf_files': bool (DFTB+ only, else True; a
                                       generic "looks like a real SK
                                       folder" proxy -- does NOT check that
                                       every element pair this system
                                       actually needs has its own .skf,
                                       only that the folder isn't empty),
            'ok'                    : bool,
        }
    """
    qc_model = getattr(system, 'qcModel', None)
    if qc_model is None:
        return None

    qc_model_class = type(qc_model).__name__
    info = _QC_ENGINE_INFO.get(qc_model_class)
    if info is None:
        return None

    scratch_path = getattr(qc_model, 'scratch', None)
    scratch_exists = bool(scratch_path) and os.path.isdir(scratch_path)
    scratch_parent_exists = bool(scratch_path) and os.path.isdir(os.path.dirname(scratch_path) or os.sep)

    executable_path = os.environ.get(info['command_env_var'])
    executable_exists = bool(executable_path) and os.path.isfile(executable_path)
    executable_executable = executable_exists and os.access(executable_path, os.X_OK)

    if qc_model_class == 'QCModelDFTB':
        skf_path = getattr(qc_model, 'skfPath', None)
        skf_path_exists = bool(skf_path) and os.path.isdir(skf_path)
        skf_path_has_skf_files = skf_path_exists and len(glob.glob(os.path.join(skf_path, '*.skf'))) > 0
    else:
        skf_path = None
        skf_path_exists = True
        skf_path_has_skf_files = True

    return {
        'engine'                 : info['label'],
        'qc_model_class'         : qc_model_class,
        'scratch_path'           : scratch_path,
        'scratch_exists'         : scratch_exists,
        'scratch_parent_exists'  : scratch_parent_exists,
        'executable_env_var'     : info['command_env_var'],
        'executable_path'        : executable_path,
        'executable_exists'      : executable_exists,
        'executable_executable'  : executable_executable,
        'skf_path'               : skf_path,
        'skf_path_exists'        : skf_path_exists,
        'skf_path_has_skf_files' : skf_path_has_skf_files,
        'ok'                     : scratch_exists and executable_executable and skf_path_has_skf_files,
    }


def _fallback_scratch_path ():
    """The one scratch location every other tool in this app already
    dumps QC job files into successfully ($PDYNAMO3_SCRATCH itself, not
    a per-engine subfolder -- confirmed by inspecting a real, long-lived
    scratch folder: XTB/ORCA job files sit directly in it, not under an
    'XTBScratch'/'orcaScratch' subfolder). Returns None if that env var
    isn't set or doesn't point to a real directory -- no fallback."""
    base = os.environ.get('PDYNAMO3_SCRATCH')
    if base and os.path.isdir(base):
        return base
    return None


def _fallback_skf_path ():
    """DFTB+'s own bundled example Slater-Koster set -- the exact same
    default src/gui/.../qc_setup/setup_dftbplus.py itself pre-fills for a
    NEW DFTB+ setup. Returns None unless it's actually there with real
    .skf files in it (never propose a fallback that would just fail the
    same way)."""
    home = os.environ.get('PDYNAMO3_HOME')
    if not home:
        return None
    candidate = os.path.join(home, 'examples', 'dftbPlus', 'data', 'skf')
    if os.path.isdir(candidate) and glob.glob(os.path.join(candidate, '*.skf')):
        return candidate
    return None


def fix_qc_engine_paths (system):
    """Runs check_qc_engine_paths(system) and, for any problem that has a
    known, REAL, working fallback on THIS machine (scratch folder, DFTB+'s
    skfPath), applies it directly to system.qcModel and records the
    change. Does NOT attempt to fix the executable path -- there's no safe
    way to guess "the right" xtb/orca/dftb+ binary if the configured one
    is missing without risking silently running a DIFFERENT version than
    the one actually intended, so that one stays report-only.

    Returns (final_report, fixes): final_report is check_qc_engine_paths()
    re-run after any fixes (so it reflects the corrected state); fixes is
    a list of {'field': 'scratch'|'skfPath', 'old': <old value>,
    'new': <new value>} for whichever problems were actually corrected
    (empty if the system has no QC model, was already ok, or nothing
    fixable was found).
    """
    report = check_qc_engine_paths(system)
    if report is None or report['ok']:
        return report, []

    qc_model = system.qcModel
    fixes = []

    if not report['scratch_exists']:
        fallback = _fallback_scratch_path()
        if fallback:
            fixes.append({'field': 'scratch', 'old': report['scratch_path'], 'new': fallback})
            qc_model.scratch = fallback

    if report['qc_model_class'] == 'QCModelDFTB' and not report['skf_path_has_skf_files']:
        fallback = _fallback_skf_path()
        if fallback:
            fixes.append({'field': 'skfPath', 'old': report['skf_path'], 'new': fallback})
            qc_model.skfPath = fallback

    if fixes:
        report = check_qc_engine_paths(system)

    return report, fixes


_FIELD_LABELS = {
    'scratch': 'scratch folder',
    'skfPath': 'Slater-Koster parameter folder',
}


def describe_qc_engine_report (system_label, report, fixes = ()):
    """Combines a (post-fix) check_qc_engine_paths() report and the fixes
    fix_qc_engine_paths() applied into a list of entries suitable for
    building a GUI report:
        [{'text': <human-readable line>,
          'highlight': <substring to visually emphasize, e.g. the new
                        redirected path -- or None>}, ...]
    Fixed items are listed first (good news, with the new path meant to
    be highlighted); remaining, un-fixable problems (if any) follow with
    no highlight. Returns [] if report is None or fully ok with no fixes
    applied (nothing worth telling the user about)."""
    if report is None:
        return []

    entries = []
    fixed_fields = {fix['field'] for fix in fixes}

    for fix in fixes:
        entries.append({
            'text': 'System "{}" ({} QC Model): {} did not exist ({}) -- redirected to: {}'.format(
                system_label, report['engine'], _FIELD_LABELS[fix['field']], fix['old'], fix['new']),
            'highlight': fix['new'],
        })

    if report['ok']:
        return entries

    if not report['scratch_exists'] and 'scratch' not in fixed_fields:
        if not report['scratch_path']:
            entries.append({'highlight': None, 'text':
                'System "{}" ({} QC Model): no scratch folder configured.'.format(
                    system_label, report['engine'])})
        elif not report['scratch_parent_exists']:
            entries.append({'highlight': None, 'text':
                'System "{}" ({} QC Model): scratch folder does not exist and cannot be '
                'auto-created (its parent directory is also missing): {}'.format(
                    system_label, report['engine'], report['scratch_path'])})
        else:
            entries.append({'highlight': None, 'text':
                'System "{}" ({} QC Model): scratch folder does not exist yet (will be '
                'created automatically when needed): {}'.format(
                    system_label, report['engine'], report['scratch_path'])})

    if not report['executable_executable']:
        if not report['executable_path']:
            entries.append({'highlight': None, 'text':
                'System "{}" ({} QC Model): {} is not set -- the {} executable path is unknown.'.format(
                    system_label, report['engine'], report['executable_env_var'], report['engine'])})
        elif not report['executable_exists']:
            entries.append({'highlight': None, 'text':
                'System "{}" ({} QC Model): {} points to a file that does not exist: {}'.format(
                    system_label, report['engine'], report['executable_env_var'], report['executable_path'])})
        else:
            entries.append({'highlight': None, 'text':
                'System "{}" ({} QC Model): {} is not executable: {}'.format(
                    system_label, report['engine'], report['executable_env_var'], report['executable_path'])})

    if (report['qc_model_class'] == 'QCModelDFTB' and not report['skf_path_has_skf_files']
            and 'skfPath' not in fixed_fields):
        if not report['skf_path']:
            entries.append({'highlight': None, 'text':
                'System "{}" ({} QC Model): no Slater-Koster parameter folder (skfPath) configured.'.format(
                    system_label, report['engine'])})
        elif not report['skf_path_exists']:
            entries.append({'highlight': None, 'text':
                'System "{}" ({} QC Model): Slater-Koster parameter folder does not exist: {}'.format(
                    system_label, report['engine'], report['skf_path'])})
        else:
            entries.append({'highlight': None, 'text':
                'System "{}" ({} QC Model): Slater-Koster parameter folder has no ".skf" files: {}'.format(
                    system_label, report['engine'], report['skf_path'])})

    return entries
