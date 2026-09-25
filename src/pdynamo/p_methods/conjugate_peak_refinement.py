#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Conjugate Peak Refinement (CPR) transition-state search
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
#      Runs pDynamo3's addOns.pyCPR (Conjugate Peak Refinement, Fischer &
#      Karplus 1992 / Gisdon, Culka & Ullmann 2016) transition-state
#      search: given only reactant and product structures (no initial
#      Hessian needed, unlike the Baker search), builds an interpolated
#      initial path with GrowingStringInitialPath() -- the same path
#      builder ChainOfStatesOptimizePath (NEB) already uses in
#      molecular_dynamics's sibling chain_of_states.py -- then refines
#      that path with ConjugatePeakRefinementOptimizePath() until the
#      highest-energy image on it is a genuine saddle point.
#
import os

from pBabel                    import *
from pCore                     import *
from pMolecule                 import *
from pScientific               import *
from pScientific.Arrays        import *
from pScientific.Geometry3     import *
from pSimulation                import *
from addOns.pyCPR               import ConjugatePeakRefinementOptimizePath


class ConjugatePeakRefinementRun:
    """ Class doc """

    def __init__(self):
        """ Class initialiser """
        self.trajectory = None
        self.logFile2   = None

    def run(self, parameters):
        """ Function doc """
        trajectoryPath = os.path.join(parameters['folder'], parameters['trajectory_name'] + ".ptGeo")
        system = parameters['system']

        # . Assign the reactant and product coordinates.
        reactants = parameters['reac_coordinates']
        products  = parameters['prod_coordinates']

        # . Get an initial (interpolated) path -- same builder used by
        # Nudged_Elastic_Band's chain_of_states.py, reused here so CPR
        # gets a reasonable initial guess instead of just the two raw
        # endpoints.
        GrowingStringInitialPath(system, parameters['number_of_structures'], reactants, products, trajectoryPath, log=None)

        logPath = os.path.join(trajectoryPath, 'output.log')
        parameters['logfile'] = logPath

        # [EN] BUG FIX (reported by the user: the pDynamo3 example script
        # for CPR prints far more detail than our saved log). A first
        # attempt reassigned the "logFile" NAME inside each addOns.pyCPR
        # submodule to self.logFile2, which fixed every DIRECT
        # "logFile.Paragraph(...)" call -- but several key methods
        # (CPRSaddlePointRefinement.LogHeader/LogIteration/LogTableStop,
        # ConjugatePeakRefinement.StructureToOptimizer) declare
        # "log=logFile" as a DEFAULT ARGUMENT and are always called
        # with no override anywhere in the call chain. Python evaluates
        # a default argument once, at "def" time (module import), and
        # freezes a reference to whatever OBJECT "logFile" pointed to
        # then -- reassigning the module-level NAME afterwards can never
        # reach an already-bound default (that's exactly the content
        # still missing: the "==== Starting CPRrun N ====" headers and
        # the whole per-iteration "RMS gradient tolerance ... Tau ...
        # RMSGradient" table from the saddle-point line search, i.e. the
        # single most useful diagnostic block for understanding a CPR
        # run). The only fix that reaches BOTH call shapes uniformly is
        # mutating pCore's global "logFile" SINGLETON's own .file/.path
        # in place -- every frozen default and every live lookup already
        # points at this exact same object, so redirecting where IT
        # writes, without replacing the object itself, redirects
        # everything at once. Safe here because each simulation job runs
        # in its own dedicated multiprocessing.Process (see
        # simulations_mixin.py's run_simulation()) -- this global is
        # only ever touched by the one job using it.
        defaultLog = logFile
        originalFile, originalPath, originalActive = defaultLog.file, defaultLog.path, defaultLog.isActive
        defaultLog.file     = open(logPath, "w", 1)
        defaultLog.path     = logPath
        defaultLog.isActive = True
        self.logFile2 = defaultLog

        system.Summary(log=self.logFile2)
        self.logFile2.Header()

        self.trajectory = ExportTrajectory(trajectoryPath, system, append=True)

        try:
            # [EN] CPRIterate() traps its own internal failures (e.g.
            # "Need at least two different structures!") into the
            # returned state dict's "Error" key without raising -- same
            # silent-failure shape already found and fixed for
            # SteepestDescentPath_SystemGeometry in geometry_optimization.
            # py's _run_reaction_path(). Raising here routes it through
            # the normal Geometry_Optimization-style error path instead
            # of masquerading as a successful, empty run.
            state = ConjugatePeakRefinementOptimizePath(
                system, self.trajectory,
                rmsGradientTolerance       = parameters['rmsGradientTolerance'],
                rmsdMaximumTolerance       = parameters['rmsdMaximumTolerance'],
                maxCPRrun                  = parameters['maxCPRrun'],
                stepsPerSegment            = parameters['stepsPerSegment'],
                betaType                   = parameters['betaType'],
                initialOrthogonalRuns      = parameters['initialOrthogonalRuns'],
                increaseTau                = parameters['increaseTau'],
                breakIfTauReached          = parameters['breakIfTauReached'],
                finalUnrefinableRefinement = parameters['finalUnrefinableRefinement'],
                log                        = self.logFile2,
            )
        finally:
            self.logFile2.Footer()  # . Also closes defaultLog.file (path is set).
            defaultLog.file, defaultLog.path, defaultLog.isActive = originalFile, originalPath, originalActive
            self.logFile2 = None

        if state.get('Error') is not None:
            raise ValueError("Conjugate Peak Refinement failed: {}".format(state['Error']))
