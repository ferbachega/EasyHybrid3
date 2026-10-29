#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: Baker saddle-point search with two pDynamo3 fixes
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
#      Drop-in replacement for pSimulation's BakerSaddleOptimize_SystemGeometry
#      (2026-10-19). pDynamo3's BakerOptimizer has two problems, fixed here in a
#      subclass so the local pDynamo3 installation is left untouched:
#
#      1. followMode > 0: mode tracking calls Array.AbsoluteMaximumIndex(), which
#         does not exist, so the search stops after 3 steps with an error.
#      2. analyticFrequency > 0: SystemGeometryObjectiveFunction has no
#         FunctionGradientsHessian (the base class stub returns 0 and leaves the
#         gradients and Hessian untouched), so those steps were wasted. Here the
#         Hessian is recomputed by finite differences, the same Hessian pDynamo
#         uses for the first step and for normal modes.
#
import math

from pScientific.ObjectiveFunctionIterators import BakerOptimizer
from pSimulation.GeometryOptimization       import _SetUpOptimization


class BakerOptimizerEH ( BakerOptimizer ):
    """BakerOptimizer with working mode following and Hessian recalculation."""

    def DetermineStep ( self, state ):
        """Same as BakerOptimizer.DetermineStep except for the mode-tracking index."""
        state.eigenVectors.VectorMultiply ( state.g, state.h, transpose = True )
        if self.useNewtonRaphsonStep and ( state.numberOfNegativeModes == self.numberOfNegativeModes ):
            state.h.Divide ( state.eigenValues )
            state.eigenVectors.VectorMultiply ( state.h, state.d )
        else:
            # . Follow the mode with the largest overlap with the previous step.
            if ( self.numberOfNegativeModes > 0 ) and ( state.currentMode > 0 ) and ( state.numberOfIterations > 1 ):
                state.eigenVectors.VectorMultiply ( state.d, state.w, transpose = True )
                state.currentMode = max ( range ( state.numberOfVariables ), key = lambda i: math.fabs ( state.w[i] ) )
            if ( state.rmsGradient <= self.rmsGradientTolerance ) and ( state.numberOfNegativeModes != self.numberOfNegativeModes ):
                v = max ( 0, state.currentMode )
                for i in range ( state.numberOfVariables ):
                    state.d[i] = self.resetStep * state.eigenVectors[i,v]
            else:
                for ( i, ( e, g ) ) in enumerate ( zip ( state.eigenValues, state.h ) ):
                    state.w[i] = - 2.0 * g / ( math.fabs ( e ) + math.sqrt ( e * e + 4.0 * g * g ) )
                if state.currentMode >= 0: state.w[state.currentMode] *= -1.0
                state.eigenVectors.VectorMultiply ( state.w, state.d )
        state.objectiveFunction.ApplyConstraintsToVector ( state.d )
        dSize = state.d.Norm2 ( )
        if   dSize < self.minimumStep: state.d.Scale ( self.minimumStep / dSize )
        elif dSize > self.maximumStep: state.d.Scale ( self.maximumStep / dSize )
        state.stepType = "N"
        state.x.Add ( state.d )

    def FunctionGradients ( self, state ):
        """Recompute the Hessian every analyticFrequency steps, otherwise update it."""
        # . numberOfIterations is incremented after this call; the starting call (0)
        #   already has the starting Hessian.
        if ( self.analyticFrequency > 0 ) and ( state.stepType != "I" ) and \
           ( ( state.numberOfIterations + 1 ) % self.analyticFrequency == 0 ):
            # . StartingHessian() is the finite-difference Hessian with rotations and
            #   translations projected out. It moves the coordinates, so the energy and
            #   gradients are evaluated afterwards.
            state.hessian = state.objectiveFunction.StartingHessian ( state.x )
            state.f       = state.objectiveFunction.FunctionGradients ( state.x, state.g )
            state.numberOfFunctionCalls += 1
            state.rmsGradient = state.g.RootMeanSquare ( )
            self.TreatHessian ( state )
            state.stepType += "H"
        else:
            super ( BakerOptimizerEH, self ).FunctionGradients ( state )


def BakerSaddleOptimizeEH_SystemGeometry ( system, **options ):
    """Saddle-point search; same arguments as BakerSaddleOptimize_SystemGeometry."""
    # . BOFILL by default. pDynamo's own test for this runs after the class defaults
    #   (BFGS) were merged in, so it never applied.
    options = dict ( options ) ; options.setdefault ( "hessianUpdatingOption", "BOFILL" )
    ( of, options, log ) = _SetUpOptimization ( system, BakerOptimizerEH._attributable, {}, options, True )
    options["followMode"           ] = max ( options["followMode"], 0 )
    options["numberOfNegativeModes"] = 1
    optimizer = BakerOptimizerEH.WithOptions ( **options )
    optimizer.Summary ( log = log )
    return optimizer.Iterate ( of, log = log )
