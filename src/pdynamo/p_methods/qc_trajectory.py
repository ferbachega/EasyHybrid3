#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  qc_trajectory.py
#
#  Copyright 2022-2025 Fernando Bachega <ferbachega@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
"""
qc_trajectory
=============

A drop-in replacement for pDynamo's ExportTrajectory that ALSO backs up the
external-QC-program log (ORCA / xTB / ...) of each frame, right next to the
frame's geometry.

Why this exists
---------------
pDynamo's geometry trajectory (SystemGeometryTrajectory, the '.ptGeo' format)
saves one 'frame{N}.pkl' per frame via its WriteOwnerData() method. It knows
nothing about QC program logs. By subclassing it and overriding WriteOwnerData()
we can, immediately after the geometry frame is written, copy that frame's QC log
into the SAME '.ptGeo' folder as 'frame{N}.<engine>.log'. Because self.owner is
the System (which carries the qcModel) and the frame index is known, each log is
correctly paired with its geometry and, in the umbrella-sampling / scan case,
with the reaction-coordinate value that produced it.

This makes analyses such as "Mayer/Wiberg bond orders vs. reaction-coordinate
distance" possible: every window/frame keeps its own QC log, from which the
bond orders (and other properties) can be parsed.

Scope / caveat
--------------
This pairs the log with the frame reliably when the QC calculation of that frame
is the most recent one when the frame is written -- e.g. umbrella sampling or a
scan where each window ends with a geometry optimization. In plain MD, many QC
steps run between two trajectory writes and share one scratch log, so the log
captured for an MD frame is whatever was last written to scratch, not a perfect
per-step match. That is inherent to how the scratch log is reused and is fine for
the intended (optimized-window) use case.

Usage
-----
Replace, in a simulation method:

    from pBabel import ExportTrajectory
    trajectory = ExportTrajectory ( path, system, log = None )

with:

    from pdynamo.p_methods.qc_trajectory import EasyHybridExportTrajectory
    trajectory = EasyHybridExportTrajectory ( path, system, log = None )

The returned object behaves exactly like the pDynamo trajectory, plus the log
backup. If anything about the QC backup fails, it is swallowed so the trajectory
itself never breaks.
"""

import os

from pBabel import SystemGeometryTrajectory

from pdynamo.p_methods._common import backup_qc_files
from util.debug              import dprint


class EasyHybridExportTrajectory ( SystemGeometryTrajectory ):
    """A SystemGeometryTrajectory that also backs up each frame's QC log."""

    def WriteOwnerData ( self, index = -1 ):
        """Write the geometry frame, then back up this frame's QC log.

        The geometry is written exactly as pDynamo does (via super). Then the
        QC-program log currently in scratch is copied into the trajectory folder
        as 'frame{N}.<engine>.log', with N the same index used for the geometry.
        """
        # frame number that the base class is about to use (see WriteFrame:
        # it uses self.numberOfFrames when frame is None, i.e. index < 0)
        frame = index if index >= 0 else self.numberOfFrames

        # 1) write the geometry frame exactly as pDynamo would
        super ( EasyHybridExportTrajectory, self ).WriteOwnerData ( index = index )

        # 2) back up this frame's QC log next to it (never break the trajectory)
        try:
            backup_qc_files ( system        = self.owner                       ,
                              output_folder = self.path                        ,
                              output_name   = "frame{:d}".format ( frame )      )
        except Exception as exc:
            dprint ( "QC log backup for frame {} failed: {}".format ( frame, exc ) )


def EasyHybridExportTrajectory_factory ( path, owner, append = False, log = None ):
    """Drop-in replacement for ExportTrajectory that backs up QC logs per frame.

    Mirrors pDynamo's ExportTrajectory call signature closely enough for the
    geometry-trajectory ('.ptGeo') use in EasyHybrid. 'log' is accepted and
    ignored (EasyHybrid always passes log = None here).
    """
    return EasyHybridExportTrajectory.WriterFromPathAndOwner ( path, owner, append = append )
