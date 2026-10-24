#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: keeps the original outputs of the external programs (xtb, ORCA) of a
#          parametrization (pdynamo/param_report.py)
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  The calculations run in scratch folders that are removed afterwards; every
#  runner that takes a `log_dir` copies its inputs/outputs there first.
#  log_dir None: nothing is kept (the default everywhere).
#

import os
import shutil
import threading

_LOCK = threading.Lock ( )


def keep ( work, log_dir, files = ( ), texts = None ):
    """ Copies `files` (names inside the scratch folder `work`, missing ones skipped)
        and writes `texts` ( { name: text } ) into log_dir. """
    if not log_dir: return
    os.makedirs ( log_dir, exist_ok = True )
    for name in files:
        path = os.path.join ( work, name )
        if os.path.isfile ( path ): shutil.copy2 ( path, os.path.join ( log_dir, os.path.basename ( name ) ) )
    for name, text in ( texts or { } ).items ( ):
        with open ( os.path.join ( log_dir, name ), "w" ) as handle: handle.write ( text or "" )


def append ( log_dir, name, header, text ):
    """ Appends one run's output (several single points in one file). """
    if not log_dir: return
    os.makedirs ( log_dir, exist_ok = True )
    with _LOCK, open ( os.path.join ( log_dir, name ), "a" ) as handle:
        handle.write ( "\n{0}\n# {1}\n{0}\n".format ( "#" * 78, header ) )
        handle.write ( text or "" )


def write_xyz ( path, symbols, coordinates, comment = "" ):
    os.makedirs ( os.path.dirname ( path ), exist_ok = True )
    with open ( path, "w" ) as handle:
        handle.write ( "{}\n{}\n".format ( len ( symbols ), comment ) )
        for s, p in zip ( symbols, coordinates ): handle.write ( "{:2s} {:12.6f} {:12.6f} {:12.6f}\n".format ( s, *p ) )
