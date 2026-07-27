#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Debug/verbose-output toggle
#
#  Description:
#      Centralises all the ad-hoc "print(...)" debug statements scattered
#      across the codebase behind a single ON/OFF switch, so that by
#      default EasyHybrid3 runs quietly, but the exact same messages can be
#      brought back at any time for troubleshooting.
#
#      Enable by setting the environment variable before starting
#      EasyHybrid3, e.g.:
#
#          EASYHYBRID_DEBUG=1 python easyhybrid.py
#
#      or at runtime from within the application:
#
#          from util.debug import set_debug
#          set_debug(True)
#
import os

_TRUE_VALUES = ( '1', 'true', 'yes', 'on' )

# . Read once at import time; can still be overridden at runtime via
#   set_debug() (e.g. from a future "Debug mode" menu item/checkbox).
DEBUG = os.environ.get ( 'EASYHYBRID_DEBUG', '0' ).strip ( ).lower ( ) in _TRUE_VALUES


def set_debug ( value ):
    """ Turns debug output on/off at runtime. """
    global DEBUG
    DEBUG = bool ( value )


def is_debug ( ):
    """ Returns whether debug output is currently enabled. """
    return DEBUG


def dprint ( *args, **kwargs ):
    """
    Drop-in replacement for the builtin print(), silenced unless debug
    output is enabled (see set_debug() / the EASYHYBRID_DEBUG environment
    variable). Used to gate the countless development/debug prints spread
    throughout the codebase without deleting the information they carry.
    """
    if DEBUG:
        print ( *args, **kwargs )
