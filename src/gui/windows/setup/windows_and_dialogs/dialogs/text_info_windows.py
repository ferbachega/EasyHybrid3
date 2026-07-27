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
from gi.repository import Gtk, Pango
from gi.repository import Gdk
from gi.repository import GdkPixbuf
import cairo


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

#import Pickle
from gui.widgets.custom_widgets import FolderChooserButton
from gui.widgets.custom_widgets import VismolTrajectoryFrame
from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.widgets.custom_widgets import get_colorful_square_pixel_buffer
from gui.widgets.custom_widgets import ReactionCoordinateBox

#from gui.widgets.custom_widgets import get_distance
from util.geometric_analysis            import get_distance 
from util.geometric_analysis            import get_dihedral 
from util.geometric_analysis            import get_angle 

from pdynamo.p_methods import LogFile


import util.orca_qc_keywords as orca_keys

from util.file_parser import get_file_type  
from util.file_parser import read_MOL2  
import pprint
import numpy as np
import gc
import os

import traceback


VISMOL_HOME      = os.environ.get('VISMOL_HOME')
HOME             = os.environ.get('HOME')
PDYNAMO3_SCRATCH = os.environ.get('PDYNAMO3_SCRATCH')


import re


# ─────────────────────────────────────────────────────────────────────────────
#  L O G   S Y N T A X   H I G H L I G H T I N G
#
#  Lightweight, regex-based highlighting for the plain-text pDynamo logs
#  shown in TextWindow (system summaries, energy/SCF results, geometry
#  scans, MD/process-manager logs, RMSD/PES analysis output, ...).
#  Colors are chosen to read well on TextView's default light background.
# ─────────────────────────────────────────────────────────────────────────────
_LOG_SEPARATOR_RE      = re.compile(r'^[\-=]{10,}\s*$')
_LOG_EMBEDDED_TITLE_RE = re.compile(r'^-{3,}[^\-]+-{3,}$')
_LOG_DATA_ROW_RE       = re.compile(r'^(DATA)\b')

_LOG_INLINE_RULES = (
    # . Problems -- checked before "success" so "Not Converged" wins over
    #   the "Converged" substring it contains.
    ( re.compile ( r'\b(Error|ERROR|Failed|FAILED|Not\s+Converged|WARNING|Warning)\b' ), 'log_error'   ) ,
    # . Success / good news.
    ( re.compile ( r'(?<!Not\s)(?<!Not)\b(Converged|Successful|Done|OK!)\b' ),            'log_success' ) ,
    # . Timing metadata.
    ( re.compile ( r'\b(Start Time|Stop Time|CPU Time)\b' ),                              'log_time'    ) ,
)


def _ensure_log_tags ( buffer ):
    """Creates (once per buffer) the Gtk.TextTags used for log highlighting."""
    table = buffer.get_tag_table ( )

    def ensure ( name, **props ):
        tag = table.lookup ( name )
        if tag is None:
            tag = buffer.create_tag ( name, **props )
        return tag

    ensure ( 'log_header'  , foreground = '#1B4F72', weight = Pango.Weight.BOLD )
    ensure ( 'log_data'    , foreground = '#1F618D', weight = Pango.Weight.BOLD )
    ensure ( 'log_time'    , foreground = '#7D3C98', style  = Pango.Style.ITALIC )
    ensure ( 'log_success' , foreground = '#1E8449', weight = Pango.Weight.BOLD )
    ensure ( 'log_error'   , foreground = '#B03A2E', weight = Pango.Weight.BOLD )


def apply_log_highlighting ( buffer ):
    """
    Applies syntax highlighting to the whole contents of a Gtk.TextBuffer
    holding a pDynamo-style plain-text log: section separators/titles,
    "DATA" scan/table rows, timestamps, and success/error keywords.
    Safe to call on any plain text -- lines/keywords that don't match
    anything are simply left as-is.
    """
    _ensure_log_tags ( buffer )

    start, end = buffer.get_bounds ( )
    text = buffer.get_text ( start, end, True )

    offset = 0
    for line in text.split ( '\n' ):
        stripped = line.strip ( )

        if _LOG_SEPARATOR_RE.match ( line ) or _LOG_EMBEDDED_TITLE_RE.match ( stripped ):
            s = buffer.get_iter_at_offset ( offset )
            e = buffer.get_iter_at_offset ( offset + len ( line ) )
            buffer.apply_tag_by_name ( 'log_header', s, e )
        else:
            data_match = _LOG_DATA_ROW_RE.match ( stripped )
            if data_match:
                lead = len ( line ) - len ( line.lstrip ( ) )
                s = buffer.get_iter_at_offset ( offset + lead )
                e = buffer.get_iter_at_offset ( offset + lead + data_match.end ( ) )
                buffer.apply_tag_by_name ( 'log_data', s, e )

            for pattern, tag_name in _LOG_INLINE_RULES:
                for m in pattern.finditer ( line ):
                    s = buffer.get_iter_at_offset ( offset + m.start ( ) )
                    e = buffer.get_iter_at_offset ( offset + m.end ( ) )
                    buffer.apply_tag_by_name ( tag_name, s, e )

        offset += len ( line ) + 1   # +1 for the '\n' split away


class TextWindow:
    """ Class doc """
    
    def __init__ (self, text = 'No text', title = None):
        """ Class initialiser """
        self.window = Gtk.Window(title=title)
        self.window.set_default_size(1100, 600)
        
        #if title:
        #    self.window.title = title
        
        
        self.textview = Gtk.TextView()
        self.textbuffer = self.textview.get_buffer()
        self.textbuffer.set_text(text)

        apply_log_highlighting ( self.textbuffer )

        # Create a Pango font description with the desired font family and size
        fontdesc = Pango.FontDescription()
        fontdesc.set_family("Monospace")
        fontdesc.set_size(12 * Pango.SCALE)  # 12 point size
        
        # Apply the font description to the text view
        self.textview.modify_font(fontdesc)
        
        # Set the text color to black
        style = self.textview.get_style_context()
        style.add_class("text-black")
        
        
        scrolledwindow = Gtk.ScrolledWindow()
        scrolledwindow.set_hexpand(True)
        scrolledwindow.set_vexpand(True)
        scrolledwindow.add(self.textview)
        
        self.window.add(scrolledwindow)
        self.window.show_all()


class InfoWindow:
    """ 
    Create a text window. Currently used to display system 
    information or log output. 
    """
    
    def __init__ (self, system = None, text = None):
        """ Class initialiser """
        
        if text:
            pass
        else:
            text = ''
        
        
        if system:
            log  = LogFile(system)
            path = log.path
            with open(path, "r") as f:
                header = f.read()
        else:
            header = '' 
        
        text = header+text
        textwindow = TextWindow(text)
