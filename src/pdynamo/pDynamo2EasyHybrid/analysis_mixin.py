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
from util.debug import dprint
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


class _StatusbarHoverLabel:
    """ [EN] User's own request: "quero colocar uma statusbar que 'on
    the fly' mostra os dados de x e y com a movimentacao do mouse" on
    the WHAM results windows (histograms/PMF for 1D, the PMF heatmap for
    2D). Both easyplot.ImagePlot (already had this -- see its own self.
    RC_label / on_motion()) and easyplot.XYPlot (extended alongside this
    change, same convention) already compute the live x/y[/z] value under
    the cursor on every motion event -- they just had nowhere real to
    show it (ImagePlot's own RC_label defaulted to None everywhere it was
    actually used, silently falling back to dprint(); XYPlot had no such
    hook at all until now). Both call `.set_text(str)` on whatever
    self.RC_label is set to. A real Gtk.Statusbar has no such method
    (it uses push()/pop() with a context id instead of a plain "replace
    the text" call) -- this tiny adapter is the one place that
    difference is bridged, so ImagePlot/XYPlot's own calling convention
    never needs to special-case "statusbar vs label". """

    def __init__ (self, statusbar):
        self.statusbar  = statusbar
        self.context_id = statusbar.get_context_id ( "plot_hover" )

    def set_text (self, text):
        self.statusbar.pop  ( self.context_id )
        self.statusbar.push ( self.context_id, text )


def _export_plot_or_warn (plot_widget, filepath, scale = 4):
    """ [EN] User's own request: "ao gerar o pmf, vamos automaticamente
    gerar duas figuras no disco, o plot do pmf e o plot das gaussianas."
    Thin wrapper around util.easyplot.export_utils.export_plot_to_png()
    (already existed, already used elsewhere for on-demand high-res
    export -- see that module's own docstring for why it re-runs on_draw()
    into a bigger offscreen surface rather than screenshotting the live
    widget) -- just adds a try/except so a failed PNG write (e.g. the
    widget wasn't actually allocated a real size yet) can't take down the
    "WHAM finished successfully" flow around it; the WHAM results windows
    themselves are still shown either way, this is a best-effort side
    file, not something the user is blocked on. """
    from util.easyplot.export_utils import export_plot_to_png
    try:
        export_plot_to_png ( plot_widget, filepath, scale = scale )
    except Exception as exc:
        dprint ( "WARNING: could not export plot to '{}': {}".format ( filepath, exc ) )
    return GLib.SOURCE_REMOVE   # [EN] run exactly once via GLib.idle_add -- see call sites below;
                                 # returning True (GLib.SOURCE_CONTINUE) here would make GLib keep
                                 # re-invoking this same export forever, on every idle tick.


def _pack_plot_with_statusbar (window, plot_widget):
    """ [EN] Shared by every WHAM results window below (histograms/PMF/
    heatmap): wraps `plot_widget` in a vertical GtkBox with a GtkStatusbar
    docked at the bottom, and returns a _StatusbarHoverLabel already
    wired as `plot_widget.RC_label` -- so a caller just does
    `_pack_plot_with_statusbar(window, self.plot)` instead of repeating
    this same box/statusbar/adapter boilerplate 3 times. """
    box = Gtk.Box ( orientation = Gtk.Orientation.VERTICAL )
    box.pack_start ( plot_widget, True, True, 0 )
    statusbar = Gtk.Statusbar ( )
    box.pack_start ( statusbar, False, False, 0 )
    window.add ( box )
    plot_widget.RC_label = _StatusbarHoverLabel ( statusbar )
    return plot_widget.RC_label


class pAnalysis:
    """ Class doc """

    def __init__ (self):
        """ Class initialiser """
        self.imgPlot = None

    def on_mouse_button_press (self, widget, event):
        dprint (widget, event )



    def run_analysis (self, parameters):
        """ Function doc """
        if parameters['analysis_type'] == 'wham':
            wham_analysis = WHAMAnalysis()
            TrueFalse, results = wham_analysis.run ( parameters )
            
            
            if TrueFalse:
                self.main.simple_dialog.info(msg = results['msg'] )
                
                if results['type'] == 0:
                    from util.easyplot import ImagePlot, XYPlot
                    import random
                    
                    '''                 Histograms                 '''
                    #self.plot = XYPlot(bg_color = [0,0,0])
                    self.plot = XYPlot( )
                    
                    for i , log in enumerate(results['histograms']):
                        
                        X = [] 
                        Y = [] 
                        
                        r = random.random()
                        g = random.random()
                        b = random.random()
                        rgb = [r,g,b,]
                        data = open(log, 'r')
                        for line in data:
                            line2 = line.split()
                            X.append(float(line2[0]) )
                            Y.append(float(line2[1]))
                        self.plot.add ( X = X, Y = Y,
                                        symbol = None, sym_color = [1,1,1], sym_fill = False,
                                        line = 'solid', line_color = rgb, energy_label = None,
                                        label = os.path.basename ( log ) )

                    #self.plot.Ymax_list= [100]
                    window =  Gtk.Window()
                    window.set_default_size(800, 300)
                    window.move(900, 300)
                    window.set_title('Histograms')
                    _pack_plot_with_statusbar ( window, self.plot )

                    # [EN] User's own request -- "ao clicar no grafico
                    # (sobre uma das gaussianas) ele reconhece qual
                    # gaussiana/janela foi selecionada": curve index i ==
                    # results['histograms'][i], the SAME order they were
                    # add()'ed in just above, so the clicked curve's own
                    # 'label' (its source file's basename) already
                    # identifies the umbrella window directly -- no extra
                    # index-to-window bookkeeping needed. Shown in the
                    # WINDOW TITLE rather than the hover statusbar: the
                    # statusbar already gets overwritten by the very next
                    # mouse-motion event's own x/y readout, which would
                    # make a "you selected window N" message disappear
                    # again the instant the cursor moves off the line --
                    # the title bar persists until the next click instead.
                    def _on_histogram_curve_clicked ( index, label, window = window ):
                        if index is None:
                            window.set_title ( 'Histograms' )
                        else:
                            window.set_title ( 'Histograms -- selected: {}'.format ( label ) )
                    self.plot.on_curve_click_callback = _on_histogram_curve_clicked

                    window.show_all()

                    # [EN] User's own request -- auto-save a PNG of this
                    # plot to disk the moment WHAM finishes, no manual
                    # export click needed. Same folder + logfile-name
                    # prefix WHAMAnalysis.run() itself already uses for
                    # every other output file (results['pmf']/['histograms'],
                    # see wham.py's own PMF_file/output_file naming), so
                    # the PNGs sit right next to the data they were
                    # plotted from. Deferred one GLib idle tick past
                    # show_all() -- export_plot_to_png() needs the widget's
                    # OWN get_allocated_width/height() to already be real
                    # (not 0x0), which GTK only guarantees once the normal
                    # size-allocate cycle for this brand-new window has
                    # actually run, not necessarily synchronously inside
                    # show_all() itself.
                    histograms_png = os.path.join ( parameters['folder'], parameters['logfile'] + '_histograms.png' )
                    GLib.idle_add ( _export_plot_or_warn, self.plot, histograms_png )

                    X = []
                    Y = []
                    
                    '''                   PMF plot                  '''
                    #self.plot2 = XYPlot(bg_color = [0,0,0])
                    self.plot2 = XYPlot( )
                    dprint(results['pmf'])
                    data2 = open(results['pmf'], 'r')
                    for line in data2:
                        line2 = line.split()
                        if line2[1] == 'inf':
                            pass
                        else:
                            X.append(float(line2[0]) )
                            Y.append(float(line2[1]))

                    
                    
                    self.plot2.add ( X = X, Y = Y,
                                    symbol = 'dot', sym_color = [0,0,0], sym_fill = False, 
                                    line = 'solid', line_color = [0,0,0], energy_label = None)
                    
                    window2 =  Gtk.Window()
                    window2.set_default_size(800, 300)
                    window2.move(100, 300)
                    window2.set_title(results['pmf'])
                    _pack_plot_with_statusbar ( window2, self.plot2 )
                    window2.show_all()

                    pmf_png = os.path.join ( parameters['folder'], parameters['logfile'] + '_pmf.png' )
                    GLib.idle_add ( _export_plot_or_warn, self.plot2, pmf_png )
                    
                if results['type'] == 1:
                    from util.easyplot import ImagePlot, XYPlot
                    data = open(results['pmf'], 'r')
                    self.imgPlot = ImagePlot()
                    self.imgPlot.connect("button_press_event", self.on_mouse_button_press)
                    # [EN] No external motion-notify-event connect needed --
                    # ImagePlot.__init__() already self-connects its OWN
                    # on_motion() (image_plot.py), which already pushes a
                    # live "i | j | rc1 | rc2 | E" readout into self.
                    # RC_label whenever it's set (see _pack_plot_with_
                    # statusbar() below) -- this used to ALSO connect
                    # pAnalysis.on_motion() here, a near-duplicate of that
                    # same logic that only ever dprint()'d it (RC_label was
                    # never set), so both handlers fired on every motion
                    # event for no visible benefit. Removed instead of kept
                    # alongside, per the user's own request to make this
                    # live readout actually visible in the UI.
                    X =[]
                    Y =[]
                    Z =[]
                    for line in data:
                        line2 = line.split()

                        X.append(float(line2[0]) )
                        Y.append(float(line2[1]))
                        Z.append(float(line2[2]))


                    sizex = X.count(X[0])
                    sizey = Y.count(Y[0])

                    RC1 = []
                    RC2 = []
                    _Z  = []
                    for i in range(sizex):
                        RC1.append(X[sizey*i:sizey*(i+1)])
                        _Z.append( Z[sizey*i:sizey*(i+1)] )
                            #print(i,j, X[sizey*i:sizey*(i+1)])

                    for j in range(len(RC1)):
                        RC2.append(  Y[sizey*j:sizey*(j+1)]  )

                    #print(len(RC1),len(RC1[1]) )
                    #print(len(RC2),len(RC2[1]) )
                    #
                    #print(RC1[1])
                    #print(RC2[1])
                    #print(_Z )
                    data = {
                            'name': 'output.log', 
                            'type': 'self.imgPlot2D', 
                            'RC1': RC1,
                            'RC2': RC2,
                            'Z'  : _Z,
                           }

                    self.imgPlot.show()
                    self.imgPlot.data    =  data['Z']
                    self.imgPlot.dataRC1 =  data['RC1']
                    self.imgPlot.dataRC2 =  data['RC2']
                    self.imgPlot.set_label_mode(mode = 1)

                    #self.imgPlot.set_threshold_color ( _min = 0, _max = 100, cmap = 'jet')
                    window =  Gtk.Window()
                    window.set_default_size(800, 300)
                    window.move(900, 300)
                    window.set_title(results['pmf'])
                    _pack_plot_with_statusbar ( window, self.imgPlot )
                    window.show_all()
                    
                    

            else:
                self.main.simple_dialog.error(msg = results['msg'] )
                #dialog = Gtk.MessageDialog(
                #        parent=self.main.window,
                #        flags=Gtk.DialogFlags.MODAL,
                #        type=Gtk.MessageType.ERROR,
                #        buttons=Gtk.ButtonsType.OK,
                #        message_format=msg
                #    )
                #dialog.run()
                #dialog.destroy()

            
            #print (TrueFalse,msg )
            return TrueFalse
