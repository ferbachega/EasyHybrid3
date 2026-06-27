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

class pAnalysis:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        self.imgPlot = None
        
    def on_mouse_button_press (self, widget, event):
        print (widget, event )
    
    def on_motion (self, widget, event):
        """ Function doc """
        #print(widget, event)

        x_on_plot, y_on_plot, x, y = self.imgPlot.get_i_and_j_from_click_event (event)
        i = y_on_plot
        j = x_on_plot
        
        j_size = len(self.imgPlot.norm_data)
        i_size = len(self.imgPlot.norm_data[0])
        
        if i < 0 or j < 0:
            pass
        else:
            if i < i_size and j < j_size:
                text = 'i {}    |    j {}    |    rc1 {:4.2f}    |    rc2 {:4.2f}    |    E = {:6.3f}'.format(i, j,  
                                                                                         self.imgPlot.dataRC1[j][i], 
                                                                                         self.imgPlot.dataRC2[j][i],  
                                                                                         self.imgPlot.data[j][i])
                print(text)
                    #self.RC_label.set_text(text)






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
                                        line = 'solid', line_color = rgb, energy_label = None)
                    
                    #self.plot.Ymax_list= [100]
                    window =  Gtk.Window()
                    window.set_default_size(800, 300)
                    window.move(900, 300)
                    window.set_title('Histograms')
                    window.add(self.plot)
                    window.show_all()
                    
                    X = []
                    Y = []
                    
                    '''                   PMF plot                  '''
                    #self.plot2 = XYPlot(bg_color = [0,0,0])
                    self.plot2 = XYPlot( )
                    print(results['pmf'])
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
                    window2.add(self.plot2)
                    window2.show_all()
                    
                if results['type'] == 1:
                    from util.easyplot import ImagePlot, XYPlot
                    data = open(results['pmf'], 'r')
                    self.imgPlot = ImagePlot()
                    self.imgPlot.connect("button_press_event", self.on_mouse_button_press)
                    self.imgPlot.connect("motion-notify-event", self.on_motion)
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
                    window.add(self.imgPlot)
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
