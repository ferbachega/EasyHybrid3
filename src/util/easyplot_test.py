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
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk
import cairo
import random
import numpy as np
import sys
from pprint import pprint
from  easyplot import ImagePlot, XYPlot, XYScatterPlot
import random

'''                 Histograms                 '''
#plot = XYPlot( )
plot = ImagePlot()

#data = open('/home/fernando/PROJECTS/water_trimer/output_2D_pmf.log', 'r')
#data = open('/home/fernando/PROJECTS/water_trimer/PMFs2D/output_50_35_pmf.log', 'r')
#data = open('/home/fernando/PROJECTS/CH3Cl_Br/output2d_pmf.log', 'r')
#data = open('/home/fernando/PROJECTS/CH3Cl_Br/output_pmf.log', 'r')
#data = open('/home/fernando/PROJECTS/CH3Cl_Br/output2d_50_50_pmf.log', 'r')
data = open('/home/fernando/Downloads/3qva_xtb_etapa1_2d_pmf(1).log', 'r')
#data = open('/home/fernando/Downloads/3qva_2d_pmf_pmf.log', 'r')
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
        'type': 'plot2D', 
        'RC1': RC1,
        'RC2': RC2,
        'Z'  : _Z,
       }

plot.show()
plot.data    =  data['Z']
plot.dataRC1 =  data['RC1']
plot.dataRC2 =  data['RC2']
plot.set_label_mode(mode = 1)

window =  Gtk.Window()
window.set_default_size(800, 300)
window.move(900, 300)
window.set_title('Histograms')
window.add(plot)
window.show_all()
Gtk.main()
#'''
