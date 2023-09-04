#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  relaxed_scaning.py
#  
#  Copyright 2023 Fernando Bachega <fernando@fernando-Inspiron-7537>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
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
#  

import sys

'''
This is a basic script to open an already prepared system in pDynamo. 
It simply consists of loading a pkl file and doing a quick system 
assessment.
'''

from pBabel                    import *                                     
from pCore                     import *                                     
from pMolecule                 import *                  
from pScientific               import *                                     
from pScientific.Arrays        import *                                     
from pScientific.Geometry3     import *                 
from pSimulation               import *

sys.path.insert(1, '/home/fernando/programs/EasyHybrid3/pdynamo')
import p_methods  as pMethods




# importing system
system = ImportSystem ('/home/fernando/programs/EasyHybrid3/MOLECULES.pkl')
# summary
system.Summary()



parameters = {'NmaxThreads'          : 1,                  # Use 1 for one-dimensional reaction coordinate scans 
                                                           
             'RC1': {'ATOMS'         : [0, 3],             # atom indices. For single-distance two atoms are needed
                     'ATOM_NAMES'    : ['O0', 'N3'],       # atom names.
                     
                     'dincre'        : 0.1,                # distance increment
                     'dminimum'      : 4.208733002805674,  # starting distance
                     'force_constant': 4000.0,             # Force constant used in harmonic restraint
                     'nsteps'        : 5,                  #
                     'rc_type'       : 'simple_distance',  # simple_distance / multiple-distance
                     
                     'sigma_pk1pk3'  : None,               # Required only when multiple-distance is used
                     'sigma_pk3pk1'  : None},              # Required only when multiple-distance is used
             
             'RC2': None,
             
             'folder'             : '/home/fernando/programs/EasyHybrid3/examples/scripts',
             'initial_coordinates': None,
             
             'log_frequency'      : 50,
             'maxIterations'      : 600.0,
             'optimizer'          : 'ConjugatedGradient',
             'rmsGradient'        : 0.1,
             'simulation_type'    : None,
             'system'             : system,
             'traj_folder_name'   : '1-MOL_AM1_QC10_PES_scan',
             'traj_type'          : 'pklfolder',
             'vobject_name'       : None
             }

simObj =  pMethods.RelaxedSurfaceScan()
simObj.run(parameters)



