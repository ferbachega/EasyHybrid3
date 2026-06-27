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

class Restraints:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
    
    def update_restaint_representation (self, e_id = None):
        """ Function doc """
        try:
            if e_id:
                pass
            else:
                e_id = self.active_id
            
            indexes = []
            for name, restraint in self.psystem[e_id].e_restraints_dict.items():
                _bool = restraint[0]
                name  = restraint[1]
                _type = restraint[2]
                if _type == 'distance':
                    atons = '{} / {}'.format(restraint[3][0],restraint[3][1]) 
                    dist_or_angle = '{:.4f}'.format(restraint[4])
                    force_const   = str(restraint[5])
                    e_id          =  restraint[6] 
                    if _bool:
                        indexes.append(restraint[3][0])
                        indexes.append(restraint[3][1])


            for vobject in self.main.vm_session.vm_objects_dic.values():
                
                if vobject.e_id == e_id:
                    
                    if indexes == []:
                        
                        vobject.representations["restraints"] = None
                    
                    else:
                        
                        if 'restraints' in vobject.representations.keys():
                            
                            if vobject.representations["restraints"] is not None:
                                #print('["restraints"] is not None',indexes)
                                #vobject.representations["restraints"].define_new_indexes_to_vbo(indexes)
                                vobject.representations["restraints"] = DashedLinesRepresentation(vobject, self.vm_session.vm_glcore,
                                                                                          active=True, indexes = indexes)
                            else:
                                #print('else',indexes)
                                vobject.representations["restraints"] = DashedLinesRepresentation(vobject, self.vm_session.vm_glcore,
                                                                                          active=True, indexes = indexes)
                        else:
                            vobject.representations["restraints"] = DashedLinesRepresentation(vobject, self.vm_session.vm_glcore,
                                                                      active=True, indexes = indexes)    
            self.main.vm_session.vm_glcore.queue_draw()
            #print(indexes)
    
        except:
            print('\n\n Failed when trying to represent harmonic potential constraints. This is just a representation error, the potencies are working normally.\n\n')
        
    def add_new_harmonic_restraint (self, parameters, _type = 'distance'):
        """ Function doc """
        #restraints = RestraintModel()
        #parameters['system'].DefineRestraintModel( restraints )
        
        if _type == 'distance':
            
            atom1 = parameters['atom1']
            atom2 = parameters['atom2']
           
            rest_name = str(parameters['system'].e_restraint_counter)
            parameters['system'].e_restraints_dict[rest_name] = [True, 
                                                                 rest_name ,
                                                                 'distance', 
                                                                 [parameters['atom1'],parameters['atom2']], 
                                                                 parameters['distance'],  
                                                                 parameters['force_constant'],
                                                                 parameters['system'].e_id] 
            
            
            parameters['system'].e_restraint_counter += 1
        
        else:
            pass

    def define_harmonic_restraints (self, parameters):
        """ Function doc """
        parameters['equilibriumValue'] = 0.0
        parameters['forceConstant']    = 0.0
        parameters['reference']        = None # <- Clone ( system.coordinates3 )
        parameters['period']           = None
        parameters['selection']        = None
                
        tetherEnergyModel  = RestraintEnergyModel.Harmonic ( 0.0, forceConstant )
        tethers            = RestraintModel ( )        
        tethers["Tethers"] = RestraintMultipleTether.WithOptions ( energyModel = tetherEnergyModel ,
                                                                   reference   = reference         , 
                                                                   selection   = heavies           )
        system.DefineRestraintModel ( tethers )

    def define_distance_harmonic_restraints (self, parameters):
        """ Function doc """
        restraints = RestraintModel()
        #print('\n\n\n\n',parameters, '\n\n\n\n')
        parameters['system'].DefineRestraintModel( restraints )
        
        distance = parameters['dminimum_RC1'] + ( parameters['dincre_RC1'] * float(i) )
        
        rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['force_constant_1'])
        restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom1, point2= atom2)
        restraints["RC1"] = restraint

    def define_position_harmonic_restraints (self, parameters):
        """ Function doc """
        # . parameters['reference']
        # . parameters['selection']
        # . parameters['system']
        # . parameters['force_constant']

        # . Harmonically restrain heavy atoms.
        tethers       = None
        forceConstant = parameters['force_constant']
        
        reference          = Clone ( parameters['reference'] ) # Change it later
        tetherEnergyModel  = RestraintEnergyModel.Harmonic ( 0.0, forceConstant )
        tethers            = RestraintModel ( )
        tethers["Tethers"] = RestraintMultipleTether.WithOptions ( energyModel = tetherEnergyModel ,
                                                                   reference   = reference         , 
                                                                   selection   = parameters['selection'])
        parameters['system'].DefineRestraintModel ( tethers )



    def clear_restraints (self):
        """ Function doc """
        parameters['system'].DefineRestraintModel(None)
