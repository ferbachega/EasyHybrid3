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

class Atom:
    """ Class doc """
    
    def __init__(self, vismol_object, name="Xx", index=None, residue=None,
                 chain=None, pos=None, symbol=None, atom_id=None, color=None,
                 vdw_rad=None, cov_rad=None, ball_rad=None, electronegativity=None,
                 occupancy=0.0, bfactor=0.0, charge=0.0, bonds_indexes=None):
        """ Class initializer """
        self.vm_object   = vismol_object
        self.vm_session  = vismol_object.vm_session
        
        
        self.name     = name
        self.index    = index   # - Remember that the "index" attribute refers to the numbering of atoms (it is not a zero base, it starts at 1 for the first atom)
        self.residue  = residue
        self.chain    = chain
        self.molecule = None

        self.pos = pos     # - coordinates of the first frame
        self.unique_id = None
        
        # self.symbol = self._get_symbol(self.name) if (symbol is None) else symbol
        if symbol is None:
            self.symbol = self._get_symbol()
        else:
            self.symbol = symbol
        self.atom_id = atom_id
        
        if color is None:
            self.color = self._init_color()
        else:
            self.color = color
        
        if vdw_rad is None:
            self.vdw_rad = self._init_vdw_rad()
        else:
            self.vdw_rad = vdw_rad
        
        if cov_rad is None:
            self.cov_rad = self._init_cov_rad()
        else:
            self.cov_rad = cov_rad
        
        if ball_rad is None:
            self.ball_rad = self._init_ball_rad()
        else:
            self.ball_rad = ball_rad

        #Bachega April/05/2026
        if electronegativity is None:
            self.electronegativity = self._init_electronegativity()
        else:
            self.electronegativity = electronegativity           

        self.color_id = None
        self.occupancy = occupancy
        self.bfactor = bfactor
        self.charge = charge
        if bonds_indexes is None:
            self.bonds_indexes = []
        else:
            self.bonds_indexes = bonds_indexes
        
        self.selected       = False
        self.lines          = True
        self.dots           = False
        self.nonbonded      = False
        self.impostor       = False
        self.ribbons        = False
        self.ribbon_sphere  = False
        self.dynamic        = False
        self.ball_and_stick = False
        self.sticks         = False
        self.stick_spheres  = False
        self.spheres        = False
        self.vdw_spheres    = False
        self.dash           = False
        self.surface        = False
        self.bonds          = []
        self.isfree         = True
        self.labels         = False
        self.label_text     = None
        
    def _get_symbol(self):
        """ Function doc """
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        
        name = self.name.strip()
        if name == "":
            return ""
        
        _n = name
        for char in name:
            if char.isnumeric():
                _n = _n.replace(char, "")
        name = _n
        
        if len(name) >= 3:
            name = name[:2]
        
        # This can't happen since before all numbers have been converted to strings
        # if len(name) == 2:
        #     if name[1].isnumeric():
        #         symbol = name[0]
        
        # The capitalization of hte name can solve all the next ifs
        # name = name.lower().capitalize()
        if name in ATOM_TYPES.keys():
            return name
        else:
            if name[0] == "H":
                if name[1] == "g":
                    symbol =  "Hg"
                elif name[1] =="e":
                    symbol = "He"
                else:
                    symbol =  "H"
            
            elif name[0] == "C":
                if name[1] == "a":
                    symbol = "Ca"
                elif name[1] =="l":
                    symbol = "Cl"
                elif name[1] =="L":
                    symbol = "Cl"
                elif name[1] =="d":
                    symbol = "Cd"
                elif name[1] =="u":
                    symbol = "Cu"
                elif name[1] =="U":
                    symbol = "Cu"
                else:
                    symbol = "C"
            
            elif name[0] == "N":
                if name[1] == "i" or name[1] == "I":
                    symbol = "Ni"
                elif name[1] == "a":
                    symbol = "Na"
                elif name[1] == "e":
                    symbol = "Ne"
                elif name[1] == "b":
                    symbol = "Nb"
                else:
                    symbol = "N"
            
            elif name[0] == "O":
                if name[1] == "s":
                    symbol = "Os"
                else:
                    symbol = "O"
            
            elif name[0] == "S":
                if name[1] == "I":
                    symbol = "Si"
                elif name[1] == "e":
                    symbol = "Se"
                
                elif name[1] == "O":
                    try:
                        if name[2] == "D":
                            symbol = "Na"
                        else:
                            pass
                    except:
                        symbol = "S"
                else:
                    symbol = "S"
            
            elif name[0] == "P":
                if name[1] == "d":
                    symbol = "Pd"
                elif  name[1] == "b":
                    symbol = "Pb"
                elif  name[1] == "o":
                    symbol = "Po"
                else:
                    symbol = "P" 
            
            elif name[0] == "Z":
                if name[1] == "r":
                    symbol = "Zr"
                elif  name[1] == "N":
                    symbol = "Zn"
                else:
                    symbol = "Zn" 
            
            elif name[0] == "F":
                if name[1] == "E":
                    symbol = "Fe"
                elif  name[1] == "e":
                    symbol = "Fe"
                else:
                    symbol = "F" 
            
            elif name[0] == "M":
                if name[1] == "n":
                    symbol = "Mn"
                elif name[1] == "N":
                    symbol = "Mn"
                elif name[1] == "o":
                    symbol = "Mo"
                elif name[1] == "G":
                    symbol = "Mg"
                else:
                    symbol = "X"
            else:
                symbol = "X"
        return symbol
    
    def _init_color(self):
        """ Return the color of an atom in RGB. Note that the returned
            value is in scale of 0 to 1, but you can change this in the
            index. If the atomname does not match any of the names
            given, it returns the default dummy value of atom X.
        """
        try:
            color = self.vm_object.color_palette[self.name]
        except KeyError:
            # #print(self.symbol, "Atom")
            color = self.vm_object.color_palette[self.symbol]
        return np.array(color, dtype=np.float32)
    
    def _generate_atom_unique_color_id(self):
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        """ Function doc """
        r = (self.unique_id & 0x000000FF) >> 0
        g = (self.unique_id & 0x0000FF00) >> 8
        b = (self.unique_id & 0x00FF0000) >> 16
        # pickedID = r + g * 256 + b * 256*256
        self.color_id = np.array([r/255.0, g/255.0, b/255.0], dtype=np.float32)
    
    def _init_vdw_rad(self):
        """ Function doc """
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        try:
            vdw = ATOM_TYPES[self.name][6]
        except KeyError:
            vdw = ATOM_TYPES[self.symbol][6]
        return vdw
    
    def _init_cov_rad(self):
        """ Function doc """
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        try:
            cov = ATOM_TYPES[self.name][5]
        except KeyError:
            cov = ATOM_TYPES[self.symbol][5]
        return cov
        #try:
        #    cov = ATOM_TYPES[self.name][7]
        #except KeyError:
        #    cov = ATOM_TYPES[self.symbol][7]
        #return cov
    
    def _init_ball_rad(self):
        """ Function doc """
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        try:
            ball = ATOM_TYPES[self.name][6]
        except KeyError:
            ball = ATOM_TYPES[self.symbol][6]
        return ball
    
    def coords(self, frame=None):
        """ 
        frame = int
        
        Returns the coordinates of an atom according to the specified frame. 
        If no frame is specified, the frame set by easyhybrid (probably by the 
        scale bar of the trajectory manipulation window) is used. If the object 
        (vobject) has a smaller number of frames than the one set by the 
        interface, the last frame of the object is used.
        
        return  xyz 
        """
        if frame is None:
            frame  = self.vm_object.vm_session.frame
            #print (frame, len(self.vm_object.frames))
            if len(self.vm_object.frames)-1 <= frame:
                frame = len(self.vm_object.frames)-1
            else:
                pass
        return self.vm_object.frames[frame, self.atom_id]
    
    def get_grid_position(self, gridsize=3, frame=None):
        """ Function doc """
        coords = self.coords(frame)
        gridpos = (int(coords[0]/gridsize), int(coords[1]/gridsize), int(coords[2]/gridsize))
        return gridpos
    
    def get_cov_rad(self):
        """ Function doc """
        return self.cov_rad
    
    def init_color_rgb(self, name):
        """ Return the color of an atom in RGB. Note that the returned
            value is in scale of 0 to 1, but you can change this in the
            index. If the atomname does not match any of the names
            given, it returns the default dummy value of atom X.
        """
        
        try:
            color = color =self.vm_object.color_palette[name]
            #color = ATOM_TYPES[name][1]
        except KeyError:
            if name[0] == "H":# or name in self.hydrogen:
                #color = ATOM_TYPES["H"][1]
                color = self.vm_object.color_palette["H"]
            
            elif name[0] == "C":
                #color = ATOM_TYPES["C"][1]
                color = self.vm_object.color_palette["C"]
            
            elif name[0] == "O":
                #color = ATOM_TYPES["O"][1]
                color = self.vm_object.color_palette["O"]
            
            elif name[0] == "N":
                #color = ATOM_TYPES["N"][1]
                color = self.vm_object.color_palette["N"]
                
            elif name[0] == "S":
                #color = ATOM_TYPES["S"][1]
                color = self.vm_object.color_palette["S"]
            else:
                #color = ATOM_TYPES["X"][1]
                color = self.vm_object.color_palette["X"]
                
        color = [int(color[0]*250), int(color[1]*250), int(color[2]*250)]
        return color


    def _init_electronegativity(self):
        """ Function doc """
        #return False
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        try:
            en_UFF = ATOM_TYPES[self.name][8]
        except KeyError:
            try:
                en_UFF = ATOM_TYPES[self.symbol][8]
            except:
                en_UFF = 1.0
        return en_UFF

    def init_radius(self, name):
        """
        """
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        try:
            rad = ATOM_TYPES[name][6]/5.0
        except KeyError:
            if name[0] == "H" or name in self.hydrogen:
                rad = ATOM_TYPES["H"][6]/5.0
            elif name[0] == "C":
                rad = ATOM_TYPES["C"][6]/5.0
            elif name[0] == "O":
                rad = ATOM_TYPES["O"][6]/5.0
            elif name[0] == "N":
                rad = ATOM_TYPES["N"][6]/5.0
            elif name[0] == "S":
                rad = ATOM_TYPES["S"][6]/5.0
            else:
                rad = 0.30
        return rad


def generate_random_code(length):
    # Define the character set: uppercase letters, lowercase letters, and digits
    characters = string.ascii_letters + string.digits
    # Generate the random code by randomly selecting characters from the set
    random_code = ''.join(random.choice(characters) for _ in range(length))
    return random_code


def export_special_PDB (vobject = None, frame = -1, output = 'temp.pdb'):
    """ Function doc """
    
    
    data = open(output, 'w')
    change_chain = False
    text = ''
    for index in vobject.atoms.keys():
        atom    = vobject.atoms[index]
        name    = atom.name
        
        resi    = atom.residue
        chain   = atom.chain
        chain_name = chain.name
        
        if chain_name =='':
            if change_chain:
                chain_name = 'B'
            else:
                chain_name = 'A'
        #resn    = vobject.chain[chain].
        
        #chain   = atom.chain
        
        symbol  = atom.symbol
        coords  = atom.coords(frame)
        
        resi_index =resi.index
        
        if resi_index > 9999:
            change_chain = True
            resi_index = resi_index-9999
        
        #print (atom, name, resi, chain)
        #           ATOM    120  CA  VAL A  17      36.365  31.982  42.405  1.00 53.96           C  
        #           ATOM    514 H514 SER  34        17.413  42.503  16.453  1.00  1.00      H   
        #           HETATM63640  H1  WAT  1916       7.005  19.149   8.699  0.00  0.00          H   
        #           ATOM    116 S116 CYS  9         23.989  35.368   6.299  1.00  1.00      S   
        #           HETATM63640  H1         WAT  1916       7.005  19.149   8.699  0.00  0.00          H   
        ATOMLINE = "{:<6s}{:5d} {:<4s} {:3s} {:1s}{:>4s}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:<4s}\n".format(
        #HETATM
        'ATOM  ',
        index+1,
        str(name),
        str(resi.name),
        str(chain_name),
        str(resi_index),
        coords[0],
        coords[1],
        coords[2],
        1.0,
        1.0,
        symbol
        )
        #print(ATOMLINE)
        text += ATOMLINE
    data.write(text)
