#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  simple.py
#
#  Copyright 2022 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
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

import logging
import gi, sys
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
from vismol.core.vismol_object import VismolObject
import time


from vismol.model.atom import Atom
from vismol.model.residue import Residue
from vismol.model.chain import Chain
from vismol.core.vismol_object import VismolObject


class EVismolObject(VismolObject):
    """ Class doc """
    
    def __init__ (self, vismol_session, index, name="UNK", active=False, trajectory=None,
                 color_palette=None, bonds_pair_of_indexes=None):
        """ Class initialiser """
        super().__init__(self, vismol_session, index, name="UNK", active=False, trajectory=None,
                 color_palette=None, bonds_pair_of_indexes=None)


    def _generate_topology_from_atom_list (self, atoms = []):
        """ Function doc """
        initial   = time.time()
        unique_id = len(self.vm_session.atom_dic_id)
        atom_id   = 0
        
        for _atom in atoms:
            if _atom["chain"] not in self.chains.keys():
                self.chains[_atom["chain"]] = Chain(self, name=_atom["chain"])
            _chain = self.chains[_atom["chain"]]
            
            if _atom["resi"] not in _chain.residues.keys():
                _r = Residue(self, name=_atom["resn"], index=_atom["resi"], chain=_chain)
                self.residues[_atom["resi"]] = _r
                _chain.residues[_atom["resi"]] = _r
            _residue = _chain.residues[_atom["resi"]]
            
            #atom = Atom()
            atom = Atom(vismol_object = self,#
                        name          = _atom["name"] , 
                        index         = _atom["index"],
                        residue       = _residue      , 
                        chain         = _chain, 
                        atom_id       = atom_id,
                        occupancy     = _atom["occupancy"], 
                        bfactor       = _atom["bfactor"],
                        charge        = _atom["charge"])
                        
                        
            atom.unique_id = unique_id
            atom._generate_atom_unique_color_id()
            _residue.atoms[atom_id] = atom
            self.atoms[atom_id] = atom
            atom_id += 1
            unique_id += 1
