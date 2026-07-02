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
from vismol.libgl.representations import PickingDotsRepresentation
from vismol.libgl.representations import DotsRepresentation
from vismol.libgl.representations import OneColorDotsRepresentation

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
    """Restraints mixin for pDynamoSession.

    Groups the methods that: (1) keep the visual REPRESENTATION of distance
    restraints up to date (the dashed lines drawn between restrained atoms) and
    (2) DEFINE the restraint models on the pDynamo 'system' object used by the
    simulation.

    It helps to keep the two coexisting layers apart:
      - e_restraints_dict / e_restraint_counter / e_id  -> EasyHybrid layer,
        used by the INTERFACE and the REPRESENTATION (what shows up in the 3D view).
      - RestraintModel / RestraintEnergyModel / RestraintDistance / ...  -> pDynamo
        layer, used in the actual ENERGY calculation during the simulation.
    """

    def __init__ (self):
        """Initialiser. Intentionally empty: this is a MIXIN.

        The class is never instantiated on its own; it is inherited by
        pDynamoSession, whose __init__ is what actually sets up the attributes
        (self.psystem, self.main, self.active_id, self.vm_session, etc.) that the
        methods below rely on through 'self'.
        """

    def update_restaint_representation_old (self, e_id = None):
        """Update the on-screen drawing of the distance restraints in the 3D view.

        Walks the restraints registered for a system and (re)creates the
        'dashed lines' representation connecting the pairs of atoms that are
        under an active restraint. It does not change the physics; purely visual.

        Args:
            e_id: id of the target system. If None, uses the active system
                  (self.active_id).
        """
        #try:
        # ---- Resolve which system will be represented ----------------------
        if e_id:
            # e_id explicitly provided: use the one passed in.
            pass
        else:
            # no e_id provided -> fall back to the currently active system.
            e_id = self.active_id

        # 'indexes' collects the indices of the atoms that must be joined by a
        # dashed line (two indices per active distance restraint).
        indexes = []
        pos_index = []
        # ---- Scan the restraints registered for this system ----------------
        # e_restraints_dict maps: restraint_name -> list with the fields
        # [active?, name, type, [atoms], value, force_constant, e_id].
        for name, restraint in self.psystem[e_id].e_restraints_dict.items():
            _bool = restraint[0]   # bool: is the restraint active?
            name  = restraint[1]   # restraint name/identifier
            _type = restraint[2]   # type: only 'distance' is handled here

            if _type == 'distance':
                # Fields specific to a distance restraint.
                # (atons/dist_or_angle/force_const/e_id are built here, but
                #  note that only 'e_id' is actually used below; the others
                #  look like leftovers for display/debug and do not affect
                #  the drawing.)
                atons = '{} / {}'.format(restraint[3][0], restraint[3][1])  # atom pair (text)
                dist_or_angle = '{:.4f}'.format(restraint[4])               # target distance (text)
                force_const   = str(restraint[5])                           # force constant (text)
                e_id          = restraint[6]                                # id of the system owning the restraint

                # Only draw the restraints that are ACTIVE (_bool == True).
                if _bool:
                    indexes.append(restraint[3][0])  # 1st atom of the pair
                    indexes.append(restraint[3][1])  # 2nd atom of the pair
            
            elif _type == 'position':
                pos_indexes = restraint[3]
                pass
            else:
                pass
                
        # ---- Apply the representation to the matching visual object ---------
        # Look, among the objects in the 3D scene, for the one whose e_id
        # matches the system of the restraints we just collected.
        for vobject in self.main.vm_session.vm_objects_dic.values():

            if vobject.e_id == e_id:

                if indexes == []:
                    # No active restraints -> drop the representation
                    # (there is nothing to draw).
                    
                    vobject.representations["restraints"] = None
                
                elif pos_indexes == []:
                    
                    vobject.representations["position"] = None

                else:
                    # There are atoms to connect. The three branches below all
                    # end up doing the SAME thing: creating a brand-new
                    # DashedLinesRepresentation with the collected indices. The
                    # commented-out 'define_new_indexes_to_vbo' call hints that
                    # the future intent was to merely UPDATE the indices of the
                    # existing object (cheaper) instead of recreating it, but
                    # that is currently disabled.
                    #if _type == 'restraints':
                    if 'restraints' in vobject.representations.keys():
                        
                        if vobject.representations["restraints"] is not None:
                            # A representation already existed: recreate it from scratch.
                            #vobject.representations["restraints"].define_new_indexes_to_vbo(indexes)
                            vobject.representations["restraints"] = DashedLinesRepresentation(
                                vobject, self.vm_session.vm_glcore,
                                active=True, indexes=indexes)
                        else:
                            # The key exists but was None: create it.
                            vobject.representations["restraints"] = DashedLinesRepresentation(
                                vobject, self.vm_session.vm_glcore,
                                active=True, indexes=indexes)
                    else:
                        # The 'restraints' key did not even exist yet: create it.
                        vobject.representations["restraints"] = DashedLinesRepresentation(
                            vobject, self.vm_session.vm_glcore,
                            active=True, indexes=indexes)
                    
                    
                    #if _type == 'position':
                    if 'position' in vobject.representations.keys():
                        #DotsRepresentation(self, self.vm_session.vm_glcore,
                        #                    active=True, indexes=list(self.atoms.keys()))
                        if vobject.representations["position"] is not None:
                            # A representation already existed: recreate it from scratch.
                            #vobject.representations["restraints"].define_new_indexes_to_vbo(indexes)
                            #vobject.representations["position"] = PickingDotsRepresentation(
                            vobject.representations["position"] = OneColorDotsRepresentation(
                                vobject, self.vm_session.vm_glcore,
                                active=True, indexes=pos_indexes,rgb=(1, 0, 1))
                        else:
                            # The key exists but was None: create it.
                            #vobject.representations["position"] = PickingDotsRepresentation(
                            vobject.representations["position"] = OneColorDotsRepresentation(
                                vobject, self.vm_session.vm_glcore,
                                active=True, indexes=pos_indexes,rgb=(1, 0, 1))
                    else:
                        # The 'restraints' key did not even exist yet: create it.
                        #vobject.representations["position"] = PickingDotsRepresentation(
                        vobject.representations["position"] = OneColorDotsRepresentation(
                            vobject, self.vm_session.vm_glcore,
                            active=True, indexes=pos_indexes,rgb=(1, 0, 1))

        # Request a redraw of the scene to reflect the changes.
        self.main.vm_session.vm_glcore.queue_draw()
        #print(indexes)

        #except:
        #    # Any failure here is treated as a REPRESENTATION-only error: the
        #    # message makes clear that the potentials/forces keep working normally
        #    # in the calculation; only the drawing failed.
        #    # (Note: a bare 'except' also swallows unexpected errors and makes
        #    #  debugging harder; catching specific exceptions would be safer.)
        #    print('\n\n Failed when trying to represent harmonic potential constraints. This is just a representation error, the potencies are working normally.\n\n')

    def add_new_harmonic_restraint (self, parameters, _type = 'distance'):
        """Register (in EasyHybrid) a new harmonic distance restraint.

        This function does NOT touch the pDynamo calculation; it only adds an
        entry to the system's e_restraints_dict so the restraint is tracked by
        the interface and drawn by update_restaint_representation. Converting it
        into the pDynamo model happens in other methods.

        Args:
            parameters: dict with 'system', 'atom1', 'atom2', 'distance',
                        'force_constant'.
            _type: only 'distance' is implemented for now.
        """
        #restraints = RestraintModel()
        #parameters['system'].DefineRestraintModel( restraints )

        if _type == 'distance':

            atom1 = parameters['atom1']
            atom2 = parameters['atom2']

            # Use the system's counter as a unique name for this restraint.
            rest_name = str(parameters['system'].e_restraint_counter)

            # Entry layout (same order read in update_restaint_representation):
            #   [ active?, name, type, [atom1, atom2], distance, k, e_id ]
            parameters['system'].e_restraints_dict[rest_name] = [
                True,                                            # active by default
                rest_name,                                       # name
                'distance',                                      # type
                [parameters['atom1'], parameters['atom2']],      # atom pair
                parameters['distance'],                          # target distance
                parameters['force_constant'],                    # force constant
                parameters['system'].e_id]                       # system id

            # Bump the counter so the next restraint gets a unique name.
            parameters['system'].e_restraint_counter += 1
        
        elif _type == 'position':
            # Use the system's counter as a unique name for this restraint.
            rest_name = str(parameters['system'].e_restraint_counter)

            # Cor dos "dots" desta restricao especifica (R,G,B). Opcional:
            # se nao for passada em parameters['color'], cai no magenta
            # (1, 0, 1) que era o valor fixo usado antes desta mudanca.
            color = parameters.get('color', (1, 0, 1))

            # Entry layout (same order read in update_restaint_representation):
            #   [ active?, name, type, [atomlist], reference, k, e_id, color ]
            parameters['system'].e_restraints_dict[rest_name] = [
                True,                                            # active by default
                rest_name,                                       # name
                'position',                                      # type
                parameters['atomlist']      ,                    # atom list
                parameters['reference']     ,                    # target reference: ( system.coordinates3 ) Dynamic?
                parameters['force_constant'],                    # force constant
                parameters['system'].e_id   ,                    # system id
                color]                                            # cor dos dots (R,G,B), por restricao

            # Bump the counter so the next restraint gets a unique name.
            parameters['system'].e_restraint_counter += 1
        
        
        else:
            # Types other than 'distance' are not supported yet.
            pass

    def define_harmonic_restraints (self, parameters):
        """Define a harmonic 'tether' (multiple positional restraint) in pDynamo.

        Builds a RestraintModel with a RestraintMultipleTether, which pins a set
        of atoms to reference positions through a harmonic potential, and binds
        it to the system via DefineRestraintModel.

        WARNING - as written, this function does NOT run: it references several
        names that are never defined (see the BUGs below). The 'parameters[...]'
        at the top are all overwritten with neutral/None values and then ignored,
        which suggests draft/incomplete code.
        """
        # These four fields are prepared but never used afterwards
        # (the names actually used are free variables, not 'parameters[...]').
        parameters['equilibriumValue'] = 0.0
        parameters['forceConstant']    = 0.0
        parameters['reference']        = None  # <- Clone ( system.coordinates3 )
        parameters['period']           = None
        parameters['selection']        = None

        # BUG: 'forceConstant' (a local variable) does not exist; it was only set
        #      as parameters['forceConstant']. This should be
        #      parameters['forceConstant'].
        tetherEnergyModel  = RestraintEnergyModel.Harmonic ( 0.0, forceConstant )

        tethers            = RestraintModel ( )

        # BUG: 'reference' and 'heavies' do not exist in this scope. They were
        #      probably meant to be parameters['reference'] and a selection of
        #      heavy atoms (heavies) that was never created here.
        tethers["Tethers"] = RestraintMultipleTether.WithOptions (
            energyModel = tetherEnergyModel,
            reference   = reference,
            selection   = heavies )

        # BUG: 'system' does not exist; it should be parameters['system'].
        system.DefineRestraintModel ( tethers )

    def define_distance_harmonic_restraints (self, parameters):
        """Define a harmonic DISTANCE restraint between 2 atoms in pDynamo.

        Creates a RestraintModel, computes the target distance for the current
        scan step and adds a harmonic RestraintDistance under the key 'RC1'
        (reaction coordinate 1).

        WARNING - it also does not run as written (see BUGs).
        """
        restraints = RestraintModel()
        #print('\n\n\n\n',parameters, '\n\n\n\n')

        # Bind a (still empty) restraint model to the system.
        parameters['system'].DefineRestraintModel( restraints )

        # Target distance for the step: minimum + increment * step index.
        # BUG: 'i' is not defined in this function. In scans, 'i' is usually the
        #      current step index and would need to be received in parameters
        #      (e.g. parameters['step']) or passed as an argument.
        distance = parameters['dminimum_RC1'] + ( parameters['dincre_RC1'] * float(i) )

        # Harmonic energy model centred on the target distance, with the given
        # force constant.
        rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['force_constant_1'])

        # BUG: 'atom1' and 'atom2' do not exist here; they should come from
        #      parameters['atom1'] / parameters['atom2'].
        restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom1, point2= atom2)

        # Register the distance restraint under reaction coordinate 'RC1'.
        # (Note: 'restraints' was bound to the system BEFORE receiving this entry.
        #  Depending on how pDynamo handles the reference, it may be necessary to
        #  call DefineRestraintModel again after filling it.)
        restraints["RC1"] = restraint



    def update_restaint_representation(self, e_id=None):
        if not e_id:
            e_id = self.active_id

        indexes = []          # pares de átomos das restrições de distância ativas
        position_groups = {}  # nome_da_restrição -> (indices_de_átomos, cor)

        for name, restraint in self.psystem[e_id].e_restraints_dict.items():
            _bool = restraint[0]
            _type = restraint[2]

            if _type == 'distance':
                if _bool:
                    indexes.append(restraint[3][0])
                    indexes.append(restraint[3][1])

            elif _type == 'position':
                if _bool:
                    # restraint[7] = cor (R,G,B) definida por
                    # add_new_harmonic_restraint. Restrições antigas (criadas
                    # antes desse campo existir) têm só 7 elementos: cai no
                    # magenta (1,0,1), que era o valor fixo anterior.
                    color = restraint[7] if len(restraint) > 7 else (1, 0, 1)
                    position_groups[restraint[1]] = (restraint[3], color)

        for vobject in self.main.vm_session.vm_objects_dic.values():
            if vobject.e_id != e_id:
                continue

            # --- representação das restrições de distância ---
            if indexes:
                vobject.representations["restraints"] = DashedLinesRepresentation(
                    vobject, self.vm_session.vm_glcore,
                    active=True, indexes=indexes)
            else:
                vobject.representations["restraints"] = None

            # --- representações das restrições de posição ---
            # OneColorDotsRepresentation só desenha UMA cor por instância, então
            # cada restrição de posição ganha sua própria representação (chave
            # "position_<nome>"), ao invés de uma única "position" compartilhada.
            active_keys = set()
            for rest_name, (atom_indexes, color) in position_groups.items():
                key = "position_{}".format(rest_name)
                active_keys.add(key)
                vobject.representations[key] = OneColorDotsRepresentation(
                    vobject, self.vm_session.vm_glcore,
                    active=True, indexes=atom_indexes, rgb=color)

            # Limpa representações de restrições de posição que foram
            # desativadas/removidas desde a última atualização, pra não deixar
            # "dots" órfãos desenhados na tela.
            stale_keys = [
                key for key in vobject.representations.keys()
                if key.startswith("position_") and key not in active_keys
            ]
            for key in stale_keys:
                vobject.representations[key] = None

        self.main.vm_session.vm_glcore.queue_draw()



    def define_position_harmonic_restraints (self, parameters):
        """Define a harmonic positional 'tether' for a selection in pDynamo.

        Unlike the distance version, this pins each atom of a selection to a
        reference POSITION (cloned coordinates) with a harmonic potential. This
        version is CORRECT and uses parameters[...] consistently.

        Expects in 'parameters':
            'reference'      -> reference coordinates (will be cloned)
            'selection'      -> atoms to be pinned
            'system'         -> target pDynamo system
            'force_constant' -> force constant of the potential
        """
        # . parameters['reference']
        # . parameters['selection']
        # . parameters['system']
        # . parameters['force_constant']

        # . Harmonically restrain heavy atoms.
        tethers       = None
        forceConstant = parameters['force_constant']

        # Clone the reference coordinates so the originals are not modified.
        reference          = Clone ( parameters['system'].coordinates3 )  # Change it later

        # Harmonic potential with equilibrium at 0.0 (zero displacement relative
        # to the reference) and the given force constant.
        tetherEnergyModel  = RestraintEnergyModel.Harmonic ( 0.0, forceConstant )

        tethers            = RestraintModel ( )
        tethers["Tethers"] = RestraintMultipleTether.WithOptions (
            energyModel = tetherEnergyModel,
            reference   = reference,
            selection   = parameters['selection'])

        # Bind the tether to the system -> it now applies in the energy calculation.
        parameters['system'].DefineRestraintModel ( tethers )

    def clear_restraints (self):
        """Remove any restraint model from the system (in pDynamo).

        Passing None to DefineRestraintModel clears the active restraints in the
        calculation.

        BUG: 'parameters' does not exist in this scope (the function takes no
             argument other than self). To work, it should receive the system,
             for example:
                 def clear_restraints(self, parameters):
                     parameters['system'].DefineRestraintModel(None)
             or use the active system, e.g. self.psystem[self.active_id].
        """
        parameters['system'].DefineRestraintModel(None)


#class Restraints:
#    """ Class doc """
#    
#    def __init__ (self):
#        """ Class initialiser """
#    
#    def update_restaint_representation (self, e_id = None):
#        """ Function doc """
#        try:
#            if e_id:
#                pass
#            else:
#                e_id = self.active_id
#            
#            indexes = []
#            for name, restraint in self.psystem[e_id].e_restraints_dict.items():
#                _bool = restraint[0]
#                name  = restraint[1]
#                _type = restraint[2]
#                if _type == 'distance':
#                    atons = '{} / {}'.format(restraint[3][0],restraint[3][1]) 
#                    dist_or_angle = '{:.4f}'.format(restraint[4])
#                    force_const   = str(restraint[5])
#                    e_id          =  restraint[6] 
#                    if _bool:
#                        indexes.append(restraint[3][0])
#                        indexes.append(restraint[3][1])
#
#
#            for vobject in self.main.vm_session.vm_objects_dic.values():
#                
#                if vobject.e_id == e_id:
#                    
#                    if indexes == []:
#                        
#                        vobject.representations["restraints"] = None
#                    
#                    else:
#                        
#                        if 'restraints' in vobject.representations.keys():
#                            
#                            if vobject.representations["restraints"] is not None:
#                                #print('["restraints"] is not None',indexes)
#                                #vobject.representations["restraints"].define_new_indexes_to_vbo(indexes)
#                                vobject.representations["restraints"] = DashedLinesRepresentation(vobject, self.vm_session.vm_glcore,
#                                                                                          active=True, indexes = indexes)
#                            else:
#                                #print('else',indexes)
#                                vobject.representations["restraints"] = DashedLinesRepresentation(vobject, self.vm_session.vm_glcore,
#                                                                                          active=True, indexes = indexes)
#                        else:
#                            vobject.representations["restraints"] = DashedLinesRepresentation(vobject, self.vm_session.vm_glcore,
#                                                                      active=True, indexes = indexes)    
#            self.main.vm_session.vm_glcore.queue_draw()
#            #print(indexes)
#    
#        except:
#            print('\n\n Failed when trying to represent harmonic potential constraints. This is just a representation error, the potencies are working normally.\n\n')
#        
#    def add_new_harmonic_restraint (self, parameters, _type = 'distance'):
#        """ Function doc """
#        #restraints = RestraintModel()
#        #parameters['system'].DefineRestraintModel( restraints )
#        
#        if _type == 'distance':
#            
#            atom1 = parameters['atom1']
#            atom2 = parameters['atom2']
#           
#            rest_name = str(parameters['system'].e_restraint_counter)
#            parameters['system'].e_restraints_dict[rest_name] = [True, 
#                                                                 rest_name ,
#                                                                 'distance', 
#                                                                 [parameters['atom1'],parameters['atom2']], 
#                                                                 parameters['distance'],  
#                                                                 parameters['force_constant'],
#                                                                 parameters['system'].e_id] 
#            
#            
#            parameters['system'].e_restraint_counter += 1
#        
#        else:
#            pass
#
#    def define_harmonic_restraints (self, parameters):
#        """ Function doc """
#        parameters['equilibriumValue'] = 0.0
#        parameters['forceConstant']    = 0.0
#        parameters['reference']        = None # <- Clone ( system.coordinates3 )
#        parameters['period']           = None
#        parameters['selection']        = None
#                
#        tetherEnergyModel  = RestraintEnergyModel.Harmonic ( 0.0, forceConstant )
#        tethers            = RestraintModel ( )        
#        tethers["Tethers"] = RestraintMultipleTether.WithOptions ( energyModel = tetherEnergyModel ,
#                                                                   reference   = reference         , 
#                                                                   selection   = heavies           )
#        system.DefineRestraintModel ( tethers )
#
#    def define_distance_harmonic_restraints (self, parameters):
#        """ Function doc """
#        restraints = RestraintModel()
#        #print('\n\n\n\n',parameters, '\n\n\n\n')
#        parameters['system'].DefineRestraintModel( restraints )
#        
#        distance = parameters['dminimum_RC1'] + ( parameters['dincre_RC1'] * float(i) )
#        
#        rmodel            = RestraintEnergyModel.Harmonic(distance, parameters['force_constant_1'])
#        restraint         = RestraintDistance.WithOptions(energyModel = rmodel, point1= atom1, point2= atom2)
#        restraints["RC1"] = restraint
#
#    def define_position_harmonic_restraints (self, parameters):
#        """ Function doc """
#        # . parameters['reference']
#        # . parameters['selection']
#        # . parameters['system']
#        # . parameters['force_constant']
#
#        # . Harmonically restrain heavy atoms.
#        tethers       = None
#        forceConstant = parameters['force_constant']
#        
#        reference          = Clone ( parameters['reference'] ) # Change it later
#        tetherEnergyModel  = RestraintEnergyModel.Harmonic ( 0.0, forceConstant )
#        tethers            = RestraintModel ( )
#        tethers["Tethers"] = RestraintMultipleTether.WithOptions ( energyModel = tetherEnergyModel ,
#                                                                   reference   = reference         , 
#                                                                   selection   = parameters['selection'])
#        parameters['system'].DefineRestraintModel ( tethers )
#
#
#
#    def clear_restraints (self):
#        """ Function doc """
#        parameters['system'].DefineRestraintModel(None)
#
#
#
