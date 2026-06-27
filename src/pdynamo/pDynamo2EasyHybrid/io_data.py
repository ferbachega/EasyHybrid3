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

class LoadAndSaveData:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass


    def save_easyhybrid_session (self, filename = 'session.easy', tmp = False):
        """   
        When the interface makes some modification to the session, 
        a temporary file "filename.easy~" is saved. When the session 
        is saved by the user, the interface checks if there is a 
        temporary file and deletes it.
        """
        easyhybrid_session_data = {}
        backup = {}
        
        '''- - - - - - - - - - pDynamo systems - - - - - - - - - - - '''
        #easyhybrid_session_data['psystem'] = self.psystem
        '''- - - - - - - - - - - - - - - - - - - - - - - - - - - - - '''
        self.main.bottom_notebook.get_active_system_text_from_textbuffer()
        easyhybrid_session_data['systems'] = [ ]
        
        '''- - - - - - - - - - camera/view orientation - - - - - - - - - - '''
        # Captures the state that fully determines how the scene is framed:
        # model_mat (world rotation/pan - every VismolObject inherits this
        # at creation time via set_model_matrix), the camera's view_matrix
        # (zoom/position), its near/far clipping planes (change with zoom
        # too), and zero_reference_point/dist_cam_zrp (used for pan, fog
        # and line width calculations). Restoring just these on load is
        # enough: objects loaded afterwards pick up vm_glcore.model_mat
        # automatically, they don't need their own per-object copy saved.
        glcamera = self.vm_session.vm_glcore.glcamera
        easyhybrid_session_data['camera'] = {
            'model_mat'           : np.copy(self.vm_session.vm_glcore.model_mat),
            'zero_reference_point': np.copy(self.vm_session.vm_glcore.zero_reference_point),
            'dist_cam_zrp'        : float(self.vm_session.vm_glcore.dist_cam_zrp),
            'view_matrix'         : np.copy(glcamera.view_matrix),
            'projection_matrix'   : np.copy(glcamera.projection_matrix),
            'z_near'              : float(glcamera.z_near),
            'z_far'               : float(glcamera.z_far),
        }
        
        for e_id, system in self.psystem.items():
            
            if system == None:
                pass
            
            else:
                data   = {}
    
                #backup[e_id] = []
                #backup[e_id].append(system.e_treeview_iter)
                #backup[e_id].append(system.e_liststore_iter)
    
                #system.e_treeview_iter   = None
                #system.e_liststore_iter  = None
                
                data['system'] = system
                data['vobjects'] = []
                for key, vobject in self.vm_session.vm_objects_dic.items():
                    if system.e_id == vobject.e_id:
                        #data['frames'] = vobject.frames
                        #data['color_palette'] = vobject.color_palette
                        vobj_data = {}
                        vobj_data['frames']        = vobject.frames
                        vobj_data['color_palette'] = vobject.color_palette
                        vobj_data['name']          = vobject.name
                        vobj_data['active']        = vobject.active
                        vobj_data['key6']          = vobject.key6
                        vobj_data['cell_coordinates'] = vobject.cell_coordinates
                        
                        if key in system.e_logfile_data.keys():
                            vobj_data['logfile_data'] = system.e_logfile_data[key]
                        
                        if getattr (vobject, 'idx_2D_xy', False):
                            vobj_data['idx_2D_xy'] = vobject.idx_2D_xy
                        
                        if getattr (vobject, 'is_surface', False):
                            vobj_data['is_surface'] = vobject.is_surface
                        
                        data['vobjects'].append(vobj_data)
                            
                easyhybrid_session_data['systems'].append(data)
        
        #---------------------------------------------------------------
        """   
        When the interface makes some modification to the session, 
        a temporary file "filename.easy~" is saved. When the session 
        is saved by the user, the interface checks if there is a 
        temporary file and deletes it.
        """
        if tmp:
            
            #.Saves a temporary file with the same name as the original session + "~"
            if filename:
                tmpfile = filename+'~'
                with open(tmpfile,'wb') as outfile:
                    pickle.dump(easyhybrid_session_data, outfile)
            
            #.Saves a temporary file when there is no previous session name.
            else:
                tmpfile = os.path.join(self.vm_session.vm_config.easyhybrid_tmp, 
                                       self.random_code+'.easy')
                with open(tmpfile,'wb') as outfile:
                    pickle.dump(easyhybrid_session_data, outfile)
                


        else:
            with open(filename,'wb') as outfile:
                pickle.dump(easyhybrid_session_data, outfile)
                self.changed = False
            if os.path.exists(filename+'~'):
                os.remove(filename+'~')
            else:
                pass
        #---------------------------------------------------------------
            
        for e_id, data in backup.items():
            pass
            #self.psystem[e_id]
            #self.psystem[e_id].e_treeview_iter   = data[0]
            #self.psystem[e_id].e_liststore_iter  = data[1]
        
        self.main.session_filename = filename
        
        if tmp:
            pass
        else:
            self.main.bottom_notebook.status_teeview_add_new_item(message = ':  {}  saved'.format(filename), 
                                                               system =  system )
        
        #'''- - - - - - - - - - - - vismol obejcts - - - - - - - - - - - '''
        #vobjects = {}
        #for key, vobject in self.vm_session.vm_objects_dic.items():
        #    parameters = {
        #                  'index'             : vobject.index            ,
        #                  'name'              : vobject.name             ,
        #                  'active'            : vobject.active           ,
        #                  'frames'            : vobject.frames           ,
        #                  'color_palette'     : vobject.color_palette    ,
        #                  'mass_center'       : vobject.mass_center      ,
        #                  'selected_atom_ids' : vobject.selected_atom_ids,
        #                  'index_bonds'       : vobject.index_bonds      ,
        #                                        
        #                  'colors'            : vobject.colors           ,
        #                  'color_indexes'     : vobject.color_indexes    ,
        #                 }
        #    vobjects[key] = parameters
        #easyhybrid_session_data['vobjects'] = vobjects
        #'''- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '''

    def load_easyhybrid_serialization_file (self, filename = None, tmp = False):
        self.main.restart() 
        self.main.bottom_notebook.status_teeview_add_new_item(message = ':  {}  loaded'.format(filename), 
                                                               system =  None )
        if filename is None:
            return None
        with open(filename, "rb") as f:
            # Load the object from the file
            easyhybrid_session_data = pickle.load(f)
        #print(easyhybrid_session_data)
        
        '''- - - - - - - - - - camera/view orientation - - - - - - - - - - '''
        # Older .easy files won't have this key - skip restoring and keep
        # whatever default view vm_glcore already has in that case.
        camera_data = easyhybrid_session_data.get('camera')
        if camera_data is not None:
            vm_glcore = self.vm_session.vm_glcore
            glcamera = vm_glcore.glcamera
            vm_glcore.model_mat            = np.copy(camera_data['model_mat'])
            vm_glcore.zero_reference_point = np.copy(camera_data['zero_reference_point'])
            vm_glcore.dist_cam_zrp         = camera_data['dist_cam_zrp']
            glcamera.set_view_matrix(np.copy(camera_data['view_matrix']))
            glcamera.set_projection_matrix(np.copy(camera_data['projection_matrix']))
            glcamera.z_near = camera_data['z_near']
            glcamera.z_far  = camera_data['z_far']
            glcamera.update_fog()
            vm_glcore.queue_draw()
        
        for data  in easyhybrid_session_data['systems']:
            system = data['system']
            
            #.checking e_job_history attribute
            if hasattr(system, 'e_job_history'):
                pass
            else:
                system.e_job_history = {}
            
            name   = system.label
            tag    = system.e_tag
            #print('\n\n\n\n',system, name, tag, data['system'], data['vobjects'] )
            if len(data['vobjects']) == 0:
                pass
            else:
                self.add_new_system_to_psession (system = system, name  = name, tag = tag)
                self.main.main_treeview.add_new_system_to_treeview (system)
                ff  =  getattr(system.mmModel, 'forceField', "None")
                self.main.bottom_notebook.status_teeview_add_new_item(message = 'New System:  {} ({}) - Force Field:  {}'.format(system.label, system.e_tag, ff), system = system)
                
                for vobj  in data['vobjects']:
                    frames = vobj['frames']
                    name   = vobj['name']
                    
                    #if 'key6' in vobj.keys():
                    
                    
                    if 'is_surface' in vobj.keys():
                        pass
                    else:
                        vobj['is_surface'] = False
                    
                    if vobj['is_surface']:
                        print(vobj['is_surface'])
                        pass
                    else:
                        vm_object = self._build_vobject_from_pdynamo_system ( system = system, name = name ) 
                        vm_object.frames = frames
                        vm_object.active = vobj['active']
                        
                        if 'key6' in vobj.keys():
                            vm_object.key6 = vobj['key6']
                         
                        
                        self.vm_session._add_vismol_object(vm_object, show_molecule = True,
                                                            autocenter = (camera_data is None))
                        
                        self.main.main_treeview.add_vismol_object_to_treeview(vm_object)
                        
                        self.main.add_vobject_to_vobject_liststore_dict(vm_object)
                        
                        self._apply_fixed_representation_to_vobject(vismol_object =vm_object)
                        self._apply_QC_representation_to_vobject(vismol_object =vm_object)
                        
                        self.main.refresh_widgets()
                        
                        if 'logfile_data' in vobj.keys():
                            system.e_logfile_data[vm_object.index] = vobj['logfile_data']
                        if 'idx_2D_xy' in vobj.keys():
                            vm_object.idx_2D_xy  = vobj['idx_2D_xy']
                            
                        if 'cell_coordinates' in vobj.keys():
                            vm_object.cell_coordinates  = vobj['cell_coordinates']
        if tmp:
            filename = filename.replace('~', '')
            self.main.session_filename = filename
        else:
            self.main.session_filename = filename
        self.main.process_manager_window.build_liststore_from_job_history (clear = True)
        
    
    def save_special_PDB (vObject):
        """ Function doc """
