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
from vismol.libgl.representations import SurfaceRepresentation

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

    def register_change_and_maybe_autosave (self):
        """ Ponto de entrada do criterio de autosave por CONTADOR DE EVENTOS
            (o criterio por TEMPO/timer periodico vive em main_window.py,
            MainWindow._on_autosave_timer_tick -- os dois disparam o mesmo
            _do_autosave() abaixo, o que vier primeiro).

            Chame este metodo (em vez de save_easyhybrid_session(tmp=True)
            direto) em qualquer acao que deva contar como "mudanca" para
            fins de autosave -- ver as chamadas ja existentes em
            treeview_menu.py/treeview_menu_new.py (_save_backup_file).

            Sempre marca self.changed = True (a mudanca aconteceu de
            qualquer forma, autosave ligado ou nao -- isso e' o que faz
            on_delete_event perguntar antes de fechar). O autosave de fato
            so' dispara se gl_parameters['autosave'] estiver ligado E o
            contador atingir gl_parameters['autosave_event_count'].
        """
        self.changed = True
        gl_parameters = self.vm_session.vm_config.gl_parameters
        if not gl_parameters.get('autosave', True):
            return
        self.autosave_change_counter += 1
        threshold = gl_parameters.get('autosave_event_count', 20)
        if self.autosave_change_counter >= threshold:
            self._do_autosave()

    def _do_autosave (self):
        """ Dispara o autosave de fato: grava <arquivo_da_sessao>~ (ou um
            arquivo temporario em vm_config.easyhybrid_tmp se a sessao ainda
            nao foi salva nenhuma vez -- ver save_easyhybrid_session(tmp=
            True) abaixo, que ja' cobria esse caso). Zera o contador de
            eventos. Chamado tanto pelo criterio de contador (acima) quanto
            pelo timer periodico (main_window.py). Nunca propaga excecao --
            uma falha de autosave nao deve interromper o uso normal do
            programa. """
        try:
            self.save_easyhybrid_session(filename=self.main.session_filename, tmp=True)
        except Exception as e:
            logger.warning('Autosave failed: %s', e)
        self.autosave_change_counter = 0

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
            # [NEW - OPTIONAL property] Size of the GLArea (viewport) at
            # save time. Used on load to restore the same aspect
            # ratio and avoid camera distortion (the projection_matrix saved
            # above was computed for THIS aspect ratio; if the .easy is
            # reopened in a window of different size/proportion, without this
            # the image would be stretched/compressed). NEW key: files
            # .easy saved before this change simply will not have it --
            # load_easyhybrid_serialization_file uses .get() and only tries
            # to restore the size if the key exists, so old sessions
            # keep loading normally, just without this extra adjustment.
            'glarea_width'        : float(self.vm_session.vm_glcore.width),
            'glarea_height'       : float(self.vm_session.vm_glcore.height),
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
                            # [NOVO] Dados extras necessarios pra reconstruir a
                            # surface on load (see load_easyhybrid_serialization_
                            # file / _rebuild_surface_vobject_from_saved_data
                            # below). These are all NEW keys within vobj_data
                            # -- .easy sessions saved before this change only have
                            # 'is_surface' (without the keys below), so the loader
                            # keeps handling that case as before (surface
                            # skipped, rest of the session loads normally).
                            #
                            # surface_trajectory: the mesh (vertices/normals/faces
                            # per frame) already computed, the same thing the
                            # surface uses to draw in memory -- save it
                            # directly instead of trying to recompute on load (avoids
                            # depending on the QC system/.cube file/etc. still
                            # being available and identical).
                            vobj_data['surface_trajectory'] = getattr(vobject, 'surface_trajectory', None)
                            vobj_data['surface_type']       = getattr(vobject, 'surface_type', None)
                            vobj_data['parameters']         = getattr(vobject, 'parameters', None)
                            vobj_data['model_mat']          = np.copy(vobject.model_mat)
                            vobj_data['trans_mat']          = np.copy(vobject.trans_mat)
                            # MEP-specific setup (colormap + color
                            # limits), read/written by the "Surface Setup" window
                            # (see treeview_menu.py). Only exists for
                            # surface_type == "mep"; getattr with default None
                            # covers the other types without error.
                            vobj_data['mep_cmap_name'] = getattr(vobject, 'mep_cmap_name', None)
                            vobj_data['mep_vmin']      = getattr(vobject, 'mep_vmin', None)
                            vobj_data['mep_vmax']      = getattr(vobject, 'mep_vmax', None)
                            # One or two SurfaceRepresentation per surface
                            # object (e.g. "surface1"/"surface2" -- positive/
                            # negative lobe of an orbital, or just "surface1"
                            # for simple density/potential/MEP). Only needs
                            # to save 'surf_name' (the key used in
                            # vismol_object.surface_trajectory[frame][surf_name]
                            # to fetch vertices/colors/normals/indices -- see
                            # SurfaceRepresentation.draw_representation) and
                            # 'active'. The constructor's 'iso_color' parameter
                            # is NOT saved here because the class never actually
                            # uses it (accepted in __init__ but never becomes an attribute
                            # nor referenced in the rest of the class body --
                            # the real color is embedded per vertex within
                            # surface_trajectory itself, already saved above).
                            #
                            # 'render_mode' ("surface"/"lines" -- wireframe),
                            # 'alpha' (opacity) and 'smooth_shading' (flat vs.
                            # per-vertex normal) are the common controls of the
                            # "Surface Setup" window (treeview_menu.py) --
                            # saved here to restore exactly as the
                            # user left it, instead of always reverting to the
                            # default (opaque, filled, flat shading).
                            vobj_data['surface_representations'] = {}
                            for rep_name, rep in vobject.representations.items():
                                if isinstance(rep, SurfaceRepresentation):
                                    vobj_data['surface_representations'][rep_name] = {
                                        'surf_name'     : getattr(rep, 'surf_name', rep_name),
                                        'active'        : rep.active,
                                        'render_mode'   : getattr(rep, 'render_mode', 'surface'),
                                        'alpha'         : getattr(rep, 'alpha', 1.0),
                                        'smooth_shading': getattr(rep, 'smooth_shading', False),
                                    }
                        
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
            
            # [NEW - OPTIONAL property] .easy files saved before this
            # change do not have 'glarea_width'/'glarea_height' -- .get() returns
            # None and the block below is skipped, without breaking anything (same
            # compatibility pattern as the 'camera' block above).
            #
            # When present: tries to resize the GLArea to the same size
            # as when it was saved (set_size_request -- 'best effort', the
            # GTK container may not respect it 100% depending on the layout) and,
            # more importantly to actually eliminate distortion, recomputes the
            # aspect ratio/projection_matrix via resize_window() using that
            # SAME saved size (with z_near/z_far already restored above) --
            # the raw projection_matrix restored just above was computed
            # for the aspect ratio AT SAVE TIME, which may not match
            # the current window's; resize_window ensures consistency.
            glarea_w = camera_data.get('glarea_width')
            glarea_h = camera_data.get('glarea_height')
            if glarea_w and glarea_h:
                vm_widget = getattr(self.vm_session, 'vm_widget', None)
                if vm_widget is not None:
                    try:
                        vm_widget.set_size_request(int(glarea_w), int(glarea_h))
                    except Exception as e:
                        dprint('Could not resize GLArea on session load:', e)
                vm_glcore.resize_window(glarea_w, glarea_h)
            
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
                        self._rebuild_surface_vobject_from_saved_data(system = system, vobj = vobj)
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
        
    def _rebuild_surface_vobject_from_saved_data (self, system, vobj):
        """ Reconstroi um VismolObject de superficie (orbital/densidade/
            potencial/MEP/cubo externo -- ver surface_analysis_window.py)
            a partir dos dados salvos em save_easyhybrid_session.

            [PROPRIEDADE OPCIONAL / COMPATIVEL] Sessoes .easy salvas antes
            desta mudanca tem 'is_surface' mas nao tem 'surface_trajectory'
            -- nesse caso, a superficie e' simplesmente pulada (mesmo
            comportamento de sempre), sem quebrar o carregamento do resto
            da sessao.

            Nao recalcula a malha (nao chama de volta o gerador de
            orbital/densidade/MEP/etc.) -- so' restaura a malha ja'
            computada e salva (vertices/cores/normais/indices por frame,
            em vobj['surface_trajectory']), o que evita depender do
            sistema QC, do arquivo .cube original ou de qualquer estado
            de calculo ainda estarem disponiveis/identicos no momento do
            load.
        """
        name = vobj.get('name', 'surface')
        if vobj.get('surface_trajectory') is None:
            dprint('Surface object "%s" has no saved geometry (older .easy '
                   'file, or was created before this feature) -- skipping.' % name)
            return
        
        vobject_tmp = VismolObject(name = name, index = -1,
                                   vismol_session        = self.vm_session,
                                   trajectory            = [],
                                   bonds_pair_of_indexes = [0, 1])
        
        vobject_tmp.model_mat = vobj.get('model_mat', np.identity(4, dtype=np.float32))
        vobject_tmp.trans_mat = vobj.get('trans_mat', np.identity(4, dtype=np.float32))
        vobject_tmp.surface_trajectory = vobj['surface_trajectory']
        vobject_tmp.parameters         = vobj.get('parameters')
        vobject_tmp.surface_type       = vobj.get('surface_type')
        vobject_tmp.frames             = vobj['frames']
        vobject_tmp.active             = vobj.get('active', True)
        vobject_tmp.is_surface         = True
        vobject_tmp.e_id               = system.e_id
        
        # Setup especifico de MEP (colormap + limites de cor). getattr/get
        # com default None: para outros surface_type, essas chaves existem
        # in vobj (saved as None) and simply do nothing here.
        if vobj.get('mep_cmap_name') is not None:
            vobject_tmp.mep_cmap_name = vobj['mep_cmap_name']
        if vobj.get('mep_vmin') is not None:
            vobject_tmp.mep_vmin = vobj['mep_vmin']
        if vobj.get('mep_vmax') is not None:
            vobject_tmp.mep_vmax = vobj['mep_vmax']
        
        # Recreates one SurfaceRepresentation per saved entry (normally
        # "surface1" and, for orbital/MEP, also "surface2" -- see
        # save_easyhybrid_session). Each one fetches its own mesh from
        # vobject_tmp.surface_trajectory[frame][surf_name] at draw
        # time, so we only need to recreate the object with the same
        # saved surf_name/active -- there is no mesh to pass here.
        surface_reps = vobj.get('surface_representations') or {}
        for rep_name, rep_data in surface_reps.items():
            rep = SurfaceRepresentation(
                vismol_object = vobject_tmp,
                vismol_glcore = self.vm_session.vm_glcore,
                name          = 'surface',
                active        = rep_data.get('active', True),
                indexes       = [],
                is_dynamic    = False,
                surface_name  = rep_data.get('surf_name', rep_name),
            )
            # Restores wireframe/opacity/shading exactly as the user
            # left it ("Surface Setup" window, treeview_menu.py) -- without this,
            # every reloaded surface would revert to the default (opaque,
            # filled, flat shading), losing any adjustment made.
            if hasattr(rep, 'set_render_mode'):
                rep.set_render_mode(rep_data.get('render_mode', 'surface'))
            if hasattr(rep, 'set_alpha'):
                rep.set_alpha(rep_data.get('alpha', 1.0))
            if hasattr(rep, 'set_shading_mode'):
                rep.set_shading_mode('smooth' if rep_data.get('smooth_shading', False) else 'flat')
            vobject_tmp.representations[rep_name] = rep
        
        self.vm_session._add_vismol_object(vobject_tmp, show_molecule = False, autocenter = False)
        # [KNOWN LIMITATION] When generating a surface for the first time
        # (surface_analysis_window.py), it is nested in the treeview under the
        # vobject that originated it (vobj_parent = <original vobject>.
        # e_treeview_iter). Here we only have 'system' (the pDynamo system) at
        # hand, not the specific parent vobject -- e_treeview_iter lives in the
        # VOBJECT, not the system (see main_treeview.py). Passing the wrong
        # iter would risk an AttributeError/incorrect behavior, so
        # for now the reconstructed surface enters as a TOP-level item in the
        # treeview (still functional/visible, just not visually nested
        # under the original molecule as the first time).
        self.main.main_treeview.add_vismol_object_to_treeview(vobject_tmp)
        self.main.add_vobject_to_vobject_liststore_dict(vobject_tmp)
        self.main.refresh_widgets()
        
    def save_special_PDB (vObject):
        """ Function doc """
