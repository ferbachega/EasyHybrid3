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
import os
import json

class VismolConfig                       :
    """ Class doc """
    
    def __init__ (self, vismol_session = None, home = None):
        """ Class initialiser """
        self.easyhybrid_home = home
        self.easyhybrid_tmp  = None
        

        
        
        #.EasyHybrid Default Parameters / Only used when the .config.json is not found 
        self.gl_parameters_default = {"background_color"           : [0.0, 0.0, 0.0, 1.0],#[1.0, 1.0, 1.0, 1.0],#"background_color"                       : [0.0, 0.0, 0.0, 1.0],
                                      "color_type"                 : 0,
                                      "dot_size"                   : 2,
                                      "dots_size"                  : 2,
                                      "dot_type"                   : 1,
                                      "dot_sel_size"               : 2.0,
                                      "line_width"                 : 3,
                                      "line_width_selection"       : 10,
                                      "line_type"                  : 1,
                                      "line_color"                 : 0,
                                      
                                      "ribbon_width"               : 0.4, # Now being used (defined in the shader)
                                      #"ribbon_width_selection"     : 0.4, # Now being used (defined in the shader)
                                      "ribbon_type"                : 2,
                                      "ribbon_color"               : 0,
                                      "sphere_type"                : 0,
                                      #"sphere_scale"              : 0.20,
                                                                
                                      "sphere_scale"               : 0.20,
                                                                
                                      "sphere_quality"             : 2,
                                      "impostor_type"              : 0,
                                      #"sticks_radius"             : 2.5,
                                      
                                      "sticks_radius"              : 0.16, # not being used (defined in the shader)
                                      "multiple_bonds"             : False,
                                      "sticks_color"               : 0,
                                      "sticks_type"                : 0,
                                      "antialias"                  : True,
                                      "mouse_rotation_sensibility" : 1.5,
                                      "scroll_step"                : 0.9,
                                      "field_of_view"              : 10,
                                      "light_position"             : [0, 0, 10.0],
                                      #"light_position"            : [-2.5, 2.5, 3.0],
                                      "light_color"                : [ 1.0, 1.0, 1.0, 1.0],
                                      "light_ambient_coef"         : 0.4,
                                      "light_shininess"            : 5.5,
                                      "light_intensity"            : [0.6, 0.6, 0.6],
                                      "light_specular_color"       : [1.0, 1.0, 1.0],
                                      "center_on_coord_sleep_time" : 0.01,
                                      "gridsize"                   : 0.8,
                                      "maxbond"                    : 2.6,
                                      "bond_tolerance"             : 1.4,
                                      
                                      "picking_dots_color"         : [0.0, 1.0, 1.0],
                                      "picking_dots_safe"          : True,
                                      "pk_label_color"             : [1.0, 1.0, 1.0, 1.0],
                                      "pk_dist_label_color"        : [1.0, 1.0, 1.0, 1.0],
                                      "dashed_dist_lines_color"    : [0.1, 0.1, 0.1, 1.0],

                                      # [NOVO] Dynamic bonds (regiao QC): opcao de
                                      # draw the dynamic bonds in a single
                                      # color, instead of per-atom color. Off
                                      # by default (keeps the current representation).
                                      # The default color is white.
                                      "dynamic_bonds_single_color" : False,
                                      "dynamic_bonds_color"        : [1.0, 1.0, 1.0, 1.0],
                                      
                                      #"pk_label_size"              : [1.0, 1.0, 1.0, 1.0],
                                      
                                      'startup_path'               : None,
                                      'tmp_files'                  : None,
                                      'autosave'     : True, 
                                      # Criterio de autosave: dispara pelo que
                                      # vier primeiro -- a cada N minutos
                                      # (timer, ver MainWindow._restart_
                                      # autosave_timer) OU a cada N mudancas
                                      # registradas (contador, ver
                                      # pDynamoSession.register_change_and_
                                      # maybe_autosave). Grava em
                                      # <session_file>~ (or a temporary
                                      # file if the session has not been
                                      # salva nenhuma vez).
                                      'autosave_interval_minutes' : 5,
                                      'autosave_event_count'      : 20,
                                      'askSaveUnsave': True, 
                                      #'askSaveUnsave': True, 
                                      # Main window dimensions, saved
                                      # ao fechar o EasyHybrid e restauradas
                                      # na proxima abertura (ver MainWindow.
                                      # __init__ / window_resize / on_delete_
                                      # event em main_window.py).
                                      'main_window_width'  : 1200,
                                      'main_window_height' : 600,
                                      # Posicao do divisor paned_V (area 3D/
                                      # treeview vs. notebook de baixo -- Status/
                                      # Annotations/Sequence), como PROPORCAO da
                                      # window height (not fixed pixels), to
                                      # follow the window size when it
                                      # muda. Default 400/600 (posicao antiga
                                      # fixa / altura default). Ver main_window.py:
                                      # window_resize / on_paned_v_position_changed.
                                      'main_window_paned_v_ratio' : 400.0/600.0,
                                      # Liga/desliga salvar/restaurar as dimensoes
                                      # of the main window between sessions (checkbox
                                      # ja existia no glade -- "Save window size" --
                                      # but the value was never used). See MainWindow.
                                      # __init__ / window_resize / on_delete_event.
                                      'save_window_size'   : True,
                                      # V-Sync ("vblank_mode" env var, read by the
                                      # Mesa/GLX driver -- see easyhybrid.py's
                                      # _maybe_disable_vsync_for_intel_igpu, which
                                      # must set it BEFORE the GL context is
                                      # created, i.e. before GTK is even imported).
                                      # 'auto' keeps the existing behavior (disable
                                      # V-Sync automatically on Intel iGPUs only,
                                      # to avoid a rotation/pan/zoom stutter);
                                      # 'on'/'off' force it either way regardless
                                      # of GPU vendor. Only takes effect after a
                                      # restart -- see Preferences > Startup.
                                      'vblank_mode'         : 'auto',
                                      }
                              
        self.n_proc = 2
        # self.representations_available = {"dots", "lines", "nonbonded", "dotted_lines",
        #                                   "ribbon", "sticks", "spheres", "impostor",
        #                                   "surface", "cartoon", "freetype",
        #                                   "picking_dots"}
        
        #.Rep list - Don't change this list
        # [EN] "cartoon" re-added (was only in the commented-out version
        # above). THIS is the file actually used at runtime -- easyhybrid.py
        # (the real entry point) does "from gui.config import VismolConfig"
        # and passes an instance of THIS class explicitly to
        # EasyHybridSession(vm_config=vconfig), which overrides the
        # submodule's own default vismol/core/vismol_config.py entirely
        # (see VismolSession.__init__: "if vm_config: self.vm_config =
        # vm_config"). An earlier fix mistakenly edited the submodule's
        # vismol_config.py instead -- same class name (VismolConfig),
        # same near-identical content, different file, never actually
        # instantiated by the real app. Without "cartoon" here,
        # vm_glcore.initialize()'s shader-compile loop
        # (for rep in self.vm_config.representations_available: ...)
        # simply never attempts "cartoon" at all -- no error at startup
        # (nothing failed, it just never got tried), and
        # shader_programs["cartoon"] stays unset until the first
        # attempt to actually draw a Cartoon representation, which is
        # exactly where the KeyError shows up.
        self.representations_available = {"dots", "lines", "nonbonded", "impostor",'dash', "posdot_type",
                                          "sticks", "spheres", 'ribbons', 'cartoon', #'ribbon_sphere', 
                                          'dynamic','vdw_spheres', 
                                          'picking_spheres','static_freetype', 'surface'}
    
        
        #.Checking temporary folder
        self._check_tmp_folder()
        
        self.vismol_session = vismol_session
        
        #.Checking configuration file.
        self._check_config_file()
        self._check_startup_path()
    
    def _check_startup_path (self):
        """ Function doc """
        
        if 'startup_path' in self.gl_parameters:
            
            if self.gl_parameters['startup_path']:                
                #. checking startup_path is trully a folder
                if os.path.isdir(self.gl_parameters['startup_path']):
                    pass
                else:
                    self.gl_parameters['startup_path'] = self.easyhybrid_home
            
            else:
                self.gl_parameters['startup_path'] = self.easyhybrid_home
            
            
            
        #else:
        #    if os.path.isdir( os.path.join (self.easyhybrid_home,'PROJECTS')):
        #        self.gl_parameters['startup_path'] = os.path.join (self.easyhybrid_home,'PROJECTS')
        #        #self.easyhybrid_tmp = os.path.join ( self.easyhybrid_home,'.tmp')
        #    
        #    else:
        #        #print('')
        #        projects = os.path.join (self.easyhybrid_home,'PROJECTS')
        #        os.mkdir(projects)
        #        self.gl_parameters['startup_path'] = projects
        #        #self.easyhybrid_tmp = os.path.join ( self.easyhybrid_home,'.tmp')
        #    
        #    
        #    #self.gl_parameters['startup_path'] = self.easyhybrid_home
        

    def _check_config_file (self):
        """ 
        This function checks if the config/json file exists. 
        If not, a new config/json  is created. 
        """
        if self.easyhybrid_home is not None:

            config = os.path.join ( self.easyhybrid_home,'.config.json')
            
            if os.path.exists(config):
                
                try:
                    dprint("Loading EasyHybrid config file.")
                    with open(config, 'r', encoding='utf-8') as f:
                        self.gl_parameters = json.load(f)
                
                        for key in self.gl_parameters_default.keys():
                            if key in self.gl_parameters.keys():
                                pass
                            else:
                                # [BUG FIX] Era 'self.gl_parameters_default[keys]'
                                # -- 'keys' does not exist (typo of 'key', the variable
                                # do loop). Isso disparava NameError toda vez que
                                # uma chave nova em gl_parameters_default (ex.:
                                # as que acabamos de adicionar: autosave_interval_
                                # minutes, main_window_width, etc.) did not exist
                                # ainda no .config.json salvo -- capturado pelo
                                # 'except' abaixo, que resetava TODAS as
                                # saved preferences (not just the missing key)
                                # de volta pro default, silenciosamente.
                                self.gl_parameters[key] = self.gl_parameters_default[key]
                                
                except:
                    dprint("Failed to open EasyHybrid configuration file. Loading default settings.")
                    self.gl_parameters = self.gl_parameters_default
                
            else:
                dprint("Configuration file not found. Creating a new file from default settings.")
                with open(config, 'w', encoding='utf-8') as f:
                    json.dump(self.gl_parameters_default, f, ensure_ascii=False, indent=4)
                    self.gl_parameters = self.gl_parameters_default
    
    def _check_tmp_folder (self):
        """ 
        This function checks if the temporary files directory exists. 
        If not, a new directory is created. 
        """ 
        if self.easyhybrid_home is not None:

            if os.path.isdir( os.path.join ( self.easyhybrid_home,'.tmp')):
                self.easyhybrid_tmp = os.path.join ( self.easyhybrid_home,'.tmp')
            
            else:
                #print('')
                os.mkdir(os.path.join ( self.easyhybrid_home,'.tmp'))
                self.easyhybrid_tmp = os.path.join ( self.easyhybrid_home,'.tmp')
        #self.tmp_folder = self.easyhybrid_home,'.tmp'
    def reset_parameters (self):
        """ Function doc """
        try:
            self.gl_parameters = self.gl_parameters_default
            return True
        except:
            return False
        #Pickle( os.path.join ( folder, filename+'.pkl'), 
        #        system.coordinates3 )

    def save_easyhybrid_config(self)                       :
        """ Function doc """
        
        if self.easyhybrid_home is not None:
            config = os.path.join ( self.easyhybrid_home,'.config.json')
            
            dprint("Saving configuration file.")
            
            with open(config, 'w', encoding='utf-8') as f:
                json.dump(self.gl_parameters, f, ensure_ascii=False, indent=4)        


    
    def load_easyhybrid_config(self, config_path):
        """ Carrega preferencias salvas (.config.json) por cima dos defaults.

            [BUG FIX] Esta funcao existia mas nunca era chamada em lugar
            nenhum do codigo -- ou seja, o botao "Apply and Save Changes"
            da janela de Preferences escrevia o arquivo (save_easyhybrid_
            config), mas nada NUNCA lia esse arquivo de volta. Toda
            preferencia salva era perdida a cada reinicio do EasyHybrid.
            Ver chamada em easyhybrid.py, logo apos criar o VismolConfig.

            [BUG FIX] Antes fazia 'self.gl_parameters = json.load(...)',
            substituindo o dict inteiro. Isso e' perigoso: se uma versao
            nova do EasyHybrid adicionar uma chave nova em gl_parameters_
            default (ex.: 'multiple_bonds', 'main_window_width'), carregar
            um .config.json salvo por uma versao ANTIGA (sem essa chave)
            apagaria a chave nova inteira, e qualquer codigo que espera
            gl_parameters['chave_nova'] quebraria com KeyError. Agora faz
            .update(...) por cima do dict default: chaves salvas sobre-
            escrevem o default, chaves novas que nunca foram salvas
            mantem o valor default.
        """
        if not os.path.isfile(config_path):
            config_path = os.path.join(os.environ["HOME"], ".VisMol", "VismolConfig.json")
        if not os.path.isfile(config_path):
            # First run (or deleted file): keeps the defaults,
            # sem erro.
            return False
        try:
            with open(config_path, "r") as config_file:
                loaded = json.load(config_file)
        except Exception as e:
            dprint("Could not load saved preferences (%s): %s" % (config_path, e))
            return False
        self.gl_parameters.update(loaded)
        return True
    






















