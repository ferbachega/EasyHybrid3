#!/usr/bin/env python3
# -*- coding                       : utf-8 -*-
#
#  vismol_config.py
#  
#  Copyright 2022 Fernando <fernando@winter>
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
                                      "line_type"                  : 0,
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
                                      
                                      "sticks_color"               : 0,
                                      "sticks_type"                : 0,
                                      "antialias"                  : True,
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
                                      
                                      #"pk_label_size"              : [1.0, 1.0, 1.0, 1.0],
                                      
                                      'startup_path'               : None,
                                      'tmp_files'                  : None,
                                      'autosave'     : True, 
                                      'askSaveUnsave': True, 
                                      #'askSaveUnsave': True, 
                                      }
                              
        self.n_proc = 2
        # self.representations_available = {"dots", "lines", "nonbonded", "dotted_lines",
        #                                   "ribbon", "sticks", "spheres", "impostor",
        #                                   "surface", "cartoon", "freetype",
        #                                   "picking_dots"}
        
        #.Rep list - Don't change this list
        self.representations_available = {"dots", "lines", "nonbonded", "impostor",'dash', #"cartoon",
                                          "sticks", "spheres", 'ribbons',#'ribbon_sphere', 
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
        #print('AQUIIIIIIIIIIIIIII')
        
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
                    print("Loading EasyHybrid config file.")
                    with open(config, 'r', encoding='utf-8') as f:
                        self.gl_parameters = json.load(f)
                except:
                    print("Failed to open EasyHybrid configuration file. Loading default settings.")
                    self.gl_parameters = self.gl_parameters_default
                
            else:
                print("Configuration file not found. Creating a new file from default settings.")
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
            
            print("Saving configuration file.")
            
            with open(config, 'w', encoding='utf-8') as f:
                json.dump(self.gl_parameters, f, ensure_ascii=False, indent=4)        


    
    def load_easyhybrid_config(self, config_path)                       :
        """ Function doc """
        if not os.path.isfile(config_path)                       :
          config_path = os.path.join(os.environ["HOME"], ".VisMol", "VismolConfig.json")
        with open(config_path, "r") as config_file                       :
            self.gl_parameters = json.load(config_file)
    























