#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  easyhybrid_pDynamo_selection.py
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
import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
import gc
import os
from util.easyplot import ImagePlot, XYPlot
from util.colormaps import COLOR_MAPS  

#import math
from gi.repository import Gtk, Pango
from gi.repository import Gdk
from gi.repository import GdkPixbuf
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.windows.setup.windows_and_dialogs import InfoWindow

from pprint import pprint
import numpy as np

class PotentialEnergyAnalysisWindow:

    def __init__(self, main = None ):
        """ Class initialiser """
        self.main                = main
        self.p_session           = self.main.p_session
        self.vm_session          = self.main.vm_session
        self.Visible             =  False  
        
        self.cmap_id = 0
        
        self.vobject_liststore = Gtk.ListStore(str,              # name
                                               int,              # vobj_id
                                               int,              # sys_id
                                               GdkPixbuf.Pixbuf) # PixBuf
        
        self.data_liststore    = Gtk.ListStore(str, int)
        
        # plotting attributes
        self.interpolate = True
        
        self.cmap_ref_dict = {}
        self.cmap_store = Gtk.ListStore(str)
        cmaps = COLOR_MAPS.keys()
        
        for i, cmap in enumerate(cmaps):
            self.cmap_store.append([cmap])
            self.cmap_ref_dict[i] = cmap
            #print(cmap)
        
        #----------------------- RC labels -----------------------------
        self.RC_label_list = ['frames', 'rc distances']
        self.RC_type_store = Gtk.ListStore(str)
        for label in self.RC_label_list:
            self.RC_type_store.append([label])
        #---------------------------------------------------------------

        self.menu_items = {
                          
                          #'header'                  : None                            ,
                                                    
                          '_separator'              : ''                              ,
                                                    
                          'Info'                    : self._show_info                  ,
                                                    
                          '_separator'              : ''                              ,
                          'Settings'                : self._menu_settings             ,
                          '_separator'              : ''                              ,

                          
                          'Optimize Pathway'        : self._menu_opt_pathway           ,
                          'Export Incomplete Matrix': self._menu_export_incomplete_data,
                          'Reduce Resolution'       : self._menu_change_data_resolution,

                          
                          
                          '_separator'              : ''                              ,
                          
                          'Delete'                  : self._menu_delete_system        ,

                          }
            
        
        
    def OpenWindow (self, vobject = None):
        """ Function doc """
        if self.Visible  ==  False:

            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.main.home,'src/gui/windows/analysis/PES_analysis_window.glade'))
            self.builder.connect_signals(self)

            self.vobject = vobject

            self.window = self.builder.get_object('window')
            self.window.set_title('PES Analysis Window')
            #self.window.set_keep_above(True)            
            self.window.set_default_size(700, 450)
            
            
            self.grid = self.builder.get_object('grid_setup')
            self.hbox = self.builder.get_object('hbox_plotting')
            
            self.checkbox_smooth = self.builder.get_object('checkbox_smooth')
            self.checkbox_smooth.connect('toggled', self.on_checkbox_smooth_toggle)
            
            self.cmap_box = self.builder.get_object('cmap_box')          
            self.cmap_combobox = Gtk.ComboBox.new_with_model(self.cmap_store)
            renderer_text = Gtk.CellRendererText()
            self.cmap_combobox.pack_start(renderer_text, True)
            self.cmap_combobox.add_attribute(renderer_text, "text", 0)            
            self.cmap_box.add(self.cmap_combobox)
            self.cmap_combobox.connect('changed', self.on_cmap_combobox_change)
            
            '''-------------------------------------------------------------'''
            self.RC_label = self.builder.get_object('RC1_RC2_label')

            self.plot = ImagePlot()
            self.plot.RC_label = self.RC_label
            #t = np.linspace(0, 2 * np.pi, 40)
            #data2d = np.sin(t)[:, np.newaxis] * np.cos(t)[np.newaxis, :]
            self.hbox.pack_start(self.plot, True, True, 0)
            #self.plot.data = data2d
            self.plot.connect("button_press_event", self.on_mouse_button_press)
            
            self.RC_type_box = self.builder.get_object('RC_type_box')          
            self.RC_type_combobox = Gtk.ComboBox.new_with_model(self.RC_type_store)
            renderer_text = Gtk.CellRendererText()
            self.RC_type_combobox.pack_start(renderer_text, True)
            self.RC_type_combobox.add_attribute(renderer_text, "text", 0)            
            self.RC_type_box.add(self.RC_type_combobox)
            self.RC_type_combobox.connect('changed', self.on_RC_type_combobox)
            self.RC_type_combobox.set_active(0)
            '''-------------------------------------------------------------'''


            '''-------------------------------------------------------------'''
            self.plot2 = XYPlot()

            self.plot2.x_minor_ticks = 1
            self.plot2.x_major_ticks = 5
            self.plot2.y_minor_ticks = 6
            self.plot2.y_major_ticks = 5
            self.plot2.decimal = 0
            self.hbox.pack_start(self.plot2, True, True, 0)
            '''-------------------------------------------------------------'''
            #self.plot2.connect("button_press_event", self.on_mouse_button_press)

            self.coordinates_combobox = CoordinatesComboBox(coordinates_liststore = self.vobject_liststore)
            self.coordinates_combobox.connect('changed', self.on_coordinates_combobox_change)
            self.grid.attach (self.coordinates_combobox, 1, 0, 1, 1)
            self.refresh_vobject_liststore ()


            





            #------------------------------------------------------------------------------------
            self.data_combobox = Gtk.ComboBox()
            renderer_text = Gtk.CellRendererText()
            self.data_combobox.pack_start(renderer_text, True)
            self.data_combobox.add_attribute(renderer_text, "text", 0)
            self.data_combobox.set_model(self.data_liststore)
            self.data_combobox.connect('changed', self.on_data_combobox_change)
            self.grid.attach (self.data_combobox, 3, 0, 1, 1)
            #------------------------------------------------------------------------------------


            self.threshold_entry= self.builder.get_object('threshold_entry')
            self.threshold_entry.connect('activate', self.on_threshold_entry_activate)

            #------------------------------------------------------------------------------------
            self.scale_traj = self.builder.get_object('scale_trajectory_from_PES')
            self.scale_traj.connect('change-value', self.on_scaler_frame_value_changed )
            self.adjustment     = Gtk.Adjustment(value         = 0,
                                                 lower         = 0,
                                                 upper         = 100,
                                                 step_increment= 1,
                                                 page_increment= 1,
                                                 page_size     = 1)

            self.scale_traj.set_adjustment ( self.adjustment)
            self.scale_traj.set_digits(0)
            #------------------------------------------------------------------------------------
            #self.scale_trajectory_energy_label = self.builder.get_object('scale_trajectory_energy_label')
            self.menu, x = self.build_tree_view_menu(self.menu_items)
            


            self.cmap_combobox.set_active(self.cmap_id)

            self.window.show_all()
            self.Visible  = True
            self.coordinates_combobox.set_active(0)



    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False

    def refresh_vobject_liststore (self):
        """ Function doc """       
        self.vobject_liststore.clear()
        for system_id , system in self.p_session.psystem.items():
            
            logfile_data = getattr ( system, "e_logfile_data", None )
            if logfile_data is not None:
                for vobject_id in system.e_logfile_data.keys():
                    try:
                        _vobject = self.vm_session.vm_objects_dic[vobject_id]
                        
                        pixbuf  = self.main.vobject_liststore_dict[system_id][_vobject.liststore_iter][3]
                        
                        self.vobject_liststore.append([_vobject.name, 
                                                        vobject_id  , 
                                                        system_id   , 
                                                        pixbuf     ])
                    except:
                        print('Log data not found!')
            else:
                pass

    def on_threshold_entry_activate (self, widget):
        """ Function doc """
        try:
            new_threshold = float(self.threshold_entry.get_text())
            self.plot.set_threshold_color ( _min = 0, _max = new_threshold)
            self.plot.queue_draw()
        except:
            pass
   
    def on_checkbox_smooth_toggle (self, widget):
        """ Function doc """
        if widget.get_active():
            self.plot.is_discrete = False
        else:
            self.plot.is_discrete = True
        self.plot.queue_draw()
        #self.builder.get_object('checkbox_interpolate').get_active()

    def on_RC_type_combobox (self, widget):
        self.RC_type = widget.get_active()
        print(self.RC_type)
        self.plot.set_label_mode(mode = self.RC_type)
        self.plot.queue_draw()

    def on_cmap_combobox_change (self, widget):
        """ Function doc """
        self.cmap_id = widget.get_active()
        print(self.cmap_id, self.cmap_ref_dict[self.cmap_id])
        self.plot.cmap = self.cmap_ref_dict[self.cmap_id]
        self.plot.set_cmap (cmap = self.cmap_ref_dict[self.cmap_id])
        
        self.on_threshold_entry_activate (None)
        
    def on_coordinates_combobox_change (self, widget):
        """ Function doc """
        try:
            vobject_index = self.coordinates_combobox.get_vobject_id()
            self.vobject  = self.main.vm_session.vm_objects_dic[vobject_index]
            
            self.data_liststore.clear()
            for index , data in enumerate(self.p_session.psystem[self.vobject.e_id].e_logfile_data[vobject_index]):
                self.data_liststore.append([data['name'], index])

            self.data_combobox.set_active(0)
        except:
            pass
        #print( self.vobject.idx_2D_xy)
    
    def on_data_combobox_change (self, widget):
        """ Function doc """
        _iter = self.data_combobox.get_active_iter()
        if _iter is not None:
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.data_combobox.get_model()
            _name, index = model[_iter][:2]
        else:
            return False
        self.data = self.p_session.psystem[self.vobject.e_id].e_logfile_data[self.vobject.index][index] 
        #print('data: ', self.data)
        
        if self.data['type'] == 'plot1D':
            pass
            self.plot2.data = []

            self.plot2.Xmin_list = []
            self.plot2.Ymin_list = []
            self.plot2.Xmax_list = []
            self.plot2.Ymax_list = []
           
            _min = None  
            for value in self.data['Z']:
                if _min == None:
                    _min = value
                else:
                    if value < _min:
                        _min = value
                    else:
                        pass
            
            for index, value in enumerate(self.data['Z']):
                self.data['Z'][index] = value-_min
            
            
            self.plot2.add( X = range(0, len(self.data['Z'])), Y = self.data['Z'], symbol = 'dot', sym_fill = False, sym_color = [0,1,1], line = 'solid', line_color = [0,1,1] )
            self.plot.hide()
            self.scale_traj_new_definitions(set_range = len(self.data['Z']))
        
        elif self.data['type'] == 'plot2D':
        
            minlist = []
            for line in self.data['Z']:
                minlist.append(min(line))
            
            _min = min(minlist)
            print(self.data)
            for i in range(0, len(self.data['Z'])):
                for j in range(0, len(self.data['Z'][i])):
                    self.data['Z'][i][j] = self.data['Z'][i][j]-_min 
            #print('\n\n\n', self.data,'\n\n\n')
            self.plot.show()
            self.plot.data = self.data['Z']
            self.plot.dataRC1 = self.data['RC1']
            self.plot.dataRC2 = self.data['RC2']
            
            
            data_min, data_max, delta , norm_data = self.plot.define_datanorm()
            
            self.builder.get_object('threshold_entry').set_text(str(data_max))
        
        self.plot.queue_draw()


    def on_mouse_button_press (self, widget, event):
        """ Function doc """
        if event.button == 1:
            #self.is_button_1_pressed = True # mouse left button
        
            '''get the points on the graph where the click was made'''
            i_on_plot, j_on_plot, x, y = widget.get_i_and_j_from_click_event(event)
            
            
            '''
            checks if the place where click happened is outside the graphic or not. 
            If not, the function returns False and the click points are erased.
            '''
            if widget.check_clicked_points( i_on_plot, j_on_plot):
                
                '''
                Interpolation takes the intermediate points between two 
                successive clicks.
                '''
                if self.builder.get_object('checkbox_interpolate').get_active():
                #if self.interpolate:
                    if len(widget.points) > 0 :
                        xy_list = build_chain_of_states( [widget.points[-1], [i_on_plot, j_on_plot]])
                        xy_list = xy_list[1:] # the first point is already in the yx list - should be deleted from the list of interpoated points
                    else:
                        xy_list = [[i_on_plot, j_on_plot]]
                
                #when the interpolation option is not active.
                else:
                    xy_list = [[i_on_plot, j_on_plot]]
                
                
                for xy in xy_list:#[1:]:
                    widget.points.append(xy)
                
                
                self.plot2.data = []
    
                self.plot2.Xmin_list = []
                self.plot2.Ymin_list = []
                self.plot2.Xmax_list = []
                self.plot2.Ymax_list = []
    
    
                
                
                
            
                x = []
                y = []
                for i, point in enumerate(widget.points):
                    x.append(i)
                    y.append(widget.data[point[0]][point[1]])
                
                self.plot2.add( X = x, Y = y, symbol = 'dot', sym_fill = False, sym_color = [0,0,1], line = 'solid', line_color = [0,0,0] )
                
                
                
                if len(x)-1 == 0:
                    self.plot2.x_major_ticks = 1
                else:
                    x_major_ticks = len(x)-1  
                    self.plot2.x_major_ticks = x_major_ticks
                self.plot2.Xmax   = 10 
                #self.plot2.x_major_ticks = 10
                print("Mouse clicker at:",  x, y, 
                                            int(i_on_plot), int(j_on_plot), 
                                            widget.data[int(i_on_plot)][int(j_on_plot)] )
            
            self.scale_traj_new_definitions()
            self.plot2.queue_draw()
            widget.queue_draw()
        
        if event.button ==3:
            self.menu.popup(None, None, None, None, 0, 3)
            
    def scale_traj_new_definitions(self, set_range = None ):
        #self.scale_traj
        if set_range:
            self.scale_traj.set_range(0, set_range)
        else:    
            self.scale_traj.set_range(0, len(self.plot.points))
        #self.scale_traj.set_value(self.vm_session.get_frame())
    
    def on_scaler_frame_value_changed (self, hscale, text= None,  data=None):
        """ Function doc """
        value = self.scale_traj.get_value()
        pos   = self.scale_traj.get_value_pos ()
        
        #self.plot.points
        
        if self.data['type'] == 'plot2D':
            #print(self.plot.data)
            #print(self.xy_traj[int(value)])
            
            if len(self.plot.points) >1:
                xy = self.plot.points[int(value)]
                
                self.plot.selected_dot = [xy[0],xy[1]]

                #print(xy, self.plot.data[xy[0]][xy[1]])
                x = [value]
                y = [  self.plot.data[xy[0]][xy[1]]   ]

                text = 'E = {:<15.6f}'.format(y[0])
                self.builder.get_object('label_energy').set_text(text)

                frame = self.vobject.idx_2D_xy[( xy[0],xy[1])]
            else:
                pass
        else:
            x = [int(value)]
            y = [self.data['Z'][int(value)]]
            
            text = 'E = {:<15.6f}'.format(y)
            self.builder.get_object('label_energy').set_text(str(y))
            
            frame = int(value)
        
        
        
        
        if len(self.plot2.data) > 1:
            self.plot2.data.pop(-1)
        
        
        try:
            self.plot2.add( X = x, Y = y, symbol = 'dot', sym_fill = True, sym_color = [1,0,0], line = 'solid', line_color = [0,1,1] )
            self.plot2.queue_draw()
        except:
            pass
            
        self.plot.queue_draw()
        #print('E = {}'.format( y[0]) )
        #v = y[0]
        #text = 'E = {6.3f}'.format(v)
        #self.scale_trajectory_energy_label.set_text( 'E = {:6.3f}'.format( y[0])  )
        if self.vobject:
            self.main.vm_session.set_frame(int(frame))
            #self.main.vm_session.frame = int(frame)set_frame
            self.main.vm_session.vm_glcore.updated_coords = True
            self.main.vm_session.vm_widget.queue_draw()
        
        '''
        self.ax2.cla()
        self.ax3.cla()
        self.ax2.plot(range(0, len(self.zdata)), self.zdata, '-ob')
        self.ax3.plot( [int(value)], [self.zdata[int(value)]], '-or')
        
        if self.vobject:
            frame = self.vobject.idx_2D_xy[(xy[0], xy[1])]
            self.main.vm_session.frame = int(frame)
            #self.main.vm_session.set_frame(int(frame)) 
            self.main.vm_session.vm_glcore.updated_coords = True
            self.main.vm_session.vm_widget.queue_draw()
        self.canvas2.draw()
        #self.scale_traj.set_value(int(value))
        '''


    def on_button_export_trajectory (self, widget):
        """ Function doc 
            RC1 and RC2 refer to the coordinates of reactions 1 and 2, respectively, defined to obtain the PES. In the case of a simple reaction coordinate, the RC is simply the value of the distance between atoms 1 and 2 defined. In the case of the multiple distance based reaction coordinate defined by three (3) atoms, the RC is defined as the distance between atoms 1 and 2 minus the distance between atoms 2 and 3.
        
        """

        active_id = self.main.p_session.active_id
        
        new_vismol_object = self.main.p_session.generate_new_empty_vismol_object(system_id = self.vobject.e_id , 
                                                                                 name      = 'new_coordinates' )
        print('\n\n')
        
        text = ''
#The RC1 and RC2 refer to the reaction coordinates I and II, respectively. In case of RC is a  
#simple reaction coordinate, RC is simply the distance value between atoms #1 and #2. For 
#multiple distance based reaction coordinate defined by three (3) atoms, the 
#RC is defined as the distance between atoms 1 and 2 minus the distance between atoms 2 and 3.
#        ''' 
        
        text += '\n  (i)   (j)   RC1       RC2       ENERGY(kJ/mol)' 
        for xy in self.plot.points:
            frame_number = self.vobject.idx_2D_xy[(xy[1], xy[0])]
            #print(self.vobject.frames[frame_number])
            new_vismol_object.frames = np.vstack((new_vismol_object.frames, 
                                                  [self.vobject.frames[frame_number]]))
        
            
            
            text += '\n {:3d}   {:3d}    {:3.4f}    {:3.4f}    {:3.6f}'.format(xy[1], xy[0], self.data['RC1'][xy[0]][xy[1]],self.data['RC2'][xy[0]][xy[1]], self.data['Z'][xy[0]][xy[1]] )
            #print(text)
        print (text)
        print('\n\n')

        self.main.p_session._apply_fixed_representation_to_vobject(system_id = None, vismol_object = new_vismol_object)
        self.main.p_session._apply_QC_representation_to_vobject   (system_id = None, vismol_object = new_vismol_object) 
        self.vm_session.main.main_treeview.refresh_number_of_frames()
        
        system = self.main.p_session.psystem[self.vobject.e_id]
        window = InfoWindow(system, text)
        
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        #active_id = self.main.p_session.active_id
        #new_vismol_object = self.main.p_session.generate_new_empty_vismol_object(system_id = self.vobject.e_id , name = 'new_coordinates' )
        #print('  (i)   (j)   RC1       RC2       ENERGY(kJ/mol)')
        

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    def  update (self):
        """ Function doc """
        pass

    
    def build_tree_view_menu (self, menu_items = None):
        """ Function doc """
        tree_view_menu = Gtk.Menu()
        menu_header    = None
        
        for label in menu_items:

            if menu_items[label] == None:
                # just a label
                
                mitem = Gtk.MenuItem(label = label)
                if label == 'header':
                    menu_header    = mitem
                
                
            elif  label == '_separator':
                mitem = Gtk.SeparatorMenuItem()
            
            else:
                mitem = Gtk.MenuItem(label = label)
                mitem.connect('activate', menu_items[label])
            
            tree_view_menu.append(mitem)

        tree_view_menu.show_all()
        return tree_view_menu, menu_header




    def _show_info                   (self, widget):
        """ Function doc """

    def _menu_settings               (self, widget):
        """ Function doc """

    def _menu_opt_pathway            (self, widget):
        """ Function doc """
        print('Here - _menu_opt_pathway')
        
        input_coord = self.plot.points
        e_matrix    = self.plot.data

        self.plot.points = run_surface_NEB (input_coord = input_coord, e_matrix = e_matrix )
        
        
        
        self.plot2.data = []

        self.plot2.Xmin_list = []
        self.plot2.Ymin_list = []
        self.plot2.Xmax_list = []
        self.plot2.Ymax_list = []

        x = []
        y = []
        for i, point in enumerate(self.plot.points):
            x.append(i)
            y.append(self.plot.data[point[0]][point[1]])
        
        self.plot2.add( X = x, Y = y, symbol = 'dot', sym_fill = False, sym_color = [1,1,1], line = 'solid', line_color = [0,0,0] )

        if len(x)-1 == 0:
            self.plot2.x_major_ticks = 1
        else:
            x_major_ticks = len(x)-1  
            self.plot2.x_major_ticks = x_major_ticks
        self.plot2.Xmax   = 10 


        
        
        
        
        
        
        
        self.scale_traj_new_definitions()
        self.plot.queue_draw()
        self.plot2.queue_draw()
        
    def _menu_export_incomplete_data (self, widget):
        """ Function doc """
        
        
        #folder = '/home/fernando/test'
        system = self.main.p_session.psystem[self.vobject.e_id]
        folder = system.e_working_folder
        
        folder = os.path.join(folder, 'Incomplete_matrix.ptGeo')
        
        if os.path.isdir( folder):
            pass
        else:
            os.mkdir(folder)
        
        x_disp   = 2
        y_disp   = 2
        x_offset = 1
        y_offset = 1
        
        for xy in self.plot.points:
            for i in range(-x_disp, x_disp+1, x_offset):
                for j in range(-y_disp, y_disp+1, y_offset):
                    
                    y = xy[1] + j
                    x = xy[0] + i
                    #frame_number = self.vobject.idx_2D_xy[(xy[1], xy[0])]
                    
                    try:
                        frame_number = self.vobject.idx_2D_xy[(y, x)]
                        self.main.p_session.get_coordinates_from_vobject_to_pDynamo_system(vobject   = self.vobject, 
                                                                        system_id = self.vobject.e_id, 
                                                                        frame     = frame_number)
                        filename = 'frame{}_{}.pkl'.format(y, x)
                        self.main.p_session.export_pdynamo_system_coordinates( folder, filename, system)
                        #Pickle( os.path.join ( folder, filename), system.coordinates3 )
                    except:
                        pass
                        
    def _menu_change_data_resolution (self, widget):
        """ Function doc """

    def _menu_delete_system          (self, widget):
        """ Function doc """
















def find_the_midpoint (coord1 , coord2):
    """ Function doc """
    #print (coord1, coord2)
    
    x = float(coord2[0] - coord1[0])
    x = (x/2)
    #print ('x', x)
    x = coord1[0] + x
    
    
    y = float(coord2[1] - coord1[1])
    y = (y/2)
    #print ('y', y)
    
    y = coord1[1] + y
    
    #print (x, y)
    #return [int(x), int(round(y))]
    return [int(x), int( y )]


def build_chain_of_states( input_coord):
 
    #print (input_coord)
    inset_point = True


    while inset_point == True:
        a = 0
        counter = 0

        while a == 0:

            try:
                coord1 =  input_coord[ counter  ]
                coord2 =  input_coord[ counter+1]

                inset_point = check_distance (coord1 , coord2)
                
                if inset_point == False:
                    counter += 1
                    #print ('inset_point == False')
                else:
                    midpoint = find_the_midpoint (coord1 , coord2)
                    #print counter, counter+1, midpoint, input_coord
                    input_coord.insert(counter+1, 0 )
                    input_coord[counter+1] = midpoint
            except:
                a = True
                #print input_coord
    return input_coord 


def check_distance (coord1 , coord2):
    """ Function doc """
    x = float(coord2[0] - coord1[0])
    y = float(coord2[1] - coord1[1])
    d =  (x**2 + y**2)**0.5
    if d < 1.42:
        return False
    else:
        #print 'not too close'
        return True


def get_gradient_from_matrix_line (matrix_line, position):
    """ Function doc """
    d1 =  matrix_line[position-1] - matrix_line[position]
    #print position,'d1' , d1
    if d1 <= 0: 
        return -1 

    d2 =  matrix_line[position+1] - matrix_line[position]
    #print position, 'd2' , d2
    if d2 <= 0:
        return 1 
    return 0


def run_surface_NEB (input_coord = None, e_matrix = None, k = 1.5  ):
    """ Function doc """
    
    '''
    from the older script
    data = log_parser (inputs['log_file'])
    data = data[1]
    k    = inputs['kf']
    e_matrix    = data['matrix']

    reac = [0,0]
    #ts   = [15,30]
    prod = [20,20]
    
    
    #
    Initial_positions = [reac, [4,10], prod]
    #Initial_positions = [reac, ts,  prod]
    '''
    
    print ('Matrix size:', len(e_matrix[0]), len(e_matrix))
    
    initial_size = len(input_coord)
    xy_surface_positions   = build_chain_of_states(input_coord)
    final_size   = len(xy_surface_positions)

    #------------------------------------------------------------------
    print ("\n\nInitial chain of states: \n")
    for xy_coord in xy_surface_positions:
        x = xy_coord[0]
        y = xy_coord[1]
        print (x, y, e_matrix[x][y])
    #------------------------------------------------------------------

    #------------------------------------------------------------------
    counter                     = None
    total_sum_grad_ateriror     = 0
    delta                       = None
    chain_of_states_convergence = False
    

    chain_of_states_convergence_max_interactions = 100
    chain_of_states_convergence_counter = 0

    #"""
    while  chain_of_states_convergence != True:
        '''
        Optmizing the chain of states 
        ----------------------------------------------------------------------------------------
        '''
        #print('primeiro while', chain_of_states_convergence_counter)
        chain_of_states_perturbation_max_interactions = 100
        chain_of_states_perturbation_counter = 0
        
        while delta != 0:
            #print('segundo while')

            total_sum_grad = 0
            
            for xy_coord in xy_surface_positions:
                index       = xy_surface_positions.index(xy_coord)
                final_index = len(xy_surface_positions)
            
                if index == 0 or index  == final_index-1:
                    pass
                
                else:
                    xy_coord_before = xy_surface_positions[index-1]
                    xy_coord_after  = xy_surface_positions[index+1]
                    
                    x_before   = xy_coord_before[0]
                    y_before   = xy_coord_before[1]
                    
                    x_midpoint = xy_coord[0]
                    y_midpoint = xy_coord[1]
                
                    x_after    = xy_coord_after[0]
                    y_after    = xy_coord_after[1]


                    #----------------------- X   perturbations  ---------------------------------------
                    energy_left  = k*( (x_midpoint-1) - x_before )**2 + k*( x_after - (x_midpoint-1) )**2
                    energy_midle = k*(  x_midpoint    - x_before )**2 + k*( x_after -  x_midpoint    )**2
                    energy_right = k*( (x_midpoint+1) - x_before )**2 + k*( x_after - (x_midpoint+1) )**2
                
                    print(  x_before, x_midpoint, x_after, energy_left, energy_midle , energy_right )
                    energy_left__from_matrix = e_matrix[x_midpoint-1][y_midpoint]
                    energy_midle_from_matrix = e_matrix[x_midpoint  ][y_midpoint]
                    energy_right_from_matrix = e_matrix[x_midpoint+1][y_midpoint]
                    
                    sum_l = energy_left  + energy_left__from_matrix
                    sum_m = energy_midle + energy_midle_from_matrix
                    sum_r = energy_right + energy_right_from_matrix
                    
                    #klinh = 1
                    #sum_l = klinh  + energy_left__from_matrix
                    #sum_m =       energy_midle_from_matrix
                    #sum_r = klinh + energy_right_from_matrix
                    
                    d1 = sum_l - sum_m
                    d2 = sum_r - sum_m  

                    if d2 < 0:
                        asn = 1 
                        total_sum_grad += d1

                    if d1 < 0: 
                        asn = -1 	
                        total_sum_grad += d2

                    if d1  > 0  and d2 > 0:
                        asn = 0
                    
                    xy_surface_positions[index][0] += asn
                    #-------------------------------------------------------------------------------------
            

                    #-----------------------   Y    perturbations  ---------------------------------------
                    energy_left  = k*( (y_midpoint-1) - y_before )**2 + k*( y_after - (y_midpoint-1) )**2
                    energy_midle = k*(  y_midpoint    - y_before )**2 + k*( y_after -  y_midpoint    )**2
                    energy_right = k*( (y_midpoint+1) - y_before )**2 + k*( y_after - (y_midpoint+1) )**2
                
                    
                    energy_left__from_matrix = e_matrix[x_midpoint][y_midpoint-1]
                    energy_midle_from_matrix = e_matrix[x_midpoint][y_midpoint  ]
                    energy_right_from_matrix = e_matrix[x_midpoint][y_midpoint+1]
                    
                    sum_l = energy_left  + energy_left__from_matrix
                    sum_m = energy_midle + energy_midle_from_matrix
                    sum_r = energy_right + energy_right_from_matrix
                    
                    #sum_l = klinh  + energy_left__from_matrix
                    #sum_m =       energy_midle_from_matrix
                    #sum_r = klinh + energy_right_from_matrix
                    
                    
                    d1 = sum_l - sum_m
                    d2 = sum_r - sum_m  

                    if d2 < 0:
                        asn = 1 
                        total_sum_grad += d1

                    if d1 < 0: 
                        asn = -1 
                        total_sum_grad += d2

                    if d1  > 0  and d2 > 0:
                        asn = 0
                    
                    xy_surface_positions[index][1] += asn
                    #-------------------------------------------------------------------------------------

            delta = total_sum_grad_ateriror - total_sum_grad
            total_sum_grad_ateriror = total_sum_grad
            print('delta: ' , delta)

            old = None
            for xy_coord in xy_surface_positions:
                #                    X     Y
                index  = xy_surface_positions.index(xy_coord)
                
                if xy_coord == old:
                    xy_surface_positions.pop(index)
                
                else:
                    x = xy_coord[0]
                    y = xy_coord[1]
                old = xy_coord
            
            chain_of_states_perturbation_counter += 1
            
            if chain_of_states_perturbation_counter == chain_of_states_convergence_max_interactions:
                delta =  0
                
        chain_of_states_convergence_counter += 1
        
        old = xy_surface_positions
        xy_surface_positions = build_chain_of_states(xy_surface_positions)
        
        if old == xy_surface_positions:
            chain_of_states_convergence = True
        else:
            pass
        
        
        if chain_of_states_convergence_counter == chain_of_states_convergence_max_interactions:
            chain_of_states_convergence = True
        else:
            pass
                
                
                
    #print xy_surface_positions
    print ("\n\nOptimized chain of states: \n")
    for xy_coord in xy_surface_positions:
        #                    X     Y
        x = xy_coord[0]
        y = xy_coord[1]
        print (x, y, e_matrix[x][y])                
    #"""
    return xy_surface_positions
        
                
                
                
                
                
                
                
