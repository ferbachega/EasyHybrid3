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
            print(cmap)
        
        
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
            
            for i in range(0, len(self.data['Z'])):
                for j in range(0, len(self.data['Z'][i])):
                    self.data['Z'][i][j] = self.data['Z'][i][j]-_min 
            
            self.plot.show()
            self.plot.data = self.data['Z']
        
            
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
            If not, the function returns Flase and the click points are erased.
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
                    else:
                        xy_list = [[i_on_plot, j_on_plot]]
                
                #when the interpolation option is not active.
                else:
                    xy_list = [[i_on_plot, j_on_plot]]
                
                
                for xy in xy_list:
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
                
                self.plot2.add( X = x, Y = y, symbol = 'dot', sym_fill = False, sym_color = [0,1,1], line = 'solid', line_color = [0,1,1] )
                
                
                
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

            #print(self.xy_traj[int(value)])
            xy = self.plot.points[int(value)]
            #print(xy, self.zdata[int(value)])
            #print(xy, self.plot.data[xy[0]][xy[1]])
            x = [value]
            y = [  self.plot.data[xy[0]][xy[1]]   ]
            frame = self.vobject.idx_2D_xy[(xy[1], xy[0])]
        else:
            x = [int(value)]
            y = [self.data['Z'][int(value)]]
            frame = int(value)
        
        
        
        
        if len(self.plot2.data) > 1:
            self.plot2.data.pop(-1)
        
        self.plot2.add( X = x, Y = y, symbol = 'dot', sym_fill = True, sym_color = [1,0,0], line = 'solid', line_color = [0,1,1] )
        self.plot2.queue_draw()

        
        if self.vobject:
            self.main.vm_session.frame = int(frame)
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
    return input_coord[1:]

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

