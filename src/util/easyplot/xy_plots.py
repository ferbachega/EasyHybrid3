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
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk
import cairo
import random
import numpy as np
import sys

try:
    from util.colormaps import COLOR_MAPS  
except:
    from colormaps import COLOR_MAPS
from pprint import pprint

class XYPlot(Gtk.DrawingArea):
    """ Class doc """
    
    def __init__ (self, bg_color = [1,1,1], 
                        pl_color = [0,0,0], 
                        sel_color = [1,0,0], 
                        bl_color = [0.5,0.5,0.5], 
                        mode = None ):
        """ Class initialiser """
        super().__init__( )

        
        self.data = []
        
        self.bg_color = bg_color
        self.decimal = 2

        
        self.x_minor_ticks = 6
        self.x_major_ticks = 5
        
        self.y_minor_ticks = 5
        self.y_major_ticks = 10
        
        self.Xmin_list = []
        self.Ymin_list = []
        
        self.Xmax_list = []
        self.Ymax_list = []





        self.Y_space = 10


        #self.set_bg_color (bg_color[0], bg_color[1], bg_color[2])
        self.set_events(Gdk.EventMask.ALL_EVENTS_MASK)
        
        self.connect("draw", self.on_draw)
        self.connect("button_press_event", self.on_mouse_button_press)
        self.connect("motion-notify-event", self.on_motion)


        self.line_color     = pl_color
        #self.bg_color       = [1,1,1]
        self.sel_color      = sel_color
        
        #self.line_color     = [1,1,1]
        #self.bg_color       = [0,0,0]
        #self.sel_color      = [1,0,0]
        
        self.bglines_color  = bl_color
        self.bx = 100#80
        self.by = 50#50 

        self.y_botton = -1
        self.y_top    =  1

    
    def data_update (self, data):
        """ Function doc """
        X = data['X']
        Y = data['Y']
        
        self.Xmin_list.append(min(X))
        self.Ymin_list.append(min(Y))
        self.Xmax_list.append(max(X))
        self.Ymax_list.append(max(Y))
        
        
        data['Xmin'] =  min(X)
        data['Ymin'] =  min(Y)
        data['Xmax'] =  max(X)
        data['Ymax'] =  max(Y)
        
        return data


    def add (self, X = None, Y = None, 
              symbol = 'dot', sym_color = [1,1,1], sym_fill = True, 
              line = 'solid', line_color = [1,1,1], energy_label = None):
        
        """ Function doc """

        data = {
                'X'          : X     ,
                'Y'          : Y     ,
                'Xmin'       : None,
                'Ymin'       : None,
                'Xmax'       : None,
                'Ymax'       : None,
                'sym'        : symbol,
                'sym_color'  : sym_color,
                'sym_fill'   : sym_fill,
                'line'       : line  ,
                'line_color' : line_color
               }

        data = self.data_update( data)
        

        #self.Xmin_list.append(min(X))
        #self.Ymin_list.append(min(Y))
        #self.Xmax_list.append(max(X))
        #self.Ymax_list.append(max(Y))
        #
        #
        #data = {
        #        'X'          : X     ,
        #        'Y'          : Y     ,
        #        'Xmin'       : min(X),
        #        'Ymin'       : min(Y),
        #        'Xmax'       : max(X),
        #        'Ymax'       : max(Y),
        #        'sym'        : symbol,
        #        'sym_color'  : sym_color,
        #        'sym_fill'   : sym_fill,
        #        'line'       : line  ,
        #        'line_color' : line_color
        #       }
        
        self.energy_label = energy_label
        
        self.data.append(data)
        self.define_xy_limits()
    
    
    def define_xy_limits (self, y_botton = -3, y_top = 2):
        """ Function doc """
        
        self.Ymin = min(self.Ymin_list)
        self.Ymax = max(self.Ymax_list)
        

        
        self.Ymin = self.Ymin + self.Ymin/self.Y_space
        self.Ymax = self.Ymax + self.Ymax/self.Y_space
        self.Ymin = int(self.Ymin)# + self.y_botton #- 3 #round(delta/20)
        self.Ymax = int(self.Ymax) + self.y_top    #+ 2 #round(delta/20)
        delta = self.Ymax - self.Ymin

        
        #self.Ymin = -2.5
        #self.Ymax =  2.5
        
        self.Xmax = max(self.Xmax_list)
        self.Xmin = min(self.Xmin_list)
        #self.Xmax = self.Xmax - self.Xmax/20
        #self.Xmin = self.Xmin - self.Xmin/20
        
        self.deltaY  = self.Ymax - self.Ymin
        self.deltaX  = self.Xmax - self.Xmin
        
        if self.deltaX == 0:
            self.deltaX = 1
        
        if self.deltaY == 0:
            self.deltaY = 1
        

    def draw_box (self, cr, line_width = 3.0 , color = [1, 1, 1] ):
        """ Function doc """
        cr.set_line_width (3.0)
        #----------------------------------------------------------------------
        #retangle
        cr.set_source_rgb( color[0], color[1], color[2])
        cr.rectangle(self.bx,
                     self.by, 
                     self.x_box_size,
                     self.y_box_size)
        cr.stroke ()
        #----------------------------------------------------------------------
    
    def draw_XY_axes (self, cr, line_width = 3.0 , color = [1, 1, 1], font_size = 12 ):
        
        
        cr.set_line_width (1.0)
        cr.set_source_rgb( color[0], color[1], color[2])
        #---------------------------------------------------------------------------------------------------
        # x axy
        cr.move_to (self.bx, self.y_box_size+ self.by)
        cr.line_to (self.bx +self.x_box_size , self.y_box_size+ self.by)
        cr.stroke ()
        
        # y axy
        cr.move_to (self.bx, self.y_box_size+ self.by)
        cr.line_to (self.bx,  self.by)
        cr.stroke ()
        #---------------------------------------------------------------------------------------------------

    def draw_XY_scale (self, cr, line_width = 3.0 , color = [1, 1, 1], font_size =12 ):
        """ Function doc """
        
        #---------------------------------------------------------------------------------------------------
        x_major_factor = (self.x_box_size) / float(self.x_major_ticks)
        x_minor_factor = x_major_factor / float(self.x_minor_ticks)
        counter = 1
        for i in range(0, self.x_major_ticks+1):
            # marcacoes em x
            cr.move_to (self.bx+ i*x_major_factor, self.y_box_size+self.by )
            cr.line_to (self.bx+ i*x_major_factor, self.y_box_size+self.by + 15 ) 
            if i == self.x_major_ticks:
                pass
            else:
                for j in range(0, self.x_minor_ticks ):
                    cr.move_to (self.bx+ i*x_major_factor+j*x_minor_factor, self.y_box_size+self.by )
                    cr.line_to (self.bx+ i*x_major_factor+j*x_minor_factor, self.y_box_size+self.by + 10 )
            counter +=1
        #---------------------------------------------------------------------------------------------------
            #--------------------------------------------------------------------------------------------
            # texto em X
            
            cr.set_font_size(font_size)
            #text = '{:^6.2f}'.format((i/self.x_major_ticks)*self.deltaX + self.Xmin)
            text = '{:^6.'+str(self.decimal)+'f}'
            text = text.format((i/self.x_major_ticks)*self.deltaX + self.Xmin)
            
            if len(text) >7:
                num = (i/self.x_major_ticks)*self.deltaX -self.Xmax 
                text = (f"{num:.2e}")
            
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)
            
            cr.move_to(  
                        self.bx + i*x_major_factor - x_advance/2, 
                        self.by+self.y_box_size + height + 20)
            cr.show_text(text)
            #--------------------------------------------------------------------------------------------
            counter +=1


            #cr.set_font_size(font_size)
            #
            #text = str(i)
            #x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)
            ##print(x_bearing, y_bearing, width, height, x_advance, y_advance)
            ##cr.move_to(self.bx + i*self.factor_x - width/2.0    ,  self.by -15 )   
            #cr.move_to(self.bx + i*self.factor_x - x_advance/2.0 +(self.factor_x)/2   ,  self.by -15 )   
            #cr.show_text(text)







        
        #---------------------------------------------------------------------------------------------------
        y_major_factor = (self.y_box_size) / float(self.y_major_ticks)
        y_minor_factor =  y_major_factor   / float(self.y_minor_ticks)
        
        counter = 1
        self.deltaY2 =  self.deltaY
        for i in range(0, self.y_major_ticks+1):
            # marcacoes em x
            cr.move_to (self.bx, self.by+ i*y_major_factor )
            cr.line_to (self.bx - 15, self.by+ i*y_major_factor ) 
            #--------------------------------------------------------------------------------------------
            if i == self.y_major_ticks:
                pass
            else:
                for j in range(0, self.y_minor_ticks ):
                    cr.move_to  (self.bx, 
                                 self.by+ i*y_major_factor+j*y_minor_factor)
                    cr.line_to (self.bx -  10, 
                                self.by+ i*y_major_factor+j*y_minor_factor )
        
        
            #--------------------------------------------------------------------------------------------
            # texto em y
            cr.set_font_size(font_size)
            text = '{:>6.3f}'.format(self.Ymax - (i/self.y_major_ticks)*self.deltaY )
            if len(text) >6:
                num = self.Ymax - (i/self.y_major_ticks)*self.deltaY
                text = (f"{num:.2e}")
            
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)
            
            cr.move_to(  
                        self.bx -x_advance -20  , 
                        self.by+ i*y_major_factor )
            cr.show_text(text)
            #--------------------------------------------------------------------------------------------
            counter +=1
        
        
        
        
        
        
        
        
        
        
        
        
        
            counter +=1
        cr.stroke ()
        #---------------------------------------------------------------------------------------------------
        
        
        
        
        
        ''' backgound lines '''
        cr.set_source_rgb( self.bglines_color[0], 
                           self.bglines_color[1],
                           self.bglines_color[2])
        #---------------------------------------------------------------------------------------------------
        x_major_factor = (self.x_box_size) / float(self.x_major_ticks)
        x_minor_factor = x_major_factor / float(self.x_minor_ticks)
        counter = 1
        for i in range(0, self.x_major_ticks+1):
            # marcacoes em x
            cr.move_to (self.bx+ i*x_major_factor, self.y_box_size+self.by )
            cr.line_to (self.bx+ i*x_major_factor,  self.by  ) 
        #---------------------------------------------------------------------------------------------------
        
        
        #---------------------------------------------------------------------------------------------------
        y_major_factor = (self.y_box_size) / float(self.y_major_ticks)
        y_minor_factor =  y_major_factor   / float(self.y_minor_ticks)
        counter = 1
        for i in range(0, self.y_major_ticks+1):
            # marcacoes em x
            cr.move_to (self.bx, self.by+ i*y_major_factor )
            cr.line_to (self.bx +self.x_box_size, self.by+ i*y_major_factor ) 
        cr.stroke ()
        #---------------------------------------------------------------------------------------------------
        
    def on_draw(self, widget, cr):
        self.cr =  cr
        self.width = widget.get_allocated_width()
        self.height = widget.get_allocated_height()

        line_color = self.line_color
        bg_color   = self.bg_color  
        cr.set_source_rgb( bg_color[0], bg_color[1], bg_color[2])
        #line_color = [0,0,0]
        cr.paint()

        self.x = self.width
        self.y = self.height





        self.x_box_size = self.width  - (self.bx*1.5)
        self.y_box_size = self.height - (self.by*1.8)

        self.draw_box    (cr, line_width = 3.0 , color = line_color) #color = [1, 1, 1])
        self.draw_XY_axes(cr, line_width = 3.0 , color = line_color) #color = [1, 1, 1])
            
        if self.data == []:
            return False
        else:

            if self.deltaY > 0:
                for data in self.data:
                    data['Ynorm'] = [(y-self.Ymin)/self.deltaY for y in data['Y']]
                    data['Xnorm'] = [(x-self.Xmin)/self.deltaX for x in data['X']]
            else:
                for data in self.data:
                    data['Ynorm'] = data['Y']
                    data['Xnorm'] = data['X']
                
                
                
            #self.Ynorm_max = max(self.norm_data)
            
            self.draw_XY_scale ( cr, line_width = 3.0 , color = line_color)#color = [1, 1, 1] )
            
            #---------------------------------------------------------------- 
            #                          LINES                                  
            #---------------------------------------------------------------- 
            cr.set_line_width (2.0)                                           
            #---------------------------------------------------------------- 
            
            
  
            cx = 0#self.bx - self.norm_X[0]                                   
            #---------------------------------------------------------------- 
            cr.set_source_rgb( line_color[0], line_color[1], line_color[2])
            for data in self.data:
                cx = self.bx -  data['Xnorm'][0]
                
                
                #r = random.random()
                #g = random.random()
                #b = random.random()
                #rgb = [r,g,b,]
                
                
                line_color = data['line_color']
                #line_color = rgb
                #print(line_color)
                cr.set_source_rgb( line_color[0], line_color[1], line_color[2]) 
                
                
                for i in range(len(data['Ynorm'])):           
                    x = data['Xnorm'][i]    
                    y = data['Ynorm'][i]
                    
                    #color = data['line_color']
                    #color = [0,0,0]
                    #cr.set_source_rgb( color[0], color[1], color[2])
                    
                    cr.line_to (x*self.x_box_size + cx, 
                            #(self.new_Ymax + self.y_box_size + self.by ) - (y+self.y_box_size))
                            (1*self.y_box_size + self.by) - (y*self.y_box_size))
                cr.stroke ()
            #---------------------------------------------------------------- 
            #                          DOTS                                  
            #---------------------------------------------------------------- 
            #'''
            for data in self.data:
                cx = self.bx -  data['Xnorm'][0]
                
                for i in range(len(data['Ynorm'])):           
                    x = data['Xnorm'][i]    
                    y = data['Ynorm'][i]
                    color = data['line_color']
                    #color = line_color
                    
                    if data['sym'] is not None:
                        #color = data['sym_color']
                        color = data['line_color']
                        
                        #cr.set_source_rgb( color[0], color[1], color[2])
                        cr.set_source_rgb( self.sel_color[0], self.sel_color[1], self.sel_color[2])
                        
                        cr.arc (x*self.x_box_size + cx, 
                                #(self.new_Ymax + self.y_box_size + self.by ) - (y+self.y_box_size))
                                (1*self.y_box_size + self.by) - (y*self.y_box_size), 5, 0, 2*3.14)
                    
                    if data['sym_fill']:
                        cr.set_source_rgb( self.sel_color[0], self.sel_color[1], self.sel_color[2])
                        cr.fill()
                    else:
                        color = data['line_color']
                        #cr.set_source_rgb( self.line_color[0], self.line_color[1], self.line_color[2])
                        cr.set_source_rgb( color[0], color[1], color[2])
                        cr.stroke ()
    
    
    def on_motion(self, widget, event):
        '''(i/self.x_major_ticks)*self.deltaX + self.Xmin'''
        if self.data == []:
            return False 
        
        else:
            (x, y) = int(event.x), int(event.y)
            
            print("Mouse moved to:", 'x = {:10.5f} y = {:10.5f}'.format( 
                
                ((((x-self.bx)) / self.x_box_size)  *  self.deltaX + self.Xmin),
                
                ((self.y_box_size-(y-self.by)) / self.y_box_size)  *  self.deltaY + self.Ymin)
                )
    
    
    def on_mouse_button_press (self, widget, event):
        """ Function doc """
        pass


class XYScatterPlot(Gtk.DrawingArea) :
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass
        super().__init__( )
        self.data = []
        
        self.bg_color = [0,0,0]
        self.decimal = 2

        
        self.x_minor_ticks = 6
        self.x_major_ticks = 5
        
        self.y_minor_ticks = 5
        self.y_major_ticks = 10
        
        self.Xmin_list = []
        self.Ymin_list = []
        
        self.Xmax_list = []
        self.Ymax_list = []





        self.Y_space = 10


        #self.set_bg_color (bg_color[0], bg_color[1], bg_color[2])
        self.set_events(Gdk.EventMask.ALL_EVENTS_MASK)
        
        self.connect("draw", self.on_draw)
        #self.connect("button_press_event", self.on_mouse_button_press)
        #self.connect("motion-notify-event", self.on_motion)


        self.line_color     = [0,0,0]
        #self.bg_color       = [1,1,1]
        self.sel_color      = [1,0,0]
        
        #self.line_color     = [1,1,1]
        #self.bg_color       = [0,0,0]
        #self.sel_color      = [1,0,0]
        
        self.bglines_color  = [0.5,0.5,0.5]
        self.bx = 100#80
        self.by = 50#50 
        
        self.y_botton = 0
        self.y_top    = 0     

    def define_xy_limits (self, y_botton = -3, y_top = 2):
        """ Function doc """
        
        self.Ymin = min(self.Ymin_list)
        self.Ymax = max(self.Ymax_list)
        

        
        self.Ymin = self.Ymin + self.Ymin/self.Y_space
        self.Ymax = self.Ymax + self.Ymax/self.Y_space
        self.Ymin = int(self.Ymin) + self.y_botton #- 3 #round(delta/20)
        self.Ymax = int(self.Ymax) + self.y_top    #+ 2 #round(delta/20)
        delta = self.Ymax - self.Ymin

        
        #self.Ymin = -2.5
        #self.Ymax =  2.5
        
        self.Xmax = max(self.Xmax_list)
        self.Xmin = min(self.Xmin_list)
        #self.Xmax = self.Xmax - self.Xmax/20
        #self.Xmin = self.Xmin - self.Xmin/20
        
        self.deltaY  = self.Ymax - self.Ymin
        self.deltaX  = self.Xmax - self.Xmin
        
        if self.deltaX == 0:
            self.deltaX = 1
        
        if self.deltaY == 0:
            self.deltaY = 1

    def add (self, X = None, Y = None, Z = None,
              symbol = 'dot', sym_color = [1,1,1], sym_fill = True, 
              line = 'solid', line_color = [1,1,1], energy_label = None):
        
        """ Function doc """

        data = {
                'X'          : X     ,
                'Y'          : Y     ,
                'Z'          : Z     ,
                'Xmin'       : None,
                'Ymin'       : None,
                'Xmax'       : None,
                'Ymax'       : None,
                'sym'        : symbol,
                'sym_color'  : sym_color,
                'sym_fill'   : sym_fill,
                'line'       : line  ,
                'line_color' : line_color
               }

        data = self.data_update( data)
        
        
        tempZ = Z.copy()
        n = tempZ.count(None)
        for i in range(n):
            tempZ.remove(None)
        
        data['Zmax'] = max(tempZ)
        data['Zmin'] = min(tempZ)
        
        self.c_factor =  1/(data['Zmax']-data['Zmin'])
        
        
        self.energy_label = energy_label
        self.data.append(data)
        self.define_xy_limits()
  
    def data_update (self, data):
        """ Function doc """
        X = data['X']
        Y = data['Y']
        
        self.Xmin_list.append(min(X))
        self.Ymin_list.append(min(Y))
        self.Xmax_list.append(max(X))
        self.Ymax_list.append(max(Y))
        
        
        data['Xmin'] =  min(X)
        data['Ymin'] =  min(Y)
        data['Xmax'] =  max(X)
        data['Ymax'] =  max(Y)
        
        return data

    def draw_box (self, cr, line_width = 3.0 , color = [1, 1, 1] ):
        """ Function doc """
        cr.set_line_width (3.0)
        #----------------------------------------------------------------------
        #retangle
        print('draw_box')
        cr.set_source_rgb( color[0], color[1], color[2])
        cr.rectangle(self.bx,
                     self.by, 
                     self.x_box_size,
                     self.y_box_size)
        cr.stroke ()
        #----------------------------------------------------------------------

    def on_draw(self, widget, cr):
        self.cr =  cr
        self.width = widget.get_allocated_width()
        self.height = widget.get_allocated_height()

        line_color = self.line_color
        #bg_color   = self.bg_color  
        bg_color   = [1,1,1]  
        #print(bg_color)
        cr.set_source_rgb( bg_color[0], bg_color[1], bg_color[2])
        #line_color = [0,0,0]
        cr.paint()

        self.x = self.width
        self.y = self.height

        #self.data['sym_fill'] = True



        self.x_box_size = self.width  - (self.bx*1.5)
        self.y_box_size = self.height - (self.by*1.8)

        self.draw_box    (cr, line_width = 3.0 , color = line_color) #color = [1, 1, 1])
        #self.draw_XY_axes(cr, line_width = 3.0 , color = line_color) #color = [1, 1, 1])
            
        if self.data == []:
            return False
        else:
            pass
            #print(
            #'\nself.Ymax     ', self.Ymax      ,
            #'\nself.Ymin     ', self.Ymin      ,
            #'\nself.Xmax     ', self.Xmax      ,
            #'\nself.Xmin     ', self.Xmin      ,
            #'\nself.deltaY   ', self.deltaY    ,
            #'\nself.deltaX   ', self.deltaX    ,
            #)
            
            
    
            if self.deltaY > 0:
                for data in self.data:
                    data['Ynorm'] = [(y-self.Ymin)/self.deltaY for y in data['Y']]
                    data['Xnorm'] = [(x-self.Xmin)/self.deltaX for x in data['X']]
            else:
                for data in self.data:
                    data['Ynorm'] = data['Y']
                    data['Xnorm'] = data['X']
                
                
                
            #self.Ynorm_max = max(self.norm_data)
            
            #self.draw_XY_scale ( cr, line_width = 3.0 , color = line_color)#color = [1, 1, 1] )
            
            #---------------------------------------------------------------- 
            #                          LINES                                  
            #---------------------------------------------------------------- 
            cr.set_line_width (2.0)                                           
            #---------------------------------------------------------------- 
            
            
  
            cx = 0 #self.bx - self.norm_X[0]                                   
            #---------------------------------------------------------------- 
            cr.set_source_rgb( line_color[0], line_color[1], line_color[2])
            '''
            for data in self.data:
                cx = self.bx -  data['Xnorm'][0]
                
                
                #r = random.random()
                #g = random.random()
                #b = random.random()
                #rgb = [r,g,b,]
                
                
                line_color = data['line_color']
                #line_color = rgb
                #print(line_color)
                cr.set_source_rgb( line_color[0], line_color[1], line_color[2]) 
                
                
                for i in range(len(data['Ynorm'])):           
                    x = data['Xnorm'][i]    
                    y = data['Ynorm'][i]
                    
                    #color = data['line_color']
                    #color = [0,0,0]
                    #cr.set_source_rgb( color[0], color[1], color[2])
                    
                    cr.line_to (x*self.x_box_size + cx, 
                            #(self.new_Ymax + self.y_box_size + self.by ) - (y+self.y_box_size))
                            (1*self.y_box_size + self.by) - (y*self.y_box_size))
                cr.stroke ()
            #'''
            #---------------------------------------------------------------- 
            #                          DOTS                                  
            #---------------------------------------------------------------- 
            #'''
            for data in self.data:
                cx = self.bx -  data['Xnorm'][0]
                
                for i in range(len(data['Ynorm'])):           
                    x = data['Xnorm'][i]    
                    y = data['Ynorm'][i]
                    color = data['line_color']
                    color = [0,0,0]
                    data['sym_fill'] = True
                    
                    if data['sym'] is not None:
                        color = data['sym_color']
                        
                        #cr.set_source_rgb( color[0], color[1], color[2])
                        #cr.set_source_rgb( self.sel_color[0], self.sel_color[1], self.sel_color[2])
                        
                        cr.arc (x*self.x_box_size + cx, 
                                #(self.new_Ymax + self.y_box_size + self.by ) - (y+self.y_box_size))
                                (1*self.y_box_size + self.by) - (y*self.y_box_size), 5, 0, 2*3.14)
                    
                        if data['sym_fill']:
                            
                            if data['Z'][i]:
                                #print(data['Z'][i]* 1/42.0)
                                f = data['Z'][i] * 1/42.0
                                #self.sel_color[0]
                                cr.set_source_rgb( f, 0, f)
                                cr.fill()
                            else:
                                cr.set_source_rgb( 1, 1, 1)
                            
                                cr.fill()
                    cr.stroke ()
