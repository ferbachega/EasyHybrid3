#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  plotting.py
#  
#  Copyright 2023 Fernando <fernando@winter>
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
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, cairo
import numpy as np
import sys






COLOR_MAPS = {
             'jet': {
                0.0 : [0 /255, 0 /255, 143 /255]         ,
                0.015625 : [0 /255, 0 /255, 159 /255]    ,
                0.03125 : [0 /255, 0 /255, 175 /255]     ,
                0.046875 : [0 /255, 0 /255, 191 /255]    ,
                0.0625 : [0 /255, 0 /255, 207 /255]      ,
                0.078125 : [0 /255, 0 /255, 223 /255]    ,
                0.09375 : [0 /255, 0 /255, 239 /255]     ,
                0.109375 : [0 /255, 0 /255, 255 /255]    ,
                0.125 : [0 /255, 15 /255, 255 /255]      ,
                0.140625 : [0 /255, 31 /255, 255 /255]   ,
                0.15625 : [0 /255, 47 /255, 255 /255]    ,
                0.171875 : [0 /255, 63 /255, 255 /255]   ,
                0.1875 : [0 /255, 79 /255, 255 /255]     ,
                0.203125 : [0 /255, 95 /255, 255 /255]   ,
                0.21875 : [0 /255, 111 /255, 255 /255]   ,
                0.234375 : [0 /255, 127 /255, 255 /255]  ,
                0.25 : [0 /255, 143 /255, 255 /255]      ,
                0.265625 : [0 /255, 159 /255, 255 /255]  ,
                0.28125 : [0 /255, 175 /255, 255 /255]   ,
                0.296875 : [0 /255, 191 /255, 255 /255]  ,
                0.3125 : [0 /255, 207 /255, 255 /255]    ,
                0.328125 : [0 /255, 223 /255, 255 /255]  ,
                0.34375 : [0 /255, 239 /255, 255 /255]   ,
                0.359375 : [0 /255, 255 /255, 255 /255]  ,
                0.375 : [15 /255, 255 /255, 239 /255]    ,
                0.390625 : [31 /255, 255 /255, 223 /255] ,
                0.40625 : [47 /255, 255 /255, 207 /255]  ,
                0.421875 : [63 /255, 255 /255, 191 /255] ,
                0.4375 : [79 /255, 255 /255, 175 /255]   ,
                0.453125 : [95 /255, 255 /255, 159 /255] ,
                0.46875 : [111 /255, 255 /255, 143 /255] ,
                0.484375 : [127 /255, 255 /255, 127 /255],
                0.5 : [143 /255, 255 /255, 111 /255]     ,
                0.515625 : [159 /255, 255 /255, 95 /255] ,
                0.53125 : [175 /255, 255 /255, 79 /255]  ,
                0.546875 : [191 /255, 255 /255, 63 /255] ,
                0.5625 : [207 /255, 255 /255, 47 /255]   ,
                0.578125 : [223 /255, 255 /255, 31 /255] ,
                0.59375 : [239 /255, 255 /255, 15 /255]  ,
                0.609375 : [255 /255, 255 /255, 0 /255]  ,
                0.625 : [255 /255, 239 /255, 0 /255]     ,
                0.640625 : [255 /255, 223 /255, 0 /255]  ,
                0.65625 : [255 /255, 207 /255, 0 /255]   ,
                0.671875 : [255 /255, 191 /255, 0 /255]  ,
                0.6875 : [255 /255, 175 /255, 0 /255]    ,
                0.703125 : [255 /255, 159 /255, 0 /255]  ,
                0.71875 : [255 /255, 143 /255, 0 /255]   ,
                0.734375 : [255 /255, 127 /255, 0 /255]  ,
                0.75 : [255 /255, 111 /255, 0 /255]      ,
                0.765625 : [255 /255, 95 /255, 0 /255]   ,
                0.78125 : [255 /255, 79 /255, 0 /255]    ,
                0.796875 : [255 /255, 63 /255, 0 /255]   ,
                0.8125 : [255 /255, 47 /255, 0 /255]     ,
                0.828125 : [255 /255, 31 /255, 0 /255]   ,
                0.84375 : [255 /255, 15 /255, 0 /255]    ,
                0.859375 : [255 /255, 0 /255, 0 /255]    ,
                0.875 : [239 /255, 0 /255, 0 /255]       ,
                0.890625 : [223 /255, 0 /255, 0 /255]    ,
                0.90625 : [207 /255, 0 /255, 0 /255]     ,
                0.921875 : [191 /255, 0 /255, 0 /255]    ,
                0.9375 : [175 /255, 0 /255, 0 /255]      ,
                0.953125 : [159 /255, 0 /255, 0 /255]    ,
                0.96875 : [143 /255, 0 /255, 0 /255]     ,
                0.984375 : [127 /255, 0 /255, 0 /255]    ,
                1.1 : [127 /255, 0 /255, 0 /255]      ,
                },                  

}



def interpolate_color(color1, color2, fraction):
    return [ ( (x1 + fraction * (x2 - x1))) for (x1, x2) in zip(color1, color2)]

def get_color(valor, color_map):#, factor = 2):   
    
    
    #valor = (valor+1)/factor   
    
    cutoffs = list(color_map.keys())
    color_list = list(color_map.values())
    counter = 0 
    
    #print(valor)
    
    
    
    for cutoff, color in color_map.items():
        
        if valor <= cutoff :
            
            fraction = (valor -  cutoffs[counter-1]) / ( cutoff - cutoffs[counter-1])
            
            
            color1  = color_list[counter-1]
            color2  = color_list[counter]
            #
            #
            #
            #
            #
            color = interpolate_color(color1, color2, fraction)
            #color = [0,0,0]
            
            #color[0] = color[0] 
            #color[1] = color[1] 
            #color[2] = color[2] 
            #print(valor, cutoff, fraction,color )

            return color
        else:
            pass
            #print(valor, cutoff)#, fraction,color )
        counter +=1

class Canvas(Gtk.DrawingArea):
    """ Class doc """
    
    def __init__ (self, bg_color = [1,1,1], cmap =  'jet' ):
        """ Class initialiser """
        super().__init__( )
        #Gtk.DrawingArea.__init__(main)
        
        self.color_map = COLOR_MAPS[cmap]
        
        self.bg_color = None
        self.set_bg_color (bg_color[0], bg_color[1], bg_color[2])
        
        self.set_events(Gdk.EventMask.ALL_EVENTS_MASK)
        self.connect("draw", self.on_draw)
        self.connect("button_press_event", self.on_mouse_button_press)
        self.connect("motion-notify-event", self.on_motion)
    
    def set_cmap (self, cmap = 'jet'):
        """ Function doc """
        
    def set_bg_color (self, r = 1 , g = 1 , b = 1):
        """ Function doc """
        self.bg_color = [r, g, b]
    
    def _make_bg (self, cr):
        """ Function doc """
        cr.set_source_rgb( self.bg_color[0],self.bg_color[1], self.bg_color[2])
        cr.paint()
    
    def on_mouse_button_press(self, widget, event):
        #(x, y) = int(event.x), int(event.y)
        ##x, y = device.get_position(widget)
        #print("Mouse moved to:", x, y)
        (x, y) = int(event.x), int(event.y)
        
        x_on_plot = int((x-self.bx)/self.factor_x)
        y_on_plot = int((y-self.by)/self.factor_y)
        

        
        if int(x_on_plot) < 0 or int(x_on_plot) >= self.size_x :
            print('canvas')
            self.points = []
        elif int(y_on_plot) < 0 or int(y_on_plot) >= self.size_y :
            print('canvas')
            self.points = []
        else:
            
            
            self.points.append((x_on_plot, y_on_plot))
        
        print("Mouse clicker at:",  x, y, int(x_on_plot), int(y_on_plot), 
                                    (x-self.bx)/self.factor_x, #/(self.x_final-self.bx), 
                                    (y-self.by)/self.factor_y) #/(self.y_final-self.by) ) #, x - self.x_final , y - self.y_final , self.x_final ,  self.y_final ,   (x-self.bx) /self.factor_x ,    (y-self.by)/self.factor_y )
        print(self.points)
        self.queue_draw()

    
    def on_motion(self, widget, event):
        (x, y) = int(event.x), int(event.y)
        x_on_plot = int((x-self.bx)/self.factor_x)
        y_on_plot = int((y-self.by)/self.factor_y)
        #print("Mouse moved to:", x, y, int(x_on_plot), int(y_on_plot))


class ImagePlot(Canvas):
    """ Class doc """
    
    def __init__ (self, data = None):
        """ Class initialiser """
        super().__init__( )

        self.bx = 50 
        self.by = 50 
        self.points = []
        
        self.data = data
        
    
    def set_cmap (self, cmap = 'jet'):
        """ Function doc """
        
    
    def draw_color_bar (self, cr, res, gap = 20, label_spacing = 2, font_size = 12):
        """ Function doc """
         
       
        cr.set_line_width (1.0)
        
        for j in range(0, res*10, 1):
            
            #--------------------------------------------------------------------------------------------
            #color bar
            #--------------------------------------------------------------------------------------------
            z = j  / float(res*10)
            color = get_color(z, self.color_map)
            cr.set_source_rgb( color[0], color[1], color[2])
            
            j = res*10 -j-1
            
            x = round(self.size_x * self.factor_x)+self.bx + gap
            y = self.by + j*self.factor_y/10
            
            size_x  = round(20)+1 
            size_y  = round((self.factor_y/10)+1)
            
            
            
            cr.rectangle(
                         x,
                         y, 
                         size_x, 
                         size_y)
            cr.fill()
        
        
        
        factor = int( self.size_y/10 )
        
        if factor ==  0:
            factor = 1
        else:
            pass
        
        
        _min = np.min(self.data)
        _max = np.max(self.data)
        factor2 = (_max - _min)/10
        print(_min, _max)
        counter = 0
        for j in range(0, self.size_y+1, factor ):    
            
            
            
            j = self.size_y -j-1
            
            #z = ((j  / float(self.size_y)) - 0.5)*2
            #print(z)
            cr.set_source_rgb( 0, 0,0)
            cr.move_to (round(self.size_x * self.factor_x)+self.bx*1.2 +37      , 
                              self.by + j*self.factor_y + self.factor_y  )
            
            
            cr.line_to (round(self.size_x * self.factor_x)+self.bx*1.2 +30  , 
                              self.by + j*self.factor_y + self.factor_y   ) 
            #--------------------------------------------------------------------------------------------
            
            cr.stroke ()
            #--------------------------------------------------------------------------------------------

            #--------------------------------------------------------------------------------------------
            # texto em y
            cr.set_font_size(font_size)
            text = '{:<6.3f}'.format(_min + factor2*counter )
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)
            cr.move_to(round(self.size_x * self.factor_x)+self.bx*1.2 + 45 ,  
                             self.by + j*self.factor_y+ self.factor_y + height/2)
            cr.show_text(text)
            #--------------------------------------------------------------------------------------------
            counter +=1








        cr.set_line_width (2.0)
        cr.rectangle(round(self.size_x * self.factor_x)+self.bx +gap,
                     self.by   , 
                     round(20+1), 
                     round(self.size_y * self.factor_y))
        cr.stroke ()


    def draw_image (self, cr):
        """ Function doc """
        n = 0 
        for i in range(0,self.size_x, ):
            for j in range(0,self.size_y):
                
                z = self.norm_data[i][j]

                color = get_color(z, self.color_map)

                cr.set_source_rgb( color[0], color[1], color[2]   )
                
                cr.rectangle(self.bx+i*self.factor_x -1,
                             self.by+j*self.factor_y -1, 
                             round(self.factor_x)+1, 
                             round(self.factor_y)+1)
                cr.fill()


    def draw_image_box (self, cr, line_width = 1, color = [0,0,0]):
        cr.set_line_width (line_width)
        cr.set_source_rgb(color[0], color[1], color[2])
        cr.rectangle(self.bx,
                     self.by, 
                     round(self.size_x * self.factor_x), 
                     round(self.size_y * self.factor_y))
        cr.stroke ()
        
    
    
    def draw_scale (self, cr, font_size = 12, x_major_ticks = 3 , y_major_ticks = 3 ):
        """ Function doc """        
        
        round(self.size_x/x_major_ticks)
        
        for i in range(0,self.size_x, x_major_ticks):
            #--------------------------------------------------------------------------------------------
            # marcacoes em x
            #print(self.bx+ i*self.factor_x + (self.factor_x)/2.0 , self.bx+ i*self.factor_x + (self.factor_x)/2.0)
            cr.move_to (self.bx+ i*self.factor_x + (self.factor_x)/2.0      , self.by )
            cr.line_to (self.bx+ i*self.factor_x + (self.factor_x)/2.0      , self.by -10 ) 
            #--------------------------------------------------------------------------------------------
            cr.stroke ()
            
            
            #--------------------------------------------------------------------------------------------
            cr.set_source_rgb(0, 0, 0)
            cr.set_font_size(font_size)
    
            text = str(i)
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)
            #print(x_bearing, y_bearing, width, height, x_advance, y_advance)
            #cr.move_to(self.bx + i*self.factor_x - width/2.0    ,  self.by -15 )   
            cr.move_to(self.bx + i*self.factor_x - x_advance/2.0 +(self.factor_x)/2   ,  self.by -15 )   
            cr.show_text(text)
            #--------------------------------------------------------------------------------------------
        
        for j in range(0,self.size_y, y_major_ticks):
            #--------------------------------------------------------------------------------------------
            # marcacoes em y
            cr.move_to (self.bx      , self.by + j*self.factor_y + (self.factor_y)/2.0 )
            cr.line_to (self.bx -10  , self.by + j*self.factor_y + (self.factor_y)/2.0 ) 
            #--------------------------------------------------------------------------------------------
            cr.stroke ()
            
            
            #--------------------------------------------------------------------------------------------
            # texto em y
            cr.set_source_rgb(0, 0, 0)
            cr.set_font_size(font_size)
    
            text = str(j)
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)

            cr.move_to(self.bx  -15 -width     , self.by + j*self.factor_y + height/2  + (self.factor_y)/2)
            cr.show_text(str(j))
            #--------------------------------------------------------------------------------------------
        
    def draw_dots (self, cr, color = [0,0,0]):
        """ Function doc """
        cr.set_source_rgb( color[0], color[1], color[2])
        
        for x, y in self.points:
            cr.arc(((x+0.5)*self.factor_x)+ self.bx, ((y+0.5)*self.factor_y)+self.by , 5, 0, 2 * 3.14)
            cr.fill()


    def draw_lines (self, cr, line_width = 2, color = [0,0,0]):
        """ Function doc """
        cr.set_source_rgb( color[0], color[1], color[2])
        cr.set_line_width (line_width)
        for x, y in self.points:
            #cr.line_to(x , y)
            cr.line_to(((x+0.5)*self.factor_x)+ self.bx, ((y+0.5)*self.factor_y)+self.by)
        cr.stroke ()


    def on_draw(self, widget, cr, ):
        self.cr =  cr
        self.width = widget.get_allocated_width()   #widget dimentions
        self.height = widget.get_allocated_height() #widget dimentions   
        self._make_bg(cr)

        
        self.size_x = len(self.data)
        self.size_y = len(self.data[0])

        
        self.x_plot_edge = self.width  -(self.bx+70)
        self.y_plot_edge = self.height -(self.by)
        
        '''(self.x_plot_edge - self.bx) = number of pixels the plot has at x 
        We divide by the dimension of the data matrix to know how much each 
        rectangle of the plot will be in length'''
        
        self.factor_x = (self.x_plot_edge - self.bx)/self.size_x 
        self.factor_y = (self.y_plot_edge - self.by)/self.size_y
        
        self.x_plot_edge = self.size_x*self.factor_x-self.factor_x
        self.y_plot_edge = self.size_y*self.factor_y-self.factor_y
        
        
        
        
        self.norm_data = self.data - np.min(self.data)
        
        
        self.data_min = np.min(self.norm_data)
        self.data_max = np.max(self.norm_data)
        delta = (self.data_max - self.data_min)
        #print(self.data_min, self.data_max )
        self.norm_data = self.norm_data/delta

        #print (self.norm_data)

        #print(np.min(self.norm_data), np.max(self.norm_data) )
        
        self.draw_color_bar (cr, res = self.size_y)

        self.draw_image (cr)
        
        self.draw_image_box (cr, line_width = 1, color = [0,0,0])

        self.draw_scale(cr)

        self.draw_dots (cr)

        self.draw_lines (cr, line_width = 2, color = [0,0,0])



class LinePlot(Canvas):
    """ Class doc """
    
    def __init__ (self, data = None):
        """ Class initialiser """
        super().__init__( )

        self.bx = 50 
        self.by = 50 
        self.points = []
        
        self.data = data





        
    


def main(args):
    
    window = Gtk.Window(title="Cairo and GTK Example")
    window.set_default_size(500, 400)    
    
    plot = ImagePlot()


    #from numpy import ma
    #data2d = np.random.random((50, 100))
    t = np.linspace(0, 2 * np.pi, 40)
    data2d = np.sin(t)[:, np.newaxis] * np.cos(t)[np.newaxis, :]
        
    #else:
        #self.data = data
    plot = ImagePlot(data2d)
    
    window.add(plot)
    
    
    
    
    
    window.connect("delete-event", Gtk.main_quit)
    window.show_all()
    Gtk.main()







    
    
    return 0

if __name__ == '__main__':
    
    sys.exit(main(sys.argv))
