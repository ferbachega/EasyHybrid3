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

# --- imports entre modulos adicionados na refatoracao ---
from util.easyplot.canvas import Canvas
from util.easyplot.color_utils import get_color, bilinear_interpolation, expand_array_size

class ImagePlot(Canvas):
    """ Class doc """
    
    def __init__ (self, data = None):
        """ Class initialiser 
        
        Now the class can receive an energy matrix "Z" with 
        infinite value "inf". To do this it was necessary 
        to modify the get_color method.
        
        """
        super().__init__( )

        self.bx = 50 
        self.by = 50 
        self.points = []
        self.extra_points = []
        
        self.data = data
        self.dataRC1 = None
        self.dataRC2 = None   
        
        self.cmap = 'jet'
        
        self.norm_data = None
        self.data_min  = None
        self.data_max  = None
        self.connect("motion-notify-event", self.on_motion)
        
        self.RC_label = None
        self.label_mode = 0
        self.is_discrete = True
        
        self.selected_dot = None #it is the i and j coordinates at the energey matrix
        self.sel_dot_rgb  = [1,0,0]

        # [NOVO] Cor unica para eixos, marcas de escala, rotulos (x/y e da
        # barra de cores) e a moldura ao redor da matriz -- antes cada um
        # desses tinha [0,0,0] fixo espalhado em lugares diferentes
        # (draw_scale, draw_color_bar, a chamada de draw_image_box em
        # on_draw), sem nenhuma forma de mudar tudo de uma vez so'.
        self.axis_color = [0, 0, 0]

        # --- contour lines (isolines) -------------------------------------
        self.show_contours       = False
        self.num_contour_levels  = 8
        self.contour_color       = [0, 0, 0]
        self.contour_line_width  = 1.0
        self.contour_label_font_size = 11
        self.contour_label_format    = '{:.1f}'
        self._contour_cache_key  = None
        self._contour_cache      = { }
        
        
    def define_datanorm (self):
        """ Function doc 
        When you have a 2D umbrella sampling data, 
        some points have energy estimated with inf (infinity).
        
        In python, "float('inf') is allowed - accepted as a float.
        
        When the value in the matrix is ​​given by inf, IMAGEPLOT 
        should draw a white square (showing that there is no energy 
        calculated at the point)
        
        In this case it is necessary to temporarily replace the inf 
        values ​​with 0 (so that it is possible to extract the maximum 
        and minimum values ​​of the matrix)
        """
        
        
        temp_data = self.data.copy()
        a = len(temp_data)
        b = len(temp_data[0])
        inf = float('inf')
        
        #replacing the inf values ​​of the matrix  
        inf_list = []
        for i in range(a):
            for j in range(b):
                if temp_data[i][j] == inf:
                    temp_data[i][j] = 0
                    inf_list.append([i, j])
        #print(temp_data)
        self.norm_data = temp_data - np.min(temp_data)
        #self.norm_data = self.data - np.min(self.data)

        self.data_min = np.min(self.norm_data)
        self.data_max = np.max(self.norm_data)
        delta = (self.data_max - self.data_min)
        self.norm_data = self.norm_data/delta
        
        #resetting the matrix inf values
        for inf_item in inf_list:
            i = inf_item[0]
            j = inf_item[1]
            self.norm_data[i][j] = float('inf')
        
        return self.data_min, self.data_max, delta , self.norm_data
        
        
    def set_threshold_color (self, _min = 0, _max = None, cmap = 'jet'):
        """ Function doc """
        self.color_map = COLOR_MAPS[self.cmap]
        fraction = _max / self.data_max #number between 0 and 1
        
        new_colormap = {}
        
        for i in range(101):
            value = i/100
            
            color = get_color(value, self.color_map)
            new_colormap[fraction*value] = color
        
        last_color = get_color(1, self.color_map)
        new_colormap[1.1] = last_color
        #pprint(self.color_map)
        #pprint( new_colormap)
        self.color_map = new_colormap
    
        
    def set_cmap (self, cmap = 'jet'):
        """ Function doc """
        self.color_map = COLOR_MAPS[cmap]
    
    
    def set_label_mode (self, mode = 0):
        """ Function doc """
        self.label_mode = mode
        
    
    def draw_color_bar (self, cr, res, gap = 20, label_spacing = 2, font_size = 12, num_labels = 10):
        """ Function doc """
         
       
        cr.set_line_width (1.0)
        #print('\n\n',res,'\n\n')
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
        
        
        
        factor = self.size_y/num_labels
        if factor ==  0:
            factor = 1
        else:
            pass
        
        
        
        cr.set_source_rgb(*self.axis_color)
        _min = self.data_min
        _max = self.data_max
        factor2 = (_max - _min)/(num_labels-1)
        #print(_min, _max, factor, self.size_y, self.size_y/factor, self.factor_y)
        counter = 0
        
        '''
           The position of the labels depends on the number of elements 
        in y of the data matrix (self.size_y). It depends on the factor 
        in y (self.factor_y), which in turn takes into account the 
        number of squares in the matrix and the size of the window.
        
        
        '''
        for i in range(0, num_labels):    
            
            j = (num_labels-1)-i
            
            x1 = round(self.size_x * self.factor_x)+self.bx*1.2 +37
            y1 = self.by + (j*self.factor_y*(self.size_y/(num_labels-1)))
            
            cr.move_to (x1,y1)
            
            x2 = round(self.size_x * self.factor_x)+self.bx*1.2 +30  
            cr.line_to (x2,y1) 
            #--------------------------------------------------------------------------------------------
            cr.stroke ()
            
            # . text
            cr.set_font_size(font_size)
            text = '{:<6.3f}'.format(i*factor2)
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)
            
            x3 = round(self.size_x * self.factor_x)+self.bx*1.2 + 45 
            cr.move_to(x3,y1 + height/2)
            cr.show_text(text)
            #--------------------------------------------------------------------------------------------
            counter +=1
        
        
        
        '''
        for j in range(0, self.size_y+1, factor ):    

            j = self.size_y -j-1
            
            #z = ((j  / float(self.size_y)) - 0.5)*2
            #print(z)
            cr.set_source_rgb(0,0,0)
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
            #print(_min, factor2, counter, _min + factor2*counter)
            text = '{:<6.3f}'.format(_min + factor2*counter )
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)
            cr.move_to(round(self.size_x * self.factor_x)+self.bx*1.2 + 45 ,  
                             self.by + j*self.factor_y+ self.factor_y + height/2)
            cr.show_text(text)
            #--------------------------------------------------------------------------------------------
            counter +=1
        '''







        cr.set_line_width (2.0)
        cr.rectangle(round(self.size_x * self.factor_x)+self.bx +gap,
                     self.by   , 
                     round(20+1), 
                     round(self.size_y * self.factor_y))
        cr.stroke ()


    def draw_discrete_square (self, i, j, cr):
        """ Function doc """
        #z is the nergey value
        z = self.norm_data[i][j]
        color = get_color(z, self.color_map)
        
        
        if color:
            cr.set_source_rgb( color[0], color[1], color[2]   )
        else:
            #.Associates the color white when color is null 
            #.(due to matrix points with value inf)
            cr.set_source_rgb( 1, 1, 1   )
        
        cr.rectangle(self.bx+i*self.factor_x -1,
                     self.by+j*self.factor_y -1, 
                     round(self.factor_x)+1, 
                     round(self.factor_y)+1)
        cr.fill()
    
    
    def draw_smooth_square (self, cr = None, square = None, colors = None):
        """ Function doc 
        
        
        """
        x0 , y0 = square[0] # is tuple of two coordinates
        x05, y0 = square[1] # is tuple of two coordinates
        x05,y05 = square[2] # is tuple of two coordinates
        x0 ,y05 = square[3] # is tuple of two coordinates
        
        pattern = cairo.MeshPattern()
        # Começa a definição do patch de padrão de malha
        pattern.begin_patch()
        pattern.move_to(x0 , y0)
        pattern.line_to(x05, y0)
        pattern.line_to(x05,y05)
        pattern.line_to(x0 ,y05)
        
        color1 = colors[0]
        color2 = colors[1]
        color3 = colors[2]
        color4 = colors[3]
        
        ### Define os vértices do patch
        ##pattern.move_to(i,     j)
        ##pattern.line_to(i+100, j)
        ##pattern.line_to(i+100, j+100)
        ##pattern.line_to(i    , j+100)

        #color1 = get_color(z1, self.color_map)
        #color2 = get_color(z2, self.color_map)
        #color3 = get_color(z3, self.color_map)
        #color4 = get_color(z4, self.color_map)
        #
        pattern.set_corner_color_rgb(0, color1[0], color1[1], color1[2]) 
        pattern.set_corner_color_rgb(1, color2[0], color2[1], color2[2]) 
        pattern.set_corner_color_rgb(2, color3[0], color3[1], color3[2]) 
        pattern.set_corner_color_rgb(3, color4[0], color4[1], color4[2]) 
        
        # Finaliza a definição do patch de padrão de malha
        #cr.scale(1, 1)
        pattern.end_patch()
        cr.set_source(pattern)
        cr.paint()          

    
    def draw_inter_square (self, surface, x, y, size, c1, c2, c3, c4, cr):
        """ Function doc """
        #self.ctx = cr
        #self.x, self.y, self.size = x, y, size
        #self.c1, self.c2, self.c3, self.c4 = c1, c2, c3, c4

        for i in range(size):
            for j in range(size):
                tx = j / (size - 1)
                ty = i / (size - 1)
                color = bilinear_interpolation(c1, c2, c3, c4, tx, ty)

                cr.set_source_rgba(*[comp / 255.0 for comp in color])
                cr.rectangle(x + j, y + i, 1, 1)
                cr.fill()


    def draw_image (self, cr):
        """ Function doc """
        if self.is_discrete == True:
            #discrete square
            #'''
            n = 0 
            for i in range(0,self.size_x, ):
                for j in range(0,self.size_y):
                    self.draw_discrete_square (i, j, cr)


        else:
            
            data = expand_array_size(self.norm_data)
        
            for i in range(1, self.size_x+1):
                for j in range(1,  self.size_y+1):
                    
                    x0, y0 = self.bx+(i-1)*self.factor_x -1   , self.by+(j-1)*self.factor_y -1
                    x1, y1 = self.bx+(i )*self.factor_x -1 , self.by+ (j )*self.factor_y -1
                    
                    x05   = (x1 + x0)/2
                    y05   = (y1 + y0)/2
                    '''
                    dprint(
                    x0 , y0   ,
                    x0 , y05 ,
                    x0 , y1   ,

                    x05 , y0  ,
                    x05 , y05,
                    x05 , y1  ,

                    x1 , y0   ,
                    x1 , y05 ,
                    x1 , y1   ,
                    )
                    '''
                    # calculating energies / avarege 
                    E_x05_y05 = data[i][j]
                    C_x05_y05 = get_color(E_x05_y05, self.color_map)
                    
                    E_x0_y0 = (data[i-1][j ] + data[i-1][j-1] + data[i ][j-1] +  data[i][j])/4
                    C_x0_y0 = get_color(E_x0_y0, self.color_map)
                    
                    E_x1_y0 = (data[i][j-1] + data[i+1][j] + data[i+1][j-1] +  data[i][j])/4
                    C_x1_y0 = get_color(E_x1_y0, self.color_map)
                    
                    E_x0_y1 = (data[i-1][j+1]+ data[i-1][j]+ data[i][j+1] +  data[i][j])/4
                    C_x0_y1 = get_color(E_x0_y1, self.color_map)
                    
                    E_x05_y1 = (data[i][j+1] +  data[i][j])/2
                    C_x05_y1 = get_color(E_x05_y1, self.color_map)
                    
                    E_x05_y0 = (data[i][j-1] +  data[i][j])/2
                    C_x05_y0 = get_color(E_x05_y0, self.color_map)
                    
                    E_x0_y05 = (data[i-1][j] +  data[i][j])/2
                    C_x0_y05 = get_color(E_x0_y05, self.color_map)
                    
                    
                    E_x1_y05 = (data[i+1][j] +  data[i][j])/2
                    C_x1_y05 = get_color(E_x1_y05, self.color_map)
                    
                    E_x1_y1 = (data[i+1][j+1] +data[i][j+1] +data[i+1][j] +  data[i ][j ])/4
                    C_x1_y1 = get_color(E_x1_y1, self.color_map)
                    #C_x1_y1 = [0 ,0 ,0]
                    
                    
                    
                    '''each square is actually 4 squares that will be drawn next.'''
                    #---------------------------------------------------
                    square = [(x0 , y0),   # 
                              (x05, y0),   #  o------o------o
                              (x05,y05),   #  |XXXXXX|      |
                              (x0 ,y05)]   #  |XXXXXX|      |
                                           #  o------o------o
                    colors = [C_x0_y0  ,   #  |      |      |
                              C_x05_y0 ,   #  |      |      |
                              C_x05_y05,   #  o------o------o
                              C_x0_y05 ]   #
                    
                    self.draw_smooth_square(cr = cr, square = square, colors = colors)
                    #---------------------------------------------------

            
                    #---------------------------------------------------
                    square = [(x05 , y0),   # 
                              (x1  , y0),   #  o------o------o
                              (x1  , y05),  #  |      |      |
                              (x05 , y05)]  #  |      |      |
                                            #  o------o------o
                    colors = [C_x05_y0,     #  |XXXXXX|      |
                              C_x1_y0,      #  |XXXXXX|      |
                              C_x1_y05,     #  o------o------o
                              C_x05_y05]    #

                    self.draw_smooth_square(cr = cr, square = square, colors = colors)
                    #---------------------------------------------------
                    
                    
                    #---------------------------------------------------
                    square = [(x0  , y05),  # 
                              (x05 ,y05),   #  o------o------o
                              (x05 , y1 ),  #  |      |XXXXXX|
                              (x0  , y1)]   #  |      |XXXXXX|
                                            #  o------o------o
                    colors = [C_x0_y05,     #  |      |      |
                              C_x05_y05,    #  |      |      |
                              C_x05_y1,     #  o------o------o
                              C_x0_y1]      #

                    self.draw_smooth_square(cr = cr, square = square, colors = colors)
                    #---------------------------------------------------

                    
                    #---------------------------------------------------
                    square = [(x05 , y05),  # 
                              (x1  ,y05),   #  o------o------o
                              (x1  , y1 ),  #  |      |      |
                              (x05  , y1)]  #  |      |      |
                                            #  o------o------o
                    colors = [C_x05_y05,    #  |      |XXXXXX|
                              C_x1_y05,     #  |      |XXXXXX|
                              C_x1_y1,      #  o------o------o
                              C_x05_y1]     #

                    self.draw_smooth_square(cr = cr, square = square, colors = colors)
                    #---------------------------------------------------
                             
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
            # [EN] BUG FIX: antes esta linha de marcacao nao definia sua
            # propria cor -- so' "funcionava" (saia preta) porque
            # draw_image_box() (chamada logo antes, em on_draw) tambem
            # deixava a cor preta fixa no cr. Agora que axis_color e'
            # configuravel, cada desenho precisa definir sua propria cor
            # explicitamente, sem depender de estado deixado por quem
            # desenhou antes.
            cr.set_source_rgb(*self.axis_color)
            cr.move_to (self.bx+ i*self.factor_x + (self.factor_x)/2.0      , self.by )
            cr.line_to (self.bx+ i*self.factor_x + (self.factor_x)/2.0      , self.by -10 ) 
            #--------------------------------------------------------------------------------------------
            cr.stroke ()
            
            
            #--------------------------------------------------------------------------------------------
            cr.set_source_rgb(*self.axis_color)
            cr.set_font_size(font_size)
            
            if self.label_mode == 0:
                text = str(i)
            elif self.label_mode == 1:
                '''
                   The label of reaction coordinate distance list is wrong  
                   self.dataRC2[i][0] - means frames on the vertical (index i)
                '''
                text = '{:3.2f}'.format(self.dataRC2[i][0])
            else:
                text = ''
                pass
            
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)
            #print(x_bearing, y_bearing, width, height, x_advance, y_advance)
            #cr.move_to(self.bx + i*self.factor_x - width/2.0    ,  self.by -15 )   
            cr.move_to(self.bx + i*self.factor_x - x_advance/2.0 +(self.factor_x)/2   ,  self.by -15 )   
            cr.show_text(text)
            #--------------------------------------------------------------------------------------------
        
        for j in range(0,self.size_y, y_major_ticks):
            #--------------------------------------------------------------------------------------------
            # marcacoes em y
            cr.set_source_rgb(*self.axis_color)
            cr.move_to (self.bx      , self.by + j*self.factor_y + (self.factor_y)/2.0 )
            cr.line_to (self.bx -10  , self.by + j*self.factor_y + (self.factor_y)/2.0 ) 
            #--------------------------------------------------------------------------------------------
            cr.stroke ()
            
            
            #--------------------------------------------------------------------------------------------
            # texto em y
            cr.set_source_rgb(*self.axis_color)
            cr.set_font_size(font_size)

            if self.label_mode == 0:
                text = str(j)
            elif self.label_mode == 1:
                '''
                The label of reaction coordinate distance list is wrong  
                self.dataRC1[j][0] - means frames on the vertical (index j)
                '''
                text = '{:3.2f}'.format(self.dataRC1[0][j])
            else:
                text = ''
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents(text)

            cr.move_to(self.bx  -15 -width     , self.by + j*self.factor_y + height/2  + (self.factor_y)/2)
            cr.show_text(text)
            #--------------------------------------------------------------------------------------------
        
    def draw_dots (self, cr, points, color = [0,0,0], fill = True):
        """ Function doc """
        cr.stroke ()
        cr.set_source_rgb( color[0], color[1], color[2])
        #print(self.points)
        for x,y in points:
            #else:
            cr.set_source_rgb( color[0], color[1], color[2])
            cr.arc(((x+0.5)*self.factor_x)+ self.bx, ((y+0.5)*self.factor_y)+self.by , 5, 0, 2 * 3.14)
            #cr.stroke ()
            if fill:
                cr.fill()
            else:
                pass
                cr.stroke ()
        # draws the red dot at the energey matrix
        if len(points) and self.selected_dot:
            x = self.selected_dot[0]
            y = self.selected_dot[1]
            cr.set_source_rgb( self.sel_dot_rgb[0], self.sel_dot_rgb[1], self.sel_dot_rgb[2])
            cr.arc(((x+0.5)*self.factor_x)+ self.bx, ((y+0.5)*self.factor_y)+self.by , 5, 0, 2 * 3.14)
            
            if fill:
                cr.fill()
            else:
                pass
                cr.stroke ()
    
    def draw_lines (self, cr, line_width = 2, color = [0,0,0]):
        """ Function doc """
        cr.set_source_rgb( color[0], color[1], color[2])
        cr.set_line_width (line_width)
        for x, y in self.points:
            #cr.line_to(x , y)
            cr.line_to(((x+0.5)*self.factor_x)+ self.bx, ((y+0.5)*self.factor_y)+self.by)
        cr.stroke ()


    # ─────────────────────────────────────────────────────────────────────
    #  C O N T O U R   L I N E S   ( I S O L I N E S )
    #
    #  Simple marching-squares implementation that traces energy contour
    #  lines directly over self.data (the real, un-normalized values), so
    #  requested levels are in the same physical units (kJ/mol) shown
    #  elsewhere in the plot (color bar, hover readout).
    # ─────────────────────────────────────────────────────────────────────
    @staticmethod
    def _interp (p1, p2, v1, v2, level):
        """ Linear interpolation of the point along edge p1-p2 where the
        data value crosses "level". """
        if v1 == v2:
            t = 0.5
        else:
            t = (level - v1) / (v2 - v1)
        x = p1[0] + t * (p2[0] - p1[0])
        y = p1[1] + t * (p2[1] - p1[1])
        return (x, y)

    @classmethod
    def _cell_contour_segments (cls, corners, values, level):
        """
        corners: 4 (x, y) grid points in cell order [bottom-left,
                 bottom-right, top-right, top-left].
        values:  the 4 data values at those corners, same order.

        Returns a list of 0, 1 or 2 line segments ((x0,y0),(x1,y1)) where
        the contour at "level" crosses this cell. Uses a generic
        edge-crossing test (equivalent to marching squares) instead of a
        hardcoded 16-case lookup table, including the standard "average
        value" rule to disambiguate saddle points (where the contour
        crosses all 4 edges of the cell and there are two valid ways to
        connect them).
        """
        edges = ( (0, 1), (1, 2), (2, 3), (3, 0) )   # bottom, right, top, left
        crossings = [ ]
        for a, b in edges:
            va, vb = values[a], values[b]
            if (va >= level) != (vb >= level):
                crossings.append ( cls._interp ( corners[a], corners[b], va, vb, level ) )

        if len ( crossings ) == 2:
            return [ ( crossings[0], crossings[1] ) ]

        if len ( crossings ) == 4:
            # . Saddle point -- disambiguate using the cell's average value.
            bottom, right, top, left = crossings
            average = sum ( values ) / 4.0
            if average >= level:
                return [ (left, bottom), (right, top) ]
            return [ (left, top), (bottom, right) ]

        return [ ]

    def _compute_contour_segments (self, levels):
        """ Runs marching squares over the whole data grid for every
        requested level. Returns { level: [ (p0, p1), ... ] } with points
        in fractional GRID coordinates (not yet converted to pixels).
        Cells touching an undefined value (inf/-inf/NaN -- e.g. missing
        points in an umbrella-sampling grid) are skipped, so contours
        never cross through regions with no data. """
        data    = self.data
        size_x  = len ( data )
        size_y  = len ( data[0] )

        segments = { level: [ ] for level in levels }

        for i in range ( size_x - 1 ):
            for j in range ( size_y - 1 ):
                corner_values = [ data[i][j], data[i+1][j], data[i+1][j+1], data[i][j+1] ]
                if any ( ( v in ( float('inf'), float('-inf') ) ) or ( v != v ) for v in corner_values ):
                    continue

                corner_points = [ (i, j), (i+1, j), (i+1, j+1), (i, j+1) ]

                for level in levels:
                    segments[level].extend (
                        self._cell_contour_segments ( corner_points, corner_values, level ) )

        return segments

    def get_contour_levels (self, num_levels = None):
        """ Evenly-spaced levels strictly between the data's min and max
        (excluding the extremes themselves, which would otherwise produce
        degenerate contours right at the edge of the surface). """
        if num_levels is None:
            num_levels = self.num_contour_levels

        finite_values = [ v for row in self.data for v in row
                          if v not in ( float('inf'), float('-inf') ) and v == v ]
        if len ( finite_values ) < 2:
            return [ ]

        vmin, vmax = min ( finite_values ), max ( finite_values )
        if vmax <= vmin:
            return [ ]

        return list ( np.linspace ( vmin, vmax, num_levels + 2 )[1:-1] )

    def _ensure_contour_cache (self, levels):
        """ (Re)computes self._contour_cache only when the data or the
        requested levels actually changed since the last call. """
        cache_key = ( id ( self.data ), tuple ( levels ) )
        if cache_key != self._contour_cache_key:
            self._contour_cache     = self._compute_contour_segments ( levels )
            self._contour_cache_key = cache_key

    def draw_contours (self, cr, levels = None, label_gaps = None):
        """ Draws energy contour lines (isolines) on top of the 2D PES
        surface. "levels" defaults to get_contour_levels() (evenly spaced
        between the data's min and max). "label_gaps" is an optional set
        of id(segment) values (see compute_contour_labels()) to skip --
        used to leave a small break in the line exactly where an inline
        label sits, so the label text doesn't overlap the line under it.
        Segment positions are cached and only recomputed when the data or
        the requested levels change, so redraws triggered by e.g. mouse
        movement stay cheap. """
        if not self.data:
            return

        if levels is None:
            levels = self.get_contour_levels ( )
        if not levels:
            return

        self._ensure_contour_cache ( levels )
        label_gaps = label_gaps or set ( )

        cr.save ( )
        cr.set_source_rgb ( *self.contour_color )
        cr.set_line_width ( self.contour_line_width )

        for level in levels:
            for segment in self._contour_cache.get ( level, [ ] ):
                if id ( segment ) in label_gaps:
                    continue
                p0, p1 = segment
                px0 = self.bx + (p0[0] + 0.5) * self.factor_x
                py0 = self.by + (p0[1] + 0.5) * self.factor_y
                px1 = self.bx + (p1[0] + 0.5) * self.factor_x
                py1 = self.by + (p1[1] + 0.5) * self.factor_y
                cr.move_to ( px0, py0 )
                cr.line_to ( px1, py1 )

        cr.stroke ( )
        cr.restore ( )

    @staticmethod
    def _trace_branches (segments, precision = 6):
        """
        Groups a flat list of contour-line segments into separate
        connected branches (visually distinct curves/loops), by chaining
        segments that share an endpoint. Two endpoints are considered
        "the same point" if they round to the same coordinates at
        "precision" decimal places (marching squares produces exactly
        shared endpoints between neighbouring cells; this is just a
        safety margin against floating-point noise).

        Returns a list of branches, each a list of the very same segment
        objects passed in (so the caller can later match one of them by
        id()/identity -- e.g. to leave a gap in the line under a label).
        """
        def point_key (pt):
            return ( round ( pt[0], precision ), round ( pt[1], precision ) )

        parent = { }

        def find (x):
            parent.setdefault ( x, x )
            root = x
            while parent[root] != root:
                root = parent[root]
            while parent[x] != root:
                parent[x], x = root, parent[x]
            return root

        def union (a, b):
            ra, rb = find ( a ), find ( b )
            if ra != rb:
                parent[ra] = rb

        for (p0, p1) in segments:
            union ( point_key ( p0 ), point_key ( p1 ) )

        branches = { }
        for segment in segments:
            p0, _p1 = segment
            root = find ( point_key ( p0 ) )
            branches.setdefault ( root, [ ] ).append ( segment )

        return list ( branches.values ( ) )

    def compute_contour_labels (self, levels = None, min_spacing_px = 45, min_branch_segments = 2):
        """
        Decides where to place inline contour labels: one per visually
        distinct branch/loop of each requested level (so a level that
        shows up as several separate curves -- e.g. around different
        basins -- gets more than a single label), skipping branches too
        short to fit a label legibly, and finally dropping any candidate
        that would land closer than "min_spacing_px" screen pixels to a
        label already kept -- so labels never end up cluttered/
        overlapping no matter how many levels/branches are requested.

        Returns a list of dicts with keys "level", "segment" (the exact
        (p0, p1) grid-coordinate segment the label sits on, so the line
        can be broken there -- see draw_contours()'s "label_gaps") and
        "pixel" (the label's screen position).
        """
        if not self.data:
            return [ ]
        if levels is None:
            levels = self.get_contour_levels ( )
        if not levels:
            return [ ]

        self._ensure_contour_cache ( levels )

        candidates = [ ]
        for level in levels:
            segments = self._contour_cache.get ( level, [ ] )
            if not segments:
                continue
            for branch in self._trace_branches ( segments ):
                if len ( branch ) < min_branch_segments:
                    continue
                mid_segment = branch[ len ( branch ) // 2 ]
                p0, p1 = mid_segment
                gx = (p0[0] + p1[0]) / 2.0
                gy = (p0[1] + p1[1]) / 2.0
                px = self.bx + (gx + 0.5) * self.factor_x
                py = self.by + (gy + 0.5) * self.factor_y
                candidates.append ( { 'level' : level, 'segment' : mid_segment, 'pixel' : (px, py) } )

        accepted = [ ]
        for candidate in candidates:
            px, py = candidate['pixel']
            too_close = False
            for other in accepted:
                ox, oy = other['pixel']
                if ( (px - ox) ** 2 + (py - oy) ** 2 ) ** 0.5 < min_spacing_px:
                    too_close = True
                    break
            if not too_close:
                accepted.append ( candidate )

        return accepted

    def draw_contour_labels (self, cr, labels):
        """ Draws the inline value label (e.g. "42.3") for each entry
        returned by compute_contour_labels(). The background is fully
        transparent -- meant to be used together with passing the same
        labels' segments as draw_contours()'s "label_gaps", which leaves
        a small break in the line so the text never overlaps it. """
        if not labels:
            return

        cr.save ( )
        cr.select_font_face ( "Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL )
        cr.set_font_size ( self.contour_label_font_size )
        cr.set_source_rgb ( *self.contour_color )

        for label in labels:
            px, py = label['pixel']
            text = self.contour_label_format.format ( label['level'] )
            x_bearing, y_bearing, width, height, x_advance, y_advance = cr.text_extents ( text )

            cr.move_to ( px - width/2.0 - x_bearing, py - height/2.0 - y_bearing )
            cr.show_text ( text )

        cr.restore ( )

    def on_draw(self, widget, cr, export_scale=1):
        """ [EN] export_scale: usado SO' pela exportacao em alta
        resolucao (ver util/easyplot/export_utils.export_plot_to_png()) --
        1 (default) para o desenho normal na tela.

        Por que isso precisa entrar aqui, e nao so' um cr.scale() de fora:
        este metodo desenha tudo numa superficie RASTER interna
        (self.temp_surface, criada logo abaixo) do tamanho exato de
        self.width/self.height, e so' no final "cola" (blit) essa
        superficie no cr recebido -- um cr.scale() aplicado de fora,
        sem tocar aqui dentro, so' ampliaria essa superficie ja' pronta
        (borrado, tipo dar zoom numa foto pequena), sem redesenhar nada
        em resolucao maior de verdade. Aplicando export_scale AQUI (na
        criacao de temp_surface e no proprio temp_cr), a grade de dados,
        as fontes e as margens sao desenhadas de novo, todas maiores
        proporcionalmente juntas -- exatamente como ficariam numa tela
        export_scale vezes maior, so' que nitido. """
        self.cr =  cr
        self.width = widget.get_allocated_width()   #widget dimentions
        self.height = widget.get_allocated_height() #widget dimentions   
        self._make_bg(cr)

        if self.data:
            pass
        else:
            return None
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
        
        if self.data_max is None:
            self.define_datanorm()
        
        '''
        self.norm_data = self.data - np.min(self.data)
        
        
        self.data_min = np.min(self.norm_data)
        self.data_max = np.max(self.norm_data)
        #self.data_max = 100
        delta = (self.data_max - self.data_min)
        self.norm_data = self.norm_data/delta
        
        self.set_threshold_color ( _min = 0, _max = 200)
        '''
        
        # Criar uma superfície de imagem temporária
        self.temp_surface = cairo.ImageSurface(cairo.FORMAT_ARGB32,
                                                int(self.width * export_scale),
                                                int(self.height * export_scale))
        self.temp_cr = cairo.Context(self.temp_surface)
        if export_scale != 1:
            # tudo abaixo continua desenhando em coordenadas logicas
            # normais (self.width/self.height inalterados) -- so' o
            # temp_cr (e portanto o temp_surface, fisicamente maior)
            # recebe a escala, entao margens/fontes/grade crescem juntos.
            self.temp_cr.scale(export_scale, export_scale)
        
        #print(np.min(self.norm_data), np.max(self.norm_data) )
        
        self.draw_color_bar (self.temp_cr, res = self.size_y)#(cr, res = self.size_y)

        self.draw_image (self.temp_cr)#(cr)

        if self.show_contours:
            levels = self.get_contour_levels ( )
            labels = self.compute_contour_labels ( levels = levels )
            label_gaps = { id ( label['segment'] ) for label in labels }
            self.draw_contours (self.temp_cr, levels = levels, label_gaps = label_gaps)
            self.draw_contour_labels (self.temp_cr, labels)

        self.draw_image_box (self.temp_cr, line_width = 1, color = self.axis_color)#(cr, line_width = 1, color = [0,0,0])

        self.draw_scale(self.temp_cr)#(cr)
        
        #this draws the black dots at the image matrix surface (red dot if it is selected)
        self.draw_dots (self.temp_cr, self.points, [0,0,0], True)#(cr)
        
        self.draw_dots (self.temp_cr, self.extra_points, [0.75, 0, 0.5], False)

        self.draw_lines (self.temp_cr, line_width = 2, color = [0,0,0])#(cr, line_width = 2, color = [0,0,0])
    
        cr.set_source_surface(self.temp_surface, 0, 0)
        cr.paint()
    
        #rgb = self.get_pixel_rgb_from_surface(self.temp_surface, 0, 0)
        #print(f"RGB at (0, 0): {rgb}")
    
        #width, height = 300, 300
        #surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
        #surface.write_to_png("mesh_pattern_example.png")
        
    def on_motion (self, widget, event):
        """ Function doc """
        if self.norm_data is not None:
            pass
        else:
            return False
        x_on_plot, y_on_plot, x, y = self.get_i_and_j_from_click_event (event)
        i = y_on_plot
        j = x_on_plot
        
        j_size = len(self.norm_data)
        i_size = len(self.norm_data[0])

        #rgb = self.get_pixel_rgb_from_surface(self.temp_surface, x, y)
        #print(f"RGB at ({x}, {y}): {rgb}")

        if i < 0 or j < 0:
            pass
        else:
            if i < i_size and j < j_size:
                text = 'i {}    |    j {}    |    rc1 {:4.2f}    |    rc2 {:4.2f}    |    E = {:6.3f}'.format(i, j,  self.dataRC1[j][i], self.dataRC2[j][i],  self.data[j][i])

                if self.RC_label:
                    self.RC_label.set_text(text)
                else:
                    dprint(text)


    def get_pixel_rgb_from_surface(self, surface, x, y):
        """Lê o valor RGB de um pixel dentro de uma cairo.ImageSurface."""
        
        # Pega os dados da superfície de Cairo
        buf = surface.get_data()

        # Converte os dados em uma array numpy (ARGB formato, 4 bytes por pixel)
        array = np.frombuffer(buf, dtype=np.uint8)
        array = array.reshape((surface.get_height(), surface.get_width(), 4))

        # Pega o valor do pixel na posição especificada (x, y) (formato ARGB)
        pixel = array[y, x]

        # Extraímos os valores de Red, Green e Blue
        r, g, b = pixel[1], pixel[2], pixel[3]

        return (r, g, b)
