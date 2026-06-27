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

def interpolate_color(color1, color2, fraction):
    return [ ( (x1 + fraction * (x2 - x1))) for (x1, x2) in zip(color1, color2)]


def get_color(value, color_map):#, factor = 2):   
    
    if value == float('inf'):
        return [1,1,1]
    #value = (value+1)/factor   
    
    cutoffs = list(color_map.keys())
    color_list = list(color_map.values())
    counter = 0 
    
    #print(value)
    
    
    
    for cutoff, color in color_map.items():
        
        if value <= cutoff :
            
            fraction = (value -  cutoffs[counter-1]) / ( cutoff - cutoffs[counter-1])
            
            
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
            #print(value, cutoff, fraction,color )

            return color
        else:
            pass
            #print(value, cutoff)#, fraction,color )
        counter +=1


def bilinear_interpolation(c1, c2, c3, c4, tx, ty):
    return (
        int((1 - tx) * (1 - ty) * c1[0] + tx * (1 - ty) * c2[0] + (1 - tx) * ty * c4[0] + tx * ty * c3[0]),
        int((1 - tx) * (1 - ty) * c1[1] + tx * (1 - ty) * c2[1] + (1 - tx) * ty * c4[1] + tx * ty * c3[1]),
        int((1 - tx) * (1 - ty) * c1[2] + tx * (1 - ty) * c2[2] + (1 - tx) * ty * c4[2] + tx * ty * c3[2]),
        int((1 - tx) * (1 - ty) * c1[3] + tx * (1 - ty) * c2[3] + (1 - tx) * ty * c4[3] + tx * ty * c3[3])
    )


def expand_array_size (norm_data):
    """ 
    this function receives an i x j matrix and transforms it into an 
    i+2 x j+2 matrix. the original matrix is ​​in the middle of the new 
    matrix, so the added elements are on the edges. they will be 
    important for the interpolation of points when the smooth function 
    is used in ImagePlot.
    
    
       extra column /extra row
         |
       |---|---|---|---|---|
     --| x | x | x | x | x |
       |---|---|---|---|---|
       | x |###|###|###| x |
       |---|---|---|---|---|
       | x |###|###|###| x |
       |---|---|---|---|---|
       | x | x | x | x | x |--
       |---|---|---|---|---|
                         |
                       extra column /extra row
    
    ### = elements in original matrix
    x   = new elements
    
    -It's not a very efficient algorithm, but it works well
    """
    size_x = len(norm_data)
    size_y = len(norm_data[0])
    #print('expanded matrix')
    bigmatrix = []#[None]*size_x+2
    for k in range (0, size_x+2):
        bigmatrix.append([])
        for l in range(0,  size_y+2):
            bigmatrix[k].append(0)
    #pprint (bigmatrix )    
    for k in range (0, size_x+2):    
        for l in range(0,  size_y+2):
           
            #    extra column /extra row
            #      |
            #    |---|---|---|---|---|
            #  --|XXX|   |   |   |   |
            #    |---|---|---|---|---|
            #    |   |###|###|###|   |
            #    |---|---|---|---|---|
            #    |   |###|###|###|   |
            #    |---|---|---|---|---|
            #    |   |   |   |   |   |--
            #    |---|---|---|---|---|
            #                      |
            #                    extra column /extra row
            
            if   k == 0 and l == 0:
                bigmatrix[k][l] = norm_data[0][0]
            
            elif k == 0 and l == 1:
                bigmatrix[k][l] = norm_data[0][0]
            
            elif k == 1 and l == 0:
                bigmatrix[k][l] = norm_data[0][0]

            
            #    extra column /extra row
            #      |
            #    |---|---|---|---|---|
            #  --|   |   |   |   |   |
            #    |---|---|---|---|---|
            #    |   |###|###|###|   |
            #    |---|---|---|---|---|
            #    |   |###|###|###|   |
            #    |---|---|---|---|---|
            #    |   |   |   |   |XXX|--
            #    |---|---|---|---|---|
            #                      |
            #                    extra column /extra row
            
            elif k == size_x+1   and l == size_y+1 :
                bigmatrix[k][l] = norm_data[size_x-1][size_y-1]
            
            elif k == size_x   and l == size_y+1 :
                bigmatrix[k][l] = norm_data[size_x-1][size_y-1]
            
            elif k == size_x+1    and l == size_y :
                bigmatrix[k][l] = norm_data[size_x-1][size_y-1]

                
            #    extra column /extra row
            #      |
            #    |---|---|---|---|---|
            #  --|   |XXX|XXX|XXX|   |
            #    |---|---|---|---|---|
            #    |   |###|###|###|   |
            #    |---|---|---|---|---|
            #    |   |###|###|###|   |
            #    |---|---|---|---|---|
            #    |   |   |   |   |   |--
            #    |---|---|---|---|---|
            #                      |
            #                    extra column /extra row
            
            elif k == 0 and   size_y+1 > l > 0  :
                bigmatrix[k][l] = norm_data[0][l-1]
            
            
            #    extra column /extra row
            #      |
            #    |---|---|---|---|---|
            #  --|   |   |   |   |   |
            #    |---|---|---|---|---|
            #    |XXX|###|###|###|   |
            #    |---|---|---|---|---|
            #    |XXX|###|###|###|   |
            #    |---|---|---|---|---|
            #    |   |   |   |   |   |-- 
            #    |---|---|---|---|---|
            #                      |
            #                    extra column /extra row
            
            elif size_x+1 > k > 0 and l == 0:
                bigmatrix[k][l] = norm_data[k-1][0]
            
            
            #---------------------------------------------------
            #    extra column /extra row
            #      |
            #    |---|---|---|---|---|
            #  --|   |   |   |   |   |
            #    |---|---|---|---|---|
            #    |   |###|###|###|xxx|
            #    |---|---|---|---|---|
            #    |   |###|###|###|xxx|
            #    |---|---|---|---|---|
            #    |   |   |   |   |   |-- 
            #    |---|---|---|---|---|
            #                      |
            #                    extra column /extra row
            
            elif size_x+1 > k > 0 and size_y+1 == l:
                bigmatrix[k][l] = norm_data[k-1][size_y-1]
            #---------------------------------------------------
               
            #---------------------------------------------------
            #    extra column /extra row
            #      |
            #    |---|---|---|---|---|
            #  --|   |   |   |   |   |
            #    |---|---|---|---|---|
            #    |   |###|###|###|   |
            #    |---|---|---|---|---|
            #    |   |###|###|###|   |
            #    |---|---|---|---|---|
            #    |   |xxx|xxx|xxx|   |-- 
            #    |---|---|---|---|---|
            #                      |
            #                    extra column /extra row
            
            elif k == size_x+1 and size_y  > l> 0:
                bigmatrix[k][l] = norm_data[k-2][l-1]
            #---------------------------------------------------


            elif   k == 0        and l == size_y+1:
                bigmatrix[k][l] = norm_data[0][l-2]
            
            elif k == size_x+1 and l == 0:
                bigmatrix[k][l] = norm_data[k-2][0]
            
            else:
                bigmatrix[k][l] = norm_data[k-1][l-1]
    return bigmatrix


class ColorSquareBilinearInterpolation:
    def __init__(self, surface, x, y, size, c1, c2, c3, c4):
        self.ctx = cairo.Context(surface)
        self.x, self.y, self.size = x, y, size
        self.c1, self.c2, self.c3, self.c4 = c1, c2, c3, c4

    def draw(self):
        for i in range(self.size):
            for j in range(self.size):
                tx = j / (self.size - 1)
                ty = i / (self.size - 1)
                color = bilinear_interpolation(self.c1, self.c2, self.c3, self.c4, tx, ty)

                self.ctx.set_source_rgba(*[comp / 255.0 for comp in color])
                self.ctx.rectangle(self.x + j, self.y + i, 1, 1)
                self.ctx.fill()
                
                cr.rectangle(self.bx+i*self.factor_x -1,
                             self.by+j*self.factor_y -1, 
                             round(self.factor_x)+1, 
                             round(self.factor_y)+1)
                cr.fill()
