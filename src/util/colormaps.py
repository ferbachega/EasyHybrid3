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
