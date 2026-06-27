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
# -----------------------------------------------------------------------------
#  easyplot/  --  layout apos refatoracao
#
#  src/util/
#  +-- easyplot/                  (pacote; __init__.py e a fachada, re-exporta
#      |                           tudo p/ manter os imports antigos funcionando)
#      +-- color_utils.py          interpolate_color, get_color, bilinear_interpolation,
#      |                           expand_array_size, ColorSquareBilinearInterpolation
#      +-- canvas.py               Canvas (base, Gtk.DrawingArea)
#      +-- image_plot.py           ImagePlot(Canvas)
#      +-- xy_plots.py             XYPlot, XYScatterPlot
# -----------------------------------------------------------------------------

from util.easyplot.color_utils import (interpolate_color, get_color, bilinear_interpolation, expand_array_size, ColorSquareBilinearInterpolation)
from util.easyplot.canvas import Canvas
from util.easyplot.image_plot import ImagePlot
from util.easyplot.xy_plots import XYPlot, XYScatterPlot
