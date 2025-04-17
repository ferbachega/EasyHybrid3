#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  simple.py
#
#  Copyright 2022 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
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
EASYHYBRID_VERSION = '3.0'

import os, sys, time
import logging
import gi 
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk

#               Installation is not necessary anymore.
#.This retrieves the absolute path of the script file that is currently being executed.
easy_main_file = os.path.abspath(__file__)

#.This extracts the directory (folder) from the absolute path obtained in the previous step.
EASYHYBRID_HOME = os.path.dirname(easy_main_file)

#.Adding GRAPHIC ENGINE LIB
sys.path.append(os.path.join(EASYHYBRID_HOME,"src/graphics_engine/src"))
sys.path.append(os.path.join(EASYHYBRID_HOME,"src/"))

from gui.main import MainWindow
from gui.eSession import EasyHybridSession
from gui.config   import VismolConfig


def main():
    logging.basicConfig(format="%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s",
                        datefmt="%Y-%m-%d:%H:%M:%S", level=logging.DEBUG)

    vconfig = VismolConfig(home = EASYHYBRID_HOME)

    vm_session = EasyHybridSession(vm_config = vconfig)
    vm_session.vm_widget.insert_glmenu()
    main_window = MainWindow(vm_session = vm_session,
                             home       =  EASYHYBRID_HOME,
                             version    = EASYHYBRID_VERSION)
    vm_session.main_session = main_window                  
    #main_window.window.connect('destroy', Gtk.main_quit)
    
    
    
    main_window.builder.get_object('toolbutton_monte_carlo').hide()
    main_window.builder.get_object('button_test')           .hide()
    main_window.builder.get_object('test_item')             .hide()
    main_window.builder.get_object('toolbutton_terminal')   .hide()
    
    main_window.builder.get_object('menuitem_RMSD_tool')   .hide()
    main_window.builder.get_object('menuitem_rama')   .hide()
    
    try:
        filein = sys.argv[-1]
        vm_session.load_molecule(filein)
    except:
        pass
    Gtk.main()
    return 0

if __name__ == "__main__":
    main()

