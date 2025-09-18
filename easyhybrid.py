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
EASYHYBRID_VERSION = '3.0.1'

import os, sys, time
import logging
import gi 
gi.require_version("Gtk", "3.0")
#from gi.repository import Gtk, Gdk
from gi.repository import Gtk, Gdk, GdkPixbuf, GLib

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
import time
import threading

# Splash Screen
class SplashScreen(Gtk.Window):
    def __init__(self):
        super().__init__(title="splash.png")
        self.set_decorated(False)  # Sem bordas
        self.set_position(Gtk.WindowPosition.CENTER)
        #self.set_default_size(710 ,710  )
        self.set_default_size(800 ,671  )
        try:
            # Carrega imagem do splash
            pixbuf = GdkPixbuf.Pixbuf.new_from_file_at_scale(
                filename=(os.path.join(EASYHYBRID_HOME, "splash.png")),   
                width=800,
                height=671,
                preserve_aspect_ratio=True
            )

            image = Gtk.Image.new_from_pixbuf(pixbuf)
            self.add(image)
        except:
            print('splash.png file not found!')

def load_modules(callback_final):
    def _load():
        print("Starting module loading...")
        time.sleep(1.5)   
        GLib.idle_add(callback_final)  #Call the finalize function in the main loop
    threading.Thread(target=_load).start()


def main():
    splash = SplashScreen()
    splash.show_all()
    
    def on_finalizado():
        logging.basicConfig(format="%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s",
                            datefmt="%Y-%m-%d:%H:%M:%S", level=logging.DEBUG)

        vconfig = VismolConfig(home = EASYHYBRID_HOME)

        vm_session = EasyHybridSession(vm_config = vconfig)
        vm_session.vm_widget.insert_glmenu()
        main_window = MainWindow(vm_session = vm_session,
                                 home       =  EASYHYBRID_HOME,
                                 version    = EASYHYBRID_VERSION)
        vm_session.main_session = main_window                  
        #print(vm_session.vm_config.gl_parameters)
        #main_window.window.connect('destroy', Gtk.main_quit)
        
        
        # do now show these itens:
        main_window.builder.get_object('toolbutton_monte_carlo').hide()
        #main_window.builder.get_object('button_test')           .hide()
        main_window.builder.get_object('test_item')             .hide()
        #main_window.builder.get_object('toolbutton_terminal')   .hide()
        #main_window.builder.get_object('menuitem_reimaging')   .hide()
        #main_window.builder.get_object('menuitem_RMSD_tool')   .hide()
        main_window.builder.get_object('menuitem_rama')   .hide()
        splash.destroy()
        try:
            filein = sys.argv[-1]
            vm_session.load_molecule(filein)
        except:
            pass
        #Gtk.main()
        return 0

    load_modules(on_finalizado)

    Gtk.main()


if __name__ == "__main__":
    main()

