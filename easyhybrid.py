#!/usr/bin/env python3
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
EASYHYBRID_VERSION = '3.0.3'

import os, sys, time, re

'''
O módulo re é nativo do Python (faz parte da biblioteca padrão).  
Você não precisa instalar nada.Ele implementa expressões regulares 
(regular expressions) para busca e manipulação de texto.
'''


'''
os.environ["PDYNAMO3_HOME"]         = "/home/fernando/programs/pDynamo3"
os.environ["PDYNAMO3_PARAMETERS"]   = "/home/fernando/programs/pDynamo3/parameters"
os.environ["PDYNAMO3_SCRATCH"]      = "/home/fernando/programs/pDynamo3/scratch"
os.environ["PDYNAMO3_ORCACOMMAND"]  = "/home/fernando/programs/orca_6_1_0_linux_x86-64_shared_openmpi418/orca"
os.environ["PDYNAMO3_XTBCOMMAND"]   = "/home/fernando/programs/xtb-6.6.1/bin/xtb"
os.environ["PDYNAMO3_MOPACCOMMAND"] = "/home/fernando/programs/MOPAC2016A/bin/mopac"
os.environ["PDYNAMO3_DFTBCOMMAND"]  = "/home/fernando/programs/dftbplus-24.1.x86_64-linux/bin/dftb+"
sys.path.append("/home/fernando/programs/pDynamo3")
#sys.path.append("/home/fernando/programs/amber24_src/lib/python3.12/site-packages")
'''

#filepath = '/home/fernando/programs/pDynamo3/installation/shellScripts/environment_bash.com'

def parse_bash_env_file(filepath):
    """
    Parse a bash environment file and return a dictionary
    with environment variables and values.
    """
    env_vars = {}
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()

            # Ignore empty lines and comments
            if not line or line.startswith("#"):
                continue

            # Remove inline comments
            line = line.split("#")[0].strip()

            # Remove "; export VAR"
            line = re.sub(r";\s*export\s+\w+", "", line)

            if "=" in line:
                var, value = line.split("=", 1)
                var = var.strip()
                value = value.strip()

                env_vars[var] = value

    # Expand variables
    for var, value in env_vars.items():
        # Expand variables like $HOME, $PDYNAMO3_HOME etc
        value = os.path.expandvars(value)
        # Expand ~
        value = os.path.expanduser(value)
        env_vars[var] = value
        # Add to environment
        os.environ[var] = value

    # Add pDynamo to python path
    if "PDYNAMO3_HOME" in env_vars:
        sys.path.append(env_vars["PDYNAMO3_HOME"])

    return env_vars


try:
    from paths import PDYNAMO_HOME
    shell_scripts = os.path.join(PDYNAMO_HOME, 'installation/shellScripts/environment_bash.com')
    parse_bash_env_file(shell_scripts)
    print('importing evironment variables from:', shell_scripts)
except:
    print ('paths not found')





import logging
import gi 
gi.require_version("Gtk", "3.0")
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
        self.set_default_size(605 ,605)
        try:
            # Carrega imagem do splash
            pixbuf = GdkPixbuf.Pixbuf.new_from_file_at_scale(
                filename=(os.path.join(EASYHYBRID_HOME, "splash.png")),   
                width=605,
                height=605,
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
        
        #main_window.builder.get_object('toolbutton_export_img').hide()
        #main_window.builder.get_object('button_task_list')           .hide()
        
        main_window.builder.get_object('test_item')             .hide() # IR spectrum
        
        
        
        
        #main_window.builder.get_object('toolbutton_terminal')   .hide()
        #main_window.builder.get_object('menuitem_reimaging')   .hide()
        #main_window.builder.get_object('menuitem_RMSD_tool')   .hide()
        
        main_window.builder.get_object('menuitem_rama')   .hide()
        #main_window.builder.get_object('menuitem_advanced_rc_scans').hide()
        
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

