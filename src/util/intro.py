#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  intro.py
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


#import gi
#gi.require_version('Gtk', '3.0')
#from gi.repository import Gtk
#
#class ImageWindow(Gtk.Window):
#    def __init__(self):
#        Gtk.Window.__init__(self, title="Image Example")
#        self.set_decorated(False)
#        image = Gtk.Image()
#        image.set_from_file("/home/fernando/Desktop/Screenshot from 2023-02-06 02-10-52.png")
#        
#        self.add(image)
#
#win = ImageWindow()
#win.connect("delete-event", Gtk.main_quit)
#win.show_all()
#Gtk.main()


import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, GObject
import threading
import time

# Função para simular o carregamento dos módulos (substitua isso com suas próprias tarefas de carregamento)
def load_modules():
    time.sleep(3)  # Simula o carregamento de módulos
    return True

class LoaderWindow(Gtk.Window):
    def __init__(self):
        Gtk.Window.__init__(self, title="Carregando Módulos")
        self.set_default_size(400, 300)

        # Logo
        self.logo = Gtk.Image.new_from_file("/home/fernando/programs/EasyHybrid3/gui/icons/easyhybrid_solo_100x100.png")  # Substitua "logo.png" pelo caminho do seu logo
        self.add(self.logo)

        self.load_thread = threading.Thread(target=self.load_modules_thread)
        self.load_thread.start()

    def load_modules_thread(self):
        # Realize aqui suas tarefas de carregamento
        loaded = load_modules()

        # Após o carregamento, chame a função de callback na thread principal
        GObject.idle_add(self.load_complete, loaded)

    def load_complete(self, loaded):
        if loaded:
            self.destroy()
        else:
            # Trate o caso de erro, se necessário
            pass

def main():
    Gdk.threads_init()
    window = LoaderWindow()
    window.connect("delete-event", Gtk.main_quit)
    window.show_all()
    Gtk.main()

if __name__ == "__main__":
    main()
