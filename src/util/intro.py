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
