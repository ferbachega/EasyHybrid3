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
#  main/  --  layout apos refatoracao
#
#  src/gui/
#  +-- main/                      (pacote; __init__.py e a fachada, re-exporta
#      |                           MainWindow p/ manter o import antigo)
#      +-- main_window.py          MainWindow            (janela principal)
#      +-- main_treeview.py        EasyHybridMainTreeView (arvore de sistemas)
#      +-- treeview_menu.py        TreeViewMenu          (menu de contexto)
#      +-- treeview_menu_new.py    TreeViewMenu_new      (CODIGO MORTO - sem uso)
#      +-- bottom_notebook.py      BottomNoteBook        (painel inferior)
#      +-- preferences_window.py   PreferencesWindow     (preferencias)
# -----------------------------------------------------------------------------

from gui.main.main_window import MainWindow
from gui.main.main_treeview import EasyHybridMainTreeView
from gui.main.treeview_menu import TreeViewMenu
from gui.main.treeview_menu_new import TreeViewMenu_new
from gui.main.bottom_notebook import BottomNoteBook
from gui.main.preferences_window import PreferencesWindow
