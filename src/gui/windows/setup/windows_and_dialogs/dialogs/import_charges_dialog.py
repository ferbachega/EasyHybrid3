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
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk


class ChargeTypeSelectionDialog(Gtk.Dialog):
    def __init__(self, parent, items):
        super().__init__(
            #title="Escolha um item",
            transient_for=parent,
            flags=0
        )

        self.set_default_size(350, 300)

        self.add_button("Cancel", Gtk.ResponseType.CANCEL)
        self.add_button("OK", Gtk.ResponseType.OK)

        self.selected_index = None

        box = self.get_content_area()

        # bool (radio selecionado), string (nome)
        self.store = Gtk.ListStore(bool, str)

        for item in items:
            self.store.append([False, item])

        self.treeview = Gtk.TreeView(model=self.store)

        # Coluna do radio button
        renderer_toggle = Gtk.CellRendererToggle()
        renderer_toggle.set_radio(True)
        renderer_toggle.connect("toggled", self.on_radio_toggled)

        col_toggle = Gtk.TreeViewColumn("Select", renderer_toggle, active=0)
        self.treeview.append_column(col_toggle)

        # Coluna do texto
        renderer_text = Gtk.CellRendererText()
        col_text = Gtk.TreeViewColumn("Item", renderer_text, text=1)
        self.treeview.append_column(col_text)

        scrolled = Gtk.ScrolledWindow()
        scrolled.add(self.treeview)

        box.add(scrolled)

        self.show_all()

    def on_radio_toggled(self, widget, path):
        for row in self.store:
            row[0] = False

        self.store[path][0] = True
        self.selected_index = int(path)

    def get_selected_item(self):
        if self.selected_index is None:
            return None
        return self.store[self.selected_index][1]

