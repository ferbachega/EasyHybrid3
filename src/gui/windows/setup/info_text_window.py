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
'''
import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk

class TextBufferExample(Gtk.Window):
    def __init__(self):
        Gtk.Window.__init__(self, title="Text Buffer Example")

        # Creating a text buffer
        self.text_buffer = Gtk.TextBuffer()

        # Adding text to the buffer
        self.text_buffer.set_text("Welcome to GTK3 Text Buffer Example!")

        # Creating a text view
        text_view = Gtk.TextView(buffer=self.text_buffer)

        # Adding the text view to the window
        self.add(text_view)

window = TextBufferExample()
window.connect("delete-event", Gtk.main_quit)
window.show_all()
Gtk.main()
'''

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk

class HighlightKeywordsExample(Gtk.Window):
    def __init__(self):
        Gtk.Window.__init__(self, title="Highlight Keywords Example")

        # Creating a text buffer
        self.text_buffer = Gtk.TextBuffer()

        # Adding text to the buffer
        self.text_buffer.set_text("Welcome to GTK3 Text Buffer Example!")

        # Defining the keywords to highlight
        self.keywords = ["GTK3", "Text Buffer"]

        # Creating a tag for each keyword
        for keyword in self.keywords:
            tag = self.text_buffer.create_tag(None, foreground="red")
            start_iter = self.text_buffer.get_start_iter()
            while True:
                # Searching for the keyword
                result = start_iter.forward_search(keyword, 0)
                if result is None:
                    break
                start, end = result
                self.text_buffer.apply_tag(tag, start, end)
                start_iter = end

        # Creating a text view
        text_view = Gtk.TextView(buffer=self.text_buffer)

        # Adding the text view to the window
        self.add(text_view)

window = HighlightKeywordsExample()
window.connect("delete-event", Gtk.main_quit)
window.show_all()
Gtk.main()
