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
from gi.repository import Gtk, Pango
from gi.repository import Gdk
from gi.repository import GdkPixbuf
import cairo


#---------------------------------------
from pBabel                    import*                                     
from pCore                     import*  
#---------------------------------------
from pMolecule                 import*                              
from pMolecule.MMModel         import*
from pMolecule.NBModel         import*                                     
from pMolecule.QCModel         import*
#---------------------------------------
from pScientific               import*                                     
from pScientific.Arrays        import*                                     
from pScientific.Geometry3     import*                                     
from pScientific.RandomNumbers import*                                     
from pScientific.Statistics    import*
from pScientific.Symmetry      import*
#---------------------------------------                              
from pSimulation               import*
#---------------------------------------

#import Pickle
from gui.widgets.custom_widgets import FolderChooserButton
from gui.widgets.custom_widgets import VismolTrajectoryFrame
from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.widgets.custom_widgets import get_colorful_square_pixel_buffer
from gui.widgets.custom_widgets import ReactionCoordinateBox

#from gui.widgets.custom_widgets import get_distance
from util.geometric_analysis            import get_distance 
from util.geometric_analysis            import get_dihedral 
from util.geometric_analysis            import get_angle 

from pdynamo.p_methods import LogFile


import util.orca_qc_keywords as orca_keys

from util.file_parser import get_file_type  
from util.file_parser import read_MOL2  
import pprint
import numpy as np
import gc
import os

import traceback


VISMOL_HOME      = os.environ.get('VISMOL_HOME')
HOME             = os.environ.get('HOME')
PDYNAMO3_SCRATCH = os.environ.get('PDYNAMO3_SCRATCH')


def call_message_dialog (text1 = '', text2 = '', transient_for = None):
    """ Function doc """
    dialog = Gtk.MessageDialog(
                transient_for = transient_for,
                flags = 0,
                message_type=Gtk.MessageType.INFO,
                buttons = Gtk.ButtonsType.OK,
                text = text1, #"MMModelError",
                )
    dialog.format_secondary_text(text2)#"""Total active MM charge is neither integral nor zero.""")
    dialog.run()
    #print("INFO dialog closed")
    dialog.destroy()


class SimpleDialog:
    """ A helper class to create simple dialog windows using GTK. """
    
    def __init__(self, main = None, is_on_top = True):
        """ 
        Class initializer. 
        Stores a reference to the main application instance (main), 
        which is expected to have a window attribute used as the parent for dialogs. 
        """
        self.main = main
        self.is_on_top = is_on_top
    def create_finished_dialog(self, parent=None, 
            msg1 = '',
            msg2 = '', 
            modal = False):
                
        # Criar a janela de diálogo
        dialog = Gtk.Dialog(
            title="Process finished",
            transient_for=parent,
            flags=0
        )
        dialog.add_button(Gtk.STOCK_OK, Gtk.ResponseType.OK)
        dialog.set_default_size(350, 150)
        dialog.set_border_width(5)

        # Área de conteúdo
        content_area = dialog.get_content_area()

        # Caixa vertical principal
        vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=10)
        vbox.set_halign(Gtk.Align.CENTER)
        vbox.set_valign(Gtk.Align.CENTER)
        content_area.add(vbox)

        # Texto acima (centralizado)
        top_label = Gtk.Label(label=msg1)
        top_label.set_justify(Gtk.Justification.CENTER)
        vbox.pack_start(top_label, False, False, 0)

        # Caixa horizontal para quadrado + mensagem
        hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=10)
        vbox.pack_start(hbox, False, False, 0)

        # Quadrado verde 15x15
        pixbuf = GdkPixbuf.Pixbuf.new(GdkPixbuf.Colorspace.RGB, True, 8, 15, 15)
        pixbuf.fill(0x00FF00FF)  # RGBA → verde
        image = Gtk.Image.new_from_pixbuf(pixbuf)
        hbox.pack_start(image, False, False, 0)

        # Texto ao lado do quadrado
        msg_label = Gtk.Label(label=msg2)
        msg_label.set_xalign(0)
        hbox.pack_start(msg_label, False, False, 0)

        vbox.show_all()
        return dialog

    def finished(self, msg=None, modal=True):
        """ 
        Show an information dialog with an OK button and a green square.
        """

        if modal:
            flags = Gtk.DialogFlags.MODAL
        else:
            print('not modal')
            flags = 0

        dialog = Gtk.MessageDialog(
            parent=self.main.window,
            flags=flags,
            type=Gtk.MessageType.INFO,
            buttons=Gtk.ButtonsType.OK,
            message_format=None  # vamos criar nosso próprio layout
        )

        # Criar a área de conteúdo
        content_area = dialog.get_content_area()

        # Caixa horizontal para alinhar ícone + texto
        hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
        content_area.add(hbox)

        # Criar um quadrado verde (pixbuf 15x15)
        surface = Gdk.cairo_surface_create_from_pixbuf(
            GdkPixbuf.Pixbuf.new(GdkPixbuf.Colorspace.RGB, True, 8, 15, 15),
            1, None
        )
        # Pintar o surface de verde
        cr = cairo.Context(surface)
        cr.set_source_rgb(0, 1, 0)  # verde
        cr.rectangle(0, 0, 15, 15)
        cr.fill()

        pixbuf = Gdk.pixbuf_get_from_surface(surface, 0, 0, 15, 15)
        image = Gtk.Image.new_from_pixbuf(pixbuf)
        hbox.pack_start(image, False, False, 0)

        # Adicionar texto
        label = Gtk.Label(label=msg)
        label.set_xalign(0)  # alinhado à esquerda
        hbox.pack_start(label, True, True, 0)

        hbox.show_all()
        dialog.set_border_width(15)
        dialog.show_all()

    def info(self, msg = None, modal = True, title = 'Title' ):
        """ 
        Show an information dialog with an OK button. 
        The dialog is modal, meaning the user must dismiss it before interacting with the main window. 
        """
        
        if modal:
            flags = Gtk.DialogFlags.MODAL
        else:
            print('not modal')
            flags = 0
            
        dialog = Gtk.MessageDialog(
                title=title,
                parent=self.main.window,     # Parent window for proper modality
                flags=flags,                 # Makes the dialog block input to the parent
                type=Gtk.MessageType.INFO,   # Sets the dialog type to "information"
                buttons=Gtk.ButtonsType.OK,  # Provides a single OK button
                message_format=msg           # The text message displayed to the user
            )
 
 
        if modal:
            dialog.run()   # Displays the dialog and waits for the user response
            dialog.destroy()  # Closes and frees resources after user dismisses it
        else:
            dialog.connect("response", lambda d, r: d.destroy())
            dialog.show_all()
        
        return None
    
    def error(self, msg):
        """ 
        Show an error dialog with an OK button. 
        Used to notify the user of an error condition. 
        """
        dialog = Gtk.MessageDialog(
                    parent=self.main.window,
                    flags=Gtk.DialogFlags.MODAL,
                    type=Gtk.MessageType.ERROR,     # Sets the dialog type to "error"
                    buttons=Gtk.ButtonsType.OK,
                    message_format=msg
                )
        dialog.run()
        dialog.destroy()
        return None

    def error_details (self, parent, msg, details = None, title='Error!'):
        """ Function doc """
        dialog = Gtk.MessageDialog(title=title,
            transient_for=parent,
            modal=True,
            destroy_with_parent=True,
        )
        # Main message shown to the user (short description of the error)
        dialog.format_secondary_text(msg)
        # Revealer is used to show/hide the error details (stacktrace/log)
        dialog.revealer = Gtk.Revealer()
        dialog.revealer.set_reveal_child(False)

        # Scrolled window to contain the log output (scrollable for long text)
        scrolled = Gtk.ScrolledWindow()
        scrolled.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
        scrolled.set_min_content_height(250)
        scrolled.set_hexpand(True)
        scrolled.set_vexpand(True)

        # TextView to display the traceback/log (non-editable)
        self.textview = Gtk.TextView()
        self.textview.set_editable(False)

        # Define style (monospace font for technical logs)
        css_provider = Gtk.CssProvider()
        css_provider.load_from_data(b"""
            textview {
                font-family: monospace;
                font-size: 10pt;
            }
        """)

        context = self.textview.get_style_context()
        context.add_provider(css_provider, Gtk.STYLE_PROVIDER_PRIORITY_USER)

        # Insert traceback or log into the TextView
        buffer = self.textview.get_buffer()
        buffer.set_text(details)

        scrolled.add(self.textview)
        self.revealer.add(scrolled)

        # Extra buttons
        self.btn_details = Gtk.Button(label="Show details")
        self.btn_details.connect("clicked", self.on_toggle_details)

        # Optional "Copy log" button (currently disabled)
        # self.btn_copy = Gtk.Button(label="Copy log")
        # self.btn_copy.connect("clicked", self.on_copy_log)

        # Add buttons and widgets to the dialog’s message area
        box = self.get_message_area()
        # box.pack_end(self.btn_copy, False, False, 0)
        box.pack_end(self.btn_details, False, False, 0)
        box.pack_end(self.revealer, True, True, 0)

        # Add default close button
        dialog.add_button("OK", Gtk.ResponseType.CLOSE)
        dialog.show_all()



            
    def question(self, msg):
        """ 
        Show a question dialog with YES and NO buttons. 
        Returns True if the user clicks YES, False if the user clicks NO. 
        """
        dialog = Gtk.MessageDialog(
            parent=self.main.window,
            flags=Gtk.DialogFlags.MODAL,
            type=Gtk.MessageType.QUESTION,        # Sets the dialog type to "question"
            buttons=Gtk.ButtonsType.YES_NO,       # Provides Yes and No buttons
            message_format=msg
        )
        if self.is_on_top:
            dialog.set_keep_above(True)
            
        response = dialog.run()   # Waits for the user's response
        dialog.destroy()
        if response == Gtk.ResponseType.YES:
            return True
        elif response == Gtk.ResponseType.NO:
            return False
