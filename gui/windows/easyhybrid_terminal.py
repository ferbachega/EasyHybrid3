#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  easyhybrid_pDynamo_selection.py
#  
#  Copyright 2022 Fernando <fernando@winter>
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
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
from gi.repository import Gdk

import os

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')

class CommandLineInterpreter:
    """ Class doc """
    
    def __init__(self, main = None):
        """ Class initialiser """
        self.main       = main
        self.visible               =  False        
        self.p_session             = main.p_session
        self.vm_session            = main.vm_session
        

    def cmd (self, command_line):
        """ Function doc """
        


        
        
class TerminalWindow():
    """ Class doc """
    
    def OpenWindow (self):
        """ Function doc /home/fernando/programs/VisMol/easyhybrid/gui/selection_list.glade"""
        if self.visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.main.home,'gui/windows/easyhybrid_terminal.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('window')
            #self.window.set_border_width(10)
            self.window.set_default_size(450, 250)  
            self.window.set_title('EasyHybrid Terminal')  

            self.window.connect('destroy-event', self.CloseWindow)
            self.window.connect('delete-event', self.CloseWindow)

            self.window.set_keep_above(True)
            self.entry_terminal = self.builder.get_object('entry_terminal')
            #self.entry_terminal.set_placeholder_text('>>>')
            #self.entry_terminal.do_insert_at_cursor ('>')
            self.textview = self.builder.get_object('entry_text_buffer')
            self.textview.set_buffer(self.textbuffer)
            self.window.show_all()                                               
            #self.system_names_combo.set_active(0)
            self.visible    =  True
            '''--------------------------------------------------------------------------------------------'''


    def CloseWindow (self, button, data  = None):
        """ Function doc """
        #self.BackUpWindowData()
        self.window.destroy()
        self.visible    =  False
        self.main.cmd_terminal_button.set_active(False)
        #print('self.visible',self.visible)
    
    def __init__(self, main = None):
        """ Class initialiser """
        self.main                  = main
        self.visible               =  False        
        self.p_session             = main.p_session
        self.vm_session            = main.vm_session
        
        self.cmd_history = []
        self.cmd_history_counter = 0
        self.textbuffer          = Gtk.TextBuffer()
        self.command_list= {
                          'list'   : True,
                          'show'   : True,
                          'hide'   : True,
                          'bond'   : True,
                          'unbond' : True,
                          'rename' : True,
                          'delete' : True,
                          'save'   : True,
                          }
    
    
    def run_cmd (self, cmd):
        """ Function doc """
        text  = '\n'+cmd
        end_iter = self.textbuffer.get_end_iter ()
        self.textbuffer.insert(end_iter, text)
        print(self.cmd_history)
    
    def on_entry_terminal (self, widget):
        """ Function doc """
        cmd = self.entry_terminal.get_text()
        self.cmd_history.append(cmd)
        self.entry_terminal.set_text('')
        self.run_cmd (cmd)
        print(self.cmd_history)
    
    def on_entry_terminal_backspace (self, widget):
        pass
        #print('on_entry_terminal_backspace', widget)
    
    def on_entry_terminal_move_cursor (self, widget, data, data2, data3 ):
        pass
        #print('on_entry_terminal_move_cursor', widget, data, data2, data3)

    def on_entry_terminal_change (self, widget):
        """ Function doc """
        pass
        #print('on_entry_terminal_change', widget)
    
    
    def update_window (self, system_names = True, coordinates = False,  selections = True ):
        """ Function doc """
        pass
        
    def on_key_press_event (self, widget, event):
        """ Function doc """
        size = -len(self.cmd_history)
        
        if size == 0:
            return None
        else:
            pass
        
        k_name = Gdk.keyval_name(event.keyval)
        print ('on_key_press_event', widget, event, k_name)

        if k_name in ['Down','Up' ]:

            if k_name == 'Up':
                self.cmd_history_counter += -1
                print('Up key!')
            
            if k_name == 'Down':
                self.cmd_history_counter += 1
                print('Down key!')
            
            print(self.cmd_history_counter, size)
            
            if self.cmd_history_counter >= 0:
                self.entry_terminal.set_text('')
                self.cmd_history_counter = 0
            
            elif self.cmd_history_counter < 0 and self.cmd_history_counter > size:
                
                self.entry_terminal.set_text(self.cmd_history[self.cmd_history_counter])

            else:
                self.cmd_history_counter = size
                self.entry_terminal.set_text(self.cmd_history[self.cmd_history_counter])
                #print(self.cmd_history_counter, size)
            
            


        elif k_name =='Tab' :
            text = self.entry_terminal.get_text()
            for key in self.command_list:
                if  text in key:
                    print(key)
