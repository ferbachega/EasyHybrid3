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
import sys
import io
import traceback
import time 
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
from gi.repository import Gdk, GLib

import os
import inspect

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')



class Command:
    def __init__(self, console, vm_session = None):
        self.console = console
        self.session = vm_session
        
    def help (self, obj_name="objeto"):
        '''help_cmd - commands'''
        self.console.write_output(f"[3D] help cmd", "green")
        metodos = [self for self in dir(Command) if callable(getattr(Command, self)) and not self.startswith("__")]
        
        for method in metodos:
            method_obj = getattr(self, method)
            print(method, method_obj.__doc__)
    
    
    def show(self, obj_name="objeto", rep='lines' ):
        '''show method'''
        self.console.write_output(f"[3D] Exibindo {obj_name}", "green")
        self.vm_session
    def hide(self, obj_name="objeto"):
        '''hide method'''
        self.console.write_output(f"[3D] Ocultando {obj_name}", "red")

    def rotate(self, obj_name="objeto", angle=90):
        '''rotate method'''
        self.console.write_output(f"[3D] Rotacionando {obj_name} em {angle} graus", "blue")

        
        
class TerminalWindow():
    """ Class doc """
    
    def OpenWindow (self):
        """ """
        if self.visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.main.home,'src/gui/windows/setup/easyhybrid_terminal.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('window')
            #self.window.set_border_width(10)
            self.window.set_default_size(650, 350)  
            self.window.set_title('EasyHybrid Terminal')  

            self.window.connect('destroy-event', self.CloseWindow)
            #self.window.connect('delete-event', self.CloseWindow)

            self.window.set_keep_above(True)
            
            # Tabela de tags para cores
            self.tag_table = self.textbuffer.get_tag_table()
            self.textbuffer.create_tag("green", foreground="green")
            self.textbuffer.create_tag("red", foreground="red")
            self.textbuffer.create_tag("blue", foreground="blue")
            #self.textbuffer.create_tag("normal", foreground="black")
            self.locals = {}
            # Contexto de execução
            self.cmd = Command(self, self.vm_session)
            self.locals['cmd'] = self.cmd

            
            self.entry_terminal = self.builder.get_object('entry_terminal')
            

            
            #self.entry_terminal.set_placeholder_text('>>>')
            #self.entry_terminal.do_insert_at_cursor ('>')
            self.textview = self.builder.get_object('entry_text_buffer')
            self.textview.set_buffer(self.textbuffer)
            
            self.textview.set_wrap_mode(Gtk.WrapMode.WORD)
            #--------------------------------------------------------------
            self.textview.get_style_context().add_provider(self.css_provider, Gtk.STYLE_PROVIDER_PRIORITY_USER)
            # Entrada de comandos
            self.entry_terminal.get_style_context().add_provider(self.css_provider, Gtk.STYLE_PROVIDER_PRIORITY_USER)
            #--------------------------------------------------------------
            
            
            self.window.show_all()                                               
            #self.system_names_combo.set_active(0)
            self.visible    =  True
            '''--------------------------------------------------------------------------------------------'''

        else:
            self.window.present()
            

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
        self.textbuffer = self.main.terminal_text_buffer
        #self.command_list= (dir(self.vm_session.cmd))
        self.command_list = []
        
        
        
        #adiciona o texto inicial ao terminal do EasyHybrid
        text = '''
        --------------------------------------------------------
               Welcome to the EasyHybrid Terminal. 
        --------------------------------------------------------
        
                   Created by J.F.R Bachega 
                   
        \n\n
        '''
        end_iter = self.textbuffer.get_end_iter()
        self.textbuffer.insert(end_iter, text)
        #self.methods = inspect.getmembers(self.vm_session.cmd, predicate=inspect.ismethod)
        
        # time controler
        self.last_tab_time = 0
        self.tab_timer_id = None
        
        #for method in self.methods:
        #    self.command_list.append(method[0])
        #print(self.command_list)
        
        
        #--------------------------------------------------------------
        self.css_provider = Gtk.CssProvider()
        self.css_provider.load_from_data(b"""
        textview {
            font-family: Monospace;
            font-size: 12pt;
        }
        
        entry {
            font-family: Monospace;
            font-size: 12pt;
        }
        """)
        #--------------------------------------------------------------

        
        
    def run_cmd (self, cmd):
        """ Function doc """
        text  = '\n>'+cmd
        end_iter = self.textbuffer.get_end_iter ()
        self.textbuffer.insert(end_iter, text)
        
        log = self.vm_session.cmd.run(cmd)
        if log is not None:
            self.textbuffer.insert(end_iter, log)
        #print(self.cmd_history)
    
    
    def write_output(self, text: str, color: str = "normal"):
        end_iter = self.textbuffer.get_end_iter()
        tag = self.tag_table.lookup(color)
        if not tag:
            tag = self.textbuffer.create_tag(color, foreground=color)
        self.textbuffer.insert_with_tags(end_iter, text + "\n", tag)
        self.textview.scroll_to_iter(end_iter, 0.0, False, 0, 0)
    
    
    def on_entry_terminal (self, widget):
        """ Function doc """
        #cmd = self.entry_terminal.get_text()
        #self.cmd_history.append(cmd)
        #self.entry_terminal.set_text('')
        #self.run_cmd (cmd)
        #print(self.cmd_history)
        
        command = self.entry_terminal.get_text()
        self.write_output(f">>> {command}", "normal")
        self.entry_terminal.set_text("")
        
        self.cmd_history.append(command)
        
        # Capturar stdout/erros
        old_stdout, old_stderr = sys.stdout, sys.stderr
        sys.stdout, sys.stderr = io.StringIO(), io.StringIO()

        try:
            # Primeiro tenta eval
            try:
                result = eval(command, {}, self.locals)
                if result is not None:
                    print(result)
            except SyntaxError:
                exec(command, {}, self.locals)
        except Exception:
            traceback.print_exc()

        # Mostrar saída
        output = sys.stdout.getvalue() + sys.stderr.getvalue()
        if output:
            self.write_output(output.strip(), "normal")

        sys.stdout, sys.stderr = old_stdout, old_stderr


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
        
        #if size == 0:
        #    return None
        #else:
        #    pass
        
        k_name = Gdk.keyval_name(event.keyval)
        #print ('on_key_press_event', widget, event, k_name)

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

        
        
        
        #print('\n\n')
        
        #elif k_name =='Tab' :
        #    text = self.entry_terminal.get_text()
        #    
        #    
        #    for key in self.command_list:
        #        if  text in key:
        #            print(key)

        if event.keyval == Gdk.KEY_Tab:
            now = time.time()
            if now - self.last_tab_time < 0.2:  # se for dentro de 300 ms
                # Cancela timer do Tab simples
                if self.tab_timer_id is not None:
                    GLib.source_remove(self.tab_timer_id)
                    self.tab_timer_id = None
                print("Double Tab detectado!")
                self.last_tab_time = 0
            else:
                # Agenda um Tab simples daqui 300 ms
                self.last_tab_time = now
                self.tab_timer_id = GLib.timeout_add(200, self._simple_Tab)
            return True
        return False
    
    def _simple_Tab(self):
        print("simple tab")
        self.tab_timer_id = None
        return False  # remove o timer
