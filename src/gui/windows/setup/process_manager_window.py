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
from gi.repository import Gtk, Gdk
from gi.repository import GdkPixbuf
from gui.widgets.custom_widgets  import get_colorful_square_pixel_buffer
from gui.windows.setup.windows_and_dialogs import TextWindow
import os, sys, time
from pprint import pprint

class ProcessManagerWindow(Gtk.Window):
    """
    A GTK window that provides an interface to manage simulation or computational jobs.

    This window displays a history of submitted jobs, their metadata, and allows
    interaction through a popup menu (e.g., rerun, abort, remove, clear).
    """
    def __init__(self, main = None):
        """
        Initialize the ProcessManagerWindow.

        Args:
            main (MainWindow): Reference to the main EasyHybrid window.
                Provides access to session objects and job history.
        """
        self.main                = main
        self.p_session           = main.p_session
        self.home                = main.home
        self.Visible             = False

    def open_window (self):
        """
        Open the process manager window.

        If the window is not already visible, it creates a new GTK window,
        sets up the treeview and popup menu, and shows the interface.
        If the window is already open, it is brought to the front.
        """
        if not self.Visible:
            self.window = Gtk.Window(title="Process Manager")
            self.window.set_default_size(600, 300)
            self.window.set_border_width(10)
            self.window.set_keep_above(True)

            # Build widgets
            self.build_treeview()
            self.build_popupmenu()

            # Add treeview inside a scrolled window
            scrolled = Gtk.ScrolledWindow()
            scrolled.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
            scrolled.add(self.treeview)
            self.window.add(scrolled)

            # Connect signals
            self.window.connect("destroy", self.close_window)
            self.treeview.connect("row-activated", self.on_row_activated)

            self.window.show_all()
            self.Visible = True
        else:
            # Bring existing window to front
            self.window.present()

    def close_window (self, button, data  = None):
        """
        Close the process manager window and update visibility state.

        Args:
            button (Gtk.Widget): The widget that triggered the close.
            data (Any, optional): Additional signal data.
        """
        self.window.destroy()
        self.Visible    =  False
    
    def build_popupmenu (self):
        """
        Build the context popup menu for job management.

        Provides the following options:
            - Rerun a job
            - Abort a job
            - Remove a job from the list
            - Clear the entire job list
        """
        # Criar popup menu
        self.popup_menu = Gtk.Menu()
        
        # Rerun job
        menu_item1 = Gtk.MenuItem(label="Rerun")
        menu_item1.connect("activate", self.rerun_job)
        self.popup_menu.append(menu_item1)

        # Abort job
        menu_item2 = Gtk.MenuItem(label="Abort")
        menu_item2.connect("activate", self.on_stop_activate)
        self.popup_menu.append(menu_item2)
        
        # Remove job
        menu_item3 = Gtk.MenuItem(label="Remove")
        menu_item3.connect("activate", self.on_remove_activate)
        self.popup_menu.append(menu_item3)
        
        # Clear all jobs
        menu_item4 = Gtk.MenuItem(label="Clear List")
        menu_item4.connect("activate", self.on_clear_list)
        self.popup_menu.append(menu_item4)
        
        self.popup_menu.show_all()

    def build_treeview (self):
        """
        Build the TreeView widget for displaying job history.

        The TreeView is populated using `main.job_history_liststore` and
        includes the following columns:
            - Job ID
            - System Name (with pixbuf icon)
            - Job Type
            - Job Index
            - Potential
            - Status
            - Start Time
            - End Time
        """
        
        self.treeview = Gtk.TreeView(model = self.main.job_history_liststore)
        
        #---------------------------------------------------------------
        renderer_int = Gtk.CellRendererText()
        column_int = Gtk.TreeViewColumn("id", renderer_int, text=7)
        column_int.set_sort_column_id(0)
        self.treeview.append_column(column_int)
        
        #------------------ system name --------------------------------
        renderer_pixbuf = Gtk.CellRendererPixbuf()
        renderer_text   = Gtk.CellRendererText()
        column_text     = Gtk.TreeViewColumn("System Name")#, renderer_text, text=2)
        
        column_text.pack_start(renderer_pixbuf, False)
        column_text.add_attribute(renderer_pixbuf, "pixbuf", 6)
        column_text.pack_start(renderer_text, True)
        column_text.add_attribute(renderer_text, "text", 0)
        column_text.set_sort_column_id(0)
        self.treeview.append_column(column_text)
        #---------------------------------------------------------------
        
        columns = {"Job Type"  : 1, 
                   'job'       : 8,
                   "Potential" : 2, 
                   "Status"    : 5,
                   "Started"   : 3, 
                   "Ended"     : 4  
                   #"Status"    : 
                   
                   }
        for title in columns.keys():
            renderer = Gtk.CellRendererText()
            i = columns[title]
            column = Gtk.TreeViewColumn(title, renderer, text=i)
            column.set_sort_column_id(i)
            self.treeview.append_column(column)

        # Enable context menu on right-click
        self.treeview.connect("button-press-event", self.on_button_press_event)

    def add_new_process (self,
                         system    = None, 
                         _type      = None, 
                         potential = None,
                         start     = None, 
                         end       = ' ', 
                         status    = None,
                         step      = None):
        
        """Add a new process entry to the process history list.

        This method creates a new row in the job history ListStore, associated 
        with the provided `system`. If no start time is given, the current time 
        will be automatically generated. The function also updates the system's 
        job history and assigns a step counter.

        Args:
            system (object, optional): The system object containing job history 
                and identifier information.
            _type (str, optional): Type of the process (e.g., minimization, 
                dynamics).
            potential (str, optional): The potential used in the process.
            start (str, optional): Start time of the process. If not provided, 
                the current time is used.
            end (str, optional): End time of the process. Defaults to a single 
                space (' ').
            status (str, optional): Current status of the process (e.g., running, 
                finished).
            step (int, optional): Step counter of the process. If not provided, 
                the value from `system.e_step_counter` is used.

        Returns:
            Gtk.TreeIter: An iterator pointing to the newly added row in 
            `job_history_liststore`.
        """
        
        # Prepare process name and graphical square icon.
        name  = ' '+ system.label
        sqr_color = get_colorful_square_pixel_buffer(system)
        
        # If no start time was provided, generate current formatted timestamp.
        if start is not None:
            pass
        else:
            #----------------------------------------------------------------------------
            current_time   = time.time()
            formatted_time = time.strftime("%d/%m %H:%M:%S", time.localtime(current_time))
            start     = formatted_time
            #----------------------------------------------------------------------------
        
        #when the step is external - from load EH session.
        # Use provided step counter or fallback to system default.
        if step is not None:
            step_counter = step
        else:
            step_counter = system.e_step_counter
        
        # Update system job history with the start time.
        system.e_job_history[step_counter]['started'] = start
        
        # Append new process entry to the job history ListStore.
        treeiter = self.main.job_history_liststore.append([
            name, _type, potential, start, end, status, sqr_color,  
            system.e_id, step_counter 
            ])
        
        return treeiter

    def set_status (self, treeiter = None, status = 'Queued' ):
        """ Function doc """
        self.main.job_history_liststore[treeiter][5] = status

    def set_time (self, treeiter, start = False, end = False):
        """ Function doc """
        e_id = self.main.job_history_liststore[treeiter][7]
        step_counter = self.main.job_history_liststore[treeiter][8]
        try:
            system = self.p_session.psystem[e_id]
        except:
            system = False
        
        if start:
            current_time   = time.time()
            #formatted_time = time.strftime("%Y-%m-%d   %H:%M:%S", time.localtime(current_time))
            formatted_time = time.strftime("%d/%m %H:%M:%S", time.localtime(current_time))
            self.main.job_history_liststore[treeiter][4] = formatted_time
            if system:
                system.e_job_history[step_counter]['started'] = formatted_time
            
            
        if end:
            current_time   = time.time()
            #formatted_time = time.strftime("%Y-%m-%d   %H:%M:%S", time.localtime(current_time))
            formatted_time = time.strftime("%d/%m %H:%M:%S", time.localtime(current_time))
            self.main.job_history_liststore[treeiter][4] = formatted_time
            if system:
                system.e_job_history[step_counter]['ended'] = formatted_time
    
    def set_step_counter (self, treeiter, e_step_counter):
        """ Function doc """
        self.main.job_history_liststore[treeiter][8] = e_step_counter
        
    def on_button_press_event(self, widget, event):
        if event.type == Gdk.EventType.BUTTON_PRESS and event.button == 3:
            # Seleciona a linha clicada com o botão direito
            path_info = widget.get_path_at_pos(int(event.x), int(event.y))
            if path_info is not None:
                path, col, cell_x, cell_y = path_info
                widget.grab_focus()
                widget.set_cursor(path, col, 0)
                self.popup_menu.popup_at_pointer(event)
            return True
        return False

    def on_row_activated(self, treeview, path, column):
        model = treeview.get_model()
        iter = model.get_iter(path)
        nome = model[iter][0]
        pid = model[iter][1]
        e_id = model[iter][7]
        step_counter = model[iter][8]
        print(f"Double click on: {e_id} {nome} {step_counter} (PID={pid})")
        system = self.p_session.psystem[e_id]
        system.e_job_history[step_counter]['e_id'] = e_id
        pprint(system.e_job_history[step_counter])
        
        
        #print('system:', system.e_id, system.e_job_history[step_counter]['backup_parameters']['system'].e_id)
        
        logfile = system.e_job_history[step_counter]['logfile']
        data = open(logfile, 'r')
        data = data.read()
        textwindow = TextWindow(data)
        
    def rerun_job(self, widget):
        #return False
        
        model, treeiter = self.treeview.get_selection().get_selected()
        
        #iter = model.get_iter(path)
        nome = model[treeiter][0]
        pid = model[treeiter][1]
        e_id = model[treeiter][7]
        step_counter = model[treeiter][8]
        
        self.main.main_treeview.set_active_system(e_id)
        system = self.p_session.psystem[e_id]
        
        parameters = system.e_job_history[step_counter]['backup_parameters']
        parameters['e_id'] = system.e_id
        
        if parameters['simulation_type'] == 'Geometry_Optimization':
            if self.main.geometry_optimization_window.Visible:
                self.main.geometry_optimization_window.close_window(None, None)
            self.main.geometry_optimization_window.open_window()
            self.main.geometry_optimization_window.restore_the_parameters_to_the_window (parameters)
        
        if parameters['simulation_type'] == 'Molecular_Dynamics':
            if self.main.molecular_dynamics_window.Visible:
                self.main.molecular_dynamics_window.close_window(None, None)
            self.main.molecular_dynamics_window.open_window()
            self.main.molecular_dynamics_window.restore_the_parameters_to_the_window (parameters)
            
        
        if parameters['simulation_type'] == 'Relaxed_Surface_Scan':
            if self.main.PES_scan_window.Visible:
                self.main.PES_scan_window.close_window(None, None)
            self.main.PES_scan_window.open_window()
            self.main.PES_scan_window.restore_the_parameters_to_the_window (parameters)
            
        
        if parameters['simulation_type'] == 'Umbrella_Sampling':
            if self.main.umbrella_sampling_window.Visible:
                self.main.umbrella_sampling_window.close_window(None, None)
            self.main.umbrella_sampling_window.open_window()
            self.main.umbrella_sampling_window.restore_the_parameters_to_the_window (parameters)
            
        
        
        if parameters['simulation_type'] == 'Nudged_Elastic_Band':
            if self.main.chain_of_states_opt_window.Visible:
                self.main.chain_of_states_opt_window.close_window(None, None)
            self.main.chain_of_states_opt_window.open_window()
            self.main.chain_of_states_opt_window.restore_the_parameters_to_the_window (parameters)
            
        
        
        #pprint(system.e_job_history[step_counter])
        
        
        
        
        #if treeiter:
        #    model[treeiter][5] = "Running"
        #    pprint(f"Iniciando: {model[treeiter][1]}")

    def on_stop_activate(self, widget):
        model, treeiter = self.treeview.get_selection().get_selected()
        if treeiter:
            e_id = model[treeiter][7]
            process = self.main.p_session.process_pool[e_id][1]
            
            # ... abort at some later point:
            process.terminate()    # requests immediate termination
            process.join(timeout=5)  # “awaits cleanup for up to 5 seconds”
            model[treeiter][5] = "Aborted"
            
            system       = self.p_session.psystem[e_id]
            step_counter = system.e_step_counter
            system.e_job_history[step_counter]['status'] = "Aborted"
            
            self.p_session.psystem[e_id].e_step_counter += 1
            self.set_time ( treeiter= treeiter,  end = True)
     
     
    def on_clear_list(self, widget):
        self.main.job_history_liststore.clear()
        #model, treeiter = self.treeview.get_selection().get_selected()
        #if treeiter:
        #    if model[treeiter][5] == "Running...":
        #        
        #        pass
        #    else:
        #        #print(f"Removing: {model[treeiter][1]}")
        #        model.remove(treeiter)
            
    def on_remove_activate(self, widget):
        model, treeiter = self.treeview.get_selection().get_selected()
        if treeiter:
            if model[treeiter][5] == "Running...":
                
                pass
            else:
                #print(f"Removing: {model[treeiter][1]}")
                model.remove(treeiter)

    def update (self, parameters = None):
        """ Function doc """
        pass

    def clear_liststore_job_history (self, clear = True):
        """ Function doc """
        if clear:
            self.main.job_history_liststore.clear()
    
    def build_liststore_from_job_history (self, clear = True):
        """Rebuild the job history ListStore from the stored session data.

        This method iterates over all systems in the current session (`p_session`),
        accesses their job history, validates job attributes, and reconstructs the
        corresponding entries in the `job_history_liststore`. Each valid job is 
        added to the ListStore via `add_new_process`.

        Args:
            clear (bool, optional): Whether to clear the current ListStore before 
                rebuilding. Defaults to True.

        Returns:
            bool: Returns False if an invalid job entry is found (e.g., with 
            unknown step counter). Otherwise, returns None.
        """
        # TODO: Implement optional clearing logic for job_history_liststore if needed.
        # Example:
        # if clear and hasattr(self, "treeview"):
        #     self.main.job_history_liststore.clear()
        
        # Iterate over all systems in the session.
        for e_id in self.p_session.psystem.keys():
            system    = self.p_session.psystem[e_id]
            
            
            print (len(self.p_session.psystem.keys()),
                   len(system.e_job_history.keys()),
                   system.label
                  
                  )
            # there might be several jobs in system.e_job_history dict for each system
            for job_num in system.e_job_history.keys():
                
                #.the list of parameters to reload the liststore
                keys = [
                    "status",
                    "step_counter",
                    "simulation_type",
                    "potential",
                    "status",
                    "started",
                    "ended",
                ]
                
                job = system.e_job_history[job_num]
                
                #. checking de parameters
                for key in keys:
                    if key in job.keys():
                        pass
                    else:
                        job[key] = 'Unk'
                
                status    = job['status']
                job_num   = job['step_counter']
                job_type  = job['simulation_type']
                potential = job['potential']
               
                started = job['started']
                ended   = job['ended']
                
                name  = system.label  # System label (system name) for display purposes.
                
                #'''
                #print('\n\n')
                #print(status, job_num,job_type, potential, started,  ended)
                #print('\n\n')
                #'''
                
                # Validate job number (step counter).
                if job_num == 'Unk':
                    return False
 
                # Append job to the ListStore.
                self.add_new_process ( 
                             system    = system, 
                             _type     = job_type, 
                             potential = potential,
                             start     = started, 
                             end       = ended, 
                             status    = status,
                             step      = job_num)
 
