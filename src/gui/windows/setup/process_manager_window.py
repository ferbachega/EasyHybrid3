import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
from gi.repository import GdkPixbuf
from gui.widgets.custom_widgets  import get_colorful_square_pixel_buffer
from gui.windows.setup.windows_and_dialogs import TextWindow
import os, sys, time


class ProcessManagerWindow(Gtk.Window):
    def __init__(self, main = None):
        self.main                = main
        self.p_session           = main.p_session
        self.home                = main.home
        
        self.Visible             =  False
        
        # Criar ListStore (modelo da TreeView)
        self.liststore = Gtk.ListStore(str,               #0.system name
                                       str,               #1.job type
                                       str,               #2.potential
                                       str,               #3.start
                                       str,               #4.ended
                                       str,               #5.status
                                       GdkPixbuf.Pixbuf,  #6.color 
                                       int,               #7.e_id
                                       int                #8.e_step_counter
                                       )                
        
        
        #self.liststore.append(["Sys001", "LigandA", "MD", "00:00", "01:00", "Running"])
        #self.liststore.append(["Sys002", "LigandB", "QM/MM", "01:10", "02:30", "Finished"])
        #self.liststore.append(["Sys003", "LigandC", "Docking", "02:40", "-", "Queued"])



    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:

            self.window = Gtk.Window(title="Process Manager")
            #super().__init__(title="Janela Auxiliar - Processos")
            self.window.set_default_size(600, 300)
            self.window.set_border_width(10)
            self.window.set_keep_above(True)


            self.build_treeview()
            self.build_popupmenu()
            

            # TreeView dentro de ScrolledWindow
            scrolled = Gtk.ScrolledWindow()
            scrolled.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
            scrolled.add(self.treeview)
            self.window.add(scrolled)
            self.window.connect('destroy', self.CloseWindow)
            self.treeview.connect("row-activated", self.on_row_activated)
            self.window.show_all()
            #self.methods_combo.set_active(0)
            self.Visible  = True

        else:
            self.window.present()
    
    def build_popupmenu (self):
        """ Function doc """
        # Criar popup menu
        self.popup_menu = Gtk.Menu()

        menu_item1 = Gtk.MenuItem(label="Rerun")
        menu_item1.connect("activate", self.rerun_job)
        self.popup_menu.append(menu_item1)

        menu_item2 = Gtk.MenuItem(label="Abort")
        menu_item2.connect("activate", self.on_stop_activate)
        self.popup_menu.append(menu_item2)

        menu_item3 = Gtk.MenuItem(label="Remove")
        menu_item3.connect("activate", self.on_remove_activate)
        self.popup_menu.append(menu_item3)

        self.popup_menu.show_all()

    def build_treeview (self):
        """ Function doc """
        # Criar TreeView
        self.treeview = Gtk.TreeView(model=self.liststore)
        
        # Adicionar colunas
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
            
        #titles = ["Type", "Potential", "Started", "End", "Status"]
        #for i, title in enumerate(titles):
        #    renderer = Gtk.CellRendererText()
        #    column = Gtk.TreeViewColumn(title, renderer, text=i+1)
        #    column.set_sort_column_id(i+1)
        #    self.treeview.append_column(column)

        # Adicionar suporte para menu de contexto
        self.treeview.connect("button-press-event", self.on_button_press_event)

    
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False


    def add_new_process (self,
                         system    = None, 
                         _type      = None, 
                         potential = None,
                         start     = None, 
                         end       = ' ', 
                         status    = None):
        """ Function doc """
        name  = ' '+ system.label
        sqr_color      = get_colorful_square_pixel_buffer(system)
        
        if start is not None:
            pass
        else:
            #----------------------------------------------------------------------------
            current_time   = time.time()
            formatted_time = time.strftime("%d/%m %H:%M:%S", time.localtime(current_time))
            start     = formatted_time
            #----------------------------------------------------------------------------
        
        treeiter = self.liststore.append([name, _type, potential, start, end, status, sqr_color,  system.e_id, -1 ])
        #self.set_time (treeiter, start = True, end = False)
        return treeiter


    def set_status (self, treeiter = None, status = 'Queued' ):
        """ Function doc """
        self.liststore[treeiter][5] = status

    
    def set_time (self, treeiter, start = False, end = False):
        """ Function doc """
        if start:
            current_time   = time.time()
            #formatted_time = time.strftime("%Y-%m-%d   %H:%M:%S", time.localtime(current_time))
            formatted_time = time.strftime("%d/%m %H:%M:%S", time.localtime(current_time))
            self.liststore[treeiter][4] = formatted_time
        
        if end:
            current_time   = time.time()
            #formatted_time = time.strftime("%Y-%m-%d   %H:%M:%S", time.localtime(current_time))
            formatted_time = time.strftime("%d/%m %H:%M:%S", time.localtime(current_time))
            self.liststore[treeiter][4] = formatted_time

    
    def set_step_counter (self, treeiter, e_step_counter):
        """ Function doc """
        self.liststore[treeiter][8] = e_step_counter
        
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
        print(system.e_job_history[step_counter])
        
        logfile = system.e_job_history[step_counter]['logfile']
        data = open(logfile, 'r')
        data = data.read()
        
        textwindow = TextWindow(data)
        
        
    def rerun_job(self, widget):
        return False
        model, treeiter = self.treeview.get_selection().get_selected()
        if treeiter:
            model[treeiter][5] = "Running"
            print(f"Iniciando: {model[treeiter][1]}")

    def on_stop_activate(self, widget):
        model, treeiter = self.treeview.get_selection().get_selected()
        if treeiter:
            e_id = model[treeiter][7]
            process = self.main.p_session.process_pool[e_id][1]
            
            # ... abort at some later point:
            process.terminate()    # requests immediate termination
            process.join(timeout=5)  # “awaits cleanup for up to 5 seconds”
            
            
            model[treeiter][5] = "Aborted"
            #print(f"Parando: {model[treeiter][1]}")
            
        
        
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
