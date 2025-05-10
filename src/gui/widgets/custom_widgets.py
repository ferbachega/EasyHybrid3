import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
from gi.repository import Gdk
from gi.repository import GdkPixbuf
import os

import threading
import time

from util.geometric_analysis            import get_distance 
from util.geometric_analysis            import get_dihedral 
from util.geometric_analysis            import get_angle 
#from util.periodic_table import atomic_dic 

from util.sequence_plot import GtkSequenceViewer


import numpy as np

HOME        = os.environ.get('HOME')


def getColouredPixmap( r, g, b, a=255 ):
    """ Given components, return a colour swatch pixmap """
    CHANNEL_BITS=8
    WIDTH=10
    HEIGHT=10
    swatch = GdkPixbuf.Pixbuf.new( GdkPixbuf.Colorspace.RGB, True, CHANNEL_BITS, WIDTH, HEIGHT ) 
    swatch.fill( (r<<24) | (g<<16) | (b<<8) | a ) # RGBA
    return swatch


def get_colorful_square_pixel_buffer (system = None,  atomic_symbol = 'C'):
    """ Function doc """
    if system is not None:
        color        =  system.e_color_palette[atomic_symbol]
        res_color    = [int(color[0]*255),int(color[1]*255),int(color[2]*255)] 
        pixelbuffer  =  getColouredPixmap( res_color[0], res_color[1], res_color[2] )
    
    else:
        res_color    = [255,255,255] 
        pixelbuffer  =  getColouredPixmap( res_color[0], res_color[1], res_color[2] )
    return pixelbuffer


class VismolSelectionTypeBox(Gtk.Box):
    """ Class doc """
    
    def __init__ (self, vm_session = None):
        """ Class initialiser """
        Gtk.Box.__init__(self)
        #self.set_orientation(Gtk.Orientation.VERTICAL)
        #self.set_orientation(Gtk.Orientation.HORIZONTAL)
        #self.set_spacing(5)
        self.box           = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
        self.vm_session = vm_session
        #combobox
        
        self.Vismol_selection_modes_ListStore = Gtk.ListStore(str)
        data = ['atom'    , 
                'residue' ,
                'chain'   , 
                'molecule',
                'c-alpha' ,
                'solvent' ,
                #'protein' , 
                #'C alpha' ,
                #'solvent' ,
                #'atom name',
                #'element',
                ]
        for i in data:
            #print (i)
            self.Vismol_selection_modes_ListStore.append([i])
            
        self.combobox_selection_type = Gtk.ComboBox.new_with_model(self.Vismol_selection_modes_ListStore)
        
        
        self.combobox_selection_type.set_model(self.Vismol_selection_modes_ListStore)
        
        self.renderer_text = Gtk.CellRendererText()
        self.combobox_selection_type.pack_start(self.renderer_text, True)
        self.combobox_selection_type.add_attribute(self.renderer_text, "text", 0)
        
        self.combobox_selection_type.connect('changed', self.on_combobox_selection_type)
        #
        
        #labels        
        self.label_selecting_by = Gtk.Label('selecting by: ')
        
        #toggle_button
        self.toggle_button_selecting_mode = Gtk.ToggleButton('Viewing')
        self.toggle_button_selecting_mode.connect('clicked', self.on_toggle_button_selecting_mode)
        
        
        #Atom name entry box
        self.entry_atom_names = None#Gtk.Entry()
        self.entry_elements   = None
        
        # Packing 
        self.box.pack_start(self.toggle_button_selecting_mode, False, False, 0)
        self.box.pack_start(self.label_selecting_by          , False, False, 0)
        self.box.pack_start(self.combobox_selection_type     , False, False, 0)
        
        self.combobox_selection_type.set_active(1)
        self.box.show_all()
    
    def on_entry_data_change (self, widget):
        """ Function doc """
        
        string = widget.get_text()
        
        if widget == self.entry_atom_names:
            self.vm_session.selections[self.vm_session.current_selection].selected_atom_names_list = []
            selected_atom_names_list = self.vm_session.selections[self.vm_session.current_selection].selected_atom_names_list 
        
            keys = string.split('+')
            for key in keys:
                selected_atom_names_list.append(key.strip())
            
            print ('entry_atom_names', self.vm_session.selections[self.vm_session.current_selection].selected_atom_names_list)
        
        
        
        elif widget == self.entry_elements:
            self.vm_session.selections[self.vm_session.current_selection].selected_element_list = []
            selected_element_list = self.vm_session.selections[self.vm_session.current_selection].selected_element_list  
            
            keys = string.split('+')
            for key in keys:
                selected_element_list.append(key.strip())
            
            print ('entry_elements', self.vm_session.selections[self.vm_session.current_selection].selected_element_list)

        
        else: 
            pass
               
    def on_combobox_selection_type (self, combobox):
        """ Function doc """
        self.active = combobox.get_active()
        #print('on_combobox_selection_type', self.active)
        if self.active == 0:
            self.vm_session.viewing_selection_mode(sel_type = 'atom')
            #self.show_or_hide_entries (name = 'atom_names', show = False)
            #self.show_or_hide_entries (name = 'element', show = False)

        elif self.active == 1:
            self.vm_session.viewing_selection_mode(sel_type = 'residue')
            #self.show_or_hide_entries (name = 'atom_names', show = False)
            #self.show_or_hide_entries (name = 'element', show = False)
        
        elif self.active == 2:
            self.vm_session.viewing_selection_mode(sel_type = 'chain')
            #self.show_or_hide_entries (name = 'atom_names', show = False)
            #self.show_or_hide_entries (name = 'element', show = False)
        
        elif self.active == 3:
            self.vm_session.viewing_selection_mode(sel_type = 'molecule')
            #self.show_or_hide_entries (name = 'atom_names', show = False)
            #self.show_or_hide_entries (name = 'element', show = False)
        
        elif self.active == 4:
            self.vm_session.viewing_selection_mode(sel_type = 'c-alpha')
        
        elif self.active == 5:
            self.vm_session.viewing_selection_mode(sel_type = 'solvent')
            #self.show_or_hide_entries (name = 'atom_names', show = False)
            #self.show_or_hide_entries (name = 'element', show = False)
        
        #elif self.active == 5:
        #    self.vm_session.viewing_selection_mode(sel_type = 'solvent')
        #    #self.show_or_hide_entries (name = 'atom_names', show = False)
        #    #self.show_or_hide_entries (name = 'element', show = False)
        #
        #elif self.active == 6:
        #    self.vm_session.viewing_selection_mode(sel_type = 'atom name')
        #    #self.show_or_hide_entries (name = 'atom_names', show = True)
        #    #self.show_or_hide_entries (name = 'element', show = False)
        #
        #elif self.active == 7:
        #    self.vm_session.viewing_selection_mode(sel_type = 'element')
        #    #self.show_or_hide_entries (name = 'atom_names', show = False)
        #    #self.show_or_hide_entries (name = 'element', show = True)
        #        
        else:pass
        
    def change_sel_type_in_combobox (self, sel_type):
        """ Function doc """
        #print('change_sel_type_in_combobox', sel_type)
        if sel_type == 'atom':
            self.combobox_selection_type.set_active(0)
        
        elif sel_type == 'residue':
            self.combobox_selection_type.set_active(1)
        
        elif sel_type == 'chain':
            self.combobox_selection_type.set_active(2)
        
        elif sel_type == 'molecule':
            self.combobox_selection_type.set_active(3)
        
        elif sel_type == 'C alpha':
            self.combobox_selection_type.set_active(4)
        
        elif sel_type == 'solvent':
            self.combobox_selection_type.set_active(5)
        
        elif sel_type == 'atom name':
            self.combobox_selection_type.set_active(6)
        
        elif sel_type == 'element':
            self.combobox_selection_type.set_active(7)
        
        else: pass
        
    def on_toggle_button_selecting_mode (self, button):
        """ Function doc """
        if button.get_active():
            state = "on"
            
            self.vm_session.picking_selection_mode = True
            button.set_label('Picking')
            #print(self.combobox_selection_type.get_active())
            self.vm_session._selection_function_set (None)
            self.vm_session.vm_glcore.queue_draw()
            
            #self.vm_session._picking_selection_mode = True
            
            self.combobox_selection_type.set_sensitive(False)
            self.label_selecting_by.set_sensitive(False)
            
            circle = Gdk.Cursor(Gdk.CursorType.CROSSHAIR)
            button.get_window().set_cursor(circle)
            
        else:
            state = "off"
            self.vm_session.picking_selection_mode = False
            button.set_label('Viewing')
            self.vm_session.vm_glcore.queue_draw()
            
            self.combobox_selection_type.set_sensitive(True)
            self.label_selecting_by.set_sensitive(True)
            
            
            circle = Gdk.Cursor(Gdk.CursorType.ARROW)
            button.get_window().set_cursor(circle)
            
    def change_toggle_button_selecting_mode_status (self, status = False):
        """ Function doc """
        self.toggle_button_selecting_mode.set_active(status)

    def update (self):
        """ Function doc """
        print('VismolSelectionTypeBox update')


class SaveTrajectoryBox:
    """ Class doc """
    
    def __init__ (self, parent, home):
        """ Class initialiser """
        
        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(home,'src/gui/widgets/trajectory_box.glade'))
        self.builder.connect_signals(self)
        
        self.box = self.builder.get_object('trajectory_box')
        
        self.folder_chooser_button = FolderChooserButton(main =  parent, home = home)
        self.builder.get_object('folder_chooser_box').pack_start(self.folder_chooser_button.btn, True, True, 0)
        self.folder_chooser_button.btn.set_sensitive(False)
                
        '''--------------------------------------------------------------------------------------------'''
        self.format_store = Gtk.ListStore(str)
        formats = [
                    "pDynamo / pkl" ,
                    #"amber / crd"   ,
                    #"charmm / dcd"  ,
                    #"xyz"           ,
            ]
        for format in formats:
            self.format_store.append([format])
            #print (format)
        self.formats_combo = self.builder.get_object('combobox_format')
        self.formats_combo.set_model(self.format_store)
        #self.formats_combo.connect("changed", self.on_name_combo_changed)
        self.formats_combo.set_model(self.format_store)
        
        renderer_text = Gtk.CellRendererText()
        self.formats_combo.pack_start(renderer_text, True)
        self.formats_combo.add_attribute(renderer_text, "text", 0)
        '''--------------------------------------------------------------------------------------------'''
        self.formats_combo.set_active(0)
        self.set_folder(None)
        #simParameters["trajectory_name"] = self.save_trajectory_box.builder.get_object('entry_trajectory_name').get_text()
    #====================================================================================
    def on_toggle_save_checkbox (self, widget):
        """ Function doc """
        if self.builder.get_object('checkbox_save_traj').get_active():
            self.builder.get_object('entry_trajectory_name').set_sensitive(True)
            
            self.builder.get_object('label_working_folder').set_sensitive(True)
            #self.builder.get_object('file_chooser_working_folder').set_sensitive(True)
            self.builder.get_object('label_format').set_sensitive(True)
            self.builder.get_object('combobox_format').set_sensitive(True)
            self.builder.get_object('label_trajectory_frequency').set_sensitive(True)
            self.builder.get_object('entry_trajectory_frequency').set_sensitive(True)
            self.folder_chooser_button.btn.set_sensitive(True)

        else:
            self.builder.get_object('entry_trajectory_name').set_sensitive(False)

            self.builder.get_object('label_working_folder').set_sensitive(False)
            #self.builder.get_object('file_chooser_working_folder').set_sensitive(False)
            self.builder.get_object('label_format').set_sensitive(False)
            self.builder.get_object('combobox_format').set_sensitive(False)
            self.builder.get_object('label_trajectory_frequency').set_sensitive(False)
            self.builder.get_object('entry_trajectory_frequency').set_sensitive(False)
            self.folder_chooser_button.btn.set_sensitive(False)
       
    
    def get_trajectory_frequency (self):
        """ Function doc """
        return int(self.builder.get_object('entry_trajectory_frequency').get_text())
    
    def get_format (self):
        return  self.formats_combo.get_active()
        
    def get_active (self):
        """ Function doc  """
        if self.builder.get_object('checkbox_save_traj').get_active():
            return True
        else:
            return False
    
    def set_active (self):
        """ Function doc  """
        self.builder.get_object('checkbox_save_traj').set_active(True) 

    #====================================================================================
    def get_folder (self):
        """ Function doc """
        return self.folder_chooser_button.get_folder()
    #====================================================================================    
    def set_folder (self, folder):
        """ Function doc """
        self.folder_chooser_button.set_folder(folder)
    
    def get_filename (self):
        """ Function doc """
        return self.builder.get_object('entry_trajectory_name').get_text()
    
    def set_filename (self, filename = 'filename'):
        """ Function doc """
        return self.builder.get_object('entry_trajectory_name').set_text(filename)


#=========================================================================================
class FolderChooserButton:
    """ Class doc """
    
    def __init__ (self, main = None, sel_type = 'folder', home  = None):
        """ Class initialiser """
        self.main     =  main
        self.btn      =  Gtk.Button()
        self.sel_type =  sel_type # file/folder
        
        #self.sel_type = 'file' # file/folder

        grid = Gtk.Grid ()
        grid.set_column_spacing (10)
        img = Gtk.Image()
        img.set_from_file(os.path.join(self.main.home,'src/gui/icons/icon_open.png'))
        self.label = Gtk.Label ('...')
        self.btn.connect('clicked', self.open_filechooser)

        grid.attach (img, 0, 0, 1, 1)
        grid.attach (self.label, 1, 0, 1, 1)
        grid.show_all ()    

        self.btn.add (grid)        
        #return self.btn 
        self.folder = None
        
        self.set_tooltop('The working directory is where EasyHybrid suggests saving any files generated by the simulation routines. This is a dynamic variable that can be easily changed at any time while using the interface.')
    #====================================================================================
    def set_folder (self, folder = '/home'):
        """ Function doc """
        if folder == None:
            self.folder = os.environ.get('HOME')
        else:
            self.folder = folder
        
        name = os.path.basename(self.folder )
        #print( name)
        self.label.set_text(name)
        return self.folder
        #self.main.pdynamo_session.systems[self.main.pdynamo_session.active_id]['working_folder'] = self.folder
    #====================================================================================
    def get_folder (self):
        """ Function doc """
        return self.folder
        
    #====================================================================================    
    def open_filechooser (self, parent = None):
        """ Function doc """
        
        if self.sel_type == 'folder':
            dialog = Gtk.FileChooserDialog(
                
                #title="Please choose a file", parent=window, action=Gtk.FileChooserAction.OPEN
                title="Please choose a folder", parent= self.main.window, action=Gtk.FileChooserAction.SELECT_FOLDER
            )
            dialog.set_select_multiple(True)
            dialog.add_buttons(
                Gtk.STOCK_CANCEL,
                Gtk.ResponseType.CANCEL,
                Gtk.STOCK_OPEN,
                Gtk.ResponseType.OK,
            )
        
        else:
            dialog = Gtk.FileChooserDialog(
                
                #title="Please choose a file", parent=window, action=Gtk.FileChooserAction.OPEN
                title="Please choose a file", parent= self.main.window, action=Gtk.FileChooserAction.OPEN
            )
            dialog.add_buttons(
                Gtk.STOCK_CANCEL,
                Gtk.ResponseType.CANCEL,
                Gtk.STOCK_OPEN,
                Gtk.ResponseType.OK,
            )
        
        
        
        if self.folder:
            dialog.set_current_folder(self.folder)
        else:
            dialog.set_current_folder(self.main.home)
            #self.add_filters(dialog)

        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            #print("Open clicked")
        
            #print("File selected: " + dialog.get_filename())
        
            #print(dialog.get_filename())
            folder = dialog.get_filename()
            
            #folder = dialog.get_filenames()
            #print('\n\n\n',folder,'\n\n\n')
            
            self.set_folder(folder =folder)
            #print(os.path.dirname( dialog.get_filename() ))

        
        elif response == Gtk.ResponseType.CANCEL:
            print("Cancel clicked")

        dialog.destroy()
       
    def set_tooltop (self, text = ''):
        """ Function doc """
        self.btn.set_tooltip_text(text)


class FileChooser:
    """ Class doc """
    
    def __init__ (self, main_window = None , easyhybrid = False):
        """ Class initialiser """
        self.main_window = main_window
        self.easyhybrid = easyhybrid 
    
    def open (self, select_multiple = False, filters = None):

        """ Function doc """
        #main = gtkmain
        main = None#self.main_window
        filename = None
        
        chooser = Gtk.FileChooserDialog("Open File...", main,0,
                                       (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                                        Gtk.STOCK_OPEN, Gtk.ResponseType.OK))
        #GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER
        #print (filters)
        
        if self.main_window.current_search_folder:
            chooser.set_current_folder(self.main_window.current_search_folder)
        else:
            pass
            
        
        if select_multiple:
            chooser.set_select_multiple(True)
            response = chooser.run()
            if response == Gtk.ResponseType.OK:
                filenames = chooser.get_filenames()
            chooser.destroy()
            
            # updating current_search_folder
            file_directory = os.path.dirname(filenames[0])
            self.main_window.current_search_folder = file_directory
            
            return filenames
            
            
            
        
        
        
        else:
            if filters:
                for _filter in filters:
                    chooser.add_filter(_filter)
            else:
                '''
                filter = Gtk.FileFilter()  
                filter.set_name("PKL files - *.pkl")

                filter.add_mime_type("PKL files")
                filter.add_pattern("*.pkl")
                #
                chooser.add_filter(filter)
                '''
                filter = Gtk.FileFilter()
                filter.set_name("All files")
                filter.add_pattern("*")
                #
                chooser.add_filter(filter)  
            response = chooser.run()
            if response == Gtk.ResponseType.OK:
                filename = chooser.get_filename()
                
                file_directory = os.path.dirname(filename)
                self.main_window.current_search_folder = file_directory
            
            chooser.destroy()
            return filename


class VismolTrajectoryFrame(Gtk.Frame):
    """ Class doc """
    
    def __init__ (self, vm_session = None):
        """ Class initialiser """
        self.vm_session = vm_session 
        self.frame_forward = True
        self.frame      = Gtk.Frame()
        #self.frame.set_shadow_type(Gtk.SHADOW_IN)
        self.frame.set_border_width(4)
        
        self.box        = Gtk.Box(orientation = Gtk.Orientation.VERTICAL, spacing = 6) 
        
        self.box.set_margin_top    (3)
        self.box.set_margin_bottom (3)
        self.box.set_margin_left   (3)
        self.box.set_margin_right  (3)
        
        self.value      = 0
        self.upper      = 100
        self.scale      = Gtk.Scale()
        self.running    = False
        #self.adjustment = Gtk.Adjustment(self.value, 1, 1, 0, 1, 0)
        self.adjustment     = Gtk.Adjustment(value         = self.value,
                                             lower         = 0,
                                             upper         = 100,
                                             step_increment= 1,
                                             page_increment= 1,
                                             page_size     = 1)
        
        
        
        self.scale.set_adjustment ( self.adjustment)
        self.scale.set_digits(0)
        #self.scale.connect("change_value", self.on_scaler_frame_change_change_value)
        self.scale.connect("value_changed", self.on_scaler_frame_value_changed)
        self.box.pack_start(self.scale, True, True, 0)
        
        self.vbox =  Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
        
        self.functions = [self.reverse, self.stop, self.play, self.forward ]
        c = 0
        for label in ['<<','#', '>','>>']:
            button = Gtk.Button(label)
            self.vbox.pack_start(button, True, True, 0)
            
            if self.functions[c]:
                button.connect("clicked", self.functions[c])
            c += 1 
            
        self.box.pack_start(self.vbox, True, True, 0)
        
        
        #----------------------------------------------------------------------------
        '''
        self.label2 =  Gtk.Label('Object:')
        self.combobox_vobjects = Gtk.ComboBox()
        #self.combobox_vobjects = Gtk.ComboBox.new_with_model(self.vm_session.Vismol_Objects_ListStore)
        
        self.combobox_vobjects.connect("changed", self.on_combobox_vobjects_changed)
        self.renderer_text = Gtk.CellRendererText()
        self.combobox_vobjects.pack_start(self.renderer_text, True)
        self.combobox_vobjects.add_attribute(self.renderer_text, "text", 1)
        '''
        #self.box.pack_start(self.combobox, True, True, 0)
        #----------------------------------------------------------------------------
        
        #----------------------------------------------------------------------------
        self.vbox2 = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
        self.label = Gtk.Label('FPS:')
        self.label_step = Gtk.Label('Step:')
        
        #self.entry = Gtk.Entry()
        #self.entry.set_text(str(25))
        
        self.fps_adjustment = Gtk.Adjustment(value          = 24 , 
                                             upper          = 10000, 
                                             step_increment = 1  , 
                                             page_increment = 10 )
        
        self.step_adjustment = Gtk.Adjustment(value          = 1 , 
                                             upper          = 10000, 
                                             step_increment = 1  , 
                                             page_increment = 10 )
        
        
        #self.fps_adjustment = Gtk.Adjustment(value         = float, 
        #                                     lower         = float, 
        #                                     upper         = float, 
        #                                     step_increment= float, 
        #                                     page_increment= float, 
        #                                     page_size     = float)
        self.entry = Gtk.SpinButton()
        self.entry.set_adjustment ( self.fps_adjustment)
        
        self.entry_step = Gtk.SpinButton()
        self.entry_step.set_adjustment ( self.step_adjustment)
        
        #self.vbox2.pack_start(self.label2, True, True, 0)
        #self.vbox2.pack_start(self.combobox_vobjects, True, True, 0)
        self.vbox2.pack_start(self.label, True, True, 0)
        self.vbox2.pack_start(self.entry, True, True, 0)
        
        self.vbox2.pack_start(self.label_step, True, True, 0)
        self.vbox2.pack_start(self.entry_step, True, True, 0)
        
        self.box.pack_start(self.vbox2, True, True, 0)
        #----------------------------------------------------------------------------

    def play (self, button):
        # Create and start a new thread when the button is clicked
        if self.running:
            pass
        
        else:
            #thread = threading.Thread(target=self.long_task)
            thread = threading.Thread(target=self.long_task_boomerang)
            thread.start()
            
            self.running = True
        
    def long_task_boomerang(self):
        self.stop_thread = False
        i = 0
        
        while self.stop_thread == False:
            if self.stop_thread:
                return
            
            if self.frame_forward:
            
                value = self.forward(step =  self.entry_step.get_value()) 
                #value = self.reverse(None) 

                if value >= self.upper-1:
                    #value = 0
                    self.frame_forward = False
                    self.scale.set_value(int(value))
                    self.vm_session.set_frame(int(value))
                #time.sleep(0.01)
                time.sleep(1/self.entry.get_value())
            
            else:
                value = self.reverse(step =  self.entry_step.get_value()) 
                if value <= 0:
                    self.frame_forward = True
                    self.scale.set_value(int(value))
                    self.vm_session.set_frame(int(value))
                time.sleep(1/self.entry.get_value())
    
    def long_task(self):
        # This flag is used to stop the thread
        self.stop_thread = False
        i = 0
        while self.stop_thread == False:
            if self.stop_thread:
                return
            value = self.forward(None) 
            #value = self.reverse(None) 

            if value >= self.upper-1:
                value = 0
                self.scale.set_value(int(value))
                self.vm_session.set_frame(int(value))
            #time.sleep(0.01)
            time.sleep(1/self.entry.get_value())
            
    def stop (self, button):
        """ Function doc """
        self.stop_thread = True
        self.running = False    

    def forward (self, button = None, step = 1):
        """ Function doc """
        value =  int(self.scale.get_value())
        value = value+step
        self.scale.set_value(int(value))
        self.vm_session.set_frame(int(value))
        #print(value)
        return value
        
        
    def reverse (self, button =  None, step = 1):
        """ Function doc """
        value = int(self.scale.get_value())
        #print(value)
        if value == 0:
            pass
        else:
            value = value-step
            if value < 0:
                value = 0
            else: 
                pass
                
        self.vm_session.set_frame(int(value))
        self.scale.set_value(value)
        return value
        #print(value)
    
    def get_box (self):
        """ Function doc """
        #self.add(self.box)
        return self.box
        #return self.frame
        
    def on_combobox_vobjects_changed (self, widget):
        """ Function doc """
        #print('\n\n',widget)
        #print('\n\n',widget.get_active())
        
        cb_index = widget.get_active()
        if cb_index in self.vm_session.vobjects_dic:
            self.VObj = self.vm_session.vobjects_dic[widget.get_active()]
            #self.VObj = self.vm_session.vobjects[widget.get_active()]
            number_of_frames = len(self.VObj.frames)
            self.scale.set_range(0, int(number_of_frames))
            self.scale.set_value(self.vm_session.get_frame())
        else:
            pass

    def on_scaler_frame_value_changed (self, hscale, text= None,  data=None):
        """ Function doc """
        value = hscale.get_value()
        pos   = hscale.get_value_pos ()
        self.vm_session.set_frame(int(value)) 
        self.scale.set_value(value)
        #print(value, pos)

    def update (self):
        """ Function doc """
        #print('VismolTrajectoryFrame update')
        #for index , vobject in self.vm_session.vobjects_dic.items():
        #last_obj = len(self.vm_session.vobjects) -1
        last_obj = len(self.vm_session.vobjects_dic.items()) -1
        self.combobox_vobjects.set_active(last_obj)
    
    def change_range (self, upper = 100):
        """ Function doc """
        #print('upper =', upper)
        self.adjustment     = Gtk.Adjustment(value         = self.value,
                                             lower         = 0,
                                             upper         = upper,
                                             step_increment= 1,
                                             page_increment= 1,
                                             page_size     = 1)
        self.scale.set_adjustment ( self.adjustment)
        self.scale.set_digits(0)
        self.upper = upper


class SystemComboBox(Gtk.ComboBox):
    """ Class doc """
    
    def __init__ (self, main = None, coord_combobox = False):
        """ Class initialiser """
        Gtk.ComboBox.__init__(self)
        
        self.main = main
        
        self.system_liststore = self.main.system_liststore
        self.set_model(self.system_liststore)
        
        #self.combobox_systems.connect("changed", self.on_combobox_systems_changed)
        
        renderer_pixbuf = Gtk.CellRendererPixbuf()
        self.pack_start(renderer_pixbuf, True)
        self.add_attribute(renderer_pixbuf, "pixbuf", 2)
        
        renderer_text2 = Gtk.CellRendererText()
        self.pack_start(renderer_text2, True)
        self.add_attribute(renderer_text2, "text", 0)

        #self.set_popup_fixed_width(100)
        self.connect("changed", self.on_change)
        self.coord_combobox = coord_combobox
    
    
    def get_system_id(self, widget = None):
        _, system_id, pixbuf = self._get_system()
        return system_id
    
    def get_system_name(self, widget = None):
        name, system_id, pixbuf = self._get_system()
        return name
        
        
    def _get_system(self, widget = None):
        """ Function doc """
        index = self.get_active()
        
        if index == -1:
            '''_id = -1 means no item inside the combobox'''
            return None, None, None
        
        else:    
            _, system_id, pixbuf = self.system_liststore[index]
            #print(_, system_id, pixbuf)
            return _, system_id, pixbuf


    def set_active_system (self, e_id = 0):
        """ Function doc """
        
        counter    = 0
        set_active = 0
        for key , system in self.main.p_session.psystem.items():
            if system:
                if e_id == int(key):
                    set_active = counter
                    counter += 1
                else:
                    counter += 1
            else:
                pass
        
        self.set_active(set_active)

    
    def on_change (self, widget):
        """ Function doc """
        #print('AQUI', self.coord_combobox)
        
        if self.coord_combobox:
            system_id = self.get_system_id()
            system  = self.main.p_session.get_system(system_id)
            
            if system_id is not None:
                self.coord_combobox.set_model(self.main.vobject_liststore_dict[system_id])
                #self.refresh_selection_liststore (system_id)            
                size  =  len(list(self.main.vobject_liststore_dict[system_id]))
                self.coord_combobox.set_active(size-1)
        else:
            pass

class CoordinatesComboBox(Gtk.ComboBox):
    """ Class doc """
    
    def __init__ (self, coordinates_liststore = None, system_combobox = None, pixbuf = True):
        """ Class initialiser """
        Gtk.ComboBox.__init__(self)
        
        self.coordinates_liststore = coordinates_liststore
        self.set_model(self.coordinates_liststore)
        #self.coordinates_combobox.connect("changed", self.on_self.coordinates_combobox_changed)
        
        if pixbuf:
            renderer_pixbuf = Gtk.CellRendererPixbuf()
            self.pack_start(renderer_pixbuf, True)
            self.add_attribute(renderer_pixbuf, "pixbuf", 3)
        
        renderer_text = Gtk.CellRendererText()
        self.pack_start(renderer_text, True)
        self.add_attribute(renderer_text, "text", 0)

    def get_vobject (self):
        """ Function doc """
        name, key1,  key2, pixbuf = self._get_coordinates_info()
        
        #self.vobject  = self.main.vm_session.vm_objects_dic[vobject_index]
     
        
    def get_vobject_id (self):
        """ Returns the id of the vobject (id number of the vobject - 
        generated in eSession.py / self._add_vismol_object() ) """
        name, key1,  key2, pixbuf  = self._get_coordinates_info ()
        return key1

    
    def get_system_id (self):
        """ Returns the id of the vobject (id number of the vobject - 
        generated in eSession.py / self._add_vismol_object() ) """
        name, key1,  key2, pixbuf  = self._get_coordinates_info ()
        return key2
        
        
    def _get_coordinates_info (self):
        """ Function doc """
        cb_id =  self.get_active()
        model = self.get_model()
        name, key1,  key2, pixbuf = model[cb_id]
        #print (name, key1,  key2)
        return name, key1,  key2, pixbuf


    def set_active_vobject (self, pos = 0):
        """ Function doc """
        #self.coordinates_liststore[cb_id]
        if pos ==-1:
            n = len(self.get_model())
            #print(776, n)
            self.set_active(n-1)
        else:
            self.set_active(pos)

    
    
    def _starting_coordinates_model_update (self, init = False):
        """ Function doc """
        #------------------------------------------------------------------------------------
        '''The combobox accesses, according to the id of the active system, 
        listostore of the dictionary object_list state_dict'''
        if self.Visible:

            e_id = self.main.p_session.active_id 
            self.combobox_starting_coordinates.set_model(self.main.vobject_liststore_dict[e_id])
            #------------------------------------------------------------------------------------
            size = len(self.main.vobject_liststore_dict[e_id])
            self.combobox_starting_coordinates.set_active(size-1)
            #------------------------------------------------------------------------------------
        else:
            if init:
                e_id = self.main.p_session.active_id 
                self.combobox_starting_coordinates.set_model(self.main.vobject_liststore_dict[e_id])
                #------------------------------------------------------------------------------------
                size = len(self.main.vobject_liststore_dict[e_id])
                self.combobox_starting_coordinates.set_active(size-1)
                #------------------------------------------------------------------------------------
            else:
                pass




def compute_sigma_a1_a3 (vobject, index1, index3):

    """ example:
        pk1 ---> pk2 ---> pk3
         N  ---   H  ---  O	    
         
         where H is the moving atom
         calculation only includes N and O ! 
    """
    periodic_table = vobject.vm_session.periodic_table
    
    atom1 = vobject.atoms[index1]
    atom3 = vobject.atoms[index3]
    
    symbol1 = atom1.symbol
    symbol3 = atom3.symbol    
    
    #mass1 = atomic_dic[symbol1][4]
    #mass3 = atomic_dic[symbol3][4]
    mass1 = periodic_table.get_atomic_mass(symbol1)
    mass3 = periodic_table.get_atomic_mass(symbol3)
    
    print(atom1.name, symbol1, mass1)
    print(atom3.name, symbol3, mass3)

    ##pk1_name
    ##pk3_name
    #mass1 = atomic_dic[pk1_name][4]
    #mass3 = atomic_dic[pk3_name][4]
    #
    sigma_pk1_pk3 =  mass1/(mass1+mass3)
    print ("sigma_pk1_pk3: ",sigma_pk1_pk3)
    #
    sigma_pk3_pk1 =  mass3/(mass1+mass3)
    sigma_pk3_pk1 = sigma_pk3_pk1*-1
    #
    print ("sigma_pk3_pk1: ", sigma_pk3_pk1)
    return(sigma_pk1_pk3, sigma_pk3_pk1)


texto_d1   = "\n\n                       -- simple-distance --\n\nFor simple-distance, select two using picking mode\nfollowing the diagram:\n\n   R                    R\n    \                  /\n     A1--A2  . . . . A3\n    /                  \ \n   R                    R\n         ^            ^\n         |            |\n        pk1  . . . . pk2\n                d1\n"
texto_d2d1 = "\n                       -- multiple-distance --\n\nFor multiple-distance, select three atoms  using picking mode\nfollowing the diagram:\n\n   R                    R\n    \                  /\n     A1--A2  . . . . A3\n    /                  \ \n   R                    R\n     ^   ^            ^\n     |   |            |\n    pk1-pk2  . . . . pk3\n       d1       d2\n"


class ReactionCoordinateBox(Gtk.Box):
    """ Class doc """
    
    def __init__ (self, main = None, mode = 0 ):
        """ Class initialiser """
        Gtk.Box.__init__(self)
        
        
        if mode == 1:
            xml = '''
<?xml version="1.0" encoding="UTF-8"?>
<!-- Generated with glade 3.22.2 -->
<interface>
  <requires lib="gtk+" version="3.20"/>
  <object class="GtkBox" id="reaction_coordinate_box">
    <property name="visible">True</property>
    <property name="can_focus">False</property>
    <property name="orientation">vertical</property>
    <property name="spacing">4</property>
    <child>
      <object class="GtkBox">
        <property name="visible">True</property>
        <property name="can_focus">False</property>
        <property name="spacing">6</property>
        <child>
          <object class="GtkLabel" id="label_RC_type">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Coordinate Type:</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="expand">False</property>
            <property name="fill">True</property>
            <property name="position">0</property>
          </packing>
        </child>
        <child>
          <object class="GtkComboBox" id="combobox_reaction_coord1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
          </object>
          <packing>
            <property name="expand">True</property>
            <property name="fill">True</property>
            <property name="position">1</property>
          </packing>
        </child>
      </object>
      <packing>
        <property name="expand">False</property>
        <property name="fill">True</property>
        <property name="position">0</property>
      </packing>
    </child>
    <child>
      <object class="GtkGrid">
        <property name="visible">True</property>
        <property name="can_focus">False</property>
        <property name="row_spacing">5</property>
        <property name="column_spacing">10</property>
        <child>
          <object class="GtkLabel" id="label_atom4_coord1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Atom 4(index):</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">0</property>
            <property name="top_attach">3</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_atom3_index_coord1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
          </object>
          <packing>
            <property name="left_attach">1</property>
            <property name="top_attach">2</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_atom2_coord1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Atom 2(index):</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">0</property>
            <property name="top_attach">1</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_atom1_coord1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Atom 1(index):</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">0</property>
            <property name="top_attach">0</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_atom3_coord1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Atom 3(index):</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">0</property>
            <property name="top_attach">2</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_name3_coord1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Name:</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">2</property>
            <property name="top_attach">2</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_name4_coord1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Name:</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">2</property>
            <property name="top_attach">3</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_name2_coord1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Name:</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">2</property>
            <property name="top_attach">1</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_name1_coord1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Name:</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">2</property>
            <property name="top_attach">0</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_atom3_name_coord1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
          </object>
          <packing>
            <property name="left_attach">3</property>
            <property name="top_attach">2</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_atom2_index_coord1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
          </object>
          <packing>
            <property name="left_attach">1</property>
            <property name="top_attach">1</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_atom1_index_coord1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
          </object>
          <packing>
            <property name="left_attach">1</property>
            <property name="top_attach">0</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_atom4_index_coord1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
          </object>
          <packing>
            <property name="left_attach">1</property>
            <property name="top_attach">3</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_atom4_name_coord1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
          </object>
          <packing>
            <property name="left_attach">3</property>
            <property name="top_attach">3</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_atom2_name_coord1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
          </object>
          <packing>
            <property name="left_attach">3</property>
            <property name="top_attach">1</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_atom1_name_coord1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
          </object>
          <packing>
            <property name="left_attach">3</property>
            <property name="top_attach">0</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_step_size">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Step Size:</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">4</property>
            <property name="top_attach">0</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_num_of_steps">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Number of Steps:</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">4</property>
            <property name="top_attach">1</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_force">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Force Constante:</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">4</property>
            <property name="top_attach">2</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_initial_distance_angle_dihedral">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Initial Distance:</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">4</property>
            <property name="top_attach">3</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_step_size1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
          </object>
          <packing>
            <property name="left_attach">5</property>
            <property name="top_attach">0</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_nsteps1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
            <property name="text" translatable="yes">10</property>
          </object>
          <packing>
            <property name="left_attach">5</property>
            <property name="top_attach">1</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_FORCE_coord1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
            <property name="text" translatable="yes">4000</property>
          </object>
          <packing>
            <property name="left_attach">5</property>
            <property name="top_attach">2</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_dmin_coord1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
          </object>
          <packing>
            <property name="left_attach">5</property>
            <property name="top_attach">3</property>
          </packing>
        </child>
      </object>
      <packing>
        <property name="expand">True</property>
        <property name="fill">True</property>
        <property name="position">1</property>
      </packing>
    </child>
    <child>
      <object class="GtkButton" id="import_picking_selection_button">
        <property name="label" translatable="yes">Import from Picking Selection</property>
        <property name="visible">True</property>
        <property name="can_focus">True</property>
        <property name="receives_default">True</property>
      </object>
      <packing>
        <property name="expand">False</property>
        <property name="fill">True</property>
        <property name="position">2</property>
      </packing>
    </child>
    <child>
      <object class="GtkCheckButton" id="mass_restraints1">
        <property name="label" translatable="yes">Apply Mass Weighted Restraints</property>
        <property name="visible">True</property>
        <property name="can_focus">True</property>
        <property name="receives_default">False</property>
        <property name="draw_indicator">True</property>
      </object>
      <packing>
        <property name="expand">False</property>
        <property name="fill">True</property>
        <property name="position">3</property>
      </packing>
    </child>
  </object>
</interface>

'''

        elif mode == 0:
            xml = '''
<?xml version="1.0" encoding="UTF-8"?>
<!-- Generated with glade 3.22.2 -->
<interface>
  <requires lib="gtk+" version="3.20"/>
  <object class="GtkBox" id="reaction_coordinate_box">
    <property name="visible">True</property>
    <property name="can_focus">False</property>
    <property name="orientation">vertical</property>
    <property name="spacing">4</property>
    <child>
      <object class="GtkBox">
        <property name="visible">True</property>
        <property name="can_focus">False</property>
        <property name="spacing">6</property>
        <child>
          <object class="GtkLabel" id="label_RC_type">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Coordinate Type:</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="expand">False</property>
            <property name="fill">True</property>
            <property name="position">0</property>
          </packing>
        </child>
        <child>
          <object class="GtkComboBox" id="combobox_reaction_coord1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
          </object>
          <packing>
            <property name="expand">True</property>
            <property name="fill">True</property>
            <property name="position">1</property>
          </packing>
        </child>
      </object>
      <packing>
        <property name="expand">False</property>
        <property name="fill">True</property>
        <property name="position">0</property>
      </packing>
    </child>
    <child>
      <object class="GtkGrid">
        <property name="visible">True</property>
        <property name="can_focus">False</property>
        <property name="row_spacing">5</property>
        <property name="column_spacing">10</property>
        <child>
          <object class="GtkLabel" id="label_atom4_coord1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Atom 4(index):</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">0</property>
            <property name="top_attach">3</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_atom3_index_coord1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
          </object>
          <packing>
            <property name="left_attach">1</property>
            <property name="top_attach">2</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_atom2_coord1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Atom 2(index):</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">0</property>
            <property name="top_attach">1</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_atom1_coord1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Atom 1(index):</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">0</property>
            <property name="top_attach">0</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_atom3_coord1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Atom 3(index):</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">0</property>
            <property name="top_attach">2</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_name3_coord1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Name:</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">2</property>
            <property name="top_attach">2</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_name4_coord1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Name:</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">2</property>
            <property name="top_attach">3</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_name2_coord1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Name:</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">2</property>
            <property name="top_attach">1</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_name1_coord1">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Name:</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">2</property>
            <property name="top_attach">0</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_atom3_name_coord1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
          </object>
          <packing>
            <property name="left_attach">3</property>
            <property name="top_attach">2</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_atom2_index_coord1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
          </object>
          <packing>
            <property name="left_attach">1</property>
            <property name="top_attach">1</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_atom1_index_coord1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
          </object>
          <packing>
            <property name="left_attach">1</property>
            <property name="top_attach">0</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_atom4_index_coord1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
          </object>
          <packing>
            <property name="left_attach">1</property>
            <property name="top_attach">3</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_atom4_name_coord1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
          </object>
          <packing>
            <property name="left_attach">3</property>
            <property name="top_attach">3</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_atom2_name_coord1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
          </object>
          <packing>
            <property name="left_attach">3</property>
            <property name="top_attach">1</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_atom1_name_coord1">
            <property name="width_request">25</property>
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="halign">start</property>
            <property name="valign">start</property>
            <property name="width_chars">8</property>
          </object>
          <packing>
            <property name="left_attach">3</property>
            <property name="top_attach">0</property>
          </packing>
        </child>
      </object>
      <packing>
        <property name="expand">True</property>
        <property name="fill">True</property>
        <property name="position">1</property>
      </packing>
    </child>
    <child>
      <object class="GtkButton" id="import_picking_selection_button">
        <property name="label" translatable="yes">Import from Picking Selection</property>
        <property name="visible">True</property>
        <property name="can_focus">True</property>
        <property name="receives_default">True</property>
      </object>
      <packing>
        <property name="expand">False</property>
        <property name="fill">True</property>
        <property name="position">2</property>
      </packing>
    </child>
    <child>
      <object class="GtkCheckButton" id="mass_restraints1">
        <property name="label" translatable="yes">Apply Mass Weighted Restraints</property>
        <property name="visible">True</property>
        <property name="can_focus">True</property>
        <property name="receives_default">False</property>
        <property name="draw_indicator">True</property>
      </object>
      <packing>
        <property name="expand">False</property>
        <property name="fill">True</property>
        <property name="position">3</property>
      </packing>
    </child>
    <child>
      <object class="GtkAlignment" id="rc_aligment">
        <property name="visible">True</property>
        <property name="can_focus">False</property>
        <property name="top_padding">5</property>
        <property name="bottom_padding">5</property>
        <property name="left_padding">10</property>
        <property name="right_padding">10</property>
        <child>
          <object class="GtkSeparator">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
          </object>
        </child>
      </object>
      <packing>
        <property name="expand">False</property>
        <property name="fill">True</property>
        <property name="position">4</property>
      </packing>
    </child>
    <child>
      <object class="GtkGrid" id="rc_grid">
        <property name="visible">True</property>
        <property name="can_focus">False</property>
        <property name="halign">end</property>
        <property name="row_spacing">1</property>
        <property name="column_spacing">10</property>
        <child>
          <object class="GtkLabel" id="label_force">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Force Constante:</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">0</property>
            <property name="top_attach">2</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_initial_distance_angle_dihedral">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Initial Distance:</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">0</property>
            <property name="top_attach">3</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_FORCE_coord1">
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="text" translatable="yes">4000</property>
          </object>
          <packing>
            <property name="left_attach">1</property>
            <property name="top_attach">2</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_dmin_coord1">
            <property name="visible">True</property>
            <property name="can_focus">True</property>
          </object>
          <packing>
            <property name="left_attach">1</property>
            <property name="top_attach">3</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_num_of_steps">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Number of Steps:</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">0</property>
            <property name="top_attach">1</property>
          </packing>
        </child>
        <child>
          <object class="GtkLabel" id="label_step_size">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="label" translatable="yes">Step Size:</property>
            <property name="xalign">1</property>
          </object>
          <packing>
            <property name="left_attach">0</property>
            <property name="top_attach">0</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_nsteps1">
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="text" translatable="yes">10</property>
          </object>
          <packing>
            <property name="left_attach">1</property>
            <property name="top_attach">1</property>
          </packing>
        </child>
        <child>
          <object class="GtkEntry" id="entry_step_size1">
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="text" translatable="yes">0.1</property>
          </object>
          <packing>
            <property name="left_attach">1</property>
            <property name="top_attach">0</property>
          </packing>
        </child>
      </object>
      <packing>
        <property name="expand">False</property>
        <property name="fill">True</property>
        <property name="position">5</property>
      </packing>
    </child>
  </object>
</interface>

'''
        
        self.main = main
        self.vm_session = main.vm_session
        self.p_session  = main.p_session
        
        #scan mode ON
        self.scan = True
        
        self.builder = Gtk.Builder()
        self.builder.add_from_string(xml)
        #self.builder.connect_signals(self)
        #self.glWidget = glWidget
        box = self.builder.get_object('reaction_coordinate_box')
        self.pack_start(box, False, False, 0)
        
        self.import_picking_selection_button = self.builder.get_object('import_picking_selection_button')
        self.import_picking_selection_button.connect( 'clicked',  self.import_picking_selection_data )
        
        
        self.mass_restraints1 = self.builder.get_object('mass_restraints1')
        self.mass_restraints1.connect( 'toggled',  self.toggle_mass_restraint1 )



        #-----------------------------------------------------------------------------------
        self.method_store = Gtk.ListStore(str)
        methods = ["simple distance", "multiple distance", "multiple distance *4 atoms", 'dihedral']

        for method in methods:
            self.method_store.append([method])
        self.combobox_reaction_coord1 = self.builder.get_object('combobox_reaction_coord1')
        self.combobox_reaction_coord1.set_model(self.method_store)
        renderer_text = Gtk.CellRendererText()
        self.combobox_reaction_coord1.pack_start(renderer_text, True)
        self.combobox_reaction_coord1.add_attribute(renderer_text, "text", 0)
        self.combobox_reaction_coord1.connect('changed', self.change_cb_coordType1)
        #-----------------------------------------------------------------------------------


        self.entry_step_size = self.builder.get_object('entry_step_size1')
        self.entry_nsteps    = self.builder.get_object('entry_nsteps1')
        self.entry_dmin_coord= self.builder.get_object('entry_dmin_coord1')
        
        self.label_step_size      = self.builder.get_object('label_step_size')
        self.label_nsteps         = self.builder.get_object('label_nsteps')
        self.label_force_constant = self.builder.get_object('label_force_constant')
        self.label_dmin           = self.builder.get_object('label_dmin')
        self.label_initial_dist_angle_dihedral = self.builder.get_object('label_initial_distance_angle_dihedral')
        
        
    def toggle_mass_restraint1 (self, widget):
        """ Function doc """
        self.refresh_dmininum(coord1 =  True)
   
    
    def refresh_dmininum (self, coord1 =  False, coord2 = False):
        """ Function doc """
        _type = self.combobox_reaction_coord1.get_active()
        if _type == 0:
            index1 = int(self.builder.get_object('entry_atom1_index_coord1').get_text() )
            index2 = int(self.builder.get_object('entry_atom2_index_coord1').get_text() )

            dist1 = get_distance(self.vobject, index1, index2 )
            self.builder.get_object('entry_dmin_coord1').set_text(str(dist1))
        
        elif _type == 1:
            index1 = int(self.builder.get_object('entry_atom1_index_coord1').get_text() )
            index2 = int(self.builder.get_object('entry_atom2_index_coord1').get_text() )
            index3 = int(self.builder.get_object('entry_atom3_index_coord1').get_text() )
            
            dist1 = get_distance(self.vobject, index1, index2 )
            dist2 = get_distance(self.vobject, index2, index3 )
            
            if self.builder.get_object('mass_restraints1').get_active():
                self.sigma_pk1_pk3, self.sigma_pk3_pk1  = compute_sigma_a1_a3(self.vobject, index1, index3)
                #print('distance a1 - a2:', dist1 - dist2)
                DMINIMUM =  (self.sigma_pk1_pk3 * dist1) -(self.sigma_pk3_pk1 * dist2*-1)
                self.builder.get_object('entry_dmin_coord1').set_text(str(DMINIMUM))
            else:
                DMINIMUM =  dist1- dist2
                self.builder.get_object('entry_dmin_coord1').set_text(str(DMINIMUM))
        
        elif _type == 2:
            index1 = int(self.builder.get_object('entry_atom1_index_coord1').get_text() )
            index2 = int(self.builder.get_object('entry_atom2_index_coord1').get_text() )
            index3 = int(self.builder.get_object('entry_atom3_index_coord1').get_text() )
            index4 = int(self.builder.get_object('entry_atom4_index_coord1').get_text() )
            
            dist1 = get_distance(self.vobject, index1, index2 )
            dist2 = get_distance(self.vobject, index3, index4 )
            
            #if self.builder.get_object('mass_restraints1').get_active():
            #    self.sigma_pk1_pk3, self.sigma_pk3_pk1  = compute_sigma_a1_a3(self.vobject, index1, index3)
            #    #print('distance a1 - a2:', dist1 - dist2)
            #    DMINIMUM =  (self.sigma_pk1_pk3 * dist1) -(self.sigma_pk3_pk1 * dist2*-1)
            #    self.builder.get_object('entry_dmin_coord1').set_text(str(DMINIMUM))
            #else:
            DMINIMUM =  dist1- dist2
            self.builder.get_object('entry_dmin_coord1').set_text(str(DMINIMUM))
        
        elif _type == 3:
            index1 = int(self.builder.get_object('entry_atom1_index_coord1').get_text() )
            index2 = int(self.builder.get_object('entry_atom2_index_coord1').get_text() )
            index3 = int(self.builder.get_object('entry_atom3_index_coord1').get_text() )
            index4 = int(self.builder.get_object('entry_atom4_index_coord1').get_text() )
            
            dihedral = get_dihedral(self.vobject, index1, index2, index3, index4)
            #dist2 = get_distance(self.vobject, index2, index3 )
            self.builder.get_object('entry_dmin_coord1').set_text(str(dihedral))
            
            #if self.builder.get_object('mass_restraints1').get_active():
            #    self.sigma_pk1_pk3, self.sigma_pk3_pk1  = compute_sigma_a1_a3(self.vobject, index1, index3)
            #    #print('distance a1 - a2:', dist1 - dist2)
            #    DMINIMUM =  (self.sigma_pk1_pk3 * dist1) -(self.sigma_pk3_pk1 * dist2*-1)
            #    self.builder.get_object('entry_dmin_coord1').set_text(str(DMINIMUM))
            #else:
            #    DMINIMUM =  dist1- dist2
            #    self.builder.get_object('entry_dmin_coord1').set_text(str(DMINIMUM))
        
        
        
        else:
            pass    
              

    def import_picking_selection_data (self, widget):
        """  
                   R                    R
                    \                  /
                     A1--A2  . . . . A3
                    /                  \ 
                   R                    R
                     ^   ^            ^
                     |   |            |
                    pk1-pk2  . . . . pk3
                       d1       d2	

                q1 =  1 / (mpk1 + mpk3)  =  [ mpk1 * r (pk3_pk2)  -   mpk3 * r (pk1_pk2) ]
              
          mpk1 = mass of pk1 atom  
          mpk2 = mass of pk2 atom  
          mpk3 = mass of pk3 atom  
                
        """       
        atom1 = self.vm_session.picking_selections.picking_selections_list[0]
        atom2 = self.vm_session.picking_selections.picking_selections_list[1]
        atom3 = self.vm_session.picking_selections.picking_selections_list[2]
        atom4 = self.vm_session.picking_selections.picking_selections_list[3]
        
        if atom1:
            self.vobject = atom1.vm_object
        else:
            return None
                

        if atom1:
            self.builder.get_object('entry_atom1_index_coord1').set_text(str(atom1.index-1) )
            self.builder.get_object('entry_atom1_name_coord1' ).set_text(str(atom1.name) )
        else: print('use picking selection to chose the central atom')            
        #-------
        if atom2:
            self.builder.get_object('entry_atom2_index_coord1').set_text(str(atom2.index-1) )
            self.builder.get_object('entry_atom2_name_coord1' ).set_text(str(atom2.name) )
        else: print('use picking selection to chose the central atom')
        #-------
        if atom3:
            self.builder.get_object('entry_atom3_index_coord1').set_text(str(atom3.index-1) )
            self.builder.get_object('entry_atom3_name_coord1' ).set_text(str(atom3.name) )
            
            if atom3.symbol == atom1.symbol:
                self.builder.get_object('mass_restraints1').set_active(False)
            else:
                self.builder.get_object('mass_restraints1').set_active(True)
            
        else: print('use picking selection to chose the central atom')
         #-------
        if atom4:
            self.builder.get_object('entry_atom4_index_coord1').set_text(str(atom4.index-1) )
            self.builder.get_object('entry_atom4_name_coord1' ).set_text(str(atom4.name) )
        else: print('use picking selection to chose the central atom')
            
        self.refresh_dmininum( coord1 =  True)
            
        
    def change_cb_coordType1 (self, combo_box):
        """ Function doc """
        
        _type = self.combobox_reaction_coord1.get_active()        
        
        if _type == 0:
            self.builder.get_object('label_atom3_coord1').hide()
            self.builder.get_object('entry_atom3_index_coord1').hide()
            self.builder.get_object('label_name3_coord1').hide()
            self.builder.get_object('entry_atom3_name_coord1').hide()
            
            self.builder.get_object('label_atom4_coord1').hide()
            self.builder.get_object('entry_atom4_index_coord1').hide()
            self.builder.get_object('label_name4_coord1').hide()
            self.builder.get_object('entry_atom4_name_coord1').hide()
            self.builder.get_object('mass_restraints1').set_sensitive(False)
            self.label_initial_dist_angle_dihedral.set_text('Initial Distance:')

        if _type == 1:
            self.builder.get_object('label_atom3_coord1').show()
            self.builder.get_object('entry_atom3_index_coord1').show()
            self.builder.get_object('label_name3_coord1').show()
            self.builder.get_object('entry_atom3_name_coord1').show()
            
            self.builder.get_object('label_atom4_coord1').hide()
            self.builder.get_object('entry_atom4_index_coord1').hide()
            self.builder.get_object('label_name4_coord1').hide()
            self.builder.get_object('entry_atom4_name_coord1').hide()
            self.builder.get_object('mass_restraints1').set_sensitive(True)
            self.label_initial_dist_angle_dihedral.set_text('Initial Distance:')

        
        if _type == 2:
            self.builder.get_object('label_atom3_coord1').show()
            self.builder.get_object('entry_atom3_index_coord1').show()
            self.builder.get_object('label_name3_coord1').show()
            self.builder.get_object('entry_atom3_name_coord1').show()
            
            self.builder.get_object('label_atom4_coord1').show()
            self.builder.get_object('entry_atom4_index_coord1').show()
            self.builder.get_object('label_name4_coord1').show()
            self.builder.get_object('entry_atom4_name_coord1').show()
            self.builder.get_object('mass_restraints1').set_sensitive(False)
            #self.label_initial_dist_angle_dihedral.set_text('Initial Distance:')
            self.label_initial_dist_angle_dihedral.set_text('Initial Distance:')

        if _type == 3:
            self.builder.get_object('label_atom3_coord1').show()
            self.builder.get_object('entry_atom3_index_coord1').show()
            self.builder.get_object('label_name3_coord1').show()
            self.builder.get_object('entry_atom3_name_coord1').show()
            
            self.builder.get_object('label_atom4_coord1').show()
            self.builder.get_object('entry_atom4_index_coord1').show()
            self.builder.get_object('label_name4_coord1').show()
            self.builder.get_object('entry_atom4_name_coord1').show()
            self.builder.get_object('mass_restraints1').set_sensitive(False)
            #self.label_initial_dist_angle_dihedral.set_text('Initial Distance:')
            self.label_initial_dist_angle_dihedral.set_text('Initial Angle:')
        
        if _type == 4:
            self.label_initial_dist_angle_dihedral.set_text('Initial Angle:')
        
        
        
        if self.scan:
            try:
                self.refresh_dmininum ( coord1 = True)
            except:
                print(texto_d1)
                print(texto_d2d1)

    
    def set_rc_type (self, rc_type = 0):
        """ Function doc """
        self.combobox_reaction_coord1.set_active(rc_type)
        #print('aqui')


    def set_hide_scan_parameters (self):
        """ Function doc """
        self.builder.get_object('rc_grid').hide()
        self.builder.get_object('mass_restraints1').hide()
        self.builder.get_object('rc_aligment').hide()
        self.scan = False

    def get_rc_data (self):
        """ Function doc """
        parameters = { }
        
        _type = self.combobox_reaction_coord1.get_active()
        if _type == 0:
            parameters["rc_type"]     = "simple_distance"
            index1 = int( self.builder.get_object('entry_atom1_index_coord1').get_text() )
            index2 = int( self.builder.get_object('entry_atom2_index_coord1').get_text() )
            
            name1 = self.builder.get_object('entry_atom1_name_coord1').get_text()
            name2 = self.builder.get_object('entry_atom2_name_coord1').get_text()
            
            dmin   = float( self.builder.get_object('entry_dmin_coord1').get_text( ))
            parameters["ATOMS"]       = [ index1, index2 ] 
            parameters["ATOM_NAMES"] = [ name1 ,  name2 ] 
            parameters["dminimum"]  = dmin 
            parameters["sigma_pk1pk3"] = None
            parameters["sigma_pk3pk1"] = None
            
        elif _type == 1:
            parameters["rc_type"]     = "multiple_distance"
            index1 = int( self.builder.get_object('entry_atom1_index_coord1').get_text() )
            index2 = int( self.builder.get_object('entry_atom2_index_coord1').get_text() )
            index3 = int( self.builder.get_object('entry_atom3_index_coord1').get_text() )
            
            name1 = self.builder.get_object('entry_atom1_name_coord1').get_text() 
            name2 = self.builder.get_object('entry_atom2_name_coord1').get_text() 
            name3 = self.builder.get_object('entry_atom3_name_coord1').get_text() 
            
            dmin   = float( self.builder.get_object('entry_dmin_coord1').get_text( ))
            parameters["ATOMS"]       = [ index1, index2, index3 ] 
            parameters["ATOM_NAMES"] = [ name1,  name2,  name3 ] 
            parameters["dminimum"]  = dmin  
            
            if self.builder.get_object('mass_restraints1').get_active():
                parameters["MC"] = True
                parameters["sigma_pk1pk3"] = self.sigma_pk1_pk3 
                parameters["sigma_pk3pk1"] = self.sigma_pk3_pk1 
            else:
                parameters["sigma_pk1pk3"] =  1.0 
                parameters["sigma_pk3pk1"] = -1.0 


        elif _type == 2:
            index1 = int( self.builder.get_object('entry_atom1_index_coord1').get_text() )
            index2 = int( self.builder.get_object('entry_atom2_index_coord1').get_text() )
            index3 = int( self.builder.get_object('entry_atom3_index_coord1').get_text() )
            index4 = int( self.builder.get_object('entry_atom4_index_coord1').get_text() )
            dmin   = float( self.builder.get_object('entry_dmin_coord1').get_text() )
            
            name1 = self.builder.get_object('entry_atom1_name_coord1').get_text() 
            name2 = self.builder.get_object('entry_atom2_name_coord1').get_text() 
            name3 = self.builder.get_object('entry_atom3_name_coord1').get_text() 
            name4 = self.builder.get_object('entry_atom4_name_coord1').get_text() 
            
            
            parameters["ATOMS"]      = [ index1,index2,index3,index4] 
            parameters["ATOM_NAMES"] = [ name1,  name2,  name3, name4] 

            parameters["dminimum"]  = dmin 
            parameters["rc_type"]     = "multiple_distance*4atoms"
            parameters["sigma_pk1pk3"] =   1.0
            parameters["sigma_pk3pk1"] =  -1.0
        
        elif _type == 3:
            index1 = int( self.builder.get_object('entry_atom1_index_coord1').get_text() )
            index2 = int( self.builder.get_object('entry_atom2_index_coord1').get_text() )
            index3 = int( self.builder.get_object('entry_atom3_index_coord1').get_text() )
            index4 = int( self.builder.get_object('entry_atom4_index_coord1').get_text() )
            dmin   = float( self.builder.get_object('entry_dmin_coord1').get_text() )
            parameters["ATOMS"]     = [ index1,index2,index3,index4] 
            parameters["dminimum"]  = dmin 
            parameters["rc_type"]     = "dihedral"
            parameters["sigma_pk1pk3"] = None
            parameters["sigma_pk3pk1"] = None
        
        parameters["nsteps"]         = int(self.builder.get_object('entry_nsteps1').get_text() )
        parameters["force_constant"] = float( self.builder.get_object('entry_FORCE_coord1').get_text() )
        try:
            parameters["dincre"]         = float( self.builder.get_object('entry_step_size1').get_text() )
        except:
            parameters["dincre"]         = 0.00
        return (parameters)
        
        

class SequenceViewerBox(Gtk.Box):
    """ Class doc """
    
    def __init__ (self, main = None, mode = 0 ):
        """ Class initialiser """
        Gtk.Box.__init__(self)
        
        self.main = main
        self.seqview   = GtkSequenceViewer(main = main)
        #-----------------------------------------------------------------------
        self.seqview.add_events(Gdk.EventMask.POINTER_MOTION_MASK)
        self.seqview.add_events(Gdk.EventMask.BUTTON_PRESS_MASK)
        self.seqview.add_events(Gdk.EventMask.BUTTON_RELEASE_MASK)

        self.seqview.connect("motion-notify-event" , self.seqview.on_motion)
        self.seqview.connect("button_press_event"  , self.seqview.button_press)
        self.seqview.connect("button-release-event", self.seqview.button_release)
        #-----------------------------------------------------------------------

        #self.add()
        
        #--------------------------------------------------------------#
        self.box_horizontal1 = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 10)
        
        self.label1  = Gtk.Label()
        self.label1.set_text('Coordinates:')
        self.box_horizontal1.pack_start(self.label1, False, False, 0)
        self.coordinates_combobox = CoordinatesComboBox(self.main.vobject_liststore_dict[self.main.p_session.active_id])
        
        self.box_horizontal1.pack_start(self.coordinates_combobox, False, False, 0)
        #--------------------------------------------------------------#




        #------------------------------------------------------------------#
        #                  CHAIN combobox and Label
        #------------------------------------------------------------------#
        self.label2  = Gtk.Label()
        self.label2.set_text('Chain:')
        self.box_horizontal1.pack_start(self.label2, False, False, 0)

        self.liststore_chains = Gtk.ListStore(str)
        self.combobox_chains = Gtk.ComboBox.new_with_model(self.liststore_chains)
        self.combobox_chains.connect("changed", self.on_combobox_chains_changed)
        renderer_text = Gtk.CellRendererText() 
        self.combobox_chains.pack_start(renderer_text, True)
        self.combobox_chains.add_attribute(renderer_text, "text", 0)
        #vbox.pack_start(self.combobox_chains, False, False, True)
        self.box_horizontal1.pack_start(self.combobox_chains, False, False, 0)
        
        #Gtk.show_all()



    def on_combobox_chains_changed (self, widget):
        """ Function doc """
        
        tree_iter = self.combobox_chains.get_active_iter()
        if tree_iter is not None:
            model = self.combobox_chains.get_model()
            chain = model[tree_iter][0]
            #print("Selected: country=%s" % country)
        
        self.current_filter_chain = chain
        #print("%s Chain selected!" % self.current_filter_chain)
        # we update the filter, which updates in turn the view
        self.chain_filter.refilter()
        
        
        
        
        
        
        
        
        
