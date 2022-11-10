import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
from gi.repository import Gdk

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
        data = ['atom'   , 
                'residue',
                'chain'  , 
                'protein', 
                'C alpha',
                'solvent',
                'atom name',
                'element',
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

    def show_or_hide_entries (self, name = 'atom_names', show = True):
        """ Function doc """
        if name == 'atom_names':
            if show:
                if self.entry_atom_names:
                    self.entry_atom_names.show()
                else:
                    self.entry_atom_names = Gtk.Entry()
                    self.entry_atom_names.set_width_chars(8)
                    self.entry_atom_names.connect('changed', self.on_entry_data_change)
                    self.box.pack_start(self.entry_atom_names     , False, False, 0)
                    self.entry_atom_names.show()
            else:
                if self.entry_atom_names:
                    self.entry_atom_names.hide()
                else:
                    pass
                    #self.entry_atom_names = Gtk.Entry()
                    #self.entry_atom_names.set_width_chars(8)
                    #self.entry_atom_names.connect('changed', self.on_entry_data_change)
                    #self.box.pack_start(self.entry_atom_names     , False, False, 0)
                    #self.entry_atom_names.hide()
        
        if name == 'element':
            if show:
                if self.entry_elements:
                    self.entry_elements.show()
                else:
                    self.entry_elements = Gtk.Entry()
                    self.entry_elements.set_width_chars(8)
                    self.entry_elements.connect('changed', self.on_entry_data_change)
                    self.box.pack_start(self.entry_elements     , False, False, 0)
                    self.entry_elements.show()
            else:
                if self.entry_elements:
                    self.entry_elements.hide()
                else:
                    pass
                    #self.entry_elements = Gtk.Entry()
                    #self.entry_elements.connect('changed', self.on_entry_data_change)
                    #self.entry_elements.set_width_chars(8)
                    #self.box.pack_start(self.entry_elements     , False, False, 0)
                    #self.entry_elements.hide()
                
                
    def on_combobox_selection_type (self, combobox):
        """ Function doc """
        self.active = combobox.get_active()

        if self.active == 0:
            self.vm_session.viewing_selection_mode(sel_type = 'atom')
            self.show_or_hide_entries (name = 'atom_names', show = False)
            self.show_or_hide_entries (name = 'element', show = False)

        elif self.active == 1:
            self.vm_session.viewing_selection_mode(sel_type = 'residue')
            self.show_or_hide_entries (name = 'atom_names', show = False)
            self.show_or_hide_entries (name = 'element', show = False)
        
        elif self.active == 2:
            self.vm_session.viewing_selection_mode(sel_type = 'chain')
            self.show_or_hide_entries (name = 'atom_names', show = False)
            self.show_or_hide_entries (name = 'element', show = False)
        
        elif self.active == 3:
            self.vm_session.viewing_selection_mode(sel_type = 'protein')
            self.show_or_hide_entries (name = 'atom_names', show = False)
            self.show_or_hide_entries (name = 'element', show = False)
        
        elif self.active == 4:
            self.vm_session.viewing_selection_mode(sel_type = 'C alpha')
            self.show_or_hide_entries (name = 'atom_names', show = False)
            self.show_or_hide_entries (name = 'element', show = False)
        
        elif self.active == 5:
            self.vm_session.viewing_selection_mode(sel_type = 'solvent')
            self.show_or_hide_entries (name = 'atom_names', show = False)
            self.show_or_hide_entries (name = 'element', show = False)
        
        elif self.active == 6:
            self.vm_session.viewing_selection_mode(sel_type = 'atom name')
            self.show_or_hide_entries (name = 'atom_names', show = True)
            self.show_or_hide_entries (name = 'element', show = False)
        
        elif self.active == 7:
            self.vm_session.viewing_selection_mode(sel_type = 'element')
            self.show_or_hide_entries (name = 'atom_names', show = False)
            self.show_or_hide_entries (name = 'element', show = True)
                
        else:pass
        
    def change_sel_type_in_combobox (self, sel_type):
        """ Function doc """
        
        if sel_type == 'atom':
            self.combobox_selection_type.set_active(0)
        
        elif sel_type == 'residue':
            self.combobox_selection_type.set_active(1)
        
        elif sel_type == 'chain':
            self.combobox_selection_type.set_active(2)
        
        elif sel_type == 'protein':
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
            self.vm_session._picking_selection_mode = True
            button.set_label('Picking')
            #print(self.combobox_selection_type.get_active())
            self.vm_session._selection_function (None)
            self.vm_session.glwidget.vm_widget.queue_draw()
            
            self.combobox_selection_type.set_sensitive(False)
            self.label_selecting_by.set_sensitive(False)
            
            circle = Gdk.Cursor(Gdk.CursorType.CROSSHAIR)
            button.get_window().set_cursor(circle)
            
        else:
            state = "off"
            self.vm_session._picking_selection_mode = False
            button.set_label('Viewing')
            self.vm_session.glwidget.vm_widget.queue_draw()
            
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





class FileChooser:
    """ Class doc """
    
    def __init__ (self, main_window = None , easyhybrid = False):
        """ Class initialiser """
        self.main_window = main_window
        self.easyhybrid = easyhybrid 
    
    def open (self, select_multiple = False, filters = None):

        """ Function doc """
        #main = gtkmain
        main = self.main_window
        filename = None
        
        chooser = Gtk.FileChooserDialog("Open File...", main,0,
                                       (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                                        Gtk.STOCK_OPEN, Gtk.ResponseType.OK))
        #GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER
        print (filters)

        if select_multiple:
            chooser.set_select_multiple(True)
            response = chooser.run()
            if response == Gtk.ResponseType.OK:
                filenames = chooser.get_filenames()
            chooser.destroy()
            return filenames
            
            
            
        
        
        
        else:
            if filters:
                print('\n\nfilters')
                for _filter in filters:
                    chooser.add_filter(_filter)

            else:
                print('else')

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
            chooser.destroy()
            return filename
