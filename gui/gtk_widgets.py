import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk
from gi.repository import Gdk
import os

import threading
import time

HOME        = os.environ.get('HOME')


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
            self.vm_session.viewing_selection_mode(sel_type = 'protein')
            #self.show_or_hide_entries (name = 'atom_names', show = False)
            #self.show_or_hide_entries (name = 'element', show = False)
        
        elif self.active == 4:
            self.vm_session.viewing_selection_mode(sel_type = 'C alpha')
            #self.show_or_hide_entries (name = 'atom_names', show = False)
            #self.show_or_hide_entries (name = 'element', show = False)
        
        elif self.active == 5:
            self.vm_session.viewing_selection_mode(sel_type = 'solvent')
            #self.show_or_hide_entries (name = 'atom_names', show = False)
            #self.show_or_hide_entries (name = 'element', show = False)
        
        elif self.active == 6:
            self.vm_session.viewing_selection_mode(sel_type = 'atom name')
            #self.show_or_hide_entries (name = 'atom_names', show = True)
            #self.show_or_hide_entries (name = 'element', show = False)
        
        elif self.active == 7:
            self.vm_session.viewing_selection_mode(sel_type = 'element')
            #self.show_or_hide_entries (name = 'atom_names', show = False)
            #self.show_or_hide_entries (name = 'element', show = True)
                
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



class EasyHybridGoToAtomWindow(Gtk.Window):
    def OpenWindow (self):
        """ Function doc """
        if self.visible  ==  False:
            
            #self.vm_session.Vismol_Objects_ListStore
            
            #------------------------------------------------------------------#
            #                  SYSTEM combobox and Label
            #------------------------------------------------------------------#
            self.box_vertical    = Gtk.Box(orientation = Gtk.Orientation.VERTICAL,   spacing = 10)
            self.box_horizontal1 = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 10)
             

            
            self.label1  = Gtk.Label()
            self.label1.set_text('System:')
            self.box_horizontal1.pack_start(self.label1, False, False, 0)

            
            self.combobox_systems = Gtk.ComboBox.new_with_model(self.system_liststore)
            self.combobox_systems.connect("changed", self.on_combobox_systems_changed)
            renderer_text = Gtk.CellRendererText()
            self.combobox_systems.pack_start(renderer_text, True)
            self.combobox_systems.add_attribute(renderer_text, "text", 0)
            
            self.box_horizontal1.pack_start(self.combobox_systems, False, False, 0)
            #------------------------------------------------------------------#
            
            
            
            
            
            #------------------------------------------------------------------#
            #                  COORDINATES combobox and Label
            #------------------------------------------------------------------#
            self.label1  = Gtk.Label()
            self.label1.set_text('Coordinates:')
            self.box_horizontal1.pack_start(self.label1, False, False, 0)
            self.coordinates_combobox = Gtk.ComboBox.new_with_model(self.coordinates_liststore)
            #self.coordinates_combobox.connect("changed", self.on_self.coordinates_combobox_changed)
            renderer_text = Gtk.CellRendererText()
            self.coordinates_combobox.pack_start(renderer_text, True)
            self.coordinates_combobox.add_attribute(renderer_text, "text", 0)
            self.box_horizontal1.pack_start(self.coordinates_combobox, False, False, 0)
            #------------------------------------------------------------------#
            
            


            
            #------------------------------------------------------------------#
            #                  CHAIN combobox and Label
            #------------------------------------------------------------------#
            self.box_horizontal2 = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
            
            
            self.label2  = Gtk.Label()
            self.label2.set_text('Chain:')
            self.box_horizontal2.pack_start(self.label2, False, False, 0)

            self.liststore_chains = Gtk.ListStore(str)
            
            self.combobox_chains = Gtk.ComboBox.new_with_model(self.liststore_chains)
            self.combobox_chains.connect("changed", self.on_combobox_chains_changed)
            renderer_text = Gtk.CellRendererText()
            self.combobox_chains.pack_start(renderer_text, True)
            self.combobox_chains.add_attribute(renderer_text, "text", 0)
            #vbox.pack_start(self.combobox_chains, False, False, True)
            self.box_horizontal2.pack_start(self.combobox_chains, False, False, 0)
            
            
            
            
            #------------------------------------------------------------------#
            #                  RESIDUES combobox and Label
            #------------------------------------------------------------------#
            
            self.label3  = Gtk.Label()
            self.label3.set_text('Residue type:')
            self.box_horizontal2.pack_start(self.label3, False, False, 0)

            self.liststore_residues = Gtk.ListStore(str)
            
            self.combobox_residues = Gtk.ComboBox.new_with_model(self.liststore_residues)
            self.combobox_residues.connect("changed", self.on_combobox_residues_changed)
            renderer_text = Gtk.CellRendererText()
            self.combobox_residues.pack_start(renderer_text, True)
            self.combobox_residues.add_attribute(renderer_text, "text", 0)
            #vbox.pack_start(self.combobox_chains, False, False, True)
            self.box_horizontal2.pack_start(self.combobox_residues, False, False, 0)
            
            
            #------------------------------------------------------------------#
            
            
            
            
            
            self.treeviewbox_horizontal = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
            
            #------------------------------------------------------------------------------------------
            #self.treeview = Gtk.TreeView(model =self.residue_liststore)
            
            
            #-----------------------------------------------------------------------------------------
            #                                      Chain filter
            #-----------------------------------------------------------------------------------------
            self.current_filter_chain = None
            # Creating the filter, feeding it with the liststore model
            self.chain_filter = self.residue_liststore.filter_new()
            # setting the filter function, note that we're not using the
            self.chain_filter.set_visible_func(self.chain_filter_func)
            #-----------------------------------------------------------------------------------------
            
            
            self.treeview = Gtk.TreeView(model = self.chain_filter)
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.MULTIPLE)
            
            self.treeview.connect("button-release-event", self.on_treeview_Objects_button_release_event)
            self.treeview.connect("row-activated", self.on_treeview_row_activated_event)
            
            for i, column_title in enumerate(
                ['', "index", "Residue",  "Chain", 'size']
            ):
                if i == 0:
                    cell = Gtk.CellRendererToggle()
                    cell.set_property('activatable', True)
                    cell.connect('toggled', self.on_chk_renderer_toggled, self.residue_liststore)
                    column = Gtk.TreeViewColumn(column_title, cell )
                    column.add_attribute(cell, 'active', 0)
                    self.treeview.append_column(column)
                    #print ('aqui')
                else:
                    renderer = Gtk.CellRendererText()
                    #renderer.connect('toggled', self.on_chk_renderer_toggled, self.residue_liststore)
                    column = Gtk.TreeViewColumn(column_title, renderer, text=i)
                    self.treeview.append_column(column)

            
            self.current_filter_chain = None
            

            self.scrollable_treelist = Gtk.ScrolledWindow()
            self.scrollable_treelist.set_vexpand(True)
            self.scrollable_treelist.add(self.treeview)
            #------------------------------------------------------------------------------------------
            


            
            
            #------------------------------------------------------------------------------------------
            self.treeview_atom = Gtk.TreeView(model =self.atom_liststore)
            self.treeview_atom.connect("button-release-event", self.on_treeview_atom_button_release_event)
            self.treeview_atom.connect("row-activated", self.on_treeview_atom_row_activated_event)

            for i, column_title in enumerate(
                ['', "index", "name", "MM atom", 'MM charge']
            ):
                if i == 0:
                    cell = Gtk.CellRendererToggle()
                    cell.set_property('activatable', True)
                    cell.connect('toggled', self.on_chk_renderer_toggled, self.atom_liststore)
                    column = Gtk.TreeViewColumn(column_title, cell )
                    column.add_attribute(cell, 'active', 0)
                    self.treeview_atom.append_column(column)
                    #print ('aqui')
                else:
                    renderer = Gtk.CellRendererText()
                    column = Gtk.TreeViewColumn(column_title, renderer, text=i)
                    self.treeview_atom.append_column(column)

            self.scrollable_treelist2 = Gtk.ScrolledWindow()
            self.scrollable_treelist2.set_vexpand(True)
            self.scrollable_treelist2.add(self.treeview_atom)
            #------------------------------------------------------------------------------------------
            
            
            
            
            
            
            
            
            
            
            
            self.box_vertical.pack_start(self.box_horizontal1, False, True, 0)
            self.box_vertical.pack_start(self.box_horizontal2, False, True, 0)
            self.treeviewbox_horizontal.pack_start(self.scrollable_treelist, True, True, 0)
            self.treeviewbox_horizontal.pack_start(self.scrollable_treelist2, True, True, 0)
            
            self.box_vertical.pack_start(self.treeviewbox_horizontal, False, True, 0)
            
            self.refresh_system_liststore()
            self.update_window (system_names = True, coordinates = True)
            
            self.combobox_systems.set_active(0)
            self.window =  Gtk.Window()
            self.window.set_border_width(10)
            self.window.set_default_size(600, 600)  
            self.window.add(self.box_vertical)
            self.window.connect("destroy", self.CloseWindow)
            self.window.set_title('Go to Atom Window') 
            self.window.show_all() 
                                                          
            #                                                                
            #self.builder.connect_signals(self)                                   
            
            self.visible  =  True
            #self.PutBackUpWindowData()
            #gtk.main()
            #----------------------------------------------------------------

    def CloseWindow (self, button):
        """ Function doc """
        #self.BackUpWindowData()
        self.window.destroy()
        self.visible    =  False
        #print('self.visible',self.visible)
    
    def __init__(self, main = None, system_liststore = None):
        """ Class initialiser """
        
        self.main_session  = main
        self.vm_session    = self.main_session.vm_session
        self.p_session     = self.main_session.p_session
        self.system_liststore      = system_liststore
        self.coordinates_liststore = Gtk.ListStore(str, int)
        
        
        self.residue_liststore = Gtk.ListStore(bool, int, str, str, int)
        self.atom_liststore    = Gtk.ListStore(bool, int, str, str, float, int, )
        self.residue_filter    = False
        self.visible           = False

    
    def update_window (self, system_names = False, coordinates = True,  selections = False ):
        """ Function doc """

        if self.visible:
            
            _id = self.combobox_systems.get_active()
            if _id == -1:
                '''_id = -1 means no item inside the combobox'''
                return None
            else:    
                _, system_id = self.system_liststore[_id]
            
            if system_names:
                self.refresh_system_liststore ()
                self.combobox_systems.set_active(_id)
            
            if coordinates:
                self.refresh_coordinates_liststore ()
                
            
            #if selections:
            #    _, system_id = self.system_liststore[_id]
            #    self.refresh_selection_liststore(system_id)
        else:
            if system_names:
                self.refresh_system_liststore ()
            if coordinates:
                self.refresh_coordinates_liststore ()
            
    
    def refresh_coordinates_liststore(self, system_id = None):
        """ Function doc """
        cb_id =  self.coordinates_combobox.get_active()
        if system_id:
            pass
        else:
            _id = self.combobox_systems.get_active()
            if _id == -1:
                return False
            else:
                #print('_id', _id)
                _, system_id = self.system_liststore[_id]
        
        self.coordinates_liststore.clear()
        n = 0
        for key , vobject in self.vm_session.vobjects_dic.items():
            if vobject.easyhybrid_system_id == system_id:
                self.coordinates_liststore.append([vobject.name, key])
                n += 1
        
        self.coordinates_combobox.set_active(n-1)
        
    def refresh_system_liststore (self):
        """ Function doc """
        self.main_session.refresh_system_liststore()
        #self.system_liststore     .clear()
        ##self.selection_liststore  .clear()
        #'''--------------------------------------------------------------------------------------------'''
        ##self.combobox_systems =self.builder.get_object('systems_combobox')
        #for key, system  in self.p_session.systems.items():
        #    try:
        #        self.system_liststore.append([system['name'], key])
        #    except:
        #        print(system)

    def on_combobox_residues_changed (self, widget):
        """ Function doc """
        tree_iter = widget.get_active_iter()
        if tree_iter is not None:
            model = widget.get_model()
            residue = model[tree_iter][0]
            #print("Selected: country=%s" % country)
        
            self.current_filter_residue = residue
            
            #print("%s Chain selected!" % self.current_filter_residue)
            
            # we update the filter, which updates in turn the view
            if self.residue_filter:
                self.residue_filter.refilter()
        
        
    def on_combobox_chains_changed (self, widget):
        """ Function doc """
        ##---------------------------------------------------------------
        #self.current_filter_chain = None
        ## Creating the filter, feeding it with the liststore model
        #self.chain_filter = self.residue_liststore.filter_new()
        ## setting the filter function, note that we're not using the
        #self.chain_filter.set_visible_func(self.chain_filter_func)
        ##---------------------------------------------------------------
        
        tree_iter = self.combobox_chains.get_active_iter()
        if tree_iter is not None:
            model = self.combobox_chains.get_model()
            chain = model[tree_iter][0]
            #print("Selected: country=%s" % country)
        
        self.current_filter_chain = chain
        #print("%s Chain selected!" % self.current_filter_chain)
        # we update the filter, which updates in turn the view
        self.chain_filter.refilter()
    
    
    def on_combobox_systems_changed (self, widget):
        """ Function doc """
        #print(widget)
        #print(widget.get_active())
        
        #print(widget.get_active_id())
        #print(widget.get_active_iter())
        
        #self.vm_session.vobjects_dic.items()
        
        #self.vm_session.vobjects_dic.items()
        cb_id = widget.get_active()
        if cb_id == -1:
            return None
        else:
            
            self.update_window (coordinates = True)
            
            cb_id =  self.coordinates_combobox.get_active()
            _, key = self.coordinates_liststore[cb_id]
            
            self.VObj = self.vm_session.vobjects_dic[key]
            
            
            self.liststore_chains = Gtk.ListStore(str)
            self.liststore_chains.append(['all'])
            chains = self.VObj.chains.keys()
            
            #self.chain_liststore = Gtk.ListStore(str)
            
            for chain in chains:
                self.liststore_chains.append([chain])
            self.combobox_chains.set_model(self.liststore_chains)
            self.combobox_chains.set_active(0)
            
            
            self.residue_liststore = Gtk.ListStore(bool, int, str, str, int)
            for chain in self.VObj.chains:
                for res in self.VObj.chains[chain].residues:
                    #print(res.resi, res.resn, chain,  len(res.atoms) ) 
                    
                    self.residue_liststore.append(list([True, res.resi, res.resn, chain,  len(res.atoms)]))
            #-----------------------------------------------------------------------------------------
            #                                      Chain filter
            #-----------------------------------------------------------------------------------------
            self.current_filter_chain = None
            # Creating the filter, feeding it with the liststore model
            self.chain_filter = self.residue_liststore.filter_new()
            # setting the filter function, note that we're not using the
            self.chain_filter.set_visible_func(self.chain_filter_func)
            #-----------------------------------------------------------------------------------------
            
            
            
            
            
            
            #-----------------------------------------------------------------------------------------
            #                                      Residue combobox
            #-----------------------------------------------------------------------------------------
            self.liststore_residues = Gtk.ListStore(str)
            self.liststore_residues.append(['all'])
            
            resn_labels = {}
            
            for residue in self.VObj.residues:
                resn_labels[residue.resn] = True
            
            for resn in resn_labels.keys():
                #print (resn)
                self.liststore_residues.append([resn])
            
            self.combobox_residues.set_model(self.liststore_residues)
            self.combobox_residues.set_active(0)
            
            #-----------------------------------------------------------------------------------------
            #                                      Residue filter
            #-----------------------------------------------------------------------------------------
            self.current_filter_residue = None
            # Creating the filter, feeding it with the liststore model
            self.residue_filter = self.chain_filter.filter_new()
            # setting the filter function, note that we're not using the
            self.residue_filter.set_visible_func(self.residue_filter_func)
            #-----------------------------------------------------------------------------------------        
            
            
            #self.treeview.set_model(self.residue_liststore)
            #self.treeview.set_model(self.chain_filter)
            self.treeview.set_model(self.residue_filter)
            
        
    def on_treeview_atom_button_release_event(self, tree, event):
        if event.button == 2:
            selection     = tree.get_selection()
            model         = tree.get_model()
            (model, iter) = selection.get_selected()
            
            
            if iter != None:
                self.selectedID  = int(model.get_value(iter, 1))-1  # @+
                atom = self.VObj.atoms[self.selectedID]
                self.vm_session.glwidget.vm_widget.center_on_atom(atom)
       
    def on_treeview_atom_row_activated_event (self, tree, rowline , column):
        """ Function doc """
        #print (A,B,C)
        selection     = tree.get_selection()
        model         = tree.get_model()
        
        #print(model)
        #print(rowline, list(model[rowline]))
        
        data  = list(model[rowline])
        #print(data)
        #self.selectedID  = int(data[1])  # @+
        #self.selectedObj = str(data[2])
        pickedID = data[-1]
        
        atom_picked = self.vm_session.atom_dic_id[pickedID]
        
        
        #res = self.VObj.chains[self.selectedChn].residues_by_index[self.selectedID]
        #
        ##print('Selecting by doble click')
        ##self.vm_session.selections[self.vm_session.current_selection].selecting_by_residue (res.atoms[0])
        ##self.vm_session.selections[self.vm_session.current_selection].selection_function_viewing (res.atoms[0])
       
        self.vm_session._selection_function (atom_picked, _type = 'atom')
        self.vm_session.glwidget.queue_draw()
        #
        #self.atom_liststore.clear()
        #for atom in res.atoms:
        #    self.atom_liststore.append(list([True, int(atom.index), atom.name, atom.symbol, atom.charge]))
    
    
    def on_treeview_row_activated_event(self, tree, rowline , column ):
        #print (A,B,C)
        selection     = tree.get_selection()
        model         = tree.get_model()
        
        #print(model)
        #print(rowline, list(model[rowline]))
        
        data  = list(model[rowline])
        self.selectedID  = int(data[1])  # @+
        self.selectedObj = str(data[2])
        self.selectedChn = str(data[3])
        
        cb_id =  self.coordinates_combobox.get_active()
        _, key = self.coordinates_liststore[cb_id]
        self.VObj = self.vm_session.vobjects_dic[key]
        
        res = self.VObj.chains[self.selectedChn].residues_by_index[self.selectedID]
        
        
        '''centering and selecting'''
        frame = self.vm_session.get_frame ()
        res.get_center_of_mass(frame = frame)
        self.vm_session.glwidget.vm_widget.center_on_coordinates(res.vobject, res.mass_center)
        
        self.vm_session._selection_function (res.atoms[0], _type = 'residue')
        self.vm_session.glwidget.queue_draw()
        
        self.atom_liststore.clear()
        
        for atom in res.atoms:
            self.atom_liststore.append(list([True, int(atom.index), atom.name, atom.symbol, atom.charge, atom.atom_id ]))

        #self.treeview_atom.set_model(self.atom_liststore)
  
    
    def on_treeview_Objects_button_release_event(self, tree, event):
        #print ( tree, event)
        
        if event.button == 3:
            print ("button 3")
            #selection     = tree.get_selection()
            #model         = tree.get_model()
            #(model, iter) = selection.get_selected()
            #
            #
            #
            #if iter != None:
            #    self.selectedID  = int(model.get_value(iter, 1))  # @+
            #    self.selectedObj = str(model.get_value(iter, 2))
            #    print(self.selectedID, self.selectedObj, self.VObj.residues[self.selectedID])
            #    #self.builder.get_object('TreeViewObjLabel').set_label('- ' +self.selectedObj+' -' )
            #
            #    #widget = self.builder.get_object('treeview_menu')
            #    #widget.popup(None, None, None, None, event.button, event.time)
            #    #print ('button == 3')


        if event.button == 2:
            #print ('button == 2')
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.SINGLE)

            
            selection     = tree.get_selection()
            model         = tree.get_model()
            (model, iter) = selection.get_selected()

            if iter != None:
                self.selectedID  = int(model.get_value(iter, 1))  # @+
                self.selectedObj = str(model.get_value(iter, 2))
                self.selectedChn = str(model.get_value(iter, 3))
                res = self.VObj.chains[self.selectedChn].residues_by_index[self.selectedID]
                frame = self.vm_session.get_frame ()
                res.get_center_of_mass(frame = frame)
                
                self.vm_session.glwidget.vm_widget.center_on_coordinates(res.vobject, res.mass_center)
        
                self.atom_liststore.clear()
                for atom in res.atoms:
                     self.atom_liststore.append(list([True, int(atom.index), atom.name, atom.symbol, atom.charge, atom.atom_id ]))
            
            
            
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.MULTIPLE)

        
        if event.button == 1:
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.SINGLE)

            selection     = tree.get_selection()
            model         = tree.get_model()
            (model, iter) = selection.get_selected()
            
            
            if iter != None:
                self.selectedID  = int(model.get_value(iter, 1))  # @+
                self.selectedObj = str(model.get_value(iter, 2))
                self.selectedChn = str(model.get_value(iter, 3))
                res = self.VObj.chains[self.selectedChn].residues_by_index[self.selectedID]
                
                
                self.atom_liststore.clear()
                #self.atom_liststore = Gtk.ListStore(bool, int, str, str, float)
                try:
                    #charges = list(self.p_session.systems[self.VObj.easyhybrid_system_id]['system'].AtomicCharges())
                    charges         = list(self.p_session.systems[self.VObj.easyhybrid_system_id]['system'].mmState.charges)
                    atomTypes       = self.p_session.systems[self.VObj.easyhybrid_system_id]['system'].mmState.atomTypes
                    atomTypeIndices = list(self.p_session.systems[self.VObj.easyhybrid_system_id]['system'].mmState.atomTypeIndices)

                    for atom in res.atoms:
                         #self.atom_liststore.append(list([True, int(atom.index), atom.name, atom.symbol,  charges[atom.index-1] , int(atom.atom_id)]))
                         self.atom_liststore.append(list([True, int(atom.index), atom.name, atomTypes[atomTypeIndices[atom.index-1]] ,  charges[atom.index-1] , int(atom.atom_id)]))
                except:
                    charges = [0]*len(self.VObj.atoms)
                    for atom in res.atoms:
                         self.atom_liststore.append(list([True, int(atom.index), atom.name, 'UNK',  charges[atom.index-1] , int(atom.atom_id)]))
                
                
                
                #self.treeview_atom.set_model(self.atom_liststore)
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.MULTIPLE)

        
      
    
    def on_chk_renderer_toggled(self, cell, path, model):
        print(model[path][0])

   
   
    def residue_filter_func(self, model, iter, data):
        """Tests if the language in the row is the one in the filter"""
        if (
            self.current_filter_residue is None
            or self.current_filter_residue == "all"
        ):
            return True
        else:
            return model[iter][2] == self.current_filter_residue
   
            
    def chain_filter_func(self, model, iter, data):
        """Tests if the language in the row is the one in the filter"""
        if (
            self.current_filter_chain is None
            or self.current_filter_chain == "all"
        ):
            return True
        else:
            return model[iter][3] == self.current_filter_chain

    #def on_selection_button_clicked(self, widget):
    #    """Called on any of the button clicks"""
    #    # we set the current language filter to the button's label
    #    self.current_filter_language = widget.get_label()
    #    print("%s language selected!" % self.current_filter_language)
    #    # we update the filter, which updates in turn the view
    #    self.language_filter.refilter()
    #
    #def on_button1_clicked(self, widget):
    #    print("Hello")
    #
    #def on_button2_clicked(self, widget):
    #    print("Goodbye")
    def update (self):
        """ Function doc """
        print('VismolGoToAtomWindow2 update')
        pass



class SaveTrajectoryBox:
    """ Class doc """
    
    def __init__ (self, parent, home):
        """ Class initialiser """
        
        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(home,'gui/windows/trajectory_box.glade'))
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
        img.set_from_file(os.path.join(self.main.home,'gui/icons/icon_open.png'))
        self.label = Gtk.Label ('...')
        self.btn.connect('clicked', self.open_filechooser)

        grid.attach (img, 0, 0, 1, 1)
        grid.attach (self.label, 1, 0, 1, 1)
        grid.show_all ()    

        self.btn.add (grid)        
        #return self.btn 
        self.folder = None
    #====================================================================================
    def set_folder (self, folder = '/home'):
        """ Function doc """
        #print('set_folder', folder)
        self.folder = folder
        name = os.path.basename(folder )
        #print( name)
        self.label.set_text(name)
        return folder
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
            self.set_folder(folder = dialog.get_filename())
            #print(os.path.dirname( dialog.get_filename() ))

        
        elif response == Gtk.ResponseType.CANCEL:
            print("Cancel clicked")

        dialog.destroy()
       


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


class VismolTrajectoryFrame(Gtk.Frame):
    """ Class doc """
    
    def __init__ (self, vm_session = None):
        """ Class initialiser """
        self.vm_session = vm_session 
        
        self.frame      =Gtk.Frame()
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
        
        self.functions = [self.reverse,self.stop, self.play, self.forward ]
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
        
        #self.entry = Gtk.Entry()
        #self.entry.set_text(str(25))
        
        self.fps_adjustment = Gtk.Adjustment(value          = 24 , 
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
        
        #self.vbox2.pack_start(self.label2, True, True, 0)
        #self.vbox2.pack_start(self.combobox_vobjects, True, True, 0)
        self.vbox2.pack_start(self.label, True, True, 0)
        self.vbox2.pack_start(self.entry, True, True, 0)
        self.box.pack_start(self.vbox2, True, True, 0)
        #----------------------------------------------------------------------------

    def play (self, button):
        # Create and start a new thread when the button is clicked
        if self.running:
            pass
        
        else:
            thread = threading.Thread(target=self.long_task)
            thread.start()
            
            self.running = True
        
        
    def long_task(self):
        # This flag is used to stop the thread
        self.stop_thread = False
        i = 0
        while self.stop_thread == False:
            if self.stop_thread:
                return
            value = self.forward(None) 

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

    def forward (self, button):
        """ Function doc """
        value =  int(self.scale.get_value())
        value = value+1
        self.scale.set_value(int(value))
        self.vm_session.set_frame(int(value))
        print(value)
        return value
        
        
    def reverse (self, button):
        """ Function doc """
        value = int(self.scale.get_value())
        
        if value == 0:
            pass
        else:
           value = value-1
        
        self.vm_session.set_frame(int(value))
        self.scale.set_value(value)
        print(value)
    
    def get_box (self):
        """ Function doc """
        #self.add(self.box)
        return self.box
        #return self.frame
        
    def on_combobox_vobjects_changed (self, widget):
        """ Function doc """
        print('\n\n',widget)
        print('\n\n',widget.get_active())
        
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
        print('VismolTrajectoryFrame update')
        #for index , vobject in self.vm_session.vobjects_dic.items():
        #last_obj = len(self.vm_session.vobjects) -1
        last_obj = len(self.vm_session.vobjects_dic.items()) -1
        self.combobox_vobjects.set_active(last_obj)
    def change_range (self, upper = 100):
        """ Function doc """
        print('upper =', upper)
        self.adjustment     = Gtk.Adjustment(value         = self.value,
                                             lower         = 0,
                                             upper         = upper,
                                             step_increment= 1,
                                             page_increment= 1,
                                             page_size     = 1)
        self.scale.set_adjustment ( self.adjustment)
        self.scale.set_digits(0)
        self.upper = upper
