import gi
import os
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Pango

from gui.widgets.custom_widgets import FolderChooserButton
#from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.widgets.custom_widgets import ReactionCoordinateBox

class EnergyRefinementWindow():

    def __init__(self, main = None ):
        """ Class initialiser """
        self.main       = main
        self.home       = main.home
        self.p_session  = main.p_session
        self.vm_session = main.vm_session
        self.Visible    =  False        
        
        self.vobject_liststore   = Gtk.ListStore(str, int)
        self.data_liststore      = Gtk.ListStore(str, int)
        
        
        self.input_types_liststore = Gtk.ListStore(str)
        self.input_types = {
                        0:'Vobject'    , #  single file
                        1:'pklfolder - pDynamo Trajectory'  , #  trajectory
                        2:'pklfolder2D - pDynamo 2D Trajectory', #  2d trajectory  
                      # 3:'pdbfile'    ,
                      # 4:'pdbfolder'  ,
                      # 5:'dcd',
                      # 6:'crd',
                      # 7:'xyz',
                      # 8:'mol2',
                      # 9:'netcdf',
                      #10:'log_file'  
                           }


    def OpenWindow (self, vobject = None):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'src/gui/windows/analysis/energy_refinement_window.glade'))
            self.builder.connect_signals(self)
            #self.vobject = vobject
            self.window = self.builder.get_object('window')
            self.window.set_title('Energy Refinement Window')
            self.window.connect('destroy', self.CloseWindow)
            
            
            
            
            
            
            #---------------------------------------------------------------------------------
            #                           INPUT TYPE COMBOBOX
            #---------------------------------------------------------------------------------
            self.comobobox_input_type = self.builder.get_object('comobobox_input_type')
            for key, input_type in self.input_types.items():
                self.input_types_liststore.append([input_type])
            
            self.comobobox_input_type.set_model(self.input_types_liststore)
            self.comobobox_input_type.connect("changed", self.on_input_types_changed)
            self.comobobox_input_type.set_model(self.input_types_liststore)
            
            renderer_text = Gtk.CellRendererText()
            self.comobobox_input_type.pack_start(renderer_text, True)
            self.comobobox_input_type.add_attribute(renderer_text, "text", 0)
            #---------------------------------------------------------------------------------
           
           
           
            
            #----------------------------------------------------------------------------------------------
            #                                Starting Coordinates ComboBox  
            #----------------------------------------------------------------------------------------------
            self.box_coordinates = self.builder.get_object('box_coordinates')
            self.combobox_starting_coordinates = CoordinatesComboBox() #self.builder.get_object('coordinates_combobox')
            #self.combobox_starting_coordinates.connect("changed", self.on_name_combo_changed)
            self.box_coordinates.pack_start(self.combobox_starting_coordinates, False, False, 0)
            self._starting_coordinates_model_update(init = True)
            #----------------------------------------------------------------------------------------------
            
            
            
            #----------------------------------------------------------------------------------------------
            #                                    Reaction Coordinate Boxes  
            ##'''--------------------------------------------------------------------------------------------'''
            self.RC_box1 = ReactionCoordinateBox(self.main)
            self.builder.get_object('rc1_aligment').add(self.RC_box1)
            self.RC_box2 = ReactionCoordinateBox(self.main)
            self.builder.get_object('rc2_aligment').add(self.RC_box2)
            #'''--------------------------------------------------------------------------------------------'''
            
            
            
            
            '''--------------------------------------------------------------------------------------------'''     
            self.folder_chooser_button = FolderChooserButton(main =  self.main, sel_type = 'folder', home =  self.home)
            self.builder.get_object('folder_chooser_box').pack_start(self.folder_chooser_button.btn, True, True, 0)
            system_id      = self.p_session.active_id
            
            #------------------------------------------------------------------------------------------------
            if self.main.p_session.psystem[self.p_session.active_id]:
                if self.main.p_session.psystem[self.p_session.active_id].e_working_folder == None:
                    folder = HOME
                else:
                    folder = self.main.p_session.psystem[self.p_session.active_id].e_working_folder
                self.folder_chooser_button.set_folder(folder = folder)
            #------------------------------------------------------------------------------------------------
            
            
            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system('energy_ref')
                self.builder.get_object('entry_logfile').set_text(output_name)  
            

            self.window.show_all()
            
            self.RC_box1.set_hide_scan_parameters()
            self.RC_box2.set_hide_scan_parameters()
            
            self.RC_box1.set_rc_type(0)
            self.RC_box2.set_rc_type(0)            
            
            self.comobobox_input_type.set_active(0)
           
            self.builder.get_object('button_cancel').connect('clicked', self.CloseWindow)
            self.builder.get_object('button_run').connect('clicked', self.on_button_run_clicked)
            
            self.Visible  = True
    
        else:
            pass


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


    def on_input_types_changed(self, widget):
        """ Function doc """
        _id = self.comobobox_input_type.get_active()
        #print(_id)
        
        if _id == 0:
            self.builder.get_object('label_coordinates').show()
            self.combobox_starting_coordinates.show()
            
            self.builder.get_object('label_file_folder').hide()
            self.builder.get_object('filechooser_file_folder').hide()
            self.builder.get_object('box_reaction_coordinate2').set_sensitive(False)
        
        elif _id in [1,2]:
            self.builder.get_object('label_coordinates').hide()
            self.combobox_starting_coordinates.hide()
            
            self.builder.get_object('label_file_folder').show()
            self.builder.get_object('filechooser_file_folder').show()
            
            if _id ==1:
                self.builder.get_object('box_reaction_coordinate2').set_sensitive(False)
            elif _id ==2:
                self.builder.get_object('box_reaction_coordinate2').set_sensitive(True)
            else:
                pass
            

    def on_coordinates_combobox_change (self, widget):
        """ Function doc """
        _id = self.coordinates_combobox.get_active()
        #print(_id)
        vobject_index = None
        #-----------------------------------------------------------------------------
        _iter = self.coordinates_combobox.get_active_iter()
        if _iter is not None:
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.coordinates_combobox.get_model()
            _name, vobject_index = model[_iter][:2]
            #print ('\n\n\_name, vobject_index:', _name, vobject_index, '\n\n')
        #-----------------------------------------------------------------------------
        self.vobject = self.main.vm_session.vm_objects_dic[vobject_index]

        self.data_liststore.clear()
        for index , data in enumerate(self.main.p_session.systems[self.vobject.easyhybrid_system_id]['logfile_data'][vobject_index]):
            #print(data)
            self.data_liststore.append([data['name'], index])
        
        
        #self.data_liststore.append(['all', 2])
        
        self.data_combobox.set_active(0)


    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False


    def on_button_run_clicked (self, button):
        """ Function doc """
        parameters = {"simulation_type":"Energy_Refinement"                ,
                      "dialog":True                                        ,                      
                      "system"     : self.p_session.psystem[self.p_session.active_id],
                      "system_name": self.p_session.psystem[self.p_session.active_id].label,
                      "initial_coordinates": None                          ,                       
                      "traj_type":'pklfolder'                              ,
                      'ignore_mm_charges': False                           ,
                      "NmaxThreads":1                                      ,
                      "show":False                                         }

        
        parameters['folder'] = self.folder_chooser_button.get_folder()
        
        input_type = self.comobobox_input_type.get_active()
        print('_type: ',input_type)
        #----------------------------------------------------------------------
        if input_type == 0:
            parameters["traj_type"] = 'vobject'
                             
            tree_iter = self.combobox_starting_coordinates.get_active_iter()
            if tree_iter is not None:
                '''selecting the vismol object from the content that is in the combobox '''
                model = self.combobox_starting_coordinates.get_model()
                name, vobject_id = model[tree_iter][:2]
                vobject = self.main.vm_session.vm_objects_dic[vobject_id]
                parameters["trajectory"] = vobject.frames
            
        
        elif input_type in [1,2]:
            if input_type == 1:
                parameters["traj_type"]  = 'pklfolder'
            else:
                parameters["traj_type"]  = 'pklfolder2D'
            
            
            parameters["data_path"] = self.builder.get_object('filechooser_file_folder').get_filename()
            parameters["filename"] = self.builder.get_object('entry_logfile').get_text()
            files = os.listdir( parameters['data_path'])
            pkl_files = []
            
            for _file in files:
                # Check if the file is a text file
                if _file.endswith('.pkl'):
                    pkl_files.append(_file)

            print ('pDynamo pkl folder:' , parameters['traj_type'])
            print ('Number of pkl files:', len(pkl_files))
            parameters["trajectory"] = pkl_files
        
        else:
            pass
        #----------------------------------------------------------------------
        
        if self.builder.get_object('check_box_MM_chrg_to_zero').get_active():
            parameters["ignore_mm_charges"] = True
        else:
            pass
        
        parameters["RC1"] = self.RC_box1.get_rc_data()
        
        #if self.builder.get_object('label_check_button_reaction_coordinate2').get_active():
        if parameters["traj_type"]  == 'pklfolder2D':
            parameters["RC2"] = self.RC_box2.get_rc_data()
            parameters["NmaxThreads"] =  int(self.builder.get_object('n_CPUs_spinbutton').get_value())
            #parameters["traj_type"]   = 'pklfolder2D'
        
        else:
            parameters["RC2"] = None
            parameters["NmaxThreads"] = 1
            #parameters["traj_type"]   = 'pklfolder'
        
        #----------------------------------------------------------------------
        
        #pprint (parameters)
        self.p_session.run_simulation( parameters = parameters )


    def update (self):
        """ Function doc """
        #print('VismolGoToAtomWindow2 update')
        self._starting_coordinates_model_update()
        if self.Visible:
            self.update_working_folder_chooser()
            
            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system('energy_ref')
                self.builder.get_object('entry_logfile').set_text(output_name)  

    def update_working_folder_chooser (self, folder = None):
        """ Function doc """
        if folder:
            #print('update_working_folder_chooser')
            self.folder_chooser_button.set_folder(folder = folder)
        else:
            
            folder = self.main.p_session.psystem[self.main.p_session.active_id].e_working_folder
            if folder:
                self.folder_chooser_button.set_folder(folder = folder)
            else:
                pass

