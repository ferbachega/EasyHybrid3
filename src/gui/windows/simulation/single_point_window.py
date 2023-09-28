import gi
import os
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Pango
from gui.widgets.custom_widgets import FolderChooserButton
from gui.widgets.custom_widgets import CoordinatesComboBox


class SinglePointWindow:
    """ Class doc """
    
    def __init__(self, main = None):
        """ Class initialiser """
        self.main    = main
        self.p_session    = main.p_session
        self.home    = main.home
        self.Visible =  False        
        self.starting_coords_liststore = Gtk.ListStore(str, int)
        self.is_single_frame  = True
    
    def OpenWindow (self, sys_selected = None):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'src/gui/windows/simulation/single_point_window.glade'))
            self.builder.connect_signals(self)
            
            
            self.window = self.builder.get_object('window')
            self.window.set_title('Single Point Window')
            self.window.set_keep_above(True)
            self.window.set_default_size(400, 200)
            '''--------------------------------------------------------------------------------------------'''
            
            
            '''--------------------------------------------------------------------------------------------'''     
            self.folder_chooser_button = FolderChooserButton(main =  self.main, 
                                                         sel_type = 'folder', 
                                                             home =  self.home)
                                                             
            self.builder.get_object('folder_chooser_box').pack_start(self.folder_chooser_button.btn, True, True, 0)
            system_id      = self.p_session.active_id
            
            #------------------------------------------------------------------------------------------------
            if self.main.p_session.psystem[self.p_session.active_id]:
                if self.main.p_session.psystem[self.p_session.active_id].e_working_folder == None:
                    folder = HOME
                else:
                    folder = self.main.p_session.psystem[self.p_session.active_id].e_working_folder
                self.folder_chooser_button.set_folder(folder = folder)
                #self.update_working_folder_chooser(folder = folder)            
            #------------------------------------------------------------------------------------------------
            
            
            
            # - - - - - - - - - - - - - Starting Coordinates ComboBox - - - - - - - - - - - - - - - - -
            '''--------------------------------------------------------------------------------------------'''
            #self.combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
            ##---------------------------------------------------------------------------------------------
            #self._starting_coordinates_model_update(init = True)
            #renderer_text = Gtk.CellRendererText()
            #self.combobox_starting_coordinates.pack_start(renderer_text, True)
            #self.combobox_starting_coordinates.add_attribute(renderer_text, "text", 0)
            #----------------------------------------------------------------------------------------------
            
            #----------------------------------------------------------------------------------------------
            self.box_coordinates = self.builder.get_object('box_coordinates')
            self.combobox_starting_coordinates = CoordinatesComboBox() #self.builder.get_object('coordinates_combobox')
            #self.combobox_starting_coordinates.connect("changed", self.on_name_combo_changed)
            self.box_coordinates.pack_start(self.combobox_starting_coordinates, False, False, 0)
            self._starting_coordinates_model_update(init = True)
            #----------------------------------------------------------------------------------------------
            
            
            
            
            
            self.builder.get_object('button_cancel').connect('clicked', self.CloseWindow)
            self.builder.get_object('button_run').connect('clicked', self.on_button_run_clicked)
            self.window.connect('destroy', self.CloseWindow)

            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system('single_point')
                self.builder.get_object('entry_logfile_name').set_text(output_name)


            self.window.show_all()
            self.Visible  = True
    
        else:
            self.window.present()
    
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False

    def on_button_run_clicked (self, button):
        """ Function doc """
        parameters={ "simulation_type":"Energy_Single_Point"   ,
                     "filename"       : None            , 
                     "folder"         : os.getcwd()     , 
                     }
        
        frame = int(self.builder.get_object('entry_frame').get_text())
        #----------------------------------------------------------------------------------
        tree_iter = self.combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            
            '''selecting the vismol object from the content that is in the combobox '''
            model = self.combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]
            vobject = self.main.vm_session.vm_objects_dic[vobject_id]
            
            '''This function imports the coordinates of a vobject into the dynamo system in memory.''' 
            self.main.p_session.get_coordinates_from_vobject_to_pDynamo_system(vobject = vobject, frame = frame)
        
        
        parameters["folder"]          = self.folder_chooser_button.get_folder()
        parameters["filename"] = self.builder.get_object('entry_logfile_name').get_text()

        
        
        #------------------------------------------------------------------#
        #                      RUN ENERGY CALCULATION                      #
        #------------------------------------------------------------------#
        energy = self.p_session.run_simulation( parameters = parameters )
        #------------------------------------------------------------------#
        
        
        
        #------------------------------------------------------------------#
        #                          DIALOG MESSAGE                          #
        #------------------------------------------------------------------#
        dialog = Gtk.MessageDialog(
            transient_for=None,
            flags=0,
            message_type=Gtk.MessageType.INFO,
            buttons=Gtk.ButtonsType.OK,
            text="Energy(kJ/mol): "+str(energy),
        )
        
        

        psystem = parameters['system']
            
        string = ''
        name    = psystem.label
        size    = len(psystem.atoms)
        string = '\nsystem: {}    \natoms: {}    '.format(name, size)

        if psystem.qcModel:
            hamiltonian   = getattr(psystem.qcModel, 'hamiltonian', False)
            
            if hamiltonian:
                pass
            else:
                try:
                    itens = psystem.qcModel.SummaryItems()
                    hamiltonian = itens[0][0]
                except:
                    hamiltonian = 'external'
            n_QC_atoms    = len(list(psystem.qcState.pureQCAtoms))
            summary_items = psystem.electronicState.SummaryItems()
            
            string += '\nhamiltonian: {}    \nQC atoms: {}    \nQC charge: {}    \nspin multiplicity {}    '.format(  hamiltonian, 
                                                                                                              n_QC_atoms,
                                                                                                              summary_items[1][1],
                                                                                                              summary_items[2][1],
                                                                                                                 )
                
        n_fixed_atoms = len(psystem.e_fixed_table)
        string += '\nfixed atoms: {}    '.format(n_fixed_atoms)
            
        if psystem.mmModel:
            forceField = psystem.mmModel.forceField
            string += '\nforceField: {}    '.format(forceField)
        
            if psystem.nbModel:
                nbmodel = psystem.mmModel.forceField
                string += 'nbModel: True    '
                
                summary_items = psystem.nbModel.SummaryItems()
                
            
            else:
                string += 'nbModel: False    '
                
                
            if psystem.symmetry:
                #nbmodel = psystem.mmModel.forceField
                string += '\nsymmetry: {}    '.format( psystem.symmetry.crystalSystem.label)
                #print(psystem.symmetry)
                #print(psystem.symmetryParameters)
                
                #summary_items = psystem.nbModel.SummaryItems()
                    
                
            else:
                string += ''
            
        
        string += 'frame: {}'.format (frame)
        string += '\n\nfile: {}'.format(os.path.join(parameters["folder"],parameters["filename"]+'.log' ))
        
        
        dialog.format_secondary_text(
            string
        )
        dialog.run()
        dialog.destroy()
        #--------------------------------------------------------------------
    
    
    
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

    def update (self, parameters = None):
        """ Function doc """
        self._starting_coordinates_model_update()
        if self.Visible:
            self.update_working_folder_chooser()
            
            if  self.p_session.psystem[self.p_session.active_id]:
                output_name = self.p_session.get_output_filename_from_system('single_point')
                self.builder.get_object('entry_logfile_name').set_text(output_name)

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

