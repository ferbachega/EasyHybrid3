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
#from GTKGUI.gtkWidgets.filechooser import FileChooser
#from easyhybrid.pDynamoMethods.pDynamo2Vismol import *
import gc
import os

VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')


class ImportANewSystemWindow(Gtk.Window):
    """ Class doc """
    
    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'gui/windows/windows_and_dialogs.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('ImportNewSystemWindow')
            self.window.set_border_width(10)
            self.window.set_default_size(500, 370)  

            
            '''--------------------------------------------------------------------------------------------'''
            self.system_type_store = Gtk.ListStore(str)
            system_types = [
                "AMBER",
                "CHARMM",
                "OPLS",
                "pdynamo files (*.pkl, *.yaml)",
                "other (*.pdb, *.xyz, *.mol2)",
                ]
            for system_type in system_types:
                self.system_type_store.append([system_type])
                print (system_type)
            
            self.system_types_combo = Gtk.ComboBox.new_with_model(self.system_type_store)
            #self.system_types_combo = self.builder.get_object('system_type_combox_from_import_a_new_system')
            self.box_combo = self.builder.get_object('box')
            self.box_combo.pack_start(self.system_types_combo, True, True, 0)
            
            self.system_types_combo.connect("changed", self.on_name_combo_changed)
            self.system_types_combo.set_model(self.system_type_store)
            
            renderer_text = Gtk.CellRendererText()
            self.system_types_combo.pack_start(renderer_text, True)
            self.system_types_combo.add_attribute(renderer_text, "text", 0)
            '''--------------------------------------------------------------------------------------------'''
           
            
            self.treeview = self.builder.get_object('gtktreeview_import_system')
            for i, column_title in enumerate(['file', "type", "number of atoms"]):
                renderer = Gtk.CellRendererText()
                column = Gtk.TreeViewColumn(column_title, renderer, text=i)
                self.treeview.append_column(column)

            
            self.window.show_all()                                               
            self.builder.connect_signals(self)                                   
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()

            self.Visible  =  True
            
            self.files    = {
                            'amber_prmtop': None,
                            'charmm_par'  : [],
                            'charmm_psf'  : None,
                            'charmm_extra': None, 
                            'opls_folder' : [],
                            'coordinates' : None,
                            }
            self.system_types_combo.set_active(0)

            #----------------------------------------------------------------

    def CloseWindow (self, button, data  = None):
        """ Function doc """
        #self.BackUpWindowData()
        self.window.destroy()
        self.Visible    =  False
        #print('self.Visible',self.Visible)
    
    def __init__(self, main = None, home = None):
        """ Class initialiser """
        
        if home:
            self.home = home
        else:
            self.home = os.environ.get('VISMOL_HOME')
        
        
        
        self.easyhybrid_main     = main
        self.Visible             =  False        
        
        self.residue_liststore = Gtk.ListStore(str, str, str)
        #self.atom_liststore    = Gtk.ListStore(bool, int, str, str, int, int)
        #self.residue_filter    = False
        
        self.charmm_txt = '''For systems prepared in the traditional CHARMM / PSF format, the required files are: parameters (prm / par), topologies (psf) and coordinates (chm, crd, pdb, xyz, ...)

NOTE: You can include more than one parameter file if needed.'''
        self.amber_txt = '''Systems prepared using the AMBER force field require  two files: topologies (top / prmtop) and coordinates (crd, pdb, xyz, ...).  '''
        self.OPLS_txt = '''Systems prepared natively in pDynamo using the OPLS force field require: A parameter folder and a topology/coordinate file (pdb, mol2, mol, ...) .  '''
        self.gmx_txt = '''Systems prepared natively in GROMACS using the CHARMM(or AMBER) force field require: A parameter/topology (top) file and coordinate file (pdb, mol2, mol, ...) .  '''
        
    
    def on_name_combo_changed(self, widget):
        """ Function doc """
        fftype = self.system_types_combo.get_active()
        print (fftype)
        self.files    = {
            'amber_prmtop': None,
            'charmm_par'  : [],
            'charmm_psf'  : None,
            'charmm_extra': None, 
            'opls_folder' : [],
            'coordinates' : None,
            }
        self.residue_liststore = Gtk.ListStore(str, str, str)
        self.treeview.set_model(self.residue_liststore)
            
        if fftype == 0: #AMBER
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()
            self.builder.get_object('gtk_label_fftype').set_text(self.amber_txt)
 
            
        if fftype == 1: # "CHARMM":
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()
            self.builder.get_object('gtk_label_fftype').set_text(self.charmm_txt)

            
        if fftype == 10: #"GROMACS":
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()
            self.builder.get_object('gtk_label_fftype').set_text(self.gmx_txt)

            
        if fftype == 2:#"OPLS":
            self.builder.get_object('gtkbox_OPLS_folderchooser').show()
            self.builder.get_object('gtk_label_fftype').set_text(self.OPLS_txt)
            
        if fftype == 3: #"pDynamo files(*.pkl,*.yaml)":
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()
            
        if fftype == 4: #"Other(*.pdb,*.xyz,*.mol2...)":
            self.builder.get_object('gtkbox_OPLS_folderchooser').hide()


    def filetype_parser(self, filein, systemtype):
        filetype = self.GetFileType(filein)
        print (filetype, systemtype)
        if filetype in ['top', 'prmtop', 'TOP', 'PRMTOP', ]:
            if systemtype == 0:
                self.files['amber_prmtop'] = filein
                return 'amber parameters/topologies'
                
   
        
        elif filetype in ['par', 'prmtop', 'prm', 'PAR', 'PRM', 'str']:
            if systemtype == 1:
                self.files['charmm_par'].append(filein)
                return 'charmm parameters'
        
        
        
        
        elif filetype in ['psf', 'PSF','psfx', 'PSFX']:
            if systemtype == 1:
                self.files['charmm_psf'] = filein
                return 'charmm topologies'
        
        
        
        elif filetype in ['pdb', 'PDB','mol','MOL','mol2','MOL2', 'xyz', 'XYZ', 'crd', 'inpcrd', 'chm']:
            #if systemtype == 1:
            self.files['coordinates'] = filein
            return 'coordinates'
        
        
        
        elif filetype in ['pkl', 'PKL']:
            #if systemtype == 1:
            self.files['coordinates'] = filein
            return 'pDynamo coordinates'
        else:
            return 'unknow'
        

    
    
    def GetFileType(self, filename):
        file_type = filename.split('.')
        return file_type[-1]



    def on_delete_files_button_clicked (self, button):
        """ Function doc """
        self.residue_liststore = Gtk.ListStore(str, str, str)
        self.treeview.set_model(self.residue_liststore)
    
    def on_import_files_button_clicked (self, button):
        """ Function doc """
        files = self.easyhybrid_main.filechooser.open(select_multiple = True)
        #print(files)

        for _file in files:
            #for res in self.VObj.chains[chain].residues:
                ##print(res.resi, res.resn, chain,  len(res.atoms) ) 
            systemtype = self.system_types_combo.get_active()
            filetype = self.filetype_parser( _file, systemtype)
            self.residue_liststore.append(list([_file, filetype, '10' ]))
        self.treeview.set_model(self.residue_liststore)
        self.files['opls_folder'] =  self.builder.get_object('OPLS_folderchooserbutton').get_filename()
        #print(self.files)


    def on_button_import_a_new_system_clicked (self, button):
        """ Function doc """
        
        if button == self.builder.get_object('ok_button_import_a_new_system'):
            print('ok_button_import_a_new_system')
            #self.on_button1_clicked_create_new_project(button)
            #self.dialog.hide()
        if button == self.builder.get_object('cancel_button_import_a_new_system'):
            print('cancel_button_import_a_new_system')
            self.dialog.hide()
            
    def on_button4_import_system_clicked (self, button):
        #print('ok_button_import_a_new_system')
        systemtype = self.system_types_combo.get_active()
        
        name =  self.builder.get_object('entry_system_name').get_text()
        tag  =  self.builder.get_object('entry_system_tag').get_text()
        print(self.files, systemtype)
        
        #'''
        self.easyhybrid_main.p_session.load_a_new_pDynamo_system_from_dict(input_files = self.files, 
                                                                           system_type = systemtype, 
                                                                           name        = name      ,
                                                                           tag         = tag
                                                                           )
        #'''
        #if systemtype == 2:
        #    self.files['opls_folder'] =  self.builder.get_object('OPLS_folderchooserbutton').get_filename()
        
        #print ('systemtype',systemtype, self.files )
        ##self.easyhybrid_main.p_session.get_bonds_from_pDynamo_system()
        #
        #vobject = self.easyhybrid_main.p_session.build_vobject_from_pDynamo_system (name = name)
        self.CloseWindow(button, data  = None)

