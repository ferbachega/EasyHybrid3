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
from gui.gtk_widgets import FolderChooserButton

#from gui.windows.geometry_optimization_window import  FolderChooserButton
from pprint import pprint
VISMOL_HOME = os.environ.get('VISMOL_HOME')
HOME        = os.environ.get('HOME')



atomic_dic = {#Symbol     name         number    Cov(r)     VdW(r)     Mass
                "H"  : ["Hydrogen"     , 1   ,  0.330000 , 1.200000,  1.007940   ],
                "He" : ["Helium"       , 2   ,  0.700000 , 1.400000,  4.002602   ],
                "Li" : ["Lithium"      , 3   ,  1.230000 , 1.820000,  6.941000   ],
                "Be" : ["Beryllium"    , 4   ,  0.900000 , 1.700000,  9.012182   ],
                "B"  : ["Boron"        , 5   ,  0.820000 , 2.080000,  10.811000  ],
                "C"  : ["Carbon"       , 6   ,  0.770000 , 1.950000,  12.010700  ],
                "N"  : ["Nitrogen"     , 7   ,  0.700000 , 1.850000,  14.006700  ],
                "O"  : ["Oxygen"       , 8   ,  0.660000 , 1.700000,  15.999400  ],
                "F"  : ["Fluorine"     , 9   ,  0.611000 , 1.730000,  18.998404  ],
                "Ne" : ["Neon"         , 10  ,  0.700000 , 1.540000,  20.179701  ],
                "Na" : ["Sodium"       , 11  ,  3.06     , 2.270000,  22.989771  ],
                "Mg" : ["Magnesium"    , 12  ,  1.360000 , 1.730000,  24.305000  ],
                "Al" : ["Aluminium"    , 13  ,  1.180000 , 2.050000,  26.981539  ],
                "Si" : ["Silicon"      , 14  ,  0.937000 , 2.100000,  28.085501  ],
                "P"  : ["Phosphorus"   , 15  ,  0.890000 , 2.080000,  30.973761  ],
                "S"  : ["Sulphur"      , 16  ,  1.040000 , 2.000000,  32.064999  ],
                "Cl" : ["Chlorine"     , 17  ,  0.997000 , 1.970000,  35.452999  ],
                "Ar" : ["Argon"        , 18  ,  1.740000 , 1.880000,  39.948002  ],
                "K"  : ["Potassium"    , 19  ,  2.030000 , 2.750000,  39.098301  ],
                "Ca" : ["Calcium"      , 20  ,  1.740000 , 1.973000,  40.077999  ],
                "Sc" : ["Scandium"     , 21  ,  1.440000 , 1.700000,  44.955910  ],
                "Ti" : ["Titanium"     , 22  ,  1.320000 , 1.700000,  47.867001  ],
                "V"  : ["Vanadium"     , 23  ,  1.220000 , 1.700000,  50.941502  ],
                "Cr" : ["Chromium"     , 24  ,  1.180000 , 1.700000,  51.996101  ],
                "Mn" : ["Manganese"    , 25  ,  1.170000 , 1.700000,  54.938049  ],
                "Fe" : ["Iron"         , 26  ,  1.170000 , 1.700000,  55.845001  ],
                "Co" : ["Cobalt"       , 27  ,  1.160000 , 1.700000,  58.933201  ],
                "Ni" : ["Nickel"       , 28  ,  1.150000 , 1.630000,  58.693401  ],
                "Cu" : ["Copper"       , 29  ,  1.170000 , 1.400000,  63.546001  ],
                "Zn" : ["Zinc"         , 30  ,  1.250000 , 1.390000,  65.408997  ],
                "Ga" : ["Gallium"      , 31  ,  1.260000 , 1.870000,  69.723000  ],
                "Ge" : ["Germanium"    , 32  ,  1.188000 , 1.700000,  72.639999  ],
                "As" : ["Arsenic"      , 33  ,  1.200000 , 1.850000,  74.921600  ],
                "Se" : ["Selenium"     , 34  ,  1.170000 , 1.900000,  78.959999  ],
                "Br" : ["Bromine"      , 35  ,  1.167000 , 2.100000,  79.903999  ],
                "Kr" : ["Krypton"      , 36  ,  1.910000 , 2.020000,  83.797997  ],
                "Rb" : ["Rubidium"     , 37  ,  2.160000 , 1.700000,  85.467796  ],
                "Sr" : ["Strontium"    , 38  ,  1.910000 , 1.700000,  87.620003  ],
                "Y"  : ["Yttrium"      , 39  ,  1.620000 , 1.700000,  88.905853  ],
                "Zr" : ["Zirconium"    , 40  ,  1.450000 , 1.700000,  91.223999  ],
                "Nb" : ["Niobium"      , 41  ,  1.340000 , 1.700000,  92.906380  ],
                "Mo" : ["Molybdenum"   , 42  ,  1.300000 , 1.700000,  95.940002  ],
                "Tc" : ["Technetium"   , 43  ,  1.270000 , 1.700000,  98.000000  ],
                "Ru" : ["Ruthenium"    , 44  ,  1.250000 , 1.700000,  101.070000 ],
                "Rh" : ["Rhodium"      , 45  ,  1.250000 , 1.700000,  102.905502 ],
                "Pd" : ["Palladium"    , 46  ,  1.280000 , 1.630000,  106.419998 ],
                "Ag" : ["Silver"       , 47  ,  1.340000 , 1.720000,  107.868202 ],
                "Cd" : ["Cadmium"      , 48  ,  1.480000 , 1.580000,  112.411003 ],
                "In" : ["Indium"       , 49  ,  1.440000 , 1.930000,  114.818001 ],
                "Sn" : ["Tin"          , 50  ,  1.385000 , 2.170000,  118.709999 ],
                "Sb" : ["Antimony"     , 51  ,  1.400000 , 2.200000,  121.760002 ],
                "Te" : ["Tellurium"    , 52  ,  1.378000 , 2.060000,  127.599998 ],
                "I"  : ["Iodine"       , 53  ,  1.387000 , 2.150000,  126.904472 ],
                "Xe" : ["Xenon"        , 54  ,  1.980000 , 2.160000,  131.292999 ],
                "Cs" : ["Cesium"       , 55  ,  2.350000 , 1.700000,  132.905457 ],
                "Ba" : ["Barium"       , 56  ,  1.980000 , 1.700000,  137.326996 ],
                "La" : ["Lanthanum"    , 57  ,  1.690000 , 1.700000,  138.905502 ],
                "Ce" : ["Cerium"       , 58  ,  1.830000 , 1.700000,  140.115997 ],
                "Pr" : ["Praseodymium" , 59  ,  1.820000 , 1.700000,  140.907654 ],
                "Nd" : ["Neodymium"    , 60  ,  1.810000 , 1.700000,  144.240005 ],
                "Pm" : ["Promethium"   , 61  ,  1.800000 , 1.700000,  145.000000 ],
                "Sm" : ["Samarium"     , 62  ,  1.800000 , 1.700000,  150.360001 ],
                "Eu" : ["Europium"     , 63  ,  1.990000 , 1.700000,  151.964005 ],
                "Gd" : ["Gadolinium"   , 64  ,  1.790000 , 1.700000,  157.250000 ],
                "Tb" : ["Terbium"      , 65  ,  1.760000 , 1.700000,  158.925339 ],
                "Dy" : ["Dysprosium"   , 66  ,  1.750000 , 1.700000,  162.500000 ],
                "Ho" : ["Holmium"      , 67  ,  1.740000 , 1.700000,  164.930313 ],
                "Er" : ["Erbium"       , 68  ,  1.730000 , 1.700000,  167.259003 ],
                "Tm" : ["Thulium"      , 69  ,  1.720000 , 1.700000,  168.934204 ],
                "Yb" : ["Ytterbium"    , 70  ,  1.940000 , 1.700000,  173.039993 ],
                "Lu" : ["Lutetium"     , 71  ,  1.720000 , 1.700000,  174.966995 ],
                "Hf" : ["Hafnium"      , 72  ,  1.440000 , 1.700000,  178.490005 ],
                "Ta" : ["Tantalum"     , 73  ,  1.340000 , 1.700000,  180.947906 ],
                "W"  : ["Tungsten"     , 74  ,  1.300000 , 1.700000,  183.839996 ],
                "Re" : ["Rhenium"      , 75  ,  1.280000 , 1.700000,  186.207001 ],
                "Os" : ["Osmium"       , 76  ,  1.260000 , 1.700000,  190.229996 ],
                "Ir" : ["Iridium"      , 77  ,  1.270000 , 1.700000,  192.216995 ],
                "Pt" : ["Platinum"     , 78  ,  1.300000 , 1.720000,  195.078003 ],
                "Au" : ["Gold"         , 79  ,  1.340000 , 1.660000,  196.966553 ],
                "Hg" : ["Mercury"      , 80  ,  1.490000 , 1.550000,  200.589996 ],
                "Tl" : ["Thallium"     , 81  ,  1.480000 , 1.960000,  204.383301 ],
                "Pb" : ["Lead"         , 82  ,  1.480000 , 2.020000,  207.199997 ],
                "Bi" : ["Bismuth"      , 83  ,  1.450000 , 1.700000,  208.980377 ],
                "Po" : ["Polonium"     , 84  ,  1.460000 , 1.700000,  209.000000 ],
                "At" : ["Astatine"     , 85  ,  1.450000 , 1.700000,  210.000000 ],
                "Rn" : ["Radon"        , 86  ,  2.400000 , 1.700000,  222.000000 ],
                "Fr" : ["Francium"     , 87  ,  2.000000 , 1.700000,  223.000000 ],
                "Ra" : ["Radium"       , 88  ,  1.900000 , 1.700000,  226.000000 ],
                "Ac" : ["Actinium"     , 89  ,  1.880000 , 1.700000,  227.000000 ],
                "Th" : ["Thorium"      , 90  ,  1.790000 , 1.700000,  232.038101 ],
                "Pa" : ["Protactinium" , 91  ,  1.610000 , 1.700000,  231.035873 ],
                "U"  : ["Uranium"      , 92  ,  1.580000 , 1.860000,  238.028915 ],
                "Np" : ["Neptunium"    , 93  ,  1.550000 , 1.700000,  237.000000 ],
                "Pu" : ["Plutionium"   , 94  ,  1.530000 , 1.700000,  244.000000 ],
                "Am" : ["Americium"    , 95  ,  1.070000 , 1.700000,  243.000000 ],
                "Cm" : ["Curium"       , 96  ,  0.000000 , 1.700000,  247.000000 ],
                "Bk" : ["Berkelium"    , 97  ,  0.000000 , 1.700000,  247.000000 ],
                "Cf" : ["Californium"  , 98  ,  0.000000 , 1.700000,  251.000000 ],
                "Es" : ["Einsteinium"  , 99  ,  0.000000 , 1.700000,  252.000000 ],
                "Fm" : ["Fermium"      , 100 ,  0.000000 , 1.700000,  257.000000 ],
                "Md" : ["Mendelevium"  , 101 ,  0.000000 , 1.700000,  258.000000 ],
                "No" : ["Nobelium"     , 102 ,  0.000000 , 1.700000,  259.000000 ],
                "Lr" : ["Lawrencium"   , 103 ,  0.000000 , 1.700000,  262.000000 ],
                "Rf" : ["Rutherfordiu" , 104 ,  0.000000 , 1.700000,  261.000000 ],
                "Db" : ["Dubnium"      , 105 ,  0.000000 , 1.700000,  262.000000 ],
                "Sg" : ["Seaborgium"   , 106 ,  0.000000 , 1.700000,  263.000000 ],
                "Bh" : ["Bohrium"      , 107 ,  0.000000 , 1.700000,  264.000000 ],
                "Hs" : ["Hassium"      , 108 ,  0.000000 , 1.700000,  265.000000 ],
                "Mt" : ["Meitnerium"   , 109 ,  0.000000 , 1.700000,  268.000000 ],
                "Xx" : ["Dummy"        , 0   ,  0.000000 , 0.000000,  0.000000   ],
                "X"  : ["Dummy"        , 0   ,  0.000000 , 0.000000,  0.000000   ]
              }

texto_d1   = "\n\n                       -- simple-distance --\n\nFor simple-distance, select two atoms in pymol using the editing mode\nfollowing the diagram:\n\n   R                    R\n    \                  /\n     A1--A2  . . . . A3\n    /                  \ \n   R                    R\n         ^            ^\n         |            |\n        pk1  . . . . pk2\n                d1\n"
texto_d2d1 = "\n                       -- multiple-distance --\n\nFor multiple-distance, select three atoms in pymol using the editing mode\nfollowing the diagram:\n\n   R                    R\n    \                  /\n     A1--A2  . . . . A3\n    /                  \ \n   R                    R\n     ^   ^            ^\n     |   |            |\n    pk1-pk2  . . . . pk3\n       d1       d2\n"



class PotentialEnergyScanWindow():

    
    def __init__(self, main = None):
        """ Class initialiser """
        self.main                = main
        self.home                = main.home
        self.p_session           = self.main.p_session
        self.vm_session          = main.vm_session
        self.Visible             =  False        
        self.residue_liststore   = Gtk.ListStore(str, str, str)
        self.opt_methods = { 0 : 'ConjugatedGradient',
                             1 : 'SteepestDescent'   ,
                             2 : 'LFBGS'             ,
                             3 : 'QuasiNewton'       ,
                             4 : 'FIRE'              }
        
        self.sym_tag = 'PES_scan'

    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home,'gui/windows/PES_scan_window.glade'))
            self.builder.connect_signals(self)
            #
            self.window = self.builder.get_object('pes_scan_window')
            self.window.set_title('PES Scan Window')
            self.window.set_keep_above(True)            
            
            #self.box_reaction_coordinate2 =  self.builder.get_object('box_reaction_coordinate2')           
            self.box_reaction_coordinate2 =  self.builder.get_object('box_reaction_coordinate22')           
            #'''--------------------------------------------------------------------------------------------'''
            #'''
            self.method_store = Gtk.ListStore(str)
            
            methods = ["simple distance", "multiple distance", 'dihedral']
            
            for method in methods:
                self.method_store.append([method])
                print (method)
            
            self.combobox_reaction_coord1 = self.builder.get_object('combobox_reaction_coord1')
            self.combobox_reaction_coord1.set_model(self.method_store)
            #self.combobox_reaction_coord1.connect("changed", self.on_name_combo_changed)
            
            self.combobox_reaction_coord2 = self.builder.get_object('combobox_reaction_coord2')
            self.combobox_reaction_coord2.set_model(self.method_store)
            #
            renderer_text = Gtk.CellRendererText()
            self.combobox_reaction_coord1.pack_start(renderer_text, True)
            self.combobox_reaction_coord1.add_attribute(renderer_text, "text", 0)
            
            self.combobox_reaction_coord2.pack_start(renderer_text, True)
            self.combobox_reaction_coord2.add_attribute(renderer_text, "text", 0)
            

            #'''--------------------------------------------------------------------------------------------'''

            self.method_store = Gtk.ListStore(str)
            
            methods = [ "Conjugate Gradient" ,
                        "FIRE"               ,
                        "L-BFGS"             ,
                        "Steepest Descent"   ]
            
            for method in methods:
                self.method_store.append([method])
                print (method)
            
            self.methods_combo = self.builder.get_object('combobox_methods')
            self.methods_combo.set_model(self.method_store)
            #self.methods_combo.connect("changed", self.on_name_combo_changed)
            self.methods_combo.set_model(self.method_store)
            #
            renderer_text = Gtk.CellRendererText()
            self.methods_combo.pack_start(renderer_text, True)
            self.methods_combo.add_attribute(renderer_text, "text", 0)
            #'''--------------------------------------------------------------------------------------------'''
            self.methods_combo.set_active(0)                     
            
            '''--------------------------------------------------------------------------------------------'''
            self.combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
            #self.system_liststore = self.main.active_system_liststore
            self.combobox_starting_coordinates.set_model(self.main.vobject_liststore_dict[self.main.p_session.active_id])
            # old
            #self.combobox_starting_coordinates.set_model( self.main.active_system_liststore)
            
            renderer_text = Gtk.CellRendererText()
            self.combobox_starting_coordinates.pack_start(renderer_text, True)
            self.combobox_starting_coordinates.add_attribute(renderer_text, "text", 0)
            
            size = len(self.main.vobject_liststore_dict[self.main.p_session.active_id])
            self.combobox_starting_coordinates.set_active(size-1)
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
                #self.update_working_folder_chooser(folder = folder)            
            #------------------------------------------------------------------------------------------------
            
            
            #------------------------------------------------------------------------------------------------
            if  self.p_session.psystem[self.p_session.active_id]:
                system = self.p_session.psystem[self.p_session.active_id]
                self.builder.get_object('traj_name').set_text(str(system.e_step_counter)+'-'+system.e_tag +'_'+ self.sym_tag)
                #self.save_trajectory_box.set_filename ( )
            else:
                pass
            #------------------------------------------------------------------------------------------------
            
            #working_folder = self.p_session.systems[system_id]['working_folder']
            #self.folder_chooser_button.set_folder(folder = working_folder) 
            
            ''' 
            self.builder.get_object('label_atom4_coord1').hide()
            self.builder.get_object('entry_atom4_index_coord1').hide()
            self.builder.get_object('label_name4_coord1').hide()
            self.builder.get_object('entry_atom4_name_coord1').hide()
            
            self.builder.get_object('label_atom4_coord2').hide()
            self.builder.get_object('entry_atom4_index_coord2').hide()
            self.builder.get_object('label_name4_coord2').hide()
            self.builder.get_object('entry_atom4_name_coord2').hide()            
            ''' 
            self.window.show_all()
            
            self.builder.get_object('label_atom4_coord1').hide()
            self.builder.get_object('entry_atom4_index_coord1').hide()
            self.builder.get_object('label_name4_coord1').hide()
            self.builder.get_object('entry_atom4_name_coord1').hide()
            
            self.builder.get_object('label_atom4_coord2').hide()
            self.builder.get_object('entry_atom4_index_coord2').hide()
            self.builder.get_object('label_name4_coord2').hide()
            self.builder.get_object('entry_atom4_name_coord2').hide()            
            self.change_check_button_reaction_coordinate (None)
            #self.box_reaction_coordinate2.set_sensitive(False)
            
            

            
            self.builder.get_object('button_cancel').connect('clicked', self.CloseWindow)
            #tag  = self.p_session.psystem[self.p_session.active_id]['tag']
            #step = str(self.p_session.psystem[self.p_session.active_id]['step_counter'])
            #tag  = step+'_'+tag+'_reaction_coord_scan'  
            #self.save_trajectory_box.builder.get_object('entry_trajectory_name').set_text(tag)
            
            #self.builder.get_object('traj_name').set_text(tag)     
            self.combobox_reaction_coord1.set_active(0)
            self.combobox_reaction_coord2.set_active(0)
            self.Visible  = True

        else:
            self.window.present()
            
    def CloseWindow (self, button, data  = None):
        """ Function doc """
        self.window.destroy()
        self.Visible    =  False
    
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

        try:
            self.refresh_dmininum ( coord1 = True)
        except:
            print(texto_d1)
            print(texto_d2d1)
                    
    def change_cb_coordType2 (self, combo_box):
        """ Function doc """
        
        _type = self.combobox_reaction_coord2.get_active()        
        
        if _type == 0:
            self.builder.get_object('label_atom3_coord2').hide()
            self.builder.get_object('entry_atom3_index_coord2').hide()
            self.builder.get_object('label_name3_coord2').hide()
            self.builder.get_object('entry_atom3_name_coord2').hide()
            
            self.builder.get_object('label_atom4_coord2').hide()
            self.builder.get_object('entry_atom4_index_coord2').hide()
            self.builder.get_object('label_name4_coord2').hide()
            self.builder.get_object('entry_atom4_name_coord2').hide()
            self.builder.get_object('mass_restraints2').set_sensitive(False)

        if _type == 1:
            self.builder.get_object('label_atom3_coord2').show()
            self.builder.get_object('entry_atom3_index_coord2').show()
            self.builder.get_object('label_name3_coord2').show()
            self.builder.get_object('entry_atom3_name_coord2').show()
            
            self.builder.get_object('label_atom4_coord2').hide()
            self.builder.get_object('entry_atom4_index_coord2').hide()
            self.builder.get_object('label_name4_coord2').hide()
            self.builder.get_object('entry_atom4_name_coord2').hide()
            self.builder.get_object('mass_restraints2').set_sensitive(True)

        if _type == 2:
            self.builder.get_object('label_atom3_coord2').show()
            self.builder.get_object('entry_atom3_index_coord2').show()
            self.builder.get_object('label_name3_coord2').show()
            self.builder.get_object('entry_atom3_name_coord2').show()
            
            self.builder.get_object('label_atom4_coord2').show()
            self.builder.get_object('entry_atom4_index_coord2').show()
            self.builder.get_object('label_name4_coord2').show()
            self.builder.get_object('entry_atom4_name_coord2').show()
            self.builder.get_object('mass_restraints2').set_sensitive(False)

        #try:
        self.refresh_dmininum ( coord2 = True)

    def toggle_mass_restraint1 (self, widget):
        """ Function doc """
        self.refresh_dmininum(coord1 =  True)
    
    def toggle_mass_restraint2 (self, widget):
        """ Function doc """
        self.refresh_dmininum(coord2 =  True)
    
    def refresh_dmininum (self, coord1 =  False, coord2 = False):
        """ Function doc """
        
        if coord1:
            _type = self.combobox_reaction_coord1.get_active()
            print('_type', _type)
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
        else:
            pass    
              
        if coord2:
            _type = self.combobox_reaction_coord2.get_active()
            try:
                if _type == 0:
                    index1 = int(self.builder.get_object('entry_atom1_index_coord2').get_text() )
                    index2 = int(self.builder.get_object('entry_atom2_index_coord2').get_text() )

                    dist1 = get_distance(self.vobject, index1, index2 )
                    self.builder.get_object('entry_dmin_coord2').set_text(str(dist1))
                if _type == 1:
                    index1 = int(self.builder.get_object('entry_atom1_index_coord2').get_text() )
                    index2 = int(self.builder.get_object('entry_atom2_index_coord2').get_text() )
                    index3 = int(self.builder.get_object('entry_atom3_index_coord2').get_text() )
                    
                    dist1 = get_distance(self.vobject, index1, index2 )
                    dist2 = get_distance(self.vobject, index2, index3 )
                    
                    if self.builder.get_object('mass_restraints2').get_active():
                        self.sigma_pk1_pk3_rc2, self.sigma_pk3_pk1_rc2  = compute_sigma_a1_a3(self.vobject, index1, index3)
                        #print('distance a1 - a2:', dist1 - dist2)
                        DMINIMUM =  (self.sigma_pk1_pk3_rc2 * dist1) -(self.sigma_pk3_pk1_rc2 * dist2*-1)
                        self.builder.get_object('entry_dmin_coord2').set_text(str(DMINIMUM))
                    else:
                        DMINIMUM =  dist1- dist2
                        self.builder.get_object('entry_dmin_coord2').set_text(str(DMINIMUM))
            except:
                print(texto_d1)
                print(texto_d2d1)
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
            
        if widget == self.builder.get_object('import_picking_selection_button1'):
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
            
        else:
            if atom1:
                self.builder.get_object('entry_atom1_index_coord2').set_text(str(atom1.index-1) )
                self.builder.get_object('entry_atom1_name_coord2' ).set_text(str(atom1.name) )
            else: print('use picking selection to chose the central atom')
            #-------
            if atom2:
                self.builder.get_object('entry_atom2_index_coord2').set_text(str(atom2.index-1) )
                self.builder.get_object('entry_atom2_name_coord2' ).set_text(str(atom2.name) )
            else: print('use picking selection to chose the central atom')            
            #-------
            if atom3:
                self.builder.get_object('entry_atom3_index_coord2').set_text(str(atom3.index-1) )
                self.builder.get_object('entry_atom3_name_coord2' ).set_text(str(atom3.name) )
                
                if atom3.symbol == atom1.symbol:
                    self.builder.get_object('mass_restraints2').set_active(False)
                else:
                    self.builder.get_object('mass_restraints2').set_active(True)
            
            else: print('use picking selection to chose the central atom')           
            #-------
            if atom4:
                self.builder.get_object('entry_atom4_index_coord2').set_text(str(atom4.index-1) )
                self.builder.get_object('entry_atom4_name_coord2' ).set_text(str(atom4.name) )
            else: print('use picking selection to chose the central atom')
    
            self.refresh_dmininum( coord2 =  True, )

    def change_check_button_reaction_coordinate (self, widget):
        """ Function doc """
        #radiobutton_bidimensional = self.builder.get_object('radiobutton_bidimensional')
        if self.builder.get_object('label_check_button_reaction_coordinate2').get_active():
            self.box_reaction_coordinate2.set_sensitive(True)
            self.builder.get_object('n_CPUs_spinbutton').set_sensitive(True)
            self.builder.get_object('n_CPUs_label').set_sensitive(True)
        else:
            self.box_reaction_coordinate2.set_sensitive(False)
            self.builder.get_object('n_CPUs_spinbutton').set_sensitive(False)
            self.builder.get_object('n_CPUs_label')     .set_sensitive(False)
        #print(widget)

    #def update_working_folder_chooser (self, folder = None):
    #    """ Function doc """
    #    if folder:
    #        #print('update_working_folder_chooser')
    #        self.save_trajectory_box.set_folder(folder = folder)
    #    else:
    #        self.save_trajectory_box.set_folder(folder = self.p_session.psystem[self.p_session.active_id]['working_folder'])
   
    def run_dialog (self, text = None, secondary_text = None):
        """ Function doc """
        
        if text == None:
            text="Folder not found"
        if secondary_text == None:
            secondary_text = "The folder you have selected does not appear to be valid. Please select a different folder or create a new one."
        
        dialog = Gtk.MessageDialog(
                    transient_for=self.main.window,
                    message_type=Gtk.MessageType.ERROR,
                    buttons=Gtk.ButtonsType.OK,
                    text=text,
                    )
        dialog.format_secondary_text(secondary_text)
        dialog.run()
        dialog.destroy()

    #======================================================================================
    def run_scan(self,button):
        '''
        Get infotmation and run the simulation
        '''         
        parameters = {"simulation_type":"Relaxed_Surface_Scan"             ,
                      "ndim":1                                             ,
                      "dialog":True                                        ,                      
                      "system":self.p_session.psystem[self.p_session.active_id],
                      "system_name":self.p_session.psystem[self.p_session.active_id].label,
                      "initial_coordinates": None                          ,                       
                      "traj_type":'pklfolder'                              ,
                      "second_coordinate": False                           ,
                      "ATOMS_RC1":None                                     ,
                      "ATOMS_RC2":None                                     ,
                      "ATOMS_RC1_NAMES":None                               ,
                      "ATOMS_RC2_NAMES":None                               ,                      
                      "nsteps_RC1":0                                       ,
                      "nsteps_RC2":0                                       ,
                      "force_constant_1":4000.0                            ,
                      "force_constant_2":4000.0                            ,
                      "maxIterations":1000                                 ,
                      "dincre_RC1":0.1                                     ,
                      "dincre_RC2":0.1                                     ,
                      "dminimum_RC1":0.0                                   ,
                      "dminimum_RC2":0.0                                   ,
                      "sigma_pk1pk3_rc1":1.0                               ,
                      "sigma_pk3pk1_rc1":-1.0                              ,
                      "sigma_pk1pk3_rc2":1.0                               ,
                      "sigma_pk3pk1_rc2":-1.0                              ,
                      "rc_type_1":"simple_distance"                        ,
                      "rc_type_2":"simple_distance"                        ,
                      "adaptative":False                                   , 
                      "save_format":".dcd"                                 ,
                      "rmsGradient":0.1                                    ,
                      "optimizer":"ConjugatedGradient"                     ,
                      "MC_RC1":False                                       ,
                      "MC_RC2":False                                       ,
                      "log_frequency":50                                   ,
                      "contour_lines":10                                   ,
                      "NmaxThreads":1                                      ,
                      "show":False                                         }
        
        parameters["optimizer"]        = self.opt_methods[self.methods_combo.get_active()]
        parameters["folder"]           = self.folder_chooser_button.get_folder()        
        parameters["maxIterations"]    = float(self.builder.get_object('entry_max_int').get_text() )
        parameters["rmsGradient"]      = float(self.builder.get_object('entry_rmsd_tol').get_text() )
        parameters["traj_folder_name"] = self.builder.get_object('traj_name').get_text()        
        parameters["vobject_name"]     = self.builder.get_object('traj_name').get_text()        
        
        combobox_starting_coordinates = self.builder.get_object('combobox_starting_coordinates')
        print(combobox_starting_coordinates)
        tree_iter = combobox_starting_coordinates.get_active_iter()
        if tree_iter is not None:
            '''selecting the vismol object from the content that is in the combobox '''
            model = combobox_starting_coordinates.get_model()
            name, vobject_id = model[tree_iter][:2]
            vobject = self.main.vm_session.vm_objects_dic[vobject_id]
            
            '''This function imports the coordinates of a vobject into the dynamo system in memory.''' 
            #print('vobject:', vobject.name, len(vobject.frames) )
            self.main.p_session.get_coordinates_from_vobject_to_pDynamo_system(vobject)
        
        #----------------------------------------------------------------------------------               
        parameters["initial_coordinates"] = vobject.name
        #----------------------------------------------------------------------------------               
        
        _type = self.combobox_reaction_coord1.get_active()
        print('_type', _type)
        if _type == 0:
            parameters["rc_type_1"]     = "simple_distance"
            index1 = int( self.builder.get_object('entry_atom1_index_coord1').get_text() )
            index2 = int( self.builder.get_object('entry_atom2_index_coord1').get_text() )
            
            name1 = self.builder.get_object('entry_atom1_name_coord1').get_text()
            name2 = self.builder.get_object('entry_atom2_name_coord1').get_text()
            
            dmin   = float( self.builder.get_object('entry_dmin_coord1').get_text( ))
            parameters["ATOMS_RC1"]       = [ index1, index2 ] 
            parameters["ATOMS_RC1_NAMES"] = [ name1 ,  name2 ] 
            parameters["dminimum_RC1"]  = dmin 
        
        elif _type == 1:
            parameters["rc_type_1"]     = "multiple_distance"
            index1 = int( self.builder.get_object('entry_atom1_index_coord1').get_text() )
            index2 = int( self.builder.get_object('entry_atom2_index_coord1').get_text() )
            index3 = int( self.builder.get_object('entry_atom3_index_coord1').get_text() )
            
            name1 = self.builder.get_object('entry_atom1_name_coord1').get_text() 
            name2 = self.builder.get_object('entry_atom2_name_coord1').get_text() 
            name3 = self.builder.get_object('entry_atom3_name_coord1').get_text() 
            
            dmin   = float( self.builder.get_object('entry_dmin_coord1').get_text( ))
            parameters["ATOMS_RC1"]       = [ index1, index2, index3 ] 
            parameters["ATOMS_RC1_NAMES"] = [ name1,  name2,  name3 ] 
            parameters["dminimum_RC1"]  = dmin  
            
            if self.builder.get_object('mass_restraints1').get_active():
                parameters["MC_RC1"] = True
                parameters["sigma_pk1pk3_rc1"] = self.sigma_pk1_pk3 
                parameters["sigma_pk3pk1_rc1"] = self.sigma_pk3_pk1 
        
        
        elif _type == 2:
            index1 = int( self.builder.get_object('entry_atom1_index_coord1').get_text() )
            index2 = int( self.builder.get_object('entry_atom2_index_coord1').get_text() )
            index3 = int( self.builder.get_object('entry_atom3_index_coord1').get_text() )
            index4 = int( self.builder.get_object('entry_atom4_index_coord1').get_text() )
            dmin   = float( self.builder.get_object('entry_dmin_coord1').get_text() )
            parameters["ATOMS_RC1"]     = [ index1,index2,index3,index4] 
            parameters["dminimum_RC1"]  = dmin 
            parameters["rc_type_1"]     = "dihedral"
        
        parameters["nsteps_RC1"]        = int( self.builder.get_object('entry_nsteps1').get_text() )
        parameters["force_constant_1"]  = float( self.builder.get_object('entry_FORCE_coord1').get_text() )
        parameters["dincre_RC1"]        = float( self.builder.get_object('entry_step_size1').get_text() )
        #----------------------------------------------------------------------------------
        '''
        coordenada  de reacao 2
        '''        
        if self.builder.get_object('label_check_button_reaction_coordinate2').get_active():
            parameters["second_coordinate"] = True
            parameters["NmaxThreads"] =  int(self.builder.get_object('n_CPUs_spinbutton').get_value())
            #parameters["nprocs"] =  int(self.builder.get_object('n_CPUs_spinbutton').get_value())
            self.is_scan2d = True
            parameters["ndim"] = 2
            parameters["traj_type"] = 'pklfolder2D'
            _type = self.combobox_reaction_coord2.get_active()            
            entry_FORCE_coord2 = int(self.builder.get_object('entry_FORCE_coord2').get_text() )
            
            
            #--------------------------------------------------------------------------
            if _type == 0: # simple
                parameters["rc_type_2"]     = "simple_distance"
                index1 = int(self.builder.get_object('entry_atom1_index_coord2').get_text() )
                index2 = int(self.builder.get_object('entry_atom2_index_coord2').get_text() )
                
                name1 = self.builder.get_object('entry_atom1_name_coord2').get_text() 
                name2 = self.builder.get_object('entry_atom2_name_coord2').get_text()
                
                dmin2  = float(self.builder.get_object('entry_dmin_coord2').get_text( ))
                parameters["ATOMS_RC2"]       = [ index1, index2 ]
                parameters["ATOMS_RC2_NAMES"] = [ name1,  name2] 

                parameters["dminimum_RC2"]  = dmin2 
            
            
            #--------------------------------------------------------------------------
            elif _type == 1: # multiple
                parameters["rc_type_2"]     = "multiple_distance"
                index1 = int(self.builder.get_object('entry_atom1_index_coord2').get_text() )
                index2 = int(self.builder.get_object('entry_atom2_index_coord2').get_text() )
                index3 = int(self.builder.get_object('entry_atom3_index_coord2').get_text() )
                
                name1 = self.builder.get_object('entry_atom1_name_coord2').get_text() 
                name2 = self.builder.get_object('entry_atom2_name_coord2').get_text() 
                name3 = self.builder.get_object('entry_atom3_name_coord2').get_text() 
                
                dmin2  = float(self.builder.get_object('entry_dmin_coord2').get_text( ))
                parameters["ATOMS_RC2"]     = [ index1, index2, index3 ]
                parameters["ATOMS_RC2_NAMES"] = [ name1,  name2,  name3 ] 

                parameters["dminimum_RC2"]  = dmin2 
                if self.builder.get_object('mass_restraints2').get_active():
                    parameters["MC_RC1"] = True
                    parameters["sigma_pk1pk3_rc2"] = self.sigma_pk1_pk3_rc2 
                    parameters["sigma_pk3pk1_rc2"] = self.sigma_pk3_pk1_rc2
                else: mass_weighted_2 = False 
            
            #--------------------------------------------------------------------------
            elif _type == 2: # dihedral
                parameters["rc_type_2"]     = "dihedral"
                index1 = int(self.builder.get_object('entry_atom1_index_coord2').get_text() )
                index2 = int(self.builder.get_object('entry_atom2_index_coord2').get_text() )
                index3 = int(self.builder.get_object('entry_atom3_index_coord2').get_text() )
                index4 = int(self.builder.get_object('entry_atom4_index_coord2').get_text() )
                dmin2   = float(self.builder.get_object('entry_dmin_coord2').get_text( ))
                parameters["ATOMS_RC2"]     = [ index1, index2, index3, index4 ] 
                parameters["dminimum_RC2"]  = dmin2 
        #_-----------------------------------------------------------------------------------
        parameters["nsteps_RC2"]        = int( self.builder.get_object('entry_nstep2').get_text() )
        parameters["force_constant_2"]  = float( self.builder.get_object('entry_FORCE_coord2').get_text() )
        parameters["dincre_RC2"]        = float( self.builder.get_object('entry_step_size2').get_text() )
        #-----------------------------------------------------------------------------------
        pprint(parameters)
        #a = input('')
        
        
        isExist = os.path.exists(parameters['folder'])
        if isExist:
            pass
        else:
            self.run_dialog()
            return None
        
        
        isExist = os.path.exists(os.path.join( parameters['folder'], parameters['traj_folder_name']+".ptGeo"))
        if isExist:

            self.run_dialog(text = 'Trajectory name not valid!', secondary_text ='Please provide a new trajectory name as the specified folder already exists' )
            return None
        else:
            pass
        
        #here!
        self.p_session.run_simulation( parameters = parameters )
        

    def update (self):
        """ Function doc """
        pass
    





def get_distance (vobject, index1, index2):
    """ Function doc """
    print( index1, index2)
    atom1 = vobject.atoms[index1]
    atom2 = vobject.atoms[index2]
    a1_coord = atom1.coords()
    a2_coord = atom2.coords()
    
    dx = a1_coord[0] - a2_coord[0]
    dy = a1_coord[1] - a2_coord[1]
    dz = a1_coord[2] - a2_coord[2]
    dist = (dx**2+dy**2+dz**2)**0.5
    print('distance a1 - a2:', dist)
    return dist
#===========================================================
def get_angle (vobject, index1, index2, index3):
    """ Function doc """
#===========================================================
def get_dihedral (vobject, index1, index2, index3, index4):
    """ Function doc """
#===========================================================
def compute_sigma_a1_a3 (vobject, index1, index3):

    """ example:
        pk1 ---> pk2 ---> pk3
         N  ---   H  ---  O	    
         
         where H is the moving atom
         calculation only includes N and O ! 
    """
    atom1 = vobject.atoms[index1]
    atom3 = vobject.atoms[index3]
    
    symbol1 = atom1.symbol
    symbol3 = atom3.symbol    
    
    mass1 = atomic_dic[symbol1][4]
    mass3 = atomic_dic[symbol3][4]
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
