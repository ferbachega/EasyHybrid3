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
import os, sys, time
import gi 
import signal
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk, Pango
from gi.repository import GdkPixbuf


from gui.widgets.custom_widgets  import VismolSelectionTypeBox
from gui.widgets.custom_widgets  import FileChooser
from gui.widgets.custom_widgets  import get_colorful_square_pixel_buffer
from gui.widgets.custom_widgets  import ReactionCoordinateBox
from gui.widgets.custom_widgets  import SequenceViewerBox

from gui.windows.setup.windows_and_dialogs import ImportANewSystemWindow
from gui.windows.setup.windows_and_dialogs import EasyHybridDialogSetQCAtoms
from gui.windows.setup.windows_and_dialogs import EasyHybridSetupQCModelWindow
from gui.windows.setup.windows_and_dialogs import EasyHybridGoToAtomWindow
#from gui.windows.setup.windows_and_dialogs import PDynamoSelectionWindow
from gui.windows.setup.windows_and_dialogs import EasyHybridSelectionWindow
from gui.windows.setup.windows_and_dialogs import ExportDataWindow
from gui.windows.setup.windows_and_dialogs import EasyHybridDialogPrune
from gui.windows.setup.windows_and_dialogs import MakeSolventBoxWindow


from gui.windows.setup.windows_and_dialogs import ImportTrajectoryWindow
from gui.windows.setup.windows_and_dialogs import TrajectoryPlayerWindow
from gui.windows.setup.windows_and_dialogs import InfoWindow
from gui.windows.setup.windows_and_dialogs import MergeSystemWindow
from gui.windows.setup.windows_and_dialogs import SolvateSystemWindow
from gui.windows.setup.windows_and_dialogs import SimpleDialog
from gui.windows.setup.edit_frames_dialog import EditFrameDialog

from gui.windows.setup.easyhybrid_terminal    import TerminalWindow
from gui.windows.setup.selection_list_window  import *
from gui.windows.setup.setup_interface        import EasyHybridPreferencesWindow
from gui.windows.setup.process_manager_window import ProcessManagerWindow

from gui.windows.simulation.single_point_window          import SinglePointWindow
from gui.windows.simulation.geometry_optimization_window import GeometryOptimization
from gui.windows.simulation.PES_scan_window              import PotentialEnergyScanWindow 
from gui.windows.simulation.PES_advanced_scan_window     import AdvancedPotentialEnergyScanWindow 
from gui.windows.simulation.molecular_dynamics_window    import MolecularDynamicsWindow 
from gui.windows.simulation.umbrella_sampling_window     import UmbrellaSamplingWindow 
from gui.windows.simulation.chain_of_states_opt_window   import ChainOfStatesOptWindow 
from gui.windows.simulation.normal_modes_window          import NormalModesWindow 

from gui.windows.analysis.WHAM_analysis_window                    import WHAMWindow 
from gui.windows.analysis.normal_modes_analysis_window            import NormalModesAnalysisWindow 
from gui.windows.analysis.surface_analysis_window                 import SurfaceAnalysisWindow 
#from gui.windows.analysis.surface_list_window                     import SurfaceListWindow 
from gui.windows.analysis.energy_refinement_window                import EnergyRefinementWindow
from gui.windows.analysis.PES_analysis_window                     import PotentialEnergyAnalysisWindow
from gui.windows.analysis.distance_angle_dihedral_analysis_window import DistanceAngleDihedralAnalysisWindow
from gui.windows.analysis.RMSD_tool                               import RMSDToolWindow
from gui.windows.analysis.RMSD_analysis_window                    import RMSDAnalysisWindow #/home/fernando/programs/EasyHybrid3/src/gui/windows/analysis/RMSD_analysis_window.py
from gui.windows.analysis.align_trajectory                        import AlignTrajectoryWindow #/home/fernando/programs/EasyHybrid3/src/gui/windows/analysis/RMSD_analysis_window.py
from gui.windows.analysis.reimaging_trajectory                    import ReimagingTrajectoryWindow #/home/fernando/programs/EasyHybrid3/src/gui/windows/analysis/RMSD_analysis_window.py

from util.geometric_analysis import get_simple_distance
from util.sequence_plot import GtkSequenceViewer
from util.rama_plot import RamachandranWindow


from pdynamo.pDynamo2EasyHybrid import pDynamoSession
import numpy as np



from pCore                     import Align                                        , \
                                      Clone                                        , \
                                      logFile                                      , \
                                      Selection                                    , \
                                      TestScript_InputDataPath                     , \
                                      TestScript_OutputDataPath                    , \
                                      XHTMLLogFileWriter

class PreferencesWindow:
    """ Class doc """
    
    def __init__ (self, main = None, e_id = None , v_id = None):
        """ Class initialiser """
        
        self.xml='''
<?xml version="1.0" encoding="UTF-8"?>
<!-- Generated with glade 3.22.2 -->
<interface>
  <requires lib="gtk+" version="3.20"/>
  <object class="GtkWindow" id="window">
    <property name="can_focus">False</property>
    <property name="default_width">300</property>
    <child type="titlebar">
      <placeholder/>
    </child>
    <child>
      <object class="GtkAlignment">
        <property name="visible">True</property>
        <property name="can_focus">False</property>
        <property name="top_padding">5</property>
        <property name="bottom_padding">5</property>
        <property name="left_padding">5</property>
        <property name="right_padding">5</property>
        <child>
          <object class="GtkBox">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="orientation">vertical</property>
            <property name="spacing">5</property>
            <child>
              <object class="GtkGrid">
                <property name="visible">True</property>
                <property name="can_focus">False</property>
                <property name="row_spacing">5</property>
                <property name="column_spacing">5</property>
                <child>
                  <object class="GtkEntry" id="entry_name">
                    <property name="visible">True</property>
                    <property name="can_focus">True</property>
                    <property name="hexpand">True</property>
                    <property name="text" translatable="yes">new pDynamo system</property>
                  </object>
                  <packing>
                    <property name="left_attach">1</property>
                    <property name="top_attach">0</property>
                  </packing>
                </child>
                <child>
                  <object class="GtkLabel" id="label_name">
                    <property name="visible">True</property>
                    <property name="can_focus">False</property>
                    <property name="halign">end</property>
                    <property name="label" translatable="yes">Name:</property>
                  </object>
                  <packing>
                    <property name="left_attach">0</property>
                    <property name="top_attach">0</property>
                  </packing>
                </child>
                <child>
                  <object class="GtkLabel" id="label_tag">
                    <property name="visible">True</property>
                    <property name="can_focus">False</property>
                    <property name="halign">end</property>
                    <property name="label" translatable="yes">Tag: </property>
                  </object>
                  <packing>
                    <property name="left_attach">0</property>
                    <property name="top_attach">1</property>
                  </packing>
                </child>
                <child>
                  <object class="GtkEntry" id="entry_tag">
                    <property name="visible">True</property>
                    <property name="can_focus">True</property>
                    <property name="max_length">15</property>
                    <property name="width_chars">10</property>
                    <property name="text" translatable="yes">molsys</property>
                  </object>
                  <packing>
                    <property name="left_attach">1</property>
                    <property name="top_attach">1</property>
                  </packing>
                </child>
                <child>
                  <object class="GtkLabel" id="label_color">
                    <property name="visible">True</property>
                    <property name="can_focus">False</property>
                    <property name="label" translatable="yes">Color:</property>
                  </object>
                  <packing>
                    <property name="left_attach">0</property>
                    <property name="top_attach">2</property>
                  </packing>
                </child>
                <child>
                  <object class="GtkColorButton" id="button_color">
                    <property name="visible">True</property>
                    <property name="can_focus">True</property>
                    <property name="receives_default">True</property>
                    <property name="rgba">rgb(138,226,52)</property>
                  </object>
                  <packing>
                    <property name="left_attach">1</property>
                    <property name="top_attach">2</property>
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
              <object class="GtkButtonBox" id="dialog-action_area2">
                <property name="visible">True</property>
                <property name="can_focus">False</property>
                <property name="layout_style">end</property>
                <child>
                  <object class="GtkButton" id="button_cancel">
                    <property name="label" translatable="yes">Cancel</property>
                    <property name="visible">True</property>
                    <property name="can_focus">True</property>
                    <property name="receives_default">True</property>
                  </object>
                  <packing>
                    <property name="expand">False</property>
                    <property name="fill">False</property>
                    <property name="position">0</property>
                  </packing>
                </child>
                <child>
                  <object class="GtkButton" id="button_apply">
                    <property name="label" translatable="yes">Apply</property>
                    <property name="visible">True</property>
                    <property name="can_focus">True</property>
                    <property name="receives_default">True</property>
                  </object>
                  <packing>
                    <property name="expand">False</property>
                    <property name="fill">False</property>
                    <property name="position">1</property>
                  </packing>
                </child>
              </object>
              <packing>
                <property name="expand">False</property>
                <property name="fill">False</property>
                <property name="position">2</property>
              </packing>
            </child>
          </object>
        </child>
      </object>
    </child>
  </object>
</interface>
'''        
        
        self.main       = main
        self.vm_session = main.vm_session
        self.p_session  = main.p_session
        
        self.e_id  = e_id 
        self.v_id  = v_id 
         
        
        self.builder = Gtk.Builder()
        self.builder.add_from_string(self.xml)

        self.window = self.builder.get_object('window')

        self.entry_name   = self.builder.get_object('entry_name')
        self.entry_tag    = self.builder.get_object('entry_tag')
        self.button_color = self.builder.get_object('button_color')
        
        self.button_cancel = self.builder.get_object('button_cancel')
        self.button_apply  = self.builder.get_object('button_apply')
        
        self.button_apply.connect('clicked', self.on_button_apply )
        self.button_cancel.connect('clicked', self.on_button_cancel )
        self.main.main_treeview.treeview_menu.rename_window_visible = True
        
        #self.window.connect("destroy", self.on_button_cancel)
        
        self.window.set_resizable(False)
        self.window.show_all()
        self.button_color.hide()
        self.builder.get_object('label_color').hide()
    
    def set_names (self, name, tag):
        """ Function doc """
        
        self.entry_name.set_text(name)
        self.entry_tag .set_text(tag)
        
    
    def on_button_apply (self, button):
        """ Function doc """
        name  = self.entry_name.get_text()
        tag   = self.entry_tag .get_text()
    
        self.rename_window_visible = False
        self.main.rename_tag(e_id = self.e_id, tag = tag)
        self.main.rename(e_id  = self.e_id, 
                         v_id  = self.v_id, 
                         name  = name)
        self.window.destroy()
        self.main.main_treeview.treeview_menu.rename_window_visible = False
    
    def on_button_cancel (self, button):
        """ Function doc """
        self.window.destroy()
        self.main.main_treeview.treeview_menu.rename_window_visible = False
