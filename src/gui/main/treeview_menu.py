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
from util.debug import dprint
import os, sys, time
import gi 
import signal
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk, Pango
from gi.repository import GdkPixbuf


from gui.widgets.custom_widgets  import VismolSelectionTypeBox
import vismol.utils.mesh_decimation as mesh_decimation
import vismol.utils.mesh_smoothing as mesh_smoothing
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
from gui.windows.setup.windows_and_dialogs import TextWindow, TabbedLogWindow
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
# [EN] Used by _surf_setup() below (the per-Vobject surface setup
# dialog): recolor_surface_lobe/recolor_mep_surface do the cheap,
# recompute-free recolouring described in their own docstrings in
# surface_analysis_window.py; COLOR_MAPS is the same colormap dictionary
# already used to populate the MEP colormap combobox in the main
# surface-generation window, reused here for consistency (same names,
# same actual gradients).
from gui.windows.analysis.surface_analysis_window                 import recolor_surface_lobe, recolor_mep_surface
from util.colormaps                                                import COLOR_MAPS
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

# --- imports entre modulos adicionados na refatoracao ---
from gui.main.preferences_window import PreferencesWindow

class TreeViewMenu:
    """ Class doc """
    
    def __init__ (self, treeview):
        """ Class initialiser """
        self.treeview = treeview
        self.main     = treeview.main 
        self.filechooser   = FileChooser()
        self.rename_window_visible =  False
        
        #menu_items = {
        #    "Abrir": self._menu_rename,
        #    "Salvar": self._menu_rename,
        #    "_separator": None,
        #
        #    "Exportar": {
        #        "Como PNG": self._menu_rename,
        #        "Como JPG": self._menu_rename,
        #        "_separator": None,
        #        "Avançado": {
        #            "Alta Qualidade": self._menu_rename,
        #            "Baixa Qualidade": self._menu_rename,
        #        }
        #    },
        #
        #    "header": None
        #}
        
        
        
        surf_menu_items = {
                            'header'                : None    ,
                            '_separator'            : ''      ,
                            'Rename'                : self._menu_rename ,
                            '_separator'            : ''      ,
                            'Setup'                 : self._surf_setup ,
                            }
        
        vobject_menu_items = {
                                'header'                : None    ,

                                
                                '_separator'            : ''      ,
                                'Rename'                : self._menu_rename ,
                                '_separator'            : ''      ,
                                'Show / Hide Cell'      : self.show_or_hide_cell,
                                
                                # [EN] User request: add Cartoon as a
                                # selectable representation option in
                                # EasyHybrid. No existing menu here let
                                # you pick a representation type at all
                                # (only Rename/Frames/Go To Atom/Export/
                                # Delete existed) -- this submenu covers
                                # the representations already wired into
                                # VismolObject.create_representation()
                                # (see vismol_object.py), Cartoon
                                # included (now that its secondary-
                                # structure bug is fixed -- see
                                # cartoon_BCK.py). Each entry toggles
                                # that representation on/off for the
                                # object that was right-clicked -- see
                                # _menu_toggle_representation() below.
                                
                                #'''
                                #'Representation': {
                                #                  'Lines'      : lambda mi, rep='lines'     : self._menu_toggle_representation(mi, rep),
                                #                  'Sticks'     : lambda mi, rep='sticks'    : self._menu_toggle_representation(mi, rep),
                                #                  'Nonbonded'  : lambda mi, rep='nonbonded' : self._menu_toggle_representation(mi, rep),
                                #                  'Dots'       : lambda mi, rep='dots'      : self._menu_toggle_representation(mi, rep),
                                #                  '_separator' : None,
                                #                  'Cartoon'    : lambda mi, rep='cartoon'   : self._menu_toggle_representation(mi, rep),
                                #                  },
                                ##'''
                                
                                '_separator'            : ''      ,
                                'Frames': {
                                        'Edit': {
                                             "Size": self.call_editframe_window,
                                             "Interolate": self.call_interpolate,
                                             },
                                        "_separator": None,
                                        'Current': {
                                             "Delete": self.call_delete_current_frame,
                                             "_separator": None,
                                             "Extract": self.call_extract_current_frame,
                                             "_separator": None,
                                             "Copy": self.call_copy_current_frame,
                                             },
                                        },
                                
                                
                                #'Edit Frames'           : self.call_editframe_window,
                                'Go To Atom'            : self._menu_go_to_atom ,
                                '_separator'            : ''      ,
                                # [NOVO] Mostrar/esconder a parte MM (tudo que
                                # nao esta' na lista QC). Util para focar na
                                # regiao QC sem perder o resto do sistema.
                                
                                
                                #'MM region': {
                                #        'Hide MM atoms': self._menu_hide_mm,
                                #        'Show MM atoms': self._menu_show_mm,
                                #        },
                                #'_separator'            : ''      ,
                                
                                
                                # [EN] User request: link a simulation's
                                # log file to the vismol_object IT
                                # created, right here in the main
                                # treeview's own right-click menu -- not
                                # just in the Process Manager's own
                                # treeview (see process_manager_window.
                                # py's _open_log_for_row(), added for the
                                # exact same underlying need). Made
                                # possible by an existing hook that
                                # already stores the WHOLE job results
                                # dict (including 'logfile') directly on
                                # the object simulations_mixin.
                                # _handle_result() creates
                                # (`vobject.results = results`) -- see
                                # _menu_view_log() below for how this is
                                # read back out, and why objects with no
                                # such .results (loaded from a file,
                                # drawn in the Builder, ...) get a plain,
                                # friendly "no log" message instead of a
                                # crash.
                                'View Log'              : self._menu_view_log ,
                                '_separator'            : ''      ,
                                'Export As...'          : self._menu_export_data_window ,
                                '_separator'            : ''      ,
                                'Delete'                : self._menu_delete_vm_object ,

                                }

        system_menu_items = {
                                
                                'header'                  : None                            ,
                                                          
                                '_separator'              : ''                              ,
                                                          
                                'Info'                    : self._show_info                  ,
                                                          
                                '_separator'              : ''                              ,
                                                          
                                'Rename'                  : self._menu_rename               ,
                                'Import Data...'          : self._menu_load_data_to_system  ,
                                'Reference Color'         : self._menu_change_color_palette ,
                                #'Edit Parameters'        : self.f2                         ,
                                'Export As...'            : self._menu_export_data_window    ,
                                
                                '_separator'              : ''                              ,
                                'Merge With...'           : self._menu_merge_system         ,
                                'Clone System'            : self._menu_clone_system         ,
                                
                                '_separator'              : ''                              ,
                                
                                'Delete'                  : self._menu_delete_system        ,
                                #'test'  : self.f1 ,
                                #'f1'    : self.f1 ,
                                #'f2'    : self.f2 ,
                                #'gordão': self.f3 ,
                                #'delete': self.f3 ,
                                }
                    
                    
                    
        self.tree_view_surf_menu  , self.tree_header_surf_menu    = self.build_tree_view_menu(surf_menu_items)
        self.tree_view_vobj_menu  , self.tree_header_vobj_menu    = self.build_tree_view_menu(vobject_menu_items)
        self.tree_view_sys_menu   , self.tree_header_sys_menu     = self.build_tree_view_menu(system_menu_items)

    def show_or_hide_cell (self, widget):
        """ Function doc """
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()

        old_name = model.get_value(iter, 2)
        v_id     = model.get_value(iter, 1)
        e_id     = model.get_value(iter, 0)
        #----------------------------------------------------------------------        
        #system = self.main.p_session.psystem[e_id]
        
        vobject = self.main.vm_session.vm_objects_dic[v_id]
        if "cell_lines" in vobject.representations.keys():
            if vobject.representations["cell_lines"].active:
                self.main.vm_session.hide_cell (vobject)
            
            else:
                self.main.vm_session.show_cell (vobject)
        else:
            self.main.vm_session.show_cell (vobject)
            
    def call_editframe_window (self, widget):
        """ Function doc """
        selection        = self.treeview.get_selection()
        (model, iter)    = selection.get_selected()
        e_id             = int(model.get_value(iter, 0))
        vm_object_index  = int(model.get_value(iter, 1))
        self.main.edit_frames_dialog.open_window (vm_object_index)
    
    def call_interpolate (self, widget):
        """
        Interpolate trajectory frames by inserting midpoint frames between
        each pair of consecutive frames in the selected VisMol object.

        This method retrieves the currently selected molecular object from the
        trajectory tree view, accesses its coordinate frames, and generates a
        new trajectory with interpolated frames. For every pair of consecutive
        frames, a new intermediate frame is computed using linear interpolation
        (t = 0.5), effectively doubling the temporal resolution of the trajectory.

        The resulting interpolated trajectory is then stored in a newly created
        VisMol object, which is added to the current EasyHybrid session using
        the name ``'edited_coords'``.

        After the new object is created, the method applies fixed-atom and
        quantum chemistry (QC) visual representations, then refreshes the main
        GUI widgets to reflect the updated number of frames.

        Parameters
        ----------
        widget : Gtk.Widget
            GTK widget that triggered the callback (typically a button).

        Returns
        -------
        None
            This method does not return any value. It modifies the application
            state by creating a new interpolated trajectory object.

        Notes
        -----
        - Interpolation is performed using simple linear interpolation:

          ``C = (1 - t) * frame1 + t * frame2``

          where ``t = 0.5`` (midpoint interpolation).

        - The first frame is copied unchanged into the new trajectory.

        - If the original trajectory contains `N` frames, the resulting
          trajectory will contain:

          ``2 * N - 1`` frames

          since one interpolated frame is inserted between every pair of
          consecutive frames.

        Side Effects
        ------------
        - Creates a new VisMol object in the current session.
        - Updates molecular representations.
        - Refreshes GUI components.
        - Prints progress information to stdout.

        Examples
        --------
        When a trajectory contains 3 frames:

            Frame0 ---- Frame1 ---- Frame2

        The interpolated trajectory becomes:

            Frame0 -- Mid01 -- Frame1 -- Mid12 -- Frame2

        Resulting in 5 total frames.
        """
        
        selection        = self.treeview.get_selection()
        (model, iter)    = selection.get_selected()
        e_id             = int(model.get_value(iter, 0))
        vm_object_index  = int(model.get_value(iter, 1))
        #self.main.edit_frames_dialog.open_window (vm_object_index)
        
        #vobject = self.main_session.vm_session.vm_objects_dic[vm_object_index]
        vobject = self.main.vm_session.vm_objects_dic[vm_object_index]
        frames  = vobject.frames

        
        atom_qtty = len(vobject.atoms.items())
        size = len(vobject.frames)
        dprint(atom_qtty, size)
        dprint('Interpolating, wait a second…')
        
        #coords
        
        
        init_frame = frames[0]
        
        
        new_traje  = np.empty([1, int(atom_qtty), 3], dtype=np.float32)
        
        for j, xyz in enumerate(init_frame):
            new_traje[0][j][0] = init_frame[j][0]
            new_traje[0][j][1] = init_frame[j][1]
            new_traje[0][j][2] = init_frame[j][2]
        
        dprint('adding:', 0)
        
        for i in range(0, len(frames)-1):
            frame1 = frames[i]
            frame2 = frames[i+1]
            
            new_frame  = np.empty([2, int(atom_qtty), 3], dtype=np.float32)
            t = 0.5
            
            C = (1 - t) * frame1 + t * frame2
            
            #'''
            for j, xyz in enumerate(frame1):
                dx  = (frame1[j][0] - frame2[j][0])/2
                dy  = (frame1[j][1] - frame2[j][1])/2
                dz  = (frame1[j][2] - frame2[j][2])/2

                
                new_frame[0][j][0] = C[j][0]  
                new_frame[0][j][1] = C[j][1]  
                new_frame[0][j][2] = C[j][2]  
                
                new_frame[1][j][0] = frame2[j][0] 
                new_frame[1][j][1] = frame2[j][1] 
                new_frame[1][j][2] = frame2[j][2] 
            #'''

            
            new_traje = np.vstack((new_traje, new_frame))
       
        #system = self.main_session.p_session.psystem[vobject.e_id]
        system = self.main.p_session.psystem[vobject.e_id]
        #vobject = self.main_session.p_session._add_vismol_object_to_easyhybrid_session (system, show_molecule=True, name = 'edited_coords')
        vobject = self.main.p_session._add_vismol_object_to_easyhybrid_session (system, show_molecule=True, name = 'edited_coords')
        vobject.frames = new_traje

        # Apply fixed representation to the VisMol object
        #self.main_session.p_session._apply_fixed_representation_to_vobject(vismol_object =vobject)
        self.main.p_session._apply_fixed_representation_to_vobject(vismol_object =vobject)
        
        # Apply QC representation to the VisMol object
        #self.main_session.p_session._apply_QC_representation_to_vobject(vismol_object =vobject)
        self.main.p_session._apply_QC_representation_to_vobject(vismol_object =vobject)
        
        # Refresh the widgets in the main window
        #self.main_session.main_treeview.refresh_number_of_frames()
        self.main.main_treeview.refresh_number_of_frames()
        #self.main_session.p_session.main.refresh_widgets()        
        #self.close_window(None, None)

    def call_delete_current_frame(self, widget):
        """
        Delete the currently selected molecular frame from the trajectory
        of the selected VisMol object.

        The function removes the active frame (or the last frame if the current
        index is out of bounds) from the object's trajectory and updates the UI
        and OpenGL view.
        """

        selection = self.treeview.get_selection()
        model, iter_ = selection.get_selected()

        vm_object_index = int(model.get_value(iter_, 1))
        vobject = self.main.vm_session.vm_objects_dic[vm_object_index]

        frames = vobject.frames
        frame_state = self.main.vm_session.get_frame()

        size = len(frames)

        # Debug (optional)
        atom_qtty = len(vobject.atoms)
        #print(atom_qtty, size, type(frames))

        # Safe frame index
        frame_index = min(frame_state, size - 1)

        # NumPy deletion (clean + explicit copy)
        vobject.frames = np.delete(frames, frame_index, axis=0)

        # UI updates
        self.treeview.refresh_number_of_frames()
        self.treeview.refresh_trajectory_scalebar()

        # OpenGL refresh
        self.main.vm_session.vm_glcore.queue_draw()

    def call_extract_current_frame (self, widget):
        """ Function doc """
        dprint('extrac_current_frame / Under construction')
        
        self.call_copy_current_frame(None)
        self.call_delete_current_frame(None)

    def call_copy_current_frame(self, widget):
        """
        Copy the currently selected molecular frame and create a new VisMol object
        containing only that frame.
        """

        selection = self.treeview.get_selection()
        model, iter_ = selection.get_selected()

        vm_object_index = int(model.get_value(iter_, 1))
        vobject = self.main.vm_session.vm_objects_dic[vm_object_index]

        frames = vobject.frames
        frame_state = self.main.vm_session.get_frame()

        atom_qtty = len(vobject.atoms)

        # Clamp frame index safely
        frame_index = min(frame_state, len(frames) - 1)
        init_frame = frames[frame_index]

        # FAST: vectorized copy instead of Python loop
        new_traje = np.empty((1, atom_qtty, 3), dtype=np.float32)
        new_traje[0] = np.asarray(init_frame, dtype=np.float32)

        system = self.main.p_session.psystem[vobject.e_id]

        new_vobject = self.main.p_session._add_vismol_object_to_easyhybrid_session(
            system,
            show_molecule=True,
            name='edited_coords'
        )

        new_vobject.frames = new_traje

        # Apply representations
        self.main.p_session._apply_fixed_representation_to_vobject(
            vismol_object=new_vobject
        )

        self.main.p_session._apply_QC_representation_to_vobject(
            vismol_object=new_vobject
        )

        self.main.main_treeview.refresh_number_of_frames()
    
    def _show_info (self, widget):
        """ Function doc """
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        e_id          = int(model.get_value(iter, 0))  # @+
        #----------------------------------------------------------------------        
        system = self.main.p_session.psystem[e_id]
        window = InfoWindow(system)
        
    def _menu_export_data_window (self,vobject = None ):
        """ Function doc """
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        e_id          = int(model.get_value(iter, 0)) 
        #print(e_id)
        self.treeview.main.export_data_window.open_window(sys_selected = e_id)
        
    def _menu_load_data_to_system (self, vobject = None ):
        """ Function doc """
        selection        = self.treeview.get_selection()
        model, iter      = selection.get_selected()
        #print (list(model))
        self.main.import_trajectory_window.open_window(sys_selected = model.get_value(iter, 0))

    def _menu_change_color_palette (self, widget):
        """ Function doc """
        #selection               = self.selections[self.current_selection]
        self.colorchooserdialog = Gtk.ColorChooserDialog()
        
        if self.colorchooserdialog.run() == Gtk.ResponseType.OK:
            color = self.colorchooserdialog.get_rgba()
            #print(color.red,color.green, color.blue )
            new_color = [color.red, color.green, color.blue]

        else:
            new_color = False
        
        
    
        if new_color:
            self.colorchooserdialog.destroy()

            #----------------------------------------------------------------------
            selection     = self.treeview.get_selection()
            (model, iter) = selection.get_selected()
            self.selectedID  = int(model.get_value(iter, 0))  # @+
            #----------------------------------------------------------------------
            
            system = self.main.p_session.psystem[self.selectedID]
            
            self.main.change_reference_color(system, new_color)

    def _menu_merge_system (self, widget):
        """ Function doc """
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        e_id          = int(model.get_value(iter, 0)) 
        self.main.merge_system_window.selected_system_id = e_id
        self.main.merge_system_window.open_window()
    
    def _menu_clone_system (self, widget):
        """ Function doc """
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        e_id          = int(model.get_value(iter, 0))
        name          = model.get_value(iter, 2) 
        #tag           = self.main.p_session.psystem[e_id].tag
        #--------------------------------------------------------------                                                            
        dialog = EasyHybridDialogPrune(main         = self.main,    
                                       num_of_atoms = 'all'    ,    
                                       name         = name     ,    
                                       tag          = 'UNK'    ,    
                                       e_id         = e_id     ,    
                                       _type        = 1        )    
        name         = dialog.name        
        tag          = dialog.tag  
        color        = dialog.color 
        vobject_id   = dialog.vobject_id
        #--------------------------------------------------------------                                                            

        
        vobject      = self.main.vm_session.vm_objects_dic[vobject_id]
        self.main.p_session.set_psystem_coordinates_from_vobject(vobject   = vobject,
                                                                           system_id =  e_id  )
        #print(e_id)
        self.main.p_session.clone_system( e_id    = e_id, 
                                          vobject = vobject, 
                                          name    = name, 
                                          tag     = tag, 
                                          color   = color)
        self._save_backup_file()
    
    def _menu_go_to_atom (self, vobject = None):
        """ Function doc """
        ##print('f2')
        #self._show_lines(vobject = self.vobjects[0], indices = [0,1,2,3,4] )
        self.treeview.main.go_to_atom_window.open_window()
        #self.treeview.vm_session.go_to_atom_window.open_window()

    # ----------------------------------------------------------------------- #
    #  [NOVO] Mostrar / esconder a regiao MM (tudo que nao esta' na lista QC)   #
    # ----------------------------------------------------------------------- #
    def _get_mm_indexes(self, vismol_object, system):
        """Indices dos atomos da regiao MM = todos os atomos do objeto que NAO
        estao na lista QC (system.qcState.pureQCAtoms / e_qc_table).

        Retorna [] se o sistema nao tem modelo QC (nesse caso 'MM vs QC' nao
        se aplica -- tudo e' MM)."""
        try:
            if getattr(system, "qcModel", None):
                qc_table = set(system.qcState.pureQCAtoms)
            else:
                qc_table = set()
        except Exception:
            qc_table = set()
        all_indexes = list(vismol_object.atoms.keys())
        return [i for i in all_indexes if i not in qc_table]

    def _menu_set_mm_visibility(self, show):
        """Mostra (show=True) ou esconde (show=False) os atomos da regiao MM do
        objeto que foi clicado com o botao direito na treeview.

        Usa o mesmo mecanismo de representacao por selecao usado no resto do
        codigo (create_new_selection -> selecting_by_indexes -> show_or_hide).
        A regiao QC nao e' tocada.
        """
        try:
            main = self.treeview.main
            vobject_index = getattr(self, "vobject_index", None)
            if vobject_index is None or vobject_index == -1:
                return
            vismol_object = main.vm_session.vm_objects_dic[vobject_index]
            system = main.p_session.psystem[vismol_object.e_id]

            mm_indexes = self._get_mm_indexes(vismol_object, system)
            if not mm_indexes:
                return  # nada a fazer (sem regiao MM distinta)

            selection = main.vm_session.create_new_selection()
            selection.selecting_by_indexes(vismol_object, mm_indexes, clear=True)

            # Aplica a mesma visibilidade as representacoes de "corpo" tipicas
            # da parte MM. 'lines' e 'sticks' cobrem o caso comum; se alguma
            # nao existir para o objeto, show_or_hide simplesmente ignora.
            for rep in ("lines"):#, "sticks", "nonbonded", "dots"):
                try:
                    main.vm_session.show_or_hide(rep_type=rep, selection=selection, show=show)
                except Exception:
                    pass

            main.vm_session.vm_glcore.queue_draw()
        except Exception as e:
            print("MM visibility toggle failed:", e)

    def _menu_hide_mm(self, vobject = None):
        """Esconde os atomos da parte MM (mantem so' a regiao QC visivel)."""
        self._menu_set_mm_visibility(show=False)

    def _menu_show_mm(self, vobject = None):
        """Mostra novamente os atomos da parte MM (QC + MM)."""
        self._menu_set_mm_visibility(show=True)
    
    def f3 (self, vobject = None):
        """ Function doc """
        
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()


        self.selectedID  = int(model.get_value(iter, 1))  # @+
        
        
        
        del self.treeview.vm_session.vobjects_dic[self.selectedID]
        '''
        vobject = self.treeview.vm_session.vobjects_dic.pop(self.selectedID)
        del vobject
        '''
        self.treeview.store.clear()
        for vobj_index ,vis_object in self.treeview.vm_session.vobjects_dic.items():
            data = [vis_object.active          , 
                    str(vobj_index),
                    vis_object.name            , 
                    str(len(vis_object.atoms)) , 
                    str(len(vis_object.frames)),
                   ]
            model.append(data)
        self.treeview.vm_session.glwidget.queue_draw()

    def _surf_setup (self, menu_item = None ):
        """ [EN] Opens a small setup dialog scoped to the ONE surface
        object that was right-clicked (self.vobject_index, set by
        open_menu() just above) -- addresses the user's report directly
        ("quando optamos entre lines ou triangulos... altera a
        representacao de TODAS as superficies"): every control here acts
        on this single vismol_object only, never on
        vm_session.vm_objects_dic as a whole.

        Detects vismol_object.surface_type (set when the surface was
        first generated -- see on_render_button() in
        surface_analysis_window.py) and only shows the controls that
        actually apply to that type:
          - always: Wireframe/Solid, Opacity, Smooth shading (these three
            are plain rendering flags on SurfaceRepresentation -- see
            representations.py -- so applying them needs no recompute at
            all, for any surface type).
          - orbital / density / potential (the types with a flat colour
            per lobe): Color (+) / Color (-) buttons, applied via
            recolor_surface_lobe() (cheap -- no recompute).
          - mep: Colormap / vmin / vmax, applied via recolor_mep_surface()
            (also cheap -- reuses the cached, pre-colormap potential
            values instead of re-running GridDensity/Isosurface/
            GridPotential/build_potential_interpolator; see that
            function's own docstring in surface_analysis_window.py, and
            the user's original request: "quero mexer no vmin e vmax do
            MEP... todo o calculo tem que ser refeito, isso nao e bom"). """
        vismol_object = self.treeview.main.vm_session.vm_objects_dic[self.vobject_index]
        vm_glcore     = self.treeview.main.vm_session.vm_glcore
        surface_type  = getattr ( vismol_object, "surface_type", None )

        def get_surface_reps ( ):
            """ Every representation of THIS object that actually supports
            the render-mode/alpha/shading setters (i.e. every
            SurfaceRepresentation it has -- normally exactly one, or two
            for a plus/minus pair, since both lobes are separate
            representation instances sharing the same vismol_object). """
            return [ r for r in vismol_object.representations.values ( )
                     if hasattr ( r, "set_render_mode" ) ]

        reps = get_surface_reps ( )
        # Estado atual (pra a janela abrir refletindo o que ja esta na
        # tela, nao sempre nos valores padrao) -- lido da primeira
        # representacao encontrada; todas as representacoes do mesmo
        # objeto devem estar em sincronia, ja que so podem ter sido
        # mudadas por esta mesma janela (uma instancia por objeto).
        current_render_mode    = reps[0].render_mode       if reps else "surface"
        current_alpha          = reps[0].alpha             if reps else 1.0
        current_smooth_shading = reps[0].smooth_shading    if reps else False

        window = Gtk.Window ( title = "Surface Setup -- {}".format ( vismol_object.name ) )
        window.set_border_width ( 10 )
        window.set_default_size ( 260, -1 )
        window.set_keep_above ( True )
        self._surf_setup_window = window   # mantem referencia viva (padrao ja usado por self.preferences etc. neste arquivo)

        vbox = Gtk.Box ( orientation = Gtk.Orientation.VERTICAL, spacing = 8 )
        window.add ( vbox )

        label_type = Gtk.Label ( label = "Type: {}".format ( surface_type or "unknown" ) )
        label_type.set_xalign ( 0 )
        vbox.pack_start ( label_type, False, False, 0 )
        vbox.pack_start ( Gtk.Separator ( orientation = Gtk.Orientation.HORIZONTAL ), False, False, 4 )

        # --- controles comuns a QUALQUER tipo de superficie ---
        chk_wireframe = Gtk.CheckButton ( label = "Wireframe" )
        chk_wireframe.set_active ( current_render_mode == "lines" )
        def on_wireframe_toggled ( w ):
            mode = "lines" if w.get_active ( ) else "surface"
            for rep in get_surface_reps ( ):
                rep.set_render_mode ( mode )
            vm_glcore.queue_draw ( )
        chk_wireframe.connect ( "toggled", on_wireframe_toggled )
        vbox.pack_start ( chk_wireframe, False, False, 0 )

        label_opacity = Gtk.Label ( label = "Opacity:" )
        label_opacity.set_xalign ( 0 )
        scale_opacity = Gtk.Scale.new_with_range ( Gtk.Orientation.HORIZONTAL, 0, 100, 1 )
        scale_opacity.set_value ( current_alpha * 100.0 )
        scale_opacity.set_value_pos ( Gtk.PositionType.RIGHT )
        def on_opacity_changed ( w ):
            alpha = w.get_value ( ) / 100.0
            for rep in get_surface_reps ( ):
                rep.set_alpha ( alpha )
            vm_glcore.queue_draw ( )
        scale_opacity.connect ( "value-changed", on_opacity_changed )
        vbox.pack_start ( label_opacity, False, False, 0 )
        vbox.pack_start ( scale_opacity, False, False, 0 )

        chk_smooth = Gtk.CheckButton ( label = "Smooth shading" )
        chk_smooth.set_active ( current_smooth_shading )
        def on_smooth_toggled ( w ):
            mode = "smooth" if w.get_active ( ) else "flat"
            for rep in get_surface_reps ( ):
                rep.set_shading_mode ( mode )
            vm_glcore.queue_draw ( )
        chk_smooth.connect ( "toggled", on_smooth_toggled )
        vbox.pack_start ( chk_smooth, False, False, 0 )

        vbox.pack_start ( Gtk.Separator ( orientation = Gtk.Orientation.HORIZONTAL ), False, False, 4 )

        # --- decimacao (Vertex Clustering -- ver vismol/utils/mesh_decimation.pyx) ---
        # Simplifica a malha agrupando vertices dentro de uma grade espacial
        # de tamanho 'cell_size' (mesma unidade da molecula, tipicamente
        # Angstrom). Aplicado direto em vismol_object.surface_trajectory
        # (todos os frames/lobulos), substituindo a malha original -- gerar
        # a superficie de novo (menu do sistema QC/analysis) desfaz, ja que
        # a malha "cheia" original nao e' mantida em paralelo.
        label_decimate = Gtk.Label ( label = "Decimate (merge vertices within, in \u00c5):" )
        label_decimate.set_xalign ( 0 )
        vbox.pack_start ( label_decimate, False, False, 0 )

        hbox_decimate = Gtk.Box ( orientation = Gtk.Orientation.HORIZONTAL, spacing = 6 )
        entry_cell_size = Gtk.Entry ( )
        entry_cell_size.set_placeholder_text ( "e.g. 0.3" )
        entry_cell_size.set_width_chars ( 8 )
        btn_decimate = Gtk.Button ( label = "Apply" )
        label_decimate_status = Gtk.Label ( label = "" )
        label_decimate_status.set_xalign ( 0 )

        def on_decimate_clicked ( w ):
            text = entry_cell_size.get_text ( ).strip ( ).replace ( ',', '.' )
            try:
                cell_size = float ( text )
                if cell_size <= 0.0:
                    raise ValueError ( "cell_size deve ser positivo" )
            except ValueError:
                label_decimate_status.set_text ( "Invalid cell size." )
                return
            traj = getattr ( vismol_object, "surface_trajectory", None )
            if not traj:
                label_decimate_status.set_text ( "No surface data to decimate." )
                return
            tris_before = sum (
                len ( frame_data[name][2] ) // 3
                for frame_data in traj for name in frame_data.keys ( )
            )
            try:
                vismol_object.surface_trajectory = mesh_decimation.decimate_surface_trajectory (
                    traj, cell_size
                )
            except Exception as e:
                label_decimate_status.set_text ( "Decimation failed: {}".format ( e ) )
                return
            tris_after = sum (
                len ( frame_data[name][2] ) // 3
                for frame_data in vismol_object.surface_trajectory for name in frame_data.keys ( )
            )
            label_decimate_status.set_text (
                "Triangles: {} -> {} ({:.0f}% reduction)".format (
                    tris_before, tris_after,
                    100.0 * ( 1.0 - tris_after / tris_before ) if tris_before else 0.0
                )
            )
            vm_glcore.queue_draw ( )

        btn_decimate.connect ( "clicked", on_decimate_clicked )
        hbox_decimate.pack_start ( entry_cell_size, False, False, 0 )
        hbox_decimate.pack_start ( btn_decimate, False, False, 0 )
        vbox.pack_start ( hbox_decimate, False, False, 0 )
        vbox.pack_start ( label_decimate_status, False, False, 0 )

        vbox.pack_start ( Gtk.Separator ( orientation = Gtk.Orientation.HORIZONTAL ), False, False, 4 )

        # --- suavizacao (Taubin -- ver vismol/utils/mesh_smoothing.pyx) ---
        # Suaviza o aspecto "escada" do marching cubes sem encolher a
        # superficie (ao contrario de um Laplacian puro). Aplicado direto
        # em vismol_object.surface_trajectory (todos os frames/lobulos),
        # substituindo a malha original -- assim como a decimacao, gerar a
        # superficie de novo desfaz.
        label_smooth = Gtk.Label ( label = "Smooth (Taubin) -- iterations:" )
        label_smooth.set_xalign ( 0 )
        vbox.pack_start ( label_smooth, False, False, 0 )

        hbox_smooth = Gtk.Box ( orientation = Gtk.Orientation.HORIZONTAL, spacing = 6 )
        entry_smooth_iters = Gtk.Entry ( )
        entry_smooth_iters.set_placeholder_text ( "e.g. 15" )
        entry_smooth_iters.set_width_chars ( 8 )
        btn_smooth = Gtk.Button ( label = "Apply" )
        label_smooth_status = Gtk.Label ( label = "" )
        label_smooth_status.set_xalign ( 0 )

        def on_smooth_clicked ( w ):
            text = entry_smooth_iters.get_text ( ).strip ( ).replace ( ',', '.' )
            try:
                n_iterations = int ( float ( text ) )
                if n_iterations <= 0:
                    raise ValueError ( "iterations deve ser positivo" )
            except ValueError:
                label_smooth_status.set_text ( "Invalid iteration count." )
                return
            traj = getattr ( vismol_object, "surface_trajectory", None )
            if not traj:
                label_smooth_status.set_text ( "No surface data to smooth." )
                return
            try:
                vismol_object.surface_trajectory = mesh_smoothing.smooth_surface_trajectory (
                    traj, iterations = n_iterations
                )
            except Exception as e:
                label_smooth_status.set_text ( "Smoothing failed: {}".format ( e ) )
                return
            label_smooth_status.set_text ( "Smoothed ({} iterations).".format ( n_iterations ) )
            vm_glcore.queue_draw ( )

        btn_smooth.connect ( "clicked", on_smooth_clicked )
        hbox_smooth.pack_start ( entry_smooth_iters, False, False, 0 )
        hbox_smooth.pack_start ( btn_smooth, False, False, 0 )
        vbox.pack_start ( hbox_smooth, False, False, 0 )
        vbox.pack_start ( label_smooth_status, False, False, 0 )

        # --- controles especificos do tipo ---
        if surface_type in ( "orbital", "density", "potential" ):
            vbox.pack_start ( Gtk.Separator ( orientation = Gtk.Orientation.HORIZONTAL ), False, False, 4 )

            def _current_lobe_rgba ( surf_name, fallback ):
                """ Le a cor ATUAL desse lobulo direto do cache (primeiro
                frame), pra pre-preencher o ColorButton -- fallback se o
                lobulo nao existir (ex: alguns tipos podem nao ter
                'obital_minus'). """
                try:
                    frame0 = vismol_object.surface_trajectory[0]
                    colors = frame0[surf_name][1]
                    return float ( colors[0] ), float ( colors[1] ), float ( colors[2] )
                except Exception:
                    return fallback

            hbox_colors = Gtk.Box ( orientation = Gtk.Orientation.HORIZONTAL, spacing = 10 )

            label_plus = Gtk.Label ( label = "Color (+):" )
            btn_color_plus = Gtk.ColorButton ( )
            r, g, b = _current_lobe_rgba ( "obital_plus", (1.0, 0.0, 0.0) )
            btn_color_plus.set_rgba ( Gdk.RGBA ( r, g, b, 1.0 ) )
            def on_color_plus_set ( w ):
                rgba = w.get_rgba ( )
                recolor_surface_lobe ( vismol_object, "obital_plus", (rgba.red, rgba.green, rgba.blue) )
                vm_glcore.queue_draw ( )
            btn_color_plus.connect ( "color-set", on_color_plus_set )

            label_minus = Gtk.Label ( label = "Color (-):" )
            btn_color_minus = Gtk.ColorButton ( )
            r, g, b = _current_lobe_rgba ( "obital_minus", (0.0, 0.0, 1.0) )
            btn_color_minus.set_rgba ( Gdk.RGBA ( r, g, b, 1.0 ) )
            def on_color_minus_set ( w ):
                rgba = w.get_rgba ( )
                recolor_surface_lobe ( vismol_object, "obital_minus", (rgba.red, rgba.green, rgba.blue) )
                vm_glcore.queue_draw ( )
            btn_color_minus.connect ( "color-set", on_color_minus_set )

            hbox_colors.pack_start ( label_plus, False, False, 0 )
            hbox_colors.pack_start ( btn_color_plus, False, False, 0 )
            hbox_colors.pack_start ( label_minus, False, False, 0 )
            hbox_colors.pack_start ( btn_color_minus, False, False, 0 )
            vbox.pack_start ( hbox_colors, False, False, 0 )

        elif surface_type == "mep":
            vbox.pack_start ( Gtk.Separator ( orientation = Gtk.Orientation.HORIZONTAL ), False, False, 4 )

            cmap_names = sorted ( COLOR_MAPS.keys ( ) )
            current_cmap = getattr ( vismol_object, "mep_cmap_name", "coolwarm" )
            current_vmin = getattr ( vismol_object, "mep_vmin", None )
            current_vmax = getattr ( vismol_object, "mep_vmax", None )

            label_cmap = Gtk.Label ( label = "Colormap:" )
            label_cmap.set_xalign ( 0 )
            cbx_cmap = Gtk.ComboBoxText ( )
            for cname in cmap_names:
                cbx_cmap.append ( cname, cname )
            cbx_cmap.set_active_id ( current_cmap if current_cmap in cmap_names else ( cmap_names[0] if cmap_names else None ) )

            label_vmin = Gtk.Label ( label = "vmin:" )
            label_vmin.set_xalign ( 0 )
            entry_vmin = Gtk.Entry ( )
            entry_vmin.set_placeholder_text ( "auto" )
            if current_vmin is not None:
                entry_vmin.set_text ( str ( current_vmin ) )

            label_vmax = Gtk.Label ( label = "vmax:" )
            label_vmax.set_xalign ( 0 )
            entry_vmax = Gtk.Entry ( )
            entry_vmax.set_placeholder_text ( "auto" )
            if current_vmax is not None:
                entry_vmax.set_text ( str ( current_vmax ) )

            label_status = Gtk.Label ( label = "" )
            label_status.set_xalign ( 0 )

            def _parse_optional_float ( entry ):
                text = entry.get_text ( ).strip ( )
                if text == "":
                    return None
                try:
                    return float ( text )
                except ValueError:
                    return None   # texto invalido -- cai no automatico (percentil)

            btn_apply = Gtk.Button ( label = "Apply" )
            def on_apply_clicked ( w ):
                vmin      = _parse_optional_float ( entry_vmin )
                vmax      = _parse_optional_float ( entry_vmax )
                cmap_name = cbx_cmap.get_active_id ( ) or "coolwarm"
                try:
                    recolor_mep_surface ( vismol_object, vmin = vmin, vmax = vmax, cmap_name = cmap_name )
                    vm_glcore.queue_draw ( )
                    label_status.set_text ( "Updated." )
                except ValueError as e:
                    label_status.set_text ( str ( e ) )
            btn_apply.connect ( "clicked", on_apply_clicked )

            vbox.pack_start ( label_cmap, False, False, 0 )
            vbox.pack_start ( cbx_cmap, False, False, 0 )
            vbox.pack_start ( label_vmin, False, False, 0 )
            vbox.pack_start ( entry_vmin, False, False, 0 )
            vbox.pack_start ( label_vmax, False, False, 0 )
            vbox.pack_start ( entry_vmax, False, False, 0 )
            vbox.pack_start ( btn_apply, False, False, 0 )
            vbox.pack_start ( label_status, False, False, 0 )

        window.show_all ( )
    
    def _menu_toggle_representation (self, menu_item, rep_type):
        """ [EN] Toggles representation `rep_type` on/off for the object
        that was right-clicked (self.vobject_index, set by open_menu()
        above). Reuses VismolObject.create_representation() (already
        wired up for every representation type shown in the
        'Representation' submenu, see vobject_menu_items above) --
        creating it fresh if it doesn't exist yet on this object (a new
        representation is active=True by construction), or just flipping
        its existing .active flag if it does, rather than tearing down
        and recreating an already-correctly-built one every time
        (particularly relevant for 'cartoon': recreating it re-runs the
        whole secondary-structure calculation and spline generation for
        no reason if it's already there and merely needs to be shown
        again). """
        vismol_object = self.treeview.main.vm_session.vm_objects_dic[self.vobject_index]
        rep = vismol_object.representations.get(rep_type)
        if rep is None:
            vismol_object.create_representation(rep_type=rep_type)
        else:
            rep.active = not rep.active
        self.treeview.main.vm_session.vm_glcore.queue_draw()

    def _menu_rename (self, menu_item = None ):
        """  
        menu_item = Gtk.MenuItem object at 0x7fbdcc035700 (GtkMenuItem at 0x37cf6c0)
        
        """
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()

        old_name = model.get_value(iter, 2)
        v_id     = model.get_value(iter, 1)
        e_id     = model.get_value(iter, 0)
        tag      = self.main.p_session.psystem[e_id].e_tag 
        
        old_name = old_name.split("- ")
        old_name = old_name[-1]
        
        if self.rename_window_visible:
            self.preferences.set_names (old_name, tag)
            pass
        else:
            
            self.preferences = PreferencesWindow(main = self.main, 
                                                 e_id = e_id     ,
                                                 v_id = v_id     )
            self.preferences.set_names (old_name, tag)
        self._save_backup_file()
        
    def destroy (self, widget):
        """ Function doc """
        self.rename_window_visible = False
    
    def _menu_delete_vm_object (self, widget):
        """ Function doc """
        self.main.delete_vm_object ( vm_object_index = self.vobject_index)
        self._save_backup_file()

    def _menu_view_log (self, widget):
        """ [EN] Opens the log file for whatever simulation CREATED the
        right-clicked object -- see the 'View Log' entry added to
        vobject_menu_items above for the full context.

        Reads vismol_object.results['logfile'] -- results is set once,
        directly on the object, by simulations_mixin._handle_result()
        (`vobject.results = results`) the moment a simulation produces a
        new vismol_object; results is the SAME dict
        _target_process()/run_simulation() build for the Process
        Manager's own job history (see process_manager_window.py's
        _open_log_for_row(), fixed for the exact same underlying gap in
        the same conversation this was added in: several simulation
        runner classes not correcting their OWN real log path -- already
        fixed at the source, in pdynamo/p_methods/*.py, so this reads
        the same, now-correct value).

        Objects with no .results at all (loaded from a file, drawn in
        the Builder, or any object predating this feature) get a plain,
        friendly message instead of a crash -- there's genuinely no log
        to show for those, not a bug to report. """
        vismol_object = self.treeview.main.vm_session.vm_objects_dic.get(self.vobject_index)
        if vismol_object is None:
            self.main.simple_dialog.error(msg="Could not find this object anymore.")
            return

        results = getattr(vismol_object, 'results', None)
        if not results:
            self.main.simple_dialog.error(
                msg="'{}' has no simulation log associated with it.\n\n"
                    "(Only objects CREATED by a simulation -- geometry "
                    "optimization, MD, a scan, ... -- have one.)".format(vismol_object.name))
            return

        logfile = results.get('logfile')
        if not logfile:
            self.main.simple_dialog.error(
                msg="'{}' has simulation data, but no log file path was recorded.".format(
                    vismol_object.name))
            return

        if not os.path.isfile(logfile):
            self.main.simple_dialog.error(
                msg="The log file for '{}' could not be found on disk "
                    "(it may have been moved or deleted):\n\n{}".format(vismol_object.name, logfile))
            return

        # ---- pDynamo log ----
        try:
            with open(logfile, 'r') as f:
                pdynamo_text = f.read()
        except Exception as exc:
            self.main.simple_dialog.error(
                msg="Could not read the log file:\n\n{}\n\nError: {}".format(logfile, exc))
            return

        # ---- QC program log (ORCA / xTB), if present ----
        # The QC log is copied to a permanent folder by the simulation
        # (backup_orca_files/backup_xtb_files in p_methods/_common.py), with the
        # SAME basename as the pDynamo log. We match by name + by the system's
        # actual QC model so we never show an ORCA log for an xTB job (or pick up
        # an unrelated log left in the folder by a previous run).
        system = None
        try:
            system = self.main.p_session.psystem[vismol_object.e_id]
        except Exception:
            system = None
        qc_label, qc_text = self._find_qc_log(logfile, system)

        tabs = [("pDynamo", pdynamo_text)]
        if qc_text:
            tabs.append((qc_label, qc_text))

        TabbedLogWindow(tabs, title="Log: {}".format(vismol_object.name))

    def _find_qc_log(self, pdynamo_logfile, system=None):
        """Locate and read the QC-program log that belongs to THIS job.

        Robust matching (avoids picking up an unrelated ORCA/xTB log left in the
        same folder by a previous run):
          1. by NAME: the QC log has the same basename as the pDynamo log, only
             swapping '.log' for '.orca.log'/'.xtb.log' (see backup_*_files).
          2. by QC MODEL: if the by-name match fails, scan the folder but only
             for the program the system's qcModel actually is -- never return an
             ORCA log for an xTB job.

        Returns (label, text); text is None when nothing suitable is found.
        """
        import glob
        folder = os.path.dirname(pdynamo_logfile)
        base = os.path.basename(pdynamo_logfile)
        if base.endswith(".log"):
            base = base[:-4]

        engine = None
        try:
            if system is not None and getattr(system, "qcModel", None):
                tag = system.qcModel.SummaryItems()[0][0]
                if "ORCA" in tag.upper():
                    engine = "ORCA"
                elif "XTB" in tag.upper():
                    engine = "XTB"
        except Exception:
            engine = None

        # 1) exact name match
        name_candidates = []
        if engine == "ORCA" or engine is None:
            name_candidates.append((os.path.join(folder, base + ".orca.log"), "ORCA"))
        if engine == "XTB" or engine is None:
            name_candidates.append((os.path.join(folder, base + ".xtb.log"), "xTB"))
        for path, label in name_candidates:
            if os.path.isfile(path):
                text = self._read_text(path)
                if text is not None:
                    return label, text

        # 2) fallback: scan folder, restricted to the engine this job used
        if folder and os.path.isdir(folder) and engine is not None:
            pattern = "*.orca.log" if engine == "ORCA" else "*.xtb.log"
            label = "ORCA" if engine == "ORCA" else "xTB"
            for path in sorted(glob.glob(os.path.join(folder, pattern))):
                if path == pdynamo_logfile:
                    continue
                text = self._read_text(path)
                if text is not None:
                    return label, text

        return "QC", None

    def _read_text(self, path):
        try:
            with open(path, 'r', errors='replace') as f:
                return f.read()
        except Exception:
            return None

    def _menu_delete_system (self, widget):
        """ Function doc """
        self.main.delete_system (system_e_id = self.system_e_id )
        self._save_backup_file()
        #self.save_easyhybrid_session( filename = self.main.session_filename, tmp = True)
    def build_tree_view_menu_old (self, menu_items = None):
        """ Function doc """
        tree_view_menu = Gtk.Menu()
        menu_header    = None
        
        for label in menu_items:
            
            
            if menu_items[label] == None:
                # just a label
                
                mitem = Gtk.MenuItem(label = label)
                if label == 'header':
                    menu_header    = mitem
                
                
            elif  label == '_separator':
                mitem = Gtk.SeparatorMenuItem()
            
            else:
                mitem = Gtk.MenuItem(label = label)
                mitem.connect('activate', menu_items[label])
            
            tree_view_menu.append(mitem)
            #mitem = Gtk.SeparatorMenuItem()
            #self.tree_view_menu.append(mitem)

        tree_view_menu.show_all()
        return tree_view_menu, menu_header

    def build_tree_view_menu (self, menu_items):
        """Cria menus e submenus a partir de um dicionário."""
        menu = Gtk.Menu()
        menu_header = None

        for label, value in menu_items.items():

            # --- Separador ---
            if label == "_separator":
                mitem = Gtk.SeparatorMenuItem()

            # --- Header (item destacado) ---
            elif label == "header":
                mitem = Gtk.MenuItem(label=label)
                mitem.set_sensitive(False)      # desabilita
                menu_header = mitem

            # --- Submenu (value é um dicionário) ---
            elif isinstance(value, dict):
                mitem = Gtk.MenuItem(label=label)
                # cria o submenu recursivamente
                submenu, _ = self.build_tree_view_menu(value)
                mitem.set_submenu(submenu)

            # --- Item simples sem callback ---
            elif value is None:
                mitem = Gtk.MenuItem(label=label)

            # --- Item com callback ---
            else:
                mitem = Gtk.MenuItem(label=label)
                mitem.connect("activate", value)

            menu.append(mitem)

        menu.show_all()
        return menu, menu_header

        def _build_submenus_from_dicts(self, menu_dict):
            """ Function doc """
            menu = Gtk.Menu()
            for key in menu_dict:
                mitem = Gtk.MenuItem(key)
                
                
                if menu_dict[key][0] == "submenu":
                    menu2 = self._build_submenus_from_dicts(menu_dict[key][1])
                    mitem.set_submenu(menu2)
                
                elif menu_dict[key][0] == "separator":
                    mitem = Gtk.SeparatorMenuItem()
                
                else:
                    if menu_dict[key][1] != None:
                        mitem.connect("activate", menu_dict[key][1])
                    else:
                        pass
                menu.append(mitem)
            return menu
        
        def _build_treemenu_from_dicts(self, menu_dict):
            """ Function doc """
            tree_view_menu = Gtk.Menu()
            menu_header    = None
            
            for key in menu_dict:
                mitem = Gtk.MenuItem(label=key)
                
                if menu_dict[key][0] == "submenu":
                    menu2 = self._build_submenus_from_dicts(menu_dict[key][1])
                    mitem.set_submenu(menu2)
                
                elif key == 'header':
                        menu_header    = mitem
                
                elif menu_dict[key][0] == "separator":
                    mitem = Gtk.SeparatorMenuItem()
                else:
                    if menu_dict[key][1] != None:
                        mitem.connect("activate", menu_dict[key][1])
                    else:
                        pass
                tree_view_menu.append(mitem)
            return tree_view_menu, menu_header


        def open_menu (self, system_e_id = None , vobject_index = None):
            """ Function doc """
            self.system_e_id     = system_e_id    
            self.vobject_index = vobject_index
            #print(system_e_id, vobject_index)
            
            system = self.treeview.main.p_session.psystem[self.system_e_id] 
            self.tree_header_sys_menu.set_label(system.label)
            
            if vobject_index == -1:
                
                self.tree_view_sys_menu.popup(None, None, None, None, 0, 0)

            if vobject_index != None and vobject_index != -1:
                
                vismol_object = self.treeview.main.vm_session.vm_objects_dic[vobject_index]
                self.tree_header_vobj_menu.set_label(vismol_object.name)
                
                self.tree_view_vobj_menu.popup(None, None, None, None, 0, 0)
                    
        def _save_backup_file (self):
            """ Function doc """
            # [ATUALIZACAO] Antes chamava save_easyhybrid_session(tmp=True)
            # direto, incondicional -- ignorava completamente o toggle
            # gl_parameters['autosave'] e nao contribuia pro criterio de
            # contador de eventos. Agora passa por register_change_and_
            # maybe_autosave, que respeita o toggle e so' salva de fato ao
            # atingir gl_parameters['autosave_event_count'] (ou via o timer
            # periodico em main_window.py).
            self.main.p_session.register_change_and_maybe_autosave()

    def open_menu (self, system_e_id = None , vobject_index = None):
        """ Function doc """
        self.system_e_id     = system_e_id    
        self.vobject_index = vobject_index
        dprint(system_e_id, vobject_index)
        
        system = self.treeview.main.p_session.psystem[self.system_e_id] 
        self.tree_header_sys_menu.set_label(system.label)
        
        if vobject_index == -1:
            
            self.tree_view_sys_menu.popup(None, None, None, None, 0, 0)

        
        if vobject_index != None and vobject_index != -1:
            
            vismol_object = self.treeview.main.vm_session.vm_objects_dic[vobject_index]
            self.tree_header_vobj_menu.set_label(vismol_object.name)
            is_surface  =  getattr(vismol_object, 'is_surface', None)
            #vismol object mgiht be a surface or a struture 
            if is_surface:
                dprint('is_surface:', vismol_object.is_surface)
                self.tree_header_surf_menu.set_label(vismol_object.name)
                self.tree_view_surf_menu.popup(None, None, None, None, 0, 0)
            else:
                self.tree_view_vobj_menu.popup(None, None, None, None, 0, 0)
                
    def _save_backup_file (self):
        """ Function doc """
        # Ver comentario na outra definicao de _save_backup_file, acima.
        self.main.p_session.register_change_and_maybe_autosave()
        
    def call_copy_current_frame_old_not_used (self, widget):
        """ Function doc """
        selection        = self.treeview.get_selection()
        (model, iter)    = selection.get_selected()
        e_id             = int(model.get_value(iter, 0))
        vm_object_index  = int(model.get_value(iter, 1))
        vobject = self.main.vm_session.vm_objects_dic[vm_object_index]
        
        frames  = vobject.frames
        frame_state = self.main.vm_session.get_frame()
        #print(vobject, frames ,frame_state)
        
        atom_qtty = len(vobject.atoms.items())
        size = len(vobject.frames)
        dprint(atom_qtty, size)
        
        if frame_state > size-1:
            init_frame = frames[-1]
        else:
            init_frame = frames[frame_state]


        new_traje  = np.empty([1, int(atom_qtty), 3], dtype=np.float32)
        for j, xyz in enumerate(init_frame):
            new_traje[0][j][0] = init_frame[j][0]
            new_traje[0][j][1] = init_frame[j][1]
            new_traje[0][j][2] = init_frame[j][2]

        #system = self.main_session.p_session.psystem[vobject.e_id]
        system = self.main.p_session.psystem[vobject.e_id]
        vobject = self.main.p_session._add_vismol_object_to_easyhybrid_session (system, show_molecule=True, name = 'edited_coords')
        vobject.frames = new_traje

        # Apply fixed representation to the VisMol object
        self.main.p_session._apply_fixed_representation_to_vobject(vismol_object =vobject)
        
        # Apply QC representation to the VisMol object
        self.main.p_session._apply_QC_representation_to_vobject(vismol_object =vobject)
        
        # Refresh the widgets in the main window
        self.main.main_treeview.refresh_number_of_frames()

    def call_delete_current_frame_old_not_used (self, widget):
        """ Function doc """
        dprint('delete_current_frame / Under construction')
        selection        = self.treeview.get_selection()
        (model, iter)    = selection.get_selected()
        e_id             = int(model.get_value(iter, 0))
        vm_object_index  = int(model.get_value(iter, 1))

        vobject = self.main.vm_session.vm_objects_dic[vm_object_index]
        frames  = vobject.frames
        
        frame_state = self.main.vm_session.get_frame()
        atom_qtty = len(vobject.atoms.items())
        size = len(vobject.frames)
        dprint(atom_qtty, size, type(frames))
        
        if frame_state > size-1:
            vobject.frames = np.delete(frames,-1, axis=0)
        else:
            vobject.frames = np.delete(frames, frame_state, axis=0)
        
        self.treeview.refresh_number_of_frames()
        self.treeview.refresh_trajectory_scalebar()
        self.main.vm_session.vm_glcore.queue_draw()
