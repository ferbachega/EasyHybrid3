#!/usr/bin/env python
# -*- coding: utf-8 -*-
#  
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
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
#      to facilitate QC/MM partitioning and molecular simulations.
#
from util.debug import dprint
import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Pango
from gi.repository import Gdk
from gi.repository import GdkPixbuf
import cairo


#---------------------------------------
from pBabel                    import*                                     
from pCore                     import*  
#---------------------------------------
from pMolecule                 import*                              
from pMolecule.MMModel         import*
from pMolecule.NBModel         import*                                     
from pMolecule.QCModel         import*
#---------------------------------------
from pScientific               import*                                     
from pScientific.Arrays        import*                                     
from pScientific.Geometry3     import*                                     
from pScientific.RandomNumbers import*                                     
from pScientific.Statistics    import*
from pScientific.Symmetry      import*
#---------------------------------------                              
from pSimulation               import*
#---------------------------------------

#import Pickle
from gui.widgets.custom_widgets import FolderChooserButton
from gui.widgets.custom_widgets import VismolTrajectoryFrame
from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.widgets.custom_widgets import get_colorful_square_pixel_buffer
from gui.widgets.custom_widgets import ReactionCoordinateBox

#from gui.widgets.custom_widgets import get_distance
from util.geometric_analysis            import get_distance 
from util.geometric_analysis            import get_dihedral 
from util.geometric_analysis            import get_angle 

from pdynamo.p_methods import LogFile


import util.orca_qc_keywords as orca_keys

from util.file_parser import get_file_type  
from util.file_parser import read_MOL2  
import pprint
import numpy as np
import gc
import os

import traceback


VISMOL_HOME      = os.environ.get('VISMOL_HOME')
HOME             = os.environ.get('HOME')
PDYNAMO3_SCRATCH = os.environ.get('PDYNAMO3_SCRATCH')


class SetupORCAWindow:
    """ Class doc """
    
    def __init__ (self, main, setup_QC_model_window):
        """ Class initialiser """
        self.main_session     = main
        self.home             = main.home
        self.Visible          = False
        self.vismol_object    = None
        self.window           = None

        self.setup_QC_model_window = setup_QC_model_window

        # . A generous but not unbounded CPU ceiling: os.cpu_count() can be
        #   None (containers/odd platforms), so fall back to something
        #   still sane rather than crashing on window build.
        ncpu_max = max ( 64, ( os.cpu_count ( ) or 8 ) * 2 )
        self.adjustment_nCPU = Gtk.Adjustment(value=1,
                                              lower=1,
                                              upper=ncpu_max,
                                              step_increment=1,
                                              page_increment=1,
                                              page_size=0)


        self.orca_methods_dict = {}
        self.restrited_label = ''

        # . Persisted across window close/reopen (unlike the comboboxes
        #   themselves, which are rebuilt from scratch every open_window()):
        #   without this, closing the window to tweak one setting used to
        #   reset method/basis/SCF/CPUs/scratch/keywords back to defaults.
        #   'extra_lines' holds whatever the user typed in textview_orca
        #   PAST line 1 (block directives, extra "!" lines, ...) -- line 1
        #   itself is always regenerated from the fields and never stored
        #   here.
        self.parameters = {
            'type_id'        : 0    ,
            'method_id'      : 0    ,
            'excited_id'     : 0    ,
            'scf_id'         : 0    ,
            'basis_group_id' : 0    ,
            'basis_id'       : 0    ,
            'ncpu'           : 1    ,
            'scratch'        : None ,
            'random_scratch' : False,
            'extra_lines'    : ''   ,
        }
        
        #---------------------------------------------------------------------
        # ORCA
        # this liststore contains the types of methods (HF, local and gradient corrections, ...)
        self.liststore_orca_type_of_method = Gtk.ListStore(str)
        for key, _type in orca_keys.type_of_method.items():
            self.liststore_orca_type_of_method.append([_type])

            # this liststore contains the methods (HF, B1LYP, B3LYP ...)
            self.orca_methods_dict[key] = Gtk.ListStore(str)
            
            for method in orca_keys.methods[key]:
                self.orca_methods_dict[key].append([method])
        #---------------------------------------------------------------------
        
        
        #---------------------------------------------------------------------
        #excited states
        self.excited_states_liststore_dict = {}
        self.excited_states_liststore_dict['semiemp'] = Gtk.ListStore(str)
        for excited_state in orca_keys.list_semiemp_excited.values():
            self.excited_states_liststore_dict['semiemp'].append([excited_state])
        
        self.excited_states_liststore_dict['HF'] = Gtk.ListStore(str)
        for excited_state in orca_keys.list_HF_excited.values():
            self.excited_states_liststore_dict['HF'].append([excited_state])
        
        
        self.excited_states_liststore_dict['DFT'] = Gtk.ListStore(str)
        for excited_state in orca_keys.list_DFT_excited.values():
            self.excited_states_liststore_dict['DFT'].append([excited_state])
        
        #---------------------------------------------------------------------
        
        
        
        #---------------------------------------------------------------------
        self.liststore_orca_scf_convergence = Gtk.ListStore(str)
        for scf_convergence in orca_keys.list_scf_convergence.values():
            self.liststore_orca_scf_convergence.append([scf_convergence])
        #---------------------------------------------------------------------
        
        
        #---------------------------------------------------------------------
        #basis set and basis set group
        self.orca_type_basis_dict = {}                       # contain  several liststores
        self.liststore_orca_type_basis = Gtk.ListStore(str)  # 
        
        for key , type_basis in  orca_keys.basis_set_group_dict.items():

            self.liststore_orca_type_basis.append([type_basis])
            self.orca_type_basis_dict[key] =  Gtk.ListStore(str)
            for basis in orca_keys.type_basis_dict[key]:
                self.orca_type_basis_dict[key].append([basis])   
        #---------------------------------------------------------------------
        self.comboboxes = {}



    def open_window (self, vismol_object = None):
        """ Function doc """
        if self.Visible and self.window is not None:
            # . Already open: bring it to the front instead of silently
            #   doing nothing (the old behaviour if the window ended up
            #   behind another one).
            self.window.present ( )
            return

        if self.Visible  ==  False:
            self.builder = Gtk.Builder()
            self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/setup/setup_qc_model_window.glade'))
            self.builder.connect_signals(self)
            
            self.window = self.builder.get_object('window_setup_orca')
            self.window.set_title('ORCA Setup')
            self.window.set_keep_above(True)

            #-------------------------------------------------------------------------------
            self.orca_type_of_method_combo = self.builder.get_object('combobox_orca_type')
            self.orca_type_of_method_combo.connect("changed", self.on_combobox_orca_type_changed)
            self.orca_type_of_method_combo.set_model(self.liststore_orca_type_of_method)
            renderer_text = Gtk.CellRendererText()
            self.orca_type_of_method_combo.pack_start(renderer_text, True)
            self.orca_type_of_method_combo.add_attribute(renderer_text, "text", 0)
            
            #-------------------------------------------------------------------------------
            
            #-------------------------------------------------------------------------------
            self.orca_methods_combo = self.builder.get_object('combobox_orca_method')
            self.orca_methods_combo.connect("changed", self.on_combobox_orca_method_changed)
            self.orca_methods_combo.set_model( self.orca_methods_dict[0] )
            renderer_text = Gtk.CellRendererText()
            self.orca_methods_combo.pack_start(renderer_text, True)
            self.orca_methods_combo.add_attribute(renderer_text, "text", 0)
            #
            #-------------------------------------------------------------------------------
            
            
            #-------------------------------------------------------------------------------
            self.combobox_orca_excited_states = self.builder.get_object('combobox_orca_excited_states')
            #self.combobox_orca_excited_states.connect("changed", self.on_combobox_orca_method_changed)
            self.combobox_orca_excited_states.set_model( self.excited_states_liststore_dict['HF'] )
            renderer_text = Gtk.CellRendererText()
            self.combobox_orca_excited_states.pack_start(renderer_text, True)
            self.combobox_orca_excited_states.add_attribute(renderer_text, "text", 0)
            #
            #-------------------------------------------------------------------------------
            
            
            #-------------------------------------------------------------------------------
            self.orca_scf_conv_combo = self.builder.get_object('combobox_orca_scf_convergence')
            #self.orca_scf_conv_combo.connect("changed", self.on_combobox_orca_method_changed)
            self.orca_scf_conv_combo.set_model( self.liststore_orca_scf_convergence )
            renderer_text = Gtk.CellRendererText()
            self.orca_scf_conv_combo.pack_start(renderer_text, True)
            self.orca_scf_conv_combo.add_attribute(renderer_text, "text", 0)
            #
            #-------------------------------------------------------------------------------
            
            
            
            #-------------------------------------------------------------------------------
            self.orca_basis_group_combo = self.builder.get_object('combobox_orca_basis_set_group')
            self.orca_basis_group_combo.connect("changed", self.on_orca_basis_group_combo_changed)
            self.orca_basis_group_combo.set_model( self.liststore_orca_type_basis )
            renderer_text = Gtk.CellRendererText()
            self.orca_basis_group_combo.pack_start(renderer_text, True)
            self.orca_basis_group_combo.add_attribute(renderer_text, "text", 0)
            #-------------------------------------------------------------------------------
            
            
            #-------------------------------------------------------------------------------
            self.orca_basis_set_combo = self.builder.get_object('combobox_orca_basis_set')
            #self.orca_basis_set_combo.connect("changed", self.on_orca_basis_set_combo_changed)
            self.orca_basis_set_combo.set_model( self.orca_type_basis_dict[0] )
            renderer_text = Gtk.CellRendererText()
            self.orca_basis_set_combo.pack_start(renderer_text, True)
            self.orca_basis_set_combo.add_attribute(renderer_text, "text", 0)
            #-------------------------------------------------------------------------------            
            
            
            #-------------------------------------------------------------------------------
            self.spinbutton_orca_ncpus = self.builder.get_object('spin_button_orca_ncpus')
            self.spinbutton_orca_ncpus.set_adjustment(self.adjustment_nCPU)
            #-------------------------------------------------------------------------------


            #-------------------------------------------------------------------------------
            self.entry_orca_scratch   = self.builder.get_object('entry_orca_scratch_folder')
            self.button_orca_browse   = self.builder.get_object('button_orca_scratch_browse')
            self.button_orca_browse.connect('clicked', self.on_button_browse_scratch)

            if self.parameters['scratch']:
                # . A previous OK already set this: restore it verbatim
                #   rather than recomputing the per-system default below.
                self.entry_orca_scratch.set_text(self.parameters['scratch'])
            else:
                self.entry_orca_scratch.set_text(self._default_scratch_folder())
            #-------------------------------------------------------------------------------


            #-------------------------------------------------------------------------------
            # . "ORCA: /path/to/orca" or a warning if PDYNAMO3_ORCACOMMAND is
            #   unset/invalid -- found out BEFORE the user builds a whole
            #   calculation around it, not deep in a NotInstalledError once
            #   the job actually runs.
            self.label_orca_engine_status = self.builder.get_object('label_orca_engine_status')
            self.refresh_orca_engine_status()
            #-------------------------------------------------------------------------------


            '''-----------------------------------------------------------------------------------'''

            self.combobox_list = [self.orca_type_of_method_combo,
                                  self.orca_methods_combo,
                                  self.combobox_orca_excited_states,
                                  self.orca_scf_conv_combo,
                                  self.orca_basis_group_combo,
                                  self.orca_basis_set_combo,
                                 ]

            self.window.show_all()

            # . Restore the persisted selections (defaults to all-zero on the
            #   first open). Order matters: setting type/basis-group fires
            #   on_combobox_orca_type_changed()/on_orca_basis_group_combo_changed(),
            #   which reset method/excited-states/basis-set back to index 0
            #   -- so those three are restored explicitly AFTER, overriding
            #   that reset.
            self.orca_type_of_method_combo.set_active(self.parameters['type_id'])
            self.orca_methods_combo.set_active(self.parameters['method_id'])
            self.combobox_orca_excited_states.set_active(self.parameters['excited_id'])
            self.orca_scf_conv_combo.set_active(self.parameters['scf_id'])
            self.orca_basis_group_combo.set_active(self.parameters['basis_group_id'])
            self.orca_basis_set_combo.set_active(self.parameters['basis_id'])
            self.spinbutton_orca_ncpus.set_value(self.parameters['ncpu'])


            self.spinbutton_orca_ncpus.connect('value-changed', self.refresh_orca_parameters)

            self.orca_type_of_method_combo.connect('changed', self.refresh_orca_parameters)
            self.orca_methods_combo.connect('changed', self.refresh_orca_parameters)
            self.combobox_orca_excited_states.connect('changed', self.refresh_orca_parameters)
            self.orca_scf_conv_combo.connect('changed', self.refresh_orca_parameters)
            self.orca_basis_group_combo.connect('changed', self.refresh_orca_parameters)
            self.orca_basis_set_combo.connect('changed', self.refresh_orca_parameters)


            self.checkbox_orca_random_scratch = self.builder.get_object('checkbox_orca_random_scratch')
            self.checkbox_orca_random_scratch.set_active(self.parameters['random_scratch'])

            self.button_ok         = self.builder.get_object('orca_button_ok')
            self.button_cancel     = self.builder.get_object('orca_button_cancel')

            self.button_ok.connect("clicked", self.on_button_ok)
            self.button_cancel.connect("clicked", self.close_window)

            self.window.connect("destroy", self.close_window)
            # . refresh_orca_parameters() restores self.parameters['extra_lines']
            #   (if any) past the EXTRA_MARKER the first time it runs, since
            #   the buffer has no marker yet at this point.
            self.refresh_orca_parameters (None)

            # . Watch for the user's own edits past the marker, so the
            #   "doesn't start with ! or %" warning below the keywords box
            #   stays live as they type.
            textbuffer = self.builder.get_object('textview_orca').get_buffer()
            textbuffer.connect('changed', self.on_textview_orca_changed)
            self.on_textview_orca_changed(textbuffer)

            self.Visible  =  True
            
    def on_orca_basis_group_combo_changed (self, widget):
        """ Function doc """
        _id = widget.get_active()
        self.orca_basis_set_combo.set_model( self.orca_type_basis_dict[_id])
        self.orca_basis_set_combo.set_active( 0 )        
    
    def on_combobox_orca_type_changed (self, widget):
        """ Function doc """
        _id = widget.get_active()
        self.orca_methods_combo.set_model( self.orca_methods_dict[_id] )
        self.orca_methods_combo.set_active( 0 )

        self.restrited = self.setup_QC_model_window.builder.get_object('radio_button_restricted').get_active()

        if _id in [0]:
            excited_model = self.excited_states_liststore_dict['HF']
            self.restrited_label = 'RHF' if self.restrited else 'UHF'

        elif _id in [1,2,3,4]:
            excited_model = self.excited_states_liststore_dict['DFT']
            self.restrited_label = 'RKS' if self.restrited else 'UKS'

        elif _id == 7:
            excited_model = self.excited_states_liststore_dict['semiemp']
            self.restrited_label = 'RHF' if self.restrited else 'UHF'

        else:
            # . No excited-states keyword list makes sense for this method
            #   type (e.g. force fields/semi-empirical variants not listed
            #   above): leave the combo empty rather than giving it a model
            #   from an unrelated type.
            excited_model = None
            self.restrited_label = 'RHF' if self.restrited else 'UHF'

        self.combobox_orca_excited_states.set_model( excited_model )
        if excited_model is not None:
            # . set_active() on a combo with no model logs a GLib warning
            #   for no benefit -- only call it when there is something to
            #   select.
            self.combobox_orca_excited_states.set_active(0)

    def on_combobox_orca_method_changed (self, widget):
        """ Function doc """

    def close_window (self, button, data  = None):
        """ Function doc """
        if self.window:
            self.window.destroy()
        self.Visible    =  False

    def _default_scratch_folder (self):
        """ PDYNAMO3_SCRATCH/<system tag>, mirroring SetupXTBWindow's own
        per-system scratch (so two systems open at once don't share ORCA
        job files); falls back to $HOME if PDYNAMO3_SCRATCH is unset. """
        base = os.environ.get('PDYNAMO3_SCRATCH') or os.environ.get('HOME') or ''
        try:
            _id = self.main_session.p_session.active_id
            system = self.main_session.p_session.psystem[_id]
            tag = getattr(system, 'e_tag', None)
        except Exception:
            tag = None
        return os.path.join(base, tag) if tag else base

    def refresh_orca_engine_status (self):
        """ Updates label_orca_engine_status with the configured
        PDYNAMO3_ORCACOMMAND, or a warning if it is unset/invalid -- found
        out here, before the user builds a whole calculation around it,
        rather than much later as a bare NotInstalledError once a job
        actually runs. """
        command = os.environ.get('PDYNAMO3_ORCACOMMAND')
        if command and os.path.isfile(command) and os.access(command, os.X_OK):
            self.label_orca_engine_status.set_text('ORCA executable: {}'.format(command))
        else:
            self.label_orca_engine_status.set_text(
                '⚠ ORCA executable not found (PDYNAMO3_ORCACOMMAND). Calculations with this '
                'setup will fail until it is configured -- see Preferences / install.py.')

    def on_button_browse_scratch (self, button):
        """ Function doc """
        dialog = Gtk.FileChooserDialog(
            title='Choose the ORCA scratch folder', parent=self.window,
            action=Gtk.FileChooserAction.SELECT_FOLDER)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                           Gtk.STOCK_OPEN,   Gtk.ResponseType.OK)
        current = self.entry_orca_scratch.get_text()
        if current and os.path.isdir(current):
            dialog.set_current_folder(current)
        if dialog.run() == Gtk.ResponseType.OK:
            self.entry_orca_scratch.set_text(dialog.get_filename())
        dialog.destroy()

    def on_button_ok (self, button):
        """ Function doc """
        textbuffer = self.builder.get_object('textview_orca').get_buffer ()
        text = textbuffer.get_text (textbuffer.get_start_iter(), textbuffer.get_end_iter(), True)

        scratch = self.entry_orca_scratch.get_text().strip()
        if not scratch:
            scratch = self._default_scratch_folder()
            self.entry_orca_scratch.set_text(scratch)

        if not os.path.isdir(scratch):
            try:
                os.makedirs(scratch, exist_ok=True)
            except Exception as e:
                self._show_message(Gtk.MessageType.ERROR,
                    'Could not create the scratch folder "{}": {}\n\nFix the path and try again.'.format(scratch, e))
                return  # . keep the window open -- OK must not "succeed" with a broken scratch path
        if not os.access(scratch, os.W_OK):
            self._show_message(Gtk.MessageType.ERROR,
                'The scratch folder "{}" is not writable.\n\nFix its permissions (or pick another folder) and try again.'.format(scratch))
            return

        ncpu = self.spinbutton_orca_ncpus.get_value_as_int()
        cpu_count = os.cpu_count()
        if cpu_count and ncpu > cpu_count:
            proceed = self._ask_yes_no(
                'Number of CPUs ({}) is higher than this machine\'s {} available cores. '
                'ORCA may run slower than with fewer, oversubscribed processes.\n\nUse it anyway?'.format(ncpu, cpu_count))
            if not proceed:
                return  # . let them lower it instead of silently capping or guessing

        command = os.environ.get('PDYNAMO3_ORCACOMMAND')
        if not (command and os.path.isfile(command) and os.access(command, os.X_OK)):
            # . A warning, not a block: this setup may be meant for a
            #   remote machine where ORCA is configured, even if it is not
            #   here.
            self._show_message(Gtk.MessageType.WARNING,
                'No valid ORCA executable is configured (PDYNAMO3_ORCACOMMAND). '
                'This setup will be saved, but running it will fail until that is fixed.')

        # . Persist every selection so reopening this window restores it
        #   instead of resetting to defaults (comboboxes are rebuilt from
        #   scratch on every open_window(), self.parameters is not).
        self.parameters['type_id']        = self.orca_type_of_method_combo.get_active()
        self.parameters['method_id']      = self.orca_methods_combo.get_active()
        self.parameters['excited_id']     = self.combobox_orca_excited_states.get_active()
        self.parameters['scf_id']         = self.orca_scf_conv_combo.get_active()
        self.parameters['basis_group_id'] = self.orca_basis_group_combo.get_active()
        self.parameters['basis_id']       = self.orca_basis_set_combo.get_active()
        self.parameters['ncpu']           = ncpu
        self.parameters['scratch']        = scratch
        self.parameters['random_scratch'] = self.checkbox_orca_random_scratch.get_active()
        self.parameters['extra_lines']    = (
            text.split(self.EXTRA_MARKER, 1)[1].lstrip('\n') if self.EXTRA_MARKER in text else '')

        self.setup_QC_model_window.orca_options        = text
        self.setup_QC_model_window.orca_scratch        = scratch
        self.setup_QC_model_window.orca_random_scratch = self.parameters['random_scratch']

        dprint (text)
        dprint ('checkbox_random_scratch', self.parameters['random_scratch'])
        self.close_window (None)

    def _show_message (self, message_type, text):
        """ A simple modal message dialog -- validation feedback for
        on_button_ok() (scratch folder, CPU count, ...). """
        dialog = Gtk.MessageDialog(transient_for=self.window, flags=0, message_type=message_type,
                                   buttons=Gtk.ButtonsType.OK, text=text)
        dialog.run()
        dialog.destroy()

    def _ask_yes_no (self, text):
        """ A simple modal Yes/No dialog; returns True on Yes. """
        dialog = Gtk.MessageDialog(transient_for=self.window, flags=0, message_type=Gtk.MessageType.QUESTION,
                                   buttons=Gtk.ButtonsType.YES_NO, text=text)
        response = dialog.run()
        dialog.destroy()
        return response == Gtk.ResponseType.YES

    # . Marks where pDynamo/ORCA-generated keywords end and the user's own
    #   free-form extra input begins (refresh_orca_parameters() regenerates
    #   everything before this line on every field change; everything from
    #   this line on is left untouched). "#" is ORCA's own comment
    #   character, so the marker is itself valid, harmless ORCA input even
    #   if left in place.
    EXTRA_MARKER = '# --- extra ORCA input below (kept as-is; use "!" for a keyword line, "%...end" for a block) ---'

    def on_textview_orca_changed (self, textbuffer):
        """ Live check of the user's own extra input (below EXTRA_MARKER):
        warns (via the preview label) about a non-blank line that starts
        with neither "!" (another simple-keyword line), "%" (a block
        directive) nor "#" (a comment) nor is a bare "end" (closing a
        block) -- ORCA will likely not parse it as input. Advisory only,
        never blocks typing or OK. """
        full_text = textbuffer.get_text(textbuffer.get_start_iter(), textbuffer.get_end_iter(), True)
        extra = full_text.split(self.EXTRA_MARKER, 1)[1] if self.EXTRA_MARKER in full_text else ''
        bad = [ line.strip() for line in extra.splitlines()
                if line.strip() and not line.lstrip().startswith(('!', '%', '#'))
                and line.strip().lower() != 'end' ]
        self._orca_extra_warning = (
            '⚠ line(s) below the marker don\'t start with "!", "%" or "#": {}{} -- ORCA may not read '
            'them as input.'.format(', '.join(repr(l) for l in bad[:3]), ', ...' if len(bad) > 3 else '')
        ) if bad else ''
        self._render_orca_preview()

    def _render_orca_preview (self):
        """ Combines the cached "line 1 pDynamo will actually send" preview
        with the live extra-input warning (if any) into label_orca_preview. """
        label = self.builder.get_object('label_orca_preview')
        text = ('Sent to ORCA as line 1: ! ENERGY BOHRS {}  '
                '(ENGRAD instead of ENERGY when gradients are needed)').format(
                    getattr(self, '_orca_preview_first_line', ''))
        if getattr(self, '_orca_extra_warning', ''):
            text += '\n' + self._orca_extra_warning
        label.set_text(text)

    def refresh_orca_parameters (self, widget):
        """ Function doc """
        _iter = self.orca_methods_combo.get_active_iter()
        model = self.orca_methods_combo.get_model()
        method = model[_iter][0].split()
        method = method[0]


        _iter =self.orca_basis_set_combo.get_active_iter()
        model =self.orca_basis_set_combo.get_model()
        basis = model[_iter][0].split()
        basis = basis[0]


        _iter =self.orca_scf_conv_combo.get_active_iter()
        model =self.orca_scf_conv_combo.get_model()
        scf_conv = model[_iter][0].split()
        scf_conv = scf_conv[0]
        if scf_conv == 'Default':
            scf_conv = ''


        _iter =self.combobox_orca_excited_states.get_active_iter()
        model =self.combobox_orca_excited_states.get_model()
        e_states = model[_iter][0].split()
        e_states = e_states[0]

        if e_states == 'TD-DFT':
            e_states = '''
%tddft
     nroots 8 # the number of excited states to be calculated.
     maxdim 30 # the maximum dimension of the expansion space in the Davidson procedure.
 end #tddft'''
        else:
            e_states = ''


        #-----------------------------------------------------------------------
        ncpu = self.spinbutton_orca_ncpus.get_value_as_int()
        if ncpu > 1:
            ncpu = '\n%pal nprocs {} \n end'.format(ncpu)
        else:
            ncpu = ''

        #-----------------------------------------------------------------------

        new_text = '{} {} {} {} {} {}'.format(self.restrited_label,  method, basis, scf_conv, e_states, ncpu )

        # . Preserve whatever the user has typed past EXTRA_MARKER instead
        #   of wiping it out on every field change (the original behaviour:
        #   textbuffer.set_text(new_text), unconditionally, destroying any
        #   extra line the moment ANY combo changed). The marker is kept in
        #   the buffer EVEN WHEN NOTHING FOLLOWS IT YET -- it must already
        #   be there the first time the user types below the generated
        #   line, not only once something worth preserving exists, or that
        #   very first addition would still be lost on the next refresh.
        textbuffer = self.builder.get_object('textview_orca').get_buffer ()
        current_text = textbuffer.get_text (textbuffer.get_start_iter(), textbuffer.get_end_iter(), True)
        if self.EXTRA_MARKER in current_text:
            extra = current_text.split(self.EXTRA_MARKER, 1)[1]
        elif self.parameters.get('extra_lines'):
            # . First regeneration after open_window() restored a previous
            #   session's extra lines: the marker is not in the (still
            #   blank) buffer yet, so pull them from self.parameters instead.
            extra = '\n' + self.parameters['extra_lines']
        else:
            extra = ''

        full_text = new_text.rstrip() + '\n\n' + self.EXTRA_MARKER + extra
        textbuffer.set_text(full_text)

        self._orca_preview_first_line = new_text.split('\n', 1)[0].strip()
        self._render_orca_preview()
