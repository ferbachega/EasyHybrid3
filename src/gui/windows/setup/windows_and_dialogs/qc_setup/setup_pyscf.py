#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: System > Add > Force Field/QC Model > PySCF Setup
#
#  Copyright 2022-2026 Fernando Bachega
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
#      Setup window for the PySCF QC model (pdynamo3-eh's addOns/pySCF/
#      QCModelPySCF.py): method (HF/DFT, combined with the main window's
#      Restricted/Unrestricted choice), functional, basis set, ECP, OMP
#      thread count and scratch folder. Modelled on SetupORCAWindow's
#      pattern (self.parameters persisted across window close/reopen, an
#      engine status line, scratch auto-create + writability check, a
#      CPU-count warning) -- see that file's own history for why each of
#      those exists. The curated functional/basis dropdowns
#      (util/pyscf_qc_keywords.py) cover common choices; "Custom
#      Functional"/"Custom Basis Set" let anyone type any PySCF/libxc-valid
#      string beyond that list (checked, when pyscf is importable, against
#      pyscf.dft.libxc.parse_xc()/pyscf.gto.basis.load() before accepting).
#
#      QCModelPySCF quirk (pdynamo3-eh/addOns/pySCF/QCModelPySCF.py
#      _CheckOptions): the scratch folder is NOT a WithOptions kwarg --
#      QCModelPySCF reads it from the PYSCF_TMPDIR environment variable
#      (falling back to PDYNAMO3_SCRATCH) every time a model is built,
#      overwriting whatever was passed in. So the scratch folder this
#      window collects must be applied by setting os.environ['PYSCF_TMPDIR']
#      before QCModelPySCF.WithOptions() runs (session.py's PySCF branch),
#      not passed as a parameter the model itself understands. Same story
#      for OMP_NUM_THREADS (PySCF's own parallelism, read from the
#      environment, not a constructor option).
#
from util.debug import dprint
import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

import util.pyscf_qc_keywords as pyscf_keys

import os


class SetupPySCFWindow:
    """ Class doc """

    def __init__ (self, main, setup_QC_model_window):
        """ Class initialiser """
        self.main_session     = main
        self.home             = main.home
        self.Visible          = False
        self.vismol_object    = None
        self.window           = None

        self.setup_QC_model_window = setup_QC_model_window

        ncpu_max = max ( 64, ( os.cpu_count ( ) or 8 ) * 2 )
        self.adjustment_nCPU = Gtk.Adjustment(value=1,
                                              lower=1,
                                              upper=ncpu_max,
                                              step_increment=1,
                                              page_increment=1,
                                              page_size=0)

        self.restrited_label = ''

        # . Persisted across window close/reopen (the comboboxes
        #   themselves are rebuilt from scratch every open_window()).
        # . 'custom_functional'/'custom_basis' override the dropdown
        #   selections when non-empty -- PySCF/libxc accept far more
        #   functional and basis strings than any curated dropdown could
        #   list (see util/pyscf_qc_keywords.py's own module docstring).
        # . 'ecp_id' == 0 means "None (all-electron)", PySCF's own default.
        self.parameters = {
            'type_id'           : 0    ,
            'functional_id'     : 0    ,
            'basis_group_id'    : 0    ,
            'basis_id'          : 0    ,
            'ecp_id'            : 0    ,
            'custom_functional' : ''   ,
            'custom_basis'      : ''   ,
            'ncpu'              : 1    ,
            'scratch'           : None ,
            'delete_job_files'  : True ,
        }

        #---------------------------------------------------------------------
        self.liststore_pyscf_type_of_method = Gtk.ListStore(str)
        for _type in pyscf_keys.type_of_method.values():
            self.liststore_pyscf_type_of_method.append([_type])

        self.liststore_pyscf_functionals = Gtk.ListStore(str)
        for functional in pyscf_keys.functionals.values():
            self.liststore_pyscf_functionals.append([functional])

        self.liststore_pyscf_basis_group = Gtk.ListStore(str)
        self.pyscf_basis_dict = {}
        for key, group_name in pyscf_keys.basis_set_group_dict.items():
            self.liststore_pyscf_basis_group.append([group_name])
            self.pyscf_basis_dict[key] = Gtk.ListStore(str)
            for basis in pyscf_keys.type_basis_dict[key]:
                self.pyscf_basis_dict[key].append([basis])

        self.liststore_pyscf_ecp = Gtk.ListStore(str)
        for ecp in pyscf_keys.ecp_dict.values():
            self.liststore_pyscf_ecp.append([ecp])
        #---------------------------------------------------------------------

    def open_window (self, vismol_object = None):
        """ Function doc """
        if self.Visible and self.window is not None:
            self.window.present ( )
            return

        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/setup/setup_qc_model_window.glade'))

        self.window = self.builder.get_object('window_setup_pyscf')
        self.window.set_title('PySCF Setup')
        self.window.set_keep_above(True)

        #-------------------------------------------------------------------------------
        self.pyscf_type_combo = self.builder.get_object('combobox_pyscf_type')
        self.pyscf_type_combo.connect("changed", self.on_combobox_pyscf_type_changed)
        self.pyscf_type_combo.set_model(self.liststore_pyscf_type_of_method)
        renderer_text = Gtk.CellRendererText()
        self.pyscf_type_combo.pack_start(renderer_text, True)
        self.pyscf_type_combo.add_attribute(renderer_text, "text", 0)

        self.pyscf_functional_combo = self.builder.get_object('combobox_pyscf_functional')
        self.pyscf_functional_combo.set_model(self.liststore_pyscf_functionals)
        renderer_text = Gtk.CellRendererText()
        self.pyscf_functional_combo.pack_start(renderer_text, True)
        self.pyscf_functional_combo.add_attribute(renderer_text, "text", 0)

        self.pyscf_basis_group_combo = self.builder.get_object('combobox_pyscf_basis_group')
        self.pyscf_basis_group_combo.connect("changed", self.on_pyscf_basis_group_combo_changed)
        self.pyscf_basis_group_combo.set_model(self.liststore_pyscf_basis_group)
        renderer_text = Gtk.CellRendererText()
        self.pyscf_basis_group_combo.pack_start(renderer_text, True)
        self.pyscf_basis_group_combo.add_attribute(renderer_text, "text", 0)

        self.pyscf_basis_combo = self.builder.get_object('combobox_pyscf_basis')
        self.pyscf_basis_combo.set_model(self.pyscf_basis_dict[0])
        renderer_text = Gtk.CellRendererText()
        self.pyscf_basis_combo.pack_start(renderer_text, True)
        self.pyscf_basis_combo.add_attribute(renderer_text, "text", 0)

        self.pyscf_ecp_combo = self.builder.get_object('combobox_pyscf_ecp')
        self.pyscf_ecp_combo.set_model(self.liststore_pyscf_ecp)
        renderer_text = Gtk.CellRendererText()
        self.pyscf_ecp_combo.pack_start(renderer_text, True)
        self.pyscf_ecp_combo.add_attribute(renderer_text, "text", 0)
        #-------------------------------------------------------------------------------

        self.entry_pyscf_custom_functional = self.builder.get_object('entry_pyscf_custom_functional')
        self.entry_pyscf_custom_functional.set_text(self.parameters['custom_functional'])
        self.entry_pyscf_custom_basis      = self.builder.get_object('entry_pyscf_custom_basis')
        self.entry_pyscf_custom_basis.set_text(self.parameters['custom_basis'])

        self.spinbutton_pyscf_ncpus = self.builder.get_object('spin_button_pyscf_ncpus')
        self.spinbutton_pyscf_ncpus.set_adjustment(self.adjustment_nCPU)

        #-------------------------------------------------------------------------------
        self.entry_pyscf_scratch = self.builder.get_object('entry_pyscf_scratch_folder')
        self.button_pyscf_browse = self.builder.get_object('button_pyscf_scratch_browse')
        self.button_pyscf_browse.connect('clicked', self.on_button_browse_scratch)

        if self.parameters['scratch']:
            self.entry_pyscf_scratch.set_text(self.parameters['scratch'])
        else:
            self.entry_pyscf_scratch.set_text(self._default_scratch_folder())
        #-------------------------------------------------------------------------------

        #-------------------------------------------------------------------------------
        # . "pyscf x.y.z" or a warning if the package is not importable --
        #   found out here, before the user builds a calculation around
        #   it, rather than as a bare NotInstalledError once a job runs.
        self.label_pyscf_engine_status = self.builder.get_object('label_pyscf_engine_status')
        self.refresh_pyscf_engine_status()
        #-------------------------------------------------------------------------------

        self.checkbox_pyscf_delete_job_files = self.builder.get_object('checkbox_pyscf_delete_job_files')

        self.window.show_all()

        # . Restore the persisted selections (all-zero/defaults on the
        #   first open). Setting type fires on_combobox_pyscf_type_changed(),
        #   which does not reset the functional combo's index (unlike
        #   ORCA's cascading type->method->excited-states resets) -- both
        #   restores below are independent, order does not matter here.
        self.pyscf_type_combo.set_active(self.parameters['type_id'])
        self.pyscf_functional_combo.set_active(self.parameters['functional_id'])
        self.pyscf_basis_group_combo.set_active(self.parameters['basis_group_id'])
        self.pyscf_basis_combo.set_active(self.parameters['basis_id'])
        self.pyscf_ecp_combo.set_active(self.parameters['ecp_id'])
        self.spinbutton_pyscf_ncpus.set_value(self.parameters['ncpu'])
        self.checkbox_pyscf_delete_job_files.set_active(self.parameters['delete_job_files'])

        self.button_ok     = self.builder.get_object('pyscf_button_ok')
        self.button_cancel = self.builder.get_object('pyscf_button_cancel')
        self.button_ok.connect("clicked", self.on_button_ok)
        self.button_cancel.connect("clicked", self.close_window)

        self.window.connect("destroy", self.close_window)
        self.Visible = True

    def on_pyscf_basis_group_combo_changed (self, widget):
        """ Function doc """
        _id = widget.get_active()
        self.pyscf_basis_combo.set_model( self.pyscf_basis_dict[_id] )
        self.pyscf_basis_combo.set_active( 0 )

    def on_combobox_pyscf_type_changed (self, widget):
        """ Function doc """
        _id = widget.get_active()
        restrited = self.setup_QC_model_window.builder.get_object('radio_button_restricted').get_active()
        is_dft = ( _id == 1 )
        self.pyscf_functional_combo.set_sensitive(is_dft)
        if is_dft:
            self.restrited_label = 'RKS' if restrited else 'UKS'
        else:
            self.restrited_label = 'RHF' if restrited else 'UHF'

    def close_window (self, button, data = None):
        """ Function doc """
        if self.window:
            self.window.destroy()
        self.Visible = False

    def _default_scratch_folder (self):
        """ PDYNAMO3_SCRATCH/<system tag>, mirroring SetupXTBWindow/
        SetupORCAWindow's own per-system scratch (so two systems open at
        once don't share PySCF job files); falls back to $HOME if
        PDYNAMO3_SCRATCH is unset. """
        base = os.environ.get('PDYNAMO3_SCRATCH') or os.environ.get('HOME') or ''
        try:
            _id = self.main_session.p_session.active_id
            system = self.main_session.p_session.psystem[_id]
            tag = getattr(system, 'e_tag', None)
        except Exception:
            tag = None
        return os.path.join(base, tag) if tag else base

    def refresh_pyscf_engine_status (self):
        """ Updates label_pyscf_engine_status: pyscf has no external
        executable to locate (it is a Python package), so the check here
        is simply whether it is importable. """
        try:
            import pyscf
            version = getattr(pyscf, '__version__', '?')
            self.label_pyscf_engine_status.set_text('pyscf {} available.'.format(version))
        except ImportError:
            self.label_pyscf_engine_status.set_text(
                '⚠ pyscf is not installed in this Python environment. Calculations with this '
                'setup will fail until it is (pip install pyscf).')

    def on_button_browse_scratch (self, button):
        """ Function doc """
        dialog = Gtk.FileChooserDialog(
            title='Choose the PySCF scratch folder', parent=self.window,
            action=Gtk.FileChooserAction.SELECT_FOLDER)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                           Gtk.STOCK_OPEN,   Gtk.ResponseType.OK)
        current = self.entry_pyscf_scratch.get_text()
        if current and os.path.isdir(current):
            dialog.set_current_folder(current)
        if dialog.run() == Gtk.ResponseType.OK:
            self.entry_pyscf_scratch.set_text(dialog.get_filename())
        dialog.destroy()

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

    def on_button_ok (self, button):
        """ Function doc """
        scratch = self.entry_pyscf_scratch.get_text().strip()
        if not scratch:
            scratch = self._default_scratch_folder()
            self.entry_pyscf_scratch.set_text(scratch)

        if not os.path.isdir(scratch):
            try:
                os.makedirs(scratch, exist_ok=True)
            except Exception as e:
                self._show_message(Gtk.MessageType.ERROR,
                    'Could not create the scratch folder "{}": {}\n\nFix the path and try again.'.format(scratch, e))
                return
        if not os.access(scratch, os.W_OK):
            self._show_message(Gtk.MessageType.ERROR,
                'The scratch folder "{}" is not writable.\n\nFix its permissions (or pick another folder) and try again.'.format(scratch))
            return

        ncpu = self.spinbutton_pyscf_ncpus.get_value_as_int()
        cpu_count = os.cpu_count()
        if cpu_count and ncpu > cpu_count:
            proceed = self._ask_yes_no(
                'Number of CPUs ({}) is higher than this machine\'s {} available cores.\n\nUse it anyway?'.format(ncpu, cpu_count))
            if not proceed:
                return

        try:
            import pyscf
            pyscf_available = True
        except ImportError:
            pyscf_available = False
            self._show_message(Gtk.MessageType.WARNING,
                'pyscf is not installed in this Python environment. '
                'This setup will be saved, but running it will fail until it is (pip install pyscf).')

        type_id = self.pyscf_type_combo.get_active()
        functional_id = self.pyscf_functional_combo.get_active()
        basis_group_id = self.pyscf_basis_group_combo.get_active()
        basis_id = self.pyscf_basis_combo.get_active()
        ecp_id = self.pyscf_ecp_combo.get_active()

        functional_model = self.pyscf_functional_combo.get_model()
        functional = functional_model[functional_id][0].split()[0] if functional_model is not None else ''

        basis_model = self.pyscf_basis_combo.get_model()
        basis = basis_model[basis_id][0].split()[0] if basis_model is not None else ''

        ecp_model = self.pyscf_ecp_combo.get_model()
        ecp = ecp_model[ecp_id][0].split()[0] if ecp_model is not None else 'None'

        custom_functional = self.entry_pyscf_custom_functional.get_text().strip()
        custom_basis      = self.entry_pyscf_custom_basis.get_text().strip()
        if custom_functional:
            functional = custom_functional
        if custom_basis:
            basis = custom_basis

        # . A light sanity check against the real PySCF/libxc tables, so a
        #   typo in a custom functional/basis (or one PySCF genuinely
        #   doesn't know) is caught here instead of deep inside a QC/MM
        #   job's traceback. Advisory only (some legitimate basis sets
        #   only cover certain elements, e.g. heavy-element-only ECP
        #   basis sets, so a single test element can't be a hard gate).
        if pyscf_available:
            problems = []
            if type_id == 1:  # DFT: functional actually matters
                try:
                    pyscf.dft.libxc.parse_xc(functional)
                except Exception as e:
                    problems.append('Functional "{}": {}'.format(functional, e))
            try:
                pyscf.gto.basis.load(basis, 'C')
            except Exception as e:
                problems.append('Basis set "{}" (tested on carbon): {}'.format(basis, e))
            if problems:
                proceed = self._ask_yes_no(
                    "PySCF didn't recognise:\n\n{}\n\nThis may be fine if your system doesn't contain "
                    "carbon and the basis is simply element-limited (e.g. an ECP/-PP set). "
                    "Save this setup anyway?".format('\n'.join(problems)))
                if not proceed:
                    return

        self.parameters['type_id']           = type_id
        self.parameters['functional_id']     = functional_id
        self.parameters['basis_group_id']    = basis_group_id
        self.parameters['basis_id']          = basis_id
        self.parameters['ecp_id']            = ecp_id
        self.parameters['custom_functional'] = custom_functional
        self.parameters['custom_basis']      = custom_basis
        self.parameters['ncpu']              = ncpu
        self.parameters['scratch']           = scratch
        self.parameters['delete_job_files']  = self.checkbox_pyscf_delete_job_files.get_active()

        # . What session.py actually needs to build QCModelPySCF.WithOptions().
        self.parameters['method']       = self.restrited_label  # RHF/UHF/RKS/UKS
        self.parameters['functional']   = functional.lower() if type_id == 1 else None
        self.parameters['orbitalBasis'] = basis
        self.parameters['ecp']          = None if ecp_id == 0 else ecp

        dprint(self.parameters)
        self.close_window(None)
