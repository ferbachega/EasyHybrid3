#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Run NAMD window
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
#      "Run NAMD" window -- configures and runs NAMD molecular dynamics
#      jobs directly from AMBER prmtop/inpcrd or CHARMM psf/parameters
#      input (a file, or the current frame of a vobject already loaded
#      in EasyHybrid), either a single custom run or a chained
#      "Custom Pipeline" of steps built in the treeview (restart-
#      chained the same way NAMD_protocol/AMBER_protocol's own 9-step
#      real-world equilibration protocol is -- offered there as a
#      loadable preset, see _PRESETS below). Any pipeline row can also
#      be hand-edited as raw .namd text via "Setup step..." (see
#      on_button_setup_step_clicked()) when the treeview's own five
#      columns aren't enough -- a "Custom" row is used AS-IS at run
#      time instead of being rebuilt from those columns. Pipelines can
#      be saved to / loaded from a JSON file (including any "Setup
#      step..." overrides). NAMD runs as a real subprocess in the
#      background (never blocking the GTK main loop), polled via
#      GLib.timeout_add the same way process_manager_window.py already
#      polls job status; Stop reuses that same module's psutil-based
#      process-tree kill helpers.
#

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk, GLib, Pango

import json
import os
import re
import time
import traceback

from gui.widgets.custom_widgets import FolderChooserButton
from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.windows.setup.process_manager_window import _terminate_process_tree, _force_kill_pids, _pid_alive
from pdynamo.pDynamo2EasyHybrid.helpers import export_special_PDB
from util import namd_runner


# Display <-> internal-id mappings shared by the "Add current settings as
# step" button and the pipeline treeview's own editable combo cells (see
# liststore_pipeline_ensemble_options/liststore_pipeline_restraint_options
# in the glade -- their row text IS these display strings, so a
# GtkCellRendererCombo edit always hands back one of them verbatim).
_ENSEMBLE_DISPLAY = {'minimization': 'Minimization', 'nve': 'NVE', 'nvt': 'NVT', 'npt': 'NPT'}
_ENSEMBLE_FROM_DISPLAY = {v: k for k, v in _ENSEMBLE_DISPLAY.items()}
_RESTRAINT_DISPLAY = {'none': 'None', 'solute': 'Solute', 'backbone': 'Backbone'}
_RESTRAINT_FROM_DISPLAY = {v: k for k, v in _RESTRAINT_DISPLAY.items()}

# Built-in pipeline starting points, offered via combo_pipeline_presets +
# button_load_pipeline_preset -- each value is a list of (ensemble,
# temperature, steps, restraint_mode, scaling) tuples, the same 5 fields
# a pipeline row holds (see on_button_add_pipeline_step_clicked()).
#
# "AMBER Equilibration Protocol" replicates the real 9-step protocol at
# ~/programs/NAMD_protocol/AMBER_protocol/ (already cross-checked
# elsewhere in this codebase -- see namd_runner.py's module docstring)
# AS CLOSELY AS THIS GENERIC PIPELINE MODEL CAN: 8 of its 9 steps map
# over exactly (ensemble, constraintScaling schedule 100/100/100/10/10/
# 10/1/0.1/off, conskcol O for steps 1-4 vs. B for 5-8 -- here "Solute"
# vs. "Backbone"). The one step this model genuinely cannot reproduce
# is the original step 2's single continuous 100-increment Tcl
# `langevinTemp`/`run` loop heating 100K -> 298K WITHIN one NAMD job --
# a pipeline row is always a separate NAMD run. Approximated here as 5
# separate NVT rows at increasing target temperatures instead (150K,
# 194K, 238K, 282K, 298K -- evenly spaced, restart-chained normally),
# which needs no special-casing anywhere in the pipeline execution
# engine, at the cost of a coarser (not perfectly continuous) ramp --
# and, since every row here inherits binVelocities from the row before
# it as a matter of course (unlike the original protocol's real step 2,
# which deliberately does NOT inherit velocities so it can start a
# fresh Maxwell-Boltzmann assignment at 100K), velocities simply carry
# forward from the minimization step instead of being reinitialized.
_PRESETS = {
    'AMBER Equilibration Protocol (9 steps, approximated)': [
        ('minimization', 298.0, 2000, 'solute', 100.0),
        ('nvt', 150.0, 4000, 'solute', 100.0),
        ('nvt', 194.0, 4000, 'solute', 100.0),
        ('nvt', 238.0, 4000, 'solute', 100.0),
        ('nvt', 282.0, 4000, 'solute', 100.0),
        ('nvt', 298.0, 4000, 'solute', 100.0),
        ('npt', 298.0, 1000, 'solute', 100.0),
        ('npt', 298.0, 1000, 'solute', 10.0),
        ('minimization', 298.0, 1000, 'backbone', 10.0),
        ('npt', 298.0, 1000, 'backbone', 10.0),
        ('npt', 298.0, 1000, 'backbone', 1.0),
        ('npt', 298.0, 1000, 'backbone', 0.1),
        ('npt', 298.0, 1000, 'none', 0.0),
    ],
    'Quick Minimization': [
        ('minimization', 298.0, 5000, 'none', 0.0),
    ],
    'Quick NPT Equilibration': [
        ('minimization', 298.0, 2000, 'none', 0.0),
        ('npt', 298.0, 5000, 'none', 0.0),
    ],
}

# Per-line highlight rules for the "Log" tab -- checked in this order
# (first match wins) against every line NAMD prints. Patterns matched
# against real NAMD output seen this session (both success and FATAL
# ERROR runs) -- see namd_runner.py's check_namd_log_success() for the
# same "FATAL"/"ERROR:" markers used there to decide pass/fail.
_LOG_ERROR_RE = re.compile(r'FATAL|ERROR:|Error!', re.IGNORECASE)
_LOG_WARNING_RE = re.compile(r'Warning', re.IGNORECASE)
_LOG_SUCCESS_RE = re.compile(r'WRITING COORDINATES|End of program|WallClock:')

# Per-line highlight rules for the "Setup step..." editor (see
# on_button_setup_step_clicked()/_highlight_step_setup_text()) -- every
# NAMD keyword build_namd_config() ever emits, split into two buckets:
#
# "Safe" -- ordinary tuning knobs. Changing one of these on a single
# pipeline step (temperature, timestep, cutoff, output frequencies,
# restraint scaling, Langevin/piston tuning, ensemble switches) affects
# only that step's own physics -- it can't break the chain to the next
# step or crash NAMD outright.
#
# "Critical" -- structural/chaining/force-field-intrinsic values.
# Getting one of these wrong either breaks continuity with the next
# step silently (outputName must keep matching the basename this
# window itself expects for THIS step -- see
# _build_pipeline_steps()'s own docstring; binCoordinates/binVelocities/
# extendedSystem must keep pointing at the PREVIOUS step's real files),
# or can crash NAMD the same way the real "PMEGridSizeX too small"
# FATAL ERROR did this session (PME*/cellBasisVector*/cellOrigin), or
# silently changes which force field convention is in effect
# (amber/readexclusions/scnb/oneFourScaling/exclude/switching/
# paraTypeCharmm are all intrinsic to AMBER-vs-CHARMM mode, not
# simulation choices -- see build_namd_config()'s own docstring).
_NAMD_SAFE_KEYWORDS = (
    'temperature', 'timeStep', 'cutoff', 'pairListDist',
    'langevinDamping', 'langevinTemp', 'langevinHydrogen',
    'langevinPistonTarget', 'langevinPistonPeriod', 'langevinPistonDecay',
    'langevinPistonTemp', 'restartfreq', 'dcdfreq', 'outputEnergies',
    'outputPressure', 'minimize', 'run', 'constraintScaling',
    'useFlexibleCell', 'useConstantArea', 'useGroupPressure', 'watermodel',
    'langevin', 'langevinPiston',
)
_NAMD_CRITICAL_KEYWORDS = (
    'amber', 'parmfile', 'ambercoor', 'structure', 'parameters',
    'paraTypeCharmm', 'coordinates', 'binCoordinates', 'binVelocities',
    'extendedSystem', 'outputName', 'cellBasisVector1', 'cellBasisVector2',
    'cellBasisVector3', 'cellOrigin', 'wrapAll', 'PME', 'PMEGridSpacing',
    'PMEGridSizeX', 'PMEGridSizeY', 'PMEGridSizeZ', 'PMETolerance',
    'PMEInterpOrder', 'readexclusions', 'scnb', 'oneFourScaling', 'exclude',
    'switching', 'switchdist', 'rigidBonds', 'rigidTolerance',
    'rigidIterations', 'useSettle', 'fullElectFrequency', 'nonBondedFreq',
    'stepspercycle', 'consref', 'conskfile', 'conskcol', 'constraints',
    'firsttimestep',
)
_NAMD_SAFE_RE = re.compile(
    r'^\s*(?:' + '|'.join(re.escape(k) for k in _NAMD_SAFE_KEYWORDS) + r')\b')
_NAMD_CRITICAL_RE = re.compile(
    r'^\s*(?:' + '|'.join(re.escape(k) for k in _NAMD_CRITICAL_KEYWORDS) + r')\b')


class PrepareNamdRunWindow:
    """ "Run NAMD" window. """

    def __init__(self, main=None):
        """ Class initialiser """
        self.main       = main
        self.vm_session = main.vm_session
        self.p_session  = main.p_session
        self.home       = main.home
        self.Visible    = False

        # . Job-in-progress state -- only one job (a single run, or one
        #   step of a pipeline) is ever running at a time from this
        #   window, whichever mode it's in.
        self.current_process   = None   # subprocess.Popen or None
        self.current_log_path  = None
        self._log_read_offset  = 0
        self._poll_timeout_id  = None
        self._protocol_steps   = None   # list of dicts built by on_button_run_pipeline_clicked(), or None (single-run mode)
        self._protocol_index   = None   # index into self._protocol_steps of the step currently running
        self._input_files      = None   # dict from _validate_input_files() -- describes the AMBER or CHARMM input the CURRENT run started from
        self._protocol_system_name = None
        self._run_mode         = None   # 'single' or 'pipeline' -- which import-back path to use when the job ends

        # . "Setup step..." editor state -- which pipeline row (by index
        #   into liststore_pipeline) window_step_setup is currently
        #   editing, and that step's freshly auto-generated text (for
        #   "Reset to auto-generated" -- see on_button_setup_step_clicked()).
        self._step_setup_row_index = None
        self._step_setup_auto_text = None

    #-------------------------------------------------------------------------
    #  W I N D O W   L I F E C Y C L E
    #-------------------------------------------------------------------------
    def open_window(self):
        """ Function doc """
        if self.Visible:
            self.window.present()
            return

        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/setup/prepare_namd_run_window.glade'))
        self.builder.connect_signals(self)

        self.window = self.builder.get_object('window')
        self.window.set_title('Run NAMD')

        # -------------------- widget shortcuts --------------------
        self.radio_input_amber  = self.builder.get_object('radio_input_amber')
        self.radio_input_charmm = self.builder.get_object('radio_input_charmm')
        self.box_amber_input    = self.builder.get_object('box_amber_input')
        self.box_charmm_input   = self.builder.get_object('box_charmm_input')
        self.entry_prmtop_path = self.builder.get_object('entry_prmtop_path')
        self.entry_inpcrd_path = self.builder.get_object('entry_inpcrd_path')
        self.combo_amber_coords_source = self.builder.get_object('combo_amber_coords_source')
        # . No single wrapper box for either state any more (each is its
        #   own row of individually-shown/hidden widgets instead) -- see
        #   on_amber_coords_source_changed()'s own comment for why.
        self.label_amber_inp = self.builder.get_object('label_amber_inp')
        self.button_browse_inpcrd = self.builder.get_object('button_browse_inpcrd')
        self.label_system = self.builder.get_object('label_system')
        self.label_object = self.builder.get_object('label_object')
        self.box_amber_coords_system = self.builder.get_object('box_amber_coords_system')
        self.box_amber_coords_object = self.builder.get_object('box_amber_coords_object')
        self.entry_charmm_psf_path = self.builder.get_object('entry_charmm_psf_path')
        self.entry_charmm_coordinates_path = self.builder.get_object('entry_charmm_coordinates_path')
        self.combo_charmm_coords_source = self.builder.get_object('combo_charmm_coords_source')
        self.box_charmm_coords_file = self.builder.get_object('box_charmm_coords_file')
        self.box_charmm_coords_vobject = self.builder.get_object('box_charmm_coords_vobject')
        self.treeview_charmm_parameters = self.builder.get_object('treeview_charmm_parameters')
        self.liststore_charmm_parameters = self.builder.get_object('liststore_charmm_parameters')
        self.checkbox_charmm_xplor = self.builder.get_object('checkbox_charmm_xplor')
        self.entry_namd_command = self.builder.get_object('entry_namd_command')
        self.spinbtn_n_procs   = self.builder.get_object('spinbtn_n_procs')
        self.entry_system_name = self.builder.get_object('entry_system_name')
        self.label_info        = self.builder.get_object('label_info')
        self.textview_log      = self.builder.get_object('textview_log')
        log_buffer = self.textview_log.get_buffer()
        self._log_tag_error = log_buffer.create_tag('namd_log_error', foreground='#e01b24', weight=Pango.Weight.BOLD)
        self._log_tag_warning = log_buffer.create_tag('namd_log_warning', foreground='#c47a00')
        self._log_tag_success = log_buffer.create_tag('namd_log_success', foreground='#26a269', weight=Pango.Weight.BOLD)
        self.button_stop        = self.builder.get_object('button_stop')

        self.radio_mode_single = self.builder.get_object('radio_mode_single')
        self.radio_mode_pipeline = self.builder.get_object('radio_mode_pipeline')
        self.frame_restart_single = self.builder.get_object('frame_restart_single')
        self.frame_custom_pipeline = self.builder.get_object('frame_custom_pipeline')
        self.button_add_pipeline_step = self.builder.get_object('button_add_pipeline_step')

        self.spinbtn_temperature = self.builder.get_object('spinbtn_temperature')
        self.spinbtn_timestep    = self.builder.get_object('spinbtn_timestep')
        self.spinbtn_cutoff      = self.builder.get_object('spinbtn_cutoff')
        self.combo_watermodel    = self.builder.get_object('combo_watermodel')
        self.radio_ensemble_minimization = self.builder.get_object('radio_ensemble_minimization')
        self.radio_ensemble_nve  = self.builder.get_object('radio_ensemble_nve')
        self.radio_ensemble_nvt  = self.builder.get_object('radio_ensemble_nvt')
        self.radio_ensemble_npt  = self.builder.get_object('radio_ensemble_npt')
        self.spinbtn_minimize_steps = self.builder.get_object('spinbtn_minimize_steps')
        self.spinbtn_run_steps      = self.builder.get_object('spinbtn_run_steps')
        self.combo_restraint_mode      = self.builder.get_object('combo_restraint_mode')
        self.spinbtn_constraint_scaling = self.builder.get_object('spinbtn_constraint_scaling')
        self.label_constraint_scaling   = self.builder.get_object('label_constraint_scaling')
        self.checkbox_restart_single = self.builder.get_object('checkbox_restart_single')
        self.grid_restart_single     = self.builder.get_object('grid_restart_single')
        self.entry_restart_coor = self.builder.get_object('entry_restart_coor')
        self.entry_restart_vel  = self.builder.get_object('entry_restart_vel')
        self.entry_restart_xsc  = self.builder.get_object('entry_restart_xsc')
        self.button_run_single  = self.builder.get_object('button_run_single')

        self.treeview_pipeline = self.builder.get_object('treeview_pipeline')
        self.liststore_pipeline = self.builder.get_object('liststore_pipeline')
        self.button_run_pipeline = self.builder.get_object('button_run_pipeline')
        self.combo_pipeline_presets = self.builder.get_object('combo_pipeline_presets')

        # . Everything that EDITS the pipeline (as opposed to running it,
        #   or loading/saving a whole pipeline -- both stay as regular
        #   buttons) used to be its own row of buttons below the table;
        #   moved into a right-click context menu instead, same pattern
        #   process_manager_window.py's own job-history treeview already
        #   uses (build the Gtk.Menu here, catch right-clicks via
        #   "button-press-event", see _on_pipeline_treeview_button_press()).
        #   Reuses the very same handler methods the old buttons called
        #   -- a GtkMenuItem's "activate" signal hands back just the
        #   widget, same shape as GtkButton's "clicked".
        self._pipeline_popup_menu = Gtk.Menu()
        for label, handler in (
            ('Move up', self.on_button_pipeline_move_up_clicked),
            ('Move down', self.on_button_pipeline_move_down_clicked),
            ('Remove step', self.on_button_remove_pipeline_step_clicked),
            ('Clear pipeline', self.on_button_clear_pipeline_clicked),
            (None, None),
            ('Setup step...', self.on_button_setup_step_clicked),
        ):
            if label is None:
                self._pipeline_popup_menu.append(Gtk.SeparatorMenuItem())
                continue
            menu_item = Gtk.MenuItem(label=label)
            menu_item.connect('activate', handler)
            self._pipeline_popup_menu.append(menu_item)
        self._pipeline_popup_menu.show_all()
        self.treeview_pipeline.connect('button-press-event', self._on_pipeline_treeview_button_press)

        # -------------------- "Setup step..." editor --------------------
        self.window_step_setup = self.builder.get_object('window_step_setup')
        self.label_step_setup_title = self.builder.get_object('label_step_setup_title')
        self.textview_step_setup = self.builder.get_object('textview_step_setup')
        step_setup_buffer = self.textview_step_setup.get_buffer()
        self._step_setup_tag_safe = step_setup_buffer.create_tag('namd_var_safe', foreground='#26a269')
        self._step_setup_tag_critical = step_setup_buffer.create_tag(
            'namd_var_critical', foreground='#c47a00', weight=Pango.Weight.BOLD)

        self.checkbox_pbc = self.builder.get_object('checkbox_pbc')
        self.frame_periodic_cell = self.builder.get_object('frame_periodic_cell')
        self.frame_implicit_solvent = self.builder.get_object('frame_implicit_solvent')
        self.checkbox_implicit_solvent = self.builder.get_object('checkbox_implicit_solvent')
        self.grid_implicit_solvent = self.builder.get_object('grid_implicit_solvent')
        self.spinbtn_solvent_dielectric = self.builder.get_object('spinbtn_solvent_dielectric')
        self.spinbtn_ion_concentration = self.builder.get_object('spinbtn_ion_concentration')
        self.spinbtn_alpha_cutoff = self.builder.get_object('spinbtn_alpha_cutoff')
        self.checkbox_sasa = self.builder.get_object('checkbox_sasa')
        self.spinbtn_surface_tension = self.builder.get_object('spinbtn_surface_tension')

        self.checkbox_manual_cell = self.builder.get_object('checkbox_manual_cell')
        self.grid_manual_cell     = self.builder.get_object('grid_manual_cell')
        self.spinbtn_cell_a = self.builder.get_object('spinbtn_cell_a')
        self.spinbtn_cell_b = self.builder.get_object('spinbtn_cell_b')
        self.spinbtn_cell_c = self.builder.get_object('spinbtn_cell_c')
        self.spinbtn_cell_origin_x = self.builder.get_object('spinbtn_cell_origin_x')
        self.spinbtn_cell_origin_y = self.builder.get_object('spinbtn_cell_origin_y')
        self.spinbtn_cell_origin_z = self.builder.get_object('spinbtn_cell_origin_z')
        self.spinbtn_pme_grid_x = self.builder.get_object('spinbtn_pme_grid_x')
        self.spinbtn_pme_grid_y = self.builder.get_object('spinbtn_pme_grid_y')
        self.spinbtn_pme_grid_z = self.builder.get_object('spinbtn_pme_grid_z')
        self.checkbox_flexible_cell = self.builder.get_object('checkbox_flexible_cell')
        self.checkbox_wrap_all = self.builder.get_object('checkbox_wrap_all')

        self.spinbtn_restartfreq = self.builder.get_object('spinbtn_restartfreq')
        self.spinbtn_dcdfreq = self.builder.get_object('spinbtn_dcdfreq')
        self.spinbtn_outputenergies = self.builder.get_object('spinbtn_outputenergies')
        self.spinbtn_outputpressure = self.builder.get_object('spinbtn_outputpressure')

        self.spinbtn_langevin_damping = self.builder.get_object('spinbtn_langevin_damping')
        self.checkbox_langevin_hydrogen = self.builder.get_object('checkbox_langevin_hydrogen')
        self.checkbox_override_langevin_temp = self.builder.get_object('checkbox_override_langevin_temp')
        self.spinbtn_langevin_temp = self.builder.get_object('spinbtn_langevin_temp')

        # -------------------- Force Field (AMBER/CHARMM, NAMD ug/node13.html) --------------------
        self.combo_exclude = self.builder.get_object('combo_exclude')
        self.spinbtn_one_four_scaling = self.builder.get_object('spinbtn_one_four_scaling')
        self.checkbox_switching = self.builder.get_object('checkbox_switching')
        self.label_switchdist = self.builder.get_object('label_switchdist')
        self.spinbtn_switchdist = self.builder.get_object('spinbtn_switchdist')
        self.combo_rigidbonds = self.builder.get_object('combo_rigidbonds')
        self.spinbtn_rigidtolerance = self.builder.get_object('spinbtn_rigidtolerance')
        self.checkbox_vdw_force_switching = self.builder.get_object('checkbox_vdw_force_switching')
        self.checkbox_readexclusions = self.builder.get_object('checkbox_readexclusions')
        self.spinbtn_scnb = self.builder.get_object('spinbtn_scnb')

        # -------------------- coordinates-from-vobject pickers (AMBER + CHARMM) --------------------
        self.amber_coords_system_combo = SystemComboBox(self.main)
        self.amber_coords_system_combo.connect("changed", self._on_amber_coords_system_changed)
        self.box_amber_coords_system.pack_start(self.amber_coords_system_combo, False, False, 0)
        self.amber_coords_object_combo = CoordinatesComboBox()
        self.box_amber_coords_object.pack_start(self.amber_coords_object_combo, False, False, 0)

        self.charmm_coords_system_combo = SystemComboBox(self.main)
        self.charmm_coords_system_combo.connect("changed", self._on_charmm_coords_system_changed)
        self.builder.get_object('box_charmm_coords_system').pack_start(self.charmm_coords_system_combo, False, False, 0)
        self.charmm_coords_object_combo = CoordinatesComboBox()
        self.builder.get_object('box_charmm_coords_object').pack_start(self.charmm_coords_object_combo, False, False, 0)

        # -------------------- working-folder chooser --------------------
        self.box_work_folder = self.builder.get_object('box_work_folder')
        self.work_folder_chooser = FolderChooserButton(self.main, 'folder', self.home)
        self.box_work_folder.pack_start(self.work_folder_chooser.btn, False, False, 0)
        scratch = os.environ.get('PDYNAMO3_SCRATCH')
        if scratch and os.path.isdir(scratch):
            self.work_folder_chooser.set_folder(folder=scratch)

        # -------------------- namd executable --------------------
        namd_command = self.main.vm_session.vm_config.gl_parameters.get('namd_command')
        if not namd_command or not os.path.isfile(namd_command):
            namd_command = namd_runner.find_namd_executable()
        self.entry_namd_command.set_text(namd_command or '')
        if not namd_command:
            self.label_info.set_text(
                'namd3/namd2 was not found automatically on PATH -- please enter the full path below.')

        # -------------------- CHARMM parameters (kept in Python, mirrored into the treeview) --------------------
        self.charmm_parameters = []

        # -------------------- pipeline presets --------------------
        self.combo_pipeline_presets.remove_all()
        for name in _PRESETS:
            self.combo_pipeline_presets.append_text(name)
        self.combo_pipeline_presets.set_active(0)

        self.window.show_all()
        # . Gtk.Window.show_all() forces EVERY descendant widget
        #   visible, silently overriding any visible="False" set in the
        #   glade file (confirmed the hard way: grid_restart_single's
        #   own visible="False" was being ignored) -- every
        #   conditionally-shown widget has to be re-synced to its
        #   actual starting state here, same as box_charmm_input
        #   already was.
        self.on_input_type_radio_toggled(self.radio_input_amber)
        self.on_amber_coords_source_changed(self.combo_amber_coords_source)
        self.on_namd_mode_radio_toggled(self.radio_mode_single)
        self.on_checkbox_switching_toggled(self.checkbox_switching)
        self.on_checkbox_pbc_toggled(self.checkbox_pbc)
        self.on_checkbox_sasa_toggled(self.checkbox_sasa)
        self.on_checkbox_implicit_solvent_toggled(self.checkbox_implicit_solvent)
        self.box_charmm_coords_file.set_visible(self.combo_charmm_coords_source.get_active_id() != 'vobject')
        self.box_charmm_coords_vobject.set_visible(self.combo_charmm_coords_source.get_active_id() == 'vobject')
        self.grid_restart_single.set_visible(self.checkbox_restart_single.get_active())
        has_restraints = self.combo_restraint_mode.get_active_id() != 'none'
        self.spinbtn_constraint_scaling.set_visible(has_restraints)
        self.label_constraint_scaling.set_visible(has_restraints)
        # . show_all() above forced every hidden widget momentarily
        #   visible (CHARMM box, restart grid, ...), which sizes the
        #   window to fit ALL of that -- GTK does not auto-shrink a
        #   toplevel back down once it's grown, even after the widgets
        #   are re-hidden right above. reshow_with_initial_size() resets
        #   the window back to its glade default-width/default-height
        #   instead of leaving a huge dead area below the real content.
        self.window.reshow_with_initial_size()
        self.window.connect('destroy', self.close_window)
        self.Visible = True

    def close_window(self, button=None, data=None):
        """ Function doc """
        if not self.Visible:
            return
        self.window.destroy()
        self.Visible = False

    #-------------------------------------------------------------------------
    #  H E L P E R S
    #-------------------------------------------------------------------------
    def _browse_file(self, entry, title, patterns, pattern_name):
        dialog = Gtk.FileChooserDialog(
            title=title, parent=self.window, action=Gtk.FileChooserAction.OPEN,
        )
        dialog.add_buttons(
            Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
            Gtk.STOCK_OPEN, Gtk.ResponseType.OK,
        )
        if patterns:
            file_filter = Gtk.FileFilter()
            file_filter.set_name(pattern_name)
            for pattern in patterns:
                file_filter.add_pattern(pattern)
            dialog.add_filter(file_filter)
        if dialog.run() == Gtk.ResponseType.OK:
            entry.set_text(dialog.get_filename())
        dialog.destroy()

    def _append_log(self, text):
        """ Appends `text` to the Log tab, tagging each LINE (not
            substring -- simpler, and reads fine for whole-line NAMD
            status/error messages) with a colour per _LOG_ERROR_RE/
            _LOG_WARNING_RE/_LOG_SUCCESS_RE, first match wins, plain
            for everything else. `text` is whatever _tail_log() just
            read off the log file -- not necessarily line-aligned, so
            a keyword that happens to straddle two consecutive reads
            may go untagged just that once; purely cosmetic, every
            byte still gets inserted either way.
        """
        buf = self.textview_log.get_buffer()
        for line in text.splitlines(keepends=True):
            if _LOG_ERROR_RE.search(line):
                tag = self._log_tag_error
            elif _LOG_WARNING_RE.search(line):
                tag = self._log_tag_warning
            elif _LOG_SUCCESS_RE.search(line):
                tag = self._log_tag_success
            else:
                tag = None
            end_iter = buf.get_end_iter()
            if tag is not None:
                buf.insert_with_tags(end_iter, line, tag)
            else:
                buf.insert(end_iter, line)
        self.textview_log.scroll_to_iter(buf.get_end_iter(), 0.0, False, 0.0, 0.0)

    def _highlight_step_setup_text(self, text):
        """ Fills the "Setup step..." editor with `text`, tagging each
            LINE green (_NAMD_SAFE_RE) or orange-bold (_NAMD_CRITICAL_RE)
            by its leading NAMD keyword -- see those two regexes' own
            comment for the safe/critical split and why each keyword
            landed where it did. First match wins per line (a line can
            only start with one keyword); no tag at all for lines this
            module doesn't recognise (Tcl `set`/comments/blanks/the
            heating loop's `for`/`langevinTemp $current_temp`/`run`
            inside it, etc.) -- same "best-effort, not a real parser"
            spirit as the Log tab's own highlighting.
        """
        buf = self.textview_step_setup.get_buffer()
        buf.set_text(text)
        for line_no, line in enumerate(text.splitlines()):
            if _NAMD_SAFE_RE.match(line):
                tag = self._step_setup_tag_safe
            elif _NAMD_CRITICAL_RE.match(line):
                tag = self._step_setup_tag_critical
            else:
                continue
            line_start = buf.get_iter_at_line(line_no)
            line_end = line_start.copy()
            line_end.forward_to_line_end()
            buf.apply_tag(tag, line_start, line_end)

    def _job_running(self):
        return self.current_process is not None

    def _set_controls_sensitive(self, sensitive):
        self.button_run_single.set_sensitive(sensitive)
        self.button_run_pipeline.set_sensitive(sensitive)
        self.button_stop.set_sensitive(not sensitive)

    #-------------------------------------------------------------------------
    #  S I G N A L S :  A M B E R   I N P U T
    #-------------------------------------------------------------------------
    def on_button_browse_prmtop_clicked(self, widget):
        self._browse_file(self.entry_prmtop_path, "Select AMBER topology",
                           ["*.prmtop", "*.top", "*.parm7"], "AMBER topology (*.prmtop, *.top, *.parm7)")

    def on_entry_prmtop_path_changed(self, widget):
        """ Prefills the Advanced tab's manual-cell fields (Box a/b/c,
            origin, PME grid) the moment a prmtop is picked (Browse, or
            typed/pasted directly -- GtkEntry.set_text() itself fires
            "changed", so Browse needs no separate call) FROM the
            prmtop's own "%FLAG BOX_DIMENSIONS" section -- see
            namd_runner.cell_from_prmtop_box_dimensions()'s docstring
            for the format and why cellOrigin here is only an
            approximation. Purely a convenience default: it does NOT
            check "Manually set periodic cell" for the user, so a Run
            still uses the coordinate-accurate compute_cell_from_coordinates()
            unless they opt in by checking that box themselves.
        """
        self._autofill_cell_from_prmtop(widget.get_text().strip())

    def _autofill_cell_from_prmtop(self, prmtop_path):
        if not prmtop_path or not os.path.isfile(prmtop_path):
            return
        try:
            cell = namd_runner.cell_from_prmtop_box_dimensions(prmtop_path)
        except (OSError, ValueError, IndexError):
            return
        if cell is None:
            return

        self.spinbtn_cell_a.set_value(cell['cell_basis_vector1'][0])
        self.spinbtn_cell_b.set_value(cell['cell_basis_vector2'][1])
        self.spinbtn_cell_c.set_value(cell['cell_basis_vector3'][2])
        self.spinbtn_cell_origin_x.set_value(cell['cell_origin'][0])
        self.spinbtn_cell_origin_y.set_value(cell['cell_origin'][1])
        self.spinbtn_cell_origin_z.set_value(cell['cell_origin'][2])
        self.spinbtn_pme_grid_x.set_value(cell['pme_grid_size'][0])
        self.spinbtn_pme_grid_y.set_value(cell['pme_grid_size'][1])
        self.spinbtn_pme_grid_z.set_value(cell['pme_grid_size'][2])

    def on_button_browse_inpcrd_clicked(self, widget):
        self._browse_file(self.entry_inpcrd_path, "Select AMBER coordinates",
                           ["*.inpcrd", "*.crd", "*.rst7"], "AMBER coordinates (*.inpcrd, *.crd, *.rst7)")

    #-------------------------------------------------------------------------
    #  S I G N A L S :  C H A R M M   I N P U T
    #-------------------------------------------------------------------------
    def on_input_type_radio_toggled(self, widget):
        if not widget.get_active():
            return   # fires for both the button losing AND gaining active state
        is_amber = self.radio_input_amber.get_active()
        self.box_amber_input.set_visible(is_amber)
        self.box_charmm_input.set_visible(self.radio_input_charmm.get_active())
        # . readexclusions/scnb (Setup tab's Force Field section) only
        #   mean anything in `amber yes` mode -- greyed out, not hidden,
        #   in CHARMM mode as a hint they won't apply (build_namd_config()
        #   itself already omits both entirely for CHARMM regardless of
        #   what these are set to).
        self.checkbox_readexclusions.set_sensitive(is_amber)
        self.spinbtn_scnb.set_sensitive(is_amber)

    def on_amber_coords_source_changed(self, widget):
        """ Unlike the CHARMM equivalent (on_charmm_coords_source_changed()
            below), AMBER's file-vs-vobject rows have no single wrapper
            box to toggle any more -- each widget in each row is shown/
            hidden individually. """
        is_vobject = widget.get_active_id() == 'vobject'
        self.label_amber_inp.set_visible(not is_vobject)
        self.entry_inpcrd_path.set_visible(not is_vobject)
        self.button_browse_inpcrd.set_visible(not is_vobject)

        self.label_system.set_visible(is_vobject)
        self.label_object.set_visible(is_vobject)
        self.box_amber_coords_system.set_visible(is_vobject)
        self.box_amber_coords_object.set_visible(is_vobject)

    def on_charmm_coords_source_changed(self, widget):
        is_vobject = widget.get_active_id() == 'vobject'
        self.box_charmm_coords_file.set_visible(not is_vobject)
        self.box_charmm_coords_vobject.set_visible(is_vobject)

    def _on_amber_coords_system_changed(self, widget):
        """ Mirrors prepare_amber_system.py's own
            on_combobox_systemsbox_changed() -- populates the "Object /
            frame" combobox with the vobjects belonging to whichever
            system was just selected, defaulting to the last one. """
        system_id = self.amber_coords_system_combo.get_system_id()
        if system_id is not None:
            self.amber_coords_object_combo.set_model(self.main.vobject_liststore_dict[system_id])
            size = len(list(self.main.vobject_liststore_dict[system_id]))
            self.amber_coords_object_combo.set_active(size - 1)

    def _on_charmm_coords_system_changed(self, widget):
        system_id = self.charmm_coords_system_combo.get_system_id()
        if system_id is not None:
            self.charmm_coords_object_combo.set_model(self.main.vobject_liststore_dict[system_id])
            size = len(list(self.main.vobject_liststore_dict[system_id]))
            self.charmm_coords_object_combo.set_active(size - 1)

    def on_button_browse_charmm_psf_clicked(self, widget):
        self._browse_file(self.entry_charmm_psf_path, "Select CHARMM psf", ["*.psf"], "CHARMM topology (*.psf)")

    def on_button_browse_charmm_coordinates_clicked(self, widget):
        self._browse_file(self.entry_charmm_coordinates_path, "Select starting coordinates",
                           ["*.pdb"], "Coordinates (*.pdb)")

    def on_button_add_charmm_parameter_clicked(self, widget):
        dialog = Gtk.FileChooserDialog(
            title="Add CHARMM parameter file(s) (.prm/.par/.str/.xplor)",
            parent=self.window, action=Gtk.FileChooserAction.OPEN,
        )
        dialog.set_select_multiple(True)
        dialog.add_buttons(
            Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
            Gtk.STOCK_OPEN, Gtk.ResponseType.OK,
        )
        file_filter = Gtk.FileFilter()
        file_filter.set_name("CHARMM parameter files (*.prm, *.par, *.str, *.xplor, *.rtf)")
        for pattern in ("*.prm", "*.par", "*.str", "*.xplor", "*.rtf"):
            file_filter.add_pattern(pattern)
        dialog.add_filter(file_filter)

        if dialog.run() == Gtk.ResponseType.OK:
            for path in dialog.get_filenames():
                self.charmm_parameters.append(path)
                self.liststore_charmm_parameters.append([os.path.basename(path)])
        dialog.destroy()

    def on_button_clear_charmm_parameters_clicked(self, widget):
        self.charmm_parameters = []
        self.liststore_charmm_parameters.clear()

    #-------------------------------------------------------------------------
    #  S I G N A L S :  S I N G L E   R U N   T A B
    #-------------------------------------------------------------------------
    def on_namd_mode_radio_toggled(self, widget):
        """ "Single run" vs "Pipeline" -- the Temperature/Timestep/
            Cutoff/Watermodel/Ensemble/Steps/Restraints grid above stays
            visible either way (it's shared: used directly for a single
            run, or as the "current settings" staged by "Add current
            settings as step" into the pipeline table). Only the parts
            that are meaningless in the OTHER mode are hidden: the
            restart-from-previous-run frame and "Run single" only make
            sense for a standalone run; "Add current settings as step"
            and the whole Custom Pipeline table only make sense when
            actually building a chain.
        """
        if not widget.get_active():
            return   # fires for both the button losing AND gaining active state
        is_pipeline = self.radio_mode_pipeline.get_active()
        self.frame_restart_single.set_visible(not is_pipeline)
        self.button_run_single.set_visible(not is_pipeline)
        self.button_add_pipeline_step.set_visible(is_pipeline)
        self.frame_custom_pipeline.set_visible(is_pipeline)

    def on_combo_restraint_mode_changed(self, widget):
        # . Hidden, not just greyed out, when unused -- "Scaling" only
        #   means anything once a restraint mode is actually picked
        #   (see the window-wide compacting pass this followed).
        has_restraints = widget.get_active_id() != 'none'
        self.spinbtn_constraint_scaling.set_visible(has_restraints)
        self.label_constraint_scaling.set_visible(has_restraints)

    def on_ensemble_radio_toggled(self, widget):
        """ "Minimization" has no dynamics at all -- run_steps is
            meaningless (and, if left non-zero, NAMD would still run
            that many dynamics steps AFTER the minimization, which is
            not what "Minimization" is supposed to mean here), so it's
            forced to 0 and disabled while that radio is selected.
        """
        if not widget.get_active():
            return   # fires for both the button losing AND gaining active state
        is_minimization = self.radio_ensemble_minimization.get_active()
        self.spinbtn_run_steps.set_sensitive(not is_minimization)
        if is_minimization:
            self.spinbtn_run_steps.set_value(0)

    def _get_ensemble_mode(self):
        """ Returns 'minimization', 'nve', 'nvt' or 'npt'. """
        if self.radio_ensemble_minimization.get_active():
            return 'minimization'
        if self.radio_ensemble_nve.get_active():
            return 'nve'
        if self.radio_ensemble_nvt.get_active():
            return 'nvt'
        return 'npt'

    def on_checkbox_restart_single_toggled(self, widget):
        # . Hidden, not just greyed out, when off -- the common case
        #   (most runs don't restart), so keeping 3 file-picker rows
        #   permanently on screen for it cost more space than it was
        #   worth.
        self.grid_restart_single.set_visible(widget.get_active())
        self.grid_restart_single.set_sensitive(widget.get_active())

    def on_button_browse_restart_coor_clicked(self, widget):
        self._browse_file(self.entry_restart_coor, "Select restart coordinates (.coor)", ["*.coor"], "NAMD binary coordinates (*.coor)")

    def on_button_browse_restart_vel_clicked(self, widget):
        self._browse_file(self.entry_restart_vel, "Select restart velocities (.vel)", ["*.vel"], "NAMD binary velocities (*.vel)")

    def on_button_browse_restart_xsc_clicked(self, widget):
        self._browse_file(self.entry_restart_xsc, "Select extended system (.xsc)", ["*.xsc"], "NAMD extended system (*.xsc)")

    #-------------------------------------------------------------------------
    #  S I G N A L S :  C U S T O M   P I P E L I N E
    #-------------------------------------------------------------------------
    def _renumber_pipeline_steps(self):
        for i, row in enumerate(self.liststore_pipeline):
            row[0] = str(i + 1)

    def _clear_pipeline_customization(self, treeiter):
        """ Clears one row's "Custom" override (columns 7/8), if it has
            one. A customized step's saved text has its PREDECESSOR's
            restart-file paths (binCoordinates/binVelocities/
            extendedSystem) and its own outputName baked in from
            whatever position it had at the moment "Setup step..." was
            saved (see on_button_setup_step_clicked()'s docstring) --
            once that row's position actually changes (Move up/down, or
            Remove step shifting everything after it), those baked-in
            paths silently point at the wrong step. Cheaper and safer
            to just drop back to auto-generated (which always reflects
            the CURRENT position) than to try to rewrite the saved text.
            Returns True if there was a customization to clear.
        """
        if not self.liststore_pipeline[treeiter][7]:
            return False
        self.liststore_pipeline[treeiter][7] = ''
        self.liststore_pipeline[treeiter][8] = ''
        return True

    def _clear_customization_after_reorder(self, treeiters):
        """ Runs _clear_pipeline_customization() over every row whose
            position just changed (Move up/down's two swapped rows, or
            every row from the removed index onward after Remove step)
            and reports which step number(s) it actually affected via
            label_info -- non-blocking (no confirmation dialog: this
            can happen often during ordinary pipeline editing), but
            still visible so a customization isn't silently lost
            without the user noticing.
        """
        cleared = []
        for treeiter in treeiters:
            if self._clear_pipeline_customization(treeiter):
                index = self.liststore_pipeline.get_path(treeiter).get_indices()[0]
                cleared.append(index + 1)
        if not cleared:
            return
        cleared.sort()
        numbers = ', '.join(str(n) for n in cleared)
        plural = 's' if len(cleared) > 1 else ''
        self.label_info.set_text(
            'Cleared the "Custom" override on step{} {} (position changed -- a customized '
            "step's text has its predecessor's restart files baked in).".format(plural, numbers))

    def on_button_add_pipeline_step_clicked(self, widget):
        """ Appends a new step built from whatever the ensemble/
            temperature/minimize-or-run-steps/restraints fields above
            are CURRENTLY set to (NOT temperature/timestep/cutoff/
            output-frequency-type fields, which stay shared/global
            across the whole pipeline -- only these five things vary
            per step, matching what actually varies between stages of
            the reference AMBER equilibration protocol this pipeline
            can reproduce).
        """
        ensemble_mode = self._get_ensemble_mode()
        steps = (int(self.spinbtn_minimize_steps.get_value()) if ensemble_mode == 'minimization'
                 else int(self.spinbtn_run_steps.get_value()))
        if steps <= 0:
            self.main.simple_dialog.info(msg='Set at least 1 step (minimize or run, depending on the ensemble) before adding.')
            return

        self.liststore_pipeline.append([
            str(len(self.liststore_pipeline) + 1),
            _ENSEMBLE_DISPLAY[ensemble_mode],
            '{:.1f}'.format(self.spinbtn_temperature.get_value()),
            str(steps),
            _RESTRAINT_DISPLAY[self.combo_restraint_mode.get_active_id()],
            '{:.1f}'.format(self.spinbtn_constraint_scaling.get_value()),
            'Queued', '', '',
        ])

    def on_pipeline_ensemble_edited(self, cell, path, new_text):
        if new_text in _ENSEMBLE_FROM_DISPLAY:
            self.liststore_pipeline[path][1] = new_text

    def on_pipeline_temperature_edited(self, cell, path, new_text):
        try:
            value = float(new_text)
        except ValueError:
            return
        self.liststore_pipeline[path][2] = '{:.1f}'.format(value)

    def on_pipeline_steps_edited(self, cell, path, new_text):
        try:
            value = int(new_text)
        except ValueError:
            return
        if value > 0:
            self.liststore_pipeline[path][3] = str(value)

    def on_pipeline_restraints_edited(self, cell, path, new_text):
        if new_text in _RESTRAINT_FROM_DISPLAY:
            self.liststore_pipeline[path][4] = new_text

    def on_pipeline_scaling_edited(self, cell, path, new_text):
        try:
            value = float(new_text)
        except ValueError:
            return
        self.liststore_pipeline[path][5] = '{:.1f}'.format(value)

    def _get_selected_pipeline_iter(self):
        model, treeiter = self.treeview_pipeline.get_selection().get_selected()
        return treeiter   # None if nothing selected

    def _on_pipeline_treeview_button_press(self, widget, event):
        """ Right-click context menu for the pipeline table -- same
            pattern as process_manager_window.py's own job-history
            treeview. Selects the row under the pointer FIRST (a right-
            click on a row the user hadn't already left-clicked should
            still act on THAT row, not whatever was selected before),
            then pops the menu up there; a right-click on empty space
            leaves the selection alone (Move up/down/Remove/Setup step
            need a selection and just no-op without one -- only "Clear
            pipeline" works regardless).
        """
        if event.type == Gdk.EventType.BUTTON_PRESS and event.button == 3:
            path_info = widget.get_path_at_pos(int(event.x), int(event.y))
            if path_info is not None:
                path, column, cell_x, cell_y = path_info
                widget.grab_focus()
                widget.set_cursor(path, column, 0)
            self._pipeline_popup_menu.popup_at_pointer(event)
            return True
        return False

    def on_button_remove_pipeline_step_clicked(self, widget):
        treeiter = self._get_selected_pipeline_iter()
        if treeiter is None:
            return
        removed_index = self.liststore_pipeline.get_path(treeiter).get_indices()[0]
        self.liststore_pipeline.remove(treeiter)

        # . Every row after the removed one just shifted position by
        #   one -- see _clear_pipeline_customization()'s own docstring.
        shifted_iters = [self.liststore_pipeline.get_iter(i)
                          for i in range(removed_index, len(self.liststore_pipeline))]
        self._renumber_pipeline_steps()
        self._clear_customization_after_reorder(shifted_iters)

    def on_button_clear_pipeline_clicked(self, widget):
        self.liststore_pipeline.clear()

    def on_button_pipeline_move_up_clicked(self, widget):
        treeiter = self._get_selected_pipeline_iter()
        if treeiter is None:
            return
        path = self.liststore_pipeline.get_path(treeiter)
        index = path.get_indices()[0]
        if index == 0:
            return
        previous_iter = self.liststore_pipeline.get_iter(index - 1)
        self.liststore_pipeline.swap(treeiter, previous_iter)
        self._renumber_pipeline_steps()
        self._clear_customization_after_reorder([treeiter, previous_iter])

    def on_button_pipeline_move_down_clicked(self, widget):
        treeiter = self._get_selected_pipeline_iter()
        if treeiter is None:
            return
        path = self.liststore_pipeline.get_path(treeiter)
        index = path.get_indices()[0]
        if index >= len(self.liststore_pipeline) - 1:
            return
        next_iter = self.liststore_pipeline.get_iter(index + 1)
        self.liststore_pipeline.swap(treeiter, next_iter)
        self._renumber_pipeline_steps()
        self._clear_customization_after_reorder([treeiter, next_iter])

    def _confirm_replace_pipeline(self):
        """ True if it's OK to overwrite the current pipeline contents
            -- always True when it's already empty, otherwise asks
            first (same confirmation pattern as on_button_stop_clicked).
        """
        if len(self.liststore_pipeline) == 0:
            return True
        return self.main.simple_dialog.question(
            'Replace the {} step(s) currently in the pipeline?'.format(len(self.liststore_pipeline)))

    def _populate_pipeline_from_rows(self, rows):
        """ Clears the pipeline and fills it from `rows` -- an iterable
            of (ensemble, temperature, steps, restraint_mode, scaling)
            5-tuples, OR (..., custom_text) 6-tuples, using the INTERNAL
            ensemble/restraint ids ('minimization'/'nve'/'nvt'/'npt',
            'none'/'solute'/'backbone'), same shape _PRESETS entries and
            the JSON save/load format use. `custom_text` (the "Setup
            step..." override, see on_button_setup_step_clicked()) is
            always '' for _PRESETS rows -- presets are pure auto-
            generated pipelines -- and only ever non-empty when
            reloading a previously saved pipeline that had one. """
        self.liststore_pipeline.clear()
        for row in rows:
            ensemble, temperature, steps, restraint_mode, scaling = row[:5]
            custom_text = row[5] if len(row) > 5 else ''
            self.liststore_pipeline.append([
                str(len(self.liststore_pipeline) + 1),
                _ENSEMBLE_DISPLAY[ensemble],
                '{:.1f}'.format(temperature),
                str(steps),
                _RESTRAINT_DISPLAY[restraint_mode],
                '{:.1f}'.format(scaling),
                'Queued', custom_text, 'Custom' if custom_text else '',
            ])

    def on_button_load_pipeline_preset_clicked(self, widget):
        name = self.combo_pipeline_presets.get_active_text()
        if name is None or not self._confirm_replace_pipeline():
            return
        self._populate_pipeline_from_rows(_PRESETS[name])

    def on_button_save_pipeline_clicked(self, widget):
        if len(self.liststore_pipeline) == 0:
            self.main.simple_dialog.info(msg='The pipeline is empty -- add at least one step first.')
            return

        dialog = Gtk.FileChooserDialog(
            title="Save pipeline", parent=self.window, action=Gtk.FileChooserAction.SAVE,
        )
        dialog.set_do_overwrite_confirmation(True)
        dialog.set_current_name('pipeline.json')
        dialog.add_buttons(
            Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
            Gtk.STOCK_SAVE, Gtk.ResponseType.OK,
        )
        if dialog.run() == Gtk.ResponseType.OK:
            path = dialog.get_filename()
            if not path.endswith('.json'):
                path += '.json'
            data = {
                'pipeline': [
                    {
                        'ensemble': _ENSEMBLE_FROM_DISPLAY[row[1]],
                        'temperature': float(row[2]),
                        'steps': int(row[3]),
                        'restraints': _RESTRAINT_FROM_DISPLAY[row[4]],
                        'scaling': float(row[5]),
                        'custom_text': row[7],
                    }
                    for row in self.liststore_pipeline
                ],
            }
            try:
                with open(path, 'w') as f:
                    json.dump(data, f, indent=2)
            except Exception as error:
                traceback.print_exc()
                self.main.simple_dialog.info(msg='Could not save the pipeline:\n{}'.format(error))
        dialog.destroy()

    def on_button_load_pipeline_file_clicked(self, widget):
        dialog = Gtk.FileChooserDialog(
            title="Load pipeline", parent=self.window, action=Gtk.FileChooserAction.OPEN,
        )
        dialog.add_buttons(
            Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
            Gtk.STOCK_OPEN, Gtk.ResponseType.OK,
        )
        file_filter = Gtk.FileFilter()
        file_filter.set_name("Pipeline files (*.json)")
        file_filter.add_pattern("*.json")
        dialog.add_filter(file_filter)

        path = None
        if dialog.run() == Gtk.ResponseType.OK:
            path = dialog.get_filename()
        dialog.destroy()
        if path is None:
            return

        try:
            with open(path) as f:
                data = json.load(f)
            rows = [
                (step['ensemble'], float(step['temperature']), int(step['steps']),
                 step['restraints'], float(step['scaling']), step.get('custom_text', ''))
                for step in data['pipeline']
            ]
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg='Could not load "{}":\n{}'.format(path, error))
            return

        if not self._confirm_replace_pipeline():
            return
        self._populate_pipeline_from_rows(rows)

    #-------------------------------------------------------------------------
    #  S I G N A L S :  A D V A N C E D   T A B
    #-------------------------------------------------------------------------
    def on_checkbox_pbc_toggled(self, widget):
        """ "Periodic Cell" (explicit box + PME) vs "Implicit Solvent"
            (NAMD's GBIS, no periodic cell at all) -- mutually exclusive,
            see namd_runner.build_namd_config()'s own `pbc` docstring.
        """
        is_pbc = widget.get_active()
        self.frame_periodic_cell.set_visible(is_pbc)
        self.frame_implicit_solvent.set_visible(not is_pbc)

    def on_checkbox_sasa_toggled(self, widget):
        # . Hidden, not just greyed out, when off -- surfaceTension
        #   means nothing without SASA on (same compacting pattern as
        #   switchdist/the restraint-scaling field elsewhere).
        self.spinbtn_surface_tension.set_visible(widget.get_active())

    def on_checkbox_implicit_solvent_toggled(self, widget):
        # . Greyed out, not hidden, when off -- same pattern as
        #   checkbox_manual_cell/grid_manual_cell: the fields become
        #   irrelevant (a plain-vacuum run with no solvent model at
        #   all), not gone, since this checkbox IS this frame's own
        #   title/label (only ever visible to begin with when PBC is
        #   off -- see on_checkbox_pbc_toggled()).
        self.grid_implicit_solvent.set_sensitive(widget.get_active())

    def _get_implicit_solvent_kwargs(self):
        """ Reads the "Implicit Solvent" section into the kwargs
            namd_runner.build_namd_config() expects -- only meaningful
            (and only ever passed by the callers below) when PBC is
            off; see that function's own docstring for what each one
            means and which NAMD keyword it becomes.
        """
        return {
            'implicit_solvent': self.checkbox_implicit_solvent.get_active(),
            'solvent_dielectric': self.spinbtn_solvent_dielectric.get_value(),
            'ion_concentration': self.spinbtn_ion_concentration.get_value(),
            'alpha_cutoff': self.spinbtn_alpha_cutoff.get_value(),
            'sasa': self.checkbox_sasa.get_active(),
            'surface_tension': self.spinbtn_surface_tension.get_value(),
        }

    def on_checkbox_manual_cell_toggled(self, widget):
        self.grid_manual_cell.set_sensitive(widget.get_active())

    def on_checkbox_override_langevin_temp_toggled(self, widget):
        self.spinbtn_langevin_temp.set_sensitive(widget.get_active())

    def _get_cell_override(self):
        """ Returns a cell dict (same shape as
            namd_runner.compute_cell_from_coordinates()'s return value)
            built from the Advanced tab's manual fields, or None if
            "Manually set periodic cell" is unchecked (the normal case
            -- the cell is then computed automatically from the input
            coordinates, or inherited from a previous step's .xsc).
        """
        if not self.checkbox_manual_cell.get_active():
            return None
        return {
            'cell_basis_vector1': (self.spinbtn_cell_a.get_value(), 0.0, 0.0),
            'cell_basis_vector2': (0.0, self.spinbtn_cell_b.get_value(), 0.0),
            'cell_basis_vector3': (0.0, 0.0, self.spinbtn_cell_c.get_value()),
            'cell_origin': (self.spinbtn_cell_origin_x.get_value(),
                             self.spinbtn_cell_origin_y.get_value(),
                             self.spinbtn_cell_origin_z.get_value()),
            'pme_grid_size': (int(self.spinbtn_pme_grid_x.get_value()),
                               int(self.spinbtn_pme_grid_y.get_value()),
                               int(self.spinbtn_pme_grid_z.get_value())),
        }

    def _get_output_freq_kwargs(self):
        """ restartfreq/dcdfreq/outputEnergies/outputPressure, read from
            the Advanced tab -- shared by both Single Run and the
            equilibration protocol (every step of the latter uses the
            same values, matching build_equilibration_protocol()'s own
            "forwarded to every step" contract). """
        return {
            'restartfreq': int(self.spinbtn_restartfreq.get_value()),
            'dcdfreq': int(self.spinbtn_dcdfreq.get_value()),
            'outputenergies': int(self.spinbtn_outputenergies.get_value()),
            'outputpressure': int(self.spinbtn_outputpressure.get_value()),
        }

    def _get_langevin_temp_override(self):
        if self.checkbox_override_langevin_temp.get_active():
            return self.spinbtn_langevin_temp.get_value()
        return None

    def on_checkbox_switching_toggled(self, widget):
        # . Hidden, not just greyed out, when off -- switchdist means
        #   nothing without switching on (matches the compacting
        #   pattern the "Manually set periodic cell"/restraints fields
        #   already use elsewhere in this window).
        self.label_switchdist.set_visible(widget.get_active())
        self.spinbtn_switchdist.set_visible(widget.get_active())

    def _get_force_field_kwargs(self):
        """ Reads the "Force Field" section (Setup tab) into the kwargs
            namd_runner.build_namd_config() expects for the values NAMD's
            own AMBER-support docs (ug/node13.html) call out as
            characterizing AMBER-vs-CHARMM usage -- see that function's
            own docstring for exactly what each keyword means and which
            NAMD keyword name it ends up under. readexclusions/scnb are
            read regardless of AMBER-vs-CHARMM here (harmless -- build_
            namd_config() itself already omits both entirely in CHARMM
            mode, see its own `if not is_charmm` guard); this window's
            own on_input_type_radio_toggled() just greys the two widgets
            out when CHARMM is selected, as a hint they won't apply.
        """
        return {
            'exclude': self.combo_exclude.get_active_id(),
            'switching': self.checkbox_switching.get_active(),
            'switchdist': self.spinbtn_switchdist.get_value(),
            'one_four_scaling': self.spinbtn_one_four_scaling.get_value(),
            'readexclusions': self.checkbox_readexclusions.get_active(),
            'scnb': self.spinbtn_scnb.get_value(),
            'rigidbonds': self.combo_rigidbonds.get_active_id(),
            'rigidtolerance': self.spinbtn_rigidtolerance.get_value(),
            'vdw_force_switching': self.checkbox_vdw_force_switching.get_active(),
        }

    #-------------------------------------------------------------------------
    #  S I G N A L S :  B U T T O N S
    #-------------------------------------------------------------------------
    def on_button_cancel_clicked(self, widget):
        """ Function doc """
        self.close_window()

    #-------------------------------------------------------------------------
    #  C O M M O N   V A L I D A T I O N
    #-------------------------------------------------------------------------
    def _validate_common_inputs(self):
        """ Reads/validates the fields shared by both tabs REGARDLESS of
            input type (executable/working folder/system name -- NOT
            the AMBER-vs-CHARMM structure files themselves, see
            _validate_input_files() for those). Returns (namd_command,
            work_folder, system_name) or None (and shows a dialog) if
            something is missing/invalid.
        """
        if self._job_running():
            self.main.simple_dialog.info(msg='A job is already running -- stop it first.')
            return None

        namd_command = self.entry_namd_command.get_text().strip()
        if not namd_command or not (os.path.isfile(namd_command) and os.access(namd_command, os.X_OK)):
            self.main.simple_dialog.info(msg='Please provide a valid path to the namd3/namd2 executable.')
            return None
        self.main.vm_session.vm_config.gl_parameters['namd_command'] = namd_command

        work_folder = self.work_folder_chooser.get_folder()
        if not work_folder:
            self.main.simple_dialog.info(msg='Please choose a working folder.')
            return None

        system_name = self.entry_system_name.get_text().strip() or 'namd_run'

        return namd_command, work_folder, system_name

    def _get_input_type(self):
        return 'charmm' if self.radio_input_charmm.get_active() else 'amber'

    def _get_selected_vobject(self, system_combo, object_combo):
        """ Reads a vobject out of a SystemComboBox/CoordinatesComboBox
            pair (same widgets/API prepare_amber_system.py already
            uses) -- returns the VismolObject, or None (nothing valid
            selected). """
        vobject_id = object_combo.get_vobject_id()
        if vobject_id is None:
            return None
        return self.main.vm_session.vm_objects_dic.get(vobject_id)

    def _write_amber_crd_from_vobject(self, vobject, output_path):
        """ Writes the SELECTED vobject's current frame as a plain
            AMBER coordinate file (see namd_runner.write_amber_crd()'s
            own docstring for the format/verification). Atom order
            follows vobject.atoms' own iteration order (index-sorted,
            same convention export_special_PDB() already relies on
            elsewhere in this codebase) -- this MUST match the
            accompanying prmtop's atom order for the run to be
            physically meaningful (same requirement any restart file
            already has; not independently re-validated here, just
            like a manually-picked restart .coor isn't either).
        """
        positions = [vobject.atoms[index].coords(-1) for index in vobject.atoms.keys()]
        namd_runner.write_amber_crd(positions, output_path, title=vobject.name)

    def _validate_input_files(self, work_folder):
        """ Reads/validates the structure-file section, branching on
            the "Input type" radio AND, independently, each mode's own
            "Coordinates source" (file vs a vobject already loaded in
            EasyHybrid). When the source is a vobject, this WRITES a
            fresh coordinates file into `work_folder` from its current
            frame (an AMBER .crd via write_amber_crd(), or a plain PDB
            via the already-existing export_special_PDB() for CHARMM)
            and returns that path exactly like a file the user had
            browsed to -- everything downstream (cell computation,
            restraint PDB, build_namd_config, ...) stays unaware of
            where the coordinates actually came from.

            Returns a dict:
              {'type': 'amber', 'prmtop': ..., 'inpcrd': ...}
            or
              {'type': 'charmm', 'psf': ..., 'parameters': [...],
               'coordinates': ..., 'xplor': bool}
            or None (and shows a dialog) if something is missing/invalid.
        """
        if self._get_input_type() == 'charmm':
            psf = self.entry_charmm_psf_path.get_text().strip()
            parameters = list(self.charmm_parameters)
            if not psf or not os.path.isfile(psf):
                self.main.simple_dialog.info(msg='Please select a valid CHARMM psf file.')
                return None
            if not parameters:
                self.main.simple_dialog.info(msg='Please add at least one CHARMM parameter file.')
                return None

            if self.combo_charmm_coords_source.get_active_id() == 'vobject':
                vobject = self._get_selected_vobject(self.charmm_coords_system_combo, self.charmm_coords_object_combo)
                if vobject is None:
                    self.main.simple_dialog.info(msg='Please select a system/object for the coordinates.')
                    return None
                coordinates = os.path.join(work_folder, '_namd_coords_from_vobject.pdb')
                try:
                    export_special_PDB(vobject=vobject, frame=-1, output=coordinates)
                except Exception as error:
                    traceback.print_exc()
                    self.main.simple_dialog.info(msg='Could not export coordinates from the selected object:\n{}'.format(error))
                    return None
            else:
                coordinates = self.entry_charmm_coordinates_path.get_text().strip()
                if not coordinates or not os.path.isfile(coordinates):
                    self.main.simple_dialog.info(msg='Please select a valid coordinates (pdb) file.')
                    return None

            return {
                'type': 'charmm', 'psf': psf, 'parameters': parameters,
                'coordinates': coordinates, 'xplor': self.checkbox_charmm_xplor.get_active(),
            }

        prmtop = self.entry_prmtop_path.get_text().strip()
        if not prmtop or not os.path.isfile(prmtop):
            self.main.simple_dialog.info(msg='Please select a valid AMBER prmtop file.')
            return None

        if self.combo_amber_coords_source.get_active_id() == 'vobject':
            vobject = self._get_selected_vobject(self.amber_coords_system_combo, self.amber_coords_object_combo)
            if vobject is None:
                self.main.simple_dialog.info(msg='Please select a system/object for the coordinates.')
                return None
            inpcrd = os.path.join(work_folder, '_namd_coords_from_vobject.crd')
            try:
                self._write_amber_crd_from_vobject(vobject, inpcrd)
            except Exception as error:
                traceback.print_exc()
                self.main.simple_dialog.info(msg='Could not export coordinates from the selected object:\n{}'.format(error))
                return None
        else:
            inpcrd = self.entry_inpcrd_path.get_text().strip()
            if not inpcrd or not os.path.isfile(inpcrd):
                self.main.simple_dialog.info(msg='Please select a valid AMBER inpcrd file.')
                return None

        return {'type': 'amber', 'prmtop': prmtop, 'inpcrd': inpcrd}

    #-------------------------------------------------------------------------
    #  R U N :  S I N G L E
    #-------------------------------------------------------------------------
    def on_button_run_single_clicked(self, widget):
        """ Function doc """
        validated = self._validate_common_inputs()
        if validated is None:
            return
        namd_command, work_folder, system_name = validated

        input_files = self._validate_input_files(work_folder)
        if input_files is None:
            return
        is_charmm = input_files['type'] == 'charmm'

        restraint_mode = self.combo_restraint_mode.get_active_id()
        restraints = None
        if restraint_mode != 'none':
            try:
                restraint_pdb = os.path.join(work_folder, system_name + '_const.pdb')
                if is_charmm:
                    namd_runner.build_restraint_pdb_from_charmm(input_files['coordinates'], restraint_pdb)
                else:
                    namd_runner.build_restraint_pdb(input_files['prmtop'], input_files['inpcrd'], restraint_pdb)
            except Exception as error:
                traceback.print_exc()
                self.main.simple_dialog.info(msg='Could not build the restraint PDB:\n{}'.format(error))
                return
            restraints = {
                'pdb': restraint_pdb,
                'column': 'O' if restraint_mode == 'solute' else 'B',
                'scaling': self.spinbtn_constraint_scaling.get_value(),
            }

        bin_coordinates = bin_velocities = extended_system = None
        if self.checkbox_restart_single.get_active():
            bin_coordinates = self.entry_restart_coor.get_text().strip() or None
            bin_velocities = self.entry_restart_vel.get_text().strip() or None
            extended_system = self.entry_restart_xsc.get_text().strip() or None

        pbc = self.checkbox_pbc.get_active()
        cell = None
        if pbc:
            cell = self._get_cell_override()
            if cell is None:
                try:
                    if is_charmm:
                        cell = namd_runner.compute_cell_from_charmm(input_files['coordinates'])
                    else:
                        cell = namd_runner.compute_cell_from_coordinates(input_files['prmtop'], input_files['inpcrd'])
                except Exception as error:
                    traceback.print_exc()
                    self.main.simple_dialog.info(msg='Could not compute the periodic cell from the input coordinates:\n{}'.format(error))
                    return

        minimize_steps = int(self.spinbtn_minimize_steps.get_value())
        run_steps = int(self.spinbtn_run_steps.get_value())

        ensemble_mode = self._get_ensemble_mode()
        if ensemble_mode == 'minimization':
            run_steps = 0
            if minimize_steps <= 0:
                self.main.simple_dialog.info(msg='"Minimization" needs at least 1 minimize step.')
                return
        # . "Minimization" and "NVE" both mean NO Langevin thermostat --
        #   the fine-tune fields in the Advanced tab (damping/temp
        #   override/hydrogen) simply go unused in that case (see
        #   build_namd_config()'s own docstring: they're skipped
        #   entirely when langevin=False).
        langevin = ensemble_mode in ('nvt', 'npt')
        namd_ensemble = 'NPT' if ensemble_mode == 'npt' else 'NVT'   # only "NPT" turns the piston on -- see build_namd_config()'s docstring

        input_kwargs = (
            {
                'charmm_psf': input_files['psf'], 'charmm_parameters': input_files['parameters'],
                'charmm_xplor': input_files['xplor'], 'coordinates': input_files['coordinates'],
            } if is_charmm else
            {'parmfile': input_files['prmtop'], 'ambercoor': input_files['inpcrd']}
        )

        config_text = namd_runner.build_namd_config(
            output_basename=system_name, **input_kwargs,
            bin_coordinates=bin_coordinates, bin_velocities=bin_velocities,
            extended_system=extended_system, cell=cell,
            temperature=self.spinbtn_temperature.get_value(),
            timestep=self.spinbtn_timestep.get_value(),
            cutoff=self.spinbtn_cutoff.get_value(),
            watermodel=self.combo_watermodel.get_active_id() or 'tip3',
            minimize_steps=minimize_steps, run_steps=run_steps,
            ensemble=namd_ensemble,
            restraints=restraints,
            langevin=langevin,
            langevin_damping=self.spinbtn_langevin_damping.get_value(),
            langevin_temp=self._get_langevin_temp_override(),
            langevin_hydrogen=self.checkbox_langevin_hydrogen.get_active(),
            flexible_cell=self.checkbox_flexible_cell.get_active(),
            pbc=pbc,
            wrap_all=self.checkbox_wrap_all.get_active(),
            **self._get_implicit_solvent_kwargs(),
            **self._get_force_field_kwargs(),
            **self._get_output_freq_kwargs()
        )

        config_path = os.path.join(work_folder, system_name + '.namd')
        with open(config_path, 'w') as f:
            f.write(config_text)
        log_path = os.path.join(work_folder, system_name + '.log')

        self._run_mode = 'single'
        self._protocol_steps = None
        self._input_files = input_files
        self._protocol_system_name = system_name
        self._protocol_work_folder = work_folder

        self._start_process(config_path, work_folder, namd_command, log_path,
                             int(self.spinbtn_n_procs.get_value()))
        self.label_info.set_text('Running NAMD ("{}")...'.format(system_name))

    #-------------------------------------------------------------------------
    #  R U N :  C U S T O M   P I P E L I N E
    #-------------------------------------------------------------------------
    def _build_pipeline_plan(self):
        """ Validates everything and builds every pipeline step's config
            text FRESH from the treeview + Single Run/Advanced tabs --
            shared by on_button_run_pipeline_clicked() (which actually
            runs the resulting chain) and on_button_setup_step_clicked()
            (which only needs one step's up-to-date auto-generated text
            to preview/edit, before anything has run). Each row restarts
            from the previous one's own coordinates/velocities/cell
            (binCoordinates/binVelocities/extendedSystem), exactly the
            mechanism the real 9-step AMBER protocol this window's
            "AMBER Equilibration Protocol" preset is based on uses (see
            _PRESETS below). Works with either AMBER or CHARMM input --
            nothing about chaining restart files is AMBER-specific.

            Every returned step dict's 'config_text' is the freshly
            auto-generated text REGARDLESS of whether that row also has
            a saved 'custom_text' override -- callers decide which one
            to actually use (_launch_protocol_step() prefers
            'custom_text' when non-empty; the Setup editor's "Reset to
            auto-generated" needs the fresh one specifically).

            Returns None (having already shown a dialog) on any
            validation failure -- same contract as
            _validate_common_inputs()/_validate_input_files().
        """
        validated = self._validate_common_inputs()
        if validated is None:
            return None
        namd_command, work_folder, system_name = validated

        input_files = self._validate_input_files(work_folder)
        if input_files is None:
            return None

        if len(self.liststore_pipeline) == 0:
            self.main.simple_dialog.info(msg='Add at least one step to the pipeline first.')
            return None
        is_charmm = input_files['type'] == 'charmm'

        # . Parse every row up front, so a bad cell value is caught
        #   before anything is built/run, not mid-chain.
        rows = [
            {
                'ensemble': _ENSEMBLE_FROM_DISPLAY[row[1]],
                'temperature': float(row[2]),
                'steps': int(row[3]),
                'restraint_mode': _RESTRAINT_FROM_DISPLAY[row[4]],
                'scaling': float(row[5]),
                'custom_text': row[7],
            }
            for row in self.liststore_pipeline
        ]

        restraint_pdb = None
        if any(r['restraint_mode'] != 'none' for r in rows):
            try:
                restraint_pdb = os.path.join(work_folder, system_name + '_const.pdb')
                if is_charmm:
                    namd_runner.build_restraint_pdb_from_charmm(input_files['coordinates'], restraint_pdb)
                else:
                    namd_runner.build_restraint_pdb(input_files['prmtop'], input_files['inpcrd'], restraint_pdb)
            except Exception as error:
                traceback.print_exc()
                self.main.simple_dialog.info(msg='Could not build the restraint PDB:\n{}'.format(error))
                return None

        pbc = self.checkbox_pbc.get_active()
        cell = None
        if pbc:
            cell = self._get_cell_override()
            if cell is None:
                try:
                    if is_charmm:
                        cell = namd_runner.compute_cell_from_charmm(input_files['coordinates'])
                    else:
                        cell = namd_runner.compute_cell_from_coordinates(input_files['prmtop'], input_files['inpcrd'])
                except Exception as error:
                    traceback.print_exc()
                    self.main.simple_dialog.info(msg='Could not compute the periodic cell from the input coordinates:\n{}'.format(error))
                    return None

        input_kwargs = (
            {
                'charmm_psf': input_files['psf'], 'charmm_parameters': input_files['parameters'],
                'charmm_xplor': input_files['xplor'], 'coordinates': input_files['coordinates'],
            } if is_charmm else
            {'parmfile': input_files['prmtop'], 'ambercoor': input_files['inpcrd']}
        )
        # . Shared across every step -- only ensemble/temperature/steps/
        #   restraints (parsed per row above) vary between them, same as
        #   what actually varies between stages of the reference
        #   protocol (see build_equilibration_protocol()'s docstring).
        timestep = self.spinbtn_timestep.get_value()
        cutoff = self.spinbtn_cutoff.get_value()
        watermodel = self.combo_watermodel.get_active_id() or 'tip3'
        langevin_damping = self.spinbtn_langevin_damping.get_value()
        langevin_temp = self._get_langevin_temp_override()
        langevin_hydrogen = self.checkbox_langevin_hydrogen.get_active()
        flexible_cell = self.checkbox_flexible_cell.get_active()
        wrap_all = self.checkbox_wrap_all.get_active()
        implicit_solvent_kwargs = self._get_implicit_solvent_kwargs()
        force_field_kwargs = self._get_force_field_kwargs()
        output_freq_kwargs = self._get_output_freq_kwargs()

        steps_out = []
        previous_basename = None
        for i, r in enumerate(rows):
            is_minimization = r['ensemble'] == 'minimization'
            restraints = None
            if r['restraint_mode'] != 'none':
                restraints = {
                    'pdb': restraint_pdb,
                    'column': 'O' if r['restraint_mode'] == 'solute' else 'B',
                    'scaling': r['scaling'],
                }

            if previous_basename is None:
                bin_coordinates = bin_velocities = extended_system = None
            else:
                bin_coordinates = os.path.join(work_folder, previous_basename + '.coor')
                bin_velocities = os.path.join(work_folder, previous_basename + '.vel')
                extended_system = os.path.join(work_folder, previous_basename + '.xsc')
            # . `cell` is passed to EVERY step (not just the first) --
            #   PMEGridSizeX/Y/Z has to stay fixed for the whole chain
            #   (see build_equilibration_protocol()'s own docstring for
            #   why); build_namd_config() itself already suppresses the
            #   cellBasisVector*/cellOrigin lines once extended_system is
            #   set (a restart step reads its cell from the previous
            #   step's own .xsc instead), so passing the same `cell`
            #   here for every step is correct, not redundant. wrapAll
            #   is NOT suppressed for restart steps -- every step with a
            #   periodic cell gets its own trajectory wrapped.

            output_basename = '{}_step{:02d}'.format(system_name, i + 1)
            try:
                config_text = namd_runner.build_namd_config(
                    output_basename=output_basename, **input_kwargs,
                    bin_coordinates=bin_coordinates, bin_velocities=bin_velocities,
                    extended_system=extended_system, cell=cell,
                    temperature=r['temperature'], timestep=timestep,
                    cutoff=cutoff, watermodel=watermodel,
                    minimize_steps=r['steps'] if is_minimization else 0,
                    run_steps=0 if is_minimization else r['steps'],
                    ensemble='NPT' if r['ensemble'] == 'npt' else 'NVT',
                    restraints=restraints,
                    langevin=r['ensemble'] in ('nvt', 'npt'),
                    langevin_damping=langevin_damping, langevin_temp=langevin_temp,
                    langevin_hydrogen=langevin_hydrogen,
                    flexible_cell=flexible_cell,
                    pbc=pbc,
                    wrap_all=wrap_all,
                    **implicit_solvent_kwargs,
                    **force_field_kwargs,
                    **output_freq_kwargs
                )
            except Exception as error:
                traceback.print_exc()
                self.main.simple_dialog.info(msg='Could not build step {}:\n{}'.format(i + 1, error))
                return None

            steps_out.append({
                'description': 'Step {} ({})'.format(i + 1, _ENSEMBLE_DISPLAY[r['ensemble']]),
                'config_text': config_text,
                'custom_text': r['custom_text'],
                'config_path': os.path.join(work_folder, output_basename + '.namd'),
                'output_basename': output_basename,
            })
            previous_basename = output_basename

        return {
            'namd_command': namd_command, 'work_folder': work_folder,
            'system_name': system_name, 'input_files': input_files, 'steps': steps_out,
        }

    def on_button_run_pipeline_clicked(self, widget):
        """ Runs the chain built by _build_pipeline_plan() -- see that
            method's own docstring for how each step is put together.
        """
        plan = self._build_pipeline_plan()
        if plan is None:
            return

        for row in self.liststore_pipeline:
            row[6] = 'Queued'

        self._run_mode = 'pipeline'
        self._protocol_steps = plan['steps']
        self._protocol_index = 0
        self._input_files = plan['input_files']
        self._protocol_system_name = plan['system_name']
        self._protocol_namd_command = plan['namd_command']
        self._protocol_n_procs = int(self.spinbtn_n_procs.get_value())
        self._protocol_work_folder = plan['work_folder']

        self._launch_protocol_step()

    def _step_count(self):
        return len(self._protocol_steps)

    def _set_step_status(self, index, text):
        treeiter = self.liststore_pipeline.get_iter(index)
        self.liststore_pipeline[treeiter][6] = text

    #-------------------------------------------------------------------------
    #  S I G N A L S :  S T E P   S E T U P   E D I T O R
    #-------------------------------------------------------------------------
    def on_button_setup_step_clicked(self, widget):
        """ Opens the selected pipeline row's own .namd text for hand-
            editing, highlighted green ("safe", _NAMD_SAFE_KEYWORDS) /
            orange-bold ("critical", _NAMD_CRITICAL_KEYWORDS) -- see
            those two regexes' own comment for the split. Rebuilds the
            WHOLE pipeline plan first (_build_pipeline_plan() does the
            same validation a real Run would) purely to get THIS one
            step's up-to-date auto-generated text -- cheap (a dozen
            short strings) and guarantees the preview always matches
            what an actual Run would produce right now, not a stale
            snapshot from whenever the row was last built/run.
        """
        treeiter = self._get_selected_pipeline_iter()
        if treeiter is None:
            self.main.simple_dialog.info(msg='Select a pipeline step first.')
            return

        plan = self._build_pipeline_plan()
        if plan is None:
            return

        path = self.liststore_pipeline.get_path(treeiter)
        index = path.get_indices()[0]

        self._step_setup_row_index = index
        self._step_setup_auto_text = plan['steps'][index]['config_text']
        custom_text = self.liststore_pipeline[treeiter][7]
        self.label_step_setup_title.set_text('Step {} -- edit .namd text'.format(index + 1))
        self._highlight_step_setup_text(custom_text or self._step_setup_auto_text)
        self.window_step_setup.show_all()

    def on_button_step_setup_save_clicked(self, widget):
        buf = self.textview_step_setup.get_buffer()
        text = buf.get_text(buf.get_start_iter(), buf.get_end_iter(), True)
        treeiter = self.liststore_pipeline.get_iter(self._step_setup_row_index)
        self.liststore_pipeline[treeiter][7] = text
        self.liststore_pipeline[treeiter][8] = 'Custom'
        self.window_step_setup.hide()

    def on_button_step_setup_reset_clicked(self, widget):
        """ Discards the saved override (if any) and refills the editor
            with a fresh auto-generated text -- does NOT close the
            window, so the user can keep tweaking from this clean
            baseline; if they just close/Cancel from here, the row
            stays non-custom (built automatically like any other row).
        """
        treeiter = self.liststore_pipeline.get_iter(self._step_setup_row_index)
        self.liststore_pipeline[treeiter][7] = ''
        self.liststore_pipeline[treeiter][8] = ''
        self._highlight_step_setup_text(self._step_setup_auto_text)

    def on_button_step_setup_cancel_clicked(self, widget):
        self.window_step_setup.hide()

    def on_step_setup_delete_event(self, widget, event):
        self.window_step_setup.hide()
        return True   # . Hide, don't destroy -- same window instance is reused next time.

    def _launch_protocol_step(self):
        step = self._protocol_steps[self._protocol_index]
        with open(step['config_path'], 'w') as f:
            f.write(step['custom_text'] or step['config_text'])
        log_path = os.path.join(self._protocol_work_folder, '{:02d}.log'.format(self._protocol_index + 1))

        self._set_step_status(self._protocol_index, 'Running...')
        self.label_info.set_text('Step {}/{}: {}'.format(
            self._protocol_index + 1, self._step_count(), step['description']))

        self._start_process(step['config_path'], self._protocol_work_folder,
                             self._protocol_namd_command, log_path, self._protocol_n_procs)

    #-------------------------------------------------------------------------
    #  P R O C E S S   L I F E C Y C L E   ( S H A R E D )
    #-------------------------------------------------------------------------
    def _start_process(self, config_path, work_folder, namd_command, log_path, n_procs):
        self.current_process = namd_runner.run_namd(config_path, work_folder, namd_command, log_path, n_procs=n_procs)
        self.current_log_path = log_path
        self._log_read_offset = 0
        buf = self.textview_log.get_buffer()
        buf.set_text('')
        self._set_controls_sensitive(False)
        if self._poll_timeout_id is None:
            self._poll_timeout_id = GLib.timeout_add(500, self._poll_job)

    def _tail_log(self):
        if not self.current_log_path or not os.path.isfile(self.current_log_path):
            return
        with open(self.current_log_path, 'r', errors='replace') as f:
            f.seek(self._log_read_offset)
            new_text = f.read()
            self._log_read_offset = f.tell()
        if new_text:
            self._append_log(new_text)

    def _poll_job(self):
        """ GLib.timeout_add callback (every 500ms) -- same non-blocking
            polling pattern process_manager_window.py already uses for
            job tracking. Returns True to keep polling, False to stop.
        """
        self._tail_log()

        if self.current_process is None:
            return False

        returncode = self.current_process.poll()
        if returncode is None:
            return True   # still running

        ok = (returncode == 0) and namd_runner.check_namd_log_success(self.current_log_path)
        self.current_process = None

        if self._run_mode == 'single':
            self._poll_timeout_id = None
            self._set_controls_sensitive(True)
            if not ok:
                self.label_info.set_text('NAMD failed -- see the log above.')
                return False
            self.label_info.set_text('NAMD finished -- importing the result back into EasyHybrid...')
            self._import_result_back(self._input_files,
                                      [self._protocol_work_folder_dcd_single()], self._protocol_system_name)
            return False

        # . 'pipeline' mode -- mark this step, then either advance to
        #   the next one or finish/abort the whole chain.
        step_count = self._step_count()
        self._set_step_status(self._protocol_index, 'Finished' if ok else 'Failed')
        if not ok:
            self._poll_timeout_id = None
            self._set_controls_sensitive(True)
            for i in range(self._protocol_index + 1, step_count):
                self._set_step_status(i, 'Skipped')
            self.label_info.set_text(
                'Step {}/{} failed -- see the log above. The remaining steps were not run.'.format(
                    self._protocol_index + 1, step_count))
            return False

        if self._protocol_index < step_count - 1:
            self._protocol_index += 1
            self._launch_protocol_step()
            return True   # same timeout keeps polling the next step

        # . Every step finished successfully.
        self._poll_timeout_id = None
        self._set_controls_sensitive(True)
        self.label_info.set_text('Pipeline finished -- importing the result back into EasyHybrid...')
        # . Every step's own .dcd, in order -- stitched into ONE continuous
        #   trajectory by _import_result_back() below, not just the last
        #   step's (which is all that used to be imported: each step's
        #   file is genuinely separate on disk -- see build_namd_config()'s
        #   per-step `output_basename` -- but only the final one ever made
        #   it into EasyHybrid/vismol, which is what actually looked like
        #   "each step overwrites the previous one").
        dcd_paths = [os.path.join(self._protocol_work_folder, step['output_basename'] + '.dcd')
                     for step in self._protocol_steps]
        self._import_result_back(self._input_files, dcd_paths, self._protocol_system_name)
        return False

    def _protocol_work_folder_dcd_single(self):
        """ Single-run mode's own .dcd path (built the same way
            on_button_run_single_clicked() names its outputName). """
        return os.path.join(self._protocol_work_folder, self._protocol_system_name + '.dcd')

    #-------------------------------------------------------------------------
    #  I M P O R T   B A C K
    #-------------------------------------------------------------------------
    def _import_result_back(self, input_files, dcd_paths, system_name):
        """ Imports the ORIGINAL topology/coordinates back into EasyHybrid
            (same path the tleap window already uses for AMBER --
            system_type 0; system_type 1 for CHARMM, mirroring
            session.py's own CHARMM branch), then loads the produced
            trajectory/trajectories (.dcd) into the new object via
            pDynamo3's own DCDTrajectoryFileReader (already used by the
            "Import Trajectory" window -- see import_trajectory.py's
            _import_dcd_file()/import_data()) so the whole simulated
            trajectory becomes playable in vismol, not just the final
            frame.

            `dcd_paths` is a list -- a single-run job passes a
            one-element list, a pipeline job passes EVERY step's own
            .dcd IN ORDER, imported one after another into the SAME
            vobject (the first import_data() call creates it via
            new_vobj_name; every next one instead passes that same
            vobject back in and omits new_vobj_name, which
            import_trajectory.py's _import_dcd_file() reads as "append
            these frames, don't create a new object" -- see its own
            branch on `parameters['new_vobj_name']`). The result is ONE
            continuous trajectory spanning the whole chain, matching
            what the run actually did (each step restarts from the
            previous one's own binCoordinates/binVelocities/xsc) --
            previously only the LAST step's .dcd was ever imported,
            which is what made the whole pipeline look like each step
            was overwriting the one before it (the files were always
            separate on disk -- see build_namd_config()'s per-step
            `output_basename` -- they just never all made it into
            EasyHybrid). A step whose run was too short to write even
            one frame at its own dcdfreq is silently skipped here (its
            .dcd simply doesn't exist) -- the resulting trajectory has a
            gap there, but the chain itself was still built and run
            correctly; only the on-screen trajectory is missing that
            segment.

            NOTE (CHARMM only): unlike NAMD itself (which parses psf/
            parameters/coordinates with its own, unaffected reader),
            this step DOES need pDynamo3 to fully resolve the CHARMM
            parameter set -- pDynamo3's own CHARMMParameterFileReader
            was confirmed (against real files) to silently parse ".str"
            CHARMM "stream" parameter files as empty, which then makes
            this specific step fail with a KeyError even though NAMD's
            own run already succeeded. That failure is caught below and
            reported clearly rather than crashing -- the .namd config,
            log and output files are still on disk regardless.
        """
        try:
            if input_files['type'] == 'charmm':
                self.p_session.load_a_new_pDynamo_system_from_dict(
                    input_files={
                        'charmm_psf': input_files['psf'], 'charmm_par': input_files['parameters'],
                        'coordinates': input_files['coordinates'],
                    },
                    system_type=1, name=system_name, tag='NAMD',
                    working_folder=self._protocol_work_folder,
                )
            else:
                top_path, crd_path = namd_runner.ensure_amber_top_crd(
                    input_files['prmtop'], input_files['inpcrd'], self._protocol_work_folder)
                self.p_session.load_a_new_pDynamo_system_from_dict(
                    input_files={'amber_prmtop': top_path, 'coordinates': crd_path},
                    system_type=0, name=system_name, tag='NAMD',
                    working_folder=self._protocol_work_folder,
                )
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg='NAMD finished, but the system could not be imported back into '
                                              'EasyHybrid:\n{}'.format(error))
            self.label_info.set_text('NAMD finished, but importing the result failed -- see the message above.')
            return

        existing_dcd_paths = [p for p in dcd_paths if os.path.isfile(p)]
        if not existing_dcd_paths:
            self.label_info.set_text(
                'NAMD finished and "{}" was added, but no trajectory (.dcd) frame was written '
                '(run too short for the output frequency) -- nothing to import.'.format(system_name))
            return

        system_id = list(self.p_session.psystem.keys())[-1]
        vobject = None
        try:
            for dcd_path in existing_dcd_paths:
                parameters = {
                    'system_id': system_id, 'data_path': dcd_path, 'data_type': 'dcd',
                    'new_vobj_name': None if vobject else system_name + '_traj',
                    'vobject_id': vobject.index if vobject else None, 'vobject': vobject,
                    'logfile': None, 'first': None, 'last': None, 'stride': None,
                }
                self.p_session.import_data(parameters)
                vobject = parameters['vobject']
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg='The system was imported, but the trajectory could not be loaded:\n{}'.format(error))
            self.label_info.set_text('NAMD finished and "{}" was added, but loading the trajectory failed.'.format(system_name))
            return

        n_skipped = len(dcd_paths) - len(existing_dcd_paths)
        if n_skipped:
            self.label_info.set_text(
                'Done! "{}" and its trajectory were added to the treeview ({} of {} step(s) wrote no '
                'frames -- run too short for their output frequency, skipped).'.format(
                    system_name, n_skipped, len(dcd_paths)))
        else:
            self.label_info.set_text('Done! "{}" and its trajectory were added to the treeview.'.format(system_name))

    #-------------------------------------------------------------------------
    #  S T O P
    #-------------------------------------------------------------------------
    def on_button_stop_clicked(self, widget):
        """ Function doc """
        if self.current_process is None:
            return
        if not self.main.simple_dialog.question('Stop the current NAMD run?'):
            return

        if self._poll_timeout_id is not None:
            GLib.source_remove(self._poll_timeout_id)
            self._poll_timeout_id = None

        pids_to_watch = _terminate_process_tree(self.current_process.pid)
        self.label_info.set_text('Stopping...')
        GLib.timeout_add(200, self._check_stop_progress, pids_to_watch, time.time())

    def _check_stop_progress(self, pids_to_watch, started_at):
        """ Same non-blocking grace-period pattern as
            process_manager_window.py's own _check_abort_progress() --
            checked every 200ms for up to 5 seconds before escalating to
            SIGKILL, without freezing the GTK main thread while waiting.
        """
        GRACE_PERIOD_SECONDS = 5.0

        still_alive = any(_pid_alive(pid) for pid in pids_to_watch)
        if still_alive and (time.time() - started_at) < GRACE_PERIOD_SECONDS:
            return True

        if still_alive:
            _force_kill_pids(pids_to_watch)

        if self.current_process is not None:
            try:
                self.current_process.wait(timeout=1)
            except Exception:
                pass
        self.current_process = None

        if self._run_mode == 'pipeline' and self._protocol_index is not None:
            self._set_step_status(self._protocol_index, 'Stopped')
            for i in range(self._protocol_index + 1, self._step_count()):
                self._set_step_status(i, 'Skipped')

        self._set_controls_sensitive(True)
        self.label_info.set_text('Stopped by user.')
        return False
