#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: Prepare AMBER System (tLeap) window
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
#      "Prepare AMBER System" window -- builds a tleap (AmberTools) input
#      script for a system already loaded in EasyHybrid, runs the real
#      `tleap` executable, and imports the resulting prmtop/inpcrd back
#      in as a new AMBER system. Ported from a never-finished GTK2
#      mockup in the old gtkDynamo predecessor project
#      (gui/_old_gui/22_window_tleap.glade in gtkdynamo2-master) -- that
#      old repo had no Python logic behind it at all, so everything here
#      (the tleap script generation and subprocess handling in
#      util/tleap_runner.py, and this window) is a fresh implementation,
#      only the general layout/fields are inspired by that old design.
#

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

import os
import traceback

from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from gui.widgets.custom_widgets import FolderChooserButton
from util import tleap_runner
from pdynamo.pDynamo2EasyHybrid.helpers import export_special_PDB


class PrepareAmberSystemWindow:
    """ "Prepare AMBER System (tLeap)" window. """

    def __init__(self, main=None):
        """ Class initialiser """
        self.main       = main
        self.vm_session = main.vm_session
        self.p_session  = main.p_session
        self.home       = main.home
        self.Visible    = False

    #-------------------------------------------------------------------------
    #  W I N D O W   L I F E C Y C L E
    #-------------------------------------------------------------------------
    def open_window(self):
        """ Function doc """
        if self.Visible:
            self.window.present()
            return

        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(self.home, 'src/gui/windows/setup/prepare_amber_system_window.glade'))
        self.builder.connect_signals(self)

        self.window = self.builder.get_object('window')
        self.window.set_title('Prepare AMBER System (tLeap)')

        # -------------------- widget shortcuts --------------------
        self.label_info          = self.builder.get_object('label_info')
        self.combo_protein_ff    = self.builder.get_object('combo_protein_ff')
        self.combo_gaff_ff       = self.builder.get_object('combo_gaff_ff')
        self.combo_glycam_ff     = self.builder.get_object('combo_glycam_ff')
        self.checkbox_solvate    = self.builder.get_object('checkbox_solvate')
        self.frame_solvate       = self.builder.get_object('frame_solvate')
        self.combo_water_model   = self.builder.get_object('combo_water_model')
        self.combo_box_type      = self.builder.get_object('combo_box_type')
        self.spinbtn_buffer      = self.builder.get_object('spinbtn_buffer')
        self.radio_ions_none        = self.builder.get_object('radio_ions_none')
        self.radio_ions_neutralize  = self.builder.get_object('radio_ions_neutralize')
        self.radio_ions_add         = self.builder.get_object('radio_ions_add')
        self.grid_ions_add          = self.builder.get_object('grid_ions_add')
        self.combo_cation        = self.builder.get_object('combo_cation')
        self.spinbtn_n_cations   = self.builder.get_object('spinbtn_n_cations')
        self.combo_anion         = self.builder.get_object('combo_anion')
        self.spinbtn_n_anions    = self.builder.get_object('spinbtn_n_anions')
        self.entry_system_name   = self.builder.get_object('entry_system_name')
        self.entry_tleap_command = self.builder.get_object('entry_tleap_command')
        self.treeview_extra_files    = self.builder.get_object('treeview_extra_files')
        self.liststore_extra_files   = self.builder.get_object('liststore_extra_files')
        self.treeview_bonds          = self.builder.get_object('treeview_bonds')
        self.liststore_bonds         = self.builder.get_object('liststore_bonds')

        # -------------------- system / frame comboboxes --------------------
        self.box_system = self.builder.get_object('box_system')
        self.combobox_systems = SystemComboBox(self.main)
        self.combobox_systems.connect("changed", self.on_combobox_systemsbox_changed)
        self.box_system.pack_start(self.combobox_systems, False, False, 0)

        self.box_coordinates = self.builder.get_object('box_coordinates')
        self.coordinates_combobox = CoordinatesComboBox()
        self.box_coordinates.pack_start(self.coordinates_combobox, False, False, 0)

        # -------------------- working-folder chooser --------------------
        self.box_work_folder = self.builder.get_object('box_work_folder')
        self.work_folder_chooser = FolderChooserButton(self.main, 'folder', self.home)
        self.box_work_folder.pack_start(self.work_folder_chooser.btn, False, False, 0)
        scratch = os.environ.get('PDYNAMO3_SCRATCH')
        if scratch and os.path.isdir(scratch):
            self.work_folder_chooser.set_folder(folder=scratch)

        # -------------------- populate force field / water comboboxes --------------------
        self._populate_combo(self.combo_protein_ff, tleap_runner.list_leaprc_files('protein.'),
                              preferred='protein.ff14SB')
        self._populate_combo(self.combo_gaff_ff, ['None'] + tleap_runner.list_leaprc_files('gaff'))
        self._populate_combo(self.combo_glycam_ff, ['None'] + tleap_runner.list_leaprc_files('GLYCAM_'))
        self._populate_combo(self.combo_water_model, tleap_runner.list_leaprc_files('water.'),
                              preferred='water.tip3p')
        self._populate_combo(self.combo_cation, tleap_runner.COMMON_CATIONS)
        self._populate_combo(self.combo_anion, tleap_runner.COMMON_ANIONS)

        # -------------------- tleap executable --------------------
        tleap_command = self.main.vm_session.vm_config.gl_parameters.get('tleap_command')
        if not tleap_command or not os.path.isfile(tleap_command):
            tleap_command = tleap_runner.find_tleap_executable()
        self.entry_tleap_command.set_text(tleap_command or '')
        if not tleap_command:
            self.label_info.set_text(
                'tleap was not found automatically (checked PATH and $AMBERHOME/bin) -- '
                'please enter the full path to the executable below.')

        # -------------------- extra files (kept in Python, mirrored into the treeview) --------------------
        self.extra_files = []

        # -------------------- bonds/links (kept in Python, mirrored into the treeview) --------------------
        # Each entry: ((resnum1, atomname1), (resnum2, atomname2)) -- see
        # on_button_add_bond_clicked() and tleap_runner.build_tleap_script's
        # "bonds" parameter.
        self.bonds = []

        self.frame_solvate.set_sensitive(self.checkbox_solvate.get_active())
        self.grid_ions_add.set_sensitive(self.radio_ions_add.get_active())

        if self.p_session.psystem:
            self.combobox_systems.set_active_system(e_id=self.p_session.active_id)

        self.window.show_all()
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
    def _populate_combo(self, combo, items, preferred=None):
        """ Fills a GtkComboBoxText with `items` (strings), selecting
            `preferred` if present, otherwise the first item. Leaves the
            combobox with no active item if `items` is empty (e.g. no
            $AMBERHOME leaprc.* files of that kind were found).
        """
        combo.remove_all()
        active_index = 0
        for index, item in enumerate(items):
            combo.append_text(item)
            if preferred is not None and item == preferred:
                active_index = index
        if items:
            combo.set_active(active_index)

    def _selected_extra_file_paths(self):
        return list(self.extra_files)

    #-------------------------------------------------------------------------
    #  S I G N A L S :  S Y S T E M   /   O B J E C T   S E L E C T I O N
    #-------------------------------------------------------------------------
    def on_combobox_systemsbox_changed(self, widget):
        """ Function doc """
        system_id = self.combobox_systems.get_system_id()
        if system_id is not None:
            self.coordinates_combobox.set_model(self.main.vobject_liststore_dict[system_id])
            size = len(list(self.main.vobject_liststore_dict[system_id]))
            self.coordinates_combobox.set_active(size - 1)

    #-------------------------------------------------------------------------
    #  S I G N A L S :  C H E C K B O X E S   /   R A D I O S
    #-------------------------------------------------------------------------
    def on_checkbox_solvate_toggled(self, widget):
        """ Function doc """
        self.frame_solvate.set_sensitive(widget.get_active())

    def on_radio_ions_add_toggled(self, widget):
        """ Function doc """
        self.grid_ions_add.set_sensitive(widget.get_active())

    #-------------------------------------------------------------------------
    #  S I G N A L S :  A D D I T I O N A L   F I L E S
    #-------------------------------------------------------------------------
    def on_button_add_file_clicked(self, widget):
        """ Function doc """
        dialog = Gtk.FileChooserDialog(
            title="Add parameter file (.frcmod / .lib / .off / .mol2)",
            parent=self.window,
            action=Gtk.FileChooserAction.OPEN,
        )
        dialog.set_select_multiple(True)
        dialog.add_buttons(
            Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
            Gtk.STOCK_OPEN, Gtk.ResponseType.OK,
        )
        file_filter = Gtk.FileFilter()
        file_filter.set_name("tLeap parameter files (*.frcmod, *.lib, *.off, *.mol2)")
        for pattern in ("*.frcmod", "*.lib", "*.off", "*.mol2"):
            file_filter.add_pattern(pattern)
        dialog.add_filter(file_filter)

        if dialog.run() == Gtk.ResponseType.OK:
            for path in dialog.get_filenames():
                ext = os.path.splitext(path)[1].lower().lstrip('.')
                self.extra_files.append(path)
                self.liststore_extra_files.append([os.path.basename(path), ext])
        dialog.destroy()

    def on_button_clear_files_clicked(self, widget):
        """ Function doc """
        self.extra_files = []
        self.liststore_extra_files.clear()

    #-------------------------------------------------------------------------
    #  S I G N A L S :  B O N D S   /   L I N K S
    #-------------------------------------------------------------------------
    def _picked_pair_label(self, atom):
        """ Short display string for a picked atom in the bonds treeview,
            e.g. "CYS 10 / SG". """
        return "{} {} / {}".format(atom.residue.name, atom.residue.index, atom.name)

    def on_button_add_bond_clicked(self, widget):
        """ Reads whatever is currently picked in the 3D view (Picking
            mode -- see vm_session.picking_selections.picking_selections_list,
            a 4-slot [pk1, pk2, pk3, pk4] buffer) and adds any complete
            pair(s) found there (slots 0+1, and/or 2+3) to self.bonds as
            an explicit tleap `bond` -- for connections (disulfide
            bridges, metal coordination, linking a ligand to the
            protein, ...) tleap wouldn't otherwise detect by distance.

            Both atoms of a pair must belong to the SAME vismol object
            selected in the "Object / frame" combobox above -- that is
            the exact object export_special_PDB() turns into the PDB
            tleap loads, and a picked atom's .residue.index/.name are
            the very same values that end up in that PDB's residue-
            number/atom-name columns (verified against real tleap runs
            -- see tleap_runner.build_tleap_script's "bonds" docstring).
        """
        vobject_id = self.coordinates_combobox.get_vobject_id()
        if vobject_id is None:
            self.main.simple_dialog.info(msg='Please select a system/object first.')
            return
        vobject = self.main.vm_session.vm_objects_dic.get(vobject_id)
        if vobject is None:
            self.main.simple_dialog.info(msg='Please select a system/object first.')
            return

        picked = self.main.vm_session.picking_selections.picking_selections_list
        candidate_pairs = [(picked[0], picked[1]), (picked[2], picked[3])]

        added = 0
        mismatched = 0
        for atom_a, atom_b in candidate_pairs:
            if atom_a is None or atom_b is None:
                continue
            if atom_a.vm_object is not vobject or atom_b.vm_object is not vobject:
                mismatched += 1
                continue
            pair = ((atom_a.residue.index, atom_a.name), (atom_b.residue.index, atom_b.name))
            self.bonds.append(pair)
            self.liststore_bonds.append([self._picked_pair_label(atom_a), self._picked_pair_label(atom_b)])
            added += 1

        if added:
            # Clear the picking buffer so the next pick starts clean,
            # without the user having to click empty space first.
            self.main.vm_session.picking_selections.selection_function_picking(None)
            self.main.vm_session.vm_glcore.queue_draw()

        if mismatched:
            self.main.simple_dialog.info(
                msg='{} picked pair(s) belonged to a different object than the one selected above '
                    'and were skipped.'.format(mismatched))
        elif not added:
            self.main.simple_dialog.info(
                msg='No complete pair is currently picked. Switch to "Picking" mode in the toolbar '
                    'and click two atoms in the 3D view first.')

    def on_button_clear_bonds_clicked(self, widget):
        """ Function doc """
        self.bonds = []
        self.liststore_bonds.clear()

    #-------------------------------------------------------------------------
    #  S I G N A L S :  B U T T O N S
    #-------------------------------------------------------------------------
    def on_button_cancel_clicked(self, widget):
        """ Function doc """
        self.close_window()

    #-------------------------------------------------------------------------
    #  R U N
    #-------------------------------------------------------------------------
    def on_button_run_clicked(self, widget):
        """ Function doc """
        system_id = self.combobox_systems.get_system_id()
        if system_id is None:
            self.main.simple_dialog.info(msg='Please select a system.')
            return

        vobject_id = self.coordinates_combobox.get_vobject_id()
        if vobject_id is None:
            self.main.simple_dialog.info(msg='Please select an object/frame.')
            return
        vobject = self.main.vm_session.vm_objects_dic.get(vobject_id)
        if vobject is None:
            self.main.simple_dialog.info(msg='Please select an object/frame.')
            return

        tleap_command = self.entry_tleap_command.get_text().strip()
        if not tleap_command or not (os.path.isfile(tleap_command) and os.access(tleap_command, os.X_OK)):
            self.main.simple_dialog.info(msg='Please provide a valid path to the tleap executable.')
            return
        # Persisted so the next time this window opens it doesn't have to
        # re-detect it (same pattern as other per-preference settings --
        # see vm_config.gl_parameters).
        self.main.vm_session.vm_config.gl_parameters['tleap_command'] = tleap_command

        system_name = self.entry_system_name.get_text().strip() or 'prepared_system'
        work_folder = self.work_folder_chooser.get_folder()
        if not work_folder:
            self.main.simple_dialog.info(msg='Please choose a working folder.')
            return

        protein_ff = self.combo_protein_ff.get_active_text()
        gaff_ff    = self.combo_gaff_ff.get_active_text()
        glycam_ff  = self.combo_glycam_ff.get_active_text()
        gaff_ff   = None if gaff_ff   in (None, 'None') else gaff_ff
        glycam_ff = None if glycam_ff in (None, 'None') else glycam_ff

        solvate = None
        if self.checkbox_solvate.get_active():
            water_model = self.combo_water_model.get_active_text()
            if not water_model:
                self.main.simple_dialog.info(msg='Please select a water model, or uncheck "Solvate".')
                return
            solvate = {
                'water_model': water_model[len('water.'):] if water_model.startswith('water.') else water_model,
                'box_type': self.combo_box_type.get_active_id() or 'box',
                'buffer': self.spinbtn_buffer.get_value(),
            }

        ions = None
        if self.radio_ions_neutralize.get_active():
            ions = {'mode': 'neutralize'}
        elif self.radio_ions_add.get_active():
            ions = {
                'mode': 'add',
                'cation': self.combo_cation.get_active_text(),
                'n_cation': int(self.spinbtn_n_cations.get_value()),
                'anion': self.combo_anion.get_active_text(),
                'n_anion': int(self.spinbtn_n_anions.get_value()),
            }

        # . Export the currently selected object/frame to a PDB -- this
        #   is what tleap's loadpdb actually reads.
        pdb_path = os.path.join(work_folder, system_name + '_input.pdb')
        try:
            # frame=-1 -- the "Object / frame" combobox above selects
            # WHICH vobject to use (not a frame index within it, despite
            # the name -- see get_vobject_id()), so the actual frame
            # exported is just that vobject's current/last one.
            export_special_PDB(vobject=vobject, frame=-1, output=pdb_path)
            tleap_runner.insert_ter_records(pdb_path)
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg='Could not export the structure to PDB:\n{}'.format(error))
            return

        script_text = tleap_runner.build_tleap_script(
            pdb_path=pdb_path,
            protein_ff=protein_ff,
            glycam_ff=glycam_ff,
            gaff_ff=gaff_ff,
            extra_files=self._selected_extra_file_paths(),
            bonds=self.bonds,
            solvate=solvate,
            ions=ions,
            output_basename=system_name,
        )

        self.label_info.set_text('Running tleap...')
        while Gtk.events_pending():
            Gtk.main_iteration()

        result = tleap_runner.run_tleap(script_text, workdir=work_folder,
                                         tleap_command=tleap_command, output_basename=system_name)

        if not result['ok']:
            self.label_info.set_text('tleap failed -- see the error dialog for details.')
            self.main.simple_dialog.error_details(
                parent=self.window,
                msg='tleap did not produce a valid AMBER system (prmtop/inpcrd). '
                    'The generated script is at "{}".'.format(os.path.join(work_folder, 'tleap.in')),
                details=script_text + '\n\n----- tleap output -----\n\n' + result['stdout'] + result['stderr'],
                title='tLeap Error',
            )
            return

        # . Import the prepared system back into EasyHybrid, reusing the
        #   existing (and already hardened -- see define_NBModel's
        #   return-value check) AMBER import path.
        try:
            self.p_session.load_a_new_pDynamo_system_from_dict(
                input_files={'amber_prmtop': result['prmtop'], 'coordinates': result['inpcrd']},
                system_type=0,
                name=system_name,
                tag='AMBER',
                working_folder=work_folder,
            )
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg='tleap succeeded, but the resulting system could not be '
                                              'imported into EasyHybrid:\n{}'.format(error))
            self.label_info.set_text('tleap succeeded, but importing the result failed -- see the message above.')
            return

        self.label_info.set_text('Done! "{}" was added to the treeview.'.format(system_name))
        self.close_window()
