#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: "Add Missing Hydrogens" window
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
#      GUI wiring for util/hydrogen_builder.py -- see that module's own
#      docstring for the actual chemistry (pDynamo3's PDBComponentLibrary
#      residue templates + BuildHydrogenCoordinates3FromConnectivity).
#      This window: lets the user pick a loaded system and a scope
#      (entire system / residues touched by the current 3D selection),
#      runs analyze_system(), shows a review list (only residues that
#      actually need something, with a protonation-state combo for the
#      genuinely ambiguous ones -- His/Asp/Glu/Cys/Lys), then commits via
#      rebuild_system_with_added_hydrogens() and refreshes the SAME
#      VismolObject in place (eSession._refresh_vobject_from_pdynamo_system,
#      new this feature -- no existing precedent for updating an
#      already-displayed object's atom count without discarding it).
#
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

import os
import traceback

from gui.widgets.custom_widgets import SystemComboBox
from util import hydrogen_builder


class PrepareAddHydrogensWindow:

    def __init__(self, main=None):
        """ Class initialiser """
        self.main = main
        self.home = main.home
        self.p_session = main.p_session
        self.vm_session = main.vm_session
        self.Visible = False

        self._results = []
        self._system_id = None
        self._undo_snapshot = None   # (system_id, old_system) -- one level, this tool's own runs only

    # -------------------------------------------------------------------
    #  Window setup
    # -------------------------------------------------------------------

    def open_window(self):
        """ Function doc """
        if self.Visible:
            self.window.present()
            return

        self.builder = Gtk.Builder()
        self.builder.add_from_file(os.path.join(
            self.home, 'src/gui/windows/setup/prepare_add_hydrogens_window.glade'))
        self.builder.connect_signals(self)

        self.window = self.builder.get_object('window')
        self.box_system_combo = self.builder.get_object('box_system_combo')
        self.combobox_systems = SystemComboBox(self.main)
        self.box_system_combo.pack_start(self.combobox_systems, True, True, 0)
        self.radio_scope_entire_system = self.builder.get_object('radio_scope_entire_system')
        self.radio_scope_selected_atoms = self.builder.get_object('radio_scope_selected_atoms')
        self.liststore_residues = self.builder.get_object('liststore_residues')
        self.treeview_residues = self.builder.get_object('treeview_residues')
        self.cellrenderer_choice = self.builder.get_object('cellrenderer_choice')
        self.treeviewcolumn_choice = self.builder.get_object('treeviewcolumn_choice')
        self.treeviewcolumn_choice.set_cell_data_func(self.cellrenderer_choice, self._do_get_choice_model)
        self.label_status = self.builder.get_object('label_status')
        self.button_add_hydrogens = self.builder.get_object('button_add_hydrogens')
        self.button_undo = self.builder.get_object('button_undo')
        self.radio_remove_all = self.builder.get_object('radio_remove_all')
        self.radio_remove_apolar_only = self.builder.get_object('radio_remove_apolar_only')

        self.window.show_all()
        self._results = []
        self._system_id = None
        self._undo_snapshot = None
        self.button_add_hydrogens.set_sensitive(False)
        self.button_undo.set_sensitive(False)
        self.Visible = True

    def close_window(self, button=None, data=None):
        """ Function doc """
        if not self.Visible:
            return
        self.window.destroy()
        self.Visible = False

    def on_button_close_clicked(self, widget):
        self.close_window()

    # -------------------------------------------------------------------
    #  Analyze
    # -------------------------------------------------------------------

    def on_button_analyze_clicked(self, widget):
        system_id = self.combobox_systems.get_system_id()
        if system_id is None:
            self.main.simple_dialog.info(msg='Choose a system first.')
            return
        system = self.p_session.psystem.get(system_id)
        if system is None:
            self.main.simple_dialog.info(msg='Could not find the pDynamo system for the selected object.')
            return

        residue_keys = None
        if self.radio_scope_selected_atoms.get_active():
            residue_keys = self._selected_residue_keys(system_id, system)
            if not residue_keys:
                self.main.simple_dialog.info(
                    msg='No atoms of this system are currently selected in the 3D view.')
                return

        try:
            results = hydrogen_builder.analyze_system(system, residue_keys=residue_keys)
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg='Analysis failed:\n{}'.format(error))
            return

        self._system_id = system_id
        self._results = results
        self._refresh_residue_liststore()

        n_standard = sum(1 for r in results if r['is_standard_residue'])
        n_skipped = sum(1 for r in results if not r['is_standard_residue'])
        n_missing = sum(len(r.get('missing_hydrogens', [])) for r in results)
        status = '{} residue(s) need hydrogens ({} atom(s) total).'.format(n_standard, n_missing)
        if n_skipped:
            status += ' {} non-standard residue(s) skipped.'.format(n_skipped)
        if not results:
            status = 'Nothing to do -- every standard residue already has all its hydrogens.'
        self.label_status.set_text(status)
        self.button_add_hydrogens.set_sensitive(bool(results))

    def _selected_residue_keys(self, system_id, system):
        """ Maps the current 3D-view selection's atoms (only the ones
            belonging to the chosen system's own VismolObject) back to
            residue keys, via vismol Atom.index (1-based) - 1 = the
            pDynamo/global atom index -- the exact conversion
            xtb_fragment_charges_window.py's own
            _get_current_selection_indexes() already relies on.
        """
        vm_object = self.vm_session.vm_objects_dic.get(system_id)
        if vm_object is None:
            return set()
        selection = self.vm_session.selections[self.vm_session.current_selection]
        _residues, atom_index_to_key = hydrogen_builder.build_residue_index(system)
        keys = set()
        for atom in selection.selected_atoms:
            if atom.vm_object is not vm_object:
                continue
            pdynamo_index = atom.index - 1
            key = atom_index_to_key.get(pdynamo_index)
            if key is not None:
                keys.add(key)
        return keys

    def _refresh_residue_liststore(self):
        self.liststore_residues.clear()
        for r in self._results:
            chain, resName, resSeq = r['key']
            display = '{} / {} {}'.format(chain, resName, resSeq)
            if not r['is_standard_residue']:
                if r.get('is_organic_candidate'):
                    # . No PDBComponentLibrary entry, but has heavy atoms
                    # -- a ligand/cofactor/organic molecule, handled at
                    # commit time via hydrogen_builder.protonate_non_standard_residue()
                    # (OpenBabel, pH-aware) instead of being skipped.
                    status = 'organic/ligand -- {} heavy atom(s)'.format(r.get('n_heavy_atoms', 0))
                    self.liststore_residues.append([
                        display, status, '', False, '', chain, resName, resSeq, False, '7.4', True])
                else:
                    self.liststore_residues.append([
                        display, 'non-standard (skipped)', '', False, '', chain, resName, resSeq, False, '', False])
                continue
            n_missing = len(r['missing_hydrogens'])
            status_bits = ['{} missing'.format(n_missing)]
            if r.get('disulfide_bonded'):
                status_bits.append('disulfide-bonded')
            if r['unbuildable']:
                status_bits.append('{} not addable ({})'.format(
                    len(r['unbuildable']), ', '.join(r['unbuildable'])))
            ambiguous = r.get('ambiguous_choices')
            choices_csv = '|'.join(ambiguous) if ambiguous else ''
            self.liststore_residues.append([
                display, '; '.join(status_bits), r['default_choice'], bool(ambiguous),
                choices_csv, chain, resName, resSeq, True, '', False])

    def on_residue_choice_edited(self, renderer, path, new_text):
        self.liststore_residues[path][2] = new_text

    def on_residue_ph_edited(self, renderer, path, new_text):
        try:
            float(new_text)
        except ValueError:
            self.main.simple_dialog.info(msg='pH must be a number (e.g. 7.4).')
            return
        self.liststore_residues[path][9] = new_text

    def _do_get_choice_model(self, column, cell, model, treeiter, data=None):
        """ CellDataFunc: builds a per-row options GtkListStore for the
            combo renderer, since different rows (HIS vs ASP/GLU vs
            CYS/LYS vs non-ambiguous) have different valid choices --
            a single static Glade-defined combo model can't do this,
            so it's built here in code instead.
        """
        choices_csv = model[treeiter][4]
        if choices_csv:
            options = Gtk.ListStore(str)
            for choice in choices_csv.split('|'):
                options.append([choice])
            cell.set_property('model', options)

    # -------------------------------------------------------------------
    #  Add Hydrogens / Remove Hydrogens / Undo
    # -------------------------------------------------------------------

    def _commit_rebuilt_system(self, old_system, new_system, vm_object):
        """ Shared by on_button_add_hydrogens_clicked() and
            on_button_remove_hydrogens_clicked(): both build a brand-new
            pDynamo System (System.FromSequence()) -- neither the
            addition nor the removal path can mutate System.atoms/
            coordinates3 in place (confirmed: those are fixed-size
            pDynamo Cython containers) -- and a fresh System carries
            NONE of EasyHybrid's own bookkeeping attributes, ALL named
            "e_*" and attached at original registration time by
            p_session.add_new_system_to_psession() (e_working_folder,
            e_color_palette, e_selections, e_custom_colors,
            e_job_history, e_logfile_data, e_restraints_dict,
            e_annotations, e_treeview_iter, e_liststore_iter,
            e_charges_backup, ...). Confirmed real (twice, for two
            different call sites of this same problem): without this,
            "Export data" crashed with "no attribute 'e_working_folder'"
            and the QC-model window crashed with "no attribute
            'e_qc_table'" on a system that had just gone through the
            Add Hydrogens tool. These belong to the LOGICAL system, not
            the specific pDynamo object identity, so copying most of
            them forward verbatim is correct.

            EXCEPT the handful that store raw ATOM INDICES
            (e_qc_table/e_qc_residue_table/e_fixed_table -- QM/MM region
            definitions; e_color_segments; e_bonds) -- adding OR
            removing atoms shifts the index space, so an old index
            copied verbatim can silently point at a DIFFERENT atom in
            the new system. These must be explicitly RESET to the exact
            same empty defaults add_new_system_to_psession() itself
            uses -- NOT copied (stale indices), NOT left unset either
            (System.FromSequence() never sets them itself, so leaving
            them out entirely reproduces the same AttributeError crash
            the very first time any other tool reads one).

            Returns True if the OLD system had any non-empty
            index-dependent data (so the caller can warn the user it
            was cleared).
        """
        index_dependent = ('e_qc_table', 'e_qc_residue_table', 'e_fixed_table',
                            'e_color_segments', 'e_bonds')
        had_index_dependent_data = any(getattr(old_system, a, None) for a in index_dependent)
        for attr_name, value in vars(old_system).items():
            if attr_name.startswith('e_') and attr_name not in index_dependent:
                setattr(new_system, attr_name, value)
        new_system.e_qc_table = []
        new_system.e_qc_residue_table = {}
        new_system.e_fixed_table = []
        new_system.e_color_segments = {}
        new_system.e_bonds = None
        new_system.e_id = self._system_id
        self._undo_snapshot = (self._system_id, old_system)
        self.p_session.psystem[self._system_id] = new_system
        self.p_session._refresh_vobject_from_pdynamo_system(vm_object=vm_object, system=new_system)
        return had_index_dependent_data

    def on_button_add_hydrogens_clicked(self, widget):
        if self._system_id is None or not self._results:
            return
        system = self.p_session.psystem.get(self._system_id)
        vm_object = self.vm_session.vm_objects_dic.get(self._system_id)
        if system is None or vm_object is None:
            self.main.simple_dialog.info(msg='The selected system is no longer available.')
            return

        residue_variant_choices = {}
        for row in self.liststore_residues:
            chain, resName, resSeq = row[5], row[6], row[7]
            if row[8]:
                choice = row[2]
                variant_label = None if choice in ('', 'Protonated') else choice
                residue_variant_choices[(chain, resName, resSeq)] = variant_label
            elif row[10]:
                residue_variant_choices[(chain, resName, resSeq)] = (
                    hydrogen_builder.ORGANIC_CHOICE_KIND, float(row[9]))

        if not residue_variant_choices:
            self.main.simple_dialog.info(msg='Nothing to add.')
            return

        try:
            report = hydrogen_builder.rebuild_system_with_added_hydrogens(system, residue_variant_choices)
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg='Could not add hydrogens:\n{}'.format(error))
            return

        new_system = report['system']
        if new_system is None:
            self.main.simple_dialog.info(msg='Nothing to add.')
            return

        had_index_dependent_data = self._commit_rebuilt_system(system, new_system, vm_object)

        n_added = sum(len(v) for v in report['added'].values())
        msg_bits = ['{} hydrogen(s) added across {} residue(s).'.format(n_added, len(report['added']))]
        if report['skipped_non_standard']:
            msg_bits.append('{} residue(s) skipped (non-standard).'.format(len(report['skipped_non_standard'])))
        if report['unbuildable']:
            total_unbuildable = sum(len(v) for v in report['unbuildable'].values())
            msg_bits.append('{} atom(s) reported but not added (missing heavy-atom neighbor).'.format(total_unbuildable))
        if report['unbuilt_by_geometry']:
            msg_bits.append('{} hydrogen(s) could not be given a geometry: {}.'.format(
                len(report['unbuilt_by_geometry']), ', '.join(report['unbuilt_by_geometry'])))
        if report.get('charge_changes'):
            change_bits = []
            for key, changes in report['charge_changes'].items():
                chain, resName, resSeq = key
                for (atom_label, new_charge) in changes:
                    change_bits.append('{} {} {} ({:+d})'.format(resName, resSeq, atom_label, new_charge))
            msg_bits.append('Protonation state changed on: {}.'.format(', '.join(change_bits)))
        if had_index_dependent_data:
            msg_bits.append(
                'Note: any previously-defined QC/MM region, fixed-atom selection, or custom '
                'color segments were cleared -- their atom indices are no longer valid after '
                'adding hydrogens. Please redefine them.')
        self.label_status.set_text(' '.join(msg_bits))
        self.main.simple_dialog.info(msg='\n'.join(msg_bits))
        self.button_undo.set_sensitive(True)
        self.button_add_hydrogens.set_sensitive(False)
        self._results = []
        self.liststore_residues.clear()

    def on_button_remove_hydrogens_clicked(self, widget):
        """ Uses the SAME system combo + scope radios as the Add
            Hydrogens flow above (no separate UI needed -- reused
            directly), but does not depend on self._results/Analyze
            having been run first: removal doesn't need per-residue
            review, just a mode (all / nonpolar-only) and a scope.
        """
        system_id = self.combobox_systems.get_system_id()
        if system_id is None:
            self.main.simple_dialog.info(msg='Choose a system first.')
            return
        system = self.p_session.psystem.get(system_id)
        vm_object = self.vm_session.vm_objects_dic.get(system_id)
        if system is None or vm_object is None:
            self.main.simple_dialog.info(msg='Could not find the pDynamo system for the selected object.')
            return

        residue_keys = None
        if self.radio_scope_selected_atoms.get_active():
            residue_keys = self._selected_residue_keys(system_id, system)
            if not residue_keys:
                self.main.simple_dialog.info(
                    msg='No atoms of this system are currently selected in the 3D view.')
                return

        mode = 'apolar' if self.radio_remove_apolar_only.get_active() else 'all'

        try:
            atoms_to_remove = hydrogen_builder.find_removable_hydrogens(system, mode=mode, residue_keys=residue_keys)
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg='Could not analyze hydrogens:\n{}'.format(error))
            return

        if not atoms_to_remove:
            self.main.simple_dialog.info(msg='No matching hydrogens to remove.')
            return

        try:
            new_system = hydrogen_builder.remove_hydrogens_from_system(system, atoms_to_remove)
        except Exception as error:
            traceback.print_exc()
            self.main.simple_dialog.info(msg='Could not remove hydrogens:\n{}'.format(error))
            return

        self._system_id = system_id
        had_index_dependent_data = self._commit_rebuilt_system(system, new_system, vm_object)

        msg_bits = ['{} hydrogen(s) removed ({}).'.format(
            len(atoms_to_remove), 'nonpolar only' if mode == 'apolar' else 'all')]
        if had_index_dependent_data:
            msg_bits.append(
                'Note: any previously-defined QC/MM region, fixed-atom selection, or custom '
                'color segments were cleared -- their atom indices are no longer valid after '
                'removing hydrogens. Please redefine them.')
        self.label_status.set_text(' '.join(msg_bits))
        self.main.simple_dialog.info(msg='\n'.join(msg_bits))
        self.button_undo.set_sensitive(True)
        # . Any previous Analyze results are now stale -- the atom set
        # just changed underneath them.
        self._results = []
        self.liststore_residues.clear()
        self.button_add_hydrogens.set_sensitive(False)

    def on_button_undo_clicked(self, widget):
        if self._undo_snapshot is None:
            return
        system_id, old_system = self._undo_snapshot
        vm_object = self.vm_session.vm_objects_dic.get(system_id)
        if vm_object is None:
            self.main.simple_dialog.info(msg='The object is no longer available -- cannot undo.')
            return
        self.p_session.psystem[system_id] = old_system
        self.p_session._refresh_vobject_from_pdynamo_system(vm_object=vm_object, system=old_system)
        self._undo_snapshot = None
        self.button_undo.set_sensitive(False)
        self.label_status.set_text('Undone -- restored the previous state.')
