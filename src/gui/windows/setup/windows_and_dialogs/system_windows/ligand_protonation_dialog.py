#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: "Ligand protonation" review dialog (Add Missing Hydrogens)
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      [EN] 2026-10-18, user's request: inspect the protonation a ligand gets
#      and change it where needed. Lists the N/O/S atoms of one ligand with
#      the state each has in the baseline (current hydrogens, or OpenBabel's
#      prediction at a pH) and a combo to add/remove a proton; the net charge
#      follows. OK returns a plan (util/hydrogen_builder.ligand_protonation_plan)
#      that the Add Hydrogens window applies.
#
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, GLib

import traceback

from util import hydrogen_builder as HB


class LigandProtonationDialog:
    """ run() -> plan dict, or None when cancelled. """

    # liststore columns
    C_LABEL, C_GROUP, C_BONDED, C_BASE, C_STATE, C_STATES, C_WEIGHT, C_INDEX, C_RULE = range(9)

    def __init__(self, parent_window, main, system, vobject, key, atoms, ph=7.4, previous=None,
                 baseline_fn=None, has_h=None):
        """ baseline_fn(ph or None) -> baseline (hydrogen_builder.ligand_protonation_baseline
            format; None = the current hydrogens): another source of ligand chemistry
            (Prepare OPLS System, 2026-10-18). Default: the pDynamo atoms of `system`. """
        self.main, self.system, self.vobject, self.key = main, system, vobject, key
        self.atoms = list(atoms)
        self.vm_session = main.vm_session
        self.baseline_fn = baseline_fn or (lambda ph: HB.ligand_protonation_baseline(self.system, self.atoms, ph=ph))
        self.has_obabel = True       # native perception (hydrogen_builder._ph_protonation)
        self.has_h = any(a.atomicNumber == 1 for a in self.atoms) if has_h is None else has_h
        self.baseline, self.sites, self.deltas = None, [], {}
        self._recompute_id = None

        chain, resName, resSeq = key
        self.dialog = Gtk.Dialog(title='Ligand protonation -- {} {}'.format(resName, resSeq),
                                 transient_for=parent_window, modal=True)
        self.dialog.set_default_size(820, 460)
        self.dialog.add_button('Reset', 1)
        self.dialog.add_button('Cancel', Gtk.ResponseType.CANCEL)
        ok = self.dialog.add_button('Use this protonation', Gtk.ResponseType.OK)
        ok.get_style_context().add_class('suggested-action')
        self.dialog.set_default_response(Gtk.ResponseType.OK)

        box = self.dialog.get_content_area()
        box.set_spacing(6)
        for side in ('top', 'bottom', 'start', 'end'):
            getattr(box, 'set_margin_' + side)(10)

        # . start from
        bar = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
        bar.pack_start(Gtk.Label(label='Start from:'), False, False, 0)
        self.combo_source = Gtk.ComboBoxText()
        if self.has_h or baseline_fn is None:
            self.combo_source.append('current', 'Current hydrogens' if self.has_h else 'Current (no hydrogens)')
        self.combo_source.append('ph', 'pH model (group pKa rules)')
        bar.pack_start(self.combo_source, False, False, 0)
        bar.pack_start(Gtk.Label(label='pH'), False, False, 0)
        self.spin_ph = Gtk.SpinButton.new_with_range(0.0, 14.0, 0.1)
        self.spin_ph.set_digits(1)
        self.spin_ph.set_value(float(ph))
        bar.pack_start(self.spin_ph, False, False, 0)
        self.check_show = Gtk.CheckButton(label='Show in 3D')
        self.check_show.set_active(vobject is not None)
        self.check_show.set_sensitive(vobject is not None)
        bar.pack_end(self.check_show, False, False, 0)
        box.pack_start(bar, False, False, 0)

        hint = Gtk.Label(xalign=0.0)
        hint.set_line_wrap(True)
        hint.set_markup('<small>N, O and S atoms of the ligand. Pick a state in <b>Set to</b> to add or remove a '
                        'proton (formal charge +1/-1); changed atoms are shown in bold. Other atoms keep their '
                        'hydrogens. The pH model applies group pKa rules (shown in pH rule); check unusual groups.</small>')
        box.pack_start(hint, False, False, 0)

        # . sites
        self.store = Gtk.ListStore(str, str, str, str, str, object, int, int, str)
        self.view = Gtk.TreeView(model=self.store)
        for title, col, expand in (('Atom', self.C_LABEL, False), ('Group', self.C_GROUP, True),
                                   ('Bonded to', self.C_BONDED, True), ('Baseline', self.C_BASE, False)):
            cell = Gtk.CellRendererText()
            column = Gtk.TreeViewColumn(title, cell, text=col, weight=self.C_WEIGHT)
            column.set_expand(expand)
            column.set_resizable(True)
            self.view.append_column(column)
        cell = Gtk.CellRendererCombo(editable=True, has_entry=False, text_column=0)
        cell.connect('changed', self._on_state_changed)
        column = Gtk.TreeViewColumn('Set to', cell, text=self.C_STATE, weight=self.C_WEIGHT)
        column.set_cell_data_func(cell, self._state_model)
        self.view.append_column(column)
        cell = Gtk.CellRendererText()
        column = Gtk.TreeViewColumn('pH rule', cell, text=self.C_RULE)
        column.set_resizable(True)
        self.view.append_column(column)
        self.view.get_selection().connect('changed', self._on_row_selected)
        scrolled = Gtk.ScrolledWindow()
        scrolled.set_vexpand(True)
        scrolled.add(self.view)
        box.pack_start(scrolled, True, True, 0)

        self.label_charge = Gtk.Label(xalign=0.0)
        box.pack_start(self.label_charge, False, False, 0)

        default = 'ph' if (self.has_obabel and not self.has_h) else 'current'
        if previous is not None:
            default = previous['source']
            if previous.get('ph') is not None:
                self.spin_ph.set_value(previous['ph'])
        self.combo_source.set_active_id(default)
        self._recompute()
        if previous is not None:
            # . keep the edits made the last time the dialog was opened
            self.deltas = {k: v for k, v in previous.get('deltas', {}).items()
                           if any(id(s['atom']) == k for s in self.sites)}
            self._fill()
        self.combo_source.connect('changed', lambda *a: self._recompute())
        self.spin_ph.connect('value-changed', self._on_ph_changed)
        self.dialog.show_all()

    # -------------------------------------------------------------------
    def _on_ph_changed(self, spin):
        if self.combo_source.get_active_id() != 'ph':
            return
        if self._recompute_id is not None:
            GLib.source_remove(self._recompute_id)
        self._recompute_id = GLib.timeout_add(350, self._recompute_later)

    def _recompute_later(self):
        self._recompute_id = None
        self._recompute()
        return False

    def _recompute(self):
        source = self.combo_source.get_active_id()
        if source == 'ph' and not self.has_obabel:
            self.combo_source.set_active_id('current')
            return
        self.spin_ph.set_sensitive(source == 'ph')
        ph = self.spin_ph.get_value() if source == 'ph' else None
        try:
            baseline = self.baseline_fn(ph)
        except Exception:
            traceback.print_exc()
            baseline = None
        if baseline is None:
            self.label_charge.set_markup('<span foreground="#c01c28">The pH model failed for this ligand; '
                                         'showing its current hydrogens.</span>')
            baseline = self.baseline_fn(None)
        self.baseline = baseline
        self.sites = HB.ligand_protonation_sites(baseline)
        self.deltas = {}
        self._fill()

    def _fill(self):
        self.store.clear()
        for k, site in enumerate(self.sites):
            base = HB.proton_state_text(site['symbol'], site['n_h'], site['charge'])
            delta = self.deltas.get(id(site['atom']), 0)
            state = dict(site['states']).get(delta, base)
            self.store.append([site['label'], site['group'], site['bonded'], base, state, site['states'],
                               700 if delta else 400, k, site.get('rule', '')])
        self._update_charge()

    def _state_model(self, column, cell, model, treeiter, data=None):
        options = Gtk.ListStore(str)
        for _delta, text in model[treeiter][self.C_STATES]:
            options.append([text])
        cell.set_property('model', options)

    def _on_state_changed(self, combo, path, new_iter):
        text = combo.get_property('model')[new_iter][0]
        row = self.store[path]
        site = self.sites[row[self.C_INDEX]]
        delta = {t: d for d, t in site['states']}.get(text, 0)
        if delta:
            self.deltas[id(site['atom'])] = delta
        else:
            self.deltas.pop(id(site['atom']), None)
        row[self.C_STATE] = text
        row[self.C_WEIGHT] = 700 if delta else 400
        self._update_charge()

    def _update_charge(self):
        if self.baseline is None:
            return
        base = HB.ligand_net_charge(self.baseline)
        net = HB.ligand_net_charge(self.baseline, self.deltas)
        source = ('current hydrogens' if self.baseline['source'] == 'current'
                  else 'pH model at pH {:.1f}'.format(self.baseline['ph']))
        text = 'Net charge of the ligand: <b>{:+d}</b>'.format(net)
        if net != base or self.deltas:
            text += '   (baseline {:+d}, {} change(s))'.format(base, len(self.deltas))
        text += '   ·   from {}'.format(source)
        self.label_charge.set_markup(text)

    def _on_row_selected(self, selection):
        model, treeiter = selection.get_selected()
        if treeiter is None or not self.check_show.get_active() or self.vobject is None:
            return
        atom = self.sites[model[treeiter][self.C_INDEX]]['atom']
        if getattr(atom, 'index', None) is None:
            return
        vatom = self.vobject.atoms.get(atom.index)
        if vatom is None:
            return
        session = self.vm_session
        try:
            frame = getattr(session, 'selection_box_frame', None)
            if frame is not None:
                frame.change_toggle_button_selecting_mode_status(True)
            else:
                session.picking_selection_mode = True
            picking = session.picking_selections
            picking.selection_function_picking(None)
            picking.selection_function_picking(vatom)
            import numpy as np
            session.vm_glcore.center_on_coordinates(self.vobject, np.array(vatom.coords(), dtype=np.float32))
            session.vm_glcore.queue_draw()
        except Exception:
            traceback.print_exc()

    # -------------------------------------------------------------------
    def run(self):
        while True:
            response = self.dialog.run()
            if response == 1:          # Reset
                self.deltas = {}
                self._fill()
                continue
            break
        plan = None
        if response == Gtk.ResponseType.OK and self.baseline is not None:
            plan = HB.ligand_protonation_plan(self.baseline, self.deltas)
            plan['deltas'] = dict(self.deltas)
            plan['h_counts'] = {a.label: self.baseline['n_h'][id(a)] + self.deltas.get(id(a), 0)
                                for a in self.baseline['heavy']}
        if self._recompute_id is not None:
            GLib.source_remove(self._recompute_id)
        self.dialog.destroy()
        return plan
