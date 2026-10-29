#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: Trajectory Analysis window
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      Analysis > Trajectory Analysis... - one window for the trajectory
#      analyses that read an Object's frames: radius of gyration, secondary
#      structure (DSSP), SASA, hydrogen bonds, contact maps, interaction
#      energies, PCA, cross-correlations, clustering, B-factors and the
#      average structure, quasi-harmonic entropy, side-chain rotamers,
#      diffusion (MSD) and the statistics of any series. The math is in
#      util/traj_analysis.py (GTK-free); the plots in util/analysis_plot.py.
#

import math
import os
import threading
import time
import traceback

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, GLib

import numpy as np

from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from util import traj_analysis as TA
from util.analysis_plot import AnalysisPlot, PALETTE, colormap

WATER_NAMES = {"HOH", "WAT", "SOL", "TIP3", "TIP", "H2O", "T3P", "SPC", "TIP4", "T4P", "DOD", "TP3"}

# . Tien et al., PLoS One 8, e80635 (2013), theoretical maximum ASA (A^2)
MAX_ASA = {"ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225, "GLU": 223, "GLY": 104,
           "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240, "PRO": 159, "SER": 155,
           "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174}

DSSP_COLORS = {"H": (0.85, 0.15, 0.20), "G": (0.98, 0.55, 0.45), "I": (0.60, 0.10, 0.45), "E": (0.95, 0.80, 0.10),
               "B": (0.75, 0.60, 0.10), "T": (0.30, 0.65, 0.85), "S": (0.60, 0.80, 0.60), "-": (0.96, 0.96, 0.96)}

ANALYSES = [
    ("Structure", "rg", "Radius of gyration"),
    ("Structure", "dssp", "Secondary structure (DSSP)"),
    ("Structure", "sasa", "Solvent accessible surface (SASA)"),
    ("Interactions", "hbonds", "Hydrogen bonds"),
    ("Interactions", "contacts", "Contact / distance map"),
    ("Interactions", "energy", "Interaction energy"),
    ("Interactions", "et", "Electron transfer pathways"),
    ("Dynamics", "pca", "Principal components (PCA)"),
    ("Dynamics", "dccm", "Cross-correlation (DCCM)"),
    ("Dynamics", "cluster", "Clustering"),
    ("Dynamics", "bfactors", "B-factors and average structure"),
    ("Dynamics", "qh", "Quasi-harmonic entropy"),
    ("Dynamics", "rotamers", "Side-chain rotamers (chi1, chi2)"),
    ("Dynamics", "msd", "Diffusion (MSD)"),
    ("Statistics", "stats", "Statistics of a series"),
]

ATOM_SUBSETS = ["CA atoms", "Backbone (N, CA, C, O)", "Heavy atoms", "All atoms"]


class _Cancelled(Exception):
    pass


class _LazyFrames:
    """Frames of an object read one at a time (float64), optionally for some atoms."""

    def __init__(self, vobject, frames_idx, atoms=None):
        self.v, self.idx, self.atoms = vobject, list(frames_idx), atoms

    def __len__(self):
        return len(self.idx)

    def __iter__(self):
        for f in self.idx:
            X = np.asarray(self.v.frames[f], dtype=float)
            yield X if self.atoms is None else X[self.atoms]


class TrajectoryAnalysisWindow:

    def __init__(self, main=None):
        self.main = main
        self.home = main.home
        self.p_session = main.p_session
        self.vm_session = main.vm_session
        self.Visible = False
        self.window = None
        self.series_registry = {}       # . name -> (x, y, xlabel, ylabel) for Statistics
        self.result = None
        self._thread = None
        self._cancel = False
        self._progress = 0.0

    # =====================================================================
    #  Window
    # =====================================================================
    def open_window(self):
        if self.Visible:
            self.window.present()
            return
        self.window = Gtk.Window(title="Trajectory Analysis")
        self.window.set_default_size(1240, 780)
        self.window.set_transient_for(self.main.window)
        self.window.connect("destroy", self._on_destroy)

        outer = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6, margin=8)
        self.window.add(outer)
        outer.pack_start(self._build_input_frame(), False, False, 0)

        paned = Gtk.Paned(orientation=Gtk.Orientation.HORIZONTAL)
        outer.pack_start(paned, True, True, 0)

        # . left: analyses + options
        left = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
        left.set_size_request(330, -1)
        self.store = Gtk.TreeStore(str, str)
        parents = {}
        for cat, key, title in ANALYSES:
            if cat not in parents:
                parents[cat] = self.store.append(None, ["<b>{}</b>".format(cat), ""])
            self.store.append(parents[cat], [title, key])
        self.tree = Gtk.TreeView(model=self.store)
        self.tree.set_headers_visible(False)
        renderer = Gtk.CellRendererText()
        self.tree.append_column(Gtk.TreeViewColumn("Analysis", renderer, markup=0))
        self.tree.expand_all()
        self.tree.get_selection().connect("changed", self._on_analysis_changed)
        sw = Gtk.ScrolledWindow()
        sw.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)
        sw.add(self.tree)
        sw.set_size_request(-1, 300)
        left.pack_start(sw, False, False, 0)

        frame = Gtk.Frame(label=" Options ")
        self.stack = Gtk.Stack()
        self.option_widgets = {}
        for _cat, key, _title in ANALYSES:
            self.stack.add_named(self._build_options(key), key)
        frame.add(self.stack)
        left.pack_start(frame, False, False, 0)

        hb = Gtk.Box(spacing=6)
        self.btn_run = Gtk.Button(label="Run")
        self.btn_run.connect("clicked", self.on_run)
        self.btn_cancel = Gtk.Button(label="Cancel")
        self.btn_cancel.set_sensitive(False)
        self.btn_cancel.connect("clicked", lambda *_: setattr(self, "_cancel", True))
        hb.pack_start(self.btn_run, True, True, 0)
        hb.pack_start(self.btn_cancel, False, False, 0)
        left.pack_start(hb, False, False, 0)
        self.progress = Gtk.ProgressBar(show_text=True)
        self.progress.set_text("")
        left.pack_start(self.progress, False, False, 0)
        self.label_hint = Gtk.Label(xalign=0)
        self.label_hint.set_line_wrap(True)
        self.label_hint.set_max_width_chars(44)
        left.pack_start(self.label_hint, False, False, 0)
        paned.pack1(left, False, False)

        # . right: result
        right = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6, margin_start=6)
        bar = Gtk.Box(spacing=6)
        bar.pack_start(Gtk.Label(label="View"), False, False, 0)
        self.combo_view = Gtk.ComboBoxText()
        self.combo_view.connect("changed", self._on_view_changed)
        bar.pack_start(self.combo_view, False, False, 0)
        self.btn_export = Gtk.Button(label="Export Image...")
        self.btn_export.connect("clicked", self.on_export_image)
        self.btn_save = Gtk.Button(label="Save Data...")
        self.btn_save.connect("clicked", self.on_save_data)
        bar.pack_end(self.btn_save, False, False, 0)
        bar.pack_end(self.btn_export, False, False, 0)
        right.pack_start(bar, False, False, 0)

        vpaned = Gtk.Paned(orientation=Gtk.Orientation.VERTICAL)
        self.plot = AnalysisPlot()
        self.plot.on_click = self._on_plot_click
        vpaned.pack1(self.plot, True, False)

        nb = Gtk.Notebook()
        self.table = Gtk.TreeView()
        self.table.connect("row-activated", self._on_table_row)
        sw = Gtk.ScrolledWindow()
        sw.add(self.table)
        nb.append_page(sw, Gtk.Label(label="Table"))
        self.summary = Gtk.TextView(editable=False, monospace=True, wrap_mode=Gtk.WrapMode.WORD)
        sw = Gtk.ScrolledWindow()
        sw.add(self.summary)
        nb.append_page(sw, Gtk.Label(label="Summary"))
        nb.set_size_request(-1, 220)
        vpaned.pack2(nb, False, False)
        vpaned.set_position(430)
        right.pack_start(vpaned, True, True, 0)

        self.actions_box = Gtk.Box(spacing=6)
        right.pack_start(self.actions_box, False, False, 0)
        paned.pack2(right, True, False)

        self.window.show_all()
        self.Visible = True
        if self.p_session.psystem:
            self.combobox_systems.set_active_system(e_id=self.p_session.active_id)
        first = self.store.iter_children(self.store.get_iter_first())
        self.tree.get_selection().select_iter(first)
        self._set_result_widgets(None)

    def _build_input_frame(self):
        frame = Gtk.Frame(label=" Trajectory ")
        grid = Gtk.Grid(column_spacing=8, row_spacing=6, margin=8)
        frame.add(grid)
        self.combobox_systems = SystemComboBox(self.main)
        self.combobox_systems.connect("changed", self._on_system_changed)
        self.coordinates_combobox = CoordinatesComboBox()
        self.coordinates_combobox.connect("changed", self._on_object_changed)
        grid.attach(Gtk.Label(label="System", xalign=0), 0, 0, 1, 1)
        grid.attach(self.combobox_systems, 1, 0, 1, 1)
        grid.attach(Gtk.Label(label="Object", xalign=0), 2, 0, 1, 1)
        grid.attach(self.coordinates_combobox, 3, 0, 1, 1)
        self.label_frames = Gtk.Label(xalign=0)
        grid.attach(self.label_frames, 4, 0, 4, 1)

        self.spin_first = Gtk.SpinButton.new_with_range(0, 1e7, 1)
        self.spin_last = Gtk.SpinButton.new_with_range(-1, 1e7, 1)
        self.spin_last.set_value(-1)
        self.spin_last.set_tooltip_text("-1 = the last frame")
        self.spin_step = Gtk.SpinButton.new_with_range(1, 1e6, 1)
        self.entry_dt = Gtk.Entry(text="0", width_chars=7)
        self.entry_dt.set_tooltip_text("Time between frames, in ps. 0 = plot against the frame number.")
        for col, (text, w) in enumerate((("First frame", self.spin_first), ("Last", self.spin_last),
                                         ("Step", self.spin_step), ("Time per frame (ps)", self.entry_dt))):
            grid.attach(Gtk.Label(label=text, xalign=0), 2 * col, 1, 1, 1)
            grid.attach(w, 2 * col + 1, 1, 1, 1)

        self.radio_all = Gtk.RadioButton.new_with_label(None, "Whole object")
        self.radio_sel = Gtk.RadioButton.new_with_label_from_widget(self.radio_all, "Current selection")
        self.radio_sel.set_tooltip_text("The atoms selected in the 3D view (the 'atom set' of each analysis)")
        self.check_solvent = Gtk.CheckButton(label="Ignore water and ions")
        self.check_solvent.set_active(True)
        self.check_solvent.set_tooltip_text("Leave water and single-atom ions out of the atom set")
        grid.attach(Gtk.Label(label="Atom set", xalign=0), 0, 2, 1, 1)
        grid.attach(self.radio_all, 1, 2, 1, 1)
        grid.attach(self.radio_sel, 2, 2, 2, 1)
        grid.attach(self.check_solvent, 4, 2, 3, 1)
        return frame

    # ------------------------------------------------------------- options
    def _opt(self, key, name, widget):
        self.option_widgets.setdefault(key, {})[name] = widget
        return widget

    def _build_options(self, key):
        g = Gtk.Grid(column_spacing=8, row_spacing=6, margin=8)
        row = [0]

        def add(label, name, widget, tip=None):
            if tip:
                widget.set_tooltip_text(tip)
            g.attach(Gtk.Label(label=label, xalign=0), 0, row[0], 1, 1)
            g.attach(self._opt(key, name, widget), 1, row[0], 1, 1)
            row[0] += 1

        def combo(items, active=0):
            c = Gtk.ComboBoxText()
            for it in items:
                c.append_text(it)
            c.set_active(active)
            return c

        def spin(lo, hi, step, value, digits=2):
            s = Gtk.SpinButton.new_with_range(lo, hi, step)
            s.set_digits(digits)
            s.set_value(value)
            return s

        def check(label, active=True):
            c = Gtk.CheckButton(label=label)
            c.set_active(active)
            g.attach(self._opt(key, label, c), 0, row[0], 2, 1)
            row[0] += 1
            return c

        hint = ""
        if key == "rg":
            check("Mass weighted")
            hint = "Radius of gyration of the atom set in every frame."
        elif key == "dssp":
            hint = ("DSSP (Kabsch & Sander, 1983) for the protein residues of the atom set. Amide H are "
                    "placed by the DSSP rule when missing. Click the map to go to a frame.")
        elif key == "sasa":
            add("Probe radius (A)", "probe", spin(0.0, 3.0, 0.1, 1.4))
            add("Points per atom", "points", spin(20, 500, 10, 100, 0))
            hint = "Shrake-Rupley SASA of the atom set, alone (other atoms ignored)."
        elif key == "hbonds":
            add("Between", "mode", combo(["the atom set", "the atom set and the rest", "the atom set and water"]))
            add("D...A distance (A)", "dist", spin(2.5, 4.5, 0.1, 3.5))
            add("H-D...A angle (deg)", "angle", spin(10, 90, 5, 30, 0))
            add("Minimum occupancy (%)", "occ", spin(0, 100, 1, 5, 0))
            hint = ("GROMACS criterion. Double-click a row to select its atoms and go to a frame where it "
                    "exists.")
        elif key == "contacts":
            add("Map", "kind", combo(["Contact frequency", "Mean CA-CA distance", "Contact change (2nd half - 1st half)"]))
            add("Contact cutoff (A)", "cutoff", spin(2.5, 12.0, 0.5, 4.5))
            check("Heavy atoms only")
            hint = "Residues of the atom set. Click a cell to select the two residues."
        elif key == "energy":
            add("Partner group", "partner", combo(["the rest of the system", "the rest without water and ions",
                                                   "water and ions"]))
            add("Cutoff (A, 0 = none)", "cutoff", spin(0, 99, 1, 0, 0))
            hint = ("Coulomb + Lennard-Jones between the atom set and the partner group, with the MM charges "
                    "and LJ parameters of the system (kJ/mol). Needs a force field.")
        elif key == "et":
            self.combo_donor = Gtk.ComboBoxText()
            self.combo_acceptor = Gtk.ComboBoxText()
            add("Donor", "donor", self.combo_donor)
            add("Acceptor", "acceptor", self.combo_acceptor)
            add("Through-space cutoff (A)", "space", spin(3.0, 8.0, 0.5, 5.0, 1))
            hint = ("Pathways model (Beratan & Onuchic): best donor-to-acceptor path through covalent bonds (0.6), "
                    "hydrogen bonds and through-space jumps, in every frame. Donor and acceptor are saved "
                    "selections (Selection > Send to Selection List) or the current selection. The atom set "
                    "is the medium (water is used only when 'Ignore water and ions' is off).")
        elif key == "pca":
            add("Atoms", "subset", combo(ATOM_SUBSETS))
            check("Mass weighted", False)
            add("Components shown", "ncomp", spin(2, 50, 1, 10, 0))
            hint = "Frames superposed on their average first. 'Animate PC' creates an object along a component."
        elif key == "dccm":
            add("Atoms", "subset", combo(ATOM_SUBSETS))
            hint = "C_ij = <dr_i . dr_j> / sqrt(<dr_i^2><dr_j^2>) after superposition."
        elif key == "cluster":
            add("Atoms", "subset", combo(ATOM_SUBSETS))
            add("RMSD cutoff (A)", "cutoff", spin(0.1, 10.0, 0.1, 1.5))
            hint = "GROMOS clustering (Daura et al., 1999) on the pairwise RMSD after superposition."
        elif key == "bfactors":
            add("Atoms", "subset", combo(ATOM_SUBSETS, 3))
            add("Superpose on", "fit", combo(["CA atoms", "Backbone (N, CA, C, O)", "Heavy atoms", "Nothing (raw frames)"]))
            hint = ("RMSF, B = 8 pi^2/3 <dr^2> and the anisotropic U tensors (as pDynamo's "
                    "CoordinateFluctuations), plus the average structure.")
        elif key == "qh":
            add("Atoms", "subset", combo(ATOM_SUBSETS, 2))
            add("Temperature (K)", "temp", spin(1, 2000, 10, 300, 1))
            hint = ("Mass-weighted covariance -> quasi-harmonic frequencies and entropy (Andricioaei & "
                    "Karplus 2001) and Schlitter's bound. Needs many more frames than 3 x atoms.")
        elif key == "rotamers":
            hint = "chi1/chi2 of the protein residues of the atom set: populations g+, t, g- and transitions."
        elif key == "msd":
            add("Atoms", "which", combo(["Water oxygens", "Atom set"]))
            add("Fit from (% of lags)", "fit0", spin(0, 90, 5, 20, 0))
            add("Fit to (% of lags)", "fit1", spin(10, 100, 5, 80, 0))
            hint = ("MSD over time origins, periodic jumps removed; D from MSD = 6 D t. Set the time per "
                    "frame above.")
        elif key == "stats":
            self.combo_series = Gtk.ComboBoxText()
            add("Series", "series", self.combo_series)
            btn = Gtk.Button(label="Load from File...")
            btn.connect("clicked", self._on_load_series)
            g.attach(btn, 1, row[0], 1, 1)
            row[0] += 1
            add("Discard first (%)", "discard", spin(0, 90, 5, 0, 0))
            hint = ("Mean, standard deviation, block-averaged standard error (Flyvbjerg & Petersen 1989), "
                    "autocorrelation time and histogram of any series computed in this window.")
        g._hint = hint
        return g

    def _on_analysis_changed(self, selection):
        model, it = selection.get_selected()
        if it is None:
            return
        key = model[it][1]
        if not key:
            return
        self.current_key = key
        page = self.stack.get_child_by_name(key)
        self.stack.set_visible_child(page)
        self.label_hint.set_markup("<small>{}</small>".format(GLib.markup_escape_text(page._hint)))
        if key == "stats":
            self._refresh_series_combo()
        if key == "et":
            self._refresh_et_combos()

    def _ov(self, key, name):
        w = self.option_widgets[key][name]
        if isinstance(w, Gtk.SpinButton):
            return w.get_value()
        if isinstance(w, Gtk.ComboBoxText):
            return w.get_active()
        if isinstance(w, Gtk.CheckButton):
            return w.get_active()
        return w

    # ------------------------------------------------------------- inputs
    def _on_system_changed(self, _w):
        system_id = self.combobox_systems.get_system_id()
        if system_id is None:
            return
        self.coordinates_combobox.set_model(self.main.vobject_liststore_dict[system_id])
        size = len(list(self.main.vobject_liststore_dict[system_id]))
        self.coordinates_combobox.set_active(size - 1)
        if getattr(self, "current_key", None) == "et":
            self._refresh_et_combos()

    def _refresh_et_combos(self):
        system_id = self.combobox_systems.get_system_id()
        system = self.p_session.psystem.get(system_id) if system_id is not None else None
        names = ["Current selection"] + sorted((getattr(system, "e_selections", None) or {}).keys())
        for combo, default in ((self.combo_donor, 0), (self.combo_acceptor, min(1, len(names) - 1))):
            current = combo.get_active_text()
            combo.remove_all()
            for n in names:
                combo.append_text(n)
            combo.set_active(names.index(current) if current in names else default)

    def _on_object_changed(self, _w):
        v = self._vobject()
        if v is None:
            self.label_frames.set_text("")
            return
        n = len(v.frames)
        self.label_frames.set_text("{} frames, {} atoms".format(n, len(v.atoms)))
        self.spin_first.set_range(0, max(n - 1, 0))
        self.spin_last.set_range(-1, max(n - 1, 0))

    def _vobject(self):
        vid = self.coordinates_combobox.get_vobject_id()
        return None if vid is None else self.vm_session.vm_objects_dic.get(vid)

    def _context(self):
        """Everything the computation needs, read in the GTK thread."""
        v = self._vobject()
        if v is None:
            raise ValueError("Choose a System and an Object (trajectory) first.")
        n_frames = len(v.frames)
        first = int(self.spin_first.get_value())
        last = int(self.spin_last.get_value())
        last = n_frames - 1 if last < 0 else min(last, n_frames - 1)
        step = int(self.spin_step.get_value())
        frames_idx = list(range(first, last + 1, step))
        if not frames_idx:
            raise ValueError("No frames in the chosen range.")
        try:
            dt = float(self.entry_dt.get_text())
        except ValueError:
            dt = 0.0
        n = len(v.atoms)
        names, elements, resname, reskey, chain, is_protein, is_solvent = [], [], [], [], [], [], []
        res_ids, res_labels, res_objs = {}, [], []
        for i in range(n):
            a = v.atoms[i]
            r = a.residue
            names.append(a.name)
            elements.append(TA.element_of(a.symbol))
            rname = getattr(r, "name", "UNK")
            resname.append(rname)
            ch = getattr(getattr(a, "chain", None), "name", "") or ""
            chain.append(ch)
            if id(r) not in res_ids:
                res_ids[id(r)] = len(res_labels)
                res_labels.append((ch, rname, getattr(r, "index", 0)))
                res_objs.append(r)
            reskey.append(res_ids[id(r)])
            is_protein.append(bool(getattr(r, "is_protein", False)))
            n_res_atoms = len(getattr(r, "atoms", {})) if r is not None else 1
            is_solvent.append(rname.upper() in WATER_NAMES or (n_res_atoms == 1 and not is_protein[-1]))
        chains = {c for c, _, _ in res_labels}
        labels = ["{}{}{}".format(c + ":" if len(chains) > 1 and c else "", rn, ri) for c, rn, ri in res_labels]
        bonds = []
        for b in v.bonds.values() if isinstance(v.bonds, dict) else v.bonds:
            bonds.append((int(b.atom_index_i), int(b.atom_index_j)))
        if self.radio_sel.get_active():
            sel = self.vm_session.selections[self.vm_session.current_selection]
            atoms = sorted({a.atom_id for a in sel.selected_atoms if a.vm_object is v})
            if not atoms:
                raise ValueError("'Current selection' is chosen, but nothing of this object is selected.")
        else:
            atoms = list(range(n))
        if self.check_solvent.get_active():
            atoms = [i for i in atoms if not is_solvent[i]]
            if not atoms:
                raise ValueError("The atom set has only water and ions ('Ignore water and ions' is on).")
        boxes = None
        cc = getattr(v, "cell_coordinates", None)
        if cc is not None and len(cc):
            boxes = [np.asarray(cc[min(f, len(cc) - 1)][-1], dtype=float) for f in frames_idx]
        system = self.p_session.psystem.get(getattr(v, "e_id", None))
        ctx = {
            "vobject": v, "system": system, "frames_idx": frames_idx, "dt": dt, "atoms": atoms,
            "names": names, "elements": elements, "resname": resname, "reskey": np.array(reskey),
            "res_labels": labels, "res_objs": res_objs, "chain": chain, "is_protein": np.array(is_protein),
            "is_solvent": np.array(is_solvent), "bonds": bonds, "boxes": boxes, "key": self.current_key,
            "options": {name: self._ov(self.current_key, name) for name in self.option_widgets.get(self.current_key, {})},
        }
        if self.current_key == "et":
            groups = {}
            saved = getattr(system, "e_selections", None) or {}
            for role, combo in (("donor", self.combo_donor), ("acceptor", self.combo_acceptor)):
                name = combo.get_active_text()
                if name == "Current selection":
                    sel = self.vm_session.selections[self.vm_session.current_selection]
                    idx = sorted({a.atom_id for a in sel.selected_atoms if a.vm_object is v})
                else:
                    idx = sorted(int(i) for i in saved.get(name, []))
                if not idx:
                    raise ValueError("The {} group ({}) is empty.".format(role, name))
                groups[role] = idx
            if set(groups["donor"]) & set(groups["acceptor"]):
                raise ValueError("Donor and acceptor share atoms.")
            ctx["et_groups"] = groups
        if self.current_key == "stats":
            ctx["series_name"] = self.combo_series.get_active_text()
            ctx["series"] = self.series_registry.get(ctx["series_name"])
        return ctx

    # =====================================================================
    #  Running
    # =====================================================================
    def on_run(self, _button):
        if self._thread is not None:
            return
        try:
            ctx = self._context()
        except ValueError as error:
            self._message(str(error))
            return
        self._cancel = False
        self._progress = 0.0
        self.btn_run.set_sensitive(False)
        self.btn_cancel.set_sensitive(True)
        self.progress.set_fraction(0.0)
        self.progress.set_text("running...")
        self._t0 = time.time()
        self._thread = threading.Thread(target=self._worker, args=(ctx,), daemon=True)
        self._thread.start()
        GLib.timeout_add(150, self._poll_progress)

    def _report(self, fraction):
        if self._cancel:
            raise _Cancelled()
        self._progress = fraction

    def _worker(self, ctx):
        try:
            result = getattr(self, "_compute_" + ctx["key"])(ctx)
            result["ctx"] = ctx
            GLib.idle_add(self._done, result, None)
        except _Cancelled:
            GLib.idle_add(self._done, None, "cancelled")
        except Exception as error:
            traceback.print_exc()
            GLib.idle_add(self._done, None, str(error))

    def _poll_progress(self):
        if self._thread is None:
            return False
        self.progress.set_fraction(min(self._progress, 1.0))
        self.progress.set_text("running... {:.0f} s".format(time.time() - self._t0))
        return True

    def _done(self, result, error):
        self._thread = None
        self.btn_run.set_sensitive(True)
        self.btn_cancel.set_sensitive(False)
        if error:
            self.progress.set_fraction(0.0)
            self.progress.set_text(error if error == "cancelled" else "failed")
            if error != "cancelled":
                self._message(error)
            return False
        self.progress.set_fraction(1.0)
        self.progress.set_text("done in {:.1f} s".format(time.time() - self._t0))
        for name, s in result.get("series", {}).items():
            self.series_registry[name] = s
        self._set_result_widgets(result)
        return False

    # ------------------------------------------------------------- helpers
    def _xs(self, ctx):
        dt = ctx["dt"]
        if dt > 0:
            return [f * dt for f in ctx["frames_idx"]], "time (ps)"
        return list(ctx["frames_idx"]), "frame"

    def _coords(self, ctx, atoms):
        v = ctx["vobject"]
        return np.asarray(v.frames[np.ix_(ctx["frames_idx"], atoms)], dtype=float)

    def _subset(self, ctx, which, base=None):
        base = ctx["atoms"] if base is None else base
        names, el = ctx["names"], ctx["elements"]
        prot = ctx["is_protein"]
        if which == 0:
            out = [i for i in base if names[i] == "CA" and prot[i]]
        elif which == 1:
            out = [i for i in base if names[i] in ("N", "CA", "C", "O") and prot[i]]
        elif which == 2:
            out = [i for i in base if el[i] != "H"]
        else:
            out = list(base)
        if not out:
            raise ValueError("No atoms of the chosen kind in the atom set ({}).".format(
                ["CA", "backbone", "heavy", "any"][which]))
        return out

    def _atom_label(self, ctx, i):
        return "{}:{}".format(ctx["res_labels"][ctx["reskey"][i]], ctx["names"][i])

    def _residues_of(self, ctx, atoms):
        order, members = [], {}
        for i in atoms:
            r = int(ctx["reskey"][i])
            if r not in members:
                members[r] = []
                order.append(r)
            members[r].append(i)
        return order, members

    # =====================================================================
    #  Analyses
    # =====================================================================
    def _compute_rg(self, ctx):
        atoms = ctx["atoms"]
        X = self._coords(ctx, atoms)
        masses = TA.masses_of([ctx["elements"][i] for i in atoms]) if ctx["options"]["Mass weighted"] else None
        rg = TA.radius_of_gyration(X, masses)
        x, xl = self._xs(ctx)
        d = TA.describe(rg)
        return {
            "title": "Radius of gyration",
            "views": [("Rg vs time", lambda p: p.set_lines([{"x": x, "y": rg, "label": "Rg"}], xl, "Rg (A)",
                                                          "Radius of gyration, {} atoms".format(len(atoms))))],
            "frames_x": x,
            "table": (["frame", xl, "Rg (A)"], [[f, xx, "{:.3f}".format(v)] for f, xx, v in zip(ctx["frames_idx"], x, rg)]),
            "summary": "Rg = {:.3f} +/- {:.3f} A (SEM, block averaging {:.3f}); min {:.3f}, max {:.3f}; {} atoms, {} frames.".format(
                d["mean"], d["std"], d["sem_block"], d["min"], d["max"], len(atoms), len(rg)),
            "series": {"Rg (A)": (x, list(rg), xl, "Rg (A)")},
            "data": (xl + "\tRg", np.column_stack([x, rg])),
        }

    def _backbone(self, ctx):
        order, members = self._residues_of(ctx, [i for i in ctx["atoms"] if ctx["is_protein"][i]])
        bb, keys = [], []
        for r in order:
            names = {ctx["names"][i]: i for i in members[r]}
            if all(k in names for k in ("N", "CA", "C", "O")):
                entry = {"N": names["N"], "CA": names["CA"], "C": names["C"], "O": names["O"],
                         "res": str(ctx["resname"][names["CA"]])}
                h = names.get("H", names.get("HN"))
                if h is not None:
                    entry["H"] = h
                bb.append(entry)
                keys.append(r)
        if len(bb) < 3:
            raise ValueError("The atom set has fewer than three protein residues with N, CA, C and O.")
        return bb, keys

    def _compute_dssp(self, ctx):
        bb, keys = self._backbone(ctx)
        need = sorted({i for b in bb for i in b.values() if isinstance(i, (int, np.integer))})
        local = {g: k for k, g in enumerate(need)}
        X = self._coords(ctx, need)
        lbb = [{k: (local[val] if isinstance(val, (int, np.integer)) else val) for k, val in b.items()} for b in bb]
        X0 = X[0]
        breaks = [k for k in range(len(bb) - 1)
                  if ctx["chain"][bb[k]["CA"]] != ctx["chain"][bb[k + 1]["CA"]]
                  or np.linalg.norm(X0[lbb[k]["C"]] - X0[lbb[k + 1]["N"]]) > 2.5]
        codes = TA.dssp(X, lbb, breaks, progress=self._report)
        F, R = codes.shape
        x, xl = self._xs(ctx)
        frac = {name: [] for name in ("helix (H, G, I)", "strand (E, B)", "turn/bend (T, S)", "coil")}
        for f in range(F):
            row = codes[f]
            frac["helix (H, G, I)"].append(np.isin(row, ["H", "G", "I"]).mean() * 100)
            frac["strand (E, B)"].append(np.isin(row, ["E", "B"]).mean() * 100)
            frac["turn/bend (T, S)"].append(np.isin(row, ["T", "S"]).mean() * 100)
            frac["coil"].append((row == "-").mean() * 100)
        stride = max(1, F // 600)
        shown = list(range(0, F, stride))
        rl = [ctx["res_labels"][r] for r in keys]
        mat = [[codes[f, r] for f in shown] for r in range(R)]
        rows = []
        for r in range(R):
            col = codes[:, r]
            pct = {c: (col == c).mean() * 100 for c in "HGIEBTS-"}
            most = max(pct, key=pct.get)
            rows.append([rl[r]] + ["{:.0f}".format(pct[c]) for c in "HGIEBTS-"] + [TA.DSSP_NAMES[most]])
        series = [{"x": x, "y": frac[k], "label": k, "color": PALETTE[i]} for i, k in enumerate(frac)]
        summary = "{} residues, {} frames. Mean content: ".format(R, F) + ", ".join(
            "{} {:.1f}%".format(k, np.mean(vv)) for k, vv in frac.items())
        return {
            "title": "Secondary structure",
            "views": [("DSSP map", lambda p: p.set_categories(mat, DSSP_COLORS, TA.DSSP_NAMES, xl, "residue",
                                                              "Secondary structure (DSSP)", row_labels=rl,
                                                              col_labels=[_short(x[f]) for f in shown])),
                      ("Content vs time", lambda p: p.set_lines(series, xl, "% of residues", "Secondary structure content", ymin=0))],
            "frames_x": x, "map_frames": [ctx["frames_idx"][f] for f in shown],
            "table": (["residue", "H %", "G %", "I %", "E %", "B %", "T %", "S %", "coil %", "most frequent"], rows),
            "summary": summary,
            "series": {"helix content (%)": (x, frac["helix (H, G, I)"], xl, "% helix"),
                       "strand content (%)": (x, frac["strand (E, B)"], xl, "% strand")},
            "data": ("residue\t" + "\t".join(_short(v) for v in x),
                     [[rl[r]] + list(codes[:, r]) for r in range(R)]),
        }

    def _compute_sasa(self, ctx):
        atoms = ctx["atoms"]
        X = self._coords(ctx, atoms)
        radii = np.array([TA.VDW_RADII.get(ctx["elements"][i], 1.8) for i in atoms])
        probe, npts = ctx["options"]["probe"], int(ctx["options"]["points"])
        per_atom = np.zeros((len(X), len(atoms)))
        for f in range(len(X)):
            per_atom[f] = TA.sasa(X[f], radii, probe, npts)
            self._report((f + 1) / len(X))
        total = per_atom.sum(axis=1)
        order, members = self._residues_of(ctx, atoms)
        pos = {g: k for k, g in enumerate(atoms)}
        res_mean = [per_atom[:, [pos[i] for i in members[r]]].sum(axis=1).mean() for r in order]
        rl = [ctx["res_labels"][r] for r in order]
        rel = []
        for r, val in zip(order, res_mean):
            name = str(ctx["resname"][members[r][0]]).upper()[:3]
            rel.append(val / MAX_ASA[name] * 100 if name in MAX_ASA else float("nan"))
        x, xl = self._xs(ctx)
        d = TA.describe(total)
        return {
            "title": "SASA",
            "views": [("Total SASA vs time", lambda p: p.set_lines([{"x": x, "y": total, "label": "SASA"}], xl, "SASA (A^2)",
                                                                  "Solvent accessible surface area")),
                      ("SASA per residue", lambda p: p.set_lines([{"x": list(range(len(rl))), "y": res_mean, "style": "bars"}],
                                                                "residue", "mean SASA (A^2)", "SASA per residue", ymin=0, xlabels=rl)),
                      ("Relative exposure", lambda p: p.set_lines([{"x": list(range(len(rl))), "y": rel, "style": "bars",
                                                                    "color": PALETTE[2]}], "residue", "relative SASA (%)",
                                                                  "Relative exposure (Tien et al. 2013 maxima)", ymin=0, xlabels=rl))],
            "frames_x": x,
            "table": (["residue", "mean SASA (A^2)", "relative (%)"],
                      [[l, "{:.1f}".format(a), "" if b != b else "{:.0f}".format(b)] for l, a, b in zip(rl, res_mean, rel)]),
            "summary": "Total SASA {:.0f} +/- {:.0f} A^2 (SEM {:.1f}); probe {:.2f} A, {} points per atom, {} atoms.".format(
                d["mean"], d["std"], d["sem_block"], probe, npts, len(atoms)),
            "series": {"SASA (A^2)": (x, list(total), xl, "SASA (A^2)")},
            "data": (xl + "\tSASA", np.column_stack([x, total])),
            "residues": order,
        }

    def _compute_hbonds(self, ctx):
        n = len(ctx["names"])
        mode = ctx["options"]["mode"]
        aset = set(ctx["atoms"])
        if mode == 0:
            groups = [(ctx["atoms"], ctx["atoms"])]
        else:
            if mode == 1:
                rest = [i for i in range(n) if i not in aset]
            else:
                rest = [i for i in range(n) if i not in aset and ctx["is_solvent"][i]]
            groups = [(ctx["atoms"], rest), (rest, ctx["atoms"])]
        all_atoms = sorted(set(ctx["atoms"]) | {i for g in groups for i in g[1]})
        pos = {g: k for k, g in enumerate(all_atoms)}
        X = _LazyFrames(ctx["vobject"], ctx["frames_idx"], all_atoms)
        lbonds = [(pos[i], pos[j]) for i, j in ctx["bonds"] if i in pos and j in pos]
        lel = [ctx["elements"][i] for i in all_atoms]
        counts = np.zeros(len(X), dtype=int)
        present, geometry = {}, {}
        for k, (dgroup, agroup) in enumerate(groups):
            dset = [pos[i] for i in dgroup]
            aset_ = [pos[i] for i in agroup]
            donors, _ = TA.donors_and_acceptors(lel, lbonds, dset)
            _, acceptors = TA.donors_and_acceptors(lel, lbonds, aset_)
            if not donors or not acceptors:
                continue
            c, p, g = TA.hydrogen_bonds(X, donors, acceptors, ctx["boxes"], ctx["options"]["dist"],
                                        ctx["options"]["angle"],
                                        progress=lambda fr, k=k: self._report((k + fr) / len(groups)))
            counts += c
            present.update(p)
            geometry.update(g)
        F = len(X)
        rows = []
        for key, frames in sorted(present.items(), key=lambda kv: -len(kv[1])):
            occ = len(frames) / F * 100
            if occ < ctx["options"]["occ"]:
                continue
            d, h, a = (all_atoms[i] for i in key)
            dist, ang = geometry[key]
            rows.append([self._atom_label(ctx, d), ctx["names"][h], self._atom_label(ctx, a), "{:.1f}".format(occ),
                         "{:.2f}".format(dist), "{:.1f}".format(ang), str(max(TA.lifetimes(frames))),
                         (d, h, a, ctx["frames_idx"][frames[0]])])
        x, xl = self._xs(ctx)
        dsc = TA.describe(counts) if F > 1 else {"mean": counts.mean(), "std": 0}
        top = rows[:25]
        return {
            "title": "Hydrogen bonds",
            "views": [("Number vs time", lambda p: p.set_lines([{"x": x, "y": counts, "label": "H-bonds"}], xl, "hydrogen bonds",
                                                              "Hydrogen bonds per frame", ymin=0)),
                      ("Occupancy (top 25)", lambda p: p.set_lines([{"x": list(range(len(top))), "y": [float(r[3]) for r in top],
                                                                     "style": "bars"}], "hydrogen bond", "occupancy (%)",
                                                                   "Most persistent hydrogen bonds", ymin=0,
                                                                   xlabels=["{}>{}".format(r[0], r[2]) for r in top]))],
            "frames_x": x,
            "table": (["donor", "H", "acceptor", "occupancy %", "D...A (A)", "angle (deg)", "longest (frames)"], rows),
            "summary": "{:.1f} +/- {:.1f} hydrogen bonds per frame; {} distinct bonds, {} with occupancy >= {:.0f}%.".format(
                dsc["mean"], dsc["std"], len(present), len(rows), ctx["options"]["occ"]),
            "series": {"H-bond count": (x, list(counts), xl, "hydrogen bonds")},
            "data": (xl + "\tcount", np.column_stack([x, counts])),
            "row_action": "hbond",
        }

    def _compute_contacts(self, ctx):
        kind = ctx["options"]["kind"]
        atoms = ctx["atoms"]
        if ctx["options"]["Heavy atoms only"]:
            atoms = [i for i in atoms if ctx["elements"][i] != "H"]
        order, members = self._residues_of(ctx, atoms)
        rl = [ctx["res_labels"][r] for r in order]
        R = len(order)
        if R < 2:
            raise ValueError("The atom set has fewer than two residues.")
        if kind == 1:
            ca = []
            for r in order:
                c = [i for i in members[r] if ctx["names"][i] == "CA"]
                if not c:
                    raise ValueError("Residue {} has no CA atom.".format(ctx["res_labels"][r]))
                ca.append(c[0])
            X = self._coords(ctx, ca)
            mean, std = TA.ca_distance_map(X)
            M, label, div, lo, hi = mean, "A", False, None, None
            title = "Mean CA-CA distance"
            pairs = [(mean[i, j], i, j, std[i, j]) for i in range(R) for j in range(i + 1, R)]
            pairs.sort()
            rows = [[rl[i], rl[j], "{:.2f}".format(d), "{:.2f}".format(s), (order[i], order[j])] for d, i, j, s in pairs[:300]]
            cols = ["residue i", "residue j", "mean d (A)", "std (A)"]
        else:
            loc = {g: k for k, g in enumerate(atoms)}
            res_of = np.zeros(len(atoms), dtype=int)
            for k, r in enumerate(order):
                for i in members[r]:
                    res_of[loc[i]] = k
            X = self._coords(ctx, atoms)
            boxes = ctx["boxes"]
            if kind == 0:
                M = TA.residue_contact_map(X, res_of, R, ctx["options"]["cutoff"], boxes, progress=self._report)
                label, div, lo, hi, title = "frequency", False, 0.0, 1.0, "Contact frequency ({:.1f} A)".format(ctx["options"]["cutoff"])
            else:
                h = len(X) // 2
                if h < 1:
                    raise ValueError("Needs at least two frames.")
                m1 = TA.residue_contact_map(X[:h], res_of, R, ctx["options"]["cutoff"], boxes[:h] if boxes else None,
                                            progress=lambda f: self._report(f / 2))
                m2 = TA.residue_contact_map(X[h:], res_of, R, ctx["options"]["cutoff"], boxes[h:] if boxes else None,
                                            progress=lambda f: self._report(0.5 + f / 2))
                M = m2 - m1
                label, div, lo, hi, title = "change", True, -1.0, 1.0, "Contact change: 2nd half - 1st half"
            pairs = [(M[i, j], i, j) for i in range(R) for j in range(i + 1, R) if abs(j - i) >= 3 and M[i, j] != 0]
            pairs.sort(key=lambda t: -abs(t[0]))
            rows = [[rl[i], rl[j], "{:.2f}".format(val), (order[i], order[j])] for val, i, j in pairs[:500]]
            cols = ["residue i", "residue j", label + " (|i-j| >= 3)"]
        mlist = M.tolist()
        return {
            "title": title,
            "views": [(title, lambda p: p.set_heatmap(mlist, "residue", "residue", title, lo, hi, div, rl, rl, label))],
            "table": (cols, rows),
            "summary": "{} residues, {} frames.".format(R, len(ctx["frames_idx"])),
            "data": ("matrix\t" + "\t".join(rl), [[rl[i]] + ["{:.4f}".format(v) for v in M[i]] for i in range(R)]),
            "cell_action": order, "row_action": "residues",
        }

    def _compute_energy(self, ctx):
        system = ctx["system"]
        v = ctx["vobject"]
        if system is None or getattr(system, "mmState", None) is None:
            raise ValueError("The system has no force field (MM model).")
        n = len(v.atoms)
        if len(system.atoms) != n:
            raise ValueError("The object and the system have different numbers of atoms.")
        mm = system.mmState
        A = list(ctx["atoms"])
        aset = set(A)
        partner = ctx["options"]["partner"]
        if partner == 0:
            B = [i for i in range(n) if i not in aset]
        elif partner == 1:
            B = [i for i in range(n) if i not in aset and not ctx["is_solvent"][i]]
        else:
            B = [i for i in range(n) if i not in aset and ctx["is_solvent"][i]]
        if not B:
            raise ValueError("The partner group is empty.")
        q = np.array([mm.charges[i] for i in range(n)])
        lj = lj14 = None
        form = "none"
        if mm.ljParameters is not None and mm.ljTypeIndices is not None:
            data = mm.ljParameters.__getstate__()
            if "epsilons" in data:
                form = str(data.get("analyticForm", "OPLS")).split(".")[-1]
                t = [int(mm.ljTypeIndices[i]) for i in range(n)]
                lj = TA.lj_pair_function(form, [data["epsilons"][k] for k in t], [data["sigmas"][k] for k in t])
                if mm.ljParameters14 is not None:
                    d14 = mm.ljParameters14.__getstate__()
                    if "epsilons" in d14:
                        lj14 = TA.lj_pair_function(form, [d14["epsilons"][k] for k in t], [d14["sigmas"][k] for k in t])
        bset = set(B)
        linked = any((i in aset and j in bset) or (j in aset and i in bset) for i, j in ctx["bonds"])
        excluded, scaled = set(), set()
        if linked:
            for i, j in mm.exclusions:
                if (i in aset and j in bset) or (j in aset and i in bset):
                    excluded.add((min(i, j), max(i, j)))
            if mm.interactions14 is not None:
                for i, j in mm.interactions14:
                    if (i in aset and j in bset) or (j in aset and i in bset):
                        scaled.add((min(i, j), max(i, j)))
        scale14 = getattr(system.mmModel, "electrostaticScale14", 0.5)
        X = _LazyFrames(v, ctx["frames_idx"])
        groups = ctx["reskey"]
        cutoff = ctx["options"]["cutoff"] or None
        el, ljv, per = TA.interaction_energy(X, A, B, q, lj, excluded, scaled, scale14, lj14, ctx["boxes"] if cutoff else None,
                                             cutoff, groups[B], progress=self._report)
        tot = el + ljv
        x, xl = self._xs(ctx)
        top = sorted(per.items(), key=lambda kv: kv[1])[:25]
        rl = ctx["res_labels"]
        d = TA.describe(tot)
        rows = [[rl[g], "{:.2f}".format(val), (g,)] for g, val in sorted(per.items(), key=lambda kv: kv[1])]
        return {
            "title": "Interaction energy",
            "views": [("Energy vs time", lambda p: p.set_lines([{"x": x, "y": tot, "label": "total"},
                                                               {"x": x, "y": el, "label": "Coulomb"},
                                                               {"x": x, "y": ljv, "label": "Lennard-Jones"}],
                                                              xl, "E (kJ/mol)", "Interaction energy", hline=0.0)),
                      ("Per residue (25 strongest)", lambda p: p.set_lines([{"x": list(range(len(top))), "y": [t[1] for t in top],
                                                                             "style": "bars", "color": PALETTE[3]}],
                                                                           "residue", "mean E (kJ/mol)", "Residue decomposition",
                                                                           xlabels=[rl[t[0]] for t in top], hline=0.0))],
            "frames_x": x,
            "table": (["residue", "mean E (kJ/mol)"], rows),
            "summary": ("E = {:.1f} +/- {:.1f} kJ/mol (SEM {:.1f}): Coulomb {:.1f}, Lennard-Jones {:.1f}. {} vs {} atoms, "
                        "LJ form {}, {}{}.").format(d["mean"], d["std"], d["sem_block"], el.mean(), ljv.mean(), len(A), len(B),
                                                     form, "no cutoff" if not cutoff else "cutoff {:.0f} A".format(cutoff),
                                                     ", covalent exclusions applied" if linked else ""),
            "series": {"E interaction (kJ/mol)": (x, list(tot), xl, "E (kJ/mol)"),
                       "E Coulomb (kJ/mol)": (x, list(el), xl, "E (kJ/mol)"),
                       "E LJ (kJ/mol)": (x, list(ljv), xl, "E (kJ/mol)")},
            "data": (xl + "\ttotal\tcoulomb\tlj", np.column_stack([x, tot, el, ljv])),
            "row_action": "residues",
        }

    def _compute_et(self, ctx):
        from util import et_pathways as ET
        g = ctx["et_groups"]
        atoms = sorted(set(ctx["atoms"]) | set(g["donor"]) | set(g["acceptor"]))
        pos = {a: k for k, a in enumerate(atoms)}
        bonds = [(pos[i], pos[j]) for i, j in ctx["bonds"] if i in pos and j in pos]
        X = self._coords(ctx, atoms)
        res = ET.pathways(X, [ctx["elements"][i] for i in atoms], bonds, [pos[i] for i in g["donor"]],
                          [pos[i] for i in g["acceptor"]], ctx["boxes"], ctx["options"]["space"], progress=self._report)
        x, xl = self._xs(ctx)
        logc = [math.log10(c) if c > 0 else float("nan") for c, _a, _k in res]
        paths = {}
        for f, (c, path, kinds) in enumerate(res):
            if not path:
                continue
            key = tuple(atoms[k] for k in path)
            entry = paths.setdefault(key, {"frames": [], "couplings": [], "kinds": kinds})
            entry["frames"].append(ctx["frames_idx"][f])
            entry["couplings"].append(c)
        if not paths:
            raise ValueError("No path between donor and acceptor (increase the through-space cutoff).")
        rows = []
        F = len(res)
        for key, e in sorted(paths.items(), key=lambda kv: -len(kv[1]["frames"])):
            k = e["kinds"]
            text = " > ".join(self._atom_label(ctx, a) for a in key)
            best = max(e["couplings"])
            rows.append(["{:.1f}".format(len(e["frames"]) / F * 100), "{:.2f}".format(np.mean(np.log10(e["couplings"]))),
                         "{:.1f}".format(ET.sigma_length(best)), "{}C {}H {}S".format(k.count("C"), k.count("H"), k.count("S")),
                         text, (list(key), e["frames"][int(np.argmax(e["couplings"]))])])
        valid = [c for c, _a, _k in res if c > 0]
        mean_log = float(np.mean(np.log10(valid)))
        top = rows[0]
        return {
            "title": "ET pathways",
            "views": [("Coupling vs time", lambda p: p.set_lines([{"x": x, "y": logc, "label": "log10 coupling"}], xl,
                                                                "log10 (product of decays)", "Best pathway coupling"))],
            "frames_x": x,
            "table": (["frames %", "mean log10 coupling", "sigma (A)", "steps", "path"], rows),
            "summary": ("{} distinct best paths over {} frames; mean log10 coupling {:.2f} (sigma {:.1f} A). Most frequent path "
                        "({}% of the frames): {}. Double-click a row to select the path.").format(
                len(paths), F, mean_log, ET.sigma_length(10 ** mean_log), top[0], top[4]),
            "series": {"ET log10 coupling": (x, logc, xl, "log10 coupling")},
            "data": (xl + "\tlog10_coupling", np.column_stack([x, logc])),
            "row_action": "path",
        }

    def _aligned_subset(self, ctx, which, mass=False):
        sub = self._subset(ctx, which)
        X = self._coords(ctx, sub)
        if len(X) < 3:
            raise ValueError("Needs at least three frames.")
        w = TA.masses_of([ctx["elements"][i] for i in sub]) if mass else None
        aligned, _avg = TA.iterative_average(X, None, 3, w)
        return sub, aligned

    def _compute_pca(self, ctx):
        mass = ctx["options"]["Mass weighted"]
        sub, aligned = self._aligned_subset(ctx, ctx["options"]["subset"], mass)
        masses = TA.masses_of([ctx["elements"][i] for i in sub]) if mass else None
        p = TA.pca(aligned, masses)
        k = min(int(ctx["options"]["ncomp"]), len(p["eigenvalues"]))
        frac = p["fraction"][:k] * 100
        cum = np.cumsum(p["fraction"])[:k] * 100
        x, xl = self._xs(ctx)
        unit = "amu A^2" if mass else "A^2"
        pc1 = (p["eigenvectors"][:, 0].reshape(-1, 3) ** 2).sum(axis=1)
        labels = [self._atom_label(ctx, i) for i in sub]
        series = {}
        for c in range(min(k, 3)):
            series["PC{} projection".format(c + 1)] = (x, list(p["projections"][:, c]), xl, "PC{} ({})".format(c + 1, unit))
        return {
            "title": "PCA",
            "views": [("Variance per component", lambda pl: pl.set_lines([{"x": list(range(1, k + 1)), "y": frac, "style": "bars", "label": "component"},
                                                                         {"x": list(range(1, k + 1)), "y": cum, "label": "cumulative", "color": PALETTE[1]}],
                                                                        "component", "variance (%)", "Principal components", ymin=0)),
                      ("PC1 vs PC2", lambda pl: pl.set_scatter(p["projections"][:, 0], p["projections"][:, 1], x,
                                                               "PC1", "PC2", "Projection on PC1 and PC2", xl)),
                      ("PC1 vs time", lambda pl: pl.set_lines([{"x": x, "y": p["projections"][:, 0], "label": "PC1"},
                                                               {"x": x, "y": p["projections"][:, 1], "label": "PC2"}],
                                                              xl, "projection", "Projections", hline=0.0)),
                      ("PC1 per atom", lambda pl: pl.set_lines([{"x": list(range(len(sub))), "y": pc1, "style": "bars"}],
                                                               "atom", "PC1 weight", "Contribution to PC1", ymin=0, xlabels=labels))],
            "frames_x": x, "scatter_frames": ctx["frames_idx"],
            "table": (["PC", "eigenvalue (" + unit + ")", "variance %", "cumulative %"],
                      [[c + 1, "{:.4f}".format(p["eigenvalues"][c]), "{:.2f}".format(frac[c]), "{:.2f}".format(cum[c])] for c in range(k)]),
            "summary": "{} atoms, {} frames. PC1 {:.1f}%, PC1-3 {:.1f}% of the variance.".format(
                len(sub), len(aligned), frac[0], cum[min(2, k - 1)]),
            "series": series,
            "data": (xl + "\t" + "\t".join("PC{}".format(c + 1) for c in range(k)), np.column_stack([x, p["projections"][:, :k]])),
            "actions": [("Animate PC1", lambda: self._animate_pc(ctx, sub, p, 0, masses)),
                        ("Animate PC2", lambda: self._animate_pc(ctx, sub, p, 1, masses)),
                        ("Color by PC1 weight", lambda: self._color_by(ctx, sub, pc1))],
        }

    def _compute_dccm(self, ctx):
        sub, aligned = self._aligned_subset(ctx, ctx["options"]["subset"])
        C = TA.dccm(aligned)
        labels = [self._atom_label(ctx, i) if ctx["names"][i] != "CA" else ctx["res_labels"][ctx["reskey"][i]] for i in sub]
        R = len(sub)
        pairs = sorted(((C[i, j], i, j) for i in range(R) for j in range(i + 3, R)), key=lambda t: -abs(t[0]))[:300]
        Cl = C.tolist()
        return {
            "title": "DCCM",
            "views": [("Cross-correlation", lambda p: p.set_heatmap(Cl, "atom", "atom", "Dynamic cross-correlation", -1.0, 1.0,
                                                                   True, labels, labels, "C_ij"))],
            "table": (["atom i", "atom j", "C_ij"], [[labels[i], labels[j], "{:+.3f}".format(c), (int(ctx["reskey"][sub[i]]), int(ctx["reskey"][sub[j]]))] for c, i, j in pairs]),
            "summary": "{} atoms, {} frames. Positive: in-phase motion, negative: opposite.".format(R, len(aligned)),
            "data": ("C\t" + "\t".join(labels), [[labels[i]] + ["{:.4f}".format(v) for v in C[i]] for i in range(R)]),
            "row_action": "residues", "cell_action": [int(ctx["reskey"][i]) for i in sub],
        }

    def _compute_cluster(self, ctx):
        sub = self._subset(ctx, ctx["options"]["subset"])
        X = self._coords(ctx, sub)
        F = len(X)
        if F > 4000:
            raise ValueError("{} frames: use a frame step to keep at most 4000 for clustering.".format(F))
        D = TA.rmsd_matrix(X, progress=self._report)
        labels, centres = TA.gromos_clusters(D, ctx["options"]["cutoff"])
        x, xl = self._xs(ctx)
        sizes = [int((labels == c).sum()) for c in range(len(centres))]
        rows = []
        for c, centre in enumerate(centres):
            m = np.nonzero(labels == c)[0]
            intra = D[np.ix_(m, m)][np.triu_indices(len(m), 1)].mean() if len(m) > 1 else 0.0
            rows.append([c + 1, sizes[c], "{:.1f}".format(sizes[c] / F * 100), ctx["frames_idx"][centre],
                         "{:.2f}".format(intra), (ctx["frames_idx"][centre],)])
        stride = max(1, F // 300)
        shown = list(range(0, F, stride))
        Dl = [[D[i, j] for j in shown] for i in shown]
        return {
            "title": "Clustering",
            "views": [("Cluster vs time", lambda p: p.set_lines([{"x": x, "y": labels + 1, "style": "points"}], xl, "cluster",
                                                                "Cluster membership", ymin=0)),
                      ("Cluster sizes", lambda p: p.set_lines([{"x": list(range(1, len(sizes) + 1)), "y": sizes, "style": "bars"}],
                                                              "cluster", "frames", "Cluster populations", ymin=0)),
                      ("RMSD matrix", lambda p: p.set_heatmap(Dl, xl, xl, "Pairwise RMSD", 0.0, None, False,
                                                              [_short(x[i]) for i in shown], [_short(x[i]) for i in shown], "RMSD (A)"))],
            "frames_x": x,
            "table": (["cluster", "frames", "%", "representative frame", "mean RMSD in cluster (A)"], rows),
            "summary": "{} clusters at {:.2f} A over {} atoms and {} frames; the 3 largest hold {:.0f}% of the frames.".format(
                len(centres), ctx["options"]["cutoff"], len(sub), F, sum(sizes[:3]) / F * 100),
            "series": {"cluster id": (x, list(labels + 1), xl, "cluster")},
            "data": (xl + "\tcluster", np.column_stack([x, labels + 1])),
            "row_action": "frame",
            "actions": [("Create Object with Representatives", lambda: self._new_object(
                ctx, [ctx["frames_idx"][c] for c in centres], "clusters"))],
        }

    def _compute_bfactors(self, ctx):
        sub = self._subset(ctx, ctx["options"]["subset"])
        fit = ctx["options"]["fit"]
        X = self._coords(ctx, sub)
        if fit < 3:
            fit_atoms = self._subset(ctx, fit)
            pos = {g: k for k, g in enumerate(sub)}
            missing = [i for i in fit_atoms if i not in pos]
            if missing:
                fit_X = self._coords(ctx, fit_atoms)
                aligned_fit, avg = TA.iterative_average(fit_X)
                # . apply the same transforms to the analysed atoms
                X = _apply_fits(X, fit_X, avg)
            else:
                X, _ = TA.iterative_average(X, [pos[i] for i in fit_atoms])
        fl = TA.fluctuations(X)
        order, members = self._residues_of(ctx, sub)
        pos = {g: k for k, g in enumerate(sub)}
        res_b = [fl["b_factors"][[pos[i] for i in members[r]]].mean() for r in order]
        rl = [ctx["res_labels"][r] for r in order]
        rows = [[self._atom_label(ctx, i), "{:.3f}".format(fl["rmsf"][k]), "{:.2f}".format(fl["b_factors"][k])]
                for k, i in enumerate(sub)]
        self._last_fluct = (sub, fl)
        return {
            "title": "B-factors",
            "views": [("B per residue", lambda p: p.set_lines([{"x": list(range(len(rl))), "y": res_b, "style": "bars"}],
                                                             "residue", "B (A^2)", "B-factors per residue (mean)", ymin=0, xlabels=rl)),
                      ("RMSF per atom", lambda p: p.set_lines([{"x": list(range(len(sub))), "y": fl["rmsf"]}], "atom", "RMSF (A)",
                                                              "RMSF per atom", ymin=0))],
            "table": (["atom", "RMSF (A)", "B (A^2)"], rows),
            "summary": "Mean B {:.2f} A^2, mean RMSF {:.3f} A over {} atoms and {} frames ({}).".format(
                fl["b_factors"].mean(), fl["rmsf"].mean(), len(sub), len(X),
                ["fitted on CA", "fitted on the backbone", "fitted on heavy atoms", "raw frames"][fit]),
            "data": ("atom\tRMSF\tB", [[self._atom_label(ctx, i), "{:.4f}".format(fl["rmsf"][k]), "{:.3f}".format(fl["b_factors"][k])]
                                       for k, i in enumerate(sub)]),
            "actions": [("Create Average Structure", lambda: self._average_object(ctx, sub, fl)),
                        ("Color by B", lambda: self._color_by(ctx, sub, fl["b_factors"])),
                        ("Export PDB with B and ANISOU...", lambda: self._export_bfactor_pdb(ctx, sub, fl))],
        }

    def _compute_qh(self, ctx):
        sub, aligned = self._aligned_subset(ctx, ctx["options"]["subset"], mass=True)
        masses = TA.masses_of([ctx["elements"][i] for i in sub])
        T = ctx["options"]["temp"]
        r = TA.quasi_harmonic(aligned, masses, T)
        nu = np.sort(r["frequencies_cm"])
        F, dof = len(aligned), 3 * len(sub) - 6
        warn = "" if F > 3 * dof else ("  WARNING: {} frames for {} degrees of freedom - at most {} modes are "
                                       "determined and the entropy is not converged; use more frames or fewer atoms."
                                       .format(F, dof, F - 1))
        order = np.argsort(r["frequencies_cm"])
        return {
            "title": "Quasi-harmonic",
            "views": [("Entropy per mode", lambda p: p.set_lines([{"x": r["frequencies_cm"][order], "y": r["per_mode"][order], "style": "points"}],
                                                                "frequency (cm-1)", "S per mode (J/mol/K)", "Quasi-harmonic modes", ymin=0)),
                      ("Frequencies", lambda p: p.set_lines([{"x": list(range(1, len(nu) + 1)), "y": nu, "style": "points"}],
                                                           "mode", "frequency (cm-1)", "Quasi-harmonic frequencies", ymin=0))],
            "table": (["mode", "frequency (cm-1)", "S (J/mol/K)"],
                      [[k + 1, "{:.2f}".format(r["frequencies_cm"][i]), "{:.3f}".format(r["per_mode"][i])] for k, i in enumerate(order)]),
            "summary": ("S(quasi-harmonic) = {:.1f} J/(mol K), T*S = {:.1f} kJ/mol at {:.0f} K; Schlitter bound {:.1f} J/(mol K). "
                        "{} atoms, {} frames, {} modes.{}").format(r["entropy_qh"], r["entropy_qh"] * T / 1000, T,
                                                                   r["entropy_schlitter"], len(sub), F, r["modes"], warn),
            "data": ("frequency_cm\tS_J_mol_K", np.column_stack([r["frequencies_cm"][order], r["per_mode"][order]])),
        }

    def _compute_rotamers(self, ctx):
        order, members = self._residues_of(ctx, [i for i in ctx["atoms"] if ctx["is_protein"][i]])
        defs = []
        for r in order:
            rn = str(ctx["resname"][members[r][0]]).upper()
            names = {ctx["names"][i]: i for i in members[r]}
            for k, quad in enumerate(TA.CHI_ATOMS.get(rn, [])):
                if all(a in names for a in quad):
                    defs.append((r, k + 1, [names[a] for a in quad]))
        if not defs:
            raise ValueError("No residue of the atom set has side-chain dihedrals.")
        need = sorted({i for _, _, q in defs for i in q})
        pos = {g: k for k, g in enumerate(need)}
        X = self._coords(ctx, need)
        x, xl = self._xs(ctx)
        rows, angles = [], {}
        for r, k, quad in defs:
            ang = TA.dihedral(X, *[pos[i] for i in quad])
            pops, trans = TA.rotamer_summary(ang)
            mean = math.degrees(math.atan2(np.sin(np.radians(ang)).mean(), np.cos(np.radians(ang)).mean()))
            name = "{} chi{}".format(ctx["res_labels"][r], k)
            angles[name] = ang
            rows.append([ctx["res_labels"][r], "chi{}".format(k), "{:.0f}".format(mean), "{:.0f}".format(pops["g+"] * 100),
                         "{:.0f}".format(pops["t"] * 100), "{:.0f}".format(pops["g-"] * 100), trans, (name,)])
        first = rows[0][-1][0]
        self._rot_angles = (angles, x, xl)
        return {
            "title": "Rotamers",
            "views": [("Selected dihedral vs time", lambda p: self._plot_dihedral(p, first)),
                      ("Selected dihedral histogram", lambda p: self._plot_dihedral_hist(p, first))],
            "frames_x": x,
            "table": (["residue", "angle", "circular mean", "g+ %", "t %", "g- %", "transitions"], rows),
            "summary": "{} side-chain dihedrals in {} residues. Double-click a row to plot it.".format(len(defs), len(order)),
            "series": {"{} (deg)".format(n): (x, list(a), xl, "angle (deg)") for n, a in list(angles.items())[:40]},
            "data": (xl + "\t" + "\t".join(angles), np.column_stack([x] + [angles[n] for n in angles])),
            "row_action": "dihedral",
        }

    def _plot_dihedral(self, p, name):
        angles, x, xl = self._rot_angles
        p.set_lines([{"x": x, "y": angles[name], "style": "points", "label": name}], xl, "angle (deg)", name)

    def _plot_dihedral_hist(self, p, name):
        angles, _x, _xl = self._rot_angles
        h, edges = np.histogram(angles[name], bins=36, range=(-180, 180))
        centres = 0.5 * (edges[1:] + edges[:-1])
        p.set_lines([{"x": centres, "y": h / max(h.sum(), 1) * 100, "style": "bars"}], "angle (deg)", "% of frames", name, ymin=0)

    def _compute_msd(self, ctx):
        n = len(ctx["names"])
        if ctx["options"]["which"] == 0:
            atoms = [i for i in range(n) if ctx["is_solvent"][i] and ctx["elements"][i] == "O"]
            if not atoms:
                raise ValueError("There are no water oxygens in this object.")
        else:
            atoms = ctx["atoms"]
        X = self._coords(ctx, atoms)
        if len(X) < 4:
            raise ValueError("Needs at least four frames.")
        lags, msd = TA.mean_square_displacement(X, len(X) // 2, ctx["boxes"])
        dt = ctx["dt"]
        step = ctx["frames_idx"][1] - ctx["frames_idx"][0] if len(ctx["frames_idx"]) > 1 else 1
        t = lags * step * (dt if dt > 0 else 1.0)
        D, slope, icpt, i0, i1 = TA.diffusion_coefficient(t, msd, ctx["options"]["fit0"] / 100, ctx["options"]["fit1"] / 100)
        tl = "lag time (ps)" if dt > 0 else "lag (frames)"
        fit = [slope * tt + icpt for tt in t[i0:i1 + 1]]
        unitD = "1e-5 cm^2/s" if dt > 0 else "A^2/frame /6"
        return {
            "title": "MSD",
            "views": [("MSD", lambda p: p.set_lines([{"x": t, "y": msd, "label": "MSD"},
                                                     {"x": t[i0:i1 + 1], "y": fit, "label": "fit", "dash": True, "color": PALETTE[3]}],
                                                    tl, "MSD (A^2)", "Mean square displacement", ymin=0))],
            "table": ([tl, "MSD (A^2)"], [["{:.3f}".format(a), "{:.4f}".format(b)] for a, b in zip(t, msd)]),
            "summary": ("D = {:.3f} {} from the slope between {} and {} ({} atoms, boxes {}). {}").format(
                D if dt > 0 else slope / 6.0, unitD, _short(t[i0]), _short(t[i1]), len(atoms),
                "used" if ctx["boxes"] else "not available", "" if dt > 0 else "Set the time per frame to get D in cm^2/s.") +
                ("  Reference: TIP3P water ~5.8, experiment 2.3 (1e-5 cm^2/s, 298 K)." if ctx["options"]["which"] == 0 else ""),
            "series": {"MSD (A^2)": (list(t), list(msd), tl, "MSD (A^2)")},
            "data": (tl + "\tMSD", np.column_stack([t, msd])),
        }

    def _compute_stats(self, ctx):
        s = ctx.get("series")
        if s is None:
            raise ValueError("Choose a series (run another analysis first or load a file).")
        x, y, xl, yl = s
        y = np.asarray(y, dtype=float)
        x = np.asarray(x, dtype=float)
        k = int(len(y) * ctx["options"]["discard"] / 100)
        x, y = x[k:], y[k:]
        if len(y) < 8:
            raise ValueError("The series is too short (fewer than 8 points).")
        d = TA.describe(y)
        levels, _sem = TA.block_average(y)
        acf, tau = TA.autocorrelation(y)
        h, edges = np.histogram(y, bins=min(40, max(8, len(y) // 10)))
        centres = 0.5 * (edges[1:] + edges[:-1])
        run = np.cumsum(y) / np.arange(1, len(y) + 1)
        return {
            "title": "Statistics",
            "views": [("Series and running mean", lambda p: p.set_lines([{"x": x, "y": y, "label": ctx["series_name"]},
                                                                        {"x": x, "y": run, "label": "running mean", "color": PALETTE[3]}],
                                                                       xl, yl, ctx["series_name"])),
                      ("Histogram", lambda p: p.set_lines([{"x": centres, "y": h / h.sum() * 100, "style": "bars"}], yl, "% of points",
                                                          "Distribution", ymin=0)),
                      ("Autocorrelation", lambda p: p.set_lines([{"x": list(range(len(acf))), "y": acf}], "lag (points)", "ACF",
                                                                "Autocorrelation (tau = {:.1f})".format(tau), hline=0.0)),
                      ("Block averaging", lambda p: p.set_lines([{"x": list(range(len(levels))), "y": [l[1] for l in levels], "label": "SEM"}],
                                                                "blocking level (block size 2^k)", "SEM", "Flyvbjerg-Petersen blocking", ymin=0))],
            "table": (["quantity", "value"], [["points", d["n"]], ["mean", "{:.6g}".format(d["mean"])], ["standard deviation", "{:.6g}".format(d["std"])],
                                              ["SEM (block averaging)", "{:.4g}".format(d["sem_block"])],
                                              ["SEM (naive, uncorrelated)", "{:.4g}".format(d["std"] / math.sqrt(d["n"]))],
                                              ["integrated correlation time (points)", "{:.2f}".format(d["tau_frames"])],
                                              ["independent samples", "{:.0f}".format(d["n_independent"])],
                                              ["min / max", "{:.6g} / {:.6g}".format(d["min"], d["max"])]]),
            "summary": "{}: {:.6g} +/- {:.3g} (block-averaged SEM), {:.0f} independent samples of {}.".format(
                ctx["series_name"], d["mean"], d["sem_block"], d["n_independent"], d["n"]),
            "data": (xl + "\t" + yl, np.column_stack([x, y])),
        }

    # =====================================================================
    #  Result display
    # =====================================================================
    def _set_result_widgets(self, result):
        self.result = result
        self.combo_view.remove_all()
        for child in self.actions_box.get_children():
            self.actions_box.remove(child)
        has = result is not None
        self.btn_export.set_sensitive(has)
        self.btn_save.set_sensitive(has)
        if not has:
            self.plot.clear()
            self.table.set_model(None)
            self.summary.get_buffer().set_text("")
            return
        for name, _setter in result["views"]:
            self.combo_view.append_text(name)
        self.combo_view.set_active(0)
        cols, rows = result.get("table", ([], []))
        for c in self.table.get_columns():
            self.table.remove_column(c)
        store = Gtk.ListStore(*([str] * len(cols)))
        self._row_payload = []
        for r in rows:
            payload = r[len(cols)] if len(r) > len(cols) else None
            self._row_payload.append(payload)
            store.append([str(v) for v in r[:len(cols)]])
        for k, c in enumerate(cols):
            col = Gtk.TreeViewColumn(c, Gtk.CellRendererText(), text=k)
            col.set_sort_column_id(k)
            col.set_resizable(True)
            self.table.append_column(col)
        self.table.set_model(store)
        self.summary.get_buffer().set_text(result.get("summary", ""))
        for label, callback in result.get("actions", []):
            b = Gtk.Button(label=label)
            b.connect("clicked", lambda _b, cb=callback: self._safe(cb))
            self.actions_box.pack_start(b, False, False, 0)
        self.actions_box.show_all()

    def _on_view_changed(self, combo):
        k = combo.get_active()
        if self.result is None or k < 0:
            return
        self.result["views"][k][1](self.plot)

    def _on_plot_click(self, info):
        r = self.result
        if r is None:
            return
        frame = None
        if self.plot.kind == "lines" and "frames_x" in r and info.get("index") is not None:
            xs = r["frames_x"]
            k = min(range(len(xs)), key=lambda i: abs(xs[i] - info["x"]))
            if abs(xs[k] - info["x"]) <= (abs(xs[1] - xs[0]) if len(xs) > 1 else 1):
                frame = r["ctx"]["frames_idx"][k]
        elif self.plot.kind == "categories" and "map_frames" in r:
            frame = r["map_frames"][info["col"]]
            self._select_residues(r["ctx"], [self._backbone_keys(r)[info["row"]]])
        elif self.plot.kind == "scatter" and "scatter_frames" in r:
            frame = r["scatter_frames"][info["index"]]
        elif self.plot.kind == "heatmap" and "cell_action" in r:
            keys = r["cell_action"]
            if info["row"] < len(keys) and info["col"] < len(keys):
                self._select_residues(r["ctx"], [keys[info["row"]], keys[info["col"]]])
        if frame is not None:
            self.vm_session.set_frame(frame)

    def _backbone_keys(self, r):
        if "_bb_keys" not in r:
            r["_bb_keys"] = self._backbone(r["ctx"])[1]
        return r["_bb_keys"]

    def _on_table_row(self, _view, path, _column):
        r = self.result
        if r is None:
            return
        payload = self._row_payload[path.get_indices()[0]]
        if payload is None:
            return
        kind = r.get("row_action")
        ctx = r["ctx"]
        if kind == "hbond":
            d, h, a, frame = payload
            self._select_atoms(ctx, [d, h, a])
            self.vm_session.set_frame(frame)
        elif kind == "residues":
            self._select_residues(ctx, list(payload))
        elif kind == "frame":
            self.vm_session.set_frame(payload[0])
        elif kind == "path":
            atoms, frame = payload
            self._select_atoms(ctx, atoms)
            self.vm_session.set_frame(frame)
        elif kind == "dihedral":
            name = payload[0]
            k = self.combo_view.get_active()
            (self._plot_dihedral if k == 0 else self._plot_dihedral_hist)(self.plot, name)
            self.result["views"] = [("Selected dihedral vs time", lambda p: self._plot_dihedral(p, name)),
                                    ("Selected dihedral histogram", lambda p: self._plot_dihedral_hist(p, name))]

    # ------------------------------------------------------------- 3D actions
    def _select_atoms(self, ctx, indices):
        v = ctx["vobject"]
        self.vm_session._selection_function_set(None)
        sel = self.vm_session.selections[self.vm_session.current_selection]
        mode = sel.selection_mode
        sel.selection_mode = "atom"
        self.vm_session._selection_function_set({v.atoms[i] for i in indices})
        sel.selection_mode = mode
        self.vm_session.vm_glcore.queue_draw()

    def _select_residues(self, ctx, residue_keys):
        atoms = [i for i in range(len(ctx["names"])) if int(ctx["reskey"][i]) in set(residue_keys)]
        self._select_atoms(ctx, atoms)

    def _color_by(self, ctx, atoms, values):
        v = ctx["vobject"]
        values = np.asarray(values, dtype=float)
        lo, hi = float(values.min()), float(values.max())
        hi = hi if hi > lo else lo + 1
        for i, val in zip(atoms, values):
            v.atoms[i].color = list(colormap((val - lo) / (hi - lo)))
        self.vm_session.set_color_by_index(vobject=v, indexes=[], color=None)
        self.vm_session.vm_glcore.queue_draw()

    def _new_object(self, ctx, frames_full, name, coords=None):
        """New object of the same system with the given frames (indices of the
        source object, or explicit coordinates (F, N, 3))."""
        system = ctx["system"]
        v = ctx["vobject"]
        if system is None or len(system.atoms) != len(v.atoms):
            raise ValueError("The object has no matching pDynamo system.")
        X = np.asarray(v.frames[frames_full], dtype=np.float32) if coords is None else np.asarray(coords, dtype=np.float32)
        c3 = system.coordinates3
        saved = [[c3[i, 0], c3[i, 1], c3[i, 2]] for i in range(len(system.atoms))]
        try:
            for i in range(len(system.atoms)):
                system.coordinates3[i, 0], system.coordinates3[i, 1], system.coordinates3[i, 2] = (float(c) for c in X[0, i])
            new = self.p_session._add_vismol_object_to_easyhybrid_session(system=system, name=name)
            new.frames = X
        finally:
            for i in range(len(system.atoms)):
                for k in range(3):
                    system.coordinates3[i, k] = saved[i][k]
        self.vm_session.vm_glcore.queue_draw()
        return new

    def _animate_pc(self, ctx, sub, p, c, masses):
        v = ctx["vobject"]
        mode = TA.mode_frames(p["mean"], p["eigenvectors"][:, c], p["eigenvalues"][c], 21, 2.0, masses)
        disp = mode - p["mean"][None]
        # . every other atom follows the nearest analysed atom of its residue (or the nearest overall)
        base = np.asarray(v.frames[ctx["frames_idx"][0]], dtype=float)
        sub_pos = {g: k for k, g in enumerate(sub)}
        owner = np.zeros(len(v.atoms), dtype=int)
        Xs = base[sub]
        by_res = {}
        for g, k in sub_pos.items():
            by_res.setdefault(int(ctx["reskey"][g]), []).append(k)
        for i in range(len(v.atoms)):
            if i in sub_pos:
                owner[i] = sub_pos[i]
                continue
            cand = by_res.get(int(ctx["reskey"][i]))
            if cand:
                owner[i] = cand[0]
            else:
                owner[i] = int(np.argmin(((Xs - base[i]) ** 2).sum(axis=1)))
        # . displacements are in the frame of the average; rotate them onto the source frame
        R = TA._kabsch_rotation(p["mean"] - p["mean"].mean(0), Xs - Xs.mean(0))
        frames = np.array([base + (disp[f] @ R)[owner] for f in range(len(disp))])
        self._new_object(ctx, None, "PC{}_{}".format(c + 1, v.name), coords=frames)

    def _average_object(self, ctx, sub, fl):
        v = ctx["vobject"]
        X = np.asarray(v.frames[ctx["frames_idx"][0]], dtype=float).copy()
        X[sub] = fl["average"]
        self._new_object(ctx, None, "average_{}".format(v.name), coords=X[None])

    def _export_bfactor_pdb(self, ctx, sub, fl):
        dialog = Gtk.FileChooserDialog(title="Export PDB with B-factors", parent=self.window, action=Gtk.FileChooserAction.SAVE)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        dialog.set_do_overwrite_confirmation(True)
        dialog.set_current_name("{}_bfactors.pdb".format(ctx["vobject"].name))
        if dialog.run() == Gtk.ResponseType.OK:
            write_bfactor_pdb(dialog.get_filename(), ctx, sub, fl)
        dialog.destroy()

    # ------------------------------------------------------------- output
    def on_export_image(self, _b):
        dialog = Gtk.FileChooserDialog(title="Export plot", parent=self.window, action=Gtk.FileChooserAction.SAVE)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        dialog.set_do_overwrite_confirmation(True)
        dialog.set_current_name("{}.png".format((self.combo_view.get_active_text() or "plot").replace(" ", "_")))
        if dialog.run() == Gtk.ResponseType.OK:
            self.plot.export(dialog.get_filename())
        dialog.destroy()

    def on_save_data(self, _b):
        if self.result is None:
            return
        dialog = Gtk.FileChooserDialog(title="Save data", parent=self.window, action=Gtk.FileChooserAction.SAVE)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        dialog.set_do_overwrite_confirmation(True)
        dialog.set_current_name("{}.dat".format(self.result["title"].replace(" ", "_").lower()))
        if dialog.run() == Gtk.ResponseType.OK:
            self.save_data(dialog.get_filename())
        dialog.destroy()

    def save_data(self, path):
        header, data = self.result["data"]
        with open(path, "w") as handle:
            handle.write("# EasyHybrid Trajectory Analysis: {}\n# {}\n".format(
                self.result["title"], self.result.get("summary", "").replace("\n", " ")))
            handle.write("# " + header + "\n")
            if isinstance(data, np.ndarray):
                for row in data:
                    handle.write("\t".join("{:.6g}".format(float(v)) for v in row) + "\n")
            else:
                for row in data:
                    handle.write("\t".join(str(v) for v in row) + "\n")

    # ------------------------------------------------------------- statistics input
    def _refresh_series_combo(self):
        current = self.combo_series.get_active_text()
        self.combo_series.remove_all()
        names = list(self.series_registry)
        for name in names:
            self.combo_series.append_text(name)
        if names:
            self.combo_series.set_active(names.index(current) if current in names else len(names) - 1)

    def _on_load_series(self, _b):
        dialog = Gtk.FileChooserDialog(title="Load a series (columns x y)", parent=self.window, action=Gtk.FileChooserAction.OPEN)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
        if dialog.run() == Gtk.ResponseType.OK:
            path = dialog.get_filename()
            try:
                data = np.loadtxt(path, comments=("#", "@"))
                if data.ndim == 1:
                    x, ys = np.arange(len(data)), [data]
                else:
                    x, ys = data[:, 0], [data[:, k] for k in range(1, data.shape[1])]
                for k, y in enumerate(ys):
                    name = "{} [col {}]".format(os.path.basename(path), k + 2)
                    self.series_registry[name] = (list(x), list(y), "x", "y")
                self._refresh_series_combo()
            except Exception as error:
                self._message("Could not read {}:\n{}".format(path, error))
        dialog.destroy()

    # ------------------------------------------------------------- misc
    def _safe(self, callback):
        try:
            callback()
        except Exception as error:
            traceback.print_exc()
            self._message(str(error))

    def _message(self, text):
        dialog = Gtk.MessageDialog(transient_for=self.window, modal=True, message_type=Gtk.MessageType.INFO,
                                   buttons=Gtk.ButtonsType.OK, text=text)
        dialog.run()
        dialog.destroy()

    def _on_destroy(self, *_):
        self._cancel = True
        self.Visible = False
        self.window = None


# =========================================================================
def _short(v):
    return "{:g}".format(round(float(v), 4))


def _apply_fits(X, fit_X, reference):
    """Superposes every frame of X with the transform that fits fit_X on reference."""
    out = np.empty_like(X)
    rc = reference.mean(axis=0)
    for f in range(len(X)):
        P = fit_X[f]
        pc = P.mean(axis=0)
        R = TA._kabsch_rotation(P - pc, reference - rc)
        out[f] = (X[f] - pc) @ R + rc
    return out


def write_bfactor_pdb(path, ctx, sub, fl):
    """PDB of the average structure of the analysed atoms with B-factors and ANISOU records."""
    with open(path, "w") as handle:
        handle.write("REMARK   EasyHybrid Trajectory Analysis: average structure, B = 8 pi^2/3 <dr^2>\n")
        for k, i in enumerate(sub):
            ch = (ctx["chain"][i] or "A")[:1]
            rkey = int(ctx["reskey"][i])
            resi = ctx["res_objs"][rkey]
            rnum = getattr(resi, "index", 0)
            name = ctx["names"][i]
            name = (" " + name) if len(name) < 4 else name
            x, y, z = fl["average"][k]
            el = ctx["elements"][i]
            handle.write("ATOM  {:5d} {:<4s} {:>3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n".format(
                (k + 1) % 100000, name[:4], str(ctx["resname"][i])[:3], ch, int(rnum) % 10000, x, y, z, 1.0,
                min(fl["b_factors"][k], 999.99), el[:2]))
            U = fl["U"][k] * 1e4
            handle.write("ANISOU{:5d} {:<4s} {:>3s} {:1s}{:4d}  {:7d}{:7d}{:7d}{:7d}{:7d}{:7d}      {:>2s}\n".format(
                (k + 1) % 100000, name[:4], str(ctx["resname"][i])[:3], ch, int(rnum) % 10000,
                int(round(U[0, 0])), int(round(U[1, 1])), int(round(U[2, 2])), int(round(U[0, 1])),
                int(round(U[0, 2])), int(round(U[1, 2])), el[:2]))
        handle.write("END\n")
