#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QC/MM and molecular simulations using pDynamo3
#  Module: Umbrella Sampling Browser
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      Analysis > Umbrella Sampling Browser... (2026-10-23) - the windows of an
#      umbrella sampling run (the <name>_umb_sam folder, pdynamo/us_data.py):
#
#        Windows        map of the windows (a strip for 1D, a matrix for 2D),
#                       coloured by frames / deviation from the target / drift /
#                       spread / overlap / PMF; click to select windows, load
#                       them (one object each, or one trajectory), representative
#                       frames, export
#        Diagnostics    RC distributions and the overlap of neighbours, mean RC
#                       against the target, RC along the run, 2D coverage
#        Free energy    PMF with the windows (pDynamo's WHAM files or a WHAM
#                       computed here), convergence (first half / second half /
#                       blocks), lowest-barrier path on a 2D PMF -> windows ->
#                       trajectory
#
#      Same layout as the Trajectory Analysis window; plots in util/analysis_plot.py.
#

import math
import os
import threading
import traceback

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, GLib

import numpy as np

from gui.widgets.custom_widgets import SystemComboBox
from util.analysis_plot import AnalysisPlot, PALETTE, palette_names

BLACK, GREY, BLUE, RED, ORANGE, GREEN = (0.1, 0.1, 0.1), (0.55, 0.55, 0.55), (0.08, 0.40, 0.75), (0.80, 0.16, 0.16), \
                                        (0.90, 0.45, 0.05), (0.18, 0.55, 0.20)
KCAL = 4.184
LOW_OVERLAP = 0.03          # . neighbours below this share almost no samples (WHAM cannot connect them)

VIEWS = [("Windows", "map", "Window map"),
         ("Diagnostics", "hist", "RC distributions and overlap"),
         ("Diagnostics", "overlap", "Overlap between neighbours"),
         ("Diagnostics", "drift", "Mean RC against the target"),
         ("Diagnostics", "series", "RC along the run (selected windows)"),
         ("Diagnostics", "coverage", "Coverage of the RC space (2D)"),
         ("Free energy", "pmf", "PMF and windows"),
         ("Free energy", "conv", "Convergence (WHAM by halves and blocks)"),
         ("Free energy", "mep", "Lowest-barrier path (2D)")]

MAP_METRICS = [("frames", "Frames (coordinates)"), ("deviation", "|mean − target| / restraint width"),
               ("drift", "Drift (2nd half − 1st half)"), ("spread", "RC spread (std)"),
               ("overlap", "Lowest overlap with a neighbour"), ("pmf", "PMF at the target")]


def _fmt_key(key):
    return "w" + "_".join(str(k) for k in key)


class UmbrellaBrowserWindow:

    def __init__(self, main=None):
        self.main = main
        self.home = main.home
        self.p_session = main.p_session
        self.vm_session = main.vm_session
        self.Visible = False
        self.window = None
        self.run = None
        self.stats = None
        self.selected = []          # . ordered window keys
        self.result = None
        self.key = None
        self._row_payload = []
        self._own = {}              # . own WHAM results: "all", "first", "second", "blocks"
        self._pmf = None
        self._mep = {"start": None, "end": None, "path": [], "click": None}
        self._busy = False

    # =====================================================================
    #  Window
    # =====================================================================
    def open_window(self, folder=None):
        if not self.Visible:
            self._build()
        else:
            self.window.present()
        if folder:
            self.entry_folder.set_text(folder)
            self._open(folder)

    def _build(self):
        self.window = Gtk.Window(title="Umbrella Sampling Browser")
        self.window.set_default_size(1280, 820)
        self.window.set_transient_for(self.main.window)
        self.window.connect("destroy", self._on_destroy)
        outer = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6, margin=8)
        self.window.add(outer)
        outer.pack_start(self._build_input_frame(), False, False, 0)

        paned = Gtk.Paned(orientation=Gtk.Orientation.HORIZONTAL)
        outer.pack_start(paned, True, True, 0)
        left = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
        left.set_size_request(320, -1)
        self.store = Gtk.TreeStore(str, str)
        parents = {}
        for cat, key, title in VIEWS:
            if cat not in parents:
                parents[cat] = self.store.append(None, ["<b>{}</b>".format(cat), ""])
            self.store.append(parents[cat], [title, key])
        self.tree = Gtk.TreeView(model=self.store)
        self.tree.set_headers_visible(False)
        self.tree.append_column(Gtk.TreeViewColumn("View", Gtk.CellRendererText(), markup=0))
        self.tree.expand_all()
        self.tree.get_selection().connect("changed", self._on_tree_changed)
        sw = Gtk.ScrolledWindow()
        sw.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)
        sw.add(self.tree)
        sw.set_size_request(-1, 260)
        left.pack_start(sw, False, False, 0)

        frame = Gtk.Frame(label=" Selected windows ")
        box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=4, margin=6)
        frame.add(box)
        self.label_sel = Gtk.Label(xalign=0)
        self.label_sel.set_line_wrap(True)
        self.label_sel.set_max_width_chars(42)
        self.label_sel.set_selectable(True)
        box.pack_start(self.label_sel, False, False, 0)
        hb = Gtk.Box(spacing=4)
        for text, cb, tip in (("Clear", self._sel_clear, "Nothing selected"),
                              ("All", self._sel_all, "Every window, in order"),
                              ("Row", lambda: self._sel_line(0), "2D: the windows with the same RC1 index as the last selected one"),
                              ("Column", lambda: self._sel_line(1), "2D: the windows with the same RC2 index as the last selected one"),
                              ("Range", self._sel_range, "1D: every window between the first and the last selected")):
            b = Gtk.Button(label=text)
            b.set_tooltip_text(tip)
            b.connect("clicked", lambda _b, f=cb: self._safe(f))
            hb.pack_start(b, True, True, 0)
            setattr(self, "btn_sel_" + text.lower(), b)
        box.pack_start(hb, False, False, 0)
        left.pack_start(frame, False, False, 0)

        frame = Gtk.Frame(label=" Load / export the selected windows ")
        box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=4, margin=6)
        frame.add(box)
        for text, cb, tip in (("One object per window", lambda: self._load(mode="each"),
                               "Each selected window as its own object (all its frames of the chosen source and stride)"),
                              ("One trajectory (windows in order)", lambda: self._load(mode="concat"),
                               "The frames of the selected windows concatenated, in the order of the selection"),
                              ("Representative frames (one per window)", lambda: self._load(mode="rep"),
                               "Per window the frame whose RC is closest to the window target: a 'movie' along the selection"),
                              ("Export...", self._export,
                               "Representative frames or all frames of the selection to DCD, XYZ or PDB")):
            b = Gtk.Button(label=text)
            b.set_tooltip_text(tip)
            b.connect("clicked", lambda _b, f=cb: self._safe(f))
            box.pack_start(b, False, False, 0)
        left.pack_start(frame, False, False, 0)
        self.progress = Gtk.ProgressBar(show_text=True)
        self.progress.set_text("")
        left.pack_start(self.progress, False, False, 0)
        self.label_hint = Gtk.Label(xalign=0)
        self.label_hint.set_line_wrap(True)
        self.label_hint.set_max_width_chars(44)
        left.pack_start(self.label_hint, False, False, 0)
        paned.pack1(left, False, False)

        right = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6, margin_start=6)
        bar = Gtk.Box(spacing=6)
        bar.pack_start(Gtk.Label(label="View"), False, False, 0)
        self.combo_view = Gtk.ComboBoxText()
        self.combo_view.connect("changed", self._on_view_changed)
        bar.pack_start(self.combo_view, False, False, 0)
        self.label_metric = Gtk.Label(label="Color by")
        bar.pack_start(self.label_metric, False, False, 0)
        self.combo_metric = Gtk.ComboBoxText()
        for key, text in MAP_METRICS:
            self.combo_metric.append(key, text)
        self.combo_metric.set_active(1)
        self.combo_metric.connect("changed", lambda *_: self._refresh())
        bar.pack_start(self.combo_metric, False, False, 0)
        self.btn_export = Gtk.Button(label="Export Image...")
        self.btn_export.connect("clicked", self.on_export_image)
        self.btn_save = Gtk.Button(label="Save Table...")
        self.btn_save.connect("clicked", self.on_save_table)
        bar.pack_end(self.btn_save, False, False, 0)
        bar.pack_end(self.btn_export, False, False, 0)
        # . colour map of the maps (2026-10-24): the palettes of the PES Analysis window
        self.check_reverse = Gtk.CheckButton(label="Reverse")
        self.check_reverse.set_tooltip_text("Reverse the colour map")
        self.combo_palette = Gtk.ComboBoxText()
        for key, text in palette_names():
            self.combo_palette.append(key, text)
        self.combo_palette.set_active(0)
        self.combo_palette.set_tooltip_text("Colour map of the maps and colour bars. Difference maps (convergence) keep a "
                                            "diverging map centred on zero unless a diverging one is chosen (bam, berlin, "
                                            "vik, coolwarm, coolwarm2)")
        apply_palette = lambda *_: self.plot.set_palette(self.combo_palette.get_active_id(), self.check_reverse.get_active())
        self.combo_palette.connect("changed", apply_palette)
        self.check_reverse.connect("toggled", apply_palette)
        bar.pack_end(self.check_reverse, False, False, 0)
        bar.pack_end(self.combo_palette, False, False, 0)
        bar.pack_end(Gtk.Label(label="Colors"), False, False, 0)
        right.pack_start(bar, False, False, 0)
        vpaned = Gtk.Paned(orientation=Gtk.Orientation.VERTICAL)
        self.plot = AnalysisPlot()
        self.plot.on_click = self._on_plot_click
        vpaned.pack1(self.plot, True, False)
        self.notebook = Gtk.Notebook()
        self.table = Gtk.TreeView()
        self.table.connect("row-activated", self._on_table_row)
        sw = Gtk.ScrolledWindow()
        sw.add(self.table)
        self.notebook.append_page(sw, Gtk.Label(label="Table"))
        self.text = Gtk.TextView(editable=False, monospace=True, wrap_mode=Gtk.WrapMode.WORD_CHAR)
        sw = Gtk.ScrolledWindow()
        sw.add(self.text)
        self.notebook.append_page(sw, Gtk.Label(label="Details"))
        self.notebook.set_size_request(-1, 200)
        vpaned.pack2(self.notebook, False, False)
        vpaned.set_position(470)
        right.pack_start(vpaned, True, True, 0)
        self.actions_box = Gtk.Box(spacing=6)
        right.pack_start(self.actions_box, False, False, 0)
        paned.pack2(right, True, False)

        self.window.show_all()
        self.Visible = True
        if self.p_session.psystem:
            self.combobox_systems.set_active_system(e_id=self.p_session.active_id)
        self._show(None)
        self._update_selection_label()

    def _build_input_frame(self):
        frame = Gtk.Frame(label=" Umbrella sampling run ")
        grid = Gtk.Grid(column_spacing=8, row_spacing=6, margin=8)
        frame.add(grid)
        self.entry_folder = Gtk.Entry(hexpand=True)
        self.entry_folder.set_placeholder_text("the <name>_umb_sam folder written by Simulation > Umbrella Sampling")
        self.entry_folder.connect("activate", lambda e: self._safe(lambda: self._open(e.get_text())))
        b = Gtk.Button(label="Browse...")
        b.connect("clicked", self.on_browse)
        grid.attach(Gtk.Label(label="Folder", xalign=0), 0, 0, 1, 1)
        grid.attach(self.entry_folder, 1, 0, 5, 1)
        grid.attach(b, 6, 0, 1, 1)
        self.combobox_systems = SystemComboBox(self.main)
        self.combobox_systems.set_tooltip_text("The system whose coordinates the windows hold (same atoms): loaded frames become new objects of it")
        grid.attach(Gtk.Label(label="System", xalign=0), 0, 1, 1, 1)
        grid.attach(self.combobox_systems, 1, 1, 1, 1)
        self.combo_source = Gtk.ComboBoxText()
        from pdynamo.us_data import SOURCES
        for code, label, _sub in SOURCES:
            self.combo_source.append(code, label)
        self.combo_source.set_active(0)
        self.combo_source.connect("changed", lambda *_: (self._update_selection_label(), self._refresh()))
        grid.attach(Gtk.Label(label="Coordinates", xalign=0), 2, 1, 1, 1)
        grid.attach(self.combo_source, 3, 1, 1, 1)
        self.spin_stride = Gtk.SpinButton.new_with_range(1, 10000, 1)
        self.spin_stride.set_tooltip_text("Every n-th frame of each window")
        grid.attach(Gtk.Label(label="Stride", xalign=0), 4, 1, 1, 1)
        grid.attach(self.spin_stride, 5, 1, 1, 1)
        self.combo_pmf = Gtk.ComboBoxText()
        self.combo_pmf.set_tooltip_text("PMF used by the views: a file written by pDynamo's WHAM (wham/*_pmf.log) or a WHAM computed here from the restraint data")
        self.combo_pmf.connect("changed", lambda *_: self._safe(self._on_pmf_changed))
        grid.attach(Gtk.Label(label="PMF", xalign=0), 0, 2, 1, 1)
        grid.attach(self.combo_pmf, 1, 2, 3, 1)
        self.combo_unit = Gtk.ComboBoxText()
        self.combo_unit.append("kJ", "kJ/mol")
        self.combo_unit.append("kcal", "kcal/mol")
        self.combo_unit.set_active(0)
        self.combo_unit.connect("changed", lambda *_: self._refresh())
        grid.attach(Gtk.Label(label="Energy", xalign=0), 4, 2, 1, 1)
        grid.attach(self.combo_unit, 5, 2, 1, 1)
        self.label_info = Gtk.Label(xalign=0)
        self.label_info.set_line_wrap(True)
        self.label_info.set_selectable(True)
        grid.attach(self.label_info, 0, 3, 7, 1)
        return frame

    # =====================================================================
    #  Opening a run
    # =====================================================================
    def on_browse(self, _b):
        chooser = Gtk.FileChooserDialog(title="Umbrella sampling folder (<name>_umb_sam)", parent=self.window,
                                        action=Gtk.FileChooserAction.SELECT_FOLDER)
        chooser.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
        system = self.p_session.psystem.get(self.p_session.active_id)
        base = getattr(system, "e_working_folder", None) if system is not None else None
        if self.run is not None:
            chooser.set_current_folder(os.path.dirname(self.run.folder))
        elif base and os.path.isdir(base):
            chooser.set_current_folder(base)
        response, path = chooser.run(), chooser.get_filename()
        chooser.destroy()
        if response == Gtk.ResponseType.OK and path:
            self.entry_folder.set_text(path)
            self._safe(lambda: self._open(path))

    def _open(self, folder):
        from pdynamo import us_data as U
        folder = folder.strip()
        run = U.USRun(folder)
        self.run = run
        self.stats = run.stats()
        self.selected = []
        self._own = {}
        self._mep = {"start": None, "end": None, "path": [], "click": None}
        self._overlap_cache = None
        self.combo_pmf.remove_all()
        for path in run.pmf_files():
            self.combo_pmf.append(path, "pDynamo WHAM: " + os.path.relpath(path, run.folder))
        self.combo_pmf.append("own", "WHAM computed here (from the restraint data)")
        self.combo_pmf.set_active(0)
        frames = [len(run.frame_files(k, "dc")) for k in run.keys]
        n = [self.stats[k]["n"] for k in run.keys]
        rcinfo = []
        for d in range(run.dims):
            p = run.params.get("RC{}".format(d + 1)) or {}
            atoms = " / ".join("{}({})–{}({}){}".format(a[0], a[1], a[2], a[3], "" if str(a[4]) == "1" else " (×{})".format(a[4]))
                               for a in p.get("RC", []))
            try:
                step = float(p.get("dincre"))
            except (TypeError, ValueError):
                step = None
            # . 2026-10-26: a run started from a trajectory has no step (its targets are
            #   the RC values of the frames): was "from 2.059 by 0.0"
            if step:
                span = "from {:.3f} by {}".format(float(p.get("dminimum", float("nan"))), p.get("dincre"))
            else:
                span = "targets taken from the starting trajectory"
            rcinfo.append("RC{}: {} windows, {} (k = {}){}".format(
                d + 1, run.shape[d], span, p.get("force_constant", "?"), ": " + atoms if atoms else ""))
        missing = int(np.prod(run.shape)) - len(run.keys)
        self.label_info.set_text("{}  ·  {}D, {} windows{}  ·  {:.0f} K  ·  {}–{} restraint samples and {}–{} frames per window\n{}".format(
            os.path.basename(run.folder), run.dims, len(run.keys), " ({} missing)".format(missing) if missing > 0 else "",
            run.temperature, min(n), max(n), min(frames) if frames else 0, max(frames) if frames else 0, "\n".join(rcinfo)))
        self.btn_sel_row.set_sensitive(run.dims == 2)
        self.btn_sel_column.set_sensitive(run.dims == 2)
        self.btn_sel_range.set_sensitive(run.dims == 1)
        self._update_selection_label()
        it = self.store.iter_children(self.store.get_iter_first())
        self.tree.get_selection().unselect_all()
        self.tree.get_selection().select_iter(it)
        self._refresh()

    def _on_pmf_changed(self):
        from pdynamo import us_data as U
        self._pmf = None
        if self.run is None:
            return
        code = self.combo_pmf.get_active_id()
        if code and code != "own":
            dims, x, F = U.read_pmf(code)
            if dims == self.run.dims:
                self._pmf = {"x": x, "F": F, "source": "pDynamo WHAM ({})".format(os.path.basename(code))}
        elif code == "own":
            if "all" in self._own:
                x, F = self._own["all"]
                self._pmf = {"x": x, "F": F, "source": "WHAM computed here"}
            else:
                self._compute_own(["all"], then=self._on_pmf_changed)
                return
        self._mep["path"] = []
        self._refresh()

    # =====================================================================
    #  Selection
    # =====================================================================
    def _toggle(self, key):
        if key in self.selected:
            self.selected.remove(key)
        else:
            self.selected.append(key)
        self._update_selection_label()
        self._refresh(keep_view=True)

    def _sel_clear(self):
        self.selected = []
        self._update_selection_label(); self._refresh(keep_view=True)

    def _sel_all(self):
        self.selected = list(self.run.keys)
        self._update_selection_label(); self._refresh(keep_view=True)

    def _sel_line(self, axis):
        if not self.selected:
            raise ValueError("Select one window first (click it in the map).")
        last = self.selected[-1]
        line = sorted(k for k in self.run.keys if k[axis] == last[axis])
        self.selected = [k for k in self.selected if k not in line] + line
        self._update_selection_label(); self._refresh(keep_view=True)

    def _sel_range(self):
        if len(self.selected) < 2:
            raise ValueError("Select the first and the last window of the range.")
        a, b = self.selected[0][0], self.selected[-1][0]
        step = 1 if b >= a else -1
        self.selected = [(i,) for i in range(a, b + step, step) if (i,) in self.run.windows]
        self._update_selection_label(); self._refresh(keep_view=True)

    def _update_selection_label(self):
        if self.run is None:
            self.label_sel.set_text("Open a run first.")
            return
        if not self.selected:
            self.label_sel.set_text("None. Click windows in the map (or in the PMF) to select them; the order of the clicks is the order of a trajectory.")
            return
        source = self.combo_source.get_active_id() or "dc"
        frames = sum(len(self.run.frame_files(k, source)) for k in self.selected)
        names = ", ".join(_fmt_key(k) for k in self.selected[:12]) + (" ..." if len(self.selected) > 12 else "")
        self.label_sel.set_text("{} window(s), {} frames ({}): {}".format(len(self.selected), frames,
                                                                         self.combo_source.get_active_text(), names))

    # =====================================================================
    #  Views
    # =====================================================================
    def _on_tree_changed(self, selection):
        model, it = selection.get_selected()
        if it is None or not model[it][1]:
            return
        self.key = model[it][1]
        self._refresh()

    def _refresh(self, keep_view=False):
        if self.window is None:
            return
        if self.run is None or self.key is None:
            self._show(None)
            return
        k = self.combo_view.get_active() if keep_view else 0
        try:
            result = getattr(self, "_view_" + self.key)()
        except Exception as error:
            traceback.print_exc()
            result = {"views": [], "summary": "Could not show this view: {}".format(error)}
        self._show(result, view=k)
        self.label_metric.set_visible(self.key == "map")
        self.combo_metric.set_visible(self.key == "map")

    def _unit(self):
        return (1.0 / KCAL, "kcal/mol") if self.combo_unit.get_active_id() == "kcal" else (1.0, "kJ/mol")

    def _target_str(self, key):
        t = self.run.target(key)
        return ", ".join("{:.3f}".format(v) for v in t) if t else "?"

    def _overlaps(self):
        if getattr(self, "_overlap_cache", None) is None:
            if self.run.dims == 1:
                pairs, data = self.run.overlap_1d()
            else:
                pairs, data = self.run.overlap_2d(), None
            low = {}
            for a, b, o in pairs:
                low[a] = min(low.get(a, 1.0), o); low[b] = min(low.get(b, 1.0), o)
            self._overlap_cache = (pairs, data, low)
        return self._overlap_cache

    def _pmf_at(self, key):
        """ PMF at the window target (nearest grid point), in kJ/mol. """
        if self._pmf is None:
            return float("nan")
        t = self.run.target(key)
        if t is None:
            return float("nan")
        if self.run.dims == 1:
            x = np.asarray(self._pmf["x"]); F = np.asarray(self._pmf["F"])
            return float(F[int(np.argmin(np.abs(x - t[0])))])
        xs, ys = self._pmf["x"]
        return float(self._pmf["F"][int(np.argmin(np.abs(np.asarray(xs) - t[0]))), int(np.argmin(np.abs(np.asarray(ys) - t[1])))])

    def _metric(self, key, metric):
        st = self.stats.get(key) or {}
        if not st.get("n"):
            return float("nan")
        if metric == "frames":
            return float(len(self.run.frame_files(key, self.combo_source.get_active_id() or "dc")))
        if metric == "deviation":
            return st["deviation"]
        if metric == "drift":
            return float(np.linalg.norm(st["drift"]))
        if metric == "spread":
            return float(np.linalg.norm(st["std"]))
        if metric == "overlap":
            return self._overlaps()[2].get(key, float("nan"))
        if metric == "pmf":
            return self._pmf_at(key) * self._unit()[0]
        return float("nan")

    def _window_table(self):
        low = self._overlaps()[2]
        rows = []
        for key in self.run.keys:
            st = self.stats.get(key) or {}
            mean = st.get("mean"); std = st.get("std"); drift = st.get("drift")
            rows.append([_fmt_key(key), self._target_str(key),
                         ", ".join("{:.3f}".format(v) for v in mean) if mean is not None else "",
                         ", ".join("{:.3f}".format(v) for v in std) if std is not None else "",
                         ", ".join("{:+.3f}".format(v) for v in drift) if drift is not None else "",
                         "{:.2f}".format(st.get("deviation", float("nan"))), "{:.3f}".format(low.get(key, float("nan"))),
                         st.get("n", 0), len(self.run.frame_files(key, "dc")), len(self.run.frame_files(key, "eq")),
                         "yes" if key in self.selected else "", ("window", key)])
        return (["window", "target", "mean RC", "std", "drift", "deviation", "min overlap", "samples", "frames (dc)",
                 "frames (eq)", "selected"], rows)

    # ------------------------------------------------------------- map
    def _view_map(self):
        run = self.run
        metric = self.combo_metric.get_active_id() or "deviation"
        label = dict(MAP_METRICS)[metric] + (" ({})".format(self._unit()[1]) if metric == "pmf" else "")
        if run.dims == 1:
            order = sorted(run.keys, key=lambda k: (run.target(k) or (k[0],))[0])
            cols = [self._target_str(k) for k in order]
            matrix = [[self._metric(k, metric) for k in order]]
            pos = {k: (0, n) for n, k in enumerate(order)}
            xl, yl = "window target ({})".format(run.labels[0]), ""
            rl = ["windows"]
        else:
            # . rows and columns in increasing target values, as the PMF views
            ni, nj = run.shape
            def target_of(axis, n):
                vals = [run.target(k)[axis] for k in run.keys if k[axis] == n and run.target(k)]
                return sum(vals) / len(vals) if vals else float(n)
            ti = {i: target_of(0, i) for i in range(ni)}; tj = {j: target_of(1, j) for j in range(nj)}
            oi = sorted(range(ni), key=lambda i: ti[i]); oj = sorted(range(nj), key=lambda j: tj[j])
            matrix = [[self._metric((i, j), metric) if (i, j) in run.windows else float("nan") for j in oj] for i in oi]
            pos = {(i, j): (r, c) for r, i in enumerate(oi) for c, j in enumerate(oj)}
            rl = ["{:.2f}".format(ti[i]) for i in oi]
            cols = ["{:.2f}".format(tj[j]) for j in oj]
            xl, yl = "RC2 target", "RC1 target"
        self._map_pos = {v: k for k, v in pos.items()}
        overlays = []
        if self.selected:
            cells = [pos[k] for k in self.selected if k in pos]
            overlays.append({"cells": cells, "style": "box", "color": RED, "width": 2.5})
            if len(cells) > 1 and run.dims == 2:
                overlays.append({"cells": cells, "style": "path", "color": RED, "width": 1.5})
        title = "{} windows coloured by {}".format(len(run.keys), label.lower())
        views = [("Window map", lambda p: p.set_heatmap(matrix, xl, yl, title, row_labels=rl, col_labels=cols, value_label=label,
                                                       overlays=overlays))]
        n_low = sum(1 for v in self._overlaps()[2].values() if v < LOW_OVERLAP)
        bad = [k for k in run.keys if (self.stats.get(k) or {}).get("deviation", 0) > 3.0]
        text = ("Click a window to select it (again to unselect); the order of the clicks is the order used by 'One trajectory' and "
                "'Representative frames'. Double-click a row of the table to load that window.\n\n"
                "deviation = |mean RC − target| in units of the restraint width sqrt(kT/2k): above ~3 the window did not stay near its "
                "target (the restraint is too weak there or the system escaped).\n"
                "Windows with a neighbour overlap below {:.2f}: {}. Windows with deviation > 3: {}{}.").format(
                    LOW_OVERLAP, n_low, len(bad), " (" + ", ".join(_fmt_key(k) for k in bad[:10]) + (" ..." if len(bad) > 10 else "") + ")" if bad else "")
        return {"views": views, "table": self._window_table(), "row_action": "window", "summary": text, "map": True,
                "hint": "Click windows to select them; double-click a table row to load one."}

    # ------------------------------------------------------------- diagnostics
    def _view_hist(self):
        run = self.run
        pairs, data, low = self._overlaps()
        if run.dims == 2:
            return self._view_coverage(note="2D runs: the distributions of all windows are shown as the coverage of the RC plane.")
        edges = data["edges"]; centers = 0.5 * (edges[1:] + edges[:-1])
        series = []
        for n, k in enumerate(data["order"]):
            h = data["hist"][k]
            dens = h / max(1e-12, (edges[1] - edges[0]))
            sel = k in self.selected
            series.append({"x": list(centers), "y": list(dens), "label": _fmt_key(k) if sel else None,
                           "color": RED if sel else PALETTE[n % len(PALETTE)], "width": 2.4 if sel else 1.2})
        rows = [[_fmt_key(a), _fmt_key(b), self._target_str(a), self._target_str(b), "{:.3f}".format(o),
                 "LOW" if o < LOW_OVERLAP else "", ("pair", (a, b))] for a, b, o in pairs]
        worst = sorted(pairs, key=lambda t: t[2])[:5]
        text = ("Distribution of the reaction coordinate in every window (normalized). Neighbouring windows must overlap for WHAM to "
                "connect them; overlap = Σ min(p, q) over a common grid (1 = identical, 0 = disjoint).\n"
                "Lowest: " + "; ".join("{}–{} {:.3f}".format(_fmt_key(a), _fmt_key(b), o) for a, b, o in worst) +
                "\n\nGaps (overlap < {:.2f}) usually need an extra window in between or a softer force constant.".format(LOW_OVERLAP))
        return {"views": [("Distributions of the RC per window", lambda p: p.set_lines(series, run.labels[0], "probability density",
                                                                                         "RC histograms of the {} windows (selected in red)".format(len(run.keys)), ymin=0))],
                "table": (["window", "next", "target", "next target", "overlap", ""], rows), "row_action": "pair", "summary": text,
                "hint": "Double-click a pair to select both windows."}

    def _view_overlap(self):
        run = self.run
        pairs, _data, low = self._overlaps()
        rows = [[_fmt_key(a), _fmt_key(b), "{:.3f}".format(o), "LOW" if o < LOW_OVERLAP else "", ("pair", (a, b))] for a, b, o in pairs]
        text = "Overlap coefficient of neighbouring windows. Below {:.2f}: {} of {} pairs.".format(
            LOW_OVERLAP, sum(1 for p in pairs if p[2] < LOW_OVERLAP), len(pairs))
        if run.dims == 1:
            labels = ["{}–{}".format(_fmt_key(a), _fmt_key(b)) for a, b, _ in pairs]
            x = list(range(len(pairs)))
            colors = [RED if o < LOW_OVERLAP else BLUE for _a, _b, o in pairs]
            views = [("Overlap per neighbour pair", lambda p: p.set_lines([{"x": x, "y": [o for _a, _b, o in pairs], "style": "bars",
                                                                            "color": BLUE, "colors": colors}], "pair of windows", "overlap",
                                                                          "Overlap of neighbouring windows (red: below {:.2f})".format(LOW_OVERLAP),
                                                                          ymin=0, hline=LOW_OVERLAP, xlabels=labels))]
        else:
            ni, nj = run.shape
            m_i = [[float("nan")] * nj for _ in range(ni)]
            m_j = [[float("nan")] * nj for _ in range(ni)]
            for (a, b, o) in pairs:
                (m_i if b[0] != a[0] else m_j)[a[0]][a[1]] = o
            lowest = [[low.get((i, j), float("nan")) for j in range(nj)] for i in range(ni)]
            views = [("Lowest overlap of each window", lambda p: p.set_heatmap(lowest, "window j (RC2)", "window i (RC1)",
                                                                               "Lowest overlap with a neighbour", value_label="overlap",
                                                                               vmin=0.0)),
                     ("Overlap with the next window along RC1", lambda p: p.set_heatmap(m_i, "window j", "window i",
                                                                                        "Overlap of (i, j) and (i+1, j)", value_label="overlap", vmin=0.0)),
                     ("Overlap with the next window along RC2", lambda p: p.set_heatmap(m_j, "window j", "window i",
                                                                                        "Overlap of (i, j) and (i, j+1)", value_label="overlap", vmin=0.0))]
        return {"views": views, "table": (["window", "neighbour", "overlap", ""], rows), "row_action": "pair", "summary": text,
                "hint": "Double-click a pair to select both windows."}

    def _view_drift(self):
        run = self.run
        keys = [k for k in run.keys if (self.stats.get(k) or {}).get("n")]
        rows = self._window_table()
        if run.dims == 1:
            keys = sorted(keys, key=lambda k: run.target(k)[0])
            t = [run.target(k)[0] for k in keys]
            mean = [float(self.stats[k]["mean"][0]) for k in keys]
            std = [float(self.stats[k]["std"][0]) for k in keys]
            views = [("Mean RC − target", lambda p: p.set_lines(
                        [{"x": t, "y": [m - a for m, a in zip(mean, t)], "label": "mean − target", "color": BLUE, "width": 2},
                         {"x": t, "y": [m - a for m, a in zip(mean, t)], "style": "points", "color": BLUE},
                         {"x": t, "y": std, "label": "std", "color": GREY, "dash": True}],
                        "window target ({})".format(run.labels[0]), "RC − target", "Where each window actually sampled", hline=0.0)),
                     ("Drift (2nd half − 1st half)", lambda p: p.set_lines(
                        [{"x": t, "y": [float(self.stats[k]["drift"][0]) for k in keys], "style": "bars", "color": ORANGE}],
                        "window target", "drift of the mean", "Change of the mean RC between the halves of the run", hline=0.0))]
        else:
            series = []
            for k in keys:
                t = run.target(k); m = self.stats[k]["mean"]
                bad = self.stats[k]["deviation"] > 3.0
                series.append({"x": [t[1], float(m[1])], "y": [t[0], float(m[0])], "color": RED if bad else BLUE, "width": 1.6})
            series.append({"x": [run.target(k)[1] for k in keys], "y": [run.target(k)[0] for k in keys], "style": "points", "color": BLACK,
                           "label": "targets"})
            views = [("Targets → means", lambda p: p.set_lines(series, "RC2", "RC1",
                                                               "Each window: target (dot) to the mean RC sampled (red: deviation > 3)"))]
        text = ("Mean of the RC in each window against its target, and the drift of the mean between the first and the second half of "
                "the data collection (a large drift: the window was not equilibrated or moved to another basin).")
        return {"views": views, "table": rows, "row_action": "window", "summary": text,
                "hint": "Double-click a row to load that window."}

    def _view_series(self):
        run = self.run
        keys = self.selected[:8] or [run.keys[len(run.keys) // 2]]
        views = []
        for d in range(run.dims):
            series = []
            for n, k in enumerate(keys):
                x = run.rc(k)[:, d]
                step = max(1, len(x) // 2000)
                series.append({"x": list(range(0, len(x), step)), "y": list(x[::step]), "label": _fmt_key(k), "color": PALETTE[n % len(PALETTE)]})
                t = run.target(k)
                if t:
                    series.append({"x": [0, len(x) - 1], "y": [t[d], t[d]], "color": PALETTE[n % len(PALETTE)], "dash": True})
            views.append(("{} along the run".format(run.labels[d]), lambda p, s=series, d=d: p.set_lines(
                s, "step (×{} MD steps)".format(run.freq_res), run.labels[d], "{} in the selected windows (dashed: targets)".format(run.labels[d]))))
        rows = [[_fmt_key(k), self._target_str(k)] + ["{:.3f}".format(v) for v in self.stats[k]["mean"]] + [("window", k)] for k in keys]
        return {"views": views, "table": (["window", "target"] + ["mean " + l for l in run.labels], rows), "row_action": "window",
                "summary": "The RC in every step of the data collection of the selected windows (the first 8; the middle window when none "
                           "is selected).", "hint": "Select windows in the map to compare them here."}

    def _view_coverage(self, note=None):
        run = self.run
        if run.dims != 2:
            return {"views": [], "summary": "Coverage is for 2D runs; see 'RC distributions and overlap'."}
        h, ex, ey = run.coverage_2d()
        dx, dy = ex[1] - ex[0], ey[1] - ey[0]
        dots = [((run.target(k)[0] - ex[0]) / dx - 0.5, (run.target(k)[1] - ey[0]) / dy - 0.5) for k in run.keys if run.target(k)]
        sel = [((run.target(k)[0] - ex[0]) / dx - 0.5, (run.target(k)[1] - ey[0]) / dy - 0.5) for k in self.selected if run.target(k)]
        matrix = [[math.log10(v) if v > 0 else float("nan") for v in row] for row in h]
        overlays = [{"cells": dots, "style": "dots", "color": BLACK, "radius": 2.0}]
        if sel:
            overlays.append({"cells": sel, "style": "dots", "color": RED, "radius": 4.0})
        empty = int((h == 0).sum())
        text = (note + "\n\n" if note else "") + ("Samples of all windows on the RC1 × RC2 plane (log10 of the count per bin; grey: no "
                "sample), with the window targets (black dots, selected in red). Grey regions between targets are holes WHAM cannot "
                "fill. {} of {} bins empty.".format(empty, h.size))
        rl = ["{:.2f}".format(v) for v in 0.5 * (ex[1:] + ex[:-1])]
        cl = ["{:.2f}".format(v) for v in 0.5 * (ey[1:] + ey[:-1])]
        return {"views": [("Coverage", lambda p: p.set_heatmap(matrix, "RC2", "RC1", "Sampling of the RC plane", row_labels=rl, col_labels=cl,
                                                               value_label="log10 samples", overlays=overlays))],
                "summary": text, "hint": "Black dots: window targets."}

    # ------------------------------------------------------------- free energy
    def _view_pmf(self):
        run = self.run
        f, unit = self._unit()
        if self._pmf is None:
            return {"views": [], "summary": "No PMF: choose a pDynamo WHAM file or 'WHAM computed here' in the PMF box above.",
                    "hint": "Choose a PMF above."}
        P = self._pmf
        if run.dims == 1:
            x = list(map(float, P["x"])); F = [float(v) * f for v in P["F"]]
            wx = [run.target(k)[0] for k in run.keys if run.target(k)]
            wy = [self._pmf_at(k) * f for k in run.keys if run.target(k)]
            sx = [run.target(k)[0] for k in self.selected]; sy = [self._pmf_at(k) * f for k in self.selected]
            series = [{"x": x, "y": F, "label": "PMF", "color": BLUE, "width": 2.2},
                      {"x": wx, "y": wy, "style": "points", "color": BLACK, "label": "window targets"}]
            if sx:
                series.append({"x": sx, "y": sy, "style": "points", "color": RED, "label": "selected"})
            self._pmf_click = ("1d", wx)
            views = [("PMF with the windows", lambda p: p.set_lines(series, run.labels[0], "PMF ({})".format(unit), P["source"], ymin=0))]
            fin = [v for v in F if v == v]
            text = "{}: barrier (max − min over the sampled range) {:.2f} {}.".format(P["source"], max(fin) - min(fin), unit) if fin else ""
        else:
            xs, ys = P["x"]; F = P["F"] * f
            dx = (xs[-1] - xs[0]) / max(1, len(xs) - 1); dy = (ys[-1] - ys[0]) / max(1, len(ys) - 1)
            cell = lambda t: ((t[0] - xs[0]) / dx, (t[1] - ys[0]) / dy)
            dots = [cell(run.target(k)) for k in run.keys if run.target(k)]
            sel = [cell(run.target(k)) for k in self.selected if run.target(k)]
            overlays = [{"cells": dots, "style": "dots", "color": BLACK, "radius": 2.0}]
            if sel:
                overlays.append({"cells": sel, "style": "dots", "color": RED, "radius": 4.5})
                if len(sel) > 1:
                    overlays.append({"cells": sel, "style": "path", "color": RED, "width": 1.5})
            self._pmf_click = ("2d", xs, ys)
            matrix = [list(row) for row in F]
            views = [("PMF with the windows", lambda p: p.set_heatmap(matrix, "RC2", "RC1", P["source"], row_labels=["{:.2f}".format(v) for v in xs],
                                                                     col_labels=["{:.2f}".format(v) for v in ys], value_label="PMF ({})".format(unit),
                                                                     overlays=overlays))]
            text = "{}. Black dots: window targets; red: the selection (in order). Click near a window to select it.".format(P["source"])
        # . 2026-10-28: without overlap the PMF between the windows is not determined by the data
        pairs = self._overlaps()[0]
        low = [o for _a, _b, o in pairs if o < LOW_OVERLAP]
        if low:
            text = ("WARNING: {} of {} pairs of neighbouring windows do not overlap (below {:.2f}, lowest {:.3f}): "
                    "across those gaps the PMF is not determined by the data, and different WHAM programs or bin "
                    "choices give different curves. ").format(len(low), len(pairs), LOW_OVERLAP, min(low)) + text
        rows = [[_fmt_key(k), self._target_str(k), "{:.2f}".format(self._pmf_at(k) * f), ("window", k)] for k in run.keys]
        return {"views": views, "table": (["window", "target", "PMF at target ({})".format(unit)], rows), "row_action": "window",
                "summary": text, "pmf_click": True, "hint": "Click near a window on the PMF to select it."}

    def _view_conv(self):
        run = self.run
        f, unit = self._unit()
        needed = ["all", "first", "second"] + (["b0", "b1", "b2", "b3"] if run.dims == 1 else [])
        missing = [k for k in needed if k not in self._own]
        if missing:
            return {"views": [], "actions": [("Compute (WHAM on halves{})".format(" and quarters" if run.dims == 1 else ""),
                                              lambda: self._compute_own(missing, then=self._refresh))],
                    "summary": "PMF from the whole data collection, from its first half and from its second half{} (WHAM computed here "
                               "from the restraint data, so the three use the same bins and equations). Their difference shows whether "
                               "the run is converged. 1D: a few seconds; 2D: about 20 s per PMF.".format(", and from each quarter" if run.dims == 1 else ""),
                    "hint": "Press Compute."}
        if run.dims == 1:
            x, Fa = self._own["all"]; _, F1 = self._own["first"]; _, F2 = self._own["second"]
            series = [{"x": list(x), "y": list(Fa * f), "label": "all", "color": BLACK, "width": 2.4},
                      {"x": list(x), "y": list(F1 * f), "label": "1st half", "color": BLUE},
                      {"x": list(x), "y": list(F2 * f), "label": "2nd half", "color": ORANGE}]
            blocks = [{"x": list(x), "y": list(self._own["b{}".format(q)][1] * f), "label": "quarter {}".format(q + 1),
                       "color": PALETTE[q]} for q in range(4)]
            diff = np.nanmax(np.abs(F2 - F1)) * f
            ok = np.isfinite(Fa)
            bar = lambda F: (np.nanmax(F[ok]) - np.nanmin(F[ok])) * f
            rows = [["all", "{:.2f}".format(bar(Fa))], ["1st half", "{:.2f}".format(bar(F1))], ["2nd half", "{:.2f}".format(bar(F2))]] + \
                   [["quarter {}".format(q + 1), "{:.2f}".format(bar(self._own["b{}".format(q)][1]))] for q in range(4)]
            views = [("Halves", lambda p: p.set_lines(series, run.labels[0], "PMF ({})".format(unit),
                                                       "PMF by halves of the data collection · largest difference {:.2f} {}".format(diff, unit), ymin=0)),
                     ("Quarters", lambda p: p.set_lines(blocks, run.labels[0], "PMF ({})".format(unit), "PMF by quarters", ymin=0))]
            text = "Largest difference between the PMFs of the two halves: {:.2f} {}.".format(diff, unit)
            return {"views": views, "table": (["data", "PMF range ({})".format(unit)], rows), "summary": text}
        (cx, cy), Fa = self._own["all"]; _, F1 = self._own["first"]; _, F2 = self._own["second"]
        D = (F2 - F1) * f
        rl = ["{:.2f}".format(v) for v in cx]; cl = ["{:.2f}".format(v) for v in cy]
        lim = np.nanmax(np.abs(D)) if np.isfinite(D).any() else 1.0
        low = np.isfinite(Fa) & (Fa * f < np.nanmin(Fa * f) + 30.0 * (1 if unit == "kJ/mol" else 1 / KCAL))
        views = [("2nd half − 1st half", lambda p: p.set_heatmap([list(r) for r in D], "RC2", "RC1", "PMF difference between the halves",
                                                                  row_labels=rl, col_labels=cl, value_label=unit, vmin=-lim, vmax=lim, diverging=True)),
                 ("PMF, all data", lambda p: p.set_heatmap([list(r) for r in Fa * f], "RC2", "RC1", "PMF (WHAM computed here)",
                                                           row_labels=rl, col_labels=cl, value_label=unit)),
                 ("PMF, 1st half", lambda p: p.set_heatmap([list(r) for r in F1 * f], "RC2", "RC1", "PMF, first half", row_labels=rl,
                                                           col_labels=cl, value_label=unit)),
                 ("PMF, 2nd half", lambda p: p.set_heatmap([list(r) for r in F2 * f], "RC2", "RC1", "PMF, second half", row_labels=rl,
                                                           col_labels=cl, value_label=unit))]
        text = "Mean |difference| between the halves where the PMF is within 30 kJ/mol of the minimum: {:.2f} {}; largest: {:.2f} {}.".format(
            float(np.nanmean(np.abs(D[low]))) if low.any() else float("nan"), unit, lim, unit)
        return {"views": views, "summary": text}

    def _view_mep(self):
        run = self.run
        f, unit = self._unit()
        if run.dims != 2:
            return {"views": [], "summary": "The lowest-barrier path is for 2D runs. In 1D, select a range of windows in the map."}
        if self._pmf is None:
            return {"views": [], "summary": "Choose a PMF above first."}
        from pdynamo import us_data as U
        xs, ys = self._pmf["x"]; F = self._pmf["F"]
        if self._mep["start"] is None or self._mep["end"] is None:
            self._mep["start"], self._mep["end"] = U.default_endpoints(F)
        if not self._mep["path"]:
            self._mep["path"] = U.lowest_barrier_path(F, self._mep["start"], self._mep["end"])
        path = self._mep["path"]
        prof = [float(F[c]) * f for c in path]
        arc = [0.0]
        for a, b in zip(path, path[1:]):
            arc.append(arc[-1] + math.hypot(xs[b[0]] - xs[a[0]], ys[b[1]] - ys[a[1]]))
        windows = self._path_windows(path)
        overlays = [{"cells": [(float(i), float(j)) for i, j in path], "style": "path", "color": RED, "width": 2.5},
                    {"cells": [(float(self._mep["start"][0]), float(self._mep["start"][1]))], "style": "dots", "color": GREEN, "radius": 5},
                    {"cells": [(float(self._mep["end"][0]), float(self._mep["end"][1]))], "style": "dots", "color": ORANGE, "radius": 5}]
        matrix = [list(r) for r in F * f]
        self._pmf_click = ("mep", xs, ys)
        top = int(np.argmax(prof)) if prof else 0
        views = [("Path on the PMF", lambda p: p.set_heatmap(matrix, "RC2", "RC1", "Lowest-barrier path (green: start, orange: end)",
                                                             row_labels=["{:.2f}".format(v) for v in xs], col_labels=["{:.2f}".format(v) for v in ys],
                                                             value_label="PMF ({})".format(unit), overlays=overlays)),
                 ("Profile along the path", lambda p: p.set_lines([{"x": arc, "y": prof, "color": BLUE, "width": 2.2},
                                                                   {"x": arc, "y": prof, "style": "points", "color": BLUE}],
                                                                  "path length (RC units)", "PMF ({})".format(unit),
                                                                  "Free energy along the path · barrier {:.2f} {}".format(
                                                                      (max(prof) - prof[0]) if prof else 0, unit)))]
        rows = [[n + 1, _fmt_key(k), self._target_str(k), "{:.2f}".format(self._pmf_at(k) * f), ("window", k)] for n, k in enumerate(windows)]
        text = ("Path from the start to the end cell that keeps the highest PMF value met as low as possible (the route through the "
                "lowest saddle), then the lowest sum. Saddle at RC1 {:.3f}, RC2 {:.3f}: {:.2f} {} above the start.\n"
                "{} windows along it (nearest window to each point). 'Use these windows' selects them in order; then 'Representative "
                "frames' gives one structure per window along the reaction.\n"
                "To move the ends: press 'Set start' or 'Set end' and click the PMF.").format(
                    xs[path[top][0]] if path else float("nan"), ys[path[top][1]] if path else float("nan"),
                    (max(prof) - prof[0]) if prof else 0, unit, len(windows))
        actions = [("Use these windows", lambda: self._use_path(windows)),
                   ("Set start", lambda: self._set_mep_click("start")), ("Set end", lambda: self._set_mep_click("end"))]
        return {"views": views, "table": (["#", "window", "target", "PMF ({})".format(unit)], rows), "row_action": "window",
                "summary": text, "actions": actions, "pmf_click": True, "hint": "Red line: path; table: its windows."}

    def _path_windows(self, path):
        xs, ys = self._pmf["x"]
        keys = [k for k in self.run.keys if self.run.target(k)]
        T = np.array([self.run.target(k) for k in keys])
        scale = np.array([max(1e-9, np.ptp(T[:, 0])), max(1e-9, np.ptp(T[:, 1]))])
        out = []
        for i, j in path:
            d = np.linalg.norm((T - np.array([xs[i], ys[j]])) / scale, axis=1)
            k = keys[int(np.argmin(d))]
            if not out or out[-1] != k:
                if k not in out:
                    out.append(k)
        return out

    def _use_path(self, windows):
        self.selected = list(windows)
        self._update_selection_label()
        self._message("{} windows of the path selected, in order.".format(len(windows)))
        self._refresh(keep_view=True)

    def _set_mep_click(self, which):
        self._mep["click"] = which
        self._message("Click the PMF to set the {} of the path.".format(which))

    # =====================================================================
    #  Own WHAM in the background
    # =====================================================================
    def _compute_own(self, parts, then=None):
        if self._busy:
            self._message("A calculation is already running."); return
        run = self.run
        spans = {"all": None, "first": (0.0, 0.5), "second": (0.5, 1.0), "b0": (0.0, 0.25), "b1": (0.25, 0.5),
                 "b2": (0.5, 0.75), "b3": (0.75, 1.0)}
        self._busy = True
        self.progress.set_fraction(0.0)
        self.progress.set_text("WHAM 0/{}".format(len(parts)))

        def work():
            out = {}
            try:
                for n, p in enumerate(parts):
                    out[p] = run.wham(part=spans[p])
                    GLib.idle_add(lambda n=n: (self.progress.set_fraction((n + 1) / len(parts)),
                                               self.progress.set_text("WHAM {}/{}".format(n + 1, len(parts)))) and False)
                GLib.idle_add(done, out, None)
            except Exception as error:
                traceback.print_exc()
                GLib.idle_add(done, out, str(error))

        def done(out, error):
            self._busy = False
            self.progress.set_text("")
            self.progress.set_fraction(0.0)
            if run is not self.run:
                return False
            self._own.update(out)
            if error:
                self._message("WHAM failed: " + error)
            elif then is not None:
                self._safe(then)
            return False
        self._message("Computing WHAM ({}) in the background...".format(", ".join(parts)))
        threading.Thread(target=work, daemon=True).start()

    # =====================================================================
    #  Loading frames
    # =====================================================================
    def _target_system(self):
        e_id = self.combobox_systems.get_system_id()
        system = self.p_session.psystem.get(e_id) if e_id is not None else None
        if system is None:
            raise ValueError("Choose the system of the run (its atoms) in the System box.")
        return system

    def _check_atoms(self, system, X):
        if len(X) != len(system.atoms):
            raise ValueError("The frames have {} atoms but '{}' has {}: choose the system the umbrella sampling was run on.".format(
                len(X), system.label, len(system.atoms)))

    def _load(self, mode):
        if self.run is None:
            raise ValueError("Open a run first.")
        if not self.selected:
            raise ValueError("Select windows first (click them in the map).")
        if self._busy:
            raise ValueError("Wait for the running task.")
        system = self._target_system()
        run = self.run
        source = self.combo_source.get_active_id() or "dc"
        stride = int(self.spin_stride.get_value())
        keys = list(self.selected)
        first = next((f for k in keys for f in run.frame_files(k, source)[:1]), None)
        if first is None:
            raise ValueError("The selected windows have no frames of this source ({}).".format(self.combo_source.get_active_text()))
        self._check_atoms(system, run.read_frame(first))          # . before reading everything
        self._busy = True

        def work():
            try:
                blocks = []
                for n, k in enumerate(keys):
                    if mode == "rep":
                        m, _d = run.representative(k, source)
                        frames, idx = run.frames(k, source, only=[m] if m is not None else [])
                    else:
                        frames, idx = run.frames(k, source, stride)
                    blocks.append((k, frames, idx))
                    GLib.idle_add(lambda n=n: (self.progress.set_fraction((n + 1) / len(keys)),
                                               self.progress.set_text("reading {}/{}".format(n + 1, len(keys)))) and False)
                GLib.idle_add(done, blocks, None)
            except Exception as error:
                traceback.print_exc()
                GLib.idle_add(done, None, str(error))

        def done(blocks, error):
            self._busy = False
            self.progress.set_text(""); self.progress.set_fraction(0.0)
            if error:
                self._message("Loading failed: " + error); return False
            try:
                tag = {"dc": "dc", "eq": "eq", "opt": "start"}[source]
                if mode == "each":
                    made = 0
                    for k, frames, _idx in blocks:
                        if frames:
                            self._new_object(system, frames, "US {} {}".format(_fmt_key(k), tag)); made += 1
                    self._message("{} object(s) added, one per window.".format(made))
                else:
                    frames = [X for _k, fr, _i in blocks for X in fr]
                    if not frames:
                        raise ValueError("No frames in the selected windows for this source.")
                    name = "US {} {}–{} ({} windows)".format("rep" if mode == "rep" else tag, _fmt_key(keys[0]), _fmt_key(keys[-1]), len(keys))
                    self._new_object(system, frames, name)
                    self._loaded = [(k, i) for k, _fr, idx in blocks for i in idx]
                    self._message("'{}': {} frames. The mini player plays it; frame n is window {}.".format(
                        name, len(frames), "n of the selection" if mode == "rep" else "in the selection order"))
            except Exception as error:
                traceback.print_exc()
                self._message("Loading failed: {}".format(error))
            return False
        self._message("Reading {} window(s)...".format(len(keys)))
        threading.Thread(target=work, daemon=True).start()

    def _new_object(self, system, frames, name):
        self._check_atoms(system, frames[0])
        X = np.asarray(frames, dtype=np.float32)
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
        # . the object was created (treeview row, QC dynamic bonds) with one frame: both are rebuilt
        #   now that it has all its frames (dynamic_bonds holds one entry per frame)
        self.p_session._apply_QC_representation_to_vobject(vismol_object=new)
        self.main.main_treeview.refresh_number_of_frames()
        self.main.main_treeview.refresh_trajectory_scalebar()
        self.vm_session.vm_glcore.queue_draw()
        return new

    def _export(self):
        if self.run is None or not self.selected:
            raise ValueError("Select windows first.")
        system = self._target_system()
        dialog = Gtk.Dialog(title="Export the selected windows", transient_for=self.window, modal=True)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, "Export...", Gtk.ResponseType.OK)
        box = dialog.get_content_area(); box.set_spacing(6); box.set_margin_start(10); box.set_margin_end(10)
        rep = Gtk.RadioButton.new_with_label(None, "Representative frames (one per window, in the selection order)")
        allf = Gtk.RadioButton.new_with_label_from_widget(rep, "All frames (with the stride above)")
        box.pack_start(rep, False, False, 0); box.pack_start(allf, False, False, 0)
        box.pack_start(Gtk.Label(label="Format by extension: .dcd (trajectory), .xyz (multi-frame), .pdb (one file per frame)", xalign=0), False, False, 0)
        dialog.show_all()
        response = dialog.run(); only_rep = rep.get_active(); dialog.destroy()
        if response != Gtk.ResponseType.OK:
            return
        path = self._ask_save("Export frames", ["dcd", "xyz", "pdb"], "umbrella_windows.dcd")
        if not path:
            return
        from pdynamo import us_data as U
        source = self.combo_source.get_active_id() or "dc"
        frames = []
        for k in self.selected:
            if only_rep:
                m, _d = self.run.representative(k, source)
                fr, _ = self.run.frames(k, source, only=[m] if m is not None else [])
            else:
                fr, _ = self.run.frames(k, source, int(self.spin_stride.get_value()))
            frames += fr
        if not frames:
            raise ValueError("No frames to export.")
        self._check_atoms(system, frames[0])
        out = U.export_frames(path, system, frames)
        self._message("{} frame(s) written: {}{}".format(len(frames), out[0], " ... ({} files)".format(len(out)) if len(out) > 1 else ""))

    # =====================================================================
    #  Display and events
    # =====================================================================
    def _show(self, result, view=0):
        self.result = result
        self.combo_view.handler_block_by_func(self._on_view_changed)
        self.combo_view.remove_all()
        self.combo_view.handler_unblock_by_func(self._on_view_changed)
        for child in self.actions_box.get_children():
            self.actions_box.remove(child)
        for c in self.table.get_columns():
            self.table.remove_column(c)
        self.table.set_model(None)
        self._row_payload = []
        if result is None:
            self.plot.clear()
            self.text.get_buffer().set_text("Open the folder of an umbrella sampling run (<name>_umb_sam) above." if self.run is None else "")
            self.label_hint.set_text("")
            self.btn_export.set_sensitive(False); self.btn_save.set_sensitive(False)
            return
        views = result.get("views") or []
        self.plot.set_visible(bool(views))
        for name, _setter in views:
            self.combo_view.append_text(name)
        if views:
            self.combo_view.set_active(min(view, len(views) - 1))
        else:
            self.plot.clear()
        self.btn_export.set_sensitive(bool(views))
        cols, rows = result.get("table", ([], []))
        self.btn_save.set_sensitive(bool(cols))
        store = Gtk.ListStore(*([str] * len(cols))) if cols else None
        for r in rows:
            self._row_payload.append(r[len(cols)] if len(r) > len(cols) else None)
            store.append([str(v) for v in r[:len(cols)]])
        for k, c in enumerate(cols):
            col = Gtk.TreeViewColumn(c, Gtk.CellRendererText(), text=k)
            col.set_sort_column_id(k)
            col.set_resizable(True)
            self.table.append_column(col)
        self.table.set_model(store)
        self.text.get_buffer().set_text(result.get("summary", ""))
        if not views:
            self.notebook.set_current_page(1)
        self.label_hint.set_text(result.get("hint", ""))
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
        if r is None or self.run is None:
            return
        if r.get("map") and self.plot.kind == "heatmap" and "row" in info:
            key = self._map_pos.get((info["row"], info["col"]))
            if key in self.run.windows:
                self._toggle(key)
            return
        if r.get("pmf_click"):
            mode = getattr(self, "_pmf_click", None)
            if mode is None:
                return
            if mode[0] == "1d" and "x" in info:
                keys = [k for k in self.run.keys if self.run.target(k)]
                k = min(keys, key=lambda k: abs(self.run.target(k)[0] - info["x"]))
                self._toggle(k)
            elif mode[0] in ("2d", "mep") and self.plot.kind == "heatmap" and "row" in info:
                xs, ys = mode[1], mode[2]
                if mode[0] == "mep" and self._mep.get("click"):
                    self._mep[self._mep["click"]] = (info["row"], info["col"])
                    self._mep["click"] = None
                    self._mep["path"] = []
                    self._refresh(keep_view=True)
                    return
                if mode[0] == "2d":
                    p = np.array([xs[info["row"]], ys[info["col"]]])
                    keys = [k for k in self.run.keys if self.run.target(k)]
                    T = np.array([self.run.target(k) for k in keys])
                    scale = np.array([max(1e-9, np.ptp(T[:, 0])), max(1e-9, np.ptp(T[:, 1]))])
                    self._toggle(keys[int(np.argmin(np.linalg.norm((T - p) / scale, axis=1)))])

    def _table_rows(self):
        cols, rows = self.result.get("table", ([], []))
        return [r[:len(cols)] for r in rows]

    def _on_table_row(self, _view, path, _column):
        model = self.table.get_model()
        row = [model[path][k] for k in range(model.get_n_columns())]
        payload = None
        for k, r in enumerate(self._table_rows()):
            if self._row_payload[k] is not None and [str(v) for v in r] == row:
                payload = self._row_payload[k]; break
        if payload is None:
            return
        if payload[0] == "window":
            key = payload[1]
            keep = list(self.selected)
            self.selected = [key]
            try:
                self._safe(lambda: self._load("each"))
            finally:
                self.selected = keep
        elif payload[0] == "pair":
            for key in payload[1]:
                if key not in self.selected:
                    self.selected.append(key)
            self._update_selection_label()
            self._refresh(keep_view=True)

    # ------------------------------------------------------------- export
    def on_export_image(self, _b):
        path = self._ask_save("Export Image", ["png", "svg"], "umbrella.png")
        if path:
            self.plot.export(path)

    def on_save_table(self, _b):
        path = self._ask_save("Save Table", ["csv"], "umbrella.csv")
        if not path or self.result is None:
            return
        import csv
        cols, _rows = self.result.get("table", ([], []))
        with open(path, "w", newline="") as handle:
            w = csv.writer(handle)
            w.writerow(cols)
            w.writerows(self._table_rows())

    def _ask_save(self, title, extensions, name):
        chooser = Gtk.FileChooserDialog(title=title, parent=self.window, action=Gtk.FileChooserAction.SAVE)
        chooser.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        chooser.set_do_overwrite_confirmation(True)
        if self.run is not None:
            chooser.set_current_folder(self.run.folder)
        chooser.set_current_name(name)
        response, path = chooser.run(), chooser.get_filename()
        chooser.destroy()
        if response != Gtk.ResponseType.OK or not path:
            return None
        if os.path.splitext(path)[1].lower().lstrip(".") not in extensions:
            path += "." + extensions[0]
        return path

    def _safe(self, callback):
        try:
            callback()
        except Exception as error:
            traceback.print_exc()
            self._message("{}".format(error))

    def _message(self, text):
        self.label_hint.set_text(text)

    def _on_destroy(self, *_):
        self.Visible = False
        self.window = None
