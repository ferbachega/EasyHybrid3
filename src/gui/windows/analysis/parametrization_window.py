#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  EasyHybrid: Python interface for QM/MM and molecular simulations using pDynamo3
#  Module: DYFF Parametrization window
#
#  Copyright 2022-2026 Fernando Bachega
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  Description:
#      System > Force Field > Parametrization Report... (2026-10-23) - the data
#      of a DYFF modified parametrization of ligands (the dyff_param_... folder
#      of pdynamo/param_report.py: parametrization.json and the original xTB /
#      ORCA logs), in the style of the Trajectory Analysis window: a tree of
#      views on the left, a cairo plot (util/analysis_plot.py) with a table and
#      a text panel on the right, and actions on the 3D view (select the atoms
#      of a term, play a torsion scan). The same data are written as
#      report.html next to the json; "Open HTML Report" shows it.
#

import glob
import math
import os
import traceback

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

import numpy as np

from gui.widgets.custom_widgets import SystemComboBox
from util.analysis_plot import AnalysisPlot, PALETTE

BLACK, GREY, BLUE, RED, ORANGE, GREEN = (0.1, 0.1, 0.1), (0.55, 0.55, 0.55), (0.08, 0.40, 0.75), (0.80, 0.16, 0.16), \
                                        (0.90, 0.45, 0.05), (0.18, 0.55, 0.20)
KJ = 4.184


def _f(v, fmt="{:.3f}"):
    return "" if v is None else fmt.format(v)


def _wrap180(x):
    return ((x + 180.0) % 360.0) - 180.0


class ParametrizationWindow:

    def __init__(self, main=None):
        self.main = main
        self.home = main.home
        self.p_session = main.p_session
        self.vm_session = main.vm_session
        self.Visible = False
        self.window = None
        self.data = None
        self.result = None
        self._row_payload = []
        self._folders = []
        self._scan_ids = set()          # . systems made by "Play scan in 3D": never the subject of the window
        self._current_id = None
        self._matches = False

    # =====================================================================
    #  Window
    # =====================================================================
    def open_window(self, folder=None, system=None):
        if not self.Visible:
            self._build()
        else:
            self.window.present()
        if system is not None:
            self.combobox_systems.set_active_system(e_id=system.e_id)
        if folder:
            self._add_folder(folder, select=True)

    def _build(self):
        self.window = Gtk.Window(title="DYFF Parametrization")
        self.window.set_default_size(1260, 800)
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
        self.tree = Gtk.TreeView(model=self.store)
        self.tree.set_headers_visible(False)
        self.tree.append_column(Gtk.TreeViewColumn("View", Gtk.CellRendererText(), markup=0))
        self.tree.get_selection().connect("changed", self._on_tree_changed)
        sw = Gtk.ScrolledWindow()
        sw.set_policy(Gtk.PolicyType.NEVER, Gtk.PolicyType.AUTOMATIC)
        sw.add(self.tree)
        left.pack_start(sw, True, True, 0)
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
        self.btn_export = Gtk.Button(label="Export Image...")
        self.btn_export.connect("clicked", self.on_export_image)
        self.btn_save = Gtk.Button(label="Save Table...")
        self.btn_save.connect("clicked", self.on_save_table)
        bar.pack_end(self.btn_save, False, False, 0)
        bar.pack_end(self.btn_export, False, False, 0)
        right.pack_start(bar, False, False, 0)

        vpaned = Gtk.Paned(orientation=Gtk.Orientation.VERTICAL)
        self.plot = AnalysisPlot()
        self.plot.on_click = self._on_plot_click
        vpaned.pack1(self.plot, True, False)
        self.notebook = Gtk.Notebook()
        self.table = Gtk.TreeView()
        self.table.connect("row-activated", self._on_table_row)
        self.table.get_selection().connect("changed", self._on_table_select)
        sw = Gtk.ScrolledWindow()
        sw.add(self.table)
        self.notebook.append_page(sw, Gtk.Label(label="Table"))
        self.text = Gtk.TextView(editable=False, monospace=True, wrap_mode=Gtk.WrapMode.WORD_CHAR)
        sw = Gtk.ScrolledWindow()
        sw.add(self.text)
        self.notebook.append_page(sw, Gtk.Label(label="Details"))
        self.notebook.set_size_request(-1, 220)
        vpaned.pack2(self.notebook, False, False)
        vpaned.set_position(440)
        right.pack_start(vpaned, True, True, 0)
        self.actions_box = Gtk.Box(spacing=6)
        right.pack_start(self.actions_box, False, False, 0)
        paned.pack2(right, True, False)

        self.window.show_all()
        self.Visible = True
        self._show(None)
        if self.p_session.psystem:
            self.combobox_systems.set_active_system(e_id=self.p_session.active_id)
            self._on_system_changed(None)

    def _build_input_frame(self):
        frame = Gtk.Frame(label=" Parametrization ")
        grid = Gtk.Grid(column_spacing=8, row_spacing=6, margin=8)
        frame.add(grid)
        self.combobox_systems = SystemComboBox(self.main)
        self.combobox_systems.connect("changed", self._on_system_changed)
        self.combo_folder = Gtk.ComboBoxText()
        self.combo_folder.set_tooltip_text("The parametrizations of this system found in its working folder "
                                           "(dyff_param_<system>_<date>), most recent first")
        self.combo_folder.connect("changed", self._on_folder_changed)
        grid.attach(Gtk.Label(label="System", xalign=0), 0, 0, 1, 1)
        grid.attach(self.combobox_systems, 1, 0, 1, 1)
        grid.attach(Gtk.Label(label="Parametrization", xalign=0), 2, 0, 1, 1)
        grid.attach(self.combo_folder, 3, 0, 1, 1)
        b = Gtk.Button(label="Open Folder...")
        b.set_tooltip_text("Read another dyff_param_... folder (parametrization.json)")
        b.connect("clicked", self.on_open_folder)
        grid.attach(b, 4, 0, 1, 1)
        self.btn_html = Gtk.Button(label="Open HTML Report")
        self.btn_html.set_tooltip_text("report.html of this parametrization (same data, with the original logs), in the web browser")
        self.btn_html.connect("clicked", self.on_open_html)
        grid.attach(self.btn_html, 5, 0, 1, 1)
        self.btn_apply = Gtk.Button(label="Apply to System")
        self.btn_apply.set_tooltip_text("Defines DYFF on the system with these parameters, without computing them again "
                                        "(same atoms and bonds needed)")
        self.btn_apply.connect("clicked", self.on_apply)
        grid.attach(self.btn_apply, 6, 0, 1, 1)
        self.label_info = Gtk.Label(xalign=0)
        self.label_info.set_line_wrap(True)
        self.label_info.set_selectable(True)
        grid.attach(self.label_info, 0, 1, 7, 1)
        return frame

    # =====================================================================
    #  Sources
    # =====================================================================
    def _system(self):
        e_id = self.combobox_systems.get_system_id() if hasattr(self.combobox_systems, "get_system_id") else None
        if e_id is None:
            e_id = self.p_session.active_id
        return self.p_session.psystem.get(e_id)

    def _vobject(self):
        system = self._system()
        if system is None:
            return None
        return self.main.get_active_vobject(system.e_id)

    def _on_system_changed(self, _w):
        e_id = self.combobox_systems.get_system_id()
        if e_id in self._scan_ids and self._current_id is not None:
            # . a scan object was just added (and made active): stay on the parametrized system
            self.combobox_systems.handler_block_by_func(self._on_system_changed)
            self.combobox_systems.set_active_system(e_id=self._current_id)
            self.combobox_systems.handler_unblock_by_func(self._on_system_changed)
            return
        if e_id == self._current_id and self.data is not None:
            return
        self._current_id = e_id
        system = self._system()
        folders = []
        if system is not None:
            v = self._vobject()
            own = getattr(v, "dyff_param_folder", None) if v is not None else None
            if own:
                folders.append(own)
            base = getattr(system, "e_working_folder", None) or os.getcwd()
            safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in (system.label or "system"))[:60]
            folders += sorted(glob.glob(os.path.join(base, "dyff_param_{}_*".format(glob.escape(safe)))), reverse=True)
        self._set_folders(folders)

    def _set_folders(self, folders, select=None):
        seen, out = set(), []
        for f in folders:
            f = os.path.abspath(f)
            if f not in seen and os.path.isfile(os.path.join(f, "parametrization.json")):
                seen.add(f); out.append(f)
        self._folders = out
        self.combo_folder.handler_block_by_func(self._on_folder_changed)
        self.combo_folder.remove_all()
        for f in out:
            self.combo_folder.append_text(os.path.basename(f))
        self.combo_folder.handler_unblock_by_func(self._on_folder_changed)
        if out:
            self.combo_folder.set_active(out.index(os.path.abspath(select)) if select and os.path.abspath(select) in out else 0)
        else:
            self._load(None)

    def _add_folder(self, folder, select=True):
        folder = os.path.abspath(folder if os.path.isdir(folder) else os.path.dirname(folder))
        self._set_folders([folder] + self._folders, select=folder if select else None)

    def _on_folder_changed(self, combo):
        k = combo.get_active()
        self._load(self._folders[k] if 0 <= k < len(self._folders) else None)

    def on_open_folder(self, _b):
        chooser = Gtk.FileChooserDialog(title="DYFF parametrization (parametrization.json)", parent=self.window,
                                        action=Gtk.FileChooserAction.OPEN)
        chooser.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
        flt = Gtk.FileFilter(); flt.set_name("parametrization.json"); flt.add_pattern("parametrization.json")
        chooser.add_filter(flt)
        system = self._system()
        base = getattr(system, "e_working_folder", None) if system is not None else None
        if base and os.path.isdir(base):
            chooser.set_current_folder(base)
        response, path = chooser.run(), chooser.get_filename()
        chooser.destroy()
        if response == Gtk.ResponseType.OK and path:
            self._add_folder(path)

    def on_open_html(self, _b):
        if self.data is None:
            return
        path = os.path.join(self.data["_folder"], "report.html")
        if not os.path.exists(path):
            from pdynamo import param_report as PR
            data = dict(self.data); data.pop("_folder", None)
            PR.write(self.data["_folder"], data)
        self.main._open_in_browser(path)

    def on_apply(self, _b):
        if self.data is None:
            return
        try:
            self.main.load_dyff_parametrization(self.data["_folder"], self._system())
            self._message("DYFF is being defined on '{}' with these parameters; the result appears in the Status panel.".format(
                self._system().label))
        except Exception as error:
            self._message("Apply to System: {}".format(error))

    # =====================================================================
    #  Data -> tree of views
    # =====================================================================
    def _load(self, folder):
        from pdynamo import param_report as PR
        self.data = None
        self.store.clear()
        if folder is not None:
            try:
                self.data = PR.load(folder)
            except Exception as error:
                self.label_info.set_text("Could not read {}: {}".format(folder, error))
        has = self.data is not None
        self.btn_html.set_sensitive(has)
        self.btn_apply.set_sensitive(has)
        if not has:
            if folder is None:
                self.label_info.set_markup("<i>No parametrization of this system. Run System &gt; Force Field &gt; Define Force Field &gt; "
                                           "Prepare DYFF System (DYFF modified), or open a dyff_param_... folder.</i>")
            self._show(None)
            return
        d = self.data
        d["_by_index"] = {a["index"]: a for a in d["atoms"]}
        rep = d.get("report") or {}
        system = self._system()
        from pdynamo import param_report as PR
        self._matches = system is not None and PR.check_system(d, system) is None
        info = "{}  ·  {}  ·  EasyHybrid {}  ·  {} atoms outside the library of {}  ·  total charge {:+.3f}".format(
            d["system"]["name"], d["created"].replace("T", " "), d["easyhybrid_version"], len(d["rest_atoms"]),
            d["system"]["atoms"], round(d["system"]["total_charge"], 3) + 0.0)
        if d.get("total_seconds") is not None:
            info += "  ·  {:.0f} s".format(d["total_seconds"])
        if not self._matches:
            info += "\n(This parametrization is not of the selected system: selecting atoms in the 3D view is off.)"
        self.label_info.set_text(info)
        add = lambda parent, text, key: self.store.append(parent, [text, key])
        p = add(None, "<b>Overview</b>", "")
        add(p, "Summary", "summary")
        p = add(None, "<b>Charges</b>", "")
        add(p, "Charges per atom", "charges")
        if any(a.get("resp") is not None for a in d["atoms"]):
            add(p, "CM5 against RESP", "parity")
        if rep.get("ligand_charges"):
            add(p, "Molecules (RESP fits)", "resp")
        p = add(None, "<b>Bonded terms (xTB Hessian)</b>", "")
        add(p, "Bonds ({})".format(len(d["edits"].get("bond_terms") or {})), "bonds")
        add(p, "Angles ({})".format(len(d["edits"].get("angle_terms") or {})), "angles")
        p = add(None, "<b>Torsions</b>", "")
        add(p, "Library torsions ({})".format(len(rep.get("library_torsions") or [])), "library")
        details = rep.get("scan_details") or []
        if "scan_targets" in rep:
            add(p, "Scanned torsions: overview ({} of {} fitted)".format(sum(1 for x in details if x.get("accepted")),
                                                                         rep.get("scan_targets", 0)), "scans")
            for k, s in enumerate(details):
                j, kk = (int(x) for x in s["bond"].split("-"))
                mark = "" if s.get("accepted") else ("  <span foreground='#c62828'>(rule kept)</span>" if not s.get("error")
                                                     else "  <span foreground='#c62828'>(failed)</span>")
                add(p, "{} – {}{}".format(self._label(j), self._label(kk), mark), "scan:{}".format(k))
        p = add(None, "<b>Original logs</b>", "")
        add(p, "xTB / ORCA files", "logs")
        self.tree.expand_all()
        self.tree.get_selection().select_iter(self.store.iter_children(self.store.get_iter_first()))

    def _short(self, index):
        """ Element + 1-based number (plot axes). """
        a = (self.data or {}).get("_by_index", {}).get(index)
        return "{}{}".format(a["element"] if a else "", index + 1)

    def _label(self, index):
        a = (self.data or {}).get("_by_index", {}).get(index)
        return "{} ({})".format(a["name"], index + 1) if a else str(index + 1)

    def _on_tree_changed(self, selection):
        model, it = selection.get_selected()
        if it is None or self.data is None:
            return
        key = model[it][1]
        if not key:
            return
        try:
            self._show(self._view(key))
        except Exception as error:
            traceback.print_exc()
            self._show({"views": [], "summary": "Could not show this view: {}".format(error)})

    # =====================================================================
    #  Views
    # =====================================================================
    def _view(self, key):
        if key.startswith("scan:"):
            return self._view_scan(int(key.split(":")[1]))
        return getattr(self, "_view_" + key)()

    def _view_summary(self):
        d = self.data
        rep = d.get("report") or {}
        timings = rep.get("timings") or {}
        names = {"cm5": "CM5", "resp": "RESP", "hessian": "Hessians", "torsions": "torsions"}
        keys = list(timings)
        rows = [["system", d["system"]["name"]], ["created", d["created"]], ["EasyHybrid", d["easyhybrid_version"]],
                ["atoms in the system", d["system"]["atoms"]], ["atoms outside the residue library", len(d["rest_atoms"])],
                ["total charge", "{:+.3f}".format(round(d["system"]["total_charge"], 3) + 0.0)],
                ["residue library", "{} of {} residues".format((d.get("library") or {}).get("matched"), (d.get("library") or {}).get("residues"))],
                ["Hessian bonds / angles", "{} / {}".format(rep.get("bonds", 0), rep.get("angles", 0))],
                ["library torsions", rep.get("torsions", 0)],
                ["scanned torsions fitted", "{} of {}".format(sum(1 for x in rep.get("scan_details") or [] if x.get("accepted")), rep.get("scan_targets", 0))]]
        for k, v in (d.get("settings") or {}).items():
            rows.append(["setting: " + k.replace("_", " "), str(v)])
        for k, v in (rep.get("charge_settings") or {}).items():
            rows.append(["charges: " + k.replace("_", " "), str(v)])
        for k, v in (rep.get("scan_settings") or {}).items():
            rows.append(["torsion scans: " + k.replace("_", " "), str(v)])
        for k in keys:
            rows.append(["time: " + names.get(k, k), "{:.1f} s".format(timings[k])])
        views = []
        if keys:
            views.append(("Wall time per step", lambda p: p.set_lines(
                [{"x": list(range(len(keys))), "y": [timings[k] for k in keys], "style": "bars", "color": BLUE}],
                "", "wall time (s)", "Wall time of each step ({:.0f} s in all)".format(d.get("total_seconds") or sum(timings.values())),
                ymin=0, xlabels=[names.get(k, k) for k in keys])))
        return {"views": views, "table": (["item", "value"], rows),
                "summary": "\n\n".join(d.get("messages") or []) or "No messages.",
                "hint": "Everything recorded for this parametrization. The Details tab has the message shown in the Status panel."}

    # ------------------------------------------------------------- charges
    def _view_charges(self):
        atoms = self.data["atoms"]
        labels = [self._short(a["index"]) for a in atoms]
        x = list(range(len(atoms)))
        has_resp = any(a.get("resp") is not None for a in atoms)
        num = lambda v: float("nan") if v is None else v

        def bars(p, keys):
            series = [{"x": x, "y": [num(a.get(k)) for a in atoms], "style": "bars", "label": lab, "color": c} for k, lab, c in keys]
            p.set_lines(series, "atom", "charge (e)", "Partial charges of the atoms outside the library", xlabels=labels, hline=0.0)

        views = []
        if has_resp:
            views.append(("CM5 and RESP-2", lambda p: bars(p, [("cm5", "CM5", GREY), ("resp", "RESP-2", BLUE)])))
            views.append(("RESP-2 minus CM5", lambda p: p.set_lines(
                [{"x": x, "y": [num(a.get("resp")) - num(a.get("cm5")) for a in atoms], "style": "bars", "color": ORANGE}],
                "atom", "RESP-2 − CM5 (e)", "Difference between the charge models", xlabels=labels, hline=0.0)))
            views.append(("ESP (unrestrained) and RESP-2", lambda p: bars(p, [("esp", "ESP", ORANGE), ("resp", "RESP-2", BLUE)])))
        views.append(("Final charges", lambda p: bars(p, [("charge", "final", BLUE)])))
        views.append(("CM5", lambda p: bars(p, [("cm5", "CM5", GREY)])))
        rows = [[a["index"] + 1, a["name"], a["residue"], a["element"], a["formal_charge"] or "", _f(a.get("cm5"), "{:.4f}"),
                 _f(a.get("esp"), "{:.4f}"), _f(a.get("resp"), "{:.4f}"), _f(a.get("charge"), "{:.4f}"), [a["index"]]] for a in atoms]
        by_el = {}
        for a in atoms:
            by_el.setdefault(a["element"], []).append(a)
        text = ["Charge sums by element (final):"] + ["  {:2s} {:3d} atoms  sum {:+.3f}  mean {:+.3f}".format(
            e, len(v), sum(a["charge"] or 0 for a in v), sum(a["charge"] or 0 for a in v) / len(v)) for e, v in sorted(by_el.items())]
        if has_resp:
            diff = [abs(a["resp"] - a["cm5"]) for a in atoms if a.get("resp") is not None and a.get("cm5") is not None]
            if diff:
                text.append("\nRESP-2 vs CM5: mean |difference| {:.3f} e, largest {:.3f} e.".format(sum(diff) / len(diff), max(diff)))
        return {"views": views, "table": (["#", "atom", "residue", "element", "formal", "CM5", "ESP", "RESP-2", "final"], rows),
                "row_action": "atoms", "click_atoms": [[a["index"]] for a in atoms], "summary": "\n".join(text),
                "hint": "Double-click a row (or click a bar) to select the atom in the 3D view.",
                "actions": [("Color by charge", lambda: self._color_by_charge())] if self._matches else []}

    def _view_parity(self):
        atoms = [a for a in self.data["atoms"] if a.get("resp") is not None and a.get("cm5") is not None]
        x = [a["resp"] for a in atoms]; y = [a["cm5"] for a in atoms]
        n = len(x)
        mae = sum(abs(a - b) for a, b in zip(x, y)) / n if n else 0.0
        slope = float(np.polyfit(x, y, 1)[0]) if n > 1 else float("nan")
        Z = {"H": 1, "C": 6, "N": 7, "O": 8, "S": 16, "P": 15}
        c = [Z.get(a["element"], 30) for a in atoms]
        return {"views": [("CM5 against RESP-2", lambda p: p.set_scatter(x, y, c, "RESP-2 (e)", "CM5 (e)",
                                                                        "MAE {:.3f} e, slope {:.2f}".format(mae, slope), "atomic number"))],
                "table": (["#", "atom", "RESP-2", "CM5", "CM5 − RESP-2"], [[a["index"] + 1, a["name"], "{:.4f}".format(a["resp"]),
                          "{:.4f}".format(a["cm5"]), "{:+.4f}".format(a["cm5"] - a["resp"]), [a["index"]]] for a in atoms]),
                "row_action": "atoms", "click_atoms": [[a["index"]] for a in atoms],
                "summary": "{} atoms. A slope below 1 means CM5 is less polar than RESP-2.".format(n),
                "hint": "Click a point to select that atom in the 3D view."}

    def _view_resp(self):
        rep = self.data.get("report") or {}
        groups = rep.get("ligand_charges") or []
        rows = []
        for g in groups:
            first = g["members"][0] if g.get("members") else None
            rows.append([self._label(first) if first is not None else "", g["atoms"], g["status"], g.get("engine") or "",
                         g.get("charge", ""), g.get("multiplicity", ""), _f(g.get("rrms"), "{:.4f}"), _f(g.get("dipole_qm"), "{:.2f}"),
                         _f(g.get("dipole_resp"), "{:.2f}"), g.get("points") or "", g.get("seconds", ""), g.get("members") or []])
        ok = [g for g in groups if g.get("rrms") is not None]
        views = []
        if ok:
            views.append(("Fit quality (RRMS)", lambda p: p.set_lines(
                [{"x": list(range(len(ok))), "y": [g["rrms"] for g in ok], "style": "bars", "color": BLUE}],
                "", "RRMS of the RESP-2 fit", "Relative RMS of the fitted potential", ymin=0,
                xlabels=[self._short(g["members"][0]) for g in ok])))
            views.append(("Dipole: QM and RESP-2", lambda p: p.set_lines(
                [{"x": list(range(len(ok))), "y": [g.get("dipole_qm") or float("nan") for g in ok], "style": "bars", "label": "QM", "color": BLACK},
                 {"x": list(range(len(ok))), "y": [g.get("dipole_resp") or float("nan") for g in ok], "style": "bars", "label": "RESP-2", "color": BLUE}],
                "", "dipole (D)", "Dipole moment of the QM density and of the charges", ymin=0,
                xlabels=[self._short(g["members"][0]) for g in ok])))
        settings = rep.get("charge_settings") or {}
        return {"views": views, "table": (["molecule (first atom)", "atoms", "charges", "engine", "charge", "mult.", "RRMS",
                                           "dipole QM (D)", "dipole RESP (D)", "points", "time (s)"], rows),
                "row_action": "atoms", "summary": "Ligand charge model: {}. Whole molecules get two-stage RESP (AMBER) on the "
                "electrostatic potential at {} Merz–Kollman points; groups bonded to the rest of the system and metal sites keep "
                "CM5 ({}). The original xtb / ORCA outputs are under Original logs (charges/resp_N).".format(
                    settings.get("ligand_charges", "cm5"), "the", settings.get("cm5", "GFN1-xTB CM5")),
                "hint": "Double-click a row to select the molecule."}

    # ------------------------------------------------------------- bonded
    def _sorted_terms(self, kind):
        terms = self.data["edits"].get(kind) or {}
        return sorted(terms.items(), key=lambda t: [int(x) for x in t[0].split("-")])

    def _view_bonds(self):
        terms = self._sorted_terms("bond_terms")
        labels = ["–".join(self._short(int(x)) for x in k.split("-")) for k, _ in terms]
        x = list(range(len(terms)))
        rows = [[" – ".join(self._label(int(i)) for i in k.split("-")), "{:.4f}".format(v[0]), "{:.1f}".format(v[1]),
                 [int(i) for i in k.split("-")]] for k, v in terms]
        return {"views": [("Force constants", lambda p: p.set_lines([{"x": x, "y": [v[1] for _, v in terms], "style": "bars", "color": BLUE}],
                                                                    "bond", "k (kJ/mol/Å²)", "Bond force constants (modified Seminario), E = k (b − b0)²",
                                                                    ymin=0, xlabels=labels)),
                          ("Equilibrium lengths", lambda p: p.set_lines([{"x": x, "y": [v[0] for _, v in terms], "style": "points", "color": BLUE}],
                                                                        "bond", "b0 (Å)", "Bond lengths at the xTB minimum", xlabels=labels))],
                "table": (["bond", "b0 (Å)", "k (kJ/mol/Å²)"], rows), "row_action": "atoms",
                "click_atoms": [r[-1] for r in rows],
                "summary": self._hessian_summary(), "hint": "Double-click a row (or click a bar) to select the bond."}

    def _view_angles(self):
        from pdynamo.param_report import _angle_from_pairs
        terms = self._sorted_terms("angle_terms")
        values = [_angle_from_pairs(v) for _, v in terms]
        labels = ["–".join(self._short(int(x)) for x in k.split("-")) for k, _ in terms]
        x = list(range(len(terms)))
        nan = float("nan")
        rows = [[" – ".join(self._label(int(i)) for i in k.split("-")), _f(t0, "{:.2f}"), _f(kk, "{:.1f}"),
                 "  ".join("{:+.2f}×{}".format(c, int(n)) for c, n in v), [int(i) for i in k.split("-")]]
                for (k, v), (t0, kk) in zip(terms, values)]
        return {"views": [("Force constants", lambda p: p.set_lines([{"x": x, "y": [kk if kk is not None else nan for _, kk in values],
                                                                      "style": "bars", "color": GREEN}], "angle", "k (kJ/mol/rad²)",
                                                                    "Angle force constants (modified Seminario)", ymin=0, xlabels=labels)),
                          ("Equilibrium angles", lambda p: p.set_lines([{"x": x, "y": [t0 if t0 is not None else nan for t0, _ in values],
                                                                         "style": "points", "color": GREEN}], "angle", "θ0 (degrees)",
                                                                       "Angles at the xTB minimum", xlabels=labels))],
                "table": (["angle", "θ0 (°)", "k (kJ/mol/rad²)", "DYFF series c×n (kJ/mol)"], rows), "row_action": "atoms",
                "click_atoms": [r[-1] for r in rows],
                "summary": self._hessian_summary(), "hint": "Double-click a row (or click a bar) to select the angle."}

    def _hessian_summary(self):
        rep = self.data.get("report") or {}
        text = ("Bonds and angles of the atoms outside the library from the GFN2-xTB Hessian of {} capped group(s) "
                "(largest subsystem {} atoms), modified Seminario method; {} term(s) with a non-physical constant were "
                "left to the DYFF rule.").format(rep.get("groups", 0), rep.get("largest_subsystem", 0), rep.get("skipped", 0))
        if rep.get("calibrated"):
            text += "\nDFT calibration by chemical environment: {}.".format(rep["calibrated"])
        return text + "\nOriginal xtb --ohess outputs: Original logs > hessian/group_N."

    # ------------------------------------------------------------- torsions
    @staticmethod
    def _series_curve(pairs, step=2.0):
        phi = np.arange(-180.0, 180.0 + step / 2, step)
        e = sum(c * np.cos(n * np.radians(phi)) for c, n in pairs) if pairs else np.zeros_like(phi)
        e = (e - e.min()) / KJ
        return list(phi), list(e)

    def _view_library(self):
        rep = self.data.get("report") or {}
        dih = self.data["edits"].get("dihedral_bonds") or {}
        keys = rep.get("library_torsions") or []
        checks = rep.get("torsion_checks") or {}
        rows = []
        for k in keys:
            c = checks.get(k)
            rows.append([" – ".join(self._label(int(i)) for i in k.split("-")),
                         "{:.2f} / {:.2f}".format(c["rule_rmse"], c["torsion_rmse"]) if c and c.get("torsion_rmse") is not None
                         else (c.get("error", "") if c else "not checked"),
                         _f(c.get("barrier_xtb"), "{:.2f}") if c else "",
                         "  ".join("{:+.3f}×{}".format(x, int(n)) for x, n in dih.get(k, [])), ("library", k)])
        self._library_keys = keys
        actions = []
        if self._matches and keys:
            actions = [("Check selected with xTB scan", lambda: self._check_torsions(self._selected_library_keys())),
                       ("Check all with xTB scan", lambda: self._check_torsions([k for k in keys if k not in checks]))]
        done = [checks[k] for k in keys if k in checks and checks[k].get("torsion_rmse") is not None]
        text = ("{} of {} rotatable bonds have a torsion from the DFT torsion library ({}). The series replaces every DYFF dihedral "
                "term about the bond and is applied to each i-j-k-l quadruple.\n\nThe library torsions were fitted to DFT, not "
                "to xTB. 'Check with xTB scan' runs the same relaxed GFN2-xTB scan as the automatic fits (gas-phase capped "
                "fragment, ~10-60 s per bond) and compares DYFF with the DYFF rule and with the library torsion. Click a "
                "bond to see it alone.").format(rep.get("torsions", 0), rep.get("rotatable", 0),
                                                ", ".join("{} {}".format(v, k) for k, v in (rep.get("levels") or {}).items()) or "none")
        if done:
            text += "\n\nChecked: {} bond(s); mean RMSE against xTB: DYFF rule {:.2f}, library torsion {:.2f} kcal/mol.".format(
                len(done), sum(c["rule_rmse"] for c in done) / len(done), sum(c["torsion_rmse"] for c in done) / len(done))
        views = [("All library series V(φ)", self._plot_all_library)] if keys else []
        if len(done) > 1:
            views.append(("xTB check: RMSE rule and library", lambda p: self._plot_check_bars(p, keys, checks)))
        return {"views": views,
                "table": (["bond", "xTB check: RMSE rule / library", "xTB barrier", "series c×n (kJ/mol), each i-j-k-l"], rows),
                "row_action": "library", "summary": text, "actions": actions,
                "hint": "Click a bond to isolate it (with its xTB check when done). Energies in kcal/mol."}

    def _plot_all_library(self, p):
        dih = self.data["edits"].get("dihedral_bonds") or {}
        keys = self._library_keys
        series = []
        for m, k in enumerate(keys[:10]):
            phi, e = self._series_curve(dih.get(k, []))
            series.append({"x": phi, "y": e, "label": "–".join(self._short(int(i)) for i in k.split("-")), "color": PALETTE[m % len(PALETTE)]})
        p.set_lines(series, "φ (degrees)", "V − Vmin (kcal/mol)", "DFT torsion library: one series per bond" + (" (first 10)" if len(keys) > 10 else ""))

    def _plot_check_bars(self, p, keys, checks):
        done = [k for k in keys if k in checks and checks[k].get("torsion_rmse") is not None]
        x = list(range(len(done)))
        p.set_lines([{"x": x, "y": [checks[k]["rule_rmse"] for k in done], "style": "bars", "label": "DYFF rule", "color": GREY},
                     {"x": x, "y": [checks[k]["torsion_rmse"] for k in done], "style": "bars", "label": "library torsion", "color": GREEN}],
                    "bond", "RMSE vs xTB (kcal/mol)", "Library torsions against xTB scans", ymin=0,
                    xlabels=["–".join(self._short(int(i)) for i in k.split("-")) for k in done])

    def _isolate_library(self, key):
        rep = self.data.get("report") or {}
        dih = self.data["edits"].get("dihedral_bonds") or {}
        c = (rep.get("torsion_checks") or {}).get(key)
        name = " – ".join(self._label(int(i)) for i in key.split("-"))
        if c and c.get("profiles"):
            views, *_ = self._profile_views(
                "{}  ·  RMSE vs xTB: rule {:.2f} → library {:.2f} kcal/mol".format(name, c["rule_rmse"], c["torsion_rmse"]),
                c["quad"], c["profiles"]["phi"], c["profiles"]["xtb"],
                [("DYFF rule", c["profiles"]["dyff_rule"], GREY, True), ("DYFF + library torsion", c["profiles"]["dyff_torsion"], GREEN, False)],
                dih.get(key, []), "Library torsion term")
        else:
            phi, e = self._series_curve(dih.get(key, []))
            views = [("Library torsion V(φ)", lambda p: p.set_lines([{"x": phi, "y": e, "label": name, "color": GREEN}], "φ (degrees)",
                                                                   "V − Vmin (kcal/mol)", "Library torsion {} (not checked against xTB yet)".format(name), ymin=0))]
        views.append(("All library series V(φ)", self._plot_all_library))
        self._set_views(views)
        if self._matches:
            self._select_atoms(c["quad"] if c else [int(i) for i in key.split("-")])

    def _selected_library_keys(self):
        model, it = self.table.get_selection().get_selected()
        if it is None:
            raise ValueError("Select a bond in the table first.")
        payload = self._payload_of([model[it][k] for k in range(model.get_n_columns())])
        return [payload[1]] if payload else []

    def _view_scans(self):
        rep = self.data.get("report") or {}
        det = rep.get("scan_details") or []
        ok = [s for s in det if s.get("before_rmse") is not None]
        labels = ["{}–{}".format(*(self._short(int(i)) for i in s["bond"].split("-"))) for s in ok]
        x = list(range(len(ok)))
        nan = float("nan")
        rows = []
        for s in det:
            j, k = (int(i) for i in s["bond"].split("-"))
            status = ("fitted" + (" (RMSE > {} kcal/mol)".format((rep.get("scan_settings") or {}).get("poor_rmse", 2.0)) if s.get("poor") else "")
                      if s.get("accepted") else ("failed: " + s["error"] if s.get("error") else "rejected, rule kept"))
            rows.append(["{} – {}".format(self._label(j), self._label(k)), " – ".join(str(a + 1) for a in s.get("quad", [])),
                         s.get("fragment_atoms", ""), _f(s.get("barrier_xtb"), "{:.2f}"), _f(s.get("before_rmse"), "{:.2f}"),
                         _f(s.get("fit_rmse"), "{:.2f}"), status, ("scan", det.index(s))])
        ss = rep.get("scan_settings") or {}
        failed = [f for f in rep.get("scan_failed") or [] if not any(f.startswith(s["bond"] + ":") for s in det)]
        text = ("{} flexible bond(s) without a library torsion. Relaxed GFN2-xTB scan in {:g}° steps on {}; DYFF relaxed at each "
                "point with a 1000 kJ/mol/rad² restraint; self-consistent cosine fit (n = 1-4). A fit is kept only when its RMSE is "
                "below {:.0%} of the rule's.").format(rep.get("scan_targets", 0), ss.get("step", 15.0),
                                                      "a capped fragment per bond ({} bonds around it{})".format(ss.get("fragment_shell", 2),
                                                      ", gas phase" if ss.get("gas_phase") else "") if ss.get("fragments", True) else
                                                      "the subsystem of the group", ss.get("min_improvement", 0.9))
        if failed:
            text += "\nNot scanned: " + "; ".join(failed)
        fit_colors = [BLUE if s.get("accepted") else RED for s in ok]
        return {"views": [("RMSE: DYFF rule and fit", lambda p: p.set_lines(
                    [{"x": x, "y": [s["before_rmse"] for s in ok], "style": "bars", "label": "DYFF rule", "color": GREY},
                     {"x": x, "y": [s.get("fit_rmse") or nan for s in ok], "style": "bars", "label": "fit (red: rejected)",
                      "color": BLUE, "colors": fit_colors}],
                    "bond", "RMSE vs xTB (kcal/mol)", "Torsion fits (rejected: {})".format(
                        ", ".join(l for l, s in zip(labels, ok) if not s.get("accepted")) or "none"), ymin=0, xlabels=labels)),
                          ("xTB barrier", lambda p: p.set_lines([{"x": x, "y": [s.get("barrier_xtb") or nan for s in ok], "style": "bars",
                                                                  "color": BLACK}], "bond", "barrier (kcal/mol)",
                                                                "Height of the xTB scan profile", ymin=0, xlabels=labels))] if ok else [],
                "table": (["bond", "dihedral", "fragment atoms", "xTB barrier", "RMSE rule", "RMSE fit", "result"], rows),
                "row_action": "scans", "click_atoms": [list(s.get("quad", [])) for s in ok], "summary": text,
                "hint": "Energies in kcal/mol. Click a bond to isolate its fit against the xTB scan (each bond also has its own "
                        "entry in the tree)."}

    def _view_scan(self, k):
        rep = self.data.get("report") or {}
        s = (rep.get("scan_details") or [])[k]
        j, kk = (int(i) for i in s["bond"].split("-"))
        title = "{} – {}".format(self._label(j), self._label(kk))
        quad = list(s.get("quad", []))
        p = s.get("profiles") or {}
        actions = []
        if self._matches and quad:
            actions.append(("Select dihedral atoms", lambda: self._select_atoms(quad)))
        logdir = self._scan_logdir(s)
        if logdir and os.path.exists(os.path.join(logdir, "xtbscan.log")):
            actions.append(("Play scan in 3D", lambda: self._play_scan(s)))
        if logdir and os.path.exists(os.path.join(logdir, "xtb_scan.out")):
            actions.append(("Show xtb scan output", lambda: self._show_file(os.path.join(logdir, "xtb_scan.out"))))
        info = ["Bond {}  ·  dihedral {}".format(title, " – ".join(self._label(a) for a in quad)),
                "Result: " + ("fitted and kept" if s.get("accepted") else ("failed: " + s["error"] if s.get("error") else
                                                                          "the fit does not improve on the DYFF rule; the rule is kept")),
                "RMSE against xTB: rule {} / fit {} kcal/mol;  xTB barrier {} kcal/mol".format(
                    _f(s.get("before_rmse"), "{:.2f}"), _f(s.get("fit_rmse"), "{:.2f}"), _f(s.get("barrier_xtb"), "{:.2f}")),
                "Scan system: {} atoms, charge {}, {}".format(s.get("fragment_atoms", "?"), s.get("fragment_charge", "?"),
                                                               "gas phase" if not s.get("scan_solvent") else "ALPB(" + str(s["scan_solvent"]) + ")")]
        if s.get("pairs"):
            info.append("Fitted series (replaces every DYFF dihedral about the bond, each i-j-k-l):  " +
                        "  ".join("{:+.3f} cos({}φ)".format(c, int(n)) for c, n in s["pairs"]) + "  kJ/mol")
        if s.get("dropped_periods"):
            info.append("Periods dropped (they cancel over the quadruples): " + ", ".join(str(x) for x in s["dropped_periods"]))
        if logdir:
            info.append("Logs: " + logdir)
        if not p.get("phi"):
            return {"views": [], "summary": "\n".join(info), "actions": actions}
        fit_color = BLUE if s.get("accepted") else RED
        fit_label = "DYFF fitted" if s.get("accepted") else "DYFF fitted (rejected)"
        views, xs, order, xtb, (rule, fit) = self._profile_views(
            title + "  ·  RMSE rule {} → fit {} kcal/mol".format(_f(s.get("before_rmse"), "{:.2f}"), _f(s.get("fit_rmse"), "{:.2f}")),
            quad, p["phi"], p["xtb"], [("DYFF rule", p["dyff_before"], GREY, True), (fit_label, p["dyff_after"], fit_color, not s.get("accepted"))],
            s.get("pairs"), "Fitted torsion term")
        rows = [["{:.1f}".format(xs[m]), "{:.2f}".format(xtb[m]), "{:.2f}".format(rule[m]), "{:.2f}".format(fit[m]),
                 "{:+.2f}".format(fit[m] - xtb[m]), order[m] + 1, ("frame", order[m])] for m in range(len(xs))]
        self._scan_order = order
        return {"views": views,
                "table": (["φ (°)", "xTB", "DYFF rule", "DYFF fit", "fit − xTB", "scan point"], rows),
                "row_action": "scan", "summary": "\n".join(info), "actions": actions, "scan": s,
                "hint": "Relative energies (kcal/mol). 'Play scan in 3D' opens the xTB scan as a new object; then click a point "
                        "of the profile (or double-click a row) to go to that geometry."}

    def _profile_views(self, title, quad, phi_raw, xtb_raw, curves, pairs, term_title):
        """ Views of a torsion profile against the xTB scan: curves = [ ( label, values, color, dashed ) ], the
            last one is the model judged (errors, point-by-point). Returns ( views, xs, order, xtb, curves sorted ). """
        phi = [_wrap180(v) for v in phi_raw]
        order = sorted(range(len(phi)), key=lambda i: phi[i])
        xs = [phi[i] for i in order]
        xtb = [xtb_raw[i] for i in order]
        sorted_curves = [[values[i] for i in order] for _l, values, _c, _d in curves]
        label, _v, color, _d = curves[-1]
        model = sorted_curves[-1]
        rmse = lambda c: math.sqrt(sum((a - b) ** 2 for a, b in zip(c, xtb)) / len(xtb))
        low = [m for m in range(len(xtb)) if xtb[m] < 5.0]
        rmse_low = lambda c: math.sqrt(sum((c[m] - xtb[m]) ** 2 for m in low) / len(low)) if low else float("nan")
        xlabel = "dihedral {} (degrees)".format("-".join(str(a + 1) for a in quad))

        def profiles(pl):
            series = [{"x": xs, "y": xtb, "label": "GFN2-xTB scan", "color": BLACK, "width": 2.0},
                      {"x": xs, "y": xtb, "style": "points", "color": BLACK}]
            for (lab, _v, col, dash), ys in zip(curves, sorted_curves):
                series.append({"x": xs, "y": ys, "label": lab, "color": col, "dash": dash, "width": 2.2 if ys is model else 1.6})
            pl.set_lines(series, xlabel, "relative energy (kcal/mol)", title, ymin=0)

        def model_only(pl):
            pl.set_lines([{"x": xs, "y": xtb, "label": "GFN2-xTB scan", "color": BLACK, "width": 2.0},
                          {"x": xs, "y": xtb, "style": "points", "color": BLACK},
                          {"x": xs, "y": model, "label": label, "color": color, "width": 2.4}],
                         xlabel, "relative energy (kcal/mol)", "{} against the xTB scan  ·  RMSE {:.2f} kcal/mol (below 5 kcal/mol: {:.2f})".format(
                             label, rmse(model), rmse_low(model)), ymin=0)

        def residuals(pl):
            pl.set_lines([{"x": xs, "y": [a - b for a, b in zip(ys, xtb)], "label": "{} − xTB".format(lab), "color": col, "dash": dash}
                          for (lab, _v, col, dash), ys in zip(curves, sorted_curves)],
                         xlabel, "error (kcal/mol)", "Error along the scan", hline=0.0)

        def correction(pl):
            cphi, ce = self._series_curve(pairs or [])
            pl.set_lines([{"x": cphi, "y": ce, "label": "one i-j-k-l quadruple", "color": color}],
                         "φ of the quadruple (degrees)", "V − Vmin (kcal/mol)", term_title + ": " +
                         (" + ".join("{:.2f} cos({}φ)".format(c / KJ, int(n)) for c, n in (pairs or [])) or "none") + " (kcal/mol)", ymin=0)

        def scatter(pl):
            pl.set_scatter(xtb, model, None, "xTB (kcal/mol)", label + " (kcal/mol)", "Point by point against the xTB scan")

        views = [("Profiles: xTB scan and DYFF", profiles), ("Fit quality: {} vs xTB".format(label), model_only),
                 ("Errors along the scan", residuals), (term_title, correction), ("Point by point", scatter)]
        return views, xs, order, xtb, sorted_curves

    def _scan_logdir(self, s):
        name = s.get("log_dir")
        if not name or self.data is None:
            return None
        path = os.path.join(self.data["_folder"], "logs", "torsions", name)
        return path if os.path.isdir(path) else None

    # ------------------------------------------------------------- logs
    def _view_logs(self):
        root = os.path.join(self.data["_folder"], "logs")
        rows = []
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames.sort()
            for name in sorted(filenames):
                path = os.path.join(dirpath, name)
                rel = os.path.relpath(path, root)
                kind = {".out": "program output", ".xyz": "geometry", ".inp": "input", ".log": "trajectory (xyz frames)"}.get(
                    os.path.splitext(name)[1], "")
                rows.append([os.path.dirname(rel), name, kind, "{:.1f}".format(os.path.getsize(path) / 1024.0), ("file", path)])
        return {"views": [], "table": (["folder", "file", "kind", "kB"], rows), "row_action": "file",
                "summary": "Double-click a file to show it here.\n\nFolder: {}".format(root),
                "hint": "The inputs and outputs of xtb and ORCA as the programs wrote them.",
                "actions": [("Open logs folder", lambda: self.main._open_in_browser(root))]}

    def _show_file(self, path, limit=3000000):
        with open(path, errors="replace") as handle:
            text = handle.read(limit)
        if os.path.getsize(path) > limit:
            text += "\n\n[... {} kB more; open the file to see it all]".format((os.path.getsize(path) - limit) // 1024)
        self.text.get_buffer().set_text("{}\n{}\n\n{}".format(path, "=" * min(100, len(path)), text))
        self.notebook.set_current_page(1)

    # =====================================================================
    #  Display
    # =====================================================================
    def _show(self, result):
        self._showing = True
        try:
            self._show_result(result)
        finally:
            self._showing = False

    def _show_result(self, result):
        self.result = result
        self.combo_view.remove_all()
        for child in self.actions_box.get_children():
            self.actions_box.remove(child)
        for c in self.table.get_columns():
            self.table.remove_column(c)
        self.table.set_model(None)
        self._row_payload = []
        if result is None:
            self.plot.clear()
            self.text.get_buffer().set_text("")
            self.label_hint.set_text("")
            self.btn_export.set_sensitive(False); self.btn_save.set_sensitive(False)
            return
        views = result.get("views") or []
        for name, _setter in views:
            self.combo_view.append_text(name)
        # . views without a plot (the logs): the table and the text take the whole height
        self.plot.set_visible(bool(views))
        self.combo_view.set_sensitive(bool(views))
        if views:
            self.combo_view.set_active(0)
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
        self.notebook.set_current_page(0 if cols else 1)
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
        if r is None or info.get("index") is None:
            return
        if r.get("row_action") == "scan" and self.plot.kind == "lines" and getattr(self, "_scan_object", None) is not None:
            m = info["index"]
            if m < len(self._scan_order):
                self._goto_scan_frame(self._scan_order[m])
        elif r.get("click_atoms") and self._matches and self.plot.data.get("xlabels"):
            k = int(round(info.get("x", info["index"]))) if self.plot.kind == "lines" else info["index"]
            if 0 <= k < len(r["click_atoms"]):
                self._select_atoms(r["click_atoms"][k])

    def _set_views(self, views):
        self.result["views"] = views
        self.plot.set_visible(bool(views))
        self.combo_view.set_sensitive(bool(views))
        self.combo_view.handler_block_by_func(self._on_view_changed)
        self.combo_view.remove_all()
        for name, _setter in views:
            self.combo_view.append_text(name)
        self.combo_view.handler_unblock_by_func(self._on_view_changed)
        if views:
            self.combo_view.set_active(0)

    def _payload_of(self, row):
        for k, r in enumerate(self._table_rows()):
            if self._row_payload[k] is not None and [str(v) for v in r] == row:
                return self._row_payload[k]
        return None

    def _on_table_select(self, selection):
        """ Single click: the torsion views isolate the chosen bond. """
        r = self.result
        if r is None or getattr(self, "_showing", False) or r.get("row_action") not in ("library", "scans"):
            return
        model, it = selection.get_selected()
        if it is None:
            return
        payload = self._payload_of([model[it][k] for k in range(model.get_n_columns())])
        if not payload:
            return
        if payload[0] == "library":
            self._isolate_library(payload[1])
        elif payload[0] == "scan":
            s = (self.data.get("report") or {}).get("scan_details", [])[payload[1]]
            if s.get("profiles"):
                detail = self._view_scan(payload[1])
                self._set_views(detail["views"])
                self.text.get_buffer().set_text(detail["summary"])
            if self._matches and s.get("quad"):
                self._select_atoms(list(s["quad"]))

    def _on_table_row(self, _view, path, _column):
        r = self.result
        if r is None:
            return
        model = self.table.get_model()
        # . the table may be sorted: the payload follows the original row order
        row = [model[path][k] for k in range(model.get_n_columns())]
        payload = None
        for k, p in enumerate(self._row_payload):
            if p is not None and [str(v) for v in self._table_rows()[k]] == row:
                payload = p
                break
        if payload is None:
            return
        kind = r.get("row_action")
        if kind == "file":
            self._show_file(payload[1])
        elif kind == "scan":
            if getattr(self, "_scan_object", None) is not None:
                self._goto_scan_frame(payload[1])
        elif kind in ("library", "scans"):
            return
        elif self._matches and payload:
            self._select_atoms(list(payload))

    def _table_rows(self):
        cols, rows = self.result.get("table", ([], []))
        return [r[:len(cols)] for r in rows]

    # ------------------------------------------------------------- 3D
    def _select_atoms(self, indices):
        v = self._vobject()
        if v is None:
            return
        self.vm_session._selection_function_set(None)
        sel = self.vm_session.selections[self.vm_session.current_selection]
        mode = sel.selection_mode
        sel.selection_mode = "atom"
        self.vm_session._selection_function_set({v.atoms[i] for i in indices if i in v.atoms})
        sel.selection_mode = mode
        self.vm_session.vm_glcore.queue_draw()

    def _color_by_charge(self):
        v = self._vobject()
        if v is None:
            return
        for a in self.data["atoms"]:
            q = a.get("charge") or 0.0
            t = max(-1.0, min(1.0, q / 0.8))
            color = [1.0, 1.0 - t, 1.0 - t] if t > 0 else [1.0 + t, 1.0 + t, 1.0]
            if a["index"] in v.atoms:
                v.atoms[a["index"]].color = color
        self.vm_session.set_color_by_index(vobject=v, indexes=[], color=None)
        self.vm_session.vm_glcore.queue_draw()
        self._message("Atoms outside the library colored by charge: red positive, blue negative (±0.8 e saturates).")

    def _play_scan(self, s):
        """ The xTB scan geometries as a new object (a new system of the fragment, caps included). """
        from pdynamo.torsion_fit import read_xyz_frames
        from pBabel import ImportSystem
        logdir = self._scan_logdir(s)
        frames = [xyz for _comment, xyz in read_xyz_frames(os.path.join(logdir, "xtbscan.log"))]
        system = ImportSystem(os.path.join(logdir, "fragment.xyz"), log=None)
        system.BondsFromCoordinates3()
        j, k = (int(i) for i in s["bond"].split("-"))
        name = "scan {}-{}".format(self._short(j), self._short(k))
        system.label = name
        keep = self._current_id
        system = self.p_session.add_prepared_pdynamo_system(system, name=name, tag="xTB scan",
                                                            working_folder=getattr(self._system(), "e_working_folder", None))
        self._scan_ids.add(system.e_id)
        if self.combobox_systems.get_system_id() != keep:
            self._on_system_changed(None)
        v = self.main.get_active_vobject(system.e_id)
        v.frames = np.asarray(frames, dtype=np.float32)
        self._scan_object = v
        self.vm_session.set_frame(0)
        self.vm_session.vm_glcore.queue_draw()
        self._message("'{}': {} frames of the relaxed xTB scan (fragment with H caps). Click a point of the profile to go to "
                      "that geometry; the mini player plays the scan.".format(name, len(frames)))

    def _check_torsions(self, keys):
        """ xTB checks of library torsions (residue_library.check_torsion) in a worker thread;
            the results go into parametrization.json (report["torsion_checks"]) and report.html. """
        import threading
        import concurrent.futures
        from gi.repository import GLib
        from pdynamo import residue_library as RL
        if not keys:
            self._message("Nothing to check (all bonds already checked)."); return
        if getattr(self, "_checking", False):
            self._message("A check is already running."); return
        system = self._system()
        if system is None or not self._matches:
            raise ValueError("The selected system is not the parametrized one.")
        n = len(system.atoms)
        z = [a.atomicNumber for a in system.atoms]
        xyz = np.array([[system.coordinates3[i, k] for k in range(3)] for i in range(n)])
        formal = [int(getattr(a, "formalCharge", 0) or 0) for a in system.atoms]
        bonds = RL.system_bonds(system)
        edits = {kind: dict(self.data["edits"].get(kind) or {}) for kind in ("bond_terms", "angle_terms", "dihedral_bonds")}
        mult = int((self.data.get("settings") or {}).get("multiplicity_override") or 1)
        data = self.data
        root = os.path.join(data["_folder"], "logs", "torsions")
        self._checking = True
        state = {"done": 0}

        def one(key):
            j, k = (int(i) for i in key.split("-"))
            try:
                return key, RL.check_torsion(z, xyz, bonds, formal, (j, k), edits, multiplicity=mult, threads=1,
                                             log_dir=os.path.join(root, "check_{}-{}".format(j + 1, k + 1)))
            except Exception as error:
                traceback.print_exc()
                return key, {"bond": key, "error": str(error).split("\n")[0][:200]}
            finally:
                state["done"] += 1
                GLib.idle_add(lambda: self._message("xTB check of library torsions: {}/{} done...".format(state["done"], len(keys))) and False)

        def work():
            with concurrent.futures.ThreadPoolExecutor(max_workers=min(4, len(keys))) as pool:
                results = list(pool.map(one, keys))
            GLib.idle_add(finish, results)

        def finish(results):
            self._checking = False
            checks = data.setdefault("report", {}).setdefault("torsion_checks", {})
            for key, r in results:
                checks[key] = r
            self._save(data)
            bad = [r["error"] for _k, r in results if "error" in r]
            if self.data is data and self.result is not None and self.result.get("row_action") == "library":
                self._show(self._view_library())
                if len(results) == 1:
                    self._isolate_library(results[0][0])
            self._message("xTB check: {} bond(s) done{}; saved in parametrization.json and report.html.".format(
                len(results) - len(bad), ", {} failed ({})".format(len(bad), bad[0]) if bad else ""))
            return False

        self._message("xTB check of {} library torsion(s) running in the background...".format(len(keys)))
        threading.Thread(target=work, daemon=True).start()

    def _save(self, data):
        from pdynamo import param_report as PR
        clean = {k: v for k, v in data.items() if not k.startswith("_")}
        PR.write(data["_folder"], clean)

    def _goto_scan_frame(self, frame):
        self.vm_session.set_frame(int(frame))
        self.vm_session.vm_glcore.queue_draw()

    # ------------------------------------------------------------- export
    def on_export_image(self, _b):
        path = self._ask_save("Export Image", ["png", "svg"])
        if path:
            self.plot.export(path)

    def on_save_table(self, _b):
        path = self._ask_save("Save Table", ["csv"])
        if not path or self.result is None:
            return
        import csv
        cols, _rows = self.result.get("table", ([], []))
        with open(path, "w", newline="") as handle:
            w = csv.writer(handle)
            w.writerow(cols)
            w.writerows(self._table_rows())

    def _ask_save(self, title, extensions):
        chooser = Gtk.FileChooserDialog(title=title, parent=self.window, action=Gtk.FileChooserAction.SAVE)
        chooser.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        chooser.set_do_overwrite_confirmation(True)
        if self.data is not None:
            chooser.set_current_folder(self.data["_folder"])
        chooser.set_current_name("parametrization." + extensions[0])
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
            self._message("Error: {}".format(error))

    def _message(self, text):
        self.label_hint.set_text(text)

    def _on_destroy(self, *_):
        self.Visible = False
        self.window = None
