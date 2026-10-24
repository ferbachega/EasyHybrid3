"""
charge_block_inspector_window
=============================

System > Charges > Charge Block Inspector (2026-10-24): the MM charge of every
residue, its deviation from an integer and the correction, shown before
anything is changed, in the layout of the analysis windows.

  left    tolerance, correction method, filter; Preview is live; Apply, Undo
  right   deviation of each residue (bars; click to select it in the 3D view),
          and the tabs Residues (editable target), Atoms (before, after, change
          of the residue chosen) and Report (summary and every change, exportable)

The calculation is in util/charge_blocks.py. Applying is refused with a QC
region: its MM charges are already redistributed around the QC atoms.
"""

import time

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, GLib

import numpy as np

from gui.widgets.custom_widgets import SystemComboBox
from util import charge_blocks as CB
from util.analysis_plot import AnalysisPlot

FILTERS = [("all", "All residues"), ("needs", "Non-integer residues"), ("user", "Targets set by hand")]
MAX_BARS = 60


class ChargeBlockInspectorWindow:

    def __init__(self, main=None):
        self.main = main
        self.p_session = main.p_session
        self.vm_session = main.vm_session
        self.Visible = False
        self.window = None
        self.system = None
        self.blocks = []
        self.rows = []
        self.targets = {}            # . key -> integer target set by hand
        self.undo_stack = []         # . (system id, charges before, time)

    # ------------------------------------------------------------------ build
    def open_window(self):
        if self.Visible:
            self.window.present()
            return
        self.window = Gtk.Window(title="Charge Block Inspector")
        self.window.set_default_size(1080, 680)
        self.window.set_transient_for(self.main.window)
        self.window.connect("destroy", self._on_destroy)
        outer = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6, margin=8)
        self.window.add(outer)

        frame = Gtk.Frame(label=" System ")
        grid = Gtk.Grid(column_spacing=8, row_spacing=4, margin=8)
        frame.add(grid)
        self.combobox_systems = SystemComboBox(self.main)
        self.combobox_systems.connect("changed", lambda *_: self.load_system())
        self.label_system = Gtk.Label(xalign=0)
        self.label_system.set_hexpand(True)
        grid.attach(Gtk.Label(label="System", xalign=0), 0, 0, 1, 1)
        grid.attach(self.combobox_systems, 1, 0, 1, 1)
        grid.attach(self.label_system, 2, 0, 1, 1)
        outer.pack_start(frame, False, False, 0)

        paned = Gtk.Paned(orientation=Gtk.Orientation.HORIZONTAL)
        outer.pack_start(paned, True, True, 0)
        paned.pack1(self._build_left(), False, False)
        paned.pack2(self._build_right(), True, False)

        self.window.show_all()
        self.Visible = True
        if self.p_session.psystem:
            self.combobox_systems.set_active_system(e_id=self.p_session.active_id)
        self.load_system()

    def _build_left(self):
        left = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
        left.set_size_request(290, -1)
        frame = Gtk.Frame(label=" Options ")
        grid = Gtk.Grid(column_spacing=8, row_spacing=6, margin=8)
        frame.add(grid)
        self.spin_tol = Gtk.SpinButton.new_with_range(1e-6, 0.1, 1e-4)
        self.spin_tol.set_digits(6)
        self.spin_tol.set_value(1e-4)
        self.spin_tol.set_tooltip_text("Residues whose charge is within this distance of their target (e) are "
                                       "left as they are.")
        self.spin_tol.connect("value-changed", lambda *_: self.update())
        self.combo_method = Gtk.ComboBoxText()
        for key, text in CB.METHODS:
            self.combo_method.append(key, text)
        self.combo_method.set_active_id("even")
        self.combo_method.set_tooltip_text("Same shift: the difference is divided evenly over the atoms of the "
                                           "residue (the classic correction).\nProportional to |q|: each atom takes "
                                           "a share proportional to its charge; atoms without charge keep zero.")
        self.combo_method.connect("changed", lambda *_: self.update())
        self.combo_filter = Gtk.ComboBoxText()
        for key, text in FILTERS:
            self.combo_filter.append(key, text)
        self.combo_filter.set_active_id("needs")
        self.combo_filter.connect("changed", lambda *_: self._fill_residues())
        rows = [("Tolerance (e)", self.spin_tol), ("Correction", self.combo_method), ("Show", self.combo_filter)]
        for r, (text, w) in enumerate(rows):
            grid.attach(Gtk.Label(label=text, xalign=0), 0, r, 1, 1)
            grid.attach(w, 1, r, 1, 1)
        left.pack_start(frame, False, False, 0)

        box = Gtk.Box(spacing=6, homogeneous=True)
        self.btn_apply = Gtk.Button(label="Apply")
        self.btn_apply.set_tooltip_text("Write the corrected charges into the MM model of the system.")
        self.btn_apply.connect("clicked", lambda *_: self.apply())
        self.btn_undo = Gtk.Button(label="Undo")
        self.btn_undo.set_tooltip_text("Back to the charges before the last Apply.")
        self.btn_undo.connect("clicked", lambda *_: self.undo())
        box.pack_start(self.btn_apply, True, True, 0)
        box.pack_start(self.btn_undo, True, True, 0)
        left.pack_start(box, False, False, 0)
        b = Gtk.Button(label="Reset Targets")
        b.set_tooltip_text("Back to the nearest integer as the target of every residue.")
        b.connect("clicked", lambda *_: (self.targets.clear(), self.update()))
        left.pack_start(b, False, False, 0)

        self.label_state = Gtk.Label(xalign=0)
        self.label_state.set_line_wrap(True)
        self.label_state.set_max_width_chars(38)
        left.pack_start(self.label_state, False, False, 0)
        hint = Gtk.Label(xalign=0)
        hint.set_line_wrap(True)
        hint.set_max_width_chars(38)
        hint.set_markup("<small>Each residue is brought to an integer charge (the nearest one, or the target "
                        "typed in the table) by shifting the charges of its atoms. Nothing changes until "
                        "Apply. Double-click a residue to select it in the 3D view.</small>")
        left.pack_end(hint, False, False, 0)
        return left

    def _build_right(self):
        right = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6, margin_start=6)
        vpaned = Gtk.Paned(orientation=Gtk.Orientation.VERTICAL)
        self.plot = AnalysisPlot()
        self.plot.on_click = self._on_plot_click
        vpaned.pack1(self.plot, True, False)

        nb = Gtk.Notebook()
        self.notebook = nb
        # . Residues
        self.res_store = Gtk.ListStore(int, str, str, int, str, str, str, str, str)
        self.res_view = Gtk.TreeView(model=self.res_store)
        titles = ("Residue", "Chain", "Atoms", "Charge (e)", "Target", "Deviation", "Shift per atom", "Status")
        for k, title in enumerate(titles):
            renderer = Gtk.CellRendererText()
            if title == "Target":
                renderer.set_property("editable", True)
                renderer.set_property("weight", 600)
                renderer.connect("edited", self._on_target_edited)
            col = Gtk.TreeViewColumn(title, renderer, text=k + 1)
            col.set_resizable(True)
            self.res_view.append_column(col)
        self.res_view.get_selection().connect("changed", self._on_residue_selected)
        self.res_view.connect("row-activated", lambda v, p, c: self.select_residues([self.res_store[p][0]]))
        sw = Gtk.ScrolledWindow(); sw.add(self.res_view)
        nb.append_page(sw, Gtk.Label(label="Residues"))
        # . Atoms
        self.atom_store = Gtk.ListStore(int, str, str, str, str)
        self.atom_view = Gtk.TreeView(model=self.atom_store)
        for k, title in enumerate(("Atom", "Name", "Before (e)", "After (e)", "Change (e)")):
            col = Gtk.TreeViewColumn(title, Gtk.CellRendererText(), text=k)
            col.set_resizable(True)
            self.atom_view.append_column(col)
        sw = Gtk.ScrolledWindow(); sw.add(self.atom_view)
        nb.append_page(sw, Gtk.Label(label="Atoms"))
        # . Report
        self.report = Gtk.TextView(editable=False, monospace=True, wrap_mode=Gtk.WrapMode.NONE)
        sw = Gtk.ScrolledWindow(); sw.add(self.report)
        nb.append_page(sw, Gtk.Label(label="Report"))
        nb.set_size_request(-1, 250)
        vpaned.pack2(nb, False, False)
        vpaned.set_position(330)
        right.pack_start(vpaned, True, True, 0)

        actions = Gtk.Box(spacing=6)
        b = Gtk.Button(label="Select Non-integer Residues")
        b.set_tooltip_text("Select in the 3D view every residue that the correction changes.")
        b.connect("clicked", lambda *_: self.select_residues([r["key"] for r in self.rows if r["needs"]]))
        actions.pack_start(b, False, False, 0)
        b = Gtk.Button(label="Export Report...")
        b.set_tooltip_text("The report and every change of charge, atom by atom, as text.")
        b.connect("clicked", lambda *_: self.on_export())
        actions.pack_end(b, False, False, 0)
        right.pack_start(actions, False, False, 0)
        return right

    # ------------------------------------------------------------------ data
    def load_system(self):
        if not self.Visible:
            return
        sid = self.combobox_systems.get_system_id()
        self.system = self.p_session.psystem.get(sid) if sid is not None else None
        self.targets = {}
        if self.system is None or self.system.mmModel is None or self.system.mmState is None:
            self.blocks, self.rows = [], []
            self.label_system.set_text("The system has no MM charges (no force field).")
            self._clear()
            return
        self.blocks = CB.residue_blocks(self.system)
        self.update()

    def _charges(self):
        return np.array(list(self.system.mmState.charges), dtype=float)

    def update(self):
        if not self.Visible or self.system is None or not self.blocks:
            return
        q = self._charges()
        self.rows = CB.analyze(q, self.blocks, self.spin_tol.get_value(), self.targets)
        method = self.combo_method.get_active_id() or "even"
        self.new = CB.correct(q, self.rows, method)
        self.sum = CB.summary(q, self.new, self.rows)
        qc = self.system.qcModel is not None
        self.label_system.set_markup(
            "{} atoms, {} residues · total MM charge <b>{:+.6f}</b> e{}".format(
                len(q), len(self.blocks), q.sum(),
                " · <span foreground='#c0392b'>QC region: Apply disabled</span>" if qc else ""))
        n = self.sum["n_changed"]
        if qc:
            state = ("This system has a QC region: its MM charges are already redistributed around the QC atoms. "
                     "Remove the QC region (System > Remove > QC Atoms) to apply the correction.")
        elif n == 0:
            state = "Every residue has an integer charge (within the tolerance): nothing to do."
        else:
            state = "{} residue(s) will change; total charge {:+.6f} -> {:+.6f} e; largest change on one atom {:.5f} e.".format(
                n, self.sum["total_before"], self.sum["total_after"], self.sum["max_atom_shift"])
        self.label_state.set_text(state)
        self.btn_apply.set_sensitive(n > 0 and not qc)
        self.btn_undo.set_sensitive(any(u[0] is self.system for u in self.undo_stack))
        self._fill_plot()
        self._fill_residues()
        self._fill_report()

    def _clear(self):
        self.plot.clear()
        self.res_store.clear()
        self.atom_store.clear()
        self.report.get_buffer().set_text("")
        self.btn_apply.set_sensitive(False)

    def _fill_plot(self):
        shown = sorted([r for r in self.rows if r["needs"]], key=lambda r: -abs(r["deviation"]))[:MAX_BARS]
        if not shown:
            # . the residue charges as they are (all integer): their charge per residue
            shown = [r for r in self.rows if abs(r["charge"]) > 0.5][:MAX_BARS]
            self._bar_keys = [r["key"] for r in shown]
            if not shown:
                self.plot.clear()
                return
            self.plot.set_lines([{"x": list(range(len(shown))), "y": [r["charge"] for r in shown], "style": "bars",
                                  "colors": [(0.84, 0.15, 0.16) if r["charge"] > 0 else (0.12, 0.47, 0.71) for r in shown]}],
                                "residue", "charge (e)", "Charged residues (all residues are integer)",
                                xlabels=[self._res_label(r) for r in shown])
            return
        shown = sorted(shown, key=lambda r: r["key"])
        self._bar_keys = [r["key"] for r in shown]
        title = "Deviation from the target charge ({} of {} residues{})".format(
            len(shown), len(self.rows), ", the largest" if self.sum["n_changed"] > MAX_BARS else "")
        self.plot.set_lines([{"x": list(range(len(shown))), "y": [r["deviation"] for r in shown], "style": "bars",
                              "colors": [(0.84, 0.15, 0.16) if r["deviation"] > 0 else (0.12, 0.47, 0.71) for r in shown]}],
                            "residue", "charge - target (e)", title, hline=0.0,
                            xlabels=[self._res_label(r) for r in shown])

    def _res_label(self, r):
        return "{}:{}".format(r["chain"], r["label"]) if r["chain"] else r["label"]

    def _fill_residues(self):
        if not self.Visible:
            return
        key = self.combo_filter.get_active_id() or "needs"
        self.res_store.clear()
        for r in self.rows:
            if key == "needs" and not r["needs"]:
                continue
            if key == "user" and not r["user_target"]:
                continue
            status = "changes" if r["needs"] else "ok"
            if r["user_target"]:
                status += " (target set)"
            shift = "" if not r["needs"] else (
                "{:+.6f}".format(-r["deviation"] / r["n"]) if (self.combo_method.get_active_id() == "even")
                else "proportional")
            self.res_store.append([r["key"], r["label"], r["chain"], r["n"], "{:+.6f}".format(r["charge"]),
                                   "{:+g}".format(r["target"]) if r["target"] else "0", "{:+.6f}".format(r["deviation"]),
                                   shift, status])

    def _fill_report(self):
        s = self.sum
        method = dict(CB.METHODS)[self.combo_method.get_active_id() or "even"]
        lines = ["Charge Block Inspector - {}".format(getattr(self.system, "label", "")),
                 time.strftime("%Y-%m-%d %H:%M"), "",
                 "Residues                 {}".format(s["n_blocks"]),
                 "Residues to correct      {}".format(s["n_changed"]),
                 "Atoms changed            {}".format(s["n_atoms_changed"]),
                 "Total charge before      {:+.6f} e".format(s["total_before"]),
                 "Total charge after       {:+.6f} e".format(s["total_after"]),
                 "Largest change (atom)    {:.6f} e".format(s["max_atom_shift"]),
                 "Tolerance                {:g} e".format(self.spin_tol.get_value()),
                 "Correction               {}".format(method)]
        if s["worst"] is not None:
            w = s["worst"]
            lines.append("Largest deviation        {} {:+.6f} e".format(self._res_label(w), w["deviation"]))
        if s["user_targets"]:
            lines += ["", "WARNING: targets set by hand away from the nearest integer (check them):"]
            lines += ["  {} charge {:+.4f} -> target {:+g}".format(self._res_label(r), r["charge"], r["target"])
                      for r in s["user_targets"]]
        if self.system.qcModel is not None:
            lines += ["", "The system has a QC region: the correction can be inspected but not applied."]
        changed = [r for r in self.rows if r["needs"]]
        if changed:
            lines += ["", "{:<16s} {:>6s} {:>12s} {:>8s} {:>12s}".format("Residue", "Atoms", "Charge", "Target", "Deviation")]
            for r in changed:
                lines.append("{:<16s} {:>6d} {:>+12.6f} {:>+8g} {:>+12.6f}".format(
                    self._res_label(r)[:16], r["n"], r["charge"], r["target"], r["deviation"]))
        self.report.get_buffer().set_text("\n".join(lines))

    def _report_text(self):
        b = self.report.get_buffer()
        text = b.get_text(b.get_start_iter(), b.get_end_iter(), False)
        q = self._charges()
        out = [text, "", "Atom changes", "{:>7s} {:<16s} {:<6s} {:>11s} {:>11s} {:>11s}".format(
            "Atom", "Residue", "Name", "Before", "After", "Change")]
        for r in self.rows:
            if not r["needs"]:
                continue
            for i in r["indices"]:
                out.append("{:>7d} {:<16s} {:<6s} {:>+11.6f} {:>+11.6f} {:>+11.6f}".format(
                    i + 1, self._res_label(r)[:16], self.system.atoms[i].label, q[i], self.new[i], self.new[i] - q[i]))
        return "\n".join(out) + "\n"

    # ------------------------------------------------------------------ events
    def _on_residue_selected(self, sel):
        model, it = sel.get_selected()
        self.atom_store.clear()
        if it is None:
            return
        key = model[it][0]
        r = self.rows[key]
        q = self._charges()
        for i in r["indices"]:
            self.atom_store.append([i + 1, self.system.atoms[i].label, "{:+.6f}".format(q[i]),
                                    "{:+.6f}".format(self.new[i]), "{:+.6f}".format(self.new[i] - q[i])])

    def _on_target_edited(self, renderer, path, text):
        key = self.res_store[path][0]
        text = text.strip()
        if not text:
            self.targets.pop(key, None)
        else:
            try:
                value = float(text)
            except ValueError:
                return
            if abs(value - round(value)) > 1e-9:
                self.label_state.set_text("The target of a residue must be an integer charge.")
                return
            self.targets[key] = float(round(value))
        self.update()

    def _on_plot_click(self, info):
        k = info.get("index") if isinstance(info, dict) else None
        keys = getattr(self, "_bar_keys", [])
        if k is not None and 0 <= k < len(keys):
            self.select_residues([keys[k]])
            for row in self.res_store:
                if row[0] == keys[k]:
                    self.res_view.get_selection().select_iter(row.iter)
                    self.res_view.scroll_to_cell(row.path, None, True, 0.5, 0)
                    self.notebook.set_current_page(0)

    def select_residues(self, keys):
        """ The atoms of these residues as the selection of the 3D view (last object of the system). """
        if self.system is None:
            return
        objs = [v for v in self.vm_session.vm_objects_dic.values()
                if getattr(v, "e_id", None) == self.system.e_id and not getattr(v, "is_surface", False)]
        if not objs:
            return
        v = objs[-1]
        indices = [i for key in keys for i in self.rows[key]["indices"]]
        self.vm_session._selection_function_set(None)
        sel = self.vm_session.selections[self.vm_session.current_selection]
        mode = sel.selection_mode
        sel.selection_mode = "atom"
        self.vm_session._selection_function_set({v.atoms[i] for i in indices if i in v.atoms})
        sel.selection_mode = mode
        self.vm_session.vm_glcore.queue_draw()

    # ------------------------------------------------------------------ apply
    def apply(self):
        if self.system is None or self.system.qcModel is not None:
            return
        before = self._charges()
        for i, value in enumerate(self.new):
            self.system.mmState.charges[i] = float(value)
        # . the reference charges used when a QC region is defined later
        self.system.e_charges_backup = [float(x) for x in self.new]
        self.undo_stack.append((self.system, before, time.time()))
        n = self.sum["n_changed"]
        msg = "Charge Block Inspector: {} residue(s) corrected, total charge {:+.6f} -> {:+.6f} e.".format(
            n, self.sum["total_before"], self.sum["total_after"])
        self._status(msg)
        self.targets = {}
        self.update()
        self.label_state.set_text(msg + " Undo brings the previous charges back.")

    def undo(self):
        for k in range(len(self.undo_stack) - 1, -1, -1):
            system, before, _t = self.undo_stack[k]
            if system is self.system:
                for i, value in enumerate(before):
                    system.mmState.charges[i] = float(value)
                system.e_charges_backup = [float(x) for x in before]
                del self.undo_stack[k]
                self._status("Charge Block Inspector: previous charges restored.")
                self.update()
                return

    def _status(self, msg):
        try:
            self.main.bottom_notebook.status_teeview_add_new_item(message=msg, system=self.system)
        except Exception:
            pass

    def on_export(self):
        if self.system is None:
            return
        dialog = Gtk.FileChooserDialog(title="Export Report", parent=self.window, action=Gtk.FileChooserAction.SAVE)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        dialog.set_do_overwrite_confirmation(True)
        dialog.set_current_name("charge_blocks_{}.txt".format(str(getattr(self.system, "label", "system")).replace(" ", "_")))
        path = dialog.get_filename() if dialog.run() == Gtk.ResponseType.OK else None
        dialog.destroy()
        if path:
            with open(path, "w") as handle:
                handle.write(self._report_text())

    def _on_destroy(self, *_):
        self.Visible = False
        self.window = None
