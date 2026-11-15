"""
esp_charges_window
==================

System > Charges > ESP / RESP Charges...

Fits ESP, RESP-1 and RESP-2 (AMBER two-stage) charges to the electrostatic
potential of a pDynamo HF/DFT wavefunction (engine in pdynamo/esp_charges.py),
for the selected molecules or the whole system, in the gas phase. The charges
can be applied to the MM model, shown as labels, or exported.

The SCF runs in a forked process, so the interface stays responsive and the
calculation can be cancelled.
"""

import csv
import json
import os
import subprocess
import sys
import tempfile
import time
import traceback
import types

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, GLib

from pdynamo import esp_charges as EC

# . treeview columns
COL_INDEX, COL_NAME, COL_RES, COL_ELEM, COL_MM, COL_ESP, COL_R1, COL_R2 = range(8)
CHARGE_SETS = (("RESP-2 (AMBER two-stage)", "resp2"), ("RESP-1", "resp1"), ("ESP", "esp"))

_ELEMENTS = {1: "H", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 11: "Na", 12: "Mg", 14: "Si",
             15: "P", 16: "S", 17: "Cl", 19: "K", 20: "Ca", 26: "Fe", 30: "Zn", 35: "Br", 53: "I"}
MAX_ATOMS_WARNING = 80


class ESPChargesWindow:

    def __init__(self, main=None):
        self.main        = main
        self.window      = None
        self.process     = None      # . subprocess.Popen of the running calculation
        self.workdir     = None
        self.result      = None
        self.job         = None      # . system, indices, labels of the last run
        self._t0         = 0.0

    # ------------------------------------------------------------------ helpers
    @property
    def p_session(self):
        return self.main.p_session

    @property
    def vm_session(self):
        return self.main.vm_session

    def _active_system(self):
        return self.p_session.psystem.get(self.p_session.active_id)

    def _vobject_of(self, system):
        objs = [v for v in self.vm_session.vm_objects_dic.values() if getattr(v, "e_id", None) == system.e_id]
        return objs[-1] if objs else None

    def _frame_coordinates(self, vobject, indices):
        """Coordinates (angstroms) of the atoms in the frame shown in the 3D view."""
        frame = min(int(self.vm_session.get_frame() or 0), len(vobject.frames) - 1)
        coords = vobject.frames[frame]
        return [[float(coords[i][0]), float(coords[i][1]), float(coords[i][2])] for i in indices]

    def _target(self):
        """(system, vobject, indices) of the atoms to fit, or raises ValueError."""
        if self.radio_selection.get_active():
            selection = self.vm_session.selections[self.vm_session.current_selection]
            atoms = list(selection.selected_atoms)
            if not atoms:
                raise ValueError("Nothing is selected. Select the molecule(s) to fit "
                                 "(selecting by molecule is the easiest way).")
            vobjects = {id(a.vm_object): a.vm_object for a in atoms}
            if len(vobjects) > 1:
                raise ValueError("The selection spans {} objects. Select atoms of one object.".format(len(vobjects)))
            vobject = atoms[0].vm_object
            system  = self.p_session.psystem.get(vobject.e_id)
            if system is None:
                raise ValueError("The selected object has no system.")
            indices = sorted({a.index - 1 for a in atoms})
        else:
            system = self._active_system()
            if system is None:
                raise ValueError("There is no active system.")
            vobject = self._vobject_of(system)
            indices = list(range(len(system.atoms)))
        if vobject is None or len(vobject.atoms) != len(system.atoms):
            raise ValueError("The 3D object does not match the system (different number of atoms).")
        return system, vobject, indices

    # ------------------------------------------------------------------ window
    def open_window(self):
        if self.window is not None:
            self.refresh_target()
            self.window.present()
            return

        self.window = Gtk.Window(title="ESP / RESP Charges")
        self.window.set_default_size(760, 640)
        self.window.set_border_width(10)
        self.window.set_transient_for(self.main.window)
        self.window.connect("destroy", self._on_destroy)

        vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
        self.window.add(vbox)

        # --- atoms
        frame = Gtk.Frame(label=" Atoms ")
        grid = Gtk.Grid(column_spacing=10, row_spacing=6, margin=8)
        frame.add(grid)
        self.radio_selection = Gtk.RadioButton.new_with_label(None, "Selected molecule(s)")
        self.radio_system    = Gtk.RadioButton.new_with_label_from_widget(self.radio_selection, "Whole active system")
        self.radio_selection.connect("toggled", lambda *_: self.refresh_target())
        refresh = Gtk.Button(label="Refresh")
        refresh.set_tooltip_text("Read the current selection again")
        refresh.connect("clicked", lambda *_: self.refresh_target())
        self.label_target = Gtk.Label(xalign=0)
        self.label_target.set_line_wrap(True)
        grid.attach(self.radio_selection, 0, 0, 1, 1)
        grid.attach(self.radio_system,    1, 0, 1, 1)
        grid.attach(refresh,              2, 0, 1, 1)
        grid.attach(self.label_target,    0, 1, 3, 1)
        vbox.pack_start(frame, False, False, 0)

        # --- calculation
        frame = Gtk.Frame(label=" Calculation ")
        grid = Gtk.Grid(column_spacing=10, row_spacing=6, margin=8)
        frame.add(grid)
        self.combo_level = Gtk.ComboBoxText()
        for level in EC.LEVELS:
            self.combo_level.append_text(level)
        self.combo_level.set_active(list(EC.LEVELS).index(EC.DEFAULT_LEVEL))
        self.combo_level.set_tooltip_text("HF/6-31G* is the level of the AMBER force fields (GAFF, ff14SB). "
                                          "GFN2-xTB is fast but approximate: about 0.05 e from HF/6-31G* RESP "
                                          "on average, up to 0.2 e for some atoms.")
        self.spin_charge = Gtk.SpinButton.new_with_range(-20, 20, 1)
        self.spin_mult   = Gtk.SpinButton.new_with_range(1, 10, 1)
        grid.attach(Gtk.Label(label="QC level", xalign=0),     0, 0, 1, 1)
        grid.attach(self.combo_level,                           1, 0, 1, 1)
        grid.attach(Gtk.Label(label="Charge", xalign=0),       2, 0, 1, 1)
        grid.attach(self.spin_charge,                           3, 0, 1, 1)
        grid.attach(Gtk.Label(label="Multiplicity", xalign=0), 4, 0, 1, 1)
        grid.attach(self.spin_mult,                             5, 0, 1, 1)
        note = Gtk.Label(xalign=0)
        note.set_markup("<small>Gas phase, pDynamo SCF (GFN2-xTB levels: xtb --esp, on xtb's own vdW grid). "
                        "Grid: Merz-Kollman surfaces at 1.4-2.0 x vdW radii. "
                        "RESP: restraint a = 0.0005 / 0.001, b = 0.1 on heavy atoms; "
                        "CH<sub>2</sub>/CH<sub>3</sub> refitted in stage 2 (Bayly et al., 1993).</small>")
        note.set_line_wrap(True)
        grid.attach(note, 0, 1, 6, 1)
        vbox.pack_start(frame, False, False, 0)

        # --- run
        hbox = Gtk.Box(spacing=8)
        self.btn_run = Gtk.Button(label="Run")
        self.btn_run.connect("clicked", self.on_run)
        self.btn_cancel = Gtk.Button(label="Cancel")
        self.btn_cancel.set_sensitive(False)
        self.btn_cancel.connect("clicked", self.on_cancel)
        self.progress = Gtk.ProgressBar(show_text=True)
        self.progress.set_text("")
        hbox.pack_start(self.btn_run, False, False, 0)
        hbox.pack_start(self.btn_cancel, False, False, 0)
        hbox.pack_start(self.progress, True, True, 0)
        vbox.pack_start(hbox, False, False, 0)

        # --- results
        self.store = Gtk.ListStore(int, str, str, str, str, float, float, float)
        self.treeview = Gtk.TreeView(model=self.store)
        for title, col in (("#", COL_INDEX), ("Atom", COL_NAME), ("Residue", COL_RES), ("Element", COL_ELEM),
                           ("Current MM", COL_MM)):
            renderer = Gtk.CellRendererText()
            column = Gtk.TreeViewColumn(title, renderer, text=col)
            column.set_sort_column_id(col)
            column.set_min_width(60)
            self.treeview.append_column(column)
        for title, col in (("ESP", COL_ESP), ("RESP-1", COL_R1), ("RESP-2", COL_R2)):
            renderer = Gtk.CellRendererText(xalign=1.0)
            column = Gtk.TreeViewColumn(title, renderer)
            column.set_min_width(90)
            column.set_cell_data_func(renderer, self._format_charge, col)
            column.set_sort_column_id(col)
            self.treeview.append_column(column)
        scrolled = Gtk.ScrolledWindow()
        scrolled.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)
        scrolled.add(self.treeview)
        vbox.pack_start(scrolled, True, True, 0)

        self.label_summary = Gtk.Label(xalign=0)
        self.label_summary.set_selectable(True)
        self.label_summary.set_line_wrap(True)
        vbox.pack_start(self.label_summary, False, False, 0)

        # --- actions
        hbox = Gtk.Box(spacing=8)
        hbox.pack_start(Gtk.Label(label="Charge set"), False, False, 0)
        self.combo_set = Gtk.ComboBoxText()
        for text, _key in CHARGE_SETS:
            self.combo_set.append_text(text)
        self.combo_set.set_active(0)
        hbox.pack_start(self.combo_set, False, False, 0)
        self.btn_apply  = Gtk.Button(label="Apply to MM Model")
        self.btn_apply.set_tooltip_text("Replace the MM charges of these atoms by the chosen set")
        self.btn_labels = Gtk.Button(label="Show as Labels")
        self.btn_export = Gtk.Button(label="Export...")
        self.btn_apply.connect("clicked", self.on_apply)
        self.btn_labels.connect("clicked", self.on_labels)
        self.btn_export.connect("clicked", self.on_export)
        for b in (self.btn_apply, self.btn_labels, self.btn_export):
            b.set_sensitive(False)
            hbox.pack_start(b, False, False, 0)
        close = Gtk.Button(label="Close")
        close.connect("clicked", lambda *_: self.window.destroy())
        hbox.pack_end(close, False, False, 0)
        vbox.pack_start(hbox, False, False, 0)

        self.window.show_all()
        self.refresh_target()

    def _format_charge(self, _column, cell, model, it, col):
        value = model[it][col]
        cell.set_property("text", "" if self.result is None else "{:+.4f}".format(value))

    def refresh_target(self):
        """Describes the atoms that would be fitted and suggests the charge."""
        if self.window is None:
            return
        try:
            system, _vobject, indices = self._target()
        except ValueError as error:
            self.label_target.set_text(str(error))
            self.btn_run.set_sensitive(self.process is None and self.radio_system.get_active())
            return
        n = len(indices)
        heavy = sum(1 for i in indices if system.atoms[i].atomicNumber > 1)
        text = "{}: {} atoms ({} heavy).".format(system.label, n, heavy)
        if n > MAX_ATOMS_WARNING:
            text += "  Large for an HF/DFT calculation: expect a long run."
        self.label_target.set_text(text)
        if self.process is None:
            self.spin_charge.set_value(EC.suggested_charge(system, indices))
            electrons = sum(system.atoms[i].atomicNumber for i in indices) - int(self.spin_charge.get_value())
            self.spin_mult.set_value(2 if electrons % 2 else 1)
            self.btn_run.set_sensitive(True)

    # ------------------------------------------------------------------ run
    def on_run(self, _button):
        if self.process is not None:
            return
        try:
            system, vobject, indices = self._target()
            Z, _xyz, bonds = EC.molecule_from_system(system, indices)
        except ValueError as error:
            self._message(str(error), Gtk.MessageType.ERROR)
            return
        xyz = self._frame_coordinates(vobject, indices)
        charge = int(self.spin_charge.get_value())
        mult   = int(self.spin_mult.get_value())
        if (sum(Z) - charge) % 2 == (mult % 2):
            self._message("Charge {} and multiplicity {} are not compatible with {} electrons."
                          .format(charge, mult, sum(Z) - charge), Gtk.MessageType.ERROR)
            return
        level = self.combo_level.get_active_text()

        labels = []
        for i in indices:
            atom = system.atoms[i]
            residue = getattr(atom, "parent", None)
            labels.append((i, atom.label, getattr(residue, "label", ""), _ELEMENTS.get(atom.atomicNumber, str(atom.atomicNumber))))
        self.job = {"system": system, "indices": indices, "labels": labels, "xyz": xyz,
                    "vobject": vobject, "level": level}

        # . a separate Python process (not a fork of this GTK/OpenGL process)
        self.workdir = tempfile.mkdtemp(prefix="easyhybrid_esp_")
        job_path = os.path.join(self.workdir, "job.json")
        with open(job_path, "w") as handle:
            json.dump({"Z": Z, "xyz": xyz, "bonds": [list(b) for b in bonds], "charge": charge,
                       "multiplicity": mult, "level": level}, handle)
        src = os.path.dirname(os.path.dirname(os.path.abspath(EC.__file__)))
        env = dict(os.environ)
        env["PYTHONPATH"] = src + os.pathsep + env.get("PYTHONPATH", "")
        self.result_path = os.path.join(self.workdir, "result.json")
        self.process = subprocess.Popen([sys.executable, "-m", "pdynamo.esp_charges", job_path, self.result_path],
                                        cwd=self.workdir, env=env,
                                        stdout=open(os.path.join(self.workdir, "pdynamo.log"), "w"),
                                        stderr=subprocess.STDOUT)
        self._t0 = time.time()
        self.result = None
        self.store.clear()
        self.label_summary.set_text("")
        for b in (self.btn_run, self.btn_apply, self.btn_labels, self.btn_export):
            b.set_sensitive(False)
        self.btn_cancel.set_sensitive(True)
        GLib.timeout_add(200, self._poll)

    def _poll(self):
        if self.process is None:
            return False
        if self.process.poll() is None:
            self.progress.pulse()
            self.progress.set_text("{}: SCF and fit... {:.0f} s".format(self.job["level"], time.time() - self._t0))
            return True
        result = None
        if os.path.exists(self.result_path):
            with open(self.result_path) as handle:
                result = json.load(handle)
        log = os.path.join(self.workdir, "pdynamo.log")
        self._finish()
        if result is None:
            tail = open(log).read()[-1500:] if os.path.exists(log) else ""
            self.progress.set_text("failed")
            self._message("The calculation stopped without a result.\n\n" + tail, Gtk.MessageType.ERROR)
        elif "error" in result:
            self.progress.set_text("failed")
            self._message("The calculation failed:\n\n" + result["error"], Gtk.MessageType.ERROR)
        else:
            self._show_result(result)
        return False

    def _finish(self):
        self.process = None
        self.btn_cancel.set_sensitive(False)
        self.btn_run.set_sensitive(True)

    def on_cancel(self, _button):
        if self.process is not None:
            self.process.terminate()
            try:
                self.process.wait(timeout=3)
            except Exception:
                self.process.kill()
            self._finish()
            self.progress.set_fraction(0.0)
            self.progress.set_text("cancelled")

    def _show_result(self, result):
        self.result = result
        system  = self.job["system"]
        mm      = getattr(system, "mmState", None)
        self.store.clear()
        for k, (i, name, res, elem) in enumerate(self.job["labels"]):
            current = "{:+.4f}".format(mm.charges[i]) if mm is not None else "-"
            self.store.append([i + 1, name, res, elem, current,
                               result["esp"]["charges"][k], result["resp1"]["charges"][k], result["resp2"]["charges"][k]])
        qc_dip = result.get("qc_dipole")
        lines = ["{} energy {:.2f} kJ/mol, {} grid points, {:.1f} s.".format(
                     result["level"], result["energy"], result["points"], result["seconds"])]
        stats = []
        for text, key in (("ESP", "esp"), ("RESP-1", "resp1"), ("RESP-2", "resp2")):
            stats.append("{} RRMS {:.4f}, dipole {:.2f} D".format(text, result[key]["rrms"], result[key]["dipole"]))
        lines.append(";  ".join(stats) + ("" if qc_dip is None else ";  QC dipole {:.2f} D".format(qc_dip)))
        groups = result.get("methyl_methylene") or {}
        if groups:
            lines.append("RESP-2 refitted {} CH2/CH3 group(s).".format(len(groups)))
        self.label_summary.set_text("\n".join(lines))
        self.progress.set_fraction(1.0)
        self.progress.set_text("done in {:.1f} s".format(result["seconds"]))
        self.btn_apply.set_sensitive(mm is not None)
        self.btn_labels.set_sensitive(True)
        self.btn_export.set_sensitive(True)

    # ------------------------------------------------------------------ actions
    def _chosen(self):
        text, key = CHARGE_SETS[self.combo_set.get_active()]
        return text, key, self.result[key]["charges"]

    def on_apply(self, _button):
        if not self.result:
            return
        system  = self.job["system"]
        indices = self.job["indices"]
        text, _key, charges = self._chosen()
        mm = system.mmState
        old = sum(mm.charges[i] for i in indices)
        new = sum(charges)
        if abs(old - new) > 1.0e-3:
            dialog = Gtk.MessageDialog(transient_for=self.window, modal=True,
                                       message_type=Gtk.MessageType.QUESTION, buttons=Gtk.ButtonsType.YES_NO,
                                       text="The net charge of these atoms changes")
            dialog.format_secondary_text("Current MM charges sum to {:+.4f}, the {} charges to {:+.4f}. "
                                         "Apply anyway?".format(old, text, new))
            answer = dialog.run()
            dialog.destroy()
            if answer != Gtk.ResponseType.YES:
                return
        EC.apply_to_mm_model(system, indices, charges)
        for row, q in zip(self.store, charges):
            row[COL_MM] = "{:+.4f}".format(q)
        total = sum(mm.charges)
        self.progress.set_text("{} charges applied to {} atoms (system total {:+.4f})".format(text, len(indices), total))
        try:
            self.main.refresh_widgets()
        except Exception:
            pass

    def on_labels(self, _button):
        """ Labels with the chosen charges on the fitted atoms (the selection is not touched). """
        if not self.result:
            return
        _text, _key, charges = self._chosen()
        vobject = self.job["vobject"]
        atoms = [vobject.atoms[i] for i in self.job["indices"]]
        target = types.SimpleNamespace(selected_atoms=set(atoms), selected_objects=[vobject])
        self.vm_session.show_or_hide(rep_type="labels", selection=target, show=True)
        for atom, q in zip(atoms, charges):
            atom.label_text = "{:.3f}".format(q)
        self.vm_session.vm_glcore.queue_draw()

    def on_export(self, _button):
        if not self.result:
            return
        dialog = Gtk.FileChooserDialog(title="Export ESP / RESP charges", parent=self.window,
                                       action=Gtk.FileChooserAction.SAVE)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        dialog.set_do_overwrite_confirmation(True)
        dialog.set_current_name("{}_resp.csv".format(self.job["system"].label or "system"))
        for name, pattern in (("CSV table (*.csv)", "*.csv"), ("Charges, one per line (*.chg)", "*.chg")):
            f = Gtk.FileFilter()
            f.set_name(name)
            f.add_pattern(pattern)
            dialog.add_filter(f)
        if dialog.run() == Gtk.ResponseType.OK:
            path = dialog.get_filename()
            try:
                self.export(path)
                self.progress.set_text("exported to {}".format(os.path.basename(path)))
            except Exception as error:
                self._message("Could not write the file: {}".format(error), Gtk.MessageType.ERROR)
        dialog.destroy()

    def export(self, path):
        r = self.result
        if path.endswith(".chg"):
            _text, _key, charges = self._chosen()
            with open(path, "w") as handle:
                for q in charges:
                    handle.write("{:10.6f}\n".format(q))
            return
        with open(path, "w", newline="") as handle:
            handle.write("# {} ESP/RESP charges, {}, charge {}, multiplicity {}\n".format(
                self.job["system"].label, r["level"], r["charge"], r["multiplicity"]))
            handle.write("# RRMS: ESP {:.4f}, RESP-1 {:.4f}, RESP-2 {:.4f}; {} grid points\n".format(
                r["esp"]["rrms"], r["resp1"]["rrms"], r["resp2"]["rrms"], r["points"]))
            writer = csv.writer(handle)
            writer.writerow(["index", "atom", "residue", "element", "esp", "resp1", "resp2"])
            for k, (i, name, res, elem) in enumerate(self.job["labels"]):
                writer.writerow([i + 1, name, res, elem] + ["{:.6f}".format(r[key]["charges"][k])
                                                            for key in ("esp", "resp1", "resp2")])

    # ------------------------------------------------------------------ misc
    def _message(self, text, kind=Gtk.MessageType.INFO):
        dialog = Gtk.MessageDialog(transient_for=self.window, modal=True, message_type=kind,
                                   buttons=Gtk.ButtonsType.OK, text=text)
        dialog.run()
        dialog.destroy()

    def _on_destroy(self, *_):
        self.on_cancel(None)
        self.window = None
