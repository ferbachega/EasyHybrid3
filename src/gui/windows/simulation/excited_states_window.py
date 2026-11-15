"""
excited_states_window
=====================

Simulate > Excited States (MNDO-CI / MOM)...: vertical excited states of
the selected molecule(s) or the active system, in the gas phase, at the
geometry shown (engine: pdynamo/excited_states.py, run in a separate
Python process).
"""

import json
import os
import subprocess
import sys
import tempfile
import time

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, GLib

from pdynamo import esp_charges as EC
from pdynamo import excited_states as XS
from util.analysis_plot import AnalysisPlot, PALETTE


class ExcitedStatesWindow:

    def __init__(self, main):
        self.main = main
        self.p_session = main.p_session
        self.vm_session = main.vm_session
        self.window = None
        self.process = None
        self.result = None

    # ------------------------------------------------------------------ window
    def open_window(self):
        if self.window is not None:
            self.window.present()
            return
        self.window = Gtk.Window(title="Excited States (MNDO-CI / MOM)")
        self.window.set_default_size(820, 700)
        self.window.set_transient_for(self.main.window)
        self.window.connect("destroy", self._on_destroy)
        box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8, margin=10)
        self.window.add(box)

        grid = Gtk.Grid(column_spacing=8, row_spacing=6)
        self.radio_sel = Gtk.RadioButton.new_with_label(None, "Selected molecule(s)")
        self.radio_sys = Gtk.RadioButton.new_with_label_from_widget(self.radio_sel, "Whole active system")
        self.radio_sel.connect("toggled", lambda *_: self.refresh_target())
        self.label_target = Gtk.Label(xalign=0)
        self.spin_charge = Gtk.SpinButton.new_with_range(-20, 20, 1)
        self.spin_mult = Gtk.SpinButton.new_with_range(1, 10, 1)
        grid.attach(Gtk.Label(label="Atoms", xalign=0), 0, 0, 1, 1)
        grid.attach(self.radio_sel, 1, 0, 1, 1)
        grid.attach(self.radio_sys, 2, 0, 1, 1)
        refresh = Gtk.Button(label="Refresh")
        refresh.connect("clicked", lambda *_: self.refresh_target())
        grid.attach(refresh, 3, 0, 1, 1)
        grid.attach(self.label_target, 0, 1, 4, 1)
        grid.attach(Gtk.Label(label="Charge", xalign=0), 0, 2, 1, 1)
        grid.attach(self.spin_charge, 1, 2, 1, 1)
        grid.attach(Gtk.Label(label="Ground-state multiplicity", xalign=0), 2, 2, 1, 1)
        grid.attach(self.spin_mult, 3, 2, 1, 1)
        box.pack_start(grid, False, False, 0)

        self.notebook = Gtk.Notebook()
        box.pack_start(self.notebook, False, False, 0)

        # . CI page
        g = Gtk.Grid(column_spacing=8, row_spacing=6, margin=8)
        self.combo_ham = Gtk.ComboBoxText()
        for h in XS.SEMIEMPIRICAL:
            self.combo_ham.append_text(h.upper())
        self.combo_ham.set_active(0)
        self.combo_ci = Gtk.ComboBoxText()
        for m in ("Singles (CIS)", "Doubles", "Singles + doubles (CISD)", "Full CI"):
            self.combo_ci.append_text(m)
        self.combo_ci.set_active(2)
        self.spin_el = Gtk.SpinButton.new_with_range(2, 30, 1)
        self.spin_el.set_value(6)
        self.spin_orb = Gtk.SpinButton.new_with_range(2, 30, 1)
        self.spin_orb.set_value(6)
        self.spin_states = Gtk.SpinButton.new_with_range(2, 200, 1)
        self.spin_states.set_value(20)
        for k, (t, w) in enumerate((("Hamiltonian", self.combo_ham), ("CI", self.combo_ci), ("Active electrons", self.spin_el),
                                    ("Active orbitals", self.spin_orb), ("States", self.spin_states))):
            g.attach(Gtk.Label(label=t, xalign=0), 2 * (k % 3), k // 3, 1, 1)
            g.attach(w, 2 * (k % 3) + 1, k // 3, 1, 1)
        note = Gtk.Label(xalign=0)
        note.set_line_wrap(True)
        note.set_markup("<small>pDynamo QCModelMNDOCI on a reference with fractional HOMO/LUMO occupancy. "
                        "The active space is centred on the HOMO/LUMO; energies only (no intensities). "
                        "Full CI grows combinatorially: keep the active space small.</small>")
        g.attach(note, 0, 2, 6, 1)
        self.notebook.append_page(g, Gtk.Label(label="MNDO-CI"))

        # . MOM page
        g = Gtk.Grid(column_spacing=8, row_spacing=6, margin=8)
        self.combo_mom = Gtk.ComboBoxText()
        for h in XS.SEMIEMPIRICAL:
            self.combo_mom.append_text(h.upper())
        for level in EC.LEVELS:
            self.combo_mom.append_text(level)
        self.combo_mom.set_active(0)
        self.entry_swaps = Gtk.Entry(text="0,0; 0,1; 1,0; 1,1")
        self.entry_swaps.set_tooltip_text("Orbital swaps o,v = HOMO-o -> LUMO+v, separated by ';'")
        g.attach(Gtk.Label(label="Method", xalign=0), 0, 0, 1, 1)
        g.attach(self.combo_mom, 1, 0, 1, 1)
        g.attach(Gtk.Label(label="Excitations (o,v)", xalign=0), 2, 0, 1, 1)
        g.attach(self.entry_swaps, 3, 0, 1, 1)
        note = Gtk.Label(xalign=0)
        note.set_line_wrap(True)
        note.set_markup("<small>Maximum overlap method (Gilbert, Besley &amp; Gill 2008): one SCF per excited "
                        "determinant (HOMO-o -> LUMO+v, alpha spin, unrestricted), transition dipole and "
                        "oscillator strength from pDynamo TransitionDipoleMoment. Open-shell singlets are "
                        "spin-mixed determinants.</small>")
        g.attach(note, 0, 1, 4, 1)
        self.notebook.append_page(g, Gtk.Label(label="MOM (with intensities)"))

        hb = Gtk.Box(spacing=8)
        self.btn_run = Gtk.Button(label="Run")
        self.btn_run.connect("clicked", self.on_run)
        self.btn_cancel = Gtk.Button(label="Cancel")
        self.btn_cancel.set_sensitive(False)
        self.btn_cancel.connect("clicked", self.on_cancel)
        self.progress = Gtk.ProgressBar(show_text=True)
        self.progress.set_text("")
        hb.pack_start(self.btn_run, False, False, 0)
        hb.pack_start(self.btn_cancel, False, False, 0)
        hb.pack_start(self.progress, True, True, 0)
        box.pack_start(hb, False, False, 0)

        paned = Gtk.Paned(orientation=Gtk.Orientation.VERTICAL)
        self.store = Gtk.ListStore(str, str, str, str, str, str)
        self.tree = Gtk.TreeView(model=self.store)
        for k, t in enumerate(("state", "dE (eV)", "dE (kJ/mol)", "lambda (nm)", "S / f", "note")):
            self.tree.append_column(Gtk.TreeViewColumn(t, Gtk.CellRendererText(), text=k))
        sw = Gtk.ScrolledWindow()
        sw.add(self.tree)
        paned.pack1(sw, True, False)
        self.plot = AnalysisPlot()
        self.plot.set_size_request(-1, 240)
        paned.pack2(self.plot, True, False)
        box.pack_start(paned, True, True, 0)

        hb = Gtk.Box(spacing=8)
        self.btn_export = Gtk.Button(label="Export Table...")
        self.btn_export.set_sensitive(False)
        self.btn_export.connect("clicked", self.on_export)
        close = Gtk.Button(label="Close")
        close.connect("clicked", lambda *_: self.window.destroy())
        hb.pack_start(self.btn_export, False, False, 0)
        hb.pack_end(close, False, False, 0)
        box.pack_start(hb, False, False, 0)
        self.window.show_all()
        self.refresh_target()

    # ------------------------------------------------------------------ target
    def _target(self):
        if self.radio_sel.get_active():
            sel = self.vm_session.selections[self.vm_session.current_selection]
            atoms = list(sel.selected_atoms)
            if not atoms:
                raise ValueError("Nothing is selected (select the molecule, e.g. selecting by molecule).")
            if len({id(a.vm_object) for a in atoms}) > 1:
                raise ValueError("The selection spans several objects.")
            v = atoms[0].vm_object
            system = self.p_session.psystem.get(v.e_id)
            idx = sorted({a.index - 1 for a in atoms})
        else:
            system = self.p_session.psystem.get(self.p_session.active_id)
            if system is None:
                raise ValueError("There is no active system.")
            objs = [o for o in self.vm_session.vm_objects_dic.values() if getattr(o, "e_id", None) == system.e_id]
            v = objs[-1] if objs else None
            idx = list(range(len(system.atoms)))
        if v is None or len(v.atoms) != len(system.atoms):
            raise ValueError("The 3D object does not match the system.")
        return system, v, idx

    def refresh_target(self):
        try:
            system, _v, idx = self._target()
        except ValueError as error:
            self.label_target.set_text(str(error))
            return
        self.label_target.set_text("{}: {} atoms".format(system.label, len(idx)))
        if self.process is None:
            q = EC.suggested_charge(system, idx)
            self.spin_charge.set_value(q)
            electrons = sum(system.atoms[i].atomicNumber for i in idx) - q
            self.spin_mult.set_value(2 if electrons % 2 else 1)

    # ------------------------------------------------------------------ run
    def on_run(self, _b):
        if self.process is not None:
            return
        try:
            system, v, idx = self._target()
            Z, _xyz, _bonds = EC.molecule_from_system(system, idx)
        except ValueError as error:
            self._message(str(error))
            return
        frame = min(int(self.vm_session.get_frame() or 0), len(v.frames) - 1)
        xyz = [[float(c) for c in v.frames[frame][i]] for i in idx]
        job = {"Z": Z, "xyz": xyz, "charge": int(self.spin_charge.get_value()), "multiplicity": int(self.spin_mult.get_value())}
        if self.notebook.get_current_page() == 0:
            job.update({"kind": "ci", "hamiltonian": XS.SEMIEMPIRICAL[self.combo_ham.get_active()],
                        "method": XS.CI_METHODS[self.combo_ci.get_active()], "electrons": int(self.spin_el.get_value()),
                        "orbitals": int(self.spin_orb.get_value()), "states": int(self.spin_states.get_value())})
        else:
            try:
                swaps = [[int(a) for a in part.split(",")] for part in self.entry_swaps.get_text().split(";") if part.strip()]
                assert swaps and all(len(s) == 2 and s[0] >= 0 and s[1] >= 0 for s in swaps)
            except Exception:
                self._message("Write the excitations as o,v pairs separated by ';' (for example 0,0; 1,0).")
                return
            text = self.combo_mom.get_active_text()
            job.update({"kind": "mom", "swaps": swaps, "method": text.lower() if text.lower() in XS.SEMIEMPIRICAL else text})
        self.job = job
        self.workdir = tempfile.mkdtemp(prefix="easyhybrid_xs_")
        job_path = os.path.join(self.workdir, "job.json")
        with open(job_path, "w") as handle:
            json.dump(job, handle)
        self.result_path = os.path.join(self.workdir, "result.json")
        src = os.path.dirname(os.path.dirname(os.path.abspath(XS.__file__)))
        env = dict(os.environ)
        env["PYTHONPATH"] = src + os.pathsep + env.get("PYTHONPATH", "")
        self.process = subprocess.Popen([sys.executable, "-m", "pdynamo.excited_states", job_path, self.result_path],
                                        cwd=self.workdir, env=env, stdout=open(os.path.join(self.workdir, "pdynamo.log"), "w"),
                                        stderr=subprocess.STDOUT)
        self._t0 = time.time()
        self.btn_run.set_sensitive(False)
        self.btn_cancel.set_sensitive(True)
        GLib.timeout_add(200, self._poll)

    def _poll(self):
        if self.process is None:
            return False
        if self.process.poll() is None:
            self.progress.pulse()
            self.progress.set_text("running... {:.0f} s".format(time.time() - self._t0))
            return True
        self.process = None
        self.btn_run.set_sensitive(True)
        self.btn_cancel.set_sensitive(False)
        result = None
        if os.path.exists(self.result_path):
            with open(self.result_path) as handle:
                result = json.load(handle)
        if result is None or "error" in result:
            self.progress.set_text("failed")
            log = os.path.join(self.workdir, "pdynamo.log")
            detail = result["error"] if result else (open(log).read()[-1500:] if os.path.exists(log) else "")
            self._message("The calculation failed:\n\n" + detail)
            return False
        self.progress.set_fraction(1.0)
        self.progress.set_text("done in {:.1f} s".format(result.get("seconds", 0)))
        self.show_result(result)
        return False

    def on_cancel(self, _b):
        if self.process is not None:
            self.process.terminate()
            self.process = None
            self.btn_run.set_sensitive(True)
            self.btn_cancel.set_sensitive(False)
            self.progress.set_text("cancelled")

    # ------------------------------------------------------------------ result
    def show_result(self, r):
        self.result = r
        self.store.clear()
        lines = []
        if "states" in r:
            for k, s in enumerate(r["states"]):
                nm = "" if not s["nm"] else "{:.0f}".format(s["nm"])
                self.store.append([str(k), "{:.3f}".format(s["energy_ev"]), "{:.1f}".format(s["energy_kj"]), nm,
                                   "S = {:.2f} ({})".format(s["S"], {1: "singlet", 2: "doublet", 3: "triplet", 4: "quartet",
                                                                    5: "quintet"}.get(s["multiplicity"], s["multiplicity"])),
                                   "lowest" if k == 0 else ""])
            ex = [s for s in r["states"][1:]]
            series = []
            for mult, color in ((1, PALETTE[0]), (3, PALETTE[3]), (2, PALETTE[2])):
                xs = [s["energy_ev"] for s in ex if s["multiplicity"] == mult]
                if xs:
                    series.append({"x": xs, "y": [1.0] * len(xs), "style": "bars", "color": color,
                                   "label": {1: "singlets", 2: "doublets", 3: "triplets"}[mult]})
            if series:
                self.plot.set_lines(series, "excitation energy (eV)", "state", "{} CI ({} e, {} orbitals): states".format(
                    r["hamiltonian"].upper(), r["active"][0], r["active"][1]), ymin=0)
        else:
            for e in r["excited"]:
                if "error" in e:
                    self.store.append([e["label"], "", "", "", "", "failed: " + e["error"][:80]])
                    continue
                note = "collapsed to the ground state" if e.get("collapsed") else "overlap {:.3f}".format(e["overlap"])
                self.store.append([e["label"], "{:.3f}".format(e["energy_ev"]), "{:.1f}".format(e["energy_kj"]),
                                   "" if not e["nm"] else "{:.0f}".format(e["nm"]), "f = {:.4f}".format(e["f"]), note])
            good = [e for e in r["excited"] if "error" not in e and not e.get("collapsed")]
            if good:
                self.plot.set_lines([{"x": [e["energy_ev"] for e in good], "y": [e["f"] for e in good], "style": "bars"}],
                                    "excitation energy (eV)", "oscillator strength", "MOM ({})".format(r["method"]), ymin=0)
        self.btn_export.set_sensitive(True)

    def on_export(self, _b):
        dialog = Gtk.FileChooserDialog(title="Export excited states", parent=self.window, action=Gtk.FileChooserAction.SAVE)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        dialog.set_do_overwrite_confirmation(True)
        dialog.set_current_name("excited_states.txt")
        if dialog.run() == Gtk.ResponseType.OK:
            with open(dialog.get_filename(), "w") as handle:
                handle.write("# EasyHybrid excited states: {}\n".format(json.dumps({k: v for k, v in self.job.items() if k not in ("Z", "xyz")})))
                handle.write("# state\tdE_eV\tdE_kJmol\tnm\tS_or_f\tnote\n")
                for row in self.store:
                    handle.write("\t".join(row[:]) + "\n")
        dialog.destroy()

    def _message(self, text):
        d = Gtk.MessageDialog(transient_for=self.window, modal=True, message_type=Gtk.MessageType.INFO,
                              buttons=Gtk.ButtonsType.OK, text=text)
        d.run()
        d.destroy()

    def _on_destroy(self, *_):
        self.on_cancel(None)
        self.window = None
