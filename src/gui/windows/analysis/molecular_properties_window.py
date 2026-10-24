"""
molecular_properties_window
===========================

Analysis > Molecular Properties...: pDynamo tools that describe one
structure (the frame shown in the 3D view):

  Symmetry          point group and symmetry-equivalent atoms
                    (pScientific.Symmetry.Find3DGraphPointGroup)
  Stereochemistry   Cahn-Ingold-Prelog R/S and E/Z labels (pSimulation.CIPLabelFinder)
  Ion mobility      projection-approximation and hard-sphere collision cross
                    sections and mobilities in He (pSimulation.HardSphereIonMobilities,
                    Shvartsburg & Jarrold 1996)

The atoms are the active system or the selected molecules (pruned copy).
"""

import threading
import time
import traceback
import types

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, GLib

from util import thermochemistry as TC


class MolecularPropertiesWindow:

    def __init__(self, main):
        self.main = main
        self.p_session = main.p_session
        self.vm_session = main.vm_session
        self.window = None
        self._thread = None

    # ------------------------------------------------------------------ window
    def open_window(self):
        if self.window is not None:
            self.window.present()
            return
        self.window = Gtk.Window(title="Molecular Properties")
        self.window.set_default_size(720, 560)
        self.window.set_transient_for(self.main.window)
        self.window.connect("destroy", self._on_destroy)
        box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8, margin=10)
        self.window.add(box)

        top = Gtk.Box(spacing=8)
        self.radio_system = Gtk.RadioButton.new_with_label(None, "Active system")
        self.radio_sel = Gtk.RadioButton.new_with_label_from_widget(self.radio_system, "Selected molecule(s)")
        top.pack_start(Gtk.Label(label="Atoms:"), False, False, 0)
        top.pack_start(self.radio_system, False, False, 0)
        top.pack_start(self.radio_sel, False, False, 0)
        self.label_target = Gtk.Label(xalign=0)
        top.pack_start(self.label_target, True, True, 0)
        box.pack_start(top, False, False, 0)

        self.notebook = Gtk.Notebook()
        box.pack_start(self.notebook, True, True, 0)
        self.notebook.append_page(self._page_symmetry(), Gtk.Label(label="Symmetry"))
        self.notebook.append_page(self._page_cip(), Gtk.Label(label="Stereochemistry (CIP)"))
        self.notebook.append_page(self._page_mobility(), Gtk.Label(label="Ion mobility (CCS)"))

        bottom = Gtk.Box(spacing=8)
        close = Gtk.Button(label="Close")
        close.connect("clicked", lambda *_: self.window.destroy())
        bottom.pack_end(close, False, False, 0)
        box.pack_start(bottom, False, False, 0)
        self.window.show_all()

    def _text(self):
        tv = Gtk.TextView(editable=False, monospace=True, wrap_mode=Gtk.WrapMode.WORD)
        sw = Gtk.ScrolledWindow()
        sw.add(tv)
        return tv, sw

    def _page_symmetry(self):
        page = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6, margin=8)
        grid = Gtk.Grid(column_spacing=8, row_spacing=6)
        self.spin_dtol = Gtk.SpinButton.new_with_range(0.001, 1.0, 0.01)
        self.spin_dtol.set_digits(3)
        self.spin_dtol.set_value(0.1)
        self.spin_atol = Gtk.SpinButton.new_with_range(0.1, 30.0, 0.5)
        self.spin_atol.set_digits(1)
        self.spin_atol.set_value(5.7)
        grid.attach(Gtk.Label(label="Distance tolerance (A)", xalign=0), 0, 0, 1, 1)
        grid.attach(self.spin_dtol, 1, 0, 1, 1)
        grid.attach(Gtk.Label(label="Angle tolerance (deg)", xalign=0), 2, 0, 1, 1)
        grid.attach(self.spin_atol, 3, 0, 1, 1)
        run = Gtk.Button(label="Find Point Group")
        run.connect("clicked", lambda *_: self._safe(self.run_symmetry))
        grid.attach(run, 4, 0, 1, 1)
        page.pack_start(grid, False, False, 0)
        self.text_sym, sw = self._text()
        page.pack_start(sw, True, True, 0)
        return page

    def _page_cip(self):
        page = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6, margin=8)
        hb = Gtk.Box(spacing=8)
        run = Gtk.Button(label="Find R/S and E/Z")
        run.connect("clicked", lambda *_: self._safe(self.run_cip))
        self.btn_cip_labels = Gtk.Button(label="Show as Labels")
        self.btn_cip_labels.set_sensitive(False)
        self.btn_cip_labels.connect("clicked", lambda *_: self._safe(self.show_cip_labels))
        hb.pack_start(run, False, False, 0)
        hb.pack_start(self.btn_cip_labels, False, False, 0)
        page.pack_start(hb, False, False, 0)
        note = Gtk.Label(xalign=0)
        note.set_markup("<small>Needs all hydrogens and the bond orders (connectivity) of the system.</small>")
        page.pack_start(note, False, False, 0)
        self.text_cip, sw = self._text()
        page.pack_start(sw, True, True, 0)
        return page

    def _page_mobility(self):
        page = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6, margin=8)
        grid = Gtk.Grid(column_spacing=8, row_spacing=6)
        self.spin_traj = Gtk.SpinButton.new_with_range(1000, 5000000, 10000)
        self.spin_traj.set_value(100000)
        self.spin_refl = Gtk.SpinButton.new_with_range(1, 100, 1)
        self.spin_refl.set_value(30)
        self.spin_temp = Gtk.SpinButton.new_with_range(1, 1000, 1)
        self.spin_temp.set_value(298.0)
        self.spin_seed = Gtk.SpinButton.new_with_range(1, 1e9, 1)
        self.spin_seed.set_value(314159)
        for k, (text, w) in enumerate((("MC trajectories", self.spin_traj), ("Reflection limit", self.spin_refl),
                                       ("Temperature (K)", self.spin_temp), ("Random seed", self.spin_seed))):
            grid.attach(Gtk.Label(label=text, xalign=0), 2 * (k % 2), k // 2, 1, 1)
            grid.attach(w, 2 * (k % 2) + 1, k // 2, 1, 1)
        self.btn_ccs = Gtk.Button(label="Calculate")
        self.btn_ccs.connect("clicked", lambda *_: self._safe(self.run_mobility))
        grid.attach(self.btn_ccs, 4, 0, 1, 1)
        page.pack_start(grid, False, False, 0)
        self.progress_ccs = Gtk.ProgressBar(show_text=True)
        self.progress_ccs.set_text("")
        page.pack_start(self.progress_ccs, False, False, 0)
        note = Gtk.Label(xalign=0)
        note.set_line_wrap(True)
        note.set_markup("<small>Helium buffer gas, charge +1 (multiply the mobility by |z| for other charges). "
                        "pDynamo's hard-sphere radii cover H, C, N, O, Na, Si, P and Fe only.</small>")
        page.pack_start(note, False, False, 0)
        self.text_ccs, sw = self._text()
        page.pack_start(sw, True, True, 0)
        return page

    # ------------------------------------------------------------------ input
    def _target(self):
        """(system, vobject, indices) of the atoms, from the active system or the selection."""
        if self.radio_sel.get_active():
            sel = self.vm_session.selections[self.vm_session.current_selection]
            atoms = list(sel.selected_atoms)
            if not atoms:
                raise ValueError("Nothing is selected.")
            vobjects = {id(a.vm_object) for a in atoms}
            if len(vobjects) > 1:
                raise ValueError("The selection spans several objects.")
            vobject = atoms[0].vm_object
            system = self.p_session.psystem.get(vobject.e_id)
            indices = sorted({a.index - 1 for a in atoms})
        else:
            system = self.p_session.psystem.get(self.p_session.active_id)
            if system is None:
                raise ValueError("There is no active system.")
            objs = [v for v in self.vm_session.vm_objects_dic.values() if getattr(v, "e_id", None) == system.e_id]
            vobject = objs[-1] if objs else None
            indices = list(range(len(system.atoms)))
        if vobject is None or len(vobject.atoms) != len(system.atoms):
            raise ValueError("The 3D object does not match the system.")
        return system, vobject, indices

    def _molecule(self):
        """A pDynamo system with only the target atoms, at the coordinates shown."""
        from pCore import Clone, Selection
        from pSimulation import PruneByAtom
        system, vobject, indices = self._target()
        frame = min(int(self.vm_session.get_frame() or 0), len(vobject.frames) - 1)
        X = vobject.frames[frame]
        if len(indices) == len(system.atoms):
            mol = Clone(system)
        else:
            try:
                mol = PruneByAtom(system, Selection.FromIterable(indices))
            except Exception as error:
                raise ValueError("Could not make a copy of the selected atoms ({}).".format(error))
        for k, i in enumerate(indices):
            for d in range(3):
                mol.coordinates3[k, d] = float(X[i][d])
        self.label_target.set_text("{}: {} atoms".format(system.label, len(indices)))
        return mol, system, vobject, indices

    # ------------------------------------------------------------------ symmetry
    def run_symmetry(self):
        import math
        from pCore import Clone
        from pScientific.Arrays import Array
        from pScientific.Symmetry import Find3DGraphPointGroup
        mol, system, vobject, indices = self._molecule()
        if len(indices) > 500:
            raise ValueError("{} atoms: the point group search is meant for molecules (at most 500 atoms).".format(len(indices)))
        report = Find3DGraphPointGroup([a.atomicNumber for a in mol.atoms], Clone(mol.coordinates3),
                                       doCharacterSymmetryOperations=False, log=None,
                                       weights=Array.FromIterable([a.mass for a in mol.atoms]),
                                       distanceTolerance=self.spin_dtol.get_value(),
                                       angleTolerance=math.radians(self.spin_atol.get_value()))
        pg = report["Point Group"].label
        ops = report["Symmetry Operations"]
        lines = ["Point group: {}".format(pg),
                 "Rotational symmetry number: {}".format(TC.symmetry_number(pg)),
                 "Symmetry operations: {}".format(", ".join("{} ({})".format(k, len(v) if hasattr(v, "__len__") else 1)
                                                            for k, v in ops.items())),
                 ""]
        # . "Node Groups": {atomic number: [[equivalent atoms], ...]}
        classes = []
        for _element, groups in (report.get("Node Groups") or {}).items():
            for group in groups:
                if len(group) > 1:
                    classes.append(sorted(int(i) for i in group))
        if classes:
            lines.append("Symmetry-equivalent atoms:")
            for group in sorted(classes):
                lines.append("  " + ", ".join("{}{}".format(mol.atoms[i].label, indices[i] + 1) for i in group))
        else:
            lines.append("No symmetry-equivalent atoms.")
        lines.append("")
        lines.append("Tolerances {:.3f} A and {:.1f} deg (pDynamo PointGroupFinder). Raise them for slightly "
                     "distorted structures.".format(self.spin_dtol.get_value(), self.spin_atol.get_value()))
        self.text_sym.get_buffer().set_text("\n".join(lines))
        return pg

    # ------------------------------------------------------------------ CIP
    def run_cip(self):
        from pSimulation import CIPLabelFinder
        mol, system, vobject, indices = self._molecule()
        if len(mol.connectivity.bonds) == 0:
            raise ValueError("The system has no bonds (connectivity): CIP labels need the bond orders.")
        result = CIPLabelFinder(mol, log=None)
        if result is None:
            raise ValueError("No CIP result (no connectivity).")
        (tc, rc, sc, uc), (dc, ec, zc, udc) = result
        labels = {}
        def name(i):
            a = mol.atoms[i]
            res = getattr(a.parent, "label", "")
            return "{} {} (atom {})".format(res, a.label, indices[i] + 1)
        nodeIndices = mol.connectivity.nodeIndices
        def atoms_of(entry):
            # . pDynamo gives the E/Z centers as bond indices (2026-10-24: were read as atoms)
            bond = mol.connectivity.bonds[entry]
            return [nodeIndices[bond.node1], nodeIndices[bond.node2]]
        lines = ["Tetrahedral centers: {}   R: {}   S: {}   undefined: {}".format(len(tc), len(rc), len(sc), len(uc))]
        for i in rc:
            lines.append("  R  " + name(i)); labels[indices[i]] = "R"
        for i in sc:
            lines.append("  S  " + name(i)); labels[indices[i]] = "S"
        lines.append("")
        lines.append("Double-bond centers: {}   E: {}   Z: {}   undefined: {}".format(len(dc), len(ec), len(zc), len(udc)))
        for tag, group in (("E", ec), ("Z", zc)):
            for entry in group:
                at = atoms_of(entry)
                lines.append("  {}  ".format(tag) + " = ".join(name(i) for i in at))
                for i in at:
                    labels[indices[i]] = tag
        if uc:
            lines.append("")
            lines.append("Tetrahedral centers with two identical substituents (no label): {}".format(len(uc)))
        self.text_cip.get_buffer().set_text("\n".join(lines))
        self._cip = (vobject, labels)
        self.btn_cip_labels.set_sensitive(bool(labels))
        return labels

    def show_cip_labels(self):
        vobject, labels = self._cip
        atoms = [vobject.atoms[i] for i in labels]
        target = types.SimpleNamespace(selected_atoms=set(atoms), selected_objects=[vobject])
        self.vm_session.show_or_hide(rep_type="labels", selection=target, show=True)
        for i, tag in labels.items():
            vobject.atoms[i].label_text = tag
        self.vm_session.vm_glcore.queue_draw()

    # ------------------------------------------------------------------ mobility
    def run_mobility(self):
        if self._thread is not None:
            return
        mol, system, vobject, indices = self._molecule()
        known = {1, 6, 7, 8, 11, 14, 15, 26}
        bad = sorted({a.atomicNumber for a in mol.atoms if a.atomicNumber not in known})
        if bad:
            from pScientific import PeriodicTable
            raise ValueError("pDynamo has no hard-sphere radius for: {}.".format(
                ", ".join(PeriodicTable.Symbol(z) for z in bad)))
        ntraj, nrefl = int(self.spin_traj.get_value()), int(self.spin_refl.get_value())
        temp, seed = self.spin_temp.get_value(), int(self.spin_seed.get_value())
        self.btn_ccs.set_sensitive(False)
        self._t0 = time.time()

        def work():
            try:
                from pScientific.RandomNumbers import RandomNumberGenerator
                from pSimulation import HardSphereIonMobilities
                r = HardSphereIonMobilities(mol, nreflections=nrefl, ntrajectories=ntraj, temperature=temp,
                                            randomNumberGenerator=RandomNumberGenerator.WithSeed(seed), log=None)
                GLib.idle_add(self._mobility_done, r, None)
            except Exception as error:
                traceback.print_exc()
                GLib.idle_add(self._mobility_done, None, str(error))
        self._thread = threading.Thread(target=work, daemon=True)
        self._thread.start()
        GLib.timeout_add(200, self._pulse)

    def _pulse(self):
        if self._thread is None:
            return False
        self.progress_ccs.pulse()
        self.progress_ccs.set_text("running... {:.0f} s".format(time.time() - self._t0))
        return True

    def _mobility_done(self, r, error):
        self._thread = None
        self.btn_ccs.set_sensitive(True)
        if error:
            self.progress_ccs.set_text("failed")
            self._message(error)
            return False
        self.progress_ccs.set_fraction(1.0)
        self.progress_ccs.set_text("done in {:.1f} s".format(time.time() - self._t0))
        lines = ["Collision cross section, projection approximation (PA): {:9.2f} A^2".format(r["PA Cross-Section"]),
                 "Collision cross section, exact hard spheres (EHS):      {:9.2f} A^2".format(r["HS Cross-Section"]),
                 "Mobility (PA):  {:.4f} cm^2/(V s)".format(r["PA Mobility"] * 1e4),
                 "Mobility (EHS): {:.4f} cm^2/(V s)".format(r["HS Mobility"] * 1e4),
                 "",
                 "{} Monte Carlo trajectories, reflection limit {}, at most {} reflections used, {:.0f} K, He, z = +1.".format(
                     r["MC Trajectories"], r["Reflection Limit"], r["Maximum Reflections"], self.spin_temp.get_value()),
                 "Shvartsburg & Jarrold, Chem. Phys. Lett. 261, 86 (1996). Run with another seed to estimate the "
                 "Monte Carlo error."]
        self.text_ccs.get_buffer().set_text("\n".join(lines))
        self.last_mobility = r
        return False

    # ------------------------------------------------------------------ misc
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
        self.window = None
