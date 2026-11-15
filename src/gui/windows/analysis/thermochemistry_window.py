"""
thermochemistry_window
======================

Normal Modes window > Thermochemistry (RRHO)...: ideal-gas RRHO
thermochemistry (util/thermochemistry.py) from the imported frequencies,
at any temperature and pressure. The rotational symmetry number is
suggested from the point group found by pDynamo (Find3DGraphPointGroup).
"""

import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk

from util import thermochemistry as TC


def find_point_group(system, distance_tolerance=0.1, angle_tolerance=0.1):
    """Point group label of a pDynamo system (pDynamo's PointGroupFinder)."""
    from pCore import Clone
    from pScientific.Arrays import Array
    from pScientific.Symmetry import Find3DGraphPointGroup
    coordinates = Clone(system.coordinates3)
    report = Find3DGraphPointGroup([a.atomicNumber for a in system.atoms], coordinates,
                                   doCharacterSymmetryOperations=False, log=None,
                                   weights=Array.FromIterable([a.mass for a in system.atoms]),
                                   distanceTolerance=distance_tolerance, angleTolerance=angle_tolerance)
    return report["Point Group"].label, report


class ThermochemistryWindow:

    def __init__(self, main, system, frequencies, title=""):
        self.main = main
        self.system = system
        self.frequencies = list(frequencies)
        self.window = Gtk.Window(title="Thermochemistry (RRHO){}".format(" - " + title if title else ""))
        self.window.set_default_size(640, 560)
        self.window.set_transient_for(main.window)

        box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8, margin=10)
        self.window.add(box)
        grid = Gtk.Grid(column_spacing=10, row_spacing=6)
        box.pack_start(grid, False, False, 0)

        def spin(lo, hi, step, value, digits):
            s = Gtk.SpinButton.new_with_range(lo, hi, step)
            s.set_digits(digits)
            s.set_value(value)
            return s

        self.spin_t = spin(1, 5000, 1, 298.15, 2)
        self.spin_p = spin(0.001, 1000, 0.1, 1.0, 3)
        self.spin_sigma = spin(1, 60, 1, 1, 0)
        multiplicity = 1
        try:
            multiplicity = int(system.electronicState.multiplicity)
        except Exception:
            pass
        self.spin_mult = spin(1, 10, 1, multiplicity, 0)
        self.check_qrrho = Gtk.CheckButton(label="Quasi-RRHO for low modes (Grimme 2012, w0 = 100 cm-1)")
        self.label_pg = Gtk.Label(xalign=0)
        rows = (("Temperature (K)", self.spin_t), ("Pressure (atm)", self.spin_p),
                ("Symmetry number", self.spin_sigma), ("Spin multiplicity", self.spin_mult))
        for k, (text, w) in enumerate(rows):
            grid.attach(Gtk.Label(label=text, xalign=0), 0, k, 1, 1)
            grid.attach(w, 1, k, 1, 1)
        grid.attach(self.label_pg, 2, 2, 1, 1)
        grid.attach(self.check_qrrho, 0, 4, 3, 1)

        btns = Gtk.Box(spacing=6)
        calc = Gtk.Button(label="Calculate")
        calc.connect("clicked", self.on_calculate)
        btns.pack_start(calc, False, False, 0)
        close = Gtk.Button(label="Close")
        close.connect("clicked", lambda *_: self.window.destroy())
        btns.pack_end(close, False, False, 0)
        box.pack_start(btns, False, False, 0)

        self.text = Gtk.TextView(editable=False, monospace=True)
        sw = Gtk.ScrolledWindow()
        sw.add(self.text)
        box.pack_start(sw, True, True, 0)

        self.partial = system.freeAtoms is not None and len(system.freeAtoms) < len(system.atoms)
        self._suggest_symmetry()
        self.window.show_all()
        self.on_calculate(None)

    def _suggest_symmetry(self):
        if self.partial:
            self.label_pg.set_text("partial Hessian (fixed atoms): vibrations only")
            return
        if len(self.system.atoms) > 300:
            self.label_pg.set_text("point group not searched (> 300 atoms)")
            return
        try:
            label, _report = find_point_group(self.system)
            sigma = TC.symmetry_number(label)
            self.spin_sigma.set_value(sigma)
            self.label_pg.set_text("point group {} (pDynamo) -> {}".format(label, sigma))
        except Exception as error:
            self.label_pg.set_text("point group not found ({})".format(error))

    def result(self):
        system = self.system
        if self.partial:
            idx = list(system.freeAtoms)
            freqs = [f for f in self.frequencies]
            return TC.rrho(freqs, [system.atoms[i].mass for i in idx],
                           [[system.coordinates3[i, k] for k in range(3)] for i in idx],
                           self.spin_t.get_value(), self.spin_p.get_value(), int(self.spin_sigma.get_value()),
                           int(self.spin_mult.get_value()), n_rt=0 if len(freqs) != 3 * len(idx) else None,
                           quasi_rrho=self.check_qrrho.get_active())
        n = len(system.atoms)
        if len(self.frequencies) != 3 * n:
            raise ValueError("{} modes are imported, but the system has 3 x {} = {} degrees of freedom. "
                             "Import all the mode folders of the Normal Modes job.".format(len(self.frequencies), n, 3 * n))
        return TC.rrho(self.frequencies, [a.mass for a in system.atoms],
                       [[system.coordinates3[i, k] for k in range(3)] for i in range(n)],
                       self.spin_t.get_value(), self.spin_p.get_value(), int(self.spin_sigma.get_value()),
                       int(self.spin_mult.get_value()), quasi_rrho=self.check_qrrho.get_active())

    def on_calculate(self, _button):
        try:
            t = self.result()
            text = TC.format_table(t)
            if self.partial:
                text = ("Partial Hessian (the system has fixed atoms): only the vibrational part is "
                        "meaningful; translation and rotation below are those of the free atoms "
                        "as an isolated molecule and should be ignored.\n\n") + text
        except Exception as error:
            text = "Could not compute the thermochemistry:\n{}".format(error)
        self.text.get_buffer().set_text(text)
