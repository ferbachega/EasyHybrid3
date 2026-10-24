"""
orbitals_surfaces_window
========================

Analysis > Surfaces (2026-10-24): molecular orbitals and QC surfaces in the
layout of the other analysis windows (Trajectory Analysis, Umbrella Sampling
Browser).

  left    the kind of surface (orbitals, density, potential, MEP, cube file),
          its options and Render
  right   the orbitals of the chosen frame: an energy-level diagram (click a
          level), the atom contributions of the selected orbital, and the
          tables Orbitals, Composition (LCAO) and Summary

The wave function is computed on a copy of the system (util/orbital_analysis.py),
in a thread. The surfaces reuse generate_grid_parallel() of the previous window
(surface_analysis_window.py) and become objects of the Objects tree as before.
"""

import math
import multiprocessing
import os
import threading
import time
import traceback

import cairo
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk, GLib

import numpy as np

from gui.widgets.custom_widgets import SystemComboBox
from gui.widgets.custom_widgets import CoordinatesComboBox
from util import orbital_analysis as OA
from util.analysis_plot import AnalysisPlot, PALETTE
from util.colormaps import COLOR_MAPS

EV = OA.HARTREE_TO_EV
MAX_WORKERS = 12

SURFACE_TYPES = [
    ("orbital",   "Molecular orbitals"),
    ("density",   "Electron density"),
    ("potential", "Electrostatic potential"),
    ("mep",       "MEP (density + potential)"),
    ("cube",      "Cube file"),
]

HINTS = {
    "orbital":   "Click a level of the diagram or a row of the table, then Render. Isovalue in a.u. "
                 "(0.03-0.06 for orbitals).",
    "density":   "The total electron density of the QC region at the isovalue (e/bohr^3).",
    "potential": "Isosurfaces of the electrostatic potential of the QC region at +/- the isovalue (a.u.).",
    "mep":       "The density isosurface (0.002-0.02) colored by the electrostatic potential: red negative, "
                 "blue positive with the diverging colormaps.",
    "cube":      "A density or orbital from a .cube file (e.g. ORCA's orca_plot) of this system.",
}

ELEMENT_COLORS = {"H": (0.75, 0.75, 0.75), "C": (0.35, 0.35, 0.35), "N": (0.19, 0.31, 0.97),
                  "O": (0.90, 0.13, 0.13), "S": (0.85, 0.70, 0.10), "P": (1.0, 0.5, 0.0),
                  "F": (0.56, 0.88, 0.31), "Cl": (0.12, 0.94, 0.12), "Br": (0.65, 0.16, 0.16)}

OCCUPIED_COLOR = (0.12, 0.47, 0.71)
VIRTUAL_COLOR = (0.85, 0.37, 0.01)


# =============================================================================
#  Energy-level diagram
# =============================================================================
class EnergyLevelDiagram(Gtk.DrawingArea):
    """ Orbital energies as levels (occupied blue with electrons, virtual orange) and the
        HOMO-LUMO gap, on an energy range the user controls (2026-10-24):

          wheel          zoom around the pointer
          drag           move the range up or down
          double-click   back to the automatic range (span levels around HOMO/LUMO)
          set_range()    a range typed by the user

        Levels closer than MIN_GAP pixels are drawn side by side (as degenerate levels
        are), so they never overlap; a click selects one (on_select(k)) and
        on_range_changed(e0, e1) reports the range shown. """

    MIN_GAP = 7.0          # . px between two levels drawn one above the other
    MAX_PER_ROW = 12       # . above this, a row of levels is drawn as a band

    def __init__(self):
        super().__init__()
        self.energies = None          # eV
        self.occupancies = None
        self.homo = -1
        self.selected = None
        self.hover = None
        self.span = 7                 # . automatic range: levels below the HOMO and above the LUMO
        self.view = None              # . (e0, e1) in eV, None = automatic
        self.on_select = None
        self.on_range_changed = None
        self.label_of = lambda k: str(k)
        self._hit = []                # . (x0, x1, y, k)
        self._geom = None             # . (y0, h, e0, e1) of the last drawing
        self._drag = None
        self.add_events(Gdk.EventMask.POINTER_MOTION_MASK | Gdk.EventMask.LEAVE_NOTIFY_MASK |
                        Gdk.EventMask.BUTTON_PRESS_MASK | Gdk.EventMask.BUTTON_RELEASE_MASK |
                        Gdk.EventMask.SCROLL_MASK)
        self.connect("draw", self._on_draw)
        self.connect("motion-notify-event", self._on_motion)
        self.connect("leave-notify-event", lambda *_: self._set_hover(None))
        self.connect("button-press-event", self._on_button)
        self.connect("button-release-event", self._on_release)
        self.connect("scroll-event", self._on_scroll)
        self.set_size_request(250, 300)

    # ------------------------------------------------------------------ data and range
    def set_data(self, energies_ev, occupancies, homo, label_of):
        self.energies = np.asarray(energies_ev, dtype=float)
        self.occupancies = np.asarray(occupancies, dtype=float)
        self.homo = int(homo)
        self.label_of = label_of
        self.selected = self.homo if self.homo >= 0 else None
        self.view = None
        self._emit_range()
        self.queue_draw()

    def clear(self):
        self.energies = None
        self.queue_draw()

    def auto_range(self):
        n = len(self.energies)
        lo = max(0, self.homo - self.span)
        hi = min(n - 1, self.homo + 1 + self.span)
        e0, e1 = float(self.energies[lo]), float(self.energies[hi])
        if e1 - e0 < 1e-6:
            e0, e1 = e0 - 1, e1 + 1
        pad = 0.06 * (e1 - e0)
        return e0 - pad, e1 + pad

    def current_range(self):
        return self.view if self.view is not None else self.auto_range()

    def set_range(self, e0, e1):
        """ A range typed by the user (eV). """
        if e1 - e0 > 1e-3:
            self.view = (float(e0), float(e1))
            self._emit_range()
            self.queue_draw()

    def set_span(self, span):
        self.span = max(1, int(span))
        self.view = None
        self._emit_range()
        self.queue_draw()

    def reset_range(self):
        self.view = None
        self._emit_range()
        self.queue_draw()

    def _emit_range(self):
        if self.on_range_changed and self.energies is not None and len(self.energies):
            self.on_range_changed(*self.current_range())

    def select(self, k):
        self.selected = k
        if self.energies is not None and k is not None:
            e0, e1 = self.current_range()
            e = self.energies[k]
            if not e0 <= e <= e1:          # . bring the level into view, keeping the zoom
                half = (e1 - e0) / 2
                self.view = (e - half, e + half)
                self._emit_range()
        self.queue_draw()

    # ------------------------------------------------------------------ mouse
    def _energy_at(self, y):
        if self._geom is None:
            return None
        y0, h, e0, e1 = self._geom
        return e0 + (y0 + h - y) / h * (e1 - e0)

    def _on_scroll(self, _w, event):
        if self.energies is None:
            return False
        up = event.direction == Gdk.ScrollDirection.UP
        if event.direction == Gdk.ScrollDirection.SMOOTH:
            if abs(event.delta_y) < 1e-6:
                return False
            up = event.delta_y < 0
        elif event.direction not in (Gdk.ScrollDirection.UP, Gdk.ScrollDirection.DOWN):
            return False
        e0, e1 = self.current_range()
        center = self._energy_at(event.y)
        if center is None:
            center = (e0 + e1) / 2
        f = 1 / 1.25 if up else 1.25            # . up = zoom in
        span = max(0.02, min((e1 - e0) * f, 1e4))
        t = (center - e0) / (e1 - e0)
        self.view = (center - t * span, center - t * span + span)
        self._emit_range()
        self.queue_draw()
        return True

    def _on_button(self, _w, event):
        if self.energies is None:
            return False
        if event.type == Gdk.EventType._2BUTTON_PRESS and event.button == 1:
            self._drag = None
            self.reset_range()
            return True
        if event.button == 1:
            self._drag = {"y": event.y, "range": self.current_range(), "moved": False, "x": event.x}
        return True

    def _on_release(self, _w, event):
        if event.button != 1 or self._drag is None:
            return False
        drag, self._drag = self._drag, None
        if not drag["moved"]:
            k = self._locate(event.x, event.y)
            if k is not None:
                self.select(k)
                if self.on_select:
                    self.on_select(k)
        return True

    def _on_motion(self, _w, event):
        if self.energies is None:
            return False
        if self._drag is not None and (event.state & Gdk.ModifierType.BUTTON1_MASK):
            dy = event.y - self._drag["y"]
            if abs(dy) > 3 or self._drag["moved"]:
                self._drag["moved"] = True
                e0, e1 = self._drag["range"]
                h = self._geom[1] if self._geom else 300
                de = dy / h * (e1 - e0)
                self.view = (e0 + de, e1 + de)
                self._emit_range()
                self.queue_draw()
                return True
        return self._set_hover(self._locate(event.x, event.y))

    def _locate(self, ex, ey):
        best, dist = None, 6.0
        for a, b, y, k in self._hit:
            if a <= ex <= b and abs(ey - y) < dist:
                best, dist = k, abs(ey - y)
        return best

    def _set_hover(self, k):
        if k != self.hover:
            self.hover = k
            self.queue_draw()
        return False

    # ------------------------------------------------------------------ drawing
    def _on_draw(self, widget, cr, width=None, height=None):
        width = width or widget.get_allocated_width()
        height = height or widget.get_allocated_height()
        cr.set_source_rgb(1, 1, 1)
        cr.paint()
        cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        if self.energies is None or len(self.energies) == 0:
            cr.set_source_rgb(0.45, 0.45, 0.45)
            cr.set_font_size(12)
            text = "No orbitals computed"
            e = cr.text_extents(text)
            cr.move_to(width / 2 - e.width / 2, height / 2)
            cr.show_text(text)
            return False
        e0, e1 = self.current_range()
        x0, y0, w, h = 56, 30, width - 56 - 96, height - 30 - 24
        self._geom = (y0, h, e0, e1)
        Y = lambda e: y0 + h - (e - e0) / (e1 - e0) * h
        # . axis
        cr.set_font_size(10.5)
        cr.set_line_width(1)
        step = _nice_step(e1 - e0)
        v = math.ceil(e0 / step) * step
        while v <= e1:
            y = Y(v)
            cr.set_source_rgb(0.92, 0.92, 0.92); cr.move_to(x0, y); cr.line_to(x0 + w, y); cr.stroke()
            cr.set_source_rgb(0.2, 0.2, 0.2)
            t = "{:g}".format(round(v, 6)); ex = cr.text_extents(t)
            cr.move_to(x0 - ex.width - 6, y + ex.height / 2); cr.show_text(t)
            v += step
        cr.save(); cr.set_font_size(12); t = "E (eV)"; ex = cr.text_extents(t)
        cr.move_to(14, y0 + h / 2 + ex.width / 2); cr.rotate(-math.pi / 2); cr.show_text(t); cr.restore()
        cr.set_font_size(13); cr.set_source_rgb(0.1, 0.1, 0.1)
        t = "Orbital energies"; ex = cr.text_extents(t)
        cr.move_to(x0 + w / 2 - ex.width / 2 + 40, 18); cr.show_text(t)
        # . the levels in range, in rows: a level closer than MIN_GAP px to the first level
        #   of the current row goes beside it (degenerate or just too close at this scale)
        visible = [k for k in np.argsort(self.energies, kind="stable") if e0 <= self.energies[k] <= e1]
        rows = []
        for k in visible:
            y = Y(self.energies[k])
            if rows and abs(y - rows[-1][0]) < self.MIN_GAP:
                rows[-1][1].append(k)
            else:
                rows.append([y, [k]])
        self._hit = []
        cr.save()
        cr.rectangle(x0, y0 - 2, w + 4, h + 4)
        cr.clip()
        cx0, cx1 = x0 + 14, x0 + w - 8
        # . one slot width for every row (the widest row sets it), so that the levels drawn
        #   side by side line up in columns; each row is centred
        m_max = max([len(g) for _y, g in rows] or [1])
        gap = 6 if m_max <= 4 else 3 if m_max <= self.MAX_PER_ROW else 1
        lw = min((cx1 - cx0 - gap * (m_max - 1)) / m_max, 90)
        for _y, group in rows:
            m = len(group)
            total = lw * m + gap * (m - 1)
            start = (cx0 + cx1) / 2 - total / 2
            for j, k in enumerate(group):
                a = start + j * (lw + gap)
                b = a + lw
                y = Y(self.energies[k])
                occupied = self.occupancies[k] > 1e-3
                if k == self.selected or k == self.hover:
                    cr.set_source_rgba(1.0, 0.85, 0.2, 0.55 if k == self.selected else 0.3)
                    cr.rectangle(a - 3, y - 6, b - a + 6, 12); cr.fill()
                cr.set_source_rgb(*(OCCUPIED_COLOR if occupied else VIRTUAL_COLOR))
                cr.set_line_width(3.0 if k == self.selected else 2.0)
                cr.move_to(a, y); cr.line_to(b, y); cr.stroke()
                if occupied and lw >= 14 and rows_spacing_ok(rows, _y, 13):
                    self._electrons(cr, (a + b) / 2, y, self.occupancies[k])
                self._hit.append((a - 2, b + 2, y, k))
        cr.restore()
        # . labels on the right: HOMO, LUMO and the selected level (no two on top of each other)
        cr.set_font_size(10.5)
        placed = []
        for k in (self.selected, self.homo, self.homo + 1):
            if k is None or k < 0 or k >= len(self.energies) or k in [p[1] for p in placed]:
                continue
            if not e0 <= self.energies[k] <= e1:
                continue
            y = Y(self.energies[k]) + 4
            if any(abs(y - py) < 12 for py, _k in placed):
                continue
            placed.append((y, k))
            bold = k == self.selected
            cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL,
                                cairo.FONT_WEIGHT_BOLD if bold else cairo.FONT_WEIGHT_NORMAL)
            cr.set_source_rgb(0.15, 0.15, 0.15)
            cr.move_to(x0 + w + 6, y); cr.show_text("{} {:.2f}".format(self.label_of(k), self.energies[k]))
        cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        # . HOMO-LUMO gap
        if 0 <= self.homo < len(self.energies) - 1:
            eh, el = self.energies[self.homo], self.energies[self.homo + 1]
            ya, yb = Y(eh), Y(el)
            if e0 <= eh <= e1 and e0 <= el <= e1 and ya - yb > 60:
                gx = x0 + w + 14
                cr.set_source_rgb(0.45, 0.45, 0.45); cr.set_line_width(1)
                cr.set_dash([3, 3]); cr.move_to(gx, ya - 8); cr.line_to(gx, yb + 8); cr.stroke(); cr.set_dash([])
                for yy, sgn in ((yb + 8, 1), (ya - 8, -1)):
                    cr.move_to(gx, yy - sgn * 4); cr.line_to(gx - 3, yy + sgn * 1); cr.move_to(gx, yy - sgn * 4)
                    cr.line_to(gx + 3, yy + sgn * 1); cr.stroke()
                cr.set_font_size(10)
                cr.move_to(gx + 6, (ya + yb) / 2 - 2); cr.show_text("gap")
                cr.move_to(gx + 6, (ya + yb) / 2 + 11); cr.show_text("{:.2f} eV".format(el - eh))
        # . footer
        cr.set_font_size(9.5); cr.set_source_rgb(0.45, 0.45, 0.45)
        t = "{} of {} levels   wheel: zoom   drag: move   double-click: reset".format(len(visible), len(self.energies))
        cr.move_to(8, height - 6); cr.show_text(t)
        return False

    def _electrons(self, cr, x, y, occupancy):
        n = 2 if occupancy > 1.5 else 1
        cr.set_source_rgb(0.1, 0.1, 0.1)
        cr.set_line_width(1.2)
        offsets = (-4, 4) if n == 2 else (0,)
        for j, dx in enumerate(offsets):
            up = j == 0
            xx = x + dx
            ya, yb = (y + 5, y - 5) if up else (y - 5, y + 5)
            cr.move_to(xx, ya); cr.line_to(xx, yb); cr.stroke()
            s = -1 if up else 1
            cr.move_to(xx, yb); cr.line_to(xx - 2.5, yb - s * 3); cr.stroke()

    def export(self, path, width=520, height=700):
        surface = cairo.SVGSurface(path, width, height) if path.lower().endswith(".svg") else \
            cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
        cr = cairo.Context(surface)
        hover, self.hover = self.hover, None
        geom = self._geom
        self._on_draw(self, cr, width, height)
        self.hover, self._geom = hover, geom
        if path.lower().endswith(".svg"):
            surface.finish()
        else:
            surface.write_to_png(path)


def rows_spacing_ok(rows, y, minimum):
    """ True when the row at y has at least `minimum` px to its neighbours (room for the electrons). """
    for k, (ry, _g) in enumerate(rows):
        if ry == y:
            before = abs(y - rows[k - 1][0]) if k > 0 else 1e9
            after = abs(rows[k + 1][0] - y) if k + 1 < len(rows) else 1e9
            return min(before, after) >= minimum
    return True


def _nice_step(span, target=6):
    raw = span / max(target, 1)
    mag = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1
    for f in (1, 2, 2.5, 5, 10):
        if raw <= f * mag:
            return f * mag
    return 10 * mag


# =============================================================================
#  Window
# =============================================================================
class OrbitalsSurfacesWindow:

    def __init__(self, main=None):
        self.main = main
        self.p_session = main.p_session
        self.vm_session = main.vm_session
        self.Visible = False
        self.window = None
        self.wf = None                  # . wave function of the frame shown
        self.wf_key = None              # . (system_id, vobject_id, frame)
        self.selected = None            # . orbital index
        self._busy = False
        self.cube_path = None

    # ------------------------------------------------------------------ build
    def open_window(self):
        if self.Visible:
            self.window.present()
            return
        self.window = Gtk.Window(title="Orbitals and Surfaces")
        self.window.set_default_size(1180, 740)
        self.window.set_transient_for(self.main.window)
        self.window.connect("destroy", self._on_destroy)
        outer = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6, margin=8)
        self.window.add(outer)
        outer.pack_start(self._build_input_frame(), False, False, 0)

        paned = Gtk.Paned(orientation=Gtk.Orientation.HORIZONTAL)
        outer.pack_start(paned, True, True, 0)
        paned.pack1(self._build_left(), False, False)
        paned.pack2(self._build_right(), True, False)

        self.window.show_all()
        self.Visible = True
        if self.p_session.psystem:
            self.combobox_systems.set_active_system(e_id=self.p_session.active_id)
        self.type_view.get_selection().select_path(Gtk.TreePath(0))
        self._show_result_page()

    def _build_input_frame(self):
        frame = Gtk.Frame(label=" Structure ")
        grid = Gtk.Grid(column_spacing=8, row_spacing=6, margin=8)
        frame.add(grid)
        self.combobox_systems = SystemComboBox(self.main)
        self.combobox_systems.connect("changed", self._on_system_changed)
        self.coordinates_combobox = CoordinatesComboBox()
        self.coordinates_combobox.connect("changed", self._on_object_changed)
        self.spin_frame = Gtk.SpinButton.new_with_range(0, 0, 1)
        self.spin_frame.set_tooltip_text("The frame whose orbitals are listed (and rendered, with 'Current frame').")
        self.spin_frame.connect("value-changed", self._on_frame_changed)
        self.btn_compute = Gtk.Button(label="Compute Orbitals")
        self.btn_compute.set_tooltip_text("SCF at this frame with the QC model of the system, on a copy of the system.")
        self.btn_compute.connect("clicked", lambda *_: self.compute_orbitals())
        self.label_status = Gtk.Label(xalign=0)
        self.label_status.set_ellipsize(3)
        grid.attach(Gtk.Label(label="System", xalign=0), 0, 0, 1, 1)
        grid.attach(self.combobox_systems, 1, 0, 1, 1)
        grid.attach(Gtk.Label(label="Object", xalign=0), 2, 0, 1, 1)
        grid.attach(self.coordinates_combobox, 3, 0, 1, 1)
        grid.attach(Gtk.Label(label="Frame", xalign=0), 4, 0, 1, 1)
        grid.attach(self.spin_frame, 5, 0, 1, 1)
        grid.attach(self.btn_compute, 6, 0, 1, 1)
        grid.attach(self.label_status, 7, 0, 1, 1)
        self.label_status.set_hexpand(True)
        return frame

    def _build_left(self):
        left = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
        left.set_size_request(300, -1)
        store = Gtk.ListStore(str, str)
        for key, title in SURFACE_TYPES:
            store.append([title, key])
        self.type_view = Gtk.TreeView(model=store)
        self.type_view.set_headers_visible(False)
        self.type_view.append_column(Gtk.TreeViewColumn("Surface", Gtk.CellRendererText(), text=0))
        self.type_view.get_selection().connect("changed", self._on_type_changed)
        frame = Gtk.Frame(label=" Surface ")
        frame.add(self.type_view)
        left.pack_start(frame, False, False, 0)

        frame = Gtk.Frame(label=" Options ")
        box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6, margin=8)
        frame.add(box)
        self.options = Gtk.Grid(column_spacing=8, row_spacing=6)
        box.pack_start(self.options, False, False, 0)
        self.entry_iso = self._entry("0.05", "Isovalue of the surface, in atomic units.")
        self.entry_spacing = self._entry("0.4", "Grid spacing in bohr: smaller is finer and slower (0.3-0.5).")
        self.color_plus = Gtk.ColorButton.new_with_rgba(Gdk.RGBA(0.20, 0.40, 1.00, 1.0))
        self.color_minus = Gtk.ColorButton.new_with_rgba(Gdk.RGBA(0.95, 0.15, 0.15, 1.0))
        self.color_density = Gtk.ColorButton.new_with_rgba(Gdk.RGBA(0.55, 0.75, 0.95, 1.0))
        self.entry_mep_min = self._entry("", "Potential at the negative end of the color scale (a.u.). "
                                             "Empty: 2nd percentile of the values on the surface.", "auto")
        self.entry_mep_max = self._entry("", "Potential at the positive end of the color scale (a.u.). "
                                             "Empty: 98th percentile.", "auto")
        self.combo_cmap = Gtk.ComboBoxText()
        self._cmaps = list(COLOR_MAPS.keys())
        for name in self._cmaps:
            self.combo_cmap.append_text(name)
        self.combo_cmap.set_active(self._cmaps.index("coolwarm") if "coolwarm" in self._cmaps else 0)
        self.combo_cmap.set_tooltip_text("Diverging maps (coolwarm, vik, berlin, bam): red negative, blue positive.")
        self.entry_pot_spacing = self._entry("", "Grid of the potential only (bohr); it is interpolated onto the "
                                                 "density surface. Empty: 2.5 x the grid spacing.", "auto (2.5x)")
        self.cube_button = Gtk.FileChooserButton(title="Cube file", action=Gtk.FileChooserAction.OPEN)
        flt = Gtk.FileFilter(); flt.set_name("Cube files"); flt.add_pattern("*.cube"); flt.add_pattern("*.cub")
        self.cube_button.add_filter(flt)
        self.cube_button.connect("file-set", lambda b: setattr(self, "cube_path", b.get_filename()))
        self.radio_current = Gtk.RadioButton.new_with_label(None, "Current frame")
        self.radio_all = Gtk.RadioButton.new_with_label_from_widget(self.radio_current, "All frames")
        self.radio_all.set_tooltip_text("One surface per frame of the object: it follows the trajectory.")
        left.pack_start(frame, False, False, 0)

        hb = Gtk.Box(spacing=6)
        self.btn_render = Gtk.Button(label="Render")
        self.btn_render.connect("clicked", lambda *_: self.render())
        hb.pack_start(self.btn_render, True, True, 0)
        left.pack_start(hb, False, False, 0)
        self.progress = Gtk.ProgressBar(show_text=True)
        self.progress.set_text("")
        left.pack_start(self.progress, False, False, 0)
        self.label_hint = Gtk.Label(xalign=0)
        self.label_hint.set_line_wrap(True)
        self.label_hint.set_max_width_chars(40)
        left.pack_start(self.label_hint, False, False, 0)
        return left

    def _entry(self, text, tooltip, placeholder=None):
        e = Gtk.Entry(text=text, width_chars=9)
        e.set_tooltip_text(tooltip)
        if placeholder:
            e.set_placeholder_text(placeholder)
        return e

    def _build_right(self):
        right = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6, margin_start=6)
        self.right_stack = Gtk.Stack()
        right.pack_start(self.right_stack, True, True, 0)

        # . orbitals page
        page = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
        vpaned = Gtk.Paned(orientation=Gtk.Orientation.VERTICAL)
        top = Gtk.Paned(orientation=Gtk.Orientation.HORIZONTAL)
        self.diagram = EnergyLevelDiagram()
        self.diagram.on_select = self._on_diagram_select
        self.diagram.on_range_changed = self._on_diagram_range
        top.pack1(self._build_diagram_controls(), False, False)
        self.plot = AnalysisPlot()
        self.plot.on_click = self._on_plot_click
        top.pack2(self.plot, True, False)
        top.set_position(360)
        vpaned.pack1(top, True, False)

        nb = Gtk.Notebook()
        # . Orbitals
        self.orb_store = Gtk.ListStore(int, str, str, str, str, str)
        self.orb_view = Gtk.TreeView(model=self.orb_store)
        for k, title in enumerate(("#", "Orbital", "Occ.", "E (eV)", "E (hartree)", "Main contributions")):
            col = Gtk.TreeViewColumn(title, Gtk.CellRendererText(), text=k)
            col.set_resizable(True)
            self.orb_view.append_column(col)
        self.orb_view.get_selection().connect("changed", self._on_orbital_row)
        sw = Gtk.ScrolledWindow(); sw.add(self.orb_view)
        nb.append_page(sw, Gtk.Label(label="Orbitals"))
        # . Composition
        self.comp_store = Gtk.ListStore(str, str, str, str, int)
        self.comp_view = Gtk.TreeView(model=self.comp_store)
        for k, title in enumerate(("Atom", "Atomic orbital", "Coefficient", "Population (%)")):
            col = Gtk.TreeViewColumn(title, Gtk.CellRendererText(), text=k)
            col.set_resizable(True)
            self.comp_view.append_column(col)
        self.comp_view.connect("row-activated", self._on_comp_row)
        sw = Gtk.ScrolledWindow(); sw.add(self.comp_view)
        nb.append_page(sw, Gtk.Label(label="Composition (LCAO)"))
        # . Summary
        self.summary = Gtk.TextView(editable=False, monospace=True, wrap_mode=Gtk.WrapMode.WORD)
        sw = Gtk.ScrolledWindow(); sw.add(self.summary)
        nb.append_page(sw, Gtk.Label(label="Summary"))
        nb.set_size_request(-1, 210)
        vpaned.pack2(nb, False, False)
        vpaned.set_position(380)
        page.pack_start(vpaned, True, True, 0)

        actions = Gtk.Box(spacing=6)
        b = Gtk.Button(label="Render Orbital")
        b.set_tooltip_text("Render the selected orbital with the options on the left.")
        b.connect("clicked", lambda *_: self.render(kind="orbital"))
        actions.pack_start(b, False, False, 0)
        b = Gtk.Button(label="Select Main Atoms")
        b.set_tooltip_text("Select in the 3D view the atoms with at least 10% of the selected orbital.")
        b.connect("clicked", lambda *_: self.select_main_atoms())
        actions.pack_start(b, False, False, 0)
        b = Gtk.Button(label="Save Table...")
        b.set_tooltip_text("The orbitals (energies and main contributions) and the LCAO of the selected "
                           "orbital, as text.")
        b.connect("clicked", lambda *_: self.on_save_table())
        actions.pack_end(b, False, False, 0)
        b = Gtk.Button(label="Export Image...")
        b.set_tooltip_text("The energy-level diagram or the contribution plot (PNG or SVG).")
        b.connect("clicked", lambda *_: self.on_export_image())
        actions.pack_end(b, False, False, 0)
        page.pack_start(actions, False, False, 0)
        self.right_stack.add_named(page, "orbitals")

        # . other surfaces: a short report of the last surface
        self.info = Gtk.TextView(editable=False, monospace=True, wrap_mode=Gtk.WrapMode.WORD,
                                 left_margin=10, top_margin=10)
        sw = Gtk.ScrolledWindow(); sw.add(self.info)
        self.right_stack.add_named(sw, "info")
        return right

    def _build_diagram_controls(self):
        """ The energy range of the diagram: typed, or set by the mouse (kept in sync). """
        box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=4)
        grid = Gtk.Grid(column_spacing=4, row_spacing=4, margin_start=4, margin_top=2)
        self.spin_emin = Gtk.SpinButton.new_with_range(-1000, 1000, 0.5)
        self.spin_emax = Gtk.SpinButton.new_with_range(-1000, 1000, 0.5)
        for spin in (self.spin_emin, self.spin_emax):
            spin.set_digits(2)
            spin.set_width_chars(6)
            spin.connect("value-changed", self._on_range_typed)
        self.spin_emin.set_tooltip_text("Lowest energy shown (eV). The mouse wheel zooms, dragging moves the range.")
        self.spin_emax.set_tooltip_text("Highest energy shown (eV).")
        self.spin_levels = Gtk.SpinButton.new_with_range(1, 500, 1)
        self.spin_levels.set_value(self.diagram.span)
        self.spin_levels.set_tooltip_text("Automatic range: this number of levels below the HOMO and above the LUMO.")
        self.spin_levels.connect("value-changed", lambda w: self.diagram.set_span(int(w.get_value())))
        fit = Gtk.Button(label="Fit")
        fit.set_tooltip_text("Back to the automatic range around the HOMO and LUMO (also: double-click the diagram).")
        fit.connect("clicked", lambda *_: self.diagram.set_span(int(self.spin_levels.get_value())))
        grid.attach(Gtk.Label(label="E (eV)", xalign=0), 0, 0, 1, 1)
        grid.attach(self.spin_emin, 1, 0, 1, 1)
        grid.attach(Gtk.Label(label="to"), 2, 0, 1, 1)
        grid.attach(self.spin_emax, 3, 0, 1, 1)
        grid.attach(Gtk.Label(label="Levels", xalign=0), 0, 1, 1, 1)
        grid.attach(self.spin_levels, 1, 1, 1, 1)
        grid.attach(fit, 3, 1, 1, 1)
        box.pack_start(grid, False, False, 0)
        box.pack_start(self.diagram, True, True, 0)
        self._range_sync = False
        return box

    def _on_diagram_range(self, e0, e1):
        self._range_sync = True
        self.spin_emin.set_value(round(e0, 2))
        self.spin_emax.set_value(round(e1, 2))
        self._range_sync = False

    def _on_range_typed(self, _w):
        if self._range_sync:
            return
        e0, e1 = self.spin_emin.get_value(), self.spin_emax.get_value()
        if e1 > e0:
            self.diagram.set_range(e0, e1)

    # ------------------------------------------------------------------ options
    def _current_type(self):
        model, it = self.type_view.get_selection().get_selected()
        return model[it][1] if it is not None else "orbital"

    def _on_type_changed(self, _sel):
        kind = self._current_type()
        for child in self.options.get_children():
            self.options.remove(child)
        rows = []
        if kind == "orbital":
            self.entry_iso.set_text("0.05")
            rows = [("Isovalue", self.entry_iso), ("Grid spacing (bohr)", self.entry_spacing),
                    ("Positive lobe", self.color_plus), ("Negative lobe", self.color_minus)]
        elif kind == "density":
            self.entry_iso.set_text("0.05")
            rows = [("Isovalue", self.entry_iso), ("Grid spacing (bohr)", self.entry_spacing),
                    ("Color", self.color_density)]
        elif kind == "potential":
            self.entry_iso.set_text("0.05")
            rows = [("Isovalue (+/-)", self.entry_iso), ("Grid spacing (bohr)", self.entry_spacing),
                    ("Positive", self.color_plus), ("Negative", self.color_minus)]
        elif kind == "mep":
            self.entry_iso.set_text("0.02")
            rows = [("Density isovalue", self.entry_iso), ("Grid spacing (bohr)", self.entry_spacing),
                    ("MEP min (a.u.)", self.entry_mep_min), ("MEP max (a.u.)", self.entry_mep_max),
                    ("Colormap", self.combo_cmap), ("Potential spacing", self.entry_pot_spacing)]
        elif kind == "cube":
            self.entry_iso.set_text("0.05")
            rows = [("File", self.cube_button), ("Isovalue", self.entry_iso),
                    ("Positive", self.color_plus), ("Negative", self.color_minus)]
        for r, (text, widget) in enumerate(rows):
            self.options.attach(Gtk.Label(label=text, xalign=0), 0, r, 1, 1)
            self.options.attach(widget, 1, r, 1, 1)
        if kind != "cube":
            r = len(rows)
            self.options.attach(Gtk.Label(label="Frames", xalign=0), 0, r, 1, 1)
            fb = Gtk.Box(spacing=6)
            for w in (self.radio_current, self.radio_all):
                if w.get_parent() is not None:
                    w.get_parent().remove(w)
                fb.pack_start(w, False, False, 0)
            self.options.attach(fb, 1, r, 1, 1)
        self.options.show_all()
        self.label_hint.set_markup("<small>{}</small>".format(GLib.markup_escape_text(HINTS[kind])))
        self._show_result_page()

    def _show_result_page(self):
        if not self.Visible:
            return
        self.right_stack.set_visible_child_name("orbitals" if self._current_type() == "orbital" else "info")

    # ------------------------------------------------------------------ structure
    def _system(self):
        sid = self.combobox_systems.get_system_id()
        return (sid, self.p_session.psystem.get(sid)) if sid is not None else (None, None)

    def _vobject(self):
        vid = self.coordinates_combobox.get_vobject_id()
        return self.vm_session.vm_objects_dic.get(vid) if vid is not None else None

    def _on_system_changed(self, _w):
        sid, _system = self._system()
        if sid is None:
            return
        model = self.main.vobject_liststore_dict[sid]
        self.coordinates_combobox.set_model(model)
        # . the last object that is not a surface
        active = 0
        for k, row in enumerate(model):
            v = self.vm_session.vm_objects_dic.get(row[1])
            if v is not None and not getattr(v, "is_surface", False):
                active = k
        self.coordinates_combobox.set_active(active)

    def _on_object_changed(self, _w):
        v = self._vobject()
        n = len(v.frames) if v is not None else 1
        self.spin_frame.set_range(0, max(0, n - 1))
        frame = int(self.vm_session.get_frame() or 0)
        self.spin_frame.set_value(min(frame, n - 1))
        self.radio_all.set_sensitive(n > 1)
        self._on_frame_changed(None)

    def _on_frame_changed(self, _w):
        self.wf = None
        self._fill_orbitals()
        sid, system = self._system()
        if system is None:
            self._status("No system.")
            return
        if system.qcModel is None:
            self._status("This system has no QC model.")
            self.btn_compute.set_sensitive(False)
            return
        if not OA.has_orbitals(system):
            self._status("{}: no orbitals (use a cube file).".format(type(system.qcModel).__name__))
            self.btn_compute.set_sensitive(False)
            return
        self.btn_compute.set_sensitive(True)
        # . NDDO wave functions of small QC regions are computed at once
        n_qc = len(system.qcState.qcAtoms)
        if getattr(system.qcModel, "hamiltonian", None) and n_qc <= 150:
            self.compute_orbitals()
        else:
            self._status("Press Compute Orbitals ({} QC atoms).".format(n_qc))

    def _status(self, text):
        self.label_status.set_text(text)

    # ------------------------------------------------------------------ wave function
    def _copy_system(self, system, vobject, frame):
        from pCore import Clone
        from util.qc_state_repair import repair_qc_state
        copy = Clone(system)
        repair_qc_state(copy)        # . Clone goes through pickle: DFT/HF lose qcState.fitBases
        X = np.asarray(vobject.frames[frame])
        for i in range(len(copy.atoms)):
            for d in range(3):
                copy.coordinates3[i, d] = float(X[i][d])
        return copy

    def compute_orbitals(self):
        if self._busy:
            return
        sid, system = self._system()
        vobject = self._vobject()
        if system is None or vobject is None:
            return
        if len(vobject.atoms) != len(system.atoms):
            self._status("The object does not match the system.")
            return
        frame = int(self.spin_frame.get_value())
        key = (sid, vobject.index, frame)
        if self.wf is not None and self.wf_key == key:
            return
        try:
            copy = self._copy_system(system, vobject, frame)
        except Exception as error:
            self._status("Could not copy the system: {}".format(error))
            return
        self._busy = True
        self._status("Computing the wave function...")
        self.btn_compute.set_sensitive(False)
        t0 = time.time()

        def work():
            try:
                wf = OA.compute_wavefunction(copy)
                GLib.idle_add(self._orbitals_done, wf, key, time.time() - t0, None)
            except Exception as error:
                traceback.print_exc()
                GLib.idle_add(self._orbitals_done, None, key, 0, str(error))
        threading.Thread(target=work, daemon=True).start()

    def _orbitals_done(self, wf, key, seconds, error):
        self._busy = False
        if not self.Visible:
            return False
        self.btn_compute.set_sensitive(True)
        if error:
            self._status("Failed: {}".format(error))
            return False
        self.wf, self.wf_key = wf, key
        homo = wf["n_occupied"] - 1
        self._status("{}, {} basis functions, frame {}: HOMO {:.2f} eV, LUMO {:.2f} eV ({:.1f} s)".format(
            wf["method"], wf["n_basis"], key[2], wf["energies"][homo] * EV,
            wf["energies"][homo + 1] * EV if homo + 1 < len(wf["energies"]) else float("nan"), seconds))
        self._fill_orbitals()
        return False

    def _label(self, k):
        return OA.orbital_label(self.wf, k) if self.wf is not None else str(k)

    def _fill_orbitals(self):
        self.orb_store.clear()
        self.comp_store.clear()
        if self.wf is None:
            self.diagram.clear()
            self.plot.clear()
            self.summary.get_buffer().set_text("")
            return
        wf = self.wf
        E = wf["energies"]
        homo = wf["n_occupied"] - 1
        self.diagram.set_data(E * EV, wf["occupancies"], homo, self._label)
        self._rows = {}
        for k in range(len(E) - 1, -1, -1):
            it = self.orb_store.append([k, self._label(k), "{:g}".format(wf["occupancies"][k]),
                                        "{:.3f}".format(E[k] * EV), "{:.5f}".format(E[k]),
                                        OA.orbital_character(wf, k)])
            self._rows[k] = it
        self._summary()
        self.select_orbital(homo)

    def _summary(self):
        wf = self.wf
        E = wf["energies"]
        homo = wf["n_occupied"] - 1
        lumo = homo + 1
        lines = ["Method                {}".format(wf["method"]),
                 "Frame                 {}".format(self.wf_key[2]),
                 "Basis functions       {}".format(wf["n_basis"]),
                 "QC atoms              {}".format(len(wf["qc_atoms"])),
                 "Electrons             {:g}".format(float(np.sum(wf["occupancies"]))),
                 "Occupied orbitals     {}".format(wf["n_occupied"]),
                 "HOMO                  {:.3f} eV  ({:.5f} hartree)  {}".format(E[homo] * EV, E[homo],
                                                                              OA.orbital_character(wf, homo))]
        if lumo < len(E):
            lines += ["LUMO                  {:.3f} eV  ({:.5f} hartree)  {}".format(E[lumo] * EV, E[lumo],
                                                                                 OA.orbital_character(wf, lumo)),
                      "HOMO-LUMO gap         {:.3f} eV".format((E[lumo] - E[homo]) * EV)]
        lines += ["",
                  "Populations: {}.".format(wf["population_kind"]),
                  ("Gaussian basis: pDynamo numbers the shells of each atom in sequence (6-31G*: 1s core, 2s/2p "
                   "inner and 3s/3p outer valence, 4d polarization); Main contributions add them by angular "
                   "momentum (O3 p)." if wf["population_kind"] == "Mulliken" else
                   "NDDO basis: one valence shell per atom (C 2s 2p; d functions on heavier elements)."),
                  "Contribution of an atomic orbital = its population in the molecular orbital; the "
                  "populations of one orbital add up to 100%. Mulliken populations of single functions "
                  "can be slightly negative.",
                  "",
                  "Orbital energies of semiempirical methods are not quantitative; compare them within "
                  "one method. By Koopmans' theorem -E(HOMO) approximates the ionization energy (HF)."]
        self.summary.get_buffer().set_text("\n".join(lines))

    # ------------------------------------------------------------------ selection
    def select_orbital(self, k):
        if self.wf is None or k is None:
            return
        self.selected = k
        self.diagram.select(k)
        it = self._rows.get(k)
        if it is not None:
            sel = self.orb_view.get_selection()
            if sel.get_selected()[1] is None or self.orb_store[sel.get_selected()[1]][0] != k:
                self._block_row = True
                sel.select_iter(it)
                self.orb_view.scroll_to_cell(self.orb_store.get_path(it), None, False, 0, 0)
                self._block_row = False
        self._fill_composition(k)

    def _on_diagram_select(self, k):
        self.select_orbital(k)

    def _on_orbital_row(self, sel):
        if getattr(self, "_block_row", False):
            return
        model, it = sel.get_selected()
        if it is not None:
            self.select_orbital(model[it][0])

    def _fill_composition(self, k):
        wf = self.wf
        self.comp_store.clear()
        for center, ao, coef, pop in OA.orbital_composition(wf, k):
            if abs(pop) < 0.001 and abs(coef) < 0.01:
                continue
            self.comp_store.append([wf["atom_labels"][center], ao, "{:+.4f}".format(coef),
                                    "{:.1f}".format(100 * pop), center])
        contrib = OA.atom_contributions(wf, k)
        items = sorted(contrib.items(), key=lambda kv: -abs(kv[1]))
        items = [kv for kv in items if abs(kv[1]) >= 0.02][:16] or items[:8]
        self._bar_centers = [c for c, _p in items]
        labels = [wf["atom_labels"][c] for c, _p in items]
        colors = [ELEMENT_COLORS.get(''.join(ch for ch in l if ch.isalpha()), (0.55, 0.35, 0.75)) for l in labels]
        title = "{} ({:.2f} eV): contribution of each atom".format(self._label(k), wf["energies"][k] * EV)
        self.plot.set_lines([{"x": list(range(len(items))), "y": [100 * p for _c, p in items], "style": "bars",
                              "colors": colors, "label": "%"}], "atom", "contribution (%)", title,
                            ymin=0, xlabels=labels)

    def _on_plot_click(self, info):
        k = info.get("index") if isinstance(info, dict) else None
        if k is None or not hasattr(self, "_bar_centers") or not (0 <= k < len(self._bar_centers)):
            return
        self._select_centers([self._bar_centers[k]])

    def _on_comp_row(self, view, path, _col):
        self._select_centers([self.comp_store[path][4]])

    def select_main_atoms(self, threshold=0.10):
        if self.wf is None or self.selected is None:
            return
        contrib = OA.atom_contributions(self.wf, self.selected)
        self._select_centers([c for c, p in contrib.items() if p >= threshold])

    def _select_centers(self, centers):
        v = self._vobject()
        if v is None or self.wf is None:
            return
        indices = [self.wf["qc_atoms"][c] for c in centers]
        self.vm_session._selection_function_set(None)
        sel = self.vm_session.selections[self.vm_session.current_selection]
        mode = sel.selection_mode
        sel.selection_mode = "atom"
        self.vm_session._selection_function_set({v.atoms[i] for i in indices})
        sel.selection_mode = mode
        self.vm_session.vm_glcore.queue_draw()

    # ------------------------------------------------------------------ render
    def _float(self, entry, default=None):
        text = entry.get_text().strip()
        if not text:
            return default
        return float(text)

    def _rgb(self, button):
        c = button.get_rgba()
        return [c.red, c.green, c.blue]

    def render(self, kind=None):
        kind = kind or self._current_type()
        if kind == "orbital" and self._current_type() != "orbital":
            self.type_view.get_selection().select_path(Gtk.TreePath(0))
        sid, system = self._system()
        vobject = self._vobject()
        if system is None or vobject is None:
            return
        try:
            iso = self._float(self.entry_iso)
            spacing = self._float(self.entry_spacing, 0.4)
            mep_min = self._float(self.entry_mep_min)
            mep_max = self._float(self.entry_mep_max)
            pot_spacing = self._float(self.entry_pot_spacing)
        except ValueError:
            self.main.simple_dialog.info(msg="Invalid number in the options.", title="Orbitals and Surfaces")
            return
        if iso is None:
            self.main.simple_dialog.info(msg="Give the isovalue.", title="Orbitals and Surfaces")
            return
        params = {"_GridSpacing": spacing, "_isovalue": iso, "_IsosurfaceTag": "Isosurface",
                  "color_plus": self._rgb(self.color_plus), "color_minus": self._rgb(self.color_minus),
                  "orbital_key": 0}
        if kind == "orbital":
            if self.wf is None or self.selected is None:
                self.main.simple_dialog.info(msg="Compute the orbitals and choose one first.",
                                             title="Orbitals and Surfaces")
                return
            params.update({"type": "orbital", "_OrbitalTag": "Grid Orbitals", "orbital_key": int(self.selected)})
            name = self._label(self.selected)
        elif kind == "density":
            color = self._rgb(self.color_density)
            params.update({"type": "density", "_OrbitalTag": "density", "_IsosurfaceTag": "density",
                           "color_plus": color, "color_minus": color})
            name = "Density"
        elif kind == "potential":
            params.update({"type": "potential", "_OrbitalTag": "potential"})
            name = "Potential"
        elif kind == "mep":
            params.update({"type": "mep", "_OrbitalTag": "density_mep", "mep_vmin": mep_min, "mep_vmax": mep_max,
                           "mep_cmap_name": self.combo_cmap.get_active_text() or "coolwarm",
                           "mep_pot_spacing": pot_spacing})
            name = "MEP"
        else:
            if not self.cube_path:
                self.main.simple_dialog.info(msg="Choose the cube file first.", title="Orbitals and Surfaces")
                return
            params.update({"type": "external_cube", "external_density_path": self.cube_path,
                           "external_potential_path": None, "mep_vmin": None, "mep_vmax": None,
                           "mep_cmap_name": "coolwarm"})
            name = os.path.splitext(os.path.basename(self.cube_path))[0]
        if kind != "cube" and not OA.has_orbitals(system):
            self.main.simple_dialog.info(msg="This QC model gives no wave function to EasyHybrid: use a cube file.",
                                         title="Orbitals and Surfaces")
            return
        n_frames = len(vobject.frames)
        if kind == "cube":
            frames = [0]
        elif self.radio_all.get_active() and n_frames > 1:
            frames = list(range(n_frames))
        else:
            frames = [int(self.spin_frame.get_value())]
            if n_frames > 1:
                name += " (frame {})".format(frames[0])
        self._run_render(kind, name, params, system, vobject, frames)

    def _run_render(self, kind, name, params, system, vobject, frames):
        from gui.windows.analysis.surface_analysis_window import generate_grid_parallel
        self.btn_render.set_sensitive(False)
        self.progress.set_text("computing {} frame(s)...".format(len(frames)))
        self.progress.pulse()
        t0 = time.time()
        if kind == "cube":
            jobs = [[0, None, None, dict(params)]]
        else:
            jobs = []
            for f in frames:
                coords = self.p_session.get_coordinates_from_vobject(vobject=vobject, frame=f)
                jobs.append([f, system, coords, dict(params)])
        if len(jobs) == 1:
            # . one frame: in a thread, on a copy of the system
            job = jobs[0]
            if job[1] is not None:
                from pCore import Clone
                from util.qc_state_repair import repair_qc_state
                job[1] = Clone(system)
                repair_qc_state(job[1])

            def work():
                try:
                    r = generate_grid_parallel(job)
                    GLib.idle_add(self._render_done, kind, name, params, vobject, frames, [r], time.time() - t0, None)
                except Exception as error:
                    traceback.print_exc()
                    GLib.idle_add(self._render_done, kind, name, params, vobject, frames, None, 0, str(error))
            threading.Thread(target=work, daemon=True).start()
            self._pulse_id = GLib.timeout_add(150, self._pulse)
            return
        # . several frames: a pool of processes (blocks the window, as before)
        backup = (getattr(system, "e_treeview_iter", None), getattr(system, "e_liststore_iter", None))
        system.e_treeview_iter = None
        system.e_liststore_iter = None
        try:
            while Gtk.events_pending():
                Gtk.main_iteration()
            with multiprocessing.Pool(processes=max(1, min(MAX_WORKERS, multiprocessing.cpu_count(), len(jobs)))) as pool:
                results = pool.map(generate_grid_parallel, jobs)
            error = None
        except Exception as e:
            traceback.print_exc()
            results, error = None, str(e)
        finally:
            system.e_treeview_iter, system.e_liststore_iter = backup
        self._render_done(kind, name, params, vobject, frames, results, time.time() - t0, error)

    def _pulse(self):
        if self.btn_render.get_sensitive() or not self.Visible:
            return False
        self.progress.pulse()
        return True

    def _render_done(self, kind, name, params, vobject, frames, results, seconds, error):
        if not self.Visible:
            return False
        self.btn_render.set_sensitive(True)
        self.progress.set_fraction(0)
        if error:
            self.progress.set_text("failed")
            self.main.simple_dialog.info(msg="The surface could not be computed:\n\n{}".format(error),
                                         title="Orbitals and Surfaces")
            return False
        self.progress.set_text("done in {:.1f} s".format(seconds))
        obj = add_surface_object(self.main, vobject, name, kind, params, results, frames)
        self._report(kind, name, params, results, seconds)
        return False

    def _report(self, kind, name, params, results, seconds):
        r = results[0] if results else {}
        lines = ["{}  ({:.1f} s, {} frame(s))".format(name, seconds, len(results)), ""]
        for key, surf in r.items():
            if isinstance(surf, list) and len(surf) >= 3:
                try:
                    lines.append("{:<14s} {} vertices".format(key.replace("obital", "lobe"), len(surf[0]) // 3))
                except Exception:
                    pass
        if kind == "mep" and r.get("mep_raw_potential_values") is not None:
            v = np.asarray(r["mep_raw_potential_values"], dtype=float)
            lines += ["", "Potential on the surface (a.u.): min {:.4f}, max {:.4f}".format(v.min(), v.max()),
                      "                       (kJ/mol/e): min {:.1f}, max {:.1f}".format(v.min() * 2625.5, v.max() * 2625.5)]
        lines += ["", "Isovalue {:g}, grid spacing {:g} bohr.".format(params["_isovalue"], params.get("_GridSpacing", 0)),
                  "The surface is an object under its structure in the Objects tree: right-click it for "
                  "Setup (opacity, wireframe, colors...) and Delete."]
        self.info.get_buffer().set_text("\n".join(lines))

    # ------------------------------------------------------------------ files
    def on_export_image(self):
        dialog = Gtk.FileChooserDialog(title="Export Image", parent=self.window, action=Gtk.FileChooserAction.SAVE)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, "Diagram", 1, "Contributions", 2)
        dialog.set_do_overwrite_confirmation(True)
        dialog.set_current_name("orbitals.png")
        response = dialog.run()
        path = dialog.get_filename()
        dialog.destroy()
        if not path or response not in (1, 2):
            return
        if not path.lower().endswith((".png", ".svg")):
            path += ".png"
        (self.diagram if response == 1 else self.plot).export(path)

    def on_save_table(self):
        if self.wf is None:
            return
        dialog = Gtk.FileChooserDialog(title="Save Table", parent=self.window, action=Gtk.FileChooserAction.SAVE)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        dialog.set_do_overwrite_confirmation(True)
        dialog.set_current_name("orbitals.txt")
        path = dialog.get_filename() if dialog.run() == Gtk.ResponseType.OK else None
        dialog.destroy()
        if path:
            self.save_table(path)

    def save_table(self, path):
        wf = self.wf
        with open(path, "w") as handle:
            handle.write("# {}  frame {}  populations: {}\n".format(wf["method"], self.wf_key[2], wf["population_kind"]))
            handle.write("# index\torbital\toccupancy\tE_eV\tE_hartree\tmain_contributions\n")
            for k in range(len(wf["energies"])):
                handle.write("{}\t{}\t{:g}\t{:.4f}\t{:.6f}\t{}\n".format(k, self._label(k), wf["occupancies"][k],
                                                                         wf["energies"][k] * EV, wf["energies"][k],
                                                                         OA.orbital_character(wf, k)))
            if self.selected is not None:
                handle.write("\n# LCAO of {} (index {})\n# atom\tatomic_orbital\tcoefficient\tpopulation\n".format(
                    self._label(self.selected), self.selected))
                for center, ao, coef, pop in OA.orbital_composition(wf, self.selected):
                    handle.write("{}\t{}\t{:+.6f}\t{:.4f}\n".format(wf["atom_labels"][center], ao, coef, pop))

    def _on_destroy(self, *_):
        self.Visible = False
        self.window = None

    def update(self, *args, **kwargs):
        """ Called by the main window when systems or objects change. """
        if not self.Visible:
            return
        sid = self.combobox_systems.get_system_id()
        if sid is not None and sid in self.main.vobject_liststore_dict:
            self.coordinates_combobox.set_model(self.main.vobject_liststore_dict[sid])


# =============================================================================
#  Surface objects (as the previous window made them)
# =============================================================================
def add_surface_object(main, vismol_object, name, kind, params, results, frames):
    from vismol.libgl.representations import SurfaceRepresentation
    from vismol.core.vismol_object import VismolObject
    vm_session = main.vm_session
    obj = VismolObject(name=name, index=-1, vismol_session=vm_session, trajectory=[], bonds_pair_of_indexes=[0, 1])
    obj.model_mat = vismol_object.model_mat
    obj.trans_mat = vismol_object.trans_mat
    obj.surface_trajectory = results
    two_lobes = kind in ("orbital", "potential") or (kind == "cube" and results and "obital_minus" in results[0])
    obj.representations["surface1"] = SurfaceRepresentation(
        vismol_object=obj, vismol_glcore=vm_session.vm_glcore, name="surface", active=True, indexes=[],
        is_dynamic=False, iso_color=params["color_plus"], surface_name="obital_plus")
    if two_lobes:
        obj.representations["surface2"] = SurfaceRepresentation(
            vismol_object=obj, vismol_glcore=vm_session.vm_glcore, name="surface", active=True, indexes=[],
            is_dynamic=False, iso_color=params["color_minus"], surface_name="obital_minus")
    obj.parameters = params
    if len(frames) == len(vismol_object.frames):
        obj.frames = vismol_object.frames
    else:
        obj.frames = np.asarray(vismol_object.frames)[frames]
    obj.active = False
    obj.is_surface = True
    obj.surface_type = {"orbital": "orbital", "density": "density", "potential": "potential",
                        "mep": "mep", "cube": "external"}[kind]
    if kind == "mep":
        obj.mep_potential_values = [r.get("mep_raw_potential_values") for r in results]
        obj.mep_polygons = [r.get("mep_raw_polygons") for r in results]
        obj.mep_vmin = params.get("mep_vmin")
        obj.mep_vmax = params.get("mep_vmax")
        obj.mep_cmap_name = params.get("mep_cmap_name", "coolwarm")
    obj.e_id = vismol_object.e_id
    vm_session._add_vismol_object(obj, show_molecule=False, autocenter=False)
    main.main_treeview.add_vismol_object_to_treeview(obj, vismol_object.e_treeview_iter)
    main.add_vobject_to_vobject_liststore_dict(obj)
    main.refresh_widgets()
    obj.active = True
    try:
        main.main_treeview.treestore.set_value(obj.e_treeview_iter, 6, True)
    except Exception:
        pass
    vm_session.vm_glcore.queue_draw()
    return obj
