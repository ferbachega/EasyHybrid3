"""
ir_spectrum_window
==================

IR spectrum of the imported normal modes (2026-10-24), opened by Plot IR
Spectrum in Analysis > Normal Modes. Options on the left, plot on the right:

  frequency scale factor        harmonic frequencies are too high (B3LYP ~0.96)
  band shape and width          Lorentzian or Gaussian, FWHM in cm-1
  range, wavenumber decreasing  4000 -> 400 cm-1 as in experimental spectra
  y axis                        intensity (km/mol), normalized absorbance or
                                transmittance (%, peaks down)
  modes                         the discrete modes as grey lines (on/off)
  peak labels                   wavenumber of the N strongest peaks

The band height is the intensity of the mode (as before). A click near a mode
selects it in the Normal Modes window (on_mode_clicked); the value under the
mouse is shown in the corner.
"""

import math

import cairo
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk

import numpy as np

from util.ir_spectrum import broaden_spectrum

Y_MODES = [("intensity", "Intensity (km/mol)"), ("absorbance", "Absorbance (normalized)"),
           ("transmittance", "Transmittance (%)")]


def _nice_step(span, target=8):
    raw = span / max(target, 1)
    mag = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1
    for f in (1, 2, 2.5, 5, 10):
        if raw <= f * mag:
            return f * mag
    return 10 * mag


class IRSpectrumPlot(Gtk.DrawingArea):

    def __init__(self):
        super().__init__()
        self.data = None            # . dict built by IRSpectrumWindow.update()
        self.hover = None
        self.on_click = None
        self._geom = None
        self.add_events(Gdk.EventMask.POINTER_MOTION_MASK | Gdk.EventMask.LEAVE_NOTIFY_MASK |
                        Gdk.EventMask.BUTTON_PRESS_MASK)
        self.connect("draw", self._on_draw)
        self.connect("motion-notify-event", self._on_motion)
        self.connect("leave-notify-event", lambda *_: self._set_hover(None))
        self.connect("button-press-event", self._on_button)
        self.set_size_request(560, 340)

    def set_data(self, data):
        self.data = data
        self.queue_draw()

    # ------------------------------------------------------------------ geometry
    def _X(self, v):
        x0, y0, w, h, a, b, c, d = self._geom
        t = (v - a) / (b - a)
        if self.data["reversed"]:
            t = 1 - t
        return x0 + t * w

    def _Y(self, v):
        x0, y0, w, h, a, b, c, d = self._geom
        return y0 + h - (v - c) / (d - c) * h

    def _wavenumber_at(self, px):
        x0, y0, w, h, a, b, c, d = self._geom
        t = (px - x0) / w
        if self.data["reversed"]:
            t = 1 - t
        return a + t * (b - a)

    # ------------------------------------------------------------------ drawing
    def _on_draw(self, widget, cr, width=None, height=None):
        width = width or widget.get_allocated_width()
        height = height or widget.get_allocated_height()
        cr.set_source_rgb(1, 1, 1)
        cr.paint()
        cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        D = self.data
        if not D:
            return False
        x0, y0 = 70, 36
        w, h = width - x0 - 24, height - y0 - 56
        a, b = D["xmin"], D["xmax"]
        c, d = D["ymin"], D["ymax"]
        self._geom = (x0, y0, w, h, a, b, c, d)
        # . grid and ticks
        cr.set_font_size(10.5); cr.set_line_width(1)
        step = _nice_step(b - a)
        v = math.ceil(a / step) * step
        while v <= b + 1e-9:
            x = self._X(v)
            cr.set_source_rgb(0.92, 0.92, 0.92); cr.move_to(x, y0); cr.line_to(x, y0 + h); cr.stroke()
            cr.set_source_rgb(0.2, 0.2, 0.2); t = "{:g}".format(v); e = cr.text_extents(t)
            cr.move_to(x - e.width / 2, y0 + h + 15); cr.show_text(t)
            v += step
        step = _nice_step(d - c, 6)
        v = math.ceil(c / step) * step
        while v <= d + 1e-9:
            y = self._Y(v)
            cr.set_source_rgb(0.92, 0.92, 0.92); cr.move_to(x0, y); cr.line_to(x0 + w, y); cr.stroke()
            cr.set_source_rgb(0.2, 0.2, 0.2); t = "{:g}".format(round(v, 6)); e = cr.text_extents(t)
            cr.move_to(x0 - e.width - 6, y + e.height / 2); cr.show_text(t)
            v += step
        cr.set_source_rgb(0.3, 0.3, 0.3); cr.rectangle(x0, y0, w, h); cr.stroke()
        # . titles
        cr.set_source_rgb(0.1, 0.1, 0.1); cr.set_font_size(12)
        t = "wavenumber (cm-1)"; e = cr.text_extents(t)
        cr.move_to(x0 + w / 2 - e.width / 2, y0 + h + 38); cr.show_text(t)
        cr.save(); t = D["ylabel"]; e = cr.text_extents(t)
        cr.move_to(18, y0 + h / 2 + e.width / 2); cr.rotate(-math.pi / 2); cr.show_text(t); cr.restore()
        cr.set_font_size(13); t = D["title"]; e = cr.text_extents(t)
        cr.move_to(x0 + w / 2 - e.width / 2, 22); cr.show_text(t)
        # . data, clipped to the frame
        cr.save(); cr.rectangle(x0, y0, w, h); cr.clip()
        if D["show_sticks"]:
            cr.set_source_rgb(0.62, 0.62, 0.62); cr.set_line_width(1.2)
            base = self._Y(D["stick_base"])
            for f, v, _k in D["sticks"]:
                if a <= f <= b:
                    x = self._X(f); cr.move_to(x, base); cr.line_to(x, self._Y(v)); cr.stroke()
        cr.set_source_rgb(*D["color"]); cr.set_line_width(1.8)
        first = True
        for xv, yv in zip(D["x"], D["y"]):
            if first:
                cr.move_to(self._X(xv), self._Y(yv)); first = False
            else:
                cr.line_to(self._X(xv), self._Y(yv))
        cr.stroke()
        # . selected mode
        if D.get("selected") is not None:
            f = D["selected"]
            if a <= f <= b:
                cr.set_source_rgba(1.0, 0.75, 0.0, 0.8); cr.set_line_width(2.5)
                cr.move_to(self._X(f), y0); cr.line_to(self._X(f), y0 + h); cr.stroke()
        cr.restore()
        # . peak labels (no two closer than 34 px)
        cr.set_font_size(10); cr.set_source_rgb(0.15, 0.15, 0.15)
        placed = []
        for f, v in D["peaks"]:
            if not a <= f <= b:
                continue
            x = self._X(f)
            if any(abs(x - p) < 34 for p in placed):
                continue
            placed.append(x)
            y = self._Y(v)
            t = "{:.0f}".format(f); e = cr.text_extents(t)
            ty = y - 6 if not D["peaks_down"] else y + e.height + 6
            ty = min(max(ty, y0 + e.height + 2), y0 + h - 3)
            cr.move_to(x - e.width / 2, ty); cr.show_text(t)
        # . hover readout
        if self.hover:
            t = self.hover
            cr.set_font_size(11); e = cr.text_extents(t)
            bx, by = x0 + w - e.width - 14, y0 + 8
            cr.set_source_rgba(1, 1, 1, 0.92); cr.rectangle(bx - 5, by - 2, e.width + 12, e.height + 9); cr.fill()
            cr.set_source_rgb(0.3, 0.3, 0.3); cr.set_line_width(0.6); cr.rectangle(bx - 5, by - 2, e.width + 12, e.height + 9); cr.stroke()
            cr.set_source_rgb(0.1, 0.1, 0.1); cr.move_to(bx, by + e.height + 2); cr.show_text(t)
        return False

    # ------------------------------------------------------------------ mouse
    def _nearest_mode(self, px):
        D = self.data
        best, dist = None, 8.0
        for f, v, k in D["sticks"]:
            if D["xmin"] <= f <= D["xmax"]:
                dx = abs(self._X(f) - px)
                if dx < dist:
                    best, dist = (f, v, k), dx
        return best

    def _on_motion(self, _w, event):
        if not self.data or self._geom is None:
            return False
        x0, y0, w, h = self._geom[:4]
        if not (x0 <= event.x <= x0 + w and y0 <= event.y <= y0 + h):
            return self._set_hover(None)
        f = self._wavenumber_at(event.x)
        text = "{:.0f} cm-1".format(f)
        mode = self._nearest_mode(event.x)
        if mode is not None:
            text = "mode {}: {:.1f} cm-1, {:.2f} km/mol".format(mode[2], mode[0], self.data["raw"][mode[2]])
        return self._set_hover(text)

    def _set_hover(self, text):
        if text != self.hover:
            self.hover = text
            self.queue_draw()
        return False

    def _on_button(self, _w, event):
        if event.button != 1 or not self.data or self._geom is None:
            return False
        mode = self._nearest_mode(event.x)
        if mode is not None and self.on_click:
            self.on_click(mode[2], mode[0])
        return True

    def export(self, path, width=1100, height=600):
        surface = cairo.SVGSurface(path, width, height) if path.lower().endswith(".svg") else \
            cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
        cr = cairo.Context(surface)
        hover, geom = self.hover, self._geom
        self.hover = None
        self._on_draw(self, cr, width, height)
        self.hover, self._geom = hover, geom
        if path.lower().endswith(".svg"):
            surface.finish()
        else:
            surface.write_to_png(path)


class IRSpectrumWindow:
    """ modes: [(mode index, frequency cm-1, intensity km/mol), ...] """

    def __init__(self, main, modes, title="", parent=None, on_mode_clicked=None):
        self.main = main
        self.modes = [(int(k), float(f), float(i)) for k, f, i in modes]
        self.on_mode_clicked = on_mode_clicked
        self.selected = None
        self.window = Gtk.Window(title="IR Spectrum" + (" - " + title if title else ""))
        if parent is not None:
            self.window.set_transient_for(parent)
        self.window.set_default_size(1060, 560)
        outer = Gtk.Box(spacing=8, margin=8)
        self.window.add(outer)
        outer.pack_start(self._build_options(), False, False, 0)
        right = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
        self.plot = IRSpectrumPlot()
        self.plot.on_click = self._on_plot_click
        right.pack_start(self.plot, True, True, 0)
        bar = Gtk.Box(spacing=6)
        self.label_info = Gtk.Label(xalign=0)
        self.label_info.set_ellipsize(3)
        bar.pack_start(self.label_info, True, True, 0)
        b = Gtk.Button(label="Save Data...")
        b.set_tooltip_text("The spectrum curve and the modes (scaled wavenumber, intensity) as text.")
        b.connect("clicked", self.on_save_data)
        bar.pack_end(b, False, False, 0)
        b = Gtk.Button(label="Export Image...")
        b.connect("clicked", self.on_export_image)
        bar.pack_end(b, False, False, 0)
        right.pack_start(bar, False, False, 0)
        outer.pack_start(right, True, True, 0)
        self.window.show_all()
        self.update()

    # ------------------------------------------------------------------ options
    def _spin(self, lo, hi, step, value, digits=0, tip=None):
        s = Gtk.SpinButton.new_with_range(lo, hi, step)
        s.set_digits(digits)
        s.set_value(value)
        if tip:
            s.set_tooltip_text(tip)
        s.connect("value-changed", lambda *_: self.update())
        return s

    def _check(self, label, active, tip=None):
        c = Gtk.CheckButton(label=label)
        c.set_active(active)
        if tip:
            c.set_tooltip_text(tip)
        c.connect("toggled", lambda *_: self.update())
        return c

    def _build_options(self):
        box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
        box.set_size_request(270, -1)

        def group(title, rows):
            frame = Gtk.Frame(label=" {} ".format(title))
            grid = Gtk.Grid(column_spacing=8, row_spacing=6, margin=8)
            for r, row in enumerate(rows):
                if isinstance(row, tuple):
                    grid.attach(Gtk.Label(label=row[0], xalign=0), 0, r, 1, 1)
                    grid.attach(row[1], 1, r, 1, 1)
                else:
                    grid.attach(row, 0, r, 2, 1)
            frame.add(grid)
            box.pack_start(frame, False, False, 0)

        self.spin_scale = self._spin(0.5, 1.5, 0.001, 1.0, 3,
                                     "Factor applied to every frequency. Harmonic frequencies are too high: "
                                     "about 0.96 for B3LYP, 0.90 for HF; semiempirical methods vary.")
        self.combo_shape = Gtk.ComboBoxText()
        for key, text in (("lorentzian", "Lorentzian"), ("gaussian", "Gaussian")):
            self.combo_shape.append(key, text)
        self.combo_shape.set_active_id("lorentzian")
        self.combo_shape.connect("changed", lambda *_: self.update())
        self.spin_fwhm = self._spin(1, 300, 1, 20, 0, "Full width at half maximum of each band (cm-1).")
        group("Bands", [("Scale factor", self.spin_scale), ("Shape", self.combo_shape),
                        ("Width (cm-1)", self.spin_fwhm)])

        self.spin_xmin = self._spin(0, 10000, 50, 400, 0, "Lowest wavenumber shown (cm-1).")
        self.spin_xmax = self._spin(100, 10000, 50, 4000, 0, "Highest wavenumber shown (cm-1).")
        self.check_reversed = self._check("Wavenumber decreasing (high to low)", True,
                                          "The usual orientation of experimental IR spectra.")
        fit = Gtk.Button(label="Fit to Modes")
        fit.set_tooltip_text("Range from the lowest to the highest mode, with a margin.")
        fit.connect("clicked", self.on_fit_range)
        group("Wavenumber axis", [("From", self.spin_xmin), ("To", self.spin_xmax), self.check_reversed, fit])

        self.combo_y = Gtk.ComboBoxText()
        for key, text in Y_MODES:
            self.combo_y.append(key, text)
        self.combo_y.set_active_id("intensity")
        self.combo_y.set_tooltip_text("Transmittance: T = 10^(-A), with the strongest band at A = 1 (10% T).")
        self.combo_y.connect("changed", lambda *_: self.update())
        self.spin_ymax = self._spin(0, 1e6, 10, 0, 1, "Top of the intensity axis (km/mol); 0 = automatic.")
        group("Intensity axis", [("Show", self.combo_y), ("Max (0 = auto)", self.spin_ymax)])

        self.check_sticks = self._check("Modes (grey lines)", True, "One line per mode, as high as its intensity.")
        self.check_labels = self._check("Label the strongest peaks", True)
        self.spin_labels = self._spin(1, 50, 1, 6, 0, "Number of peaks labelled with their wavenumber.")
        self.color = Gtk.ColorButton.new_with_rgba(Gdk.RGBA(0.85, 0.30, 0.10, 1.0))
        self.color.connect("color-set", lambda *_: self.update())
        group("Display", [self.check_sticks, self.check_labels, ("Peaks", self.spin_labels), ("Color", self.color)])

        note = Gtk.Label(xalign=0)
        note.set_line_wrap(True)
        note.set_max_width_chars(36)
        note.set_markup("<small>Click a mode in the plot to select it in the Normal Modes window "
                        "(then Play animates it). Band height = intensity of the mode.</small>")
        box.pack_start(note, False, False, 0)
        return box

    # ------------------------------------------------------------------ spectrum
    def _spectrum(self):
        scale = self.spin_scale.get_value()
        k = [m[0] for m in self.modes]
        f = np.array([m[1] for m in self.modes]) * scale
        I = np.array([m[2] for m in self.modes])
        xmin, xmax = sorted((self.spin_xmin.get_value(), self.spin_xmax.get_value()))
        if xmax - xmin < 10:
            xmax = xmin + 10
        n = int(min(20000, max(2000, (xmax - xmin) / max(self.spin_fwhm.get_value(), 1) * 40)))
        x, y = broaden_spectrum(f, I, x_min=xmin, x_max=xmax, num_points=n, fwhm=self.spin_fwhm.get_value(),
                                lineshape=self.combo_shape.get_active_id() or "lorentzian")
        # . the discrete modes: real vibrations only (no translations/rotations or imaginary modes)
        mask = f >= 50.0 * scale
        sf, si = f[mask], I[mask]
        sticks_k = [kk for kk, keep in zip(k, mask) if keep]
        return x, y, sf, si, sticks_k, xmin, xmax

    def update(self, *_):
        if not hasattr(self, "plot"):
            return
        x, y, sf, si, sk, xmin, xmax = self._spectrum()
        mode = self.combo_y.get_active_id() or "intensity"
        top = max(float(y.max()) if len(y) else 0.0, float(si.max()) if len(si) else 0.0, 1e-12)
        if mode == "intensity":
            Y, S = y, si
            ymax = self.spin_ymax.get_value() or top * 1.12
            ymin, base, ylabel, down = 0.0, 0.0, "intensity (km/mol)", False
        elif mode == "absorbance":
            Y, S = y / top, si / top
            ymin, ymax, base, ylabel, down = 0.0, 1.12, 0.0, "absorbance (normalized)", False
        else:
            Y, S = 100 * 10 ** (-y / top), 100 * 10 ** (-si / top)
            ymin, ymax, base, ylabel, down = 0.0, 105.0, 100.0, "transmittance (%)", True
        self.spin_ymax.set_sensitive(mode == "intensity")
        peaks = []
        if self.check_labels.get_active() and len(y) > 2:
            idx = np.where((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:]) & (y[1:-1] > 0.02 * top))[0] + 1
            idx = sorted(idx, key=lambda i: -y[i])[:int(self.spin_labels.get_value())]
            peaks = [(float(x[i]), float(Y[i])) for i in idx]
        c = self.color.get_rgba()
        raw = {m[0]: m[2] for m in self.modes}
        self.plot.set_data({
            "x": x, "y": Y, "xmin": xmin, "xmax": xmax, "ymin": ymin, "ymax": ymax,
            "reversed": self.check_reversed.get_active(), "ylabel": ylabel,
            "title": "IR spectrum ({}, FWHM {:g} cm-1{})".format(
                (self.combo_shape.get_active_text() or "Lorentzian").lower(), self.spin_fwhm.get_value(),
                ", scaled x{:.3f}".format(self.spin_scale.get_value()) if abs(self.spin_scale.get_value() - 1) > 1e-6 else ""),
            "show_sticks": self.check_sticks.get_active(), "sticks": list(zip(sf, S, sk)), "stick_base": base,
            "peaks": peaks, "peaks_down": down, "color": (c.red, c.green, c.blue),
            "selected": self.selected, "raw": raw})
        strongest = max(self.modes, key=lambda m: m[2]) if self.modes else None
        self.label_info.set_text("{} modes with intensity{}".format(
            len(sf), "; strongest: mode {} at {:.0f} cm-1 ({:.1f} km/mol)".format(
                strongest[0], strongest[1] * self.spin_scale.get_value(), strongest[2]) if strongest else ""))

    def on_fit_range(self, *_):
        f = [m[1] * self.spin_scale.get_value() for m in self.modes if m[1] >= 50]
        if not f:
            return
        lo = max(0, math.floor((min(f) - 150) / 50) * 50)
        if min(f) > 450:                     # . the usual lower end of a mid-IR spectrum
            lo = 400
        hi = math.ceil((max(f) + 150) / 50) * 50
        self.spin_xmin.set_value(lo)
        self.spin_xmax.set_value(hi)

    def _on_plot_click(self, mode, wavenumber):
        self.selected = wavenumber
        self.update()
        if self.on_mode_clicked:
            self.on_mode_clicked(mode)

    def select_mode(self, mode):
        for k, f, _i in self.modes:
            if k == mode:
                self.selected = f * self.spin_scale.get_value()
                self.update()
                return

    # ------------------------------------------------------------------ files
    def _save_dialog(self, title, name):
        dialog = Gtk.FileChooserDialog(title=title, parent=self.window, action=Gtk.FileChooserAction.SAVE)
        dialog.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK)
        dialog.set_do_overwrite_confirmation(True)
        dialog.set_current_name(name)
        path = dialog.get_filename() if dialog.run() == Gtk.ResponseType.OK else None
        dialog.destroy()
        return path

    def on_export_image(self, *_):
        path = self._save_dialog("Export Image", "ir_spectrum.png")
        if path:
            if not path.lower().endswith((".png", ".svg")):
                path += ".png"
            self.plot.export(path)

    def on_save_data(self, *_):
        path = self._save_dialog("Save Data", "ir_spectrum.txt")
        if path:
            self.save_data(path)

    def save_data(self, path):
        D = self.plot.data
        with open(path, "w") as handle:
            handle.write("# {}\n# wavenumber_cm-1\t{}\n".format(D["title"], D["ylabel"].replace(" ", "_")))
            for xv, yv in zip(D["x"], D["y"]):
                handle.write("{:.2f}\t{:.6g}\n".format(xv, yv))
            handle.write("\n# modes: index\tscaled_wavenumber_cm-1\tintensity_km/mol\n")
            for k, f, i in self.modes:
                handle.write("{}\t{:.2f}\t{:.4f}\n".format(k, f * self.spin_scale.get_value(), i))
